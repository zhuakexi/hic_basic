"""Schicluster single-cell Hi-C imputation, GPU-accelerated (torch).

Ported from the production implementation in the ``hemsc`` repository
(``notebooks/5_imputation.py``), which produced the reference data for the
8-cell / TE / ICM / Morula mouse embryo runs.  The parameter surface and
semantics of every function here match the source exactly, so results
continue to agree with the GPU-generated reference data in
``tests/data/schicluster_gpu/``.

Algorithm summary (details in ``issues/schicluster_gpu/alg.md``):

1. Gaussian smoothing of the raw contact matrix (conv2d, zero padding);
2. Random Walk with Restart (RWR): solve ``(I - (1-alpha) P^T) R = alpha I``
   with a batched ``torch.linalg.solve`` on GPU (``scipy.linalg.solve`` on
   CPU), then symmetrize;
3. Post-processing: zero the main diagonal, full-sum normalization, and an
   upper-band distance constraint (keep main diagonal + ``upper_band``
   super-diagonals).

Production parameters (10 kb resolution):

``alpha=0.05, kernel_size=3, sigma=2, if_convolve=True, ignore_diag=True,
norm=True, band_width=4_000_000 -> upper_band=400, convolve_dtype=torch.int32,
kernel_dtype=torch.float32, rwr_dtype=torch.float32``

The multiprocess pipeline (``hic_imputation_multiprocess``) keeps the GPU
saturated: ``nproc_loading`` loader processes fetch matrices from cool files
into a queue, a single GPU thread imputes (CUDA cannot be used in forked
subprocesses), and ``nproc_saving`` saver processes write compressed CSR
``.npz`` files.
"""
import queue
import threading
import time
from multiprocessing import JoinableQueue, Process

import cooler
import numpy as np
import torch
from scipy import ndimage
from scipy.linalg import solve as scipy_solve
from scipy.signal import windows
from scipy.sparse import csr_matrix

# --- RWR kernels --- #


def solve_rwr_inverse(stoch_matrix, alpha=0.05):
    """Solve the Random Walk with Restart (RWR) equations by matrix inverse.

    Solves ``(I - (1-alpha) * P^T) * R = alpha * I`` for the per-source
    restart distributions and symmetrizes the result.

    Args:
        stoch_matrix (np.ndarray): column-stochastic matrix ``P``, shape
            ``(N, N)``.
        alpha (float): restart probability.

    Returns:
        np.ndarray: symmetrized RWR solution, shape ``(N, N)``.
    """
    m = stoch_matrix * (1 - alpha)
    m = m.transpose()
    y = np.eye(m.shape[0])
    A = y - m

    s = scipy_solve(A, y)
    s *= alpha
    s += s.transpose()
    return s


def solve_rwr_inverse_gpu(stoch_matrix, alpha=0.05, rwr_dtype=torch.float32, device=None):
    """GPU version of :func:`solve_rwr_inverse`, batched over ``(B, C, N, N)``.

    A single batched ``torch.linalg.solve`` replaces the per-matrix
    ``scipy.linalg.solve`` of the CPU reference; this is the main speedup of
    the GPU path.

    Args:
        stoch_matrix: column-stochastic matrix or tensor, shape
            ``(batch_size, channels, N, N)``.
        alpha (float): restart probability.
        rwr_dtype (torch.dtype): arithmetic dtype for the solve.
        device (torch.device): device to run on; defaults to CUDA if
            available, else CPU.

    Returns:
        torch.Tensor: symmetrized RWR solution, shape ``(B, C, N, N)``.
    """
    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if isinstance(stoch_matrix, torch.Tensor):
        m = stoch_matrix.to(device=device, dtype=rwr_dtype)
    else:
        m = torch.tensor(stoch_matrix, device=device, dtype=rwr_dtype)

    m = m * (1 - alpha)
    m = m.transpose(-1, -2)
    # I in shape (batch_size, channels, N, N)
    y = torch.eye(m.shape[2], device=device, dtype=m.dtype).unsqueeze(0).unsqueeze(0)

    # solve: (I - M^T) x = I
    A = y - m
    s = torch.linalg.solve(A, y)

    s = s * alpha
    s = s + s.transpose(-1, -2)
    return s


# --- single-matrix imputation --- #


def schicluster_imputation_cpu(mat, alpha=0.05, kernel_size=3, sigma=2, if_convolve=True,
                               ignore_diag=True, oldway=False, norm=True, upper_band=None,
                               convolve_dtype=np.float64, kernel_dtype=np.float64,
                               rwr_dtype=np.float64, if_rwr=True):
    """Schicluster imputation of one contact matrix, NumPy reference.

    Same parameter surface as the GPU version, so the two paths can be
    cross-validated.  Note: with ``ignore_diag=True`` the input ``mat`` is
    modified in place (main diagonal zeroed) — pass a copy if the input must
    be preserved.

    Args:
        mat (np.ndarray): input contact matrix, shape ``(H, W)``.
        alpha (float): RWR restart parameter.
        kernel_size (int): size of the Gaussian kernel.
        sigma (float): standard deviation of the Gaussian kernel.
        if_convolve (bool): whether to apply Gaussian convolution.
        ignore_diag (bool): whether to zero the main diagonal (before and
            after imputation).
        oldway (bool): legacy diagonal +1 mode (differs from the production
            mode when convolution has placed non-zero values on the
            sub/super-diagonal; do not compare against the reference data).
        norm (bool): whether to full-sum-normalize the result.
        upper_band (int): if not None, keep only the main diagonal and the
            upper ``upper_band`` super-diagonals, zero the rest.
        convolve_dtype (np.dtype): dtype of the convolution output (can be
            int for research purposes, e.g. ``np.int32`` in production).
        kernel_dtype (np.dtype): dtype of the Gaussian kernel (should be
            float for proper convolution).
        rwr_dtype (np.dtype): dtype for the RWR solve (float64 default).
        if_rwr (bool): whether to run the RWR step.

    Returns:
        np.ndarray: imputed matrix, shape ``(H, W)``.
    """
    gauss_kernel_1d = windows.gaussian(kernel_size, std=sigma)
    gauss_kernel_2d = np.outer(gauss_kernel_1d, gauss_kernel_1d)

    # ignore diagonal elements
    if ignore_diag:
        np.fill_diagonal(mat, 0)

    if if_convolve:
        # upcast to kernel_dtype (float) for proper convolution computation,
        # then cast the result to convolve_dtype afterwards (can be int)
        mat = mat.astype(kernel_dtype)
        mat = ndimage.convolve(mat, gauss_kernel_2d.astype(kernel_dtype), mode="constant", cval=0.0)
        mat = mat.astype(convolve_dtype)

    # add 1 to the first sub-diagonal and super-diagonal (minimal prior so
    # that even zero-count rows/columns get a stochastic column)
    if oldway:
        np.fill_diagonal(mat[1:, :-1], mat[1:, :-1] + 1)
        np.fill_diagonal(mat[:-1, 1:], mat[:-1, 1:] + 1)
    else:
        np.fill_diagonal(mat[1:, :-1], np.diag(mat[1:, :-1]) + 1)
        np.fill_diagonal(mat[:-1, 1:], np.diag(mat[:-1, 1:]) + 1)

    # optionally run RWR
    if if_rwr:
        # column-normalize to a stochastic matrix; zero columns stay as
        # NaN and propagate identically on the GPU path (nansum semantics)
        mat = mat / np.nansum(mat, axis=0)
        mat = mat.astype(rwr_dtype)
        mat = solve_rwr_inverse(mat, alpha)
    if ignore_diag:
        np.fill_diagonal(mat, 0)

    # normalize the result
    if norm:
        mat = mat / np.nansum(mat)

    # apply the upper band constraint if specified
    if upper_band is not None:
        n = mat.shape[0]
        # zero the lower triangle and everything beyond the upper band
        mask1 = np.triu(np.ones((n, n), dtype=bool), k=upper_band + 1)
        mask2 = np.tril(np.ones((n, n), dtype=bool), k=-1)
        mask = mask1 | mask2
        mat[mask] = 0

    return mat


def schicluster_imputation_gpu(mat, alpha=0.05, kernel_size=3, sigma=2, if_convolve=True,
                               ignore_diag=True, norm=True, upper_band=None, verbosity=0,
                               device=None, convolve_dtype=torch.float32,
                               kernel_dtype=torch.float32, rwr_dtype=torch.float32,
                               if_rwr=True):
    """Schicluster imputation of a batch of contact matrices, GPU version.

    Matches the CPU reference (:func:`schicluster_imputation_cpu`) step for
    step, but operates on a batch with ``torch.nn.functional.conv2d`` and a
    batched ``torch.linalg.solve``.

    Args:
        mat: batch of input matrices, shape ``(batch_size, 1, H, W)``
            (ndarray or tensor).
        alpha (float): RWR restart parameter.
        kernel_size (int): size of the Gaussian kernel.
        sigma (float): standard deviation of the Gaussian kernel.
        if_convolve (bool): whether to apply Gaussian convolution.
        ignore_diag (bool): whether to zero the main diagonal (before and
            after imputation).
        norm (bool): whether to full-sum-normalize the result.
        upper_band (int): if not None, keep only the main diagonal and the
            upper ``upper_band`` super-diagonals, zero the rest.
        verbosity (int): 0 = silent, >=1 = per-stage timing logs.
        device (torch.device): device to run on; defaults to CUDA if
            available, else CPU.
        convolve_dtype (torch.dtype): dtype of the convolution output
            (``torch.int32`` in production).
        kernel_dtype (torch.dtype): dtype of the Gaussian kernel (float).
        rwr_dtype (torch.dtype): dtype for the RWR solve (float32 in
            production).
        if_rwr (bool): whether to run the RWR step.

    Returns:
        torch.Tensor: imputed matrices on CPU, shape ``(B, 1, H, W)``.
    """
    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if verbosity >= 1:
        start_time = time.time()
        print(f"Transferring data to device {device}...")
    # preserve the original dtype on transfer; do not force convolve_dtype upfront
    if isinstance(mat, torch.Tensor):
        mat_tensor = mat.to(device=device)
    else:
        mat_tensor = torch.tensor(mat, device=device)
    if verbosity >= 1:
        transfer_time = time.time() - start_time
        print(f"Data transferred in {transfer_time:.2f} seconds.")
        print("Start pre-processing...")
        start_time = time.time()
    # ignore diagonal elements
    if ignore_diag:
        n = mat_tensor.shape[2]
        diag_indices = torch.arange(n, device=device)
        mat_tensor[:, :, diag_indices, diag_indices] = 0
    if if_convolve:
        gauss_kernel_1d = windows.gaussian(kernel_size, std=sigma)
        gauss_kernel_2d = np.outer(gauss_kernel_1d, gauss_kernel_1d)
        # upcast to kernel_dtype (float) for proper convolution computation
        mat_tensor = mat_tensor.to(kernel_dtype)
        gauss_kernel = torch.tensor(gauss_kernel_2d, device=device, dtype=kernel_dtype)
        gauss_kernel = gauss_kernel.unsqueeze(0).unsqueeze(0)  # (1, 1, kH, kW)
        # convolution with zero padding at the borders
        mat_tensor = torch.nn.functional.conv2d(mat_tensor, gauss_kernel, padding=kernel_size // 2)
        # cast to convolve_dtype after convolution (can be int)
        mat_tensor = mat_tensor.to(convolve_dtype)

    n = mat_tensor.shape[2]
    sub_diag_indices = torch.arange(n - 1, device=device)
    # add 1 to the first sub-diagonal and super-diagonal
    mat_tensor[:, :, sub_diag_indices + 1, sub_diag_indices] += 1
    mat_tensor[:, :, sub_diag_indices, sub_diag_indices + 1] += 1

    # column-normalize to a stochastic matrix, matching the CPU pre-RWR step
    mat_tensor = mat_tensor / torch.nansum(mat_tensor, dim=2, keepdim=True)
    if if_rwr:
        mat_tensor = mat_tensor.to(rwr_dtype)
    if verbosity >= 1:
        print(f"Pre-processing done in {time.time() - start_time:.2f} seconds.")
        print(f"Start RWR on device {device}...")
        start_time = time.time()
    # optionally run RWR
    if if_rwr:
        result = solve_rwr_inverse_gpu(mat_tensor, alpha, rwr_dtype=rwr_dtype, device=device)
    else:
        result = mat_tensor
    if verbosity >= 1:
        print(f"RWR done in {time.time() - start_time:.2f} seconds on device {device}.")
        print("Start post-processing...")
        start_time = time.time()
    # fill the diagonal with 0 again after imputation
    if ignore_diag:
        result[:, :, diag_indices, diag_indices] = 0

    if norm:
        # full-sum normalization across the spatial dims (H, W)
        result = result / torch.nansum(result, dim=(2, 3), keepdim=True)

    # apply the upper band constraint using efficient triangular operations
    if upper_band is not None:
        mask1 = torch.triu(torch.ones((n, n), device=device, dtype=torch.bool), diagonal=upper_band + 1)
        mask2 = torch.tril(torch.ones((n, n), device=device, dtype=torch.bool), diagonal=-1)
        mask = mask1 | mask2
        # expand the mask to match the batch dimension
        expanded_mask = mask.unsqueeze(0).unsqueeze(0).expand(result.shape[0], result.shape[1], -1, -1)
        result = result * (~expanded_mask).to(result.dtype)
    if verbosity >= 1:
        print(f"Post-processing done in {time.time() - start_time:.2f} seconds.")
        print("Start transferring data back to CPU...")
        start_time = time.time()
    result = result.cpu()
    if verbosity >= 1:
        print(f"Data transferred back to CPU in {time.time() - start_time:.2f} seconds.")
    return result


# --- cool loading helpers --- #


def load_single_hic_matrix(coolp, chrom="chr1", resolution=10000, band_width=4_000_000, verbosity=0):
    """Load a single Hi-C matrix from a cool file into the pipeline format.

    Args:
        coolp (str): path to the cool file.
        chrom (str): chromosome (or region) to extract.
        resolution (int): resolution of the Hi-C data (bp).
        band_width (int): maximum cis interaction distance (bp); ``None``
            disables the upper-band constraint.
        verbosity (int): 0 = silent, >=1 = log loading time.

    Returns:
        tuple: ``(mat_X, params)`` where ``mat_X`` has shape
        ``(1, 1, H, W)`` and ``params`` carries ``chrom``, ``resolution``,
        ``band_width``, ``upper_band`` and ``coolp``.
    """
    clr = cooler.Cooler(coolp)
    if verbosity >= 1:
        print(f"Loading Hi-C matrix from {coolp} for {chrom} at resolution {resolution}...")
        start_time = time.time()
    mat_X = clr.matrix(balance=False).fetch(chrom)
    mat_X = np.expand_dims(np.expand_dims(mat_X, axis=0), axis=0)  # (1, 1, H, W)

    if band_width is not None:
        upper_band = band_width // resolution
    else:
        upper_band = None
    if verbosity >= 1:
        print(f"Matrix loaded in {time.time() - start_time:.2f} seconds.")

    params = {
        "chrom": chrom,
        "resolution": resolution,
        "band_width": band_width,
        "upper_band": upper_band,
        "coolp": coolp,
    }
    return mat_X, params


# --- multiprocess pipeline --- #


def matrix_loading_worker(input_queue, output_queue, verbosity=0):
    """Loader process: fetch ``(coolp, chrom, resolution, band_width, outfile)``
    jobs and push ``(mat_X, params, outfile)`` tuples to the processing queue.

    A ``None`` item on the input queue is a poison pill that stops the
    worker.
    """
    while True:
        try:
            item = input_queue.get(timeout=1)
            if item is None:  # poison pill to stop worker
                break
            coolp, chrom, resolution, band_width, outfile = item
            mat_X, params = load_single_hic_matrix(
                coolp, chrom=chrom, resolution=resolution, band_width=band_width, verbosity=verbosity
            )
            if verbosity >= 1:
                print(f"[Loader] Matrix loaded in for {coolp}. Putting to processing queue...")
            output_queue.put((mat_X, params, outfile))  # may block if the queue is full
            input_queue.task_done()
            if verbosity >= 1:
                print(f"[Loader] Queued for processing: {outfile}")
        except queue.Empty:
            continue
        except Exception as e:
            print(f"Error in matrix loading worker: {e}")
            import traceback
            traceback.print_exc()
            input_queue.task_done()


def gpu_processing_worker(input_queue, output_queue, processing_params, verbosity=0):
    """GPU worker thread: impute matrices sequentially and hand them to the
    saver queue.

    Runs as a single thread in the parent process because CUDA cannot be
    used in forked subprocesses; GPU work is inherently sequential (one
    device, one stream of batched solves), so one worker suffices.

    ``processing_params`` keys: ``method`` ('torch_ver' or 'numpy_ver'),
    ``device`` ('gpu' or 'cpu'), ``if_convolve``, ``if_rwr``,
    ``ignore_diag``, ``norm``, ``oldway``, and the three dtypes
    (``convolve_dtype`` / ``kernel_dtype`` / ``rwr_dtype``).
    """
    method = processing_params.get("method", "torch_ver")
    device_str = processing_params.get("device", "gpu")
    device_obj = torch.device("cuda" if device_str == "gpu" and torch.cuda.is_available() else "cpu")

    if_convolve = processing_params.get("if_convolve", True)
    if_rwr = processing_params.get("if_rwr", True)
    ignore_diag = processing_params.get("ignore_diag", True)
    norm = processing_params.get("norm", True)
    oldway = processing_params.get("oldway", False)

    if method == "torch_ver":
        convolve_dtype = processing_params.get("convolve_dtype", torch.float32)
        kernel_dtype = processing_params.get("kernel_dtype", torch.float32)
        rwr_dtype = processing_params.get("rwr_dtype", torch.float32)
    else:
        convolve_dtype = processing_params.get("convolve_dtype", np.float64)
        kernel_dtype = processing_params.get("kernel_dtype", np.float64)
        rwr_dtype = processing_params.get("rwr_dtype", np.float64)

    processed_count = 0
    while True:
        try:
            item = input_queue.get(timeout=0.1)
            if item is None:  # poison pill to stop worker
                if verbosity >= 1:
                    print(f"GPU worker processed {processed_count} matrices. Shutting down.")
                break
            mat_X, params, outfile = item
            processed_count += 1
            if verbosity >= 1:
                print(f"[GPU Worker] Processing {processed_count}: {outfile} "
                      f"using {method} on {device_str}...")
                start_time = time.time()

            if method == "numpy_ver":
                # NumPy CPU implementation — squeeze to 2D first
                mat_2d = mat_X[0, 0, :, :]  # (1, 1, H, W) -> (H, W)
                imputed_mat_X = schicluster_imputation_cpu(
                    mat_2d, alpha=0.05, kernel_size=3, sigma=2,
                    if_convolve=if_convolve, ignore_diag=ignore_diag,
                    oldway=oldway, norm=norm, upper_band=params["upper_band"], if_rwr=if_rwr,
                    convolve_dtype=convolve_dtype, kernel_dtype=kernel_dtype, rwr_dtype=rwr_dtype,
                )
            else:
                imputed_mat_X = schicluster_imputation_gpu(
                    mat_X, alpha=0.05, kernel_size=3, sigma=2,
                    if_convolve=if_convolve, ignore_diag=ignore_diag,
                    norm=norm, upper_band=params["upper_band"], if_rwr=if_rwr,
                    verbosity=verbosity, device=device_obj,
                    convolve_dtype=convolve_dtype, kernel_dtype=kernel_dtype, rwr_dtype=rwr_dtype,
                )
                imputed_mat_X = imputed_mat_X[0, 0, :, :].numpy()  # squeeze back to (H, W)

            if verbosity >= 1:
                print(f"[GPU Worker] Completed {processed_count}: {outfile} in "
                      f"{time.time() - start_time:.2f}s. Putting to save queue...")
            output_queue.put((imputed_mat_X, params, outfile))
            input_queue.task_done()
            if verbosity >= 1:
                print(f"[GPU Worker] Queued for saving. Ready for next matrix.")
        except queue.Empty:
            if verbosity >= 2:
                print(f"[GPU Worker] Queue empty, waiting... (processed {processed_count} so far)")
            continue
        except Exception as e:
            print(f"Error in processing worker: {e}")
            import traceback
            traceback.print_exc()
            input_queue.task_done()


def saving_worker(input_queue, verbosity=0):
    """Saver process: convert imputed matrices to CSR and write compressed
    ``.npz`` files.

    A ``None`` item on the input queue is a poison pill that stops the
    worker.
    """
    while True:
        try:
            item = input_queue.get(timeout=0.1)
            if item is None:  # poison pill to stop worker
                break
            imputed_mat_X, params, outfile = item
            if verbosity >= 1:
                print(f"[Saver] Saving {outfile}...")
            save_single_result(imputed_mat_X, params, outfile, verbosity)
            input_queue.task_done()
            if verbosity >= 1:
                print(f"[Saver] Completed {outfile}")
        except queue.Empty:
            continue
        except Exception as e:
            print(f"Error in saving worker: {e}")
            import traceback
            traceback.print_exc()
            input_queue.task_done()


def save_single_result(imputed_mat_X, params, outfile, verbosity=0):
    """Save one imputed matrix as a compressed CSR ``.npz``.

    The stored arrays are ``data`` / ``indices`` / ``indptr`` / ``shape``;
    load with ``csr_matrix((z["data"], z["indices"], z["indptr"]),
    shape=z["shape"])``.  The stored matrix is the imputed, normalized,
    band-limited upper triangle (plus main diagonal).

    Args:
        imputed_mat_X (np.ndarray): imputed matrix, shape ``(H, W)``.
        params (dict): pipeline params (informational; not stored).
        outfile (str): output ``.npz`` path.
        verbosity (int): 0 = silent, >=1 = per-stage timing logs.
    """
    total_start = time.time()
    if verbosity >= 1:
        print(f"[Save Details] Csr transformation for {outfile} ...")
        start_time = time.time()
    # convert to sparse matrix before saving
    sparse_matrix = csr_matrix(imputed_mat_X)
    if verbosity >= 1:
        print(f"[Save Details] Csr transformation done in {time.time() - start_time:.2f}s.")
        print(f"[Save Details] Start npz compression for {outfile} ...")
        start_time = time.time()
    np.savez_compressed(
        outfile, data=sparse_matrix.data, indices=sparse_matrix.indices,
        indptr=sparse_matrix.indptr, shape=sparse_matrix.shape,
    )
    if verbosity >= 1:
        print(f"[Save Details] Npz saved in {time.time() - start_time:.2f}s. "
              f"Total save time: {time.time() - total_start:.2f}s for {outfile}.")


def hic_imputation_multiprocess(coolps, outfiles, chroms=None, resolution=10000,
                                band_width=4_000_000, ignore_diag=True, norm=True,
                                verbosity=0, nproc_loading=1, nproc_saving=5,
                                method="torch_ver", device="gpu", oldway=False,
                                if_convolve=True, if_rwr=True, convolve_dtype=None,
                                kernel_dtype=None, rwr_dtype=None):
    """Multi-process Hi-C imputation using a producer-consumer pipeline.

    Layout: ``nproc_loading`` loader processes -> 1 GPU worker thread ->
    ``nproc_saving`` saver processes.  Queues are ``JoinableQueue`` with
    capacity ``3 x n_workers`` to buffer while a slower stage catches up;
    ``None`` poison pills shut each stage down in order after all work is
    drained.

    Args:
        coolps (list): paths to cool files (one per task).
        outfiles (list): output ``.npz`` paths, one per task.
        chroms (list or str): chromosome per task, or a single chromosome
            applied to all tasks.
        resolution (int): resolution of the Hi-C data (bp).
        band_width (int): maximum cis interaction distance (bp).
        ignore_diag (bool): whether to zero the main diagonal.
        norm (bool): whether to full-sum-normalize.
        verbosity (int): 0 = silent, >=1 = pipeline progress logs.
        nproc_loading (int): number of loader processes.
        nproc_saving (int): number of saver processes.
        method (str): ``'torch_ver'`` (GPU path) or ``'numpy_ver'`` (CPU,
            in-thread).
        device (str): ``'gpu'`` or ``'cpu'`` (only for ``'torch_ver'``).
        oldway (bool): legacy diagonal +1 mode for the CPU path.
        if_convolve (bool): whether to apply Gaussian convolution.
        if_rwr (bool): whether to run the RWR step.
        convolve_dtype: convolution output dtype; torch dtype for
            ``'torch_ver'`` (float32 default), numpy dtype for
            ``'numpy_ver'`` (float64 default).
        kernel_dtype: Gaussian kernel dtype (float).
        rwr_dtype: RWR solve dtype.
    """
    # validate inputs
    if len(coolps) != len(outfiles):
        raise ValueError("Length of coolps and outfiles must be equal")

    # handle chroms — convert to a list if a single string
    if chroms is None:
        chroms = ["chr1"] * len(coolps)
    elif isinstance(chroms, str):
        chroms = [chroms] * len(coolps)
    elif len(chroms) != len(coolps):
        raise ValueError("Length of chroms must match length of coolps and outfiles")

    # default dtypes if not provided
    if method == "torch_ver":
        if convolve_dtype is None:
            convolve_dtype = torch.float32
        if kernel_dtype is None:
            kernel_dtype = torch.float32
        if rwr_dtype is None:
            rwr_dtype = torch.float32
    else:
        if convolve_dtype is None:
            convolve_dtype = np.float64
        if kernel_dtype is None:
            kernel_dtype = np.float64
        if rwr_dtype is None:
            rwr_dtype = np.float64

    processing_params = {
        "method": method,
        "device": device,
        "ignore_diag": ignore_diag,
        "norm": norm,
        "oldway": oldway,
        "if_convolve": if_convolve,
        "if_rwr": if_rwr,
        "convolve_dtype": convolve_dtype,
        "kernel_dtype": kernel_dtype,
        "rwr_dtype": rwr_dtype,
    }

    # inter-process queues; 3x the number of workers so a slower stage can
    # buffer while the GPU stays fed
    task_queue = JoinableQueue()
    loading_to_processing_queue = JoinableQueue(maxsize=max(3, nproc_loading * 3))
    processing_to_saving_queue = JoinableQueue(maxsize=max(3, nproc_saving * 3))

    if verbosity >= 1:
        print(f"Queue sizes: loading->processing={max(3, nproc_loading * 3)}, "
              f"processing->saving={max(3, nproc_saving * 3)}")
        print(f"Workers: {nproc_loading} loaders, 1 GPU processor, {nproc_saving} savers")

    # start loader processes
    loading_workers = []
    for _ in range(nproc_loading):
        p = Process(target=matrix_loading_worker,
                    args=(task_queue, loading_to_processing_queue, verbosity))
        p.start()
        loading_workers.append(p)

    # start the processing worker as a THREAD (not process) because CUDA
    # cannot be used in forked subprocesses; GPU work is sequential, so one
    # thread suffices
    processing_thread = threading.Thread(
        target=gpu_processing_worker,
        args=(loading_to_processing_queue, processing_to_saving_queue, processing_params, verbosity),
    )
    processing_thread.start()

    # start saver processes
    saving_workers = []
    for _ in range(nproc_saving):
        p = Process(target=saving_worker, args=(processing_to_saving_queue, verbosity))
        p.start()
        saving_workers.append(p)

    # add tasks to the loading queue
    if verbosity >= 1:
        print(f"Adding {len(coolps)} tasks to task queue...")
    for coolp, chrom, outfile in zip(coolps, chroms, outfiles):
        task_queue.put((coolp, chrom, resolution, band_width, outfile))
    if verbosity >= 1:
        print("All tasks added. Starting pipeline...")

    # wait for all loading tasks to be consumed
    task_queue.join()
    # send poison pills to loader processes
    for _ in range(nproc_loading):
        task_queue.put(None)
    for p in loading_workers:
        p.join()
    if verbosity >= 1:
        print("All loaders finished. Waiting for processing queue...")

    # wait for all loaded matrices to be processed
    loading_to_processing_queue.join()
    # send the poison pill to the GPU worker thread
    loading_to_processing_queue.put(None)
    processing_thread.join()
    if verbosity >= 1:
        print("GPU worker finished. Waiting for saving queue...")

    # wait for all processed matrices to be saved
    processing_to_saving_queue.join()
    # send poison pills to saver processes
    for _ in range(nproc_saving):
        processing_to_saving_queue.put(None)
    for p in saving_workers:
        p.join()

    if verbosity >= 1:
        print("All workers finished.")
    print("All tasks completed.")
