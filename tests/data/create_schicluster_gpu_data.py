"""Regenerate the schicluster_gpu test dataset: inputs + reference outputs.

The reference outputs are produced by the *production* GPU imputation
implementation from the hemsc repository
(`/share/home/ychi/dev/hemsc/notebooks/5_imputation.py`), which is the code
being integrated into hic_basic. This script exists so the reference data
can be independently reproduced; it does not depend on any hic_basic
imputation code.

Outputs (written into tests/data/schicluster_gpu/):

- synthetic_input.npz          deterministic 100x100 simulated contact matrix
                               (int64 dense, key "mat")
- synthetic_expected.npz       imputed result, production parameters
- real_region.cool             real single-cell Hi-C region (chr21 first 10 Mb
                               at 10 kb) extracted from a real human mcool
- real_expected.npz            imputed result, production parameters

Reference outputs are written only when a CUDA device is available (they are
the GPU float32 production outputs). Without a GPU the script runs the CPU
reference and reports GPU/CPU agreement.

Usage (with a GPU):

    ./run_in_env.sh --gpu 1 --env hic_basic_v096 python \\
        tests/data/create_schicluster_gpu_data.py

To regenerate the real region from a different source mcool:

    SCHICLUSTER_GPU_SOURCE_MCOOL=/path/to/sample.mcool::/resolutions/10000 \\
        SCHICLUSTER_GPU_REGION=chr21:0-10000000 python create_schicluster_gpu_data.py
"""
import importlib
import os
import sys
from pathlib import Path

import cooler
import numpy as np
import torch

# production implementation from the hemsc repository (source of the port);
# the module name starts with a digit, so importlib is required
sys.path.insert(0, "/share/home/ychi/dev/hemsc")
_hemsc = importlib.import_module("notebooks.5_imputation")
schicluster_imputation_cpu = _hemsc.schicluster_imputation_cpu
schicluster_imputation_gpu = _hemsc.schicluster_imputation_gpu

DATA_DIR = Path(__file__).resolve().parent / "schicluster_gpu"

# Production parameter set (see issues/schicluster_gpu/alg.md)
PROD = dict(
    alpha=0.05,
    kernel_size=3,
    sigma=2,
    if_convolve=True,
    ignore_diag=True,
    norm=True,
    band_width=4_000_000,
    convolve_dtype=torch.int32,
    kernel_dtype=torch.float32,
    rwr_dtype=torch.float32,
)


def generate_test_matrix(size=100, seed=42):
    """Simulated Hi-C contact matrix: diagonal 20-30, off-diagonal 0-5 with
    70% sparsity, symmetric (deterministic given seed)."""
    rng = np.random.default_rng(seed)
    mat = np.zeros((size, size), dtype=int)
    diag_values = rng.integers(20, 31, size)
    np.fill_diagonal(mat, diag_values)
    mask = rng.random((size, size)) > 0.7
    upper_rand = rng.integers(1, 6, size=(size, size))
    mat_upper = upper_rand * mask
    mat += np.triu(mat_upper, k=1)
    mat += np.triu(mat_upper, k=1).T
    return mat


def make_real_region_cool(source_mcool, region, cool_path):
    """Extract a real region from a (multi)cool into a small standalone cool."""
    clr = cooler.Cooler(source_mcool)
    bins = clr.bins().fetch(region).reset_index()  # keep global bin id as "index"
    old2new = {old: new for new, old in enumerate(bins["index"])}
    pixels = clr.pixels().fetch(region)  # columns: bin1_id, bin2_id, count
    # source single-cell mcools may contain duplicate contacts; aggregate so
    # the extracted matrix matches clr.matrix(balance=False).fetch(region)
    pixels = pixels.groupby(["bin1_id", "bin2_id"], as_index=False)["count"].sum()
    pixels = pixels.assign(
        bin1_id=pixels["bin1_id"].map(old2new),
        bin2_id=pixels["bin2_id"].map(old2new),
    )
    # pixels().fetch() returns rows whose bin1 is in the region even when bin2
    # lies outside it; the matrix fetch drops those, so drop them here too
    pixels = pixels.dropna().sort_values(["bin1_id", "bin2_id"]).reset_index(drop=True)
    pixels = pixels.astype({"bin1_id": int, "bin2_id": int})
    new_bins = bins[["chrom", "start", "end"]].reset_index(drop=True)
    metadata = clr.info.get("metadata")
    cooler.create_cooler(
        str(cool_path), new_bins, pixels, ordered=True,
        assembly=metadata.get("assembly") if metadata else None,
    )
    return cool_path


def load_region_matrix(cool_path, region):
    clr = cooler.Cooler(str(cool_path))
    mat = np.asarray(clr.matrix(balance=False).fetch(region))
    return np.expand_dims(np.expand_dims(mat, axis=0), axis=0)


def save_expected(mat_2d, out_path):
    from scipy.sparse import csr_matrix

    sparse = csr_matrix(mat_2d)
    np.savez_compressed(
        out_path,
        data=sparse.data,
        indices=sparse.indices,
        indptr=sparse.indptr,
        shape=sparse.shape,
    )
    print(f"wrote {out_path}  (shape={sparse.shape}, nnz={sparse.nnz}, sum={mat_2d.sum():.4f})")


def run_and_save(mat_X, name, device, out_dir):
    upper_band = PROD["band_width"] // 10000
    imputed_gpu = schicluster_imputation_gpu(
        mat_X, alpha=PROD["alpha"], kernel_size=PROD["kernel_size"], sigma=PROD["sigma"],
        if_convolve=PROD["if_convolve"], ignore_diag=PROD["ignore_diag"],
        norm=PROD["norm"], upper_band=upper_band,
        convolve_dtype=PROD["convolve_dtype"], kernel_dtype=PROD["kernel_dtype"],
        rwr_dtype=PROD["rwr_dtype"], device=device,
    )[0, 0].numpy()

    # CPU reference: production convolution output dtype (int32) on both paths,
    # only the RWR precision differs (float32 on GPU vs float64 here)
    imputed_cpu = schicluster_imputation_cpu(
        mat_X[0, 0].copy(), alpha=PROD["alpha"], kernel_size=PROD["kernel_size"], sigma=PROD["sigma"],
        if_convolve=PROD["if_convolve"], ignore_diag=PROD["ignore_diag"],
        norm=PROD["norm"], upper_band=upper_band,
        convolve_dtype=np.int32, kernel_dtype=np.float64, rwr_dtype=np.float64,
    )

    max_diff = float(np.abs(imputed_gpu - imputed_cpu).max())
    rel = float(np.abs(imputed_gpu - imputed_cpu).sum() / max(1e-30, np.abs(imputed_gpu).sum()))
    print(f"[{name}] GPU(device={device}) vs CPU max abs diff = {max_diff:.3e}, rel L1 = {rel:.3e}")

    if device.type == "cuda":
        save_expected(imputed_gpu, out_dir / f"{name}_expected.npz")
    else:
        print(f"[{name}] no CUDA device: reference output not written "
              f"(CPU only differs by ~{max_diff:.1e} from GPU)")


def main():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"device: {device}")

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    out_dir = DATA_DIR

    # --- synthetic input ---
    mat = generate_test_matrix(size=100, seed=42)
    np.savez_compressed(out_dir / "synthetic_input.npz", mat=mat)
    print(f"wrote {out_dir / 'synthetic_input.npz'}")
    mat_X = np.expand_dims(np.expand_dims(mat, axis=0), axis=0).astype(np.float64)
    run_and_save(mat_X, "synthetic", device, out_dir)

    # --- real region ---
    source_mcool = os.environ.get(
        "SCHICLUSTER_GPU_SOURCE_MCOOL",
        "/shareb/mliu/human/processed_data/mcool/A1a5.mcool::/resolutions/10000",
    )
    region = os.environ.get("SCHICLUSTER_GPU_REGION", "chr21:0-10000000")
    cool_path = out_dir / "real_region.cool"
    if not cool_path.exists():
        if not os.path.exists(source_mcool.split("::")[0]):
            print(f"source mcool not found: {source_mcool}; skipping real region "
                  f"(re-run with SCHICLUSTER_GPU_SOURCE_MCOOL set)")
        else:
            make_real_region_cool(source_mcool, region, cool_path)
            print(f"wrote {cool_path}")
    else:
        print(f"reusing existing {cool_path}")
    if cool_path.exists():
        mat_X = load_region_matrix(cool_path, "chr21")
        print(f"real region matrix shape: {mat_X.shape}")
        run_and_save(mat_X, "real", device, out_dir)


if __name__ == "__main__":
    main()
