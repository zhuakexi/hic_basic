# Acceptance tests for hic_basic.impute.schicluster_gpu against the
# GPU-generated reference data in tests/data/schicluster_gpu/.
#
# Acceptance criteria (issues/schicluster_gpu/alg.md §6): with the production
# parameters (alpha=0.05, kernel_size=3, sigma=2, if_convolve=True,
# ignore_diag=True, norm=True, upper_band=400, convolve_dtype=torch.int32,
# kernel_dtype=torch.float32, rwr_dtype=torch.float32):
#   1. GPU path matches the reference with allclose(rtol=1e-5, atol=1e-6)
#   2. CPU reference path matches with allclose(rtol=1e-5, atol=1e-7)
#   3. Structural invariants (exact): main diagonal 0, lower triangle 0,
#      nothing beyond diagonal 400, sums 0.500000 / 0.494426
#   4. End-to-end cool -> impute -> compressed CSR .npz reproduces
#      real_expected.npz
import os
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, "/share/home/ychi/dev/hic_basic")

import cooler
import numpy as np
import torch
from scipy.sparse import csr_matrix

from hic_basic.impute.schicluster_gpu import (
    hic_imputation_multiprocess,
    save_single_result,
    schicluster_imputation_cpu,
    schicluster_imputation_gpu,
    solve_rwr_inverse,
    solve_rwr_inverse_gpu,
)

DATA_DIR = Path(os.path.dirname(__file__)) / "data" / "schicluster_gpu"

# Production parameter set (see issues/schicluster_gpu/alg.md §4)
PROD = dict(
    alpha=0.05,
    kernel_size=3,
    sigma=2,
    if_convolve=True,
    ignore_diag=True,
    norm=True,
    upper_band=400,
    convolve_dtype=torch.int32,
    kernel_dtype=torch.float32,
    rwr_dtype=torch.float32,
)
# CPU reference path: same convolution output dtype (int32), float64 RWR
PROD_CPU = dict(
    alpha=0.05,
    kernel_size=3,
    sigma=2,
    if_convolve=True,
    ignore_diag=True,
    norm=True,
    upper_band=400,
    convolve_dtype=np.int32,
    kernel_dtype=np.float64,
    rwr_dtype=np.float64,
)


def load_csr_npz(path):
    """Load a compressed CSR .npz (data/indices/indptr/shape) as a dense array."""
    z = np.load(path)
    return csr_matrix((z["data"], z["indices"], z["indptr"]), shape=z["shape"]).toarray()


def load_synthetic_input():
    return np.load(DATA_DIR / "synthetic_input.npz")["mat"]


def load_real_matrix():
    clr = cooler.Cooler(str(DATA_DIR / "real_region.cool"))
    return np.asarray(clr.matrix(balance=False).fetch("chr21"))


class TestSchiclusterGpu(unittest.TestCase):
    """Quantitative acceptance tests vs the GPU-generated reference data."""

    # --- criterion 1: GPU path (torch; auto cuda if available, else cpu) ---

    def test_gpu_synthetic_matches_reference(self):
        mat = load_synthetic_input()
        expected = load_csr_npz(DATA_DIR / "synthetic_expected.npz")
        imputed = schicluster_imputation_gpu(mat[None, None], **PROD)[0, 0].numpy()
        np.testing.assert_allclose(imputed, expected, rtol=1e-5, atol=1e-6)

    def test_gpu_real_matches_reference(self):
        mat = load_real_matrix()
        expected = load_csr_npz(DATA_DIR / "real_expected.npz")
        imputed = schicluster_imputation_gpu(mat[None, None], **PROD)[0, 0].numpy()
        np.testing.assert_allclose(imputed, expected, rtol=1e-5, atol=1e-6)

    # --- criterion 2: CPU reference path (float64 RWR) ---

    def test_cpu_synthetic_matches_reference(self):
        mat = load_synthetic_input()
        expected = load_csr_npz(DATA_DIR / "synthetic_expected.npz")
        imputed = schicluster_imputation_cpu(mat.copy(), **PROD_CPU)
        np.testing.assert_allclose(imputed, expected, rtol=1e-5, atol=1e-7)

    def test_cpu_real_matches_reference(self):
        mat = load_real_matrix()
        expected = load_csr_npz(DATA_DIR / "real_expected.npz")
        imputed = schicluster_imputation_cpu(mat.copy(), **PROD_CPU)
        np.testing.assert_allclose(imputed, expected, rtol=1e-5, atol=1e-7)

    # --- criterion 3: structural invariants (exact) ---

    def test_structural_invariants(self):
        for name, mat, expected, ref_sum in [
            ("synthetic", load_synthetic_input(),
             load_csr_npz(DATA_DIR / "synthetic_expected.npz"), 0.500000),
            ("real", load_real_matrix(),
             load_csr_npz(DATA_DIR / "real_expected.npz"), 0.494426),
        ]:
            imputed = schicluster_imputation_gpu(mat[None, None], **PROD)[0, 0].numpy()
            # main diagonal is exactly 0 (ignore_diag)
            np.testing.assert_array_equal(np.diag(imputed), 0,
                                          err_msg=f"{name}: main diagonal not zero")
            # entire lower triangle is exactly 0 (upper_band mask)
            self.assertEqual(np.count_nonzero(np.tril(imputed, -1)), 0,
                             msg=f"{name}: lower triangle not zero")
            # nothing beyond diagonal upper_band
            self.assertEqual(
                np.count_nonzero(np.triu(imputed, PROD["upper_band"] + 1)), 0,
                msg=f"{name}: entries beyond diagonal {PROD['upper_band']}")
            # normalization precedes the band mask, so the kept upper
            # triangle holds about half the mass — match the reference sum
            self.assertAlmostEqual(imputed.sum(), ref_sum, places=6,
                                   msg=f"{name}: total sum mismatch")
            # and the result agrees with the reference itself
            np.testing.assert_allclose(imputed, expected, rtol=1e-5, atol=1e-6)

    # --- criterion 4: end-to-end cool -> impute -> compressed CSR .npz ---

    def test_end_to_end_save_single_result(self):
        mat = load_real_matrix()
        imputed = schicluster_imputation_gpu(mat[None, None], **PROD)[0, 0].numpy()
        params = {
            "chrom": "chr21", "resolution": 10000, "band_width": 4_000_000,
            "upper_band": 400, "coolp": str(DATA_DIR / "real_region.cool"),
        }
        with tempfile.TemporaryDirectory() as tmp:
            outfile = str(Path(tmp) / "real_out.npz")
            save_single_result(imputed, params, outfile)
            got = load_csr_npz(outfile)
        expected = load_csr_npz(DATA_DIR / "real_expected.npz")
        np.testing.assert_allclose(got, expected, rtol=1e-5, atol=1e-6)

    def test_hic_imputation_multiprocess_end_to_end(self):
        # full production pipeline: 1 loader -> 1 GPU thread -> 1 saver
        with tempfile.TemporaryDirectory() as tmp:
            outfile = str(Path(tmp) / "real_pipeline.npz")
            hic_imputation_multiprocess(
                [str(DATA_DIR / "real_region.cool")], [outfile],
                chroms="chr21", resolution=10000, band_width=4_000_000,
                nproc_loading=1, nproc_saving=1,
                method="torch_ver", device="gpu",
                convolve_dtype=torch.int32, kernel_dtype=torch.float32,
                rwr_dtype=torch.float32,
            )
            got = load_csr_npz(outfile)
        expected = load_csr_npz(DATA_DIR / "real_expected.npz")
        np.testing.assert_allclose(got, expected, rtol=1e-5, atol=1e-6)

    # --- kernel-level and option sanity checks ---

    def test_solve_rwr_inverse_gpu_matches_cpu(self):
        rng = np.random.default_rng(0)
        P = rng.random((50, 50))
        P /= P.sum(axis=0)  # column-stochastic
        s_cpu = solve_rwr_inverse(P, alpha=0.05)
        s_gpu = solve_rwr_inverse_gpu(
            P[None, None], alpha=0.05, device=torch.device("cpu"),
        )[0, 0].numpy()
        # float32 solve vs float64 solve
        np.testing.assert_allclose(s_gpu, s_cpu, rtol=1e-3, atol=1e-6)

    def test_no_band_normalizes_to_one(self):
        # without the band mask the full matrix is kept and sums to 1
        mat = load_synthetic_input()
        imputed = schicluster_imputation_cpu(
            mat.copy(), alpha=0.05, kernel_size=3, sigma=2,
            if_convolve=True, ignore_diag=True, norm=True, upper_band=None,
        )
        self.assertAlmostEqual(imputed.sum(), 1.0, places=6)
        # RWR symmetrization leaves mass below the diagonal
        self.assertGreater(np.count_nonzero(np.tril(imputed, -1)), 0)


if __name__ == "__main__":
    unittest.main()
