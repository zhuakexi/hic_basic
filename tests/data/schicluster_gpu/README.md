# schicluster_gpu test data

Inputs and reference outputs for `tests/test_schicluster_gpu.py`.

## Files

| File | Content |
|------|---------|
| `synthetic_input.npz` | Deterministic 100×100 simulated contact matrix (int64, key `mat`), diagonal 20–30, off-diagonal 0–5 at 70% sparsity, symmetric (seed 42) |
| `synthetic_expected.npz` | Imputed result of the synthetic matrix, GPU (torch.float32) reference, production parameters |
| `real_region.cool` | Real single-cell Hi-C region: chr21 first 10 Mb at 10 kb (1000 bins) extracted from a real human single-cell mcool (sample `A1a5`, `/shareb/mliu/human/processed_data/mcool/A1a5.mcool::/resolutions/10000`); the matrix is bit-identical to the source region fetch |
| `real_expected.npz` | Imputed result of the real region, GPU reference, production parameters |

Expected `.npz` files store a compressed CSR matrix as separate arrays
(`data`, `indices`, `indptr`, `shape`) — the same format produced by
`hic_basic.impute.schicluster_gpu.save_single_result`.

## Production parameters used for the references

```
alpha = 0.05, kernel_size = 3, sigma = 2
if_convolve = True, ignore_diag = True, norm = True
band_width = 4,000,000 bp at resolution = 10,000 bp  ->  upper_band = 400
convolve_dtype = torch.int32, kernel_dtype = torch.float32, rwr_dtype = torch.float32
```

## Regeneration

```bash
# GPU reference outputs (requires a CUDA device)
./run_in_env.sh --gpu 1 --env hic_basic_v096 python \
    tests/data/create_schicluster_gpu_data.py

# CPU-only run: recomputes and reports agreement but does not overwrite the
# reference outputs (no CUDA device)
./run_in_env.sh --env hic_basic_v096 python tests/data/create_schicluster_gpu_data.py

# different source mcool / region
SCHICLUSTER_GPU_SOURCE_MCOOL=/path/to/sample.mcool::/resolutions/10000 \
SCHICLUSTER_GPU_REGION=chr1:0-10000000 \
    python tests/data/create_schicluster_gpu_data.py
```

Notes:

- The references are **GPU float32** outputs — that is the production truth.
  The CPU implementation (float64 RWR, int32 convolution output) agrees with
  them to ~5e-11 (synthetic) / ~1.2e-10 (real) max abs diff.
- The reference generation script also documents the exact source mcool and
  region; delete `real_region.cool` before re-running to rebuild it from the
  source.

## Acceptance criteria (summary)

See `issues/schicluster_gpu/alg.md` §6 for the full criteria:

- GPU implementation vs reference: `allclose(rtol=1e-5, atol=1e-6)`;
  CPU float64 reference vs reference: `allclose(rtol=1e-5, atol=1e-7)`.
- Structural invariants (exact): main diagonal 0; lower triangle 0; nothing
  beyond diagonal 400.
- Expected total sums (normalization precedes the band mask, so ~half the
  mass is kept): `0.500000` (synthetic), `0.494426` (real), to ~1e-6.
