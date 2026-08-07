#!/usr/bin/env python
"""Visual test for schicluster_gpu imputation: pre/post-imputation heatmaps.

Generates PNG heatmaps of the raw contact matrix (pre-imputation) and the
imputed matrix (post-imputation) for both the synthetic and the real
reference datasets, so the effect of the imputation can be inspected by eye.

Usage:
    python tests/test_visual_schicluster_gpu.py      # runs both datasets

    # or collect a single dataset with pytest:
    pytest tests/test_visual_schicluster_gpu.py

Output:
    tests/output/visual/schicluster_gpu/
        synthetic_pre_imputation.png
        synthetic_post_imputation.png
        real_pre_imputation.png
        real_post_imputation.png
"""
import sys
from pathlib import Path

sys.path.insert(0, "/share/home/ychi/dev/hic_basic")

import matplotlib

matplotlib.use("Agg")  # headless

import numpy as np
import torch
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import pyplot as plt

from hic_basic.impute.schicluster_gpu import schicluster_imputation_gpu

DATA_DIR = Path(__file__).parent / "data" / "schicluster_gpu"
OUTPUT_DIR = Path(__file__).parent / "output" / "visual" / "schicluster_gpu"

# production parameters (issues/schicluster_gpu/alg.md §4)
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

# single-hue sequential ramp (blue, light -> dark) for magnitude data,
# from the validated reference palette
CMAP = LinearSegmentedColormap.from_list(
    "hic_blue",
    ["#cde2fb", "#9ec5f4", "#5598e7", "#2a78d6", "#1c5cab", "#0d366b"],
)


def plot_matrix(mat, out_path, title, cbar_label, vmin, vmax, transform="log1p"):
    """Save one heatmap of a contact matrix for visual inspection.

    Raw contact counts span orders of magnitude and are plotted on a log1p
    scale.  The imputed matrix is full-sum-normalized and heavy-tailed
    (most entries are solver noise ~1e-18..1e-7, real signal up to ~1e-4),
    so it is plotted as log10 with a 1e-7 noise floor — everything at or
    below the floor renders as the lightest ramp step and only the real
    near-diagonal signal carries color.
    """
    if transform == "log1p":
        data = np.log1p(mat)
    elif transform == "log10":
        data = np.log10(mat + 1e-7)
    else:
        data = mat
    fig, ax = plt.subplots(figsize=(6.2, 5.4))
    im = ax.imshow(data, cmap=CMAP, vmin=vmin, vmax=vmax, origin="lower")
    ax.set_title(title)
    ax.set_xlabel("bin index")
    ax.set_ylabel("bin index")
    ax.set_aspect("equal")
    cbar = fig.colorbar(im, ax=ax, shrink=0.9, pad=0.02)
    cbar.set_label(cbar_label)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {out_path}")
    return out_path


def impute(mat):
    """Production GPU path (auto cuda if available, else cpu); returns (H, W)."""
    return schicluster_imputation_gpu(mat[None, None], **PROD)[0, 0].numpy()


def test_visual_synthetic_pre_and_post_imputation():
    """Heatmaps of the synthetic 100x100 matrix before and after imputation."""
    mat = np.load(DATA_DIR / "synthetic_input.npz")["mat"]
    imputed = impute(mat)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    # raw counts: log1p, clip to the 99th percentile of non-zero entries
    vmax = np.nanpercentile(np.log1p(mat[mat > 0]), 99)
    plot_matrix(mat, OUTPUT_DIR / "synthetic_pre_imputation.png",
                title="synthetic — raw matrix (pre-imputation)",
                cbar_label="log1p(contacts)", vmin=0, vmax=vmax)
    # imputed: log10 with a 1e-7 noise floor; the ramp spans the data
    # (5th..99.9th percentile) so even small variation stays visible
    scaled = np.log10(imputed + 1e-7)
    vmin, vmax = np.nanpercentile(scaled, [5, 99.9])
    plot_matrix(imputed, OUTPUT_DIR / "synthetic_post_imputation.png",
                title="synthetic — imputed (post-imputation)",
                cbar_label="log10(imputed + 1e-7)", vmin=vmin, vmax=vmax,
                transform="log10")


def test_visual_real_pre_and_post_imputation():
    """Heatmaps of the real single-cell Hi-C region before and after
    imputation (chr21:0-10 Mb @ 10 kb, 1000 bins)."""
    import cooler

    clr = cooler.Cooler(str(DATA_DIR / "real_region.cool"))
    mat = np.asarray(clr.matrix(balance=False).fetch("chr21"))
    imputed = impute(mat)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    vmax = np.nanpercentile(np.log1p(mat[mat > 0]), 99)
    plot_matrix(mat, OUTPUT_DIR / "real_pre_imputation.png",
                title="real (chr21:0-10Mb @10kb) — raw matrix (pre-imputation)",
                cbar_label="log1p(contacts)", vmin=0, vmax=vmax)
    scaled = np.log10(imputed + 1e-7)
    vmin, vmax = np.nanpercentile(scaled, [5, 99.9])
    plot_matrix(imputed, OUTPUT_DIR / "real_post_imputation.png",
                title="real (chr21:0-10Mb @10kb) — imputed (post-imputation)",
                cbar_label="log10(imputed + 1e-7)", vmin=vmin, vmax=vmax,
                transform="log10")


if __name__ == "__main__":
    print("=" * 60)
    print("Visual Test: schicluster_gpu imputation pre/post heatmaps")
    print("=" * 60)
    if not (DATA_DIR / "synthetic_input.npz").exists():
        print(f"Error: test data not found in {DATA_DIR}")
        print("Run: ./run_in_env.sh --env hic_basic_v096 python "
              "tests/data/create_schicluster_gpu_data.py")
        sys.exit(1)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    test_visual_synthetic_pre_and_post_imputation()
    test_visual_real_pre_and_post_imputation()
    print("\n" + "=" * 60)
    print(f"All visual tests completed!")
    print(f"Output directory: {OUTPUT_DIR}")
    print("=" * 60)
