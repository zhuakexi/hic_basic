#!/usr/bin/env python
"""
Reproduce Figure 1 from sc-embryo series 6 using plot_cool_with_tracks.

This script uses real data from the sc-embryo project to debug and reproduce:
- notebooks/E_loop_analysis/7_celltype_loops/6.subset_merged_tiling_tracks/
  6.subset_merged_tiling_tracks.loop01.325348.tracks.png

See README.md for full documentation of the input data and parameters.
"""

from pathlib import Path
from hic_basic.plot.hic import plot_cool_with_tracks
from hic_basic.plot.track import BigWigTrack, GeneTrack

# Output directory
OUTPUT_DIR = Path(__file__).parent
OUTPUT_PNG = OUTPUT_DIR / "figure1_reproduced.png"

# =============================================================================
# Heatmap Inputs (tiling: lower=early, upper=late)
# =============================================================================

early_mcool_uri = (
    "/share/home/ychi/dev/sc-embryo/notebooks/E_loop_analysis/"
    "7_celltype_loops/6.subset_merged_tiling_tracks/"
    "6.subset_merged_tiling_tracks.early.matched.mcool::resolutions/20000"
)

late_mcool_uri = (
    "/share/home/ychi/dev/sc-embryo/notebooks/E_loop_analysis/"
    "7_celltype_loops/6.subset_merged_tiling_tracks/"
    "6.subset_merged_tiling_tracks.late.matched.mcool::resolutions/20000"
)

# =============================================================================
# Region and Loop
# =============================================================================

REGION = "chr10:14240000-15080000"
LOOP_ID = "325348"

# =============================================================================
# Marks (loop annotation)
# =============================================================================

marks = [
    {
        "chrom1": "chr10",
        "start1": 14640000,
        "end1": 14660000,
        "chrom2": "chr10",
        "start2": 14640000,
        "end2": 14660000,
        "color": "cyan",
        "shape": "circle",
        "width": 3,
        "label": f"loop {LOOP_ID}",
    },
    {
        "chrom1": "chr10",
        "start1": 14640000,
        "end1": 14660000,
        "chrom2": "chr10",
        "start2": 14640000,
        "end2": 14660000,
        "color": "cyan",
        "shape": "circle",
        "width": 3,
    },
]

# =============================================================================
# Track Files
# =============================================================================

# ATAC
ATAC_EARLY = "/share/home/ychi/dev/sc-embryo/mm10_bigwig/ATAC_early_2cell.mm10.bw"
ATAC_LATE = "/share/home/ychi/dev/sc-embryo/mm10_bigwig/ATAC_late_2cell.mm10.bw"

# H3K27ac (multiple replicates)
H3K27AC_EARLY_REPS = [
    "/shared/ychi/repo/sc_embryo/notebooks1/Liu2024_H3K27ac/results/bigwig/Early_2cell_rep1.bw",
    "/shared/ychi/repo/sc_embryo/notebooks1/Liu2024_H3K27ac/results/bigwig/Early_2cell_rep2.bw",
]
H3K27AC_LATE_REPS = [
    "/shared/ychi/repo/sc_embryo/notebooks1/Liu2024_H3K27ac/results/bigwig/Late_2cell_rep1_run1.bw",
    "/shared/ychi/repo/sc_embryo/notebooks1/Liu2024_H3K27ac/results/bigwig/Late_2cell_rep2_run1.bw",
]

# H3K4me3 (multiple replicates)
H3K4ME3_EARLY_REPS = [
    "/shared/ychi/repo/sc_embryo/notebooks1/Zhang2016_GSE71434/results/bigwig/Early_2cell_rep1.bw",
    "/shared/ychi/repo/sc_embryo/notebooks1/Zhang2016_GSE71434/results/bigwig/Early_2cell_rep2.bw",
]
H3K4ME3_LATE_REPS = [
    "/shared/ychi/repo/sc_embryo/notebooks1/Zhang2016_GSE71434/results/bigwig/Late_2cell_rep1.bw",
    "/shared/ychi/repo/sc_embryo/notebooks1/Zhang2016_GSE71434/results/bigwig/Late_2cell_rep2.bw",
]

# CTCF (multiple replicates)
CTCF_EARLY_REPS = [
    "/share/home/ychi/dev/sc-embryo/mm10_bigwig/GSM6031218_CTCF_Early2C_rep1_clean_100bp_rpkm.mm10.bw",
    "/share/home/ychi/dev/sc-embryo/mm10_bigwig/GSM6031219_CTCF_Early2C_rep2_clean_100bp_rpkm.mm10.bw",
]
CTCF_LATE_REPS = [
    "/share/home/ychi/dev/sc-embryo/mm10_bigwig/GSM6031220_CTCF_Late2C_rep1_clean_100bp_rpkm.mm10.bw",
    "/share/home/ychi/dev/sc-embryo/mm10_bigwig/GSM6031221_CTCF_Late2C_rep2_clean_100bp_rpkm.mm10.bw",
]

# SMC3 (multiple replicates)
SMC3_EARLY_REPS = [
    "/share/home/ychi/dev/sc-embryo/mm10_bigwig/GSM6031232_SMC3_Early2C_rep1_clean_100bp_rpkm.mm10.bw",
    "/share/home/ychi/dev/sc-embryo/mm10_bigwig/GSM6031233_SMC3_Early2C_rep2_clean_100bp_rpkm.mm10.bw",
]
SMC3_LATE_REPS = [
    "/share/home/ychi/dev/sc-embryo/mm10_bigwig/GSM9186167_SMC3_Late2C_rep1.mm10.bw",
    "/share/home/ychi/dev/sc-embryo/mm10_bigwig/GSM9186168_SMC3_Late2C_rep2.mm10.bw",
]

# Gene track
GENE_GTF = (
    "/share/home/ychi/dev/sc-embryo/notebooks/E_loop_analysis/"
    "7_celltype_loops/6.subset_merged_tiling_tracks/"
    "6.subset_merged_tiling_tracks.top_loops.genes.gtf"
)

# =============================================================================
# Build Track List
# =============================================================================

# For tracks with multiple replicates, we'll use the first replicate.
# NOTE: CTCF and SMC3 bigWig files have been fixed and are now included.

tracks = [
    # ATAC
    BigWigTrack(ATAC_EARLY, bin_size=1000, display_name="ATAC early"),
    BigWigTrack(ATAC_LATE, bin_size=1000, display_name="ATAC late"),
    # H3K27ac
    BigWigTrack(H3K27AC_EARLY_REPS[0], bin_size=1000, display_name="H3K27ac early"),
    BigWigTrack(H3K27AC_LATE_REPS[0], bin_size=1000, display_name="H3K27ac late"),
    # H3K4me3
    BigWigTrack(H3K4ME3_EARLY_REPS[0], bin_size=1000, display_name="H3K4me3 early"),
    BigWigTrack(H3K4ME3_LATE_REPS[0], bin_size=1000, display_name="H3K4me3 late"),
    # CTCF
    BigWigTrack(CTCF_EARLY_REPS[0], bin_size=1000, display_name="CTCF early"),
    BigWigTrack(CTCF_LATE_REPS[0], bin_size=1000, display_name="CTCF late"),
    # SMC3
    BigWigTrack(SMC3_EARLY_REPS[0], bin_size=1000, display_name="SMC3 early"),
    BigWigTrack(SMC3_LATE_REPS[0], bin_size=1000, display_name="SMC3 late"),
    # Gene
    GeneTrack(GENE_GTF, mode="canonical", show_labels=True, label_font_size=9),
]

# =============================================================================
# Track Heights
# =============================================================================

# Full 11 tracks: 10 BigWig tracks + 1 Gene track
# Gene track gets extra height (4.0) for staggered gene labels
# BigWig tracks reduced slightly (1.0) to give more space to heatmap
track_heights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 4.0]

# =============================================================================
# Plot
# =============================================================================

print(f"Reproducing Figure 1 (loop {LOOP_ID})...")
print(f"Region: {REGION}")
print(f"Early mcool: {early_mcool_uri}")
print(f"Late mcool: {late_mcool_uri}")
print(f"Output: {OUTPUT_PNG}")

fig = plot_cool_with_tracks(
    early_mcool_uri,
    coolp2=late_mcool_uri,
    tracks=tracks,
    region=REGION,
    track_heights=track_heights,
    heatmap_height=12,  # Increased from 9 to 12 for taller heatmap
    square_heatmap=True,  # Enable square heatmap
    title=(
        f"Loop {LOOP_ID} | lower=early matched merge | upper=late matched merge<br>"
        "Track order: assay early row then assay late row; Gene last"
    ),
    balance=False,
    cmap="viridis",
    adjust_brightness=True,
    marks=marks,
)

# Update layout to match original
# Note: when square_heatmap=True, heatmap_height is auto-calculated
# Increase figure height to accommodate taller heatmap
fig.update_layout(
    width=1100,
    height=1600,  # Increased from 1450 to 1600
    title_x=0.5,
    margin=dict(l=140, r=60, t=90, b=40),
)

# Save to PNG
fig.write_image(str(OUTPUT_PNG))
print(f"Figure saved to {OUTPUT_PNG}")
