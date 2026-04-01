# Debug Data For `plot_cool_with_tracks`: sc-embryo series 6 figure 1

This document lists the exact real-data inputs used to generate:

- `notebooks/E_loop_analysis/7_celltype_loops/6.subset_merged_tiling_tracks/6.subset_merged_tiling_tracks.loop01.325348.tracks.png`

The goal is to let `hic_basic` reproduce and debug the function with the same
inputs, not a synthetic mockup.

## Output Figure

- Figure label: `series 6 figure 1`
- Output PNG:
  `/share/home/ychi/dev/sc-embryo/notebooks/E_loop_analysis/7_celltype_loops/6.subset_merged_tiling_tracks/6.subset_merged_tiling_tracks.loop01.325348.tracks.png`

## Heatmap Inputs

The tiled heatmap uses:

- lower triangle = early matched merged mcool
- upper triangle = late matched merged mcool

Exact URIs:

- early:
  `/share/home/ychi/dev/sc-embryo/notebooks/E_loop_analysis/7_celltype_loops/6.subset_merged_tiling_tracks/6.subset_merged_tiling_tracks.early.matched.mcool::resolutions/20000`
- late:
  `/share/home/ychi/dev/sc-embryo/notebooks/E_loop_analysis/7_celltype_loops/6.subset_merged_tiling_tracks/6.subset_merged_tiling_tracks.late.matched.mcool::resolutions/20000`

Matched merge summary:

- early raw count sum: `4201729`
- late raw count sum: `4213905`
- target count sum: `4201729`
- early matched count sum: `4201729`
- late matched count sum: `4201746`

Source:

- `/share/home/ychi/dev/sc-embryo/notebooks/E_loop_analysis/7_celltype_loops/6.subset_merged_tiling_tracks/6.subset_merged_tiling_tracks.merge_manifest.tsv`

## Cell Groups Behind The Merge

This figure uses an `8 vs 8` matched-cell merge.

Selected early cells:

- `0h1c12`
- `0h1c17`
- `1h1a44b2`
- `1h1a39b2`
- `5h1a17b1`
- `4h1a23b1`
- `4h1a12b2`
- `C1_0ha36b2`

Selected late matched cells:

- `B1_0ha20b2`
- `3h1a11b2`
- `1h1a17b1`
- `3h1a6b1`
- `6h1a3b1`
- `1h1a5b2`
- `2h1a22b2`
- `0h1c22`

Source:

- `/share/home/ychi/dev/sc-embryo/notebooks/E_loop_analysis/7_celltype_loops/6.subset_merged_tiling_tracks/6.subset_merged_tiling_tracks.selection_manifest.json`
- `/share/home/ychi/dev/sc-embryo/notebooks/E_loop_analysis/7_celltype_loops/6.subset_merged_tiling_tracks/6.subset_merged_tiling_tracks.group_manifest.tsv`

## Loop Being Marked

This figure is for the top loop:

- `loop_id`: `325348`
- coordinates:
  `chr10:14640000 - chr10:14660000`
- plotted region:
  `chr10:14240000-15080000`

At 20 kb resolution, the mark spans one bin per anchor.

## Mark Inputs

The current `sc-embryo` call passes two marks so the same loop is visible in
both the lower and upper triangles:

```python
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
        "label": "loop 325348",
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
```

Note:

- For this specific loop, the mirror is identical because the plotted loop is
  on the diagonal.
- The real code path still passes two entries because the same series also
  handles non-diagonal loops with explicit mirrored marks.

## Track Order

Tracks are passed in this exact order:

1. `ATAC early`
2. `ATAC late`
3. `H3K27ac early`
4. `H3K27ac late`
5. `H3K4me3 early`
6. `H3K4me3 late`
7. `CTCF early`
8. `CTCF late`
9. `SMC3 early`
10. `SMC3 late`
11. `Gene`

Track heights:

```python
[1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 2.6]
```

## Track Files

### ATAC

- early:
  `/share/home/ychi/dev/sc-embryo/mm10_bigwig/ATAC_early_2cell.mm10.bw`
- late:
  `/share/home/ychi/dev/sc-embryo/mm10_bigwig/ATAC_late_2cell.mm10.bw`

### H3K27ac

- early replicate paths:
  `/shared/ychi/repo/sc_embryo/notebooks1/Liu2024_H3K27ac/results/bigwig/Early_2cell_rep1.bw`
  `/shared/ychi/repo/sc_embryo/notebooks1/Liu2024_H3K27ac/results/bigwig/Early_2cell_rep2.bw`
- late replicate paths:
  `/shared/ychi/repo/sc_embryo/notebooks1/Liu2024_H3K27ac/results/bigwig/Late_2cell_rep1_run1.bw`
  `/shared/ychi/repo/sc_embryo/notebooks1/Liu2024_H3K27ac/results/bigwig/Late_2cell_rep2_run1.bw`

### H3K4me3

- early replicate paths:
  `/shared/ychi/repo/sc_embryo/notebooks1/Zhang2016_GSE71434/results/bigwig/Early_2cell_rep1.bw`
  `/shared/ychi/repo/sc_embryo/notebooks1/Zhang2016_GSE71434/results/bigwig/Early_2cell_rep2.bw`
- late replicate paths:
  `/shared/ychi/repo/sc_embryo/notebooks1/Zhang2016_GSE71434/results/bigwig/Late_2cell_rep1.bw`
  `/shared/ychi/repo/sc_embryo/notebooks1/Zhang2016_GSE71434/results/bigwig/Late_2cell_rep2.bw`

### CTCF

- early replicate paths:
  `/share/home/ychi/dev/sc-embryo/mm10_bigwig/GSM6031218_CTCF_Early2C_rep1_clean_100bp_rpkm.mm10.bw`
  `/share/home/ychi/dev/sc-embryo/mm10_bigwig/GSM6031219_CTCF_Early2C_rep2_clean_100bp_rpkm.mm10.bw`
- late replicate paths:
  `/share/home/ychi/dev/sc-embryo/mm10_bigwig/GSM6031220_CTCF_Late2C_rep1_clean_100bp_rpkm.mm10.bw`
  `/share/home/ychi/dev/sc-embryo/mm10_bigwig/GSM6031221_CTCF_Late2C_rep2_clean_100bp_rpkm.mm10.bw`

### SMC3

- early replicate paths:
  `/share/home/ychi/dev/sc-embryo/mm10_bigwig/GSM6031232_SMC3_Early2C_rep1_clean_100bp_rpkm.mm10.bw`
  `/share/home/ychi/dev/sc-embryo/mm10_bigwig/GSM6031233_SMC3_Early2C_rep2_clean_100bp_rpkm.mm10.bw`
- late replicate paths:
  `/share/home/ychi/dev/sc-embryo/mm10_bigwig/GSM9186167_SMC3_Late2C_rep1.mm10.bw`
  `/share/home/ychi/dev/sc-embryo/mm10_bigwig/GSM9186168_SMC3_Late2C_rep2.mm10.bw`

### Gene Track

The original external GTF mount was not visible inside `run_in_env.sh`, so the
real plotting run used a local extracted GTF:

- `/share/home/ychi/dev/sc-embryo/notebooks/E_loop_analysis/7_celltype_loops/6.subset_merged_tiling_tracks/6.subset_merged_tiling_tracks.top_loops.genes.gtf`

`GeneTrack` parameters:

```python
GeneTrack(
    LOCAL_GENE_GTF,
    mode="canonical",
    show_labels=True,
    label_font_size=9,
)
```

## Plotting Parameters

Exact plotting call shape:

```python
fig = plot_cool_with_tracks(
    early_mcool_uri,
    coolp2=late_mcool_uri,
    tracks=tracks,
    region="chr10:14240000-15080000",
    track_heights=[1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 2.6],
    heatmap_height=9,
    title=(
        "Loop 325348 | lower=early matched merge | upper=late matched merge"
        "<br>Track order: assay early row then assay late row; Gene last"
    ),
    balance=False,
    cmap="viridis",
    adjust_brightness=True,
    marks=marks,
)
```

Layout updates after plotting:

```python
fig.update_layout(
    width=1100,
    height=1450,
    title_x=0.5,
    margin=dict(l=140, r=60, t=90, b=40),
)
```

In `sc-embryo`, `template=h.template` is also applied.

## Supporting Files In sc-embryo

Useful manifests:

- `/share/home/ychi/dev/sc-embryo/notebooks/E_loop_analysis/7_celltype_loops/6.subset_merged_tiling_tracks/6.subset_merged_tiling_tracks.figure_manifest.tsv`
- `/share/home/ychi/dev/sc-embryo/notebooks/E_loop_analysis/7_celltype_loops/6.subset_merged_tiling_tracks/6.subset_merged_tiling_tracks.track_manifest.tsv`
- `/share/home/ychi/dev/sc-embryo/notebooks/E_loop_analysis/7_celltype_loops/6.subset_merged_tiling_tracks/6.subset_merged_tiling_tracks.merge_manifest.tsv`
- `/share/home/ychi/dev/sc-embryo/notebooks/E_loop_analysis/7_celltype_loops/6.subset_merged_tiling_tracks/6.subset_merged_tiling_tracks.selection_manifest.json`

Main script:

- `/share/home/ychi/dev/sc-embryo/notebooks/E_loop_analysis/7_celltype_loops/6.subset_merged_tiling_tracks/6.subset_merged_tiling_tracks.py`

## Practical Debug Target

If `hic_basic` wants one real figure to debug:

1. Use the exact early/late matched mcool URIs above
2. Use region `chr10:14240000-15080000`
3. Use the exact 11 tracks and heights above
4. Use the two `marks` entries above
5. Render to a new test PNG and compare against the existing sc-embryo figure

This is the full real-data payload for figure 1.

## Known Data Issues (Resolved)

### CTCF and SMC3 BigWig Files - FIXED

The CTCF and SMC3 bigWig files were previously corrupted (~35KB with constant values).
They have since been fixed and now contain valid signal data:

| Track | File | Status |
|-------|------|--------|
| CTCF early/late | `GSM6031218_*_clean_100bp_rpkm.mm10.bw` | ✅ Fixed |
| SMC3 early/late | `GSM6031232_*_clean_100bp_rpkm.mm10.bw` | ✅ Fixed |

The reproduce script now includes all 11 tracks as originally specified.
