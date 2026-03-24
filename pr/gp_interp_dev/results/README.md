# GP Interpolation Grid Method Comparison Results

Generated: 2026-03-23

## Test Configuration

- **Dataset**: 50 loops × 740 cells
- **Pseudotime range**: [0.0, 1.0]
- **Target numPts**: 12 (except custom which uses 11 fixed points)

## Generated Files

| File | Description | Size |
|------|-------------|------|
| `gp_interp_all_methods_comparison.html` | Main comparison heatmap (Z-score scaled) | 5.3 MB |
| `gp_interp_all_methods_zero_one.html` | Comparison heatmap ([0,1] scaled) | 5.3 MB |
| `gp_interp_grid_distribution.html` | Grid point distribution visualization | 4.5 MB |
| `gp_interp_summary_stats.csv` | Summary statistics table | 507 B |

## Methods Tested

| Method | n_points | min_spacing | max_spacing | mean_spacing | uniformity |
|--------|----------|-------------|-------------|--------------|------------|
| **linspace** | 12 | 0.091 | 0.091 | 0.091 | 1.000 |
| **quantile** | 12 | 0.030 | 0.185 | 0.091 | 0.469 |
| **density** | 12 | 0.047 | 0.218 | 0.091 | 0.492 |
| **custom** | 11 | 0.100 | 0.100 | 0.100 | 1.000 |
| ~~optimal~~ | N/A | N/A | N/A | N/A | N/A |

**Note**: `optimal` method requires `optbinning` package (not installed).

## Method Characteristics

### linspace (uniformity = 1.0)
- **Evenly spaced** grid points
- Best for: Uniform sampling, baseline comparison
- Grid: `[0.0, 0.091, 0.182, ..., 0.909, 1.0]`

### quantile (uniformity = 0.47)
- **Equal cells per bin** - grid points at quantiles of cell distribution
- Best for: Ensuring equal representation of cells in each interval
- Places more points in sparse regions, fewer in dense regions

### density (uniformity = 0.49)
- **Density-adaptive** - more points in dense cell regions
- Best for: Capturing fine-grained dynamics where most cells are
- Uses KDE to estimate cell density, then places points accordingly

### custom (uniformity = 1.0)
- **User-defined** grid points
- Best for: Domain-specific requirements, known biological milestones
- Example: `[0.0, 0.1, 0.2, ..., 0.9, 1.0]`

## Visual Inspection Guide

### Main Comparison Heatmap (`gp_interp_all_methods_comparison.html`)

Each row shows:
- **Left panel**: Raw data (740 cells) with vertical dashed lines showing grid points
- **Right panel**: GP-smoothed data interpolated to the method's grid

**What to look for:**
1. Do the grid lines align with biological features in the raw data?
2. Does the smoothed output capture the key dynamics?
3. Are there any artifacts from the grid placement?

### Grid Distribution Figure (`gp_interp_grid_distribution.html`)

Shows grid point placements overlaid on cell density histogram.

**What to look for:**
1. **linspace**: Evenly spaced regardless of cell density
2. **quantile**: More points where cells are sparse (compensates for imbalance)
3. **density**: More points where cells are dense (captures detail)
4. **custom**: User-defined placement

### Zero-One Scaled Heatmap (`gp_interp_all_methods_zero_one.html`)

Same as main comparison but with [0,1] color scaling.

**What to look for:**
1. Relative signal intensity across methods
2. Preservation of peak locations
3. Smoothing artifacts

## Recommendations

Based on the test results:

| Use Case | Recommended Method |
|----------|-------------------|
| General purpose | `linspace` (simple, predictable) |
| Uneven cell distribution | `quantile` (equal cells per bin) |
| Dense cell clusters | `density` (capture fine detail) |
| Known biological milestones | `custom` (domain-specific) |
| Data-driven optimization | `optimal` (requires optbinning) |

## How to View Results

Open the HTML files in a web browser:

```bash
# On a server with browser access
xdg-open results/gp_interp_all_methods_comparison.html

# Or copy to local machine and open
scp user@server:/path/to/results/*.html ./
open gp_interp_all_methods_comparison.html
```

## Reproduce Results

```bash
cd /share/home/ychi/dev/hic_basic/pr/gp_interp_dev
./../../run_in_env.sh python test_all_grid_methods.py
```

## Next Steps

1. **Visual inspection**: Review the heatmaps and note any artifacts
2. **Quantitative evaluation**: Consider adding metrics like:
   - Reconstruction error (raw vs smoothed)
   - Smoothness metrics
   - Biological feature preservation
3. **Method refinement**: Adjust parameters based on findings
