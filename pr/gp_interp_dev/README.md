# GP Interpolation Development Environment

This directory contains test data and visualization tools for extending `gp_interpolate()` with alternative time point selection strategies.

## Files

- `test_data.npz` - Small-scale test dataset (50 loops, 740 cells)
- `test_loop_table.tsv` - Loop metadata for the sampled loops
- `extract_test_data.py` - Script to regenerate test data from series 1
- `visualize_gp_interp.py` - Standalone visualization script
- `requirements.txt` - Python dependencies

## Quick Start

### 1. Extract test data (if needed)

```bash
cd /share/home/ychi/dev/hic_basic/pr/gp_interp_dev
python extract_test_data.py
```

This extracts 50 loops (stratified by peak_grid_index) from series 1, keeping all 740 cells.

### 2. Visualize GP interpolation results

```bash
# Default: view current GP result with 12 grid points
python visualize_gp_interp.py

# Compare different numPts settings
python visualize_gp_interp.py --numPts 8 12 20 50

# Test quantile-based time grid
python visualize_gp_interp.py --method quantile --numPts 12

# Save to custom output file
python visualize_gp_interp.py --output my_comparison.html
```

### 3. Test new time point selection methods

Edit `visualize_gp_interp.py` and add your method to `generate_time_grid()`:

```python
def generate_time_grid(method="uniform", numPts=14, cell_pseudotime=None, custom_times=None):
    if method == "uniform":
        return np.linspace(0.0, 1.0, numPts)
    elif method == "quantile":
        quantiles = np.linspace(0.0, 1.0, numPts)
        return np.quantile(cell_pseudotime, quantiles)
    elif method == "your_new_method":
        # Your implementation here
        return your_time_grid
```

## Data Structure

**test_data.npz** contains:

| Array | Shape | Description |
|-------|-------|-------------|
| `raw_matrix` | (50, 740) | Loop × cell raw signal |
| `smooth_matrix` | (50, 12) | Loop × GP grid smoothed |
| `*_z_score_scaled` | varies | Z-score normalized versions |
| `*_zero_one_scaled` | varies | [0,1] scaled versions |
| `time_grid` | (12,) | GP grid points [0.0, 0.09, ..., 1.0] |
| `cell_pseudotime` | (740,) | Cell pseudotime values |
| `loop_ids` | (50,) | Loop identifiers |

## Extending gp_interpolate()

To add support for custom time grids in `hic_basic/interpolate.py`:

1. Add a new parameter `time_grid: Optional[np.ndarray] = None` to `gp_interpolate()`
2. Modify the grid generation logic:
   ```python
   if time_grid is not None:
       grid_values = np.asarray(time_grid)
       numPts = len(grid_values)
   else:
       grid_values = np.linspace(0.0, 1.0, numPts)
   ```
3. Validate: `assert grid_values.min() >= 0 and grid_values.max() <= 1`
4. Test with this development dataset before deploying

## Dependencies

```
numpy
pandas
plotly
scipy
```

Run with the project environment: `./run_in_env.sh python visualize_gp_interp.py`
