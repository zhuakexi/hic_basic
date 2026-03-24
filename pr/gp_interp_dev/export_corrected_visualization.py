#!/usr/bin/env python
"""
Corrected visualization with proper layout and explanations.

Fixes:
1. Each row shows the SAME method's raw (left) and smooth (right)
2. Clear labeling of what each panel shows
3. Add quantile method comparison
"""

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.io import write_image
from pathlib import Path
from scipy.stats import gaussian_kde
from scipy.interpolate import interp1d


def load_test_data(data_path):
    data = np.load(data_path, allow_pickle=True)
    return {
        "raw_matrix": data["raw_matrix"],
        "smooth_matrix": data["smooth_matrix"],
        "raw_z_score_scaled": data["raw_z_score_scaled"],
        "smooth_z_score_scaled": data["smooth_z_score_scaled"],
        "raw_zero_one_scaled": data["raw_zero_one_scaled"],
        "smooth_zero_one_scaled": data["smooth_zero_one_scaled"],
        "time_grid": data["time_grid"],
        "cell_pseudotime": data["cell_pseudotime"],
        "loop_ids": data["loop_ids"],
    }


def build_time_grid(pseudotime, numPts=14, method="linspace"):
    if method == "linspace":
        return np.linspace(0.0, 1.0, numPts)
    elif method == "quantile":
        quantiles = np.linspace(0, 1, numPts)
        return np.quantile(pseudotime, quantiles)
    elif method == "density":
        pt_sorted = np.sort(pseudotime)
        kde = gaussian_kde(pt_sorted)
        x_fine = np.linspace(0, 1, 1000)
        density = kde(x_fine)
        cdf = np.cumsum(density)
        cdf = cdf / cdf[-1]
        return np.interp(np.linspace(0, 1, numPts), cdf, x_fine)
    elif method == "custom":
        return np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    else:
        raise ValueError(f"Unknown method: {method}")


def make_corrected_comparison_figure(data, methods_data, scaling="z_score"):
    """
    CORRECTED: Each row shows ONE method with raw (left) and smooth (right).
    """
    n_methods = len(methods_data)
    
    if scaling == "z_score":
        raw_key = "raw_z_score_scaled"
        smooth_key = "smooth_z_score_scaled"
        colorscale = "RdBu_r"
        zmin, zmax = -3, 3
    else:
        raw_key = "raw_zero_one_scaled"
        smooth_key = "smooth_zero_one_scaled"
        colorscale = "inferno"
        zmin, zmax = 0, 1
    
    raw_matrix = data[raw_key]
    cell_pseudotime = data["cell_pseudotime"]
    loop_ids = data["loop_ids"]
    original_grid = data["time_grid"]
    smooth_matrix_orig = data[smooth_key] if scaling == "z_score" else data[smooth_key]
    
    cell_order = np.argsort(cell_pseudotime)
    raw_sorted = raw_matrix[:, cell_order]
    n_cells = len(cell_pseudotime)
    cell_ranks = np.arange(n_cells) / (n_cells - 1)
    pseudotime_sorted = np.sort(cell_pseudotime)
    
    fig_height = max(200, len(loop_ids) * 12 * n_methods)
    fig_width = 1400
    
    fig = make_subplots(
        rows=n_methods, cols=2,
        shared_yaxes=True,
        shared_xaxes=False,
        vertical_spacing=0.04,
        horizontal_spacing=0.06,
        column_widths=[0.55, 0.45],
        subplot_titles=[f"{name}: Raw (grid lines)" for name, _ in methods_data] + 
                       [f"{name}: GP-smoothed ({len(grid)} pts)" for name, grid in methods_data]
    )
    
    for row_idx, (method_name, grid) in enumerate(methods_data, 1):
        # LEFT: Raw data with grid lines for THIS method
        fig.add_trace(
            go.Heatmap(
                z=raw_sorted,
                x=cell_ranks,
                y=list(range(len(loop_ids))),
                colorscale=colorscale,
                zmin=zmin,
                zmax=zmax,
                showscale=False,
                hovertemplate=None,
            ),
            row=row_idx, col=1
        )
        
        # Grid lines for THIS method
        for grid_t in grid:
            rank_pos = np.searchsorted(pseudotime_sorted, grid_t) / (n_cells - 1)
            fig.add_vline(x=rank_pos, line_dash="dash", line_color="black", line_width=1.5,
                          row=row_idx, col=1)
        
        # RIGHT: GP-smoothed for THIS method
        smooth_interp = np.zeros((smooth_matrix_orig.shape[0], len(grid)))
        for i in range(smooth_matrix_orig.shape[0]):
            f = interp1d(original_grid, smooth_matrix_orig[i], kind='cubic',
                         bounds_error=False, fill_value='extrapolate')
            smooth_interp[i] = f(grid)
        
        grid_indices = np.arange(len(grid))
        
        fig.add_trace(
            go.Heatmap(
                z=smooth_interp,
                x=grid_indices,
                y=list(range(len(loop_ids))),
                colorscale=colorscale,
                zmin=zmin,
                zmax=zmax,
                showscale=row_idx == n_methods,
                colorbar=dict(title="Z-score" if scaling == "z_score" else "Scaled", x=1.02) if row_idx == n_methods else None,
                hovertemplate=None,
            ),
            row=row_idx, col=2
        )
    
    fig.update_layout(
        height=fig_height,
        width=fig_width,
        showlegend=False,
        margin=dict(l=60, r=80, t=80, b=40),
        title_text="GP Interpolation: Method Comparison (CORRECTED LAYOUT)",
        title_y=0.995,
        title_font_size=20,
    )
    
    fig.update_xaxes(title_text="Cell rank (sorted by pseudotime)", row=n_methods, col=1)
    fig.update_xaxes(
        title_text="Grid point index", 
        row=n_methods, col=2,
        tickmode='linear', tick0=0, dtick=1
    )
    fig.update_yaxes(title_text="Loop index")
    fig.update_annotations(font_size=13, font_weight="bold")
    
    return fig


def main():
    data_path = Path(__file__).parent / "test_data.npz"
    if not data_path.exists():
        print(f"Error: Test data not found at {data_path}")
        return
    
    print("=" * 70)
    print("GP Interpolation: CORRECTED Visualization")
    print("=" * 70)
    
    data = load_test_data(data_path)
    pseudotime = data["cell_pseudotime"]
    print(f"Loops: {len(data['loop_ids'])}, Cells: {len(pseudotime)}")
    
    numPts = 12
    methods = {
        "linspace": lambda: build_time_grid(pseudotime, numPts, "linspace"),
        "quantile": lambda: build_time_grid(pseudotime, numPts, "quantile"),
        "density": lambda: build_time_grid(pseudotime, numPts, "density"),
        "custom": lambda: build_time_grid(pseudotime, numPts, "custom"),
    }
    
    methods_data = [(name, fn()) for name, fn in methods.items()]
    
    print("\nGrid methods:")
    for name, grid in methods_data:
        print(f"  {name}: {len(grid)} points, range [{grid.min():.3f}, {grid.max():.3f}]")
    
    output_dir = Path(__file__).parent / "results" / "png_corrected"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Z-score comparison
    print("\nGenerating CORRECTED comparison figure (Z-score)...")
    fig = make_corrected_comparison_figure(data, methods_data, scaling="z_score")
    output_path = output_dir / "gp_interp_methods_comparison_CORRECTED.png"
    write_image(fig, output_path, scale=2, width=1400, height=fig.layout.height)
    print(f"Saved: {output_path}")
    
    # Zero-one comparison
    print("Generating CORRECTED comparison figure (0-1 scaled)...")
    fig = make_corrected_comparison_figure(data, methods_data, scaling="zero_one")
    output_path = output_dir / "gp_interp_methods_comparison_zero_one_CORRECTED.png"
    write_image(fig, output_path, scale=2, width=1400, height=fig.layout.height)
    print(f"Saved: {output_path}")
    
    print("\n" + "=" * 70)
    print("EXPLANATION OF FIXES:")
    print("=" * 70)
    print("1. Each row now shows ONE method:")
    print("   - Left: Raw data with grid lines for that method")
    print("   - Right: GP-smoothed result for that method")
    print()
    print("2. All 4 methods are shown:")
    print("   - Row 1: linspace")
    print("   - Row 2: quantile")
    print("   - Row 3: density")
    print("   - Row 4: custom")
    print()
    print("3. Grid lines show WHERE each method places its grid points")
    print("   - Lines are at rank positions, not pseudotime values")
    print("   - Endpoints (0 and 1) have few cells but are included")
    print()
    print("4. Right panel uses uniform grid indices (0,1,2,...)")
    print("   - This shows the GP output at each grid point equally spaced")
    print("=" * 70)


if __name__ == "__main__":
    main()
