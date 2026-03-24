#!/usr/bin/env python
"""
Visualize the effect of endpoint_mode='half-interval' vs 'full'.

Usage:
    python visualize_endpoint_mode_comparison.py
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.io import write_image
from pathlib import Path
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


def build_time_grid(pseudotime, numPts=14, method="linspace", endpoint_mode="half-interval"):
    """Simplified version of _build_time_grid for visualization."""
    if method == "linspace":
        if endpoint_mode == "half-interval":
            interval = 1.0 / (numPts - 1)
            half_interval = interval / 2
            return np.linspace(half_interval, 1.0 - half_interval, numPts)
        else:
            return np.linspace(0.0, 1.0, numPts)
    elif method == "quantile":
        if endpoint_mode == "half-interval":
            interval = 1.0 / (numPts - 1)
            half_interval = interval / 2
            pt_sorted = np.sort(pseudotime)
            n = len(pt_sorted)
            q_min_idx = int(np.searchsorted(pt_sorted, half_interval * pt_sorted.max()))
            q_max_idx = int(np.searchsorted(pt_sorted, (1 - half_interval) * pt_sorted.max()))
            q_min = max(0.0, q_min_idx / (n - 1))
            q_max = min(1.0, q_max_idx / (n - 1))
            quantiles = np.linspace(q_min, q_max, numPts)
            return np.quantile(pseudotime, quantiles)
        else:
            quantiles = np.linspace(0, 1, numPts)
            return np.quantile(pseudotime, quantiles)
    else:
        raise ValueError(f"Unknown method: {method}")


def make_endpoint_mode_comparison(data, methods, numPts=12, scaling="z_score"):
    """
    Compare endpoint_mode='full' vs 'half-interval' for each method.
    """
    n_methods = len(methods)
    
    if scaling == "z_score":
        raw_key = "raw_z_score_scaled"
        colorscale = "RdBu_r"
        zmin, zmax = -3, 3
    else:
        raw_key = "raw_zero_one_scaled"
        colorscale = "inferno"
        zmin, zmax = 0, 1
    
    raw_matrix = data[raw_key]
    cell_pseudotime = data["cell_pseudotime"]
    loop_ids = data["loop_ids"]
    original_grid = data["time_grid"]
    smooth_matrix_orig = data["smooth_z_score_scaled"] if scaling == "z_score" else data["smooth_zero_one_scaled"]
    
    cell_order = np.argsort(cell_pseudotime)
    raw_sorted = raw_matrix[:, cell_order]
    n_cells = len(cell_pseudotime)
    cell_ranks = np.arange(n_cells) / (n_cells - 1)
    pseudotime_sorted = np.sort(cell_pseudotime)
    
    # Create figure: 2 rows per method (full vs half-interval)
    fig = make_subplots(
        rows=n_methods * 2, cols=2,
        shared_yaxes=True,
        shared_xaxes=False,
        vertical_spacing=0.03,
        horizontal_spacing=0.06,
        column_widths=[0.55, 0.45],
        subplot_titles=[f"{m}: {ep} (raw)" for m in methods for ep in ["full", "half-int"]] +
                       [f"{m}: {ep} ({len(build_time_grid(cell_pseudotime, numPts, m, ep.replace('half-int', 'half-interval')))} pts)" 
                        for m in methods for ep in ["full", "half-int"]]
    )
    
    row_idx = 1
    for method in methods:
        for endpoint_mode in ["full", "half-interval"]:
            grid = build_time_grid(cell_pseudotime, numPts, method, endpoint_mode)
            
            # Left: Raw with grid lines
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
            
            # Grid lines
            for grid_t in grid:
                rank_pos = np.searchsorted(pseudotime_sorted, grid_t) / (n_cells - 1)
                fig.add_vline(x=rank_pos, line_dash="dash", line_color="black", line_width=1.5,
                              row=row_idx, col=1)
            
            # Right: GP-smoothed
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
                    showscale=row_idx == n_methods * 2,
                    colorbar=dict(title="Z-score" if scaling == "z_score" else "Scaled", x=1.02) if row_idx == n_methods * 2 else None,
                    hovertemplate=None,
                ),
                row=row_idx, col=2
            )
            
            row_idx += 1
    
    fig_height = max(300, len(loop_ids) * 8 * n_methods * 2)
    fig_width = 1400
    
    fig.update_layout(
        height=fig_height,
        width=fig_width,
        showlegend=False,
        margin=dict(l=60, r=80, t=80, b=40),
        title_text="Endpoint Mode Comparison: 'full' vs 'half-interval'",
        title_y=0.995,
        title_font_size=20,
    )
    
    fig.update_xaxes(title_text="Cell rank", row=n_methods * 2, col=1)
    fig.update_xaxes(
        title_text="Grid point index",
        row=n_methods * 2, col=2,
        tickmode='linear', tick0=0, dtick=1
    )
    fig.update_yaxes(title_text="Loop index")
    fig.update_annotations(font_size=12, font_weight="bold")
    
    return fig


def main():
    data_path = Path(__file__).parent / "test_data.npz"
    if not data_path.exists():
        print(f"Error: Test data not found at {data_path}")
        return
    
    print("=" * 70)
    print("Endpoint Mode Comparison Visualization")
    print("=" * 70)
    
    data = load_test_data(data_path)
    pseudotime = data["cell_pseudotime"]
    print(f"Loops: {len(data['loop_ids'])}, Cells: {len(pseudotime)}")
    
    numPts = 12
    methods = ["linspace", "quantile"]
    
    output_dir = Path(__file__).parent / "results" / "png_endpoint_mode"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Z-score comparison
    print("\nGenerating endpoint mode comparison (Z-score)...")
    fig = make_endpoint_mode_comparison(data, methods, numPts, scaling="z_score")
    output_path = output_dir / "endpoint_mode_comparison.png"
    write_image(fig, output_path, scale=2, width=1400, height=fig.layout.height)
    print(f"Saved: {output_path}")
    
    # Zero-one comparison
    print("Generating endpoint mode comparison (0-1 scaled)...")
    fig = make_endpoint_mode_comparison(data, methods, numPts, scaling="zero_one")
    output_path = output_dir / "endpoint_mode_comparison_zero_one.png"
    write_image(fig, output_path, scale=2, width=1400, height=fig.layout.height)
    print(f"Saved: {output_path}")
    
    print("\n" + "=" * 70)
    print("EXPLANATION:")
    print("=" * 70)
    print("Each method has two rows:")
    print("  Row 1 (full): endpoints at t=0 and t=1 (legacy)")
    print("  Row 2 (half-interval): endpoints at t=half_interval and t=1-half_interval")
    print()
    print("Benefits of half-interval:")
    print("  - Avoids sparse cells at endpoints (t=0 and t=1 have only 1 cell each)")
    print("  - Reduces GP extrapolation uncertainty")
    print("  - More reliable interpolation at endpoints")
    print("=" * 70)


if __name__ == "__main__":
    main()
