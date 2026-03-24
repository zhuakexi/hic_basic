#!/usr/bin/env python
"""
Export GP interpolation comparison results to PNG images.

This script converts the HTML visualization results to high-resolution PNG images
for easy sharing and publication.

Usage:
    python export_gp_interp_png.py
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.io import write_image
from pathlib import Path
from scipy.stats import gaussian_kde
from scipy.interpolate import interp1d


# =============================================================================
# Grid generation functions
# =============================================================================

def build_time_grid_linspace(pseudotime, numPts=14):
    return np.linspace(0.0, 1.0, numPts)


def build_time_grid_quantile(pseudotime, numPts=14):
    quantiles = np.linspace(0, 1, numPts)
    return np.quantile(pseudotime, quantiles)


def build_time_grid_density(pseudotime, numPts=14):
    pt_sorted = np.sort(pseudotime)
    kde = gaussian_kde(pt_sorted)
    x_fine = np.linspace(0, 1, 1000)
    density = kde(x_fine)
    cdf = np.cumsum(density)
    cdf = cdf / cdf[-1]
    quantiles = np.linspace(0, 1, numPts)
    return np.interp(quantiles, cdf, x_fine)


def build_time_grid_custom(pseudotime, custom_times):
    grid_values = np.asarray(custom_times)
    return grid_values


# =============================================================================
# Data loading
# =============================================================================

def load_test_data(data_path):
    """Load test data from .npz file."""
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
        "cell_metacell": data["cell_metacell"],
        "loop_ids": data["loop_ids"],
        "variance_matrix": data["variance_matrix"],
    }


# =============================================================================
# Visualization functions for PNG export
# =============================================================================

def make_method_comparison_png(data, methods_data, scaling="z_score", dpi=300):
    """
    Create a publication-quality comparison figure for PNG export.
    
    Each row shows: Raw data (left) + Method-specific GP-smoothed (right)
    """
    n_methods = len(methods_data)
    
    if scaling == "z_score":
        raw_key = "raw_z_score_scaled"
        colorscale = "RdBu_r"
        zmin, zmax = -3, 3
        cbar_title = "Z-score"
    else:
        raw_key = "raw_zero_one_scaled"
        colorscale = "inferno"
        zmin, zmax = 0, 1
        cbar_title = "Scaled [0, 1]"
    
    raw_matrix = data[raw_key]
    cell_pseudotime = data["cell_pseudotime"]
    loop_ids = data["loop_ids"]
    original_grid = data["time_grid"]
    smooth_matrix_orig = data["smooth_z_score_scaled"] if scaling == "z_score" else data["smooth_zero_one_scaled"]
    
    # Sort cells by pseudotime
    cell_order = np.argsort(cell_pseudotime)
    raw_sorted = raw_matrix[:, cell_order]
    
    # Calculate figure size based on number of loops and methods
    loop_height = max(8, min(15, 600 / len(loop_ids)))
    row_height = max(80, len(loop_ids) * loop_height)
    fig_height = max(200, row_height * n_methods)
    fig_width = 1400
    
    # Create figure
    fig = make_subplots(
        rows=n_methods, cols=2,
        shared_yaxes=True,
        shared_xaxes=False,
        vertical_spacing=0.04,
        horizontal_spacing=0.06,
        column_widths=[0.55, 0.45],
        subplot_titles=[f"{name} (raw)" for name, _ in methods_data] + 
                       [f"{name} ({len(grid)} pts)" for name, grid in methods_data]
    )
    
    for row_idx, (method_name, grid) in enumerate(methods_data, 1):
        # Raw panel with grid lines
        fig.add_trace(
            go.Heatmap(
                z=raw_sorted,
                x=cell_pseudotime[cell_order],
                y=list(range(len(loop_ids))),
                colorscale=colorscale,
                zmin=zmin,
                zmax=zmax,
                showscale=False,
                hovertemplate=None,
            ),
            row=row_idx, col=1
        )
        
        # Add vertical grid lines
        for t in grid:
            fig.add_vline(x=t, line_dash="dash", line_color="black", line_width=1.5,
                          row=row_idx, col=1)
        
        # Smooth panel - interpolate existing smooth_matrix to new grid
        smooth_interp = np.zeros((smooth_matrix_orig.shape[0], len(grid)))
        for i in range(smooth_matrix_orig.shape[0]):
            f = interp1d(original_grid, smooth_matrix_orig[i], kind='cubic',
                         bounds_error=False, fill_value='extrapolate')
            smooth_interp[i] = f(grid)
        
        fig.add_trace(
            go.Heatmap(
                z=smooth_interp,
                x=grid,
                y=list(range(len(loop_ids))),
                colorscale=colorscale,
                zmin=zmin,
                zmax=zmax,
                showscale=row_idx == n_methods,
                colorbar=dict(
                    title=cbar_title, 
                    x=1.02,
                    thickness=20,
                    len=0.5,
                    title_font_size=14,
                    tickfont_size=12
                ) if row_idx == n_methods else None,
                hovertemplate=None,
            ),
            row=row_idx, col=2
        )
    
    fig.update_layout(
        height=fig_height,
        width=fig_width,
        showlegend=False,
        margin=dict(l=60, r=80, t=80, b=40),
        title_text="GP Interpolation: Grid Method Comparison",
        title_y=0.995,
        title_font_size=20,
        paper_bgcolor='white',
        plot_bgcolor='white',
    )
    
    fig.update_xaxes(title_text="Pseudotime", row=n_methods, col=1, 
                     title_font_size=14, tickfont_size=11)
    fig.update_xaxes(title_text="Grid time", row=n_methods, col=2,
                     title_font_size=14, tickfont_size=11)
    fig.update_yaxes(title_text="Loop index", 
                     title_font_size=14, tickfont_size=11)
    
    # Update subplot title fonts
    fig.update_annotations(font_size=14, font_weight="bold")
    
    return fig


def make_grid_distribution_png(pseudotime, methods_data, dpi=300):
    """
    Create a publication-quality grid distribution figure for PNG export.
    """
    fig = make_subplots(
        rows=len(methods_data), cols=1,
        shared_xaxes=True,
        vertical_spacing=0.06,
        row_heights=[1] * len(methods_data),
        subplot_titles=[f"{name}: {len(grid)} points" for name, grid in methods_data]
    )
    
    # Plot pseudotime histogram as background
    hist, bin_edges = np.histogram(pseudotime, bins=50, range=(0, 1))
    hist_norm = hist / hist.max() * 0.3
    
    fig_height = max(150, 100 * len(methods_data))
    fig_width = 1200
    
    for row_idx, (method_name, grid) in enumerate(methods_data, 1):
        # Add histogram background
        fig.add_trace(
            go.Bar(
                x=bin_edges[:-1],
                y=hist_norm,
                width=bin_edges[1] - bin_edges[0],
                opacity=0.3,
                marker_color="gray",
                hoverinfo="skip",
                showlegend=False,
            ),
            row=row_idx, col=1
        )
        
        # Add grid points as vertical lines
        for t in grid:
            fig.add_vline(x=t, line_dash="solid", line_color="blue", line_width=2.5,
                          row=row_idx, col=1)
        
        # Add scatter points at top for visibility
        fig.add_trace(
            go.Scatter(
                x=grid,
                y=[0.5] * len(grid),
                mode="markers",
                marker=dict(size=15, color="red", symbol="line-ns-open", line_width=2),
                hovertemplate=None,
                showlegend=False,
            ),
            row=row_idx, col=1
        )
    
    fig.update_layout(
        height=fig_height,
        width=fig_width,
        showlegend=False,
        margin=dict(l=40, r=20, t=60, b=40),
        title_text="Grid Point Distribution by Method (background: cell density)",
        title_y=0.98,
        title_font_size=18,
        paper_bgcolor='white',
        plot_bgcolor='white',
    )
    
    fig.update_xaxes(title_text="Normalized pseudotime", row=len(methods_data), col=1,
                     title_font_size=14, tickfont_size=12, range=[0, 1])
    fig.update_yaxes(showticklabels=False, showgrid=False, range=[0, 1.2])
    
    fig.update_annotations(font_size=14, font_weight="bold")
    
    return fig


def make_individual_method_panels(data, methods_data, scaling="z_score"):
    """
    Create individual PNG files for each method (one panel per method).
    """
    figures = {}
    
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
    
    for method_name, grid in methods_data:
        fig = make_subplots(
            rows=1, cols=2,
            shared_yaxes=True,
            horizontal_spacing=0.06,
            column_widths=[0.55, 0.45],
            subplot_titles=["Raw data", f"GP-smoothed ({len(grid)} pts)"]
        )
        
        # Raw panel
        fig.add_trace(
            go.Heatmap(
                z=raw_sorted,
                x=cell_pseudotime[cell_order],
                y=list(range(len(loop_ids))),
                colorscale=colorscale,
                zmin=zmin,
                zmax=zmax,
                showscale=True,
                colorbar=dict(title="Z-score" if scaling == "z_score" else "Scaled", x=0.45),
                hovertemplate=None,
            ),
            row=1, col=1
        )
        
        # Grid lines
        for t in grid:
            fig.add_vline(x=t, line_dash="dash", line_color="black", line_width=1.5,
                          row=1, col=1)
        
        # Smooth panel
        smooth_interp = np.zeros((smooth_matrix_orig.shape[0], len(grid)))
        for i in range(smooth_matrix_orig.shape[0]):
            f = interp1d(original_grid, smooth_matrix_orig[i], kind='cubic',
                         bounds_error=False, fill_value='extrapolate')
            smooth_interp[i] = f(grid)
        
        fig.add_trace(
            go.Heatmap(
                z=smooth_interp,
                x=grid,
                y=list(range(len(loop_ids))),
                colorscale=colorscale,
                zmin=zmin,
                zmax=zmax,
                showscale=True,
                colorbar=dict(title="Z-score" if scaling == "z_score" else "Scaled", x=1.02),
                hovertemplate=None,
            ),
            row=1, col=2
        )
        
        fig.update_layout(
            height=max(400, len(loop_ids) * 10),
            width=1400,
            showlegend=False,
            margin=dict(l=60, r=80, t=60, b=40),
            title_text=f"GP Interpolation: {method_name} method",
            title_y=0.98,
            title_font_size=18,
            paper_bgcolor='white',
            plot_bgcolor='white',
        )
        
        fig.update_xaxes(title_text="Pseudotime", row=1, col=1, title_font_size=14)
        fig.update_xaxes(title_text="Grid time", row=1, col=2, title_font_size=14)
        fig.update_yaxes(title_text="Loop index", title_font_size=14)
        
        figures[method_name] = fig
    
    return figures


# =============================================================================
# Main execution
# =============================================================================

def main():
    # Find test data
    data_path = Path(__file__).parent / "test_data.npz"
    if not data_path.exists():
        print(f"Error: Test data not found at {data_path}")
        return
    
    print("=" * 60)
    print("GP Interpolation: PNG Export")
    print("=" * 60)
    
    # Load data
    print(f"\nLoading test data from {data_path}...")
    data = load_test_data(data_path)
    pseudotime = data["cell_pseudotime"]
    print(f"  Loops: {len(data['loop_ids'])}, Cells: {len(pseudotime)}")
    
    # Define test methods
    numPts = 12
    custom_times = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    
    methods = {
        "linspace": lambda: build_time_grid_linspace(pseudotime, numPts),
        "quantile": lambda: build_time_grid_quantile(pseudotime, numPts),
        "density": lambda: build_time_grid_density(pseudotime, numPts),
        "custom": lambda: build_time_grid_custom(pseudotime, custom_times),
    }
    
    # Generate grids
    print(f"\nGenerating grids (numPts={numPts})...")
    methods_data = []
    for method_name, build_fn in methods.items():
        grid = build_fn()
        methods_data.append((method_name, grid))
        print(f"  {method_name}: {len(grid)} points")
    
    # Create output directory
    output_dir = Path(__file__).parent / "results" / "png"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Check if kaleido is available
    try:
        import kaleido
        print("\n✓ kaleido available for PNG export")
    except ImportError:
        print("\n✗ kaleido not found. Installing...")
        import subprocess
        subprocess.check_call(["pip", "install", "kaleido"])
        import kaleido
    
    # Export comparison figure (Z-score)
    print("\n=== Exporting Comparison Figure (Z-score) ===")
    fig_comparison = make_method_comparison_png(data, methods_data, scaling="z_score")
    output_path = output_dir / "gp_interp_all_methods_comparison.png"
    write_image(fig_comparison, output_path, scale=2, width=1400, height=fig_comparison.layout.height)
    print(f"Saved: {output_path}")
    
    # Export comparison figure (zero_one)
    print("\n=== Exporting Comparison Figure (0-1 scaled) ===")
    fig_zero_one = make_method_comparison_png(data, methods_data, scaling="zero_one")
    output_path = output_dir / "gp_interp_all_methods_zero_one.png"
    write_image(fig_zero_one, output_path, scale=2, width=1400, height=fig_zero_one.layout.height)
    print(f"Saved: {output_path}")
    
    # Export grid distribution
    print("\n=== Exporting Grid Distribution Figure ===")
    fig_distribution = make_grid_distribution_png(pseudotime, methods_data)
    output_path = output_dir / "gp_interp_grid_distribution.png"
    write_image(fig_distribution, output_path, scale=2, width=1200, height=fig_distribution.layout.height)
    print(f"Saved: {output_path}")
    
    # Export individual method panels
    print("\n=== Exporting Individual Method Panels ===")
    individual_figs = make_individual_method_panels(data, methods_data, scaling="z_score")
    for method_name, fig in individual_figs.items():
        output_path = output_dir / f"gp_interp_method_{method_name}.png"
        write_image(fig, output_path, scale=2, width=1400, height=fig.layout.height)
        print(f"Saved: {output_path}")
    
    print("\n" + "=" * 60)
    print("PNG export complete!")
    print(f"Output directory: {output_dir}")
    print("=" * 60)


if __name__ == "__main__":
    main()
