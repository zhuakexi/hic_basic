#!/usr/bin/env python
"""
Comprehensive test of all gp_interpolate grid methods with heatmap visualization.

This script:
1. Loads test data from test_data.npz
2. Tests all 5 grid methods: linspace, quantile, density, optimal, custom
3. Generates comparison heatmaps for visual inspection
4. Saves results to HTML files

Usage:
    python test_all_grid_methods.py
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pathlib import Path
from scipy.stats import gaussian_kde


# =============================================================================
# Grid generation functions (mirroring hic_basic/interpolate.py _build_time_grid)
# =============================================================================

def build_time_grid_linspace(pseudotime, numPts=14):
    """Evenly spaced grid points."""
    return np.linspace(0.0, 1.0, numPts)


def build_time_grid_quantile(pseudotime, numPts=14):
    """Quantile-based spacing (equal cells per bin)."""
    quantiles = np.linspace(0, 1, numPts)
    return np.quantile(pseudotime, quantiles)


def build_time_grid_density(pseudotime, numPts=14):
    """Density-adaptive spacing (more points in dense regions)."""
    pt_sorted = np.sort(pseudotime)
    kde = gaussian_kde(pt_sorted)
    x_fine = np.linspace(0, 1, 1000)
    density = kde(x_fine)
    cdf = np.cumsum(density)
    cdf = cdf / cdf[-1]
    quantiles = np.linspace(0, 1, numPts)
    return np.interp(quantiles, cdf, x_fine)


def build_time_grid_optimal(pseudotime, numPts=14):
    """Optimal binning (requires optbinning package)."""
    try:
        from optbinning import ContinuousOptimalBinning
    except ImportError:
        raise ImportError("optbinning required. Install with: pip install optbinning")
    
    binning = ContinuousOptimalBinning(n_bins=20, min_bin_size=0.05, monotonic_trend=None)
    binning.fit(pseudotime)
    grid_values = binning.splits
    if grid_values[0] > 0:
        grid_values = np.concatenate([[0.0], grid_values])
    if grid_values[-1] < 1:
        grid_values = np.concatenate([grid_values, [1.0]])
    return grid_values


def build_time_grid_custom(pseudotime, custom_times):
    """User-provided custom grid."""
    grid_values = np.asarray(custom_times)
    if not np.all((grid_values >= 0) & (grid_values <= 1)):
        raise ValueError("Custom grid values must be in range [0, 1].")
    if not np.all(np.diff(grid_values) > 0):
        raise ValueError("Custom grid values must be strictly increasing.")
    return grid_values


# =============================================================================
# Visualization functions
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


def make_method_comparison_figure(data, methods_data, scaling="z_score"):
    """
    Create a multi-row figure comparing all grid methods.
    
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
    
    # Sort cells by pseudotime
    cell_order = np.argsort(cell_pseudotime)
    raw_sorted = raw_matrix[:, cell_order]
    
    # Create figure
    fig = make_subplots(
        rows=n_methods, cols=2,
        shared_yaxes=True,
        shared_xaxes=False,
        vertical_spacing=0.06,
        horizontal_spacing=0.08,
        column_widths=[0.55, 0.45],
        subplot_titles=[f"{name} (raw)" for name, _ in methods_data] + 
                       [f"{name} (grid: {len(grid)} pts)" for name, grid in methods_data]
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
                hovertemplate=f"Method: {method_name}<br>Loop: %{{y}}<br>Pseudotime: %{{x:.3f}}<br>Value: %{{z:.2f}}<extra></extra>",
            ),
            row=row_idx, col=1
        )
        
        # Add vertical grid lines
        for t in grid:
            fig.add_vline(x=t, line_dash="dash", line_color="black", line_width=1,
                          row=row_idx, col=1)
        
        # Smooth panel - interpolate existing smooth_matrix to new grid
        from scipy.interpolate import interp1d
        original_grid = data["time_grid"]
        smooth_matrix = data["smooth_z_score_scaled"] if scaling == "z_score" else data["smooth_zero_one_scaled"]
        
        smooth_interp = np.zeros((smooth_matrix.shape[0], len(grid)))
        for i in range(smooth_matrix.shape[0]):
            f = interp1d(original_grid, smooth_matrix[i], kind='cubic',
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
                colorbar=dict(title=cbar_title, x=1.02) if row_idx == n_methods else None,
                hovertemplate=f"Method: {method_name}<br>Loop: %{{y}}<br>Grid time: %{{x:.3f}}<br>Value: %{{z:.2f}}<extra></extra>",
            ),
            row=row_idx, col=2
        )
    
    fig.update_layout(
        height=max(200, len(loop_ids) * 8 * n_methods),
        showlegend=False,
        margin=dict(l=20, r=40, t=80, b=20),
        title_text="GP Interpolation: Grid Method Comparison",
        title_y=0.99
    )
    
    fig.update_xaxes(title_text="Pseudotime", row=n_methods, col=1)
    fig.update_xaxes(title_text="Grid time", row=n_methods, col=2)
    fig.update_yaxes(title_text="Loop index")
    
    return fig


def make_grid_distribution_figure(pseudotime, methods_data):
    """
    Create a figure showing the grid point distributions for each method.
    """
    fig = make_subplots(
        rows=len(methods_data), cols=1,
        shared_xaxes=True,
        vertical_spacing=0.08,
        subplot_titles=[f"{name}: {len(grid)} points" for name, grid in methods_data]
    )
    
    # Plot pseudotime histogram as background
    hist, bin_edges = np.histogram(pseudotime, bins=50, range=(0, 1))
    hist_norm = hist / hist.max() * 0.3  # Scale for background
    
    for row_idx, (method_name, grid) in enumerate(methods_data, 1):
        # Add histogram background
        fig.add_trace(
            go.Bar(
                x=bin_edges[:-1],
                y=hist_norm,
                width=bin_edges[1] - bin_edges[0],
                opacity=0.2,
                marker_color="gray",
                hoverinfo="skip",
                showlegend=False,
            ),
            row=row_idx, col=1
        )
        
        # Add grid points as vertical lines
        for t in grid:
            fig.add_vline(x=t, line_dash="solid", line_color="blue", line_width=2,
                          row=row_idx, col=1)
        
        # Add scatter points at top for visibility
        fig.add_trace(
            go.Scatter(
                x=grid,
                y=[0.5] * len(grid),
                mode="markers",
                marker=dict(size=12, color="red", symbol="line-ns-open"),
                hovertemplate=f"{method_name}: %{{x:.3f}}<extra></extra>",
                showlegend=False,
            ),
            row=row_idx, col=1
        )
    
    fig.update_layout(
        height=max(150, 80 * len(methods_data)),
        showlegend=False,
        margin=dict(l=20, r=20, t=60, b=20),
        title_text="Grid Point Distribution by Method (background: cell density)",
        title_y=0.98
    )
    
    fig.update_xaxes(title_text="Normalized pseudotime", row=len(methods_data), col=1)
    fig.update_yaxes(showticklabels=False, showgrid=False)
    
    return fig


def make_summary_statistics_table(methods_data, pseudotime):
    """
    Create a summary statistics table for all methods.
    """
    stats = []
    for method_name, grid in methods_data:
        spacings = np.diff(grid)
        stats.append({
            "Method": method_name,
            "n_points": len(grid),
            "min_spacing": float(np.min(spacings)),
            "max_spacing": float(np.max(spacings)),
            "mean_spacing": float(np.mean(spacings)),
            "std_spacing": float(np.std(spacings)),
            "uniformity": float(1.0 - np.std(spacings) / np.mean(spacings)) if np.mean(spacings) > 0 else 0,
        })
    
    df = pd.DataFrame(stats)
    return df


# =============================================================================
# Main execution
# =============================================================================

def main():
    # Find test data
    data_path = Path(__file__).parent / "test_data.npz"
    if not data_path.exists():
        print(f"Error: Test data not found at {data_path}")
        print("Run extract_test_data.py first")
        return
    
    print("=" * 60)
    print("GP Interpolation: Grid Method Comparison Test")
    print("=" * 60)
    
    # Load data
    print(f"\nLoading test data from {data_path}...")
    data = load_test_data(data_path)
    pseudotime = data["cell_pseudotime"]
    print(f"  Loops: {len(data['loop_ids'])}, Cells: {len(pseudotime)}")
    print(f"  Pseudotime range: [{pseudotime.min():.3f}, {pseudotime.max():.3f}]")
    
    # Define test methods
    numPts = 12
    custom_times = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    
    methods = {
        "linspace": lambda: build_time_grid_linspace(pseudotime, numPts),
        "quantile": lambda: build_time_grid_quantile(pseudotime, numPts),
        "density": lambda: build_time_grid_density(pseudotime, numPts),
        "custom": lambda: build_time_grid_custom(pseudotime, custom_times),
    }
    
    # Try optimal method (may not have optbinning)
    try:
        methods["optimal"] = lambda: build_time_grid_optimal(pseudotime, numPts)
    except ImportError:
        print("  Note: optbinning not installed, skipping 'optimal' method")
    
    # Generate grids
    print(f"\nGenerating grids (numPts={numPts})...")
    methods_data = []
    for method_name, build_fn in methods.items():
        try:
            grid = build_fn()
            methods_data.append((method_name, grid))
            print(f"  {method_name}: {len(grid)} points, range [{grid.min():.3f}, {grid.max():.3f}]")
        except Exception as e:
            print(f"  {method_name}: FAILED - {e}")
    
    # Summary statistics
    print("\n=== Summary Statistics ===")
    stats_df = make_summary_statistics_table(methods_data, pseudotime)
    print(stats_df.to_string(index=False))
    
    # Create output directory
    output_dir = Path(__file__).parent / "results"
    output_dir.mkdir(exist_ok=True)
    
    # Generate comparison heatmap
    print("\n=== Generating Comparison Heatmap ===")
    fig_comparison = make_method_comparison_figure(data, methods_data, scaling="z_score")
    output_comparison = output_dir / "gp_interp_all_methods_comparison.html"
    fig_comparison.write_html(output_comparison)
    print(f"Saved comparison heatmap to {output_comparison}")
    
    # Generate grid distribution figure
    print("\n=== Generating Grid Distribution Figure ===")
    fig_distribution = make_grid_distribution_figure(pseudotime, methods_data)
    output_distribution = output_dir / "gp_interp_grid_distribution.html"
    fig_distribution.write_html(output_distribution)
    print(f"Saved grid distribution to {output_distribution}")
    
    # Generate zero_one scaled version
    print("\n=== Generating Zero-One Scaled Heatmap ===")
    fig_zero_one = make_method_comparison_figure(data, methods_data, scaling="zero_one")
    output_zero_one = output_dir / "gp_interp_all_methods_zero_one.html"
    fig_zero_one.write_html(output_zero_one)
    print(f"Saved zero-one scaled heatmap to {output_zero_one}")
    
    # Save summary statistics
    output_stats = output_dir / "gp_interp_summary_stats.csv"
    stats_df.to_csv(output_stats, index=False)
    print(f"Saved summary statistics to {output_stats}")
    
    print("\n" + "=" * 60)
    print("All tests completed!")
    print(f"Results saved to: {output_dir}")
    print("=" * 60)


if __name__ == "__main__":
    main()
