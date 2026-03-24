#!/usr/bin/env python
"""
Standalone visualization script for testing gp_interpolate with different time point selection strategies.

Usage:
    python visualize_gp_interp.py                          # Use default test_data.npz
    python visualize_gp_interp.py --data /path/to/data.npz # Custom test data
    python visualize_gp_interp.py --numPts 12 20 50        # Compare different numPts settings

Supports testing:
    - Default uniform spacing (np.linspace)
    - Quantile-based spacing
    - User-specified custom time grids
"""

import argparse
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pathlib import Path


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


def make_heatmap_figure(data, scaling="z_score", title_prefix="", gp_grid_lines=None):
    """
    Create a two-panel heatmap comparing raw vs GP-smoothed data.
    
    Args:
        data: dict with matrices
        scaling: "z_score" or "zero_one"
        title_prefix: optional prefix for subplot titles
        gp_grid_lines: optional array of GP grid time points for vertical lines
    
    Returns:
        Plotly figure
    """
    if scaling == "z_score":
        raw_key = "raw_z_score_scaled"
        smooth_key = "smooth_z_score_scaled"
        colorscale = "RdBu_r"
        zmin, zmax = -3, 3
        cbar_title = "Z-score"
    else:  # zero_one
        raw_key = "raw_zero_one_scaled"
        smooth_key = "smooth_zero_one_scaled"
        colorscale = "inferno"
        zmin, zmax = 0, 1
        cbar_title = "Scaled [0, 1]"
    
    raw_matrix = data[raw_key]
    smooth_matrix = data[smooth_key]
    cell_pseudotime = data["cell_pseudotime"]
    time_grid = data["time_grid"]
    loop_ids = data["loop_ids"]
    
    # Sort cells by pseudotime
    cell_order = np.argsort(cell_pseudotime)
    raw_sorted = raw_matrix[:, cell_order]
    smooth_sorted = smooth_matrix
    
    # Create figure
    fig = make_subplots(
        rows=1, cols=2,
        shared_yaxes=True,
        horizontal_spacing=0.08,
        column_widths=[0.6, 0.4],
        subplot_titles=[f"{title_prefix}Raw ({len(cell_pseudotime)} cells)", f"{title_prefix}GP-smoothed ({len(time_grid)} grid)"]
    )
    
    # Raw heatmap
    fig.add_trace(
        go.Heatmap(
            z=raw_sorted,
            x=cell_pseudotime[cell_order],
            y=list(range(len(loop_ids))),
            colorscale=colorscale,
            zmin=zmin,
            zmax=zmax,
            colorbar=dict(title=cbar_title, x=0.45),
            hovertemplate="Loop: %{y}<br>Pseudotime: %{x:.3f}<br>Value: %{z:.2f}<extra></extra>",
        ),
        row=1, col=1
    )
    
    # Add GP grid lines
    if gp_grid_lines is not None:
        for t in gp_grid_lines:
            fig.add_vline(x=t, line_dash="dash", line_color="black", line_width=1, row=1, col=1)
    
    # Smooth heatmap
    fig.add_trace(
        go.Heatmap(
            z=smooth_sorted,
            x=time_grid,
            y=list(range(len(loop_ids))),
            colorscale=colorscale,
            zmin=zmin,
            zmax=zmax,
            colorbar=dict(title=cbar_title, x=1.02),
            hovertemplate="Loop: %{y}<br>Grid time: %{x:.3f}<br>Value: %{z:.2f}<extra></extra>",
        ),
        row=1, col=2
    )
    
    fig.update_layout(
        height=max(400, len(loop_ids) * 12),
        showlegend=False,
        margin=dict(l=20, r=20, t=60, b=20),
    )
    
    fig.update_xaxes(title_text="Pseudotime (normalized)", row=1, col=1)
    fig.update_xaxes(title_text="GP grid time", row=1, col=2)
    fig.update_yaxes(title_text="Loop index", row=1, col=1)
    
    return fig


def generate_time_grid(method="uniform", numPts=14, cell_pseudotime=None, custom_times=None):
    """
    Generate time grid with different selection strategies.
    
    Args:
        method: "uniform", "quantile", or "custom"
        numPts: number of grid points
        cell_pseudotime: cell pseudotime values (required for quantile method)
        custom_times: custom time points (required for custom method)
    
    Returns:
        np.array of time points in [0, 1]
    """
    if method == "uniform":
        return np.linspace(0.0, 1.0, numPts)
    
    elif method == "quantile":
        if cell_pseudotime is None:
            raise ValueError("cell_pseudotime required for quantile method")
        quantiles = np.linspace(0.0, 1.0, numPts)
        return np.quantile(cell_pseudotime, quantiles)
    
    elif method == "custom":
        if custom_times is None:
            raise ValueError("custom_times required for custom method")
        return np.asarray(custom_times)
    
    else:
        raise ValueError(f"Unknown method: {method}")


def main():
    parser = argparse.ArgumentParser(description="Visualize GP interpolation results")
    parser.add_argument("--data", type=str, default=None, help="Path to test_data.npz")
    parser.add_argument("--numPts", type=int, nargs="+", default=[12], help="Number of GP grid points to test")
    parser.add_argument("--method", type=str, default="uniform", choices=["uniform", "quantile", "custom"],
                        help="Time point selection method")
    parser.add_argument("--scaling", type=str, default="z_score", choices=["z_score", "zero_one"],
                        help="Heatmap scaling mode")
    parser.add_argument("--output", type=str, default=None, help="Output HTML file path")
    args = parser.parse_args()
    
    # Find test data
    if args.data:
        data_path = Path(args.data)
    else:
        data_path = Path(__file__).parent / "test_data.npz"
    
    if not data_path.exists():
        print(f"Error: Test data not found at {data_path}")
        print("Run extract_test_data.py first or specify --data")
        return
    
    print(f"Loading test data from {data_path}...")
    data = load_test_data(data_path)
    print(f"  Loops: {len(data['loop_ids'])}, Cells: {len(data['cell_pseudotime'])}")
    
    # Create comparison figure for different numPts
    if len(args.numPts) > 1:
        print(f"Comparing numPts settings: {args.numPts}")
        fig = make_subplots(
            rows=len(args.numPts), cols=2,
            shared_yaxes=True,
            shared_xaxes=False,
            vertical_spacing=0.08,
            horizontal_spacing=0.08,
            column_widths=[0.6, 0.4],
            subplot_titles=[f"numPts={n}" for n in args.numPts for suffix in [" (raw)", " (smooth)"]]
        )
        
        for row_idx, numPts in enumerate(args.numPts, 1):
            # Generate grid
            time_grid = generate_time_grid(args.method, numPts=numPts, 
                                           cell_pseudotime=data["cell_pseudotime"])
            
            # For visualization, we just show the existing smooth_matrix
            # In real testing, you would re-run gp_interpolate with new numPts
            scaling_key = f"smooth_{args.scaling}_scaled"
            smooth_matrix = data[scaling_key]
            raw_key = f"raw_{args.scaling}_scaled"
            raw_matrix = data[raw_key]
            
            cell_order = np.argsort(data["cell_pseudotime"])
            raw_sorted = raw_matrix[:, cell_order]
            
            # Raw panel
            fig.add_trace(
                go.Heatmap(
                    z=raw_sorted,
                    x=data["cell_pseudotime"][cell_order],
                    y=list(range(len(data["loop_ids"]))),
                    colorscale="inferno" if args.scaling == "zero_one" else "RdBu_r",
                    zmin=0 if args.scaling == "zero_one" else -3,
                    zmax=1 if args.scaling == "zero_one" else 3,
                    showscale=False,
                    hovertemplate="Loop: %{y}<br>Pseudotime: %{x:.3f}<br>Value: %{z:.2f}<extra></extra>",
                ),
                row=row_idx, col=1
            )
            
            # Add GP grid lines
            for t in time_grid:
                fig.add_vline(x=t, line_dash="dash", line_color="black", line_width=1, 
                              row=row_idx, col=1)
            
            # Smooth panel (downsampled for visualization)
            from scipy.interpolate import interp1d
            smooth_vis = np.zeros((smooth_matrix.shape[0], numPts))
            for i in range(smooth_matrix.shape[0]):
                f = interp1d(data["time_grid"], smooth_matrix[i], kind='cubic', 
                             bounds_error=False, fill_value='extrapolate')
                smooth_vis[i] = f(time_grid)
            
            fig.add_trace(
                go.Heatmap(
                    z=smooth_vis,
                    x=time_grid,
                    y=list(range(len(data["loop_ids"]))),
                    colorscale="inferno" if args.scaling == "zero_one" else "RdBu_r",
                    zmin=0 if args.scaling == "zero_one" else -3,
                    zmax=1 if args.scaling == "zero_one" else 3,
                    showscale=row_idx == len(args.numPts),
                    colorbar=dict(title="Z-score" if args.scaling == "z_score" else "Scaled [0,1]", 
                                  x=1.02) if row_idx == len(args.numPts) else None,
                    hovertemplate="Loop: %{y}<br>Grid time: %{x:.3f}<br>Value: %{z:.2f}<extra></extra>",
                ),
                row=row_idx, col=2
            )
        
        fig.update_layout(
            height=max(300, len(data["loop_ids"]) * 10 * len(args.numPts)),
            showlegend=False,
            margin=dict(l=20, r=30, t=60, b=20),
            title_text=f"GP Interpolation Comparison (method={args.method})",
            title_y=0.98
        )
        
    else:
        # Single numPts - standard two-panel view
        numPts = args.numPts[0]
        time_grid = generate_time_grid(args.method, numPts=numPts,
                                       cell_pseudotime=data["cell_pseudotime"])
        
        print(f"Using {args.method} time grid with {numPts} points: {time_grid}")
        
        fig = make_heatmap_figure(
            data, 
            scaling=args.scaling,
            title_prefix=f"numPts={numPts}, method={args.method} - ",
            gp_grid_lines=time_grid
        )
    
    # Save or show
    if args.output:
        output_path = Path(args.output)
    else:
        output_path = Path(__file__).parent / f"gp_interp_{args.method}_numPts{'_'.join(map(str, args.numPts))}.html"
    
    fig.write_html(output_path)
    print(f"Saved figure to {output_path}")
    
    # Also print grid comparison
    print("\n=== Time Grid Comparison ===")
    print(f"Uniform (numPts={numPts}): {np.linspace(0, 1, numPts)}")
    print(f"Quantile (numPts={numPts}): {np.quantile(data['cell_pseudotime'], np.linspace(0, 1, numPts))}")


if __name__ == "__main__":
    main()
