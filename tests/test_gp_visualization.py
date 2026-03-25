#!/usr/bin/env python
"""
Tests for GP interpolation visualization functions.

Usage:
    python test_gp_visualization.py
"""

import numpy as np
import pandas as pd
import anndata as ad
from pathlib import Path
import sys

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from hic_basic.interpolate import (
    gp_interpolate,
    plot_gp_heatmap,
    _plot_grid_comparison,
    _build_time_grid,
)


def create_test_adata(n_cells=200, n_genes=50, seed=42):
    """Create a simple test AnnData object with clear expression patterns."""
    np.random.seed(seed)

    # Create pseudotime with uniform distribution for clearer visualization
    pseudotime = np.sort(np.random.uniform(0, 1, n_cells))

    # Create expression data with clear pseudotime-dependent patterns
    X = np.zeros((n_cells, n_genes))
    for i in range(n_genes):
        # Different genes have different expression patterns with less noise
        base_pattern = np.zeros(n_cells)
        
        if i % 5 == 0:
            # Increasing pattern (monotonic)
            base_pattern = pseudotime ** 2
            noise_scale = 0.15
        elif i % 5 == 1:
            # Decreasing pattern (monotonic)
            base_pattern = (1 - pseudotime) ** 2
            noise_scale = 0.15
        elif i % 5 == 2:
            # Peak in early pseudotime
            base_pattern = np.exp(-((pseudotime - 0.25) ** 2) / 0.05)
            noise_scale = 0.2
        elif i % 5 == 3:
            # Peak in middle pseudotime
            base_pattern = np.exp(-((pseudotime - 0.5) ** 2) / 0.05)
            noise_scale = 0.2
        else:
            # Peak in late pseudotime
            base_pattern = np.exp(-((pseudotime - 0.75) ** 2) / 0.05)
            noise_scale = 0.2
        
        # Add moderate noise (not too much to blur the pattern)
        X[:, i] = base_pattern + np.random.randn(n_cells) * noise_scale

    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame({"pseudotime": pseudotime}),
        var=pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)]),
    )

    return adata


def test_plot_gp_heatmap():
    """Test the plot_gp_heatmap function."""
    print("=" * 70)
    print("Testing plot_gp_heatmap")
    print("=" * 70)

    # Create test data
    adata = create_test_adata(n_cells=200, n_genes=50)
    print(f"Created test AnnData: {adata.shape[0]} cells, {adata.shape[1]} genes")

    # Test with default settings
    output_dir = Path(__file__).parent / "test_output"
    output_dir.mkdir(exist_ok=True)

    # Test 1: Default settings with half-interval endpoint mode
    print("\nTest 1: Default settings (linspace, half-interval)")
    fig = plot_gp_heatmap(
        adata,
        pseudotime_col="pseudotime",
        numPts=12,
        grid_method="linspace",
        endpoint_mode="half-interval",
        scaling="z_score",
        figsize=(8, 6),
        save_path=output_dir / "test1_gp_heatmap_default.png",
    )
    print(f"  Saved: {output_dir / 'test1_gp_heatmap_default.png'}")

    # Test 2: Different grid method (quantile)
    print("\nTest 2: Quantile grid method")
    fig = plot_gp_heatmap(
        adata,
        pseudotime_col="pseudotime",
        numPts=12,
        grid_method="quantile",
        endpoint_mode="half-interval",
        scaling="z_score",
        figsize=(8, 6),
        save_path=output_dir / "test2_gp_heatmap_quantile.png",
    )
    print(f"  Saved: {output_dir / 'test2_gp_heatmap_quantile.png'}")

    # Test 3: Density grid method
    print("\nTest 3: Density grid method")
    fig = plot_gp_heatmap(
        adata,
        pseudotime_col="pseudotime",
        numPts=12,
        grid_method="density",
        endpoint_mode="half-interval",
        scaling="z_score",
        figsize=(8, 6),
        save_path=output_dir / "test3_gp_heatmap_density.png",
    )
    print(f"  Saved: {output_dir / 'test3_gp_heatmap_density.png'}")

    # Test 4: Zero-one scaling
    print("\nTest 4: Zero-one scaling")
    fig = plot_gp_heatmap(
        adata,
        pseudotime_col="pseudotime",
        numPts=12,
        grid_method="linspace",
        endpoint_mode="half-interval",
        scaling="zero_one",
        figsize=(8, 6),
        save_path=output_dir / "test4_gp_heatmap_zero_one.png",
    )
    print(f"  Saved: {output_dir / 'test4_gp_heatmap_zero_one.png'}")

    # Test 5: Hide grid lines
    print("\nTest 5: Hide grid lines")
    fig = plot_gp_heatmap(
        adata,
        pseudotime_col="pseudotime",
        numPts=12,
        grid_method="linspace",
        endpoint_mode="half-interval",
        scaling="z_score",
        figsize=(8, 6),
        show_grid=False,
        save_path=output_dir / "test5_gp_heatmap_no_grid.png",
    )
    print(f"  Saved: {output_dir / 'test5_gp_heatmap_no_grid.png'}")

    print("\n" + "=" * 70)
    print("All plot_gp_heatmap tests completed!")
    print(f"Output directory: {output_dir}")
    print("=" * 70)


def test_plot_grid_methods_comparison_internal():
    """Test the internal _plot_grid_comparison function."""
    print("\n" + "=" * 70)
    print("Testing _plot_grid_comparison (internal debug function)")
    print("=" * 70)

    # Create test data
    adata = create_test_adata(n_cells=200, n_genes=50)
    print(f"Created test AnnData: {adata.shape[0]} cells, {adata.shape[1]} genes")

    output_dir = Path(__file__).parent / "test_output"
    output_dir.mkdir(exist_ok=True)

    # Test internal comparison function
    print("\nTesting _plot_grid_comparison...")
    fig, _ = _plot_grid_comparison(
        adata,
        pseudotime_col="pseudotime",
        numPts=12,
        grid_methods=["linspace", "quantile", "density"],
        endpoint_mode="half-interval",
        scaling="z_score",
        figsize=(10, 10),
    )
    fig.savefig(output_dir / "debug_grid_comparison.png", dpi=100, bbox_inches="tight")
    print(f"  Saved: {output_dir / 'debug_grid_comparison.png'}")

    print("\n" + "=" * 70)
    print("Internal comparison function test completed!")
    print("=" * 70)


def test_gp_interpolate_with_visualization():
    """Test gp_interpolate and visualize the results."""
    print("\n" + "=" * 70)
    print("Testing gp_interpolate with visualization")
    print("=" * 70)

    # Create test data
    adata = create_test_adata(n_cells=150, n_genes=30)
    print(f"Created test AnnData: {adata.shape[0]} cells, {adata.shape[1]} genes")

    output_dir = Path(__file__).parent / "test_output"
    output_dir.mkdir(exist_ok=True)

    # Run gp_interpolate with different methods
    methods = ["linspace", "quantile", "density"]

    for method in methods:
        print(f"\nRunning gp_interpolate with method={method}...")
        result = gp_interpolate(
            adata,
            pseudotime_col="pseudotime",
            numPts=12,
            grid_method=method,
            endpoint_mode="half-interval",
            inplace=False,
        )

        # Check results
        grid_data = result.uns["gp_interpolate"]
        print(f"  Grid points: {len(grid_data['time_grid'])}")
        print(f"  Mean shape: {grid_data['mean'].shape}")
        print(f"  Endpoint mode: {grid_data['metadata']['endpoint_mode']}")

    # Create heatmap for each method using the new public function
    print("\nCreating heatmaps with plot_gp_heatmap...")
    for method in methods:
        fig = plot_gp_heatmap(
            adata,
            pseudotime_col="pseudotime",
            numPts=12,
            grid_method=method,
            endpoint_mode="half-interval",
            scaling="z_score",
            figsize=(8, 6),
            save_path=output_dir / f"gp_heatmap_{method}.png",
        )
        print(f"  Saved: {output_dir / f'gp_heatmap_{method}.png'}")

    print("\n" + "=" * 70)
    print("gp_interpolate visualization test completed!")
    print("=" * 70)


def main():
    """Run all tests."""
    try:
        test_plot_gp_heatmap()
        test_plot_grid_methods_comparison_internal()
        test_gp_interpolate_with_visualization()

        print("\n" + "=" * 70)
        print("ALL TESTS PASSED!")
        print("=" * 70)

    except Exception as e:
        print(f"\nTEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
