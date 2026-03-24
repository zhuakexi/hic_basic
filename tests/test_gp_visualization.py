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
    plot_grid_methods_comparison,
    _build_time_grid,
)


def create_test_adata(n_cells=200, n_genes=50, seed=42):
    """Create a simple test AnnData object."""
    np.random.seed(seed)

    # Create pseudotime with some structure
    pseudotime = np.sort(np.random.beta(2, 2, n_cells))
    pseudotime = (pseudotime - pseudotime.min()) / (pseudotime.max() - pseudotime.min())

    # Create expression data with pseudotime-dependent patterns
    X = np.zeros((n_cells, n_genes))
    for i in range(n_genes):
        # Different genes have different expression patterns
        if i % 3 == 0:
            # Increasing pattern
            X[:, i] = pseudotime * np.random.randn(n_cells) * 0.1 + np.random.randn(n_cells) * 0.3
        elif i % 3 == 1:
            # Decreasing pattern
            X[:, i] = (1 - pseudotime) * np.random.randn(n_cells) * 0.1 + np.random.randn(n_cells) * 0.3
        else:
            # Peak in middle
            X[:, i] = np.exp(-((pseudotime - 0.5) ** 2) / 0.1) * np.random.randn(n_cells) * 0.2 + np.random.randn(n_cells) * 0.3

    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame({"pseudotime": pseudotime}),
        var=pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)]),
    )

    return adata


def test_plot_grid_methods_comparison():
    """Test the plot_grid_methods_comparison function."""
    print("=" * 70)
    print("Testing plot_grid_methods_comparison")
    print("=" * 70)

    # Create test data
    adata = create_test_adata(n_cells=200, n_genes=50)
    print(f"Created test AnnData: {adata.shape[0]} cells, {adata.shape[1]} genes")

    # Test with default settings
    output_dir = Path(__file__).parent / "test_output"
    output_dir.mkdir(exist_ok=True)

    # Test 1: Default methods with half-interval endpoint mode
    print("\nTest 1: Default methods (linspace, quantile, density) with half-interval")
    fig = plot_grid_methods_comparison(
        adata,
        pseudotime_col="pseudotime",
        numPts=12,
        grid_methods=["linspace", "quantile", "density"],
        endpoint_mode="half-interval",
        scaling="z_score",
        figsize=(10, 10),
        save_path=output_dir / "test1_default_half_interval.png",
    )
    print(f"  Saved: {output_dir / 'test1_default_half_interval.png'}")

    # Test 2: Default methods with full endpoint mode
    print("\nTest 2: Default methods with full endpoint mode")
    fig = plot_grid_methods_comparison(
        adata,
        pseudotime_col="pseudotime",
        numPts=12,
        grid_methods=["linspace", "quantile", "density"],
        endpoint_mode="full",
        scaling="z_score",
        figsize=(10, 10),
        save_path=output_dir / "test2_default_full.png",
    )
    print(f"  Saved: {output_dir / 'test2_default_full.png'}")

    # Test 3: Single method
    print("\nTest 3: Single method (linspace only)")
    fig = plot_grid_methods_comparison(
        adata,
        pseudotime_col="pseudotime",
        numPts=12,
        grid_methods=["linspace"],
        endpoint_mode="half-interval",
        scaling="z_score",
        figsize=(8, 4),
        save_path=output_dir / "test3_single_method.png",
    )
    print(f"  Saved: {output_dir / 'test3_single_method.png'}")

    # Test 4: Zero-one scaling
    print("\nTest 4: Zero-one scaling")
    fig = plot_grid_methods_comparison(
        adata,
        pseudotime_col="pseudotime",
        numPts=12,
        grid_methods=["linspace", "quantile"],
        endpoint_mode="half-interval",
        scaling="zero_one",
        figsize=(10, 7),
        save_path=output_dir / "test4_zero_one_scaling.png",
    )
    print(f"  Saved: {output_dir / 'test4_zero_one_scaling.png'}")

    # Test 5: All methods including custom
    print("\nTest 5: All methods including custom")
    fig = plot_grid_methods_comparison(
        adata,
        pseudotime_col="pseudotime",
        numPts=10,
        grid_methods=["linspace", "quantile", "density"],
        endpoint_mode="half-interval",
        scaling="z_score",
        figsize=(10, 12),
        save_path=output_dir / "test5_all_methods.png",
    )
    print(f"  Saved: {output_dir / 'test5_all_methods.png'}")

    print("\n" + "=" * 70)
    print("All visualization tests completed!")
    print(f"Output directory: {output_dir}")
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

    # Create comparison plot
    print("\nCreating comparison plot...")
    fig = plot_grid_methods_comparison(
        adata,
        pseudotime_col="pseudotime",
        numPts=12,
        grid_methods=methods,
        endpoint_mode="half-interval",
        scaling="z_score",
        figsize=(10, 10),
        save_path=output_dir / "gp_interpolate_comparison.png",
    )
    print(f"  Saved: {output_dir / 'gp_interpolate_comparison.png'}")

    print("\n" + "=" * 70)
    print("gp_interpolate visualization test completed!")
    print("=" * 70)


def main():
    """Run all tests."""
    try:
        test_plot_grid_methods_comparison()
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
