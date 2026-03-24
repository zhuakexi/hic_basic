#!/usr/bin/env python
"""
Test script for endpoint_mode='half-interval' feature.

Usage:
    python test_endpoint_mode.py
"""

import numpy as np
import pandas as pd
from pathlib import Path

# Import from the module
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from hic_basic.interpolate import _build_time_grid, gp_interpolate


def test_build_time_grid_endpoint_modes():
    """Test _build_time_grid with different endpoint_modes."""
    
    # Create test pseudotime data
    np.random.seed(42)
    pseudotime = pd.Series(np.concatenate([
        np.random.beta(2, 5, 100),  # Left dense
        np.random.uniform(0.2, 0.8, 200),  # Middle uniform
        np.random.beta(5, 2, 100),  # Right dense
    ]))
    pseudotime = (pseudotime - pseudotime.min()) / (pseudotime.max() - pseudotime.min())
    
    numPts = 12
    methods = ["linspace", "quantile", "density", "custom"]
    custom_grid = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    
    print("=" * 70)
    print("Testing _build_time_grid with endpoint_mode")
    print("=" * 70)
    print(f"Pseudotime range: [{pseudotime.min():.3f}, {pseudotime.max():.3f}]")
    print(f"numPts: {numPts}")
    print()
    
    for method in methods:
        print(f"=== Method: {method} ===")
        
        # Test endpoint_mode="full" (legacy)
        grid_full, _ = _build_time_grid(
            pseudotime, numPts=numPts, method=method,
            grid_custom=custom_grid if method == "custom" else None,
            endpoint_mode="full"
        )
        
        # Test endpoint_mode="half-interval" (new default)
        grid_half, _ = _build_time_grid(
            pseudotime, numPts=numPts, method=method,
            grid_custom=custom_grid if method == "custom" else None,
            endpoint_mode="half-interval"
        )
        
        print(f"  endpoint_mode='full':     [{grid_full[0]:.4f}, {grid_full[-1]:.4f}]")
        print(f"  endpoint_mode='half-interval': [{grid_half[0]:.4f}, {grid_half[-1]:.4f}]")
        
        # Verify half-interval is inside full
        assert grid_half[0] >= grid_full[0], f"{method}: half-interval start should be >= full start"
        assert grid_half[-1] <= grid_full[-1], f"{method}: half-interval end should be <= full end"
        
        # Calculate expected half-interval
        interval = 1.0 / (numPts - 1)
        expected_half_interval = interval / 2
        print(f"  Expected half-interval: {expected_half_interval:.4f}")
        print(f"  Actual half-interval:   {grid_half[0]:.4f}")
        print()
    
    print("=" * 70)
    print("All tests passed!")
    print("=" * 70)


def test_gp_interpolate_endpoint_mode():
    """Test gp_interpolate with different endpoint_modes."""
    import anndata as ad
    
    print("\n" + "=" * 70)
    print("Testing gp_interpolate with endpoint_mode")
    print("=" * 70)
    
    # Create minimal test data
    np.random.seed(42)
    n_cells = 100
    n_genes = 10
    
    X = np.random.randn(n_cells, n_genes)
    pseudotime = np.linspace(0, 1, n_cells)
    
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame({"pseudotime": pseudotime}),
        var=pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)])
    )
    
    # Test with endpoint_mode="full"
    result_full = gp_interpolate(
        adata, pseudotime_col="pseudotime",
        numPts=8, grid_method="linspace",
        endpoint_mode="full",
        inplace=False
    )
    
    # Test with endpoint_mode="half-interval"
    result_half = gp_interpolate(
        adata, pseudotime_col="pseudotime",
        numPts=8, grid_method="linspace",
        endpoint_mode="half-interval",
        inplace=False
    )
    
    grid_full = result_full.uns["gp_interpolate"]["time_grid"].values
    grid_half = result_half.uns["gp_interpolate"]["time_grid"].values
    
    print(f"endpoint_mode='full':     [{grid_full[0]:.4f}, {grid_full[-1]:.4f}]")
    print(f"endpoint_mode='half-interval': [{grid_half[0]:.4f}, {grid_half[-1]:.4f}]")
    
    # Verify metadata includes endpoint_mode
    assert result_half.uns["gp_interpolate"]["metadata"]["endpoint_mode"] == "half-interval"
    print(f"Metadata includes endpoint_mode: ✓")
    
    print("=" * 70)
    print("gp_interpolate tests passed!")
    print("=" * 70)


if __name__ == "__main__":
    test_build_time_grid_endpoint_modes()
    test_gp_interpolate_endpoint_mode()
    print("\n" + "=" * 70)
    print("ALL TESTS COMPLETED SUCCESSFULLY!")
    print("=" * 70)
