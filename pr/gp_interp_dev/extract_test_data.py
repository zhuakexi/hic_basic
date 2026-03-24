#!/usr/bin/env python
"""
Extract a small-scale test dataset from series 1 for gp_interpolate development.

Samples 50 loops (stratified by peak_grid_index) while keeping all 740 cells.
Output: test_data.npz (< 1 MB) for rapid iteration.
"""

import numpy as np
import pandas as pd
from pathlib import Path

# Paths from series 1
SERIES1_DIR = Path("/share/home/ychi/dev/sc-embryo/notebooks/E_loop_analysis/6_time_loops/1.gp_interpolated_time_loop_heatmap")
ARRAYS_PATH = SERIES1_DIR / "1.gp_interpolated_time_loop_heatmap.arrays.npz"
LOOP_TABLE_PATH = SERIES1_DIR / "1.gp_interpolated_time_loop_heatmap.loop_table.tsv"

# Output
OUTPUT_DIR = Path(__file__).parent
OUTPUT_PATH = OUTPUT_DIR / "test_data.npz"

N_TEST_LOOPS = 50
RANDOM_STATE = 42


def main():
    print(f"Loading series 1 data from {SERIES1_DIR}...")
    
    # Load cached arrays
    data = np.load(ARRAYS_PATH, allow_pickle=True)
    loop_ids = [str(x) for x in data["loop_ids"]]
    
    # Load loop table for stratified sampling
    loop_table = pd.read_csv(LOOP_TABLE_PATH, sep="\t")
    loop_table["loop_id"] = loop_table["loop_id"].astype(str)
    print(f"  Total loops: {len(loop_table)}")
    
    # Stratified sampling by peak_grid_index (0-11)
    n_per_stratum = N_TEST_LOOPS // 12  # ~4 loops per grid index
    
    sampled_loops = []
    for grid_idx in range(12):
        stratum = loop_table[loop_table["peak_grid_index"] == grid_idx]
        n_sample = min(n_per_stratum, len(stratum))
        if n_sample > 0:
            sampled = stratum.sample(n=n_sample, random_state=RANDOM_STATE + grid_idx)
            sampled_loops.append(sampled)
    
    # Fill remaining slots if needed
    remaining = N_TEST_LOOPS - sum(len(df) for df in sampled_loops)
    if remaining > 0:
        already_sampled = pd.concat(sampled_loops)["loop_id"].tolist()
        remainder_pool = loop_table[~loop_table["loop_id"].isin(already_sampled)]
        extra = remainder_pool.sample(n=remaining, random_state=RANDOM_STATE)
        sampled_loops.append(extra)
    
    sampled_loop_table = pd.concat(sampled_loops, ignore_index=True)
    test_loop_ids = sampled_loop_table["loop_id"].tolist()
    print(f"  Sampled loops: {len(test_loop_ids)} (stratified by peak_grid_index)")
    
    # Extract indices for sampled loops
    loop_id_to_idx = {lid: idx for idx, lid in enumerate(loop_ids)}
    test_indices = [loop_id_to_idx[lid] for lid in test_loop_ids]
    
    # Extract test arrays
    print("Extracting test arrays...")
    test_arrays = {
        # Core matrices
        "raw_matrix": data["raw_matrix"][test_indices],           # (50, 740)
        "smooth_matrix": data["smooth_matrix"][test_indices],     # (50, 12)
        "raw_clipped": data["raw_clipped"][test_indices],
        "smooth_clipped": data["smooth_clipped"][test_indices],
        "raw_zero_one_scaled": data["raw_zero_one_scaled"][test_indices],
        "smooth_zero_one_scaled": data["smooth_zero_one_scaled"][test_indices],
        "raw_z_score_scaled": data["raw_z_score_scaled"][test_indices],
        "smooth_z_score_scaled": data["smooth_z_score_scaled"][test_indices],
        
        # Metadata
        "variance_matrix": data["variance_matrix"][test_indices],
        "time_grid": data["time_grid"],                           # (12,)
        "cell_pseudotime": data["cell_pseudotime"],               # (740,)
        "cell_metacell": data["cell_metacell"],                   # (740,)
        "loop_ids": np.array(test_loop_ids),                      # (50,)
        "gp_length_scale": data["gp_length_scale"][test_indices],
    }
    
    # Save
    np.savez_compressed(OUTPUT_PATH, **test_arrays)
    print(f"Saved test data to {OUTPUT_PATH}")
    
    # Print summary
    print("\n=== Test Data Summary ===")
    print(f"Loops: {len(test_loop_ids)}")
    print(f"Cells: {len(data['cell_pseudotime'])}")
    print(f"GP grid points: {len(data['time_grid'])}")
    print(f"\nArrays saved:")
    for key, arr in test_arrays.items():
        if hasattr(arr, 'shape'):
            print(f"  {key}: {arr.shape} {arr.dtype}")
        else:
            print(f"  {key}: {type(arr)}")
    
    # Save loop table for reference
    sampled_loop_table.to_csv(OUTPUT_DIR / "test_loop_table.tsv", sep="\t", index=False)
    print(f"\nSaved loop table to {OUTPUT_DIR / 'test_loop_table.tsv'}")


if __name__ == "__main__":
    main()
