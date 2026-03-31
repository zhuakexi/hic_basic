#!/usr/bin/env python
"""
Create synthetic cooler and eigenvector files for testing plot_tiling_compartment.

This script generates:
1. A synthetic cooler file with random contact matrices
2. Synthetic eigenvector bedgraph files

Usage:
    python tests/data/create_synthetic_cool.py

Output:
    tests/data/synthetic_test.cool
    tests/data/synthetic_eigs1.bedgraph
    tests/data/synthetic_eigs2.bedgraph
"""
import os
import tempfile
import numpy as np
import pandas as pd
import cooler
from pathlib import Path

# Set random seed for reproducibility
np.random.seed(42)

# Output directory
OUTDIR = Path(__file__).parent
COOL_FILE = OUTDIR / "synthetic_test.cool"
EIGS1_FILE = OUTDIR / "synthetic_eigs1.bedgraph"
EIGS2_FILE = OUTDIR / "synthetic_eigs2.bedgraph"

# Chromosome sizes (mm10 style, simplified)
CHROM_SIZES = {
    "chr1": 10_000_000,
    "chr2": 8_000_000,
}

# Bin size
BINSIZE = 100_000  # 100kb


def create_cool_file():
    """Create a synthetic cooler file with random contact matrices."""
    print(f"Creating synthetic cooler file: {COOL_FILE}")
    
    # Create bins
    bins = []
    bin_id_offset = 0
    chrom_offsets = {}  # map chrom name to bin_id offset
    
    for chrom, size in CHROM_SIZES.items():
        n_bins = size // BINSIZE
        chrom_offsets[chrom] = bin_id_offset
        for i in range(n_bins):
            bins.append({
                "chrom": chrom,
                "start": i * BINSIZE,
                "end": (i + 1) * BINSIZE
            })
        bin_id_offset += n_bins
    
    bins_df = pd.DataFrame(bins)
    
    # Create contacts (cis interactions only for simplicity)
    contacts = []
    for chrom, size in CHROM_SIZES.items():
        n_bins = size // BINSIZE
        offset = chrom_offsets[chrom]
        # Generate a distance-decay contact matrix (upper triangle only)
        for i in range(n_bins):
            for j in range(i + 1, n_bins):  # start from i+1 to skip diagonal and avoid duplicates
                # Distance-dependent contact frequency
                dist = j - i
                # Exponential decay with some noise
                expected = 1000 * np.exp(-0.05 * dist) + np.random.poisson(10)
                contacts.append({
                    "bin1_id": offset + i,
                    "bin2_id": offset + j,
                    "count": max(1, int(expected))  # ensure count >= 1
                })
    
    contacts_df = pd.DataFrame(contacts)
    
    # Remove existing file if present
    if COOL_FILE.exists():
        COOL_FILE.unlink()
    
    # Create cooler
    cooler.create_cooler(
        str(COOL_FILE),
        bins=bins_df,
        pixels=contacts_df
    )
    print(f"  Created {len(bins_df)} bins, {len(contacts_df)} contacts")


def create_eigenvector_files():
    """Create synthetic eigenvector bedgraph files with proper column names."""
    print(f"Creating synthetic eigenvector files...")
    
    for eigs_file, offset in [(EIGS1_FILE, 0.3), (EIGS2_FILE, -0.2)]:
        # Create eigenvector data with compartment-like pattern
        data = []
        for chrom, size in CHROM_SIZES.items():
            n_bins = size // BINSIZE
            for i in range(n_bins):
                # Alternating A/B compartments
                compartment_phase = np.sin(2 * np.pi * i / 10) + offset
                # Add some noise
                eig_value = compartment_phase + np.random.normal(0, 0.2)
                data.append({
                    "chrom": chrom,
                    "start": i * BINSIZE,
                    "end": (i + 1) * BINSIZE,
                    "E1": eig_value
                })
        
        eigs_df = pd.DataFrame(data)
        # Write with header for proper column names
        eigs_df.to_csv(eigs_file, sep="\t", index=False, header=True)
        print(f"  Created {eigs_file} with {len(eigs_df)} bins")


if __name__ == "__main__":
    create_cool_file()
    create_eigenvector_files()
    print("\nDone! Test files created:")
    print(f"  - {COOL_FILE}")
    print(f"  - {EIGS1_FILE}")
    print(f"  - {EIGS2_FILE}")
