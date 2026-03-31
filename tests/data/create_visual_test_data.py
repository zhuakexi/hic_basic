#!/usr/bin/env python
"""
Create comprehensive test data for trapezoid Hi-C with multi-track visualization.

This script generates:
1. A synthetic cooler file with compartment-like contact patterns
2. Synthetic eigenvector bedgraph files with A/B compartment signal
3. Synthetic bigWig files for ChIP-seq-like signal
4. Synthetic GTF file with gene annotations

Usage:
    python tests/data/create_visual_test_data.py

Output:
    tests/data/visual_test.cool
    tests/data/visual_test_eigs.bedgraph
    tests/data/visual_test_signal1.bw
    tests/data/visual_test_signal2.bw
    tests/data/visual_test.gtf
"""
import os
import tempfile
import numpy as np
import pandas as pd
import cooler
import pyBigWig
from pathlib import Path

# Set random seed for reproducibility
np.random.seed(42)

# Output directory
OUTDIR = Path(__file__).parent
COOL_FILE = OUTDIR / "visual_test.cool"
EIGS_FILE = OUTDIR / "visual_test_eigs.bedgraph"
BW_FILE1 = OUTDIR / "visual_test_signal1.bw"
BW_FILE2 = OUTDIR / "visual_test_signal2.bw"
GTF_FILE = OUTDIR / "visual_test.gtf"

# Chromosome sizes (simplified mm10 style)
CHROM_SIZES = {
    "chr1": 5_000_000,  # 5Mb for testing
}

# Bin size
BINSIZE = 50_000  # 50kb


def create_cool_file():
    """Create a synthetic cooler file with compartment-like patterns."""
    print(f"Creating synthetic cooler file: {COOL_FILE}")

    # Create bins
    bins = []
    for chrom, size in CHROM_SIZES.items():
        n_bins = size // BINSIZE
        for i in range(n_bins):
            bins.append({
                "chrom": chrom,
                "start": i * BINSIZE,
                "end": (i + 1) * BINSIZE
            })

    bins_df = pd.DataFrame(bins)
    n_bins = len(bins_df)

    # Create contacts with compartment-like pattern
    # Divide into A/B compartments (alternating 500kb blocks)
    compartment_size = 10  # 10 bins = 500kb
    contacts = []
    
    for i in range(n_bins):
        # Determine compartment for each bin
        comp_i = (i // compartment_size) % 2  # 0 for A, 1 for B
        
        for j in range(i, n_bins):
            comp_j = (j // compartment_size) % 2
            
            # Distance decay
            dist = j - i
            base_count = 1000 * np.exp(-0.03 * dist)
            
            # Compartment effect: same compartment has higher contact frequency
            if comp_i == comp_j:
                compartment_factor = 2.0  # Same compartment: 2x contacts
            else:
                compartment_factor = 0.5  # Different compartment: 0.5x contacts
            
            # Add noise
            count = max(1, int(base_count * compartment_factor + np.random.poisson(50)))
            
            contacts.append({
                "bin1_id": i,
                "bin2_id": j,
                "count": count
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


def create_eigenvector_file():
    """Create synthetic eigenvector bedgraph with A/B compartment signal."""
    print(f"Creating eigenvector file: {EIGS_FILE}")

    data = []
    for chrom, size in CHROM_SIZES.items():
        n_bins = size // BINSIZE
        compartment_size = 10  # 500kb compartments
        
        for i in range(n_bins):
            # Alternating A/B compartments
            comp = (i // compartment_size) % 2
            # E1 positive for A compartment, negative for B
            eig_value = 0.5 if comp == 0 else -0.5
            # Add some noise
            eig_value += np.random.normal(0, 0.1)
            
            data.append({
                "chrom": chrom,
                "start": i * BINSIZE,
                "end": (i + 1) * BINSIZE,
                "E1": eig_value
            })

    eigs_df = pd.DataFrame(data)
    eigs_df.to_csv(EIGS_FILE, sep="\t", index=False, header=True)
    print(f"  Created {EIGS_FILE} with {len(eigs_df)} bins")


def create_bigwig_files():
    """Create synthetic bigWig files with ChIP-seq-like signal."""
    print(f"Creating bigWig files...")

    # Signal 1: H3K27ac-like (active enhancers, peaks at compartment boundaries)
    print(f"  Creating {BW_FILE1}")
    with pyBigWig.open(str(BW_FILE1), "w") as bw:
        bw.addHeader(list(CHROM_SIZES.items()))
        
        for chrom, size in CHROM_SIZES.items():
            n_bins = size // BINSIZE
            
            starts = []
            values = []
            for i in range(n_bins):
                # Peak signal at compartment boundaries
                comp_boundary = (i % 10) == 0
                if comp_boundary:
                    signal = 50 + np.random.exponential(30)
                else:
                    signal = 5 + np.random.exponential(5)
                
                starts.append(i * BINSIZE)
                values.append(signal)
            
            # Use span mode: specify chrom, start positions, values, span, and step
            bw.addEntries(chrom, starts, values=values, span=BINSIZE, step=BINSIZE)

    # Signal 2: H3K4me3-like (promoter-associated, random peaks)
    print(f"  Creating {BW_FILE2}")
    with pyBigWig.open(str(BW_FILE2), "w") as bw:
        bw.addHeader(list(CHROM_SIZES.items()))
        
        for chrom, size in CHROM_SIZES.items():
            n_bins = size // BINSIZE
            
            # Place random peaks
            peak_positions = set(np.random.choice(n_bins, size=15, replace=False))
            
            starts = []
            values = []
            for i in range(n_bins):
                if i in peak_positions:
                    signal = 80 + np.random.exponential(40)
                else:
                    signal = 2 + np.random.exponential(3)
                
                starts.append(i * BINSIZE)
                values.append(signal)
            
            bw.addEntries(chrom, starts, values=values, span=BINSIZE, step=BINSIZE)

    print(f"  Created {BW_FILE1} and {BW_FILE2}")


def create_gtf_file():
    """Create synthetic GTF file with gene annotations."""
    print(f"Creating GTF file: {GTF_FILE}")

    genes = []
    
    # Create 8 synthetic genes with varying structures
    gene_configs = [
        # (start, end, strand, n_exons, gene_name, gene_type)
        (100000, 150000, '+', 3, 'Gene1', 'protein_coding'),
        (300000, 450000, '-', 5, 'Gene2', 'protein_coding'),
        (600000, 700000, '+', 2, 'Gene3', 'lncRNA'),
        (900000, 1100000, '+', 4, 'Gene4', 'protein_coding'),
        (1300000, 1500000, '-', 6, 'Gene5', 'protein_coding'),
        (1700000, 1850000, '+', 1, 'Gene6', 'protein_coding'),  # Single exon
        (2100000, 2400000, '+', 3, 'Gene7', 'protein_coding'),
        (2600000, 2900000, '-', 4, 'Gene8', 'lncRNA'),
    ]

    for gene_start, gene_end, strand, n_exons, gene_name, gene_type in gene_configs:
        gene_id = f"GENE{gene_name[-1]}001"
        transcript_id = f"TX{gene_name[-1]}001.1"
        
        # Calculate exon positions
        gene_length = gene_end - gene_start
        exon_length = gene_length // (n_exons * 2)
        intron_length = (gene_length - exon_length * n_exons) // (n_exons - 1) if n_exons > 1 else 0
        
        exons = []
        current_pos = gene_start
        for i in range(n_exons):
            exon_start = current_pos
            exon_end = current_pos + exon_length
            exons.append((exon_start, exon_end))
            current_pos = exon_end + intron_length
        
        # Gene feature (proper GTF format with single quotes around attributes)
        genes.append(f'{gene_start}\t{gene_end}\t+\t.\tgene_id "{gene_id}"; gene_name "{gene_name}"; gene_type "{gene_type}";')
        
        # Transcript feature
        genes.append(f'{gene_start}\t{gene_end}\t+\t.\tgene_id "{gene_id}"; gene_name "{gene_name}"; transcript_id "{transcript_id}";')
        
        # Exon features
        for i, (exon_start, exon_end) in enumerate(exons):
            genes.append(f'{exon_start}\t{exon_end}\t+\t.\tgene_id "{gene_id}"; gene_name "{gene_name}"; transcript_id "{transcript_id}"; exon_number "{i+1}";')
        
        # CDS features (for protein_coding genes)
        if gene_type == "protein_coding" and len(exons) >= 2:
            # CDS starts at second exon and ends at second-to-last
            for i, (exon_start, exon_end) in enumerate(exons[1:-1], 1):
                genes.append(f'{exon_start}\t{exon_end}\t+\t0\tgene_id "{gene_id}"; gene_name "{gene_name}"; transcript_id "{transcript_id}"; exon_number "{i+1}";')

    # Write GTF format with proper structure
    with open(GTF_FILE, 'w') as f:
        for i, (gene_start, gene_end, strand, n_exons, gene_name, gene_type) in enumerate(gene_configs):
            gene_id = f"GENE{gene_name[-1]}001"
            transcript_id = f"TX{gene_name[-1]}001.1"
            
            # Calculate exon positions
            gene_length = gene_end - gene_start
            exon_length = gene_length // (n_exons * 2)
            intron_length = (gene_length - exon_length * n_exons) // (n_exons - 1) if n_exons > 1 else 0
            
            exons = []
            current_pos = gene_start
            for j in range(n_exons):
                exon_start = current_pos
                exon_end = current_pos + exon_length
                exons.append((exon_start, exon_end))
                current_pos = exon_end + intron_length
            
            # Gene feature
            f.write(f"chr1\tsynthetic\tgene\t{gene_start}\t{gene_end}\t.\t{strand}\t.\tgene_id \"{gene_id}\"; gene_name \"{gene_name}\"; gene_type \"{gene_type}\";\n")
            
            # Transcript feature
            f.write(f"chr1\tsynthetic\ttranscript\t{gene_start}\t{gene_end}\t.\t{strand}\t.\tgene_id \"{gene_id}\"; gene_name \"{gene_name}\"; transcript_id \"{transcript_id}\";\n")
            
            # Exon features
            for j, (exon_start, exon_end) in enumerate(exons):
                f.write(f"chr1\tsynthetic\texon\t{exon_start}\t{exon_end}\t.\t{strand}\t.\tgene_id \"{gene_id}\"; gene_name \"{gene_name}\"; transcript_id \"{transcript_id}\"; exon_number \"{j+1}\";\n")
            
            # CDS features (for protein_coding genes)
            if gene_type == "protein_coding" and len(exons) >= 2:
                for j, (exon_start, exon_end) in enumerate(exons[1:-1], 1):
                    f.write(f"chr1\tsynthetic\tCDS\t{exon_start}\t{exon_end}\t.\t{strand}\t0\tgene_id \"{gene_id}\"; gene_name \"{gene_name}\"; transcript_id \"{transcript_id}\"; exon_number \"{j+1}\";\n")

    print(f"  Created {GTF_FILE}")


if __name__ == "__main__":
    print("=" * 60)
    print("Creating visual test data for trapezoid Hi-C + multi-track")
    print("=" * 60)
    
    create_cool_file()
    create_eigenvector_file()
    create_bigwig_files()
    create_gtf_file()
    
    print("\n" + "=" * 60)
    print("Done! Test files created:")
    print(f"  - {COOL_FILE}")
    print(f"  - {EIGS_FILE}")
    print(f"  - {BW_FILE1}")
    print(f"  - {BW_FILE2}")
    print(f"  - {GTF_FILE}")
    print("=" * 60)
