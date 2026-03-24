#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Visualization script for crossover detection validation.

This script generates synthetic data with known crossover positions,
runs the crossover detection module, and creates a scatter plot
showing the haplotype calls with detected crossovers marked.

Usage:
    python visualize_crossover_validation.py [--output-dir OUTPUT_DIR]

Example:
    python visualize_crossover_validation.py --output-dir test_output/
"""

import argparse
import logging
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from hic_basic.crossover import CrossoverDetector

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def generate_synthetic_data(
    n_snps: int = 1000,
    crossover_positions: list = None,
    error_rate: float = 0.05,
    seed: int = 42
) -> tuple:
    """
    Generate synthetic SNP data with known crossover positions.

    Args:
        n_snps: Number of SNPs to generate
        crossover_positions: List of positions where crossovers occur
        error_rate: Rate of genotyping errors
        seed: Random seed for reproducibility

    Returns:
        Tuple of (snp_df, phasing_df, true_haplotypes)
    """
    np.random.seed(seed)

    if crossover_positions is None:
        # Default: 3 crossovers at ~25%, 50%, 75% of chromosome
        crossover_positions = [
            int(n_snps * 0.25),
            int(n_snps * 0.50),
            int(n_snps * 0.75),
        ]

    # Generate positions
    positions = np.sort(np.random.choice(1000000, n_snps, replace=False))

    # Generate haplotype alleles (reference)
    haplotype_1 = np.random.choice(['A', 'G'], n_snps)
    haplotype_2 = np.where(haplotype_1 == 'A', 'G', 'A')

    # Generate true haplotype path with crossovers
    true_haplotypes = np.zeros(n_snps, dtype=int)  # 0 = maternal (haplotype_2), 1 = paternal (haplotype_1)

    # Start with paternal
    current_hap = 1
    for i in range(n_snps):
        if i in crossover_positions:
            current_hap = 1 - current_hap  # Switch haplotype
        true_haplotypes[i] = current_hap

    # Generate sperm genotypes based on true haplotypes
    sperm_genotypes = []
    for i in range(n_snps):
        if np.random.random() < error_rate:
            # Genotyping error
            sperm_gt = np.random.choice(['A', 'G', 'N'])
        else:
            # Correct genotype based on haplotype
            if true_haplotypes[i] == 1:
                sperm_gt = haplotype_1[i]
            else:
                sperm_gt = haplotype_2[i]
        sperm_genotypes.append(sperm_gt)

    # Create DataFrames
    snp_df = pd.DataFrame({
        'chromosome': ['chr1'] * n_snps,
        'position': positions,
        'snp_id': [f'rs{i}' for i in range(n_snps)],
        'sperm_genotype': sperm_genotypes,
        'diploid_genotype': ['A/G'] * n_snps,
    })

    phasing_df = pd.DataFrame({
        'chromosome': ['chr1'] * n_snps,
        'position': positions,
        'snp_id': [f'rs{i}' for i in range(n_snps)],
        'haplotype_1': haplotype_1,
        'haplotype_2': haplotype_2,
    })

    return snp_df, phasing_df, true_haplotypes, crossover_positions


def plot_crossover_validation(
    snp_df: pd.DataFrame,
    phasing_df: pd.DataFrame,
    true_haplotypes: np.ndarray,
    detected_crossovers: list,
    true_crossover_positions: list,
    output_file: str,
    sample_id: str = "synthetic"
) -> None:
    """
    Create a scatter plot showing haplotype calls and detected crossovers.

    Args:
        snp_df: SNP genotype DataFrame
        phasing_df: Phasing information DataFrame
        true_haplotypes: True haplotype assignments (0=maternal, 1=paternal)
        detected_crossovers: List of detected Crossover objects
        true_crossover_positions: List of true crossover SNP indices
        output_file: Output file path
        sample_id: Sample identifier
    """
    # Merge data
    merged = snp_df.merge(phasing_df, on=['chromosome', 'position', 'snp_id'])

    # Encode observed haplotypes from sperm genotypes
    observed_haplotypes = []
    for _, row in merged.iterrows():
        gt = row['sperm_genotype']
        h1, h2 = row['haplotype_1'], row['haplotype_2']

        if gt == h1:
            observed_haplotypes.append(1)  # Paternal
        elif gt == h2:
            observed_haplotypes.append(0)  # Maternal
        else:
            observed_haplotypes.append(np.nan)  # Unknown/ambiguous

    observed_haplotypes = np.array(observed_haplotypes)

    # Create figure with two panels
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10),
                                    gridspec_kw={'height_ratios': [3, 1]},
                                    sharex=True)

    positions = merged['position'].values

    # Panel 1: True and Observed haplotypes
    ax1.scatter(positions, true_haplotypes, c='lightblue', s=10, alpha=0.5,
                label='True haplotype', edgecolors='none')
    ax1.scatter(positions, observed_haplotypes, c='steelblue', s=8, alpha=0.7,
                label='Observed haplotype', edgecolors='none', marker='o')

    # Mark true crossover positions
    for i, pos_idx in enumerate(true_crossover_positions):
        if pos_idx < len(positions):
            ax1.axvline(x=positions[pos_idx], color='red', linestyle='--',
                        linewidth=2, alpha=0.7, label='True CO' if i == 0 else None)

    ax1.set_ylabel('Haplotype\n(0=Mat, 1=Pat)', fontsize=11)
    ax1.set_title(f'Crossover Detection Validation - {sample_id}', fontsize=14, fontweight='bold')
    ax1.legend(loc='upper right', fontsize=10)
    ax1.set_yticks([0, 1])
    ax1.set_yticklabels(['Maternal (0)', 'Paternal (1)'])
    ax1.grid(True, alpha=0.3, linestyle='-')
    ax1.set_ylim(-0.3, 1.3)

    # Panel 2: Haplotype difference (to show transitions)
    # Calculate running window haplotype frequency
    window_size = 20
    hap_freq = np.convolve(observed_haplotypes, np.ones(window_size)/window_size, mode='same')

    ax2.plot(positions, hap_freq, color='darkgreen', linewidth=2, label='Paternal frequency')
    ax2.axhline(y=0.5, color='gray', linestyle='-', linewidth=1, alpha=0.5)

    # Mark detected crossovers
    for i, xo in enumerate(detected_crossovers):
        ax2.axvline(x=xo.position, color='purple', linestyle='-',
                    linewidth=3, alpha=0.8, label='Detected CO' if i == 0 else None)
        # Mark with star
        ax2.scatter([xo.position], [0.5], c='purple', s=200, marker='*',
                    zorder=5, edgecolors='white', linewidths=1.5)

    ax2.set_xlabel('Genomic Position', fontsize=11)
    ax2.set_ylabel('Paternal\nFrequency', fontsize=11)
    ax2.legend(loc='upper right', fontsize=10)
    ax2.grid(True, alpha=0.3, linestyle='-')
    ax2.set_ylim(-0.1, 1.1)

    # Add statistics text box
    n_detected = len(detected_crossovers)
    n_true = len(true_crossover_positions)

    # Count how many detected are close to true
    true_positions = set(positions[p] for p in true_crossover_positions if p < len(positions))
    n_correct = sum(1 for xo in detected_crossovers
                    if any(abs(xo.position - tp) < 50000 for tp in true_positions))

    stats_text = (
        f"Statistics:\n"
        f"True crossovers: {n_true}\n"
        f"Detected crossovers: {n_detected}\n"
        f"Correctly detected: {n_correct}\n"
        f"SNPs analyzed: {len(positions):,}"
    )

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax1.text(0.02, 0.98, stats_text, transform=ax1.transAxes, fontsize=10,
             verticalalignment='top', bbox=props, family='monospace')

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved plot to {output_file}")


def run_validation(
    output_dir: str,
    n_snps: int = 1000,
    crossover_positions: list = None,
    error_rate: float = 0.05,
    seed: int = 42
) -> dict:
    """
    Run complete validation pipeline.

    Args:
        output_dir: Output directory for results
        n_snps: Number of SNPs
        crossover_positions: True crossover positions (SNP indices)
        error_rate: Genotyping error rate
        seed: Random seed

    Returns:
        Dictionary with validation results
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    logger.info(f"Generating synthetic data with {n_snps} SNPs...")
    if crossover_positions is None:
        crossover_positions = [
            int(n_snps * 0.25),
            int(n_snps * 0.50),
            int(n_snps * 0.75),
        ]
    logger.info(f"True crossover positions (SNP indices): {crossover_positions}")

    # Generate data
    snp_df, phasing_df, true_haplotypes, true_co_positions = generate_synthetic_data(
        n_snps=n_snps,
        crossover_positions=crossover_positions,
        error_rate=error_rate,
        seed=seed
    )

    # Save input data
    snp_file = output_path / "synthetic_snps.tsv"
    phasing_file = output_path / "synthetic_phasing.tsv"
    snp_df.to_csv(snp_file, sep='\t', index=False)
    phasing_df.to_csv(phasing_file, sep='\t', index=False)
    logger.info(f"Saved input data to {output_path}")

    # Run crossover detection
    logger.info("Running crossover detection...")
    detector = CrossoverDetector(
        window_size=10,
        min_fragment_size=50000,  # 50kb for synthetic data
        error_rate=error_rate
    )

    detector.load_snp_data(str(snp_file))
    detector.load_phasing_data(str(phasing_file))
    results = detector.detect_crossovers("synthetic")

    # Save results
    detector.save_results(str(output_path), sperm_id="synthetic")
    logger.info(f"Detected {len(results.crossovers)} crossovers")

    # Create visualization
    plot_file = output_path / "crossover_validation.png"
    plot_crossover_validation(
        snp_df=snp_df,
        phasing_df=phasing_df,
        true_haplotypes=true_haplotypes,
        detected_crossovers=results.crossovers,
        true_crossover_positions=true_co_positions,
        output_file=str(plot_file),
        sample_id="synthetic"
    )

    # Calculate validation metrics
    true_positions = set(snp_df.iloc[p]['position'] for p in true_co_positions if p < len(snp_df))

    n_correct = 0
    for xo in results.crossovers:
        if any(abs(xo.position - tp) < 50000 for tp in true_positions):
            n_correct += 1

    validation_results = {
        'n_snps': n_snps,
        'n_true_crossovers': len(true_co_positions),
        'n_detected_crossovers': len(results.crossovers),
        'n_correctly_detected': n_correct,
        'precision': n_correct / len(results.crossovers) if results.crossovers else 0,
        'recall': n_correct / len(true_co_positions),
        'detected_positions': [(xo.chromosome, xo.position) for xo in results.crossovers],
        'true_positions': [(snp_df.iloc[p]['chromosome'], snp_df.iloc[p]['position'])
                          for p in true_co_positions if p < len(snp_df)],
    }

    # Save validation summary
    summary_file = output_path / "validation_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("Crossover Detection Validation Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Input:\n")
        f.write(f"  Number of SNPs: {validation_results['n_snps']:,}\n")
        f.write(f"  Error rate: {error_rate}\n\n")
        f.write(f"Results:\n")
        f.write(f"  True crossovers: {validation_results['n_true_crossovers']}\n")
        f.write(f"  Detected crossovers: {validation_results['n_detected_crossovers']}\n")
        f.write(f"  Correctly detected: {validation_results['n_correctly_detected']}\n\n")
        f.write(f"Metrics:\n")
        f.write(f"  Precision: {validation_results['precision']:.2%}\n")
        f.write(f"  Recall: {validation_results['recall']:.2%}\n\n")
        f.write(f"True crossover positions:\n")
        for chrom, pos in validation_results['true_positions']:
            f.write(f"  {chrom}:{pos:,}\n")
        f.write(f"\nDetected crossover positions:\n")
        for chrom, pos in validation_results['detected_positions']:
            f.write(f"  {chrom}:{pos:,}\n")

    logger.info(f"Validation summary saved to {summary_file}")

    return validation_results


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Visualize crossover detection validation with synthetic data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run with default parameters (3 crossovers)
  python visualize_crossover_validation.py --output-dir test_output/

  # Custom number of SNPs and error rate
  python visualize_crossover_validation.py --output-dir test_output/ --n-snps 2000 --error-rate 0.03

  # Custom crossover positions
  python visualize_crossover_validation.py --output-dir test_output/ --crossover-positions 100 500 900
        """
    )

    parser.add_argument(
        '--output-dir', '-o',
        default='test_output',
        help='Output directory for results (default: test_output)'
    )
    parser.add_argument(
        '--n-snps',
        type=int,
        default=1000,
        help='Number of SNPs to generate (default: 1000)'
    )
    parser.add_argument(
        '--error-rate',
        type=float,
        default=0.05,
        help='Genotyping error rate (default: 0.05)'
    )
    parser.add_argument(
        '--seed',
        type=int,
        default=42,
        help='Random seed (default: 42)'
    )
    parser.add_argument(
        '--crossover-positions',
        type=int,
        nargs='+',
        default=None,
        help='Crossover positions as SNP indices (default: 25%%, 50%%, 75%%)'
    )

    args = parser.parse_args()

    # Run validation
    results = run_validation(
        output_dir=args.output_dir,
        n_snps=args.n_snps,
        crossover_positions=args.crossover_positions,
        error_rate=args.error_rate,
        seed=args.seed
    )

    # Print summary
    print("\n" + "=" * 60)
    print("CROSSOVER DETECTION VALIDATION COMPLETE")
    print("=" * 60)
    print(f"\nTrue crossovers:      {results['n_true_crossovers']}")
    print(f"Detected crossovers:  {results['n_detected_crossovers']}")
    print(f"Correctly detected:   {results['n_correctly_detected']}")
    print(f"\nPrecision: {results['precision']:.2%}")
    print(f"Recall:    {results['recall']:.2%}")
    print(f"\nOutput files:")
    print(f"  - {args.output_dir}/crossover_validation.png")
    print(f"  - {args.output_dir}/synthetic_snps.tsv")
    print(f"  - {args.output_dir}/synthetic_phasing.tsv")
    print(f"  - {args.output_dir}/validation_summary.txt")
    print("=" * 60 + "\n")

    return 0 if results['recall'] > 0.5 else 1


if __name__ == '__main__':
    sys.exit(main())
