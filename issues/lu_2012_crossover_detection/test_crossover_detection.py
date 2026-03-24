#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test script for Crossover Detection Algorithm
Generates synthetic data and tests the detection pipeline.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from crossover_detection import CrossoverDetector, HMMCrossoverDetector


def generate_synthetic_data(
    n_snps: int = 10000,
    n_crossovers: int = 25,
    error_rate: float = 0.05,
    seed: int = 42
) -> tuple:
    """
    Generate synthetic single sperm sequencing data with known crossovers.
    
    Args:
        n_snps: Number of SNPs to generate
        n_crossovers: Number of crossover events to simulate
        error_rate: SNP calling error rate
        seed: Random seed for reproducibility
        
    Returns:
        Tuple of (positions, sperm_genotypes, haplotype_1, haplotype_2, true_crossovers)
    """
    np.random.seed(seed)
    
    # Generate SNP positions along a chromosome (0-250 Mb)
    positions = np.sort(np.random.uniform(0, 250_000_000, n_snps)).astype(int)
    
    # Generate haplotype alleles (random A/C/G/T)
    alleles = ['A', 'C', 'G', 'T']
    haplotype_1 = np.random.choice(alleles, n_snps)
    haplotype_2 = np.random.choice(alleles, n_snps)
    
    # Ensure haplotypes are different at each position (heterozygous)
    for i in range(n_snps):
        while haplotype_2[i] == haplotype_1[i]:
            haplotype_2[i] = np.random.choice(alleles)
    
    # Generate crossover positions
    true_crossovers = np.sort(np.random.uniform(
        positions[0] + 1_000_000,  # Avoid ends
        positions[-1] - 1_000_000,
        n_crossovers
    )).astype(int)
    
    # Generate sperm genotypes with crossovers
    sperm_genotypes = []
    current_haplotype = 0  # Start with haplotype 1 (paternal)
    crossover_idx = 0
    
    for i, pos in enumerate(positions):
        # Check if we crossed a crossover point
        while crossover_idx < len(true_crossovers) and pos > true_crossovers[crossover_idx]:
            current_haplotype = 1 - current_haplotype  # Switch haplotype
            crossover_idx += 1
        
        # Get genotype from current haplotype
        if current_haplotype == 0:
            gt = haplotype_1[i]
        else:
            gt = haplotype_2[i]
        
        # Add sequencing errors
        if np.random.random() < error_rate:
            gt = np.random.choice(alleles)
        
        sperm_genotypes.append(gt)
    
    return positions, sperm_genotypes, haplotype_1.tolist(), haplotype_2.tolist(), true_crossovers


def test_hmm_detector():
    """Test the HMM crossover detector with synthetic data."""
    print("=" * 60)
    print("Testing HMM Crossover Detector")
    print("=" * 60)
    
    # Generate synthetic data
    print("\nGenerating synthetic data...")
    positions, sperm_genotypes, haplotype_1, haplotype_2, true_crossovers = generate_synthetic_data(
        n_snps=5000,
        n_crossovers=20,
        error_rate=0.05,
        seed=42
    )
    
    print(f"  SNPs: {len(positions)}")
    print(f"  True crossovers: {len(true_crossovers)}")
    print(f"  Error rate: 5%")
    
    # Initialize detector
    detector = HMMCrossoverDetector(
        window_size=10,
        min_fragment_size=300_000,
        error_rate=0.05
    )
    
    # Run detection
    print("\nRunning crossover detection...")
    states, detected_crossovers = detector.detect(
        positions=positions,
        sperm_genotypes=sperm_genotypes,
        haplotype_1=haplotype_1,
        haplotype_2=haplotype_2
    )
    
    print(f"  Detected crossovers: {len(detected_crossovers)}")
    
    # Calculate performance metrics
    print("\nPerformance Metrics:")
    
    # True positives: detected crossovers within 100kb of true crossovers
    tolerance = 100_000  # 100 kb
    true_positives = 0
    for dc in detected_crossovers:
        for tc in true_crossovers:
            if abs(dc - tc) < tolerance:
                true_positives += 1
                break
    
    false_positives = len(detected_crossovers) - true_positives
    false_negatives = len(true_crossovers) - true_positives
    
    precision = true_positives / (true_positives + false_positives) if (true_positives + false_positives) > 0 else 0
    recall = true_positives / (true_positives + false_negatives) if (true_positives + false_negatives) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    
    print(f"  True Positives: {true_positives}")
    print(f"  False Positives: {false_positives}")
    print(f"  False Negatives: {false_negatives}")
    print(f"  Precision: {precision:.3f}")
    print(f"  Recall: {recall:.3f}")
    print(f"  F1 Score: {f1:.3f}")
    
    # Print detected vs true crossover positions
    print("\nCrossover Positions Comparison:")
    print(f"  True:      {sorted(true_crossovers[:10])}...")  # First 10
    print(f"  Detected:  {sorted(detected_crossovers[:10])}...")  # First 10
    
    return {
        'precision': precision,
        'recall': recall,
        'f1': f1,
        'true_crossovers': len(true_crossovers),
        'detected_crossovers': len(detected_crossovers)
    }


def test_full_pipeline():
    """Test the full crossover detection pipeline with synthetic data files."""
    print("\n" + "=" * 60)
    print("Testing Full Pipeline")
    print("=" * 60)
    
    # Create test directory
    test_dir = Path(__file__).parent / "test_output"
    test_dir.mkdir(exist_ok=True)
    
    # Generate and save synthetic data
    print("\nGenerating test data files...")
    positions, sperm_genotypes, haplotype_1, haplotype_2, _ = generate_synthetic_data(
        n_snps=3000,
        n_crossovers=15,
        seed=123
    )
    
    # Create SNP data DataFrame
    snp_df = pd.DataFrame({
        'chromosome': ['chr1'] * len(positions),
        'position': positions,
        'snp_id': [f'rs{i}' for i in range(len(positions))],
        'ref_allele': ['A'] * len(positions),
        'alt_allele': ['G'] * len(positions),
        'sperm_genotype': sperm_genotypes,
        'diploid_genotype': ['A/G'] * len(positions),
        'sperm_id': ['TEST001'] * len(positions)
    })
    
    # Create phasing data DataFrame
    phasing_df = pd.DataFrame({
        'chromosome': ['chr1'] * len(positions),
        'position': positions,
        'snp_id': [f'rs{i}' for i in range(len(positions))],
        'haplotype_1': haplotype_1,
        'haplotype_2': haplotype_2
    })
    
    # Save files
    snp_file = test_dir / "test_snps.tsv"
    phasing_file = test_dir / "test_phasing.tsv"
    
    snp_df.to_csv(snp_file, sep='\t', index=False)
    phasing_df.to_csv(phasing_file, sep='\t', index=False)
    
    print(f"  Saved SNP data: {snp_file}")
    print(f"  Saved phasing data: {phasing_file}")
    
    # Run detection
    print("\nRunning full pipeline...")
    detector = CrossoverDetector(
        window_size=10,
        min_fragment_size=300_000,
        error_rate=0.05
    )
    
    detector.load_snp_data(str(snp_file))
    detector.load_phasing_data(str(phasing_file))
    
    results = detector.detect_crossovers(sperm_id='TEST001', chromosome='chr1')
    
    # Save results
    detector.save_results(str(test_dir))
    
    print(f"\nResults saved to: {test_dir}")
    print(f"\nStatistics:")
    for key, value in results.statistics.items():
        print(f"  {key}: {value}")
    
    return results


def main():
    """Run all tests."""
    print("\n" + "=" * 60)
    print("Crossover Detection Algorithm - Test Suite")
    print("=" * 60)
    
    # Test HMM detector
    hmm_results = test_hmm_detector()
    
    # Test full pipeline
    pipeline_results = test_full_pipeline()
    
    # Summary
    print("\n" + "=" * 60)
    print("Test Summary")
    print("=" * 60)
    print(f"\nHMM Detector:")
    print(f"  Precision: {hmm_results['precision']:.3f}")
    print(f"  Recall: {hmm_results['recall']:.3f}")
    print(f"  F1 Score: {hmm_results['f1']:.3f}")
    
    print(f"\nFull Pipeline:")
    print(f"  Total crossovers detected: {pipeline_results.statistics['total_crossovers']}")
    print(f"  Total haplotype blocks: {pipeline_results.statistics['total_blocks']}")
    print(f"  Mean confidence: {pipeline_results.statistics['mean_confidence']:.3f}")
    
    print("\n" + "=" * 60)
    print("All tests completed!")
    print("=" * 60)


if __name__ == '__main__':
    main()
