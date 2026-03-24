#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Crossover Detection Algorithm
Implementation of Lu et al. (2012) Science method for detecting meiotic 
crossover events from single sperm whole-genome sequencing data.

Reference:
    Lu S, Zong C, Fan W, et al. Probing Meiotic Recombination and Aneuploidy 
    of Single Sperm Cells by Whole-Genome Sequencing. Science. 2012;338(6114):1627-1630.

Author: Generated based on Lu et al. (2012) methods
Date: 2026-03-23
"""

import argparse
import logging
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.stats import binom

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# ============================================================================
# Data Classes
# ============================================================================

@dataclass
class SNPData:
    """Container for SNP genotype data."""
    chromosome: str
    position: int
    snp_id: str
    ref_allele: str
    alt_allele: str
    sperm_genotype: str
    diploid_genotype: str
    haplotype_1: Optional[str] = None
    haplotype_2: Optional[str] = None


@dataclass
class HaplotypeBlock:
    """Container for a haplotype block."""
    chromosome: str
    start_pos: int
    end_pos: int
    haplotype: str  # 'P' (paternal), 'M' (maternal), or 'U' (undetermined)
    confidence: float
    num_snps: int


@dataclass
class Crossover:
    """Container for a crossover event."""
    sperm_id: str
    chromosome: str
    position: int
    start_pos: int
    end_pos: int
    from_haplotype: str
    to_haplotype: str
    confidence: float
    resolution_kb: float


@dataclass
class CrossoverResults:
    """Container for crossover detection results."""
    sperm_id: str
    haplotype_blocks: List[HaplotypeBlock] = field(default_factory=list)
    crossovers: List[Crossover] = field(default_factory=list)
    statistics: Dict = field(default_factory=dict)


# ============================================================================
# HMM Implementation
# ============================================================================

class HMMCrossoverDetector:
    """
    Hidden Markov Model for crossover detection.
    
    States:
        P (0): Paternal haplotype
        M (1): Maternal haplotype
        U (2): Undetermined haplotype
    
    This implements the two-step method from Lu et al. (2012):
    1. Rough estimation using SNP ratio in windows
    2. Accurate determination using HMM with Viterbi algorithm
    """
    
    # State indices
    STATE_PATERNAL = 0
    STATE_MATERNAL = 1
    STATE_UNDETERMINED = 2
    
    STATE_NAMES = {0: 'P', 1: 'M', 2: 'U'}
    
    def __init__(
        self,
        window_size: int = 10,
        min_fragment_size: int = 300000,
        error_rate: float = 0.05,
        crossover_rate: float = 0.01
    ):
        """
        Initialize HMM detector.
        
        Args:
            window_size: Number of SNPs per window for rough estimation
            min_fragment_size: Minimum fragment size in bp (filter small fragments)
            error_rate: SNP calling error rate for emission probabilities
            crossover_rate: Prior crossover rate for transition probabilities
        """
        self.window_size = window_size
        self.min_fragment_size = min_fragment_size
        self.error_rate = error_rate
        self.crossover_rate = crossover_rate
        
        # Initialize HMM parameters
        self._init_transition_matrix()
        self._init_emission_matrix()
    
    def _init_transition_matrix(self):
        """
        Initialize transition probability matrix.
        
        Based on rough crossover frequency from Step 1.
        """
        # Transition probabilities
        # P->M and M->P transitions represent crossovers
        p_crossover = self.crossover_rate / 2
        p_stay = 1 - p_crossover
        p_undetermined = 0.01  # Small probability to enter/exit undetermined state
        
        self.transition_matrix = np.array([
            # P       M       U
            [p_stay, p_crossover, p_undetermined],  # From P
            [p_crossover, p_stay, p_undetermined],  # From M
            [p_undetermined, p_undetermined, 1 - 2 * p_undetermined],  # From U
        ])
        
        # Ensure rows sum to 1
        self.transition_matrix = self.transition_matrix / self.transition_matrix.sum(axis=1, keepdims=True)
    
    def _init_emission_matrix(self):
        """
        Initialize emission probability matrix.
        
        Based on SNP calling error rate.
        Observations:
            0: Paternal allele observed
            1: Maternal allele observed
            2: Undetermined/heterozygous
        """
        e = self.error_rate
        
        self.emission_matrix = np.array([
            # Paternal_obs  Maternal_obs  Undetermined_obs
            [1 - e, e, 0.01],  # State P (Paternal)
            [e, 1 - e, 0.01],  # State M (Maternal)
            [0.01, 0.01, 0.98],  # State U (Undetermined)
        ])
        
        # Ensure rows sum to 1
        self.emission_matrix = self.emission_matrix / self.emission_matrix.sum(axis=1, keepdims=True)
    
    def _rough_estimation(
        self,
        positions: np.ndarray,
        haplotype_alleles: np.ndarray
    ) -> List[int]:
        """
        Step 1: Rough estimation of crossover positions.
        
        Calculate SNP ratio belonging to one parental haplotype over the other
        in windows of SNPs, and identify transition points.
        
        Args:
            positions: SNP positions
            haplotype_alleles: Array of 0 (paternal) or 1 (maternal) for each SNP
            
        Returns:
            List of rough crossover positions
        """
        n_snps = len(positions)
        if n_snps < self.window_size * 2:
            return []
        
        # Calculate paternal allele ratio in sliding windows
        n_windows = n_snps - self.window_size + 1
        paternal_ratios = np.zeros(n_windows)
        window_positions = np.zeros(n_windows)
        
        for i in range(n_windows):
            window = haplotype_alleles[i:i + self.window_size]
            valid_snps = window[window >= 0]  # Exclude undetermined
            if len(valid_snps) > 0:
                paternal_ratios[i] = np.sum(valid_snps == 0) / len(valid_snps)
            else:
                paternal_ratios[i] = 0.5
            window_positions[i] = np.mean(positions[i:i + self.window_size])
        
        # Identify transition points
        crossover_positions = []
        threshold = 0.5  # Ratio change threshold
        
        for i in range(1, len(paternal_ratios)):
            ratio_change = abs(paternal_ratios[i] - paternal_ratios[i-1])
            if ratio_change > threshold:
                # Transition point detected
                crossover_pos = int((window_positions[i] + window_positions[i-1]) / 2)
                crossover_positions.append(crossover_pos)
        
        logger.debug(f"Rough estimation found {len(crossover_positions)} candidate crossovers")
        return crossover_positions
    
    def _viterbi(
        self,
        observations: np.ndarray
    ) -> Tuple[np.ndarray, float]:
        """
        Step 2: Viterbi algorithm for most probable state sequence.
        
        Args:
            observations: Array of observed alleles (0=paternal, 1=maternal, 2=undetermined)
            
        Returns:
            Tuple of (best_path, log_probability)
        """
        n_states = 3
        n_obs = len(observations)
        
        # Initialize Viterbi matrices
        # delta[t, i] = max probability of being in state i at time t
        delta = np.zeros((n_obs, n_states))
        # psi[t, i] = argmax state at t-1 leading to state i at time t
        psi = np.zeros((n_obs, n_states), dtype=int)
        
        # Initial probabilities (uniform)
        delta[0, :] = 1.0 / n_states * self.emission_matrix[:, observations[0]]
        
        # Forward pass
        for t in range(1, n_obs):
            for j in range(n_states):
                # Calculate delta[t, j] = max_i(delta[t-1, i] * a_ij) * b_j(o_t)
                probs = delta[t-1, :] * self.transition_matrix[:, j]
                psi[t, j] = np.argmax(probs)
                delta[t, j] = probs[psi[t, j]] * self.emission_matrix[j, observations[t]]
            
            # Normalize to prevent underflow
            delta[t, :] /= np.sum(delta[t, :]) + 1e-300
        
        # Backtrack to find best path
        best_path = np.zeros(n_obs, dtype=int)
        best_path[-1] = np.argmax(delta[-1, :])
        log_prob = np.log(np.max(delta[-1, :]) + 1e-300)
        
        for t in range(n_obs - 2, -1, -1):
            best_path[t] = psi[t + 1, best_path[t + 1]]
        
        return best_path, log_prob
    
    def _encode_observations(
        self,
        sperm_genotypes: List[str],
        haplotype_1: List[str],
        haplotype_2: List[str]
    ) -> np.ndarray:
        """
        Encode genotypes as observations for HMM.
        
        Args:
            sperm_genotypes: Sperm genotypes (e.g., 'A', 'G', 'N' for unknown)
            haplotype_1: Paternal haplotype alleles
            haplotype_2: Maternal haplotype alleles
            
        Returns:
            Array of encoded observations (0=paternal, 1=maternal, 2=undetermined)
        """
        observations = []
        
        for sperm_gt, h1, h2 in zip(sperm_genotypes, haplotype_1, haplotype_2):
            if sperm_gt in ['N', '.', '-'] or len(sperm_gt) == 0:
                observations.append(2)  # Undetermined
            elif sperm_gt == h1:
                observations.append(0)  # Paternal
            elif sperm_gt == h2:
                observations.append(1)  # Maternal
            else:
                observations.append(2)  # Mismatch/undetermined
        
        return np.array(observations)
    
    def _filter_small_fragments(
        self,
        positions: np.ndarray,
        states: np.ndarray
    ) -> np.ndarray:
        """
        Filter out small conflicting fragments (< min_fragment_size).
        
        Args:
            positions: SNP positions
            states: HMM state assignments
            
        Returns:
            Filtered state assignments
        """
        filtered_states = states.copy()
        n_snps = len(positions)
        
        i = 0
        while i < n_snps:
            # Find fragment boundaries
            current_state = states[i]
            start_idx = i
            
            while i < n_snps and states[i] == current_state:
                i += 1
            
            end_idx = i - 1
            
            # Calculate fragment size
            if end_idx > start_idx:
                fragment_size = positions[end_idx] - positions[start_idx]
                
                # Filter small fragments (but not P/U/M transitions at ends)
                if fragment_size < self.min_fragment_size and current_state == self.STATE_UNDETERMINED:
                    # Assign to neighboring state
                    if start_idx > 0:
                        filtered_states[start_idx:end_idx+1] = states[start_idx - 1]
                    elif end_idx < n_snps - 1:
                        filtered_states[start_idx:end_idx+1] = states[end_idx + 1]
        
        return filtered_states
    
    def detect(
        self,
        positions: np.ndarray,
        sperm_genotypes: List[str],
        haplotype_1: List[str],
        haplotype_2: List[str]
    ) -> Tuple[np.ndarray, List[int]]:
        """
        Run complete crossover detection pipeline.
        
        Args:
            positions: SNP positions (sorted)
            sperm_genotypes: Sperm genotypes at each SNP
            haplotype_1: Paternal haplotype alleles
            haplotype_2: Maternal haplotype alleles
            
        Returns:
            Tuple of (state_sequence, crossover_positions)
        """
        # Step 1: Rough estimation
        haplotype_alleles = self._encode_observations(sperm_genotypes, haplotype_1, haplotype_2)
        rough_crossovers = self._rough_estimation(positions, haplotype_alleles)
        
        # Update crossover rate based on rough estimation
        if len(positions) > 0 and len(rough_crossovers) > 0:
            genome_span = positions[-1] - positions[0]
            if genome_span > 0:
                self.crossover_rate = len(rough_crossovers) / (genome_span / 1e6)  # per Mb
                self._init_transition_matrix()
        
        # Step 2: HMM with Viterbi
        observations = haplotype_alleles
        states, log_prob = self._viterbi(observations)
        
        # Filter small fragments
        filtered_states = self._filter_small_fragments(positions, states)
        
        # Identify crossover positions from state transitions
        crossover_positions = []
        for i in range(1, len(filtered_states)):
            if filtered_states[i] != filtered_states[i-1]:
                if filtered_states[i] in [self.STATE_PATERNAL, self.STATE_MATERNAL] and \
                   filtered_states[i-1] in [self.STATE_PATERNAL, self.STATE_MATERNAL]:
                    # Crossover between P and M
                    crossover_pos = int((positions[i] + positions[i-1]) / 2)
                    crossover_positions.append(crossover_pos)
        
        logger.debug(f"HMM detected {len(crossover_positions)} crossovers")
        return filtered_states, crossover_positions


# ============================================================================
# Main Detection Class
# ============================================================================

class CrossoverDetector:
    """
    Main class for crossover detection from single sperm sequencing data.
    """
    
    def __init__(
        self,
        window_size: int = 10,
        min_fragment_size: int = 300000,
        error_rate: float = 0.05
    ):
        """
        Initialize crossover detector.
        
        Args:
            window_size: SNP window size for rough estimation
            min_fragment_size: Minimum fragment size to keep (bp)
            error_rate: SNP calling error rate
        """
        self.window_size = window_size
        self.min_fragment_size = min_fragment_size
        self.error_rate = error_rate
        
        self.hmm_detector = HMMCrossoverDetector(
            window_size=window_size,
            min_fragment_size=min_fragment_size,
            error_rate=error_rate
        )
        
        self.snp_data: Optional[pd.DataFrame] = None
        self.phasing_data: Optional[pd.DataFrame] = None
        self.results: Dict[str, CrossoverResults] = {}
    
    def load_snp_data(self, snp_file: str) -> None:
        """
        Load SNP genotype data.
        
        Args:
            snp_file: Path to TSV file with SNP genotypes
        """
        logger.info(f"Loading SNP data from {snp_file}")
        self.snp_data = pd.read_csv(snp_file, sep='\t')
        logger.info(f"Loaded {len(self.snp_data)} SNPs")
    
    def load_phasing_data(self, phasing_file: str) -> None:
        """
        Load phasing information.
        
        Args:
            phasing_file: Path to TSV file with haplotype phasing
        """
        logger.info(f"Loading phasing data from {phasing_file}")
        self.phasing_data = pd.read_csv(phasing_file, sep='\t')
        logger.info(f"Loaded phasing for {len(self.phasing_data)} SNPs")
    
    def detect_crossovers(
        self,
        sperm_id: str,
        chromosome: Optional[str] = None
    ) -> CrossoverResults:
        """
        Detect crossovers for a specific sperm sample.
        
        Args:
            sperm_id: ID of sperm sample to analyze
            chromosome: Optional chromosome to analyze (default: all)
            
        Returns:
            CrossoverResults object
        """
        logger.info(f"Detecting crossovers for sperm {sperm_id}")
        
        if self.snp_data is None:
            raise ValueError("SNP data not loaded. Call load_snp_data() first.")
        
        results = CrossoverResults(sperm_id=sperm_id)
        
        # Get unique chromosomes
        chromosomes = [chromosome] if chromosome else self.snp_data['chromosome'].unique()
        
        total_crossovers = 0
        total_blocks = 0
        
        for chrom in chromosomes:
            logger.debug(f"Processing chromosome {chrom}")
            
            # Filter data for this chromosome and sperm
            chrom_data = self.snp_data[self.snp_data['chromosome'] == chrom].copy()
            
            if len(chrom_data) == 0:
                continue
            
            # Merge with phasing data
            if self.phasing_data is not None:
                chrom_data = chrom_data.merge(
                    self.phasing_data,
                    on=['chromosome', 'position', 'snp_id'],
                    how='left'
                )
            
            # Extract arrays for HMM
            positions = chrom_data['position'].values
            sperm_genotypes = chrom_data['sperm_genotype'].tolist()
            
            # Get haplotype alleles
            if 'haplotype_1' in chrom_data.columns:
                haplotype_1 = chrom_data['haplotype_1'].fillna('N').tolist()
                haplotype_2 = chrom_data['haplotype_2'].fillna('N').tolist()
            else:
                # Infer from diploid genotype
                haplotype_1 = []
                haplotype_2 = []
                for _, row in chrom_data.iterrows():
                    gt = row['diploid_genotype']
                    if '/' in gt or '|' in gt:
                        alleles = gt.replace('|', '/').split('/')
                        haplotype_1.append(alleles[0])
                        haplotype_2.append(alleles[1] if len(alleles) > 1 else alleles[0])
                    else:
                        haplotype_1.append(gt)
                        haplotype_2.append(gt)
            
            # Run HMM detection
            states, crossover_positions = self.hmm_detector.detect(
                positions=positions,
                sperm_genotypes=sperm_genotypes,
                haplotype_1=haplotype_1,
                haplotype_2=haplotype_2
            )
            
            # Convert states to haplotype blocks
            blocks = self._states_to_blocks(
                chromosome=chrom,
                positions=positions,
                states=states
            )
            results.haplotype_blocks.extend(blocks)
            total_blocks += len(blocks)
            
            # Create crossover objects
            for i, pos in enumerate(crossover_positions):
                # Determine haplotype transition
                if i < len(blocks) - 1:
                    from_hap = blocks[i].haplotype
                    to_hap = blocks[i + 1].haplotype
                    confidence = (blocks[i].confidence + blocks[i + 1].confidence) / 2
                    resolution = (blocks[i + 1].start_pos - blocks[i].end_pos) / 1000
                else:
                    from_hap = 'U'
                    to_hap = 'U'
                    confidence = 0.5
                    resolution = 0
                
                crossover = Crossover(
                    sperm_id=sperm_id,
                    chromosome=chrom,
                    position=pos,
                    start_pos=blocks[i].end_pos if i < len(blocks) else positions[0],
                    end_pos=blocks[i + 1].start_pos if i < len(blocks) - 1 else positions[-1],
                    from_haplotype=from_hap,
                    to_haplotype=to_hap,
                    confidence=confidence,
                    resolution_kb=resolution
                )
                results.crossovers.append(crossover)
                total_crossovers += 1
        
        # Calculate statistics
        results.statistics = {
            'total_crossovers': total_crossovers,
            'total_blocks': total_blocks,
            'chromosomes_analyzed': len(chromosomes),
            'mean_resolution_kb': np.mean([c.resolution_kb for c in results.crossovers]) if results.crossovers else 0,
            'mean_confidence': np.mean([c.confidence for c in results.crossovers]) if results.crossovers else 0,
        }
        
        self.results[sperm_id] = results
        logger.info(f"Detected {total_crossovers} crossovers in {total_blocks} haplotype blocks")
        
        return results
    
    def _states_to_blocks(
        self,
        chromosome: str,
        positions: np.ndarray,
        states: np.ndarray
    ) -> List[HaplotypeBlock]:
        """
        Convert HMM state sequence to haplotype blocks.
        
        Args:
            chromosome: Chromosome name
            positions: SNP positions
            states: HMM state assignments
            
        Returns:
            List of HaplotypeBlock objects
        """
        blocks = []
        n_snps = len(positions)
        
        i = 0
        while i < n_snps:
            current_state = states[i]
            start_idx = i
            start_pos = positions[i]
            
            # Find end of current block
            while i < n_snps and states[i] == current_state:
                i += 1
            
            end_idx = i - 1
            end_pos = positions[end_idx]
            num_snps = end_idx - start_idx + 1
            
            # Calculate confidence (proportion of consistent SNPs)
            block_states = states[start_idx:end_idx + 1]
            confidence = np.sum(block_states == current_state) / len(block_states)
            
            block = HaplotypeBlock(
                chromosome=chromosome,
                start_pos=start_pos,
                end_pos=end_pos,
                haplotype=HMMCrossoverDetector.STATE_NAMES[current_state],
                confidence=confidence,
                num_snps=num_snps
            )
            blocks.append(block)
        
        return blocks
    
    def save_results(
        self,
        output_dir: str,
        sperm_id: Optional[str] = None
    ) -> None:
        """
        Save detection results to files.
        
        Args:
            output_dir: Output directory path
            sperm_id: Optional specific sperm ID (default: all)
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        sperm_ids = [sperm_id] if sperm_id else list(self.results.keys())
        
        for sid in sperm_ids:
            if sid not in self.results:
                logger.warning(f"No results for sperm {sid}")
                continue
            
            results = self.results[sid]
            
            # Save haplotype blocks
            blocks_file = output_path / f"{sid}_haplotype_blocks.tsv"
            blocks_df = pd.DataFrame([
                {
                    'chromosome': b.chromosome,
                    'start_pos': b.start_pos,
                    'end_pos': b.end_pos,
                    'haplotype': b.haplotype,
                    'confidence': b.confidence,
                    'num_snps': b.num_snps,
                }
                for b in results.haplotype_blocks
            ])
            blocks_df.to_csv(blocks_file, sep='\t', index=False)
            logger.info(f"Saved haplotype blocks to {blocks_file}")
            
            # Save crossovers
            crossovers_file = output_path / f"{sid}_crossovers.tsv"
            crossovers_df = pd.DataFrame([
                {
                    'sperm_id': c.sperm_id,
                    'chromosome': c.chromosome,
                    'position': c.position,
                    'start_pos': c.start_pos,
                    'end_pos': c.end_pos,
                    'from_haplotype': c.from_haplotype,
                    'to_haplotype': c.to_haplotype,
                    'confidence': c.confidence,
                    'resolution_kb': c.resolution_kb,
                }
                for c in results.crossovers
            ])
            crossovers_df.to_csv(crossovers_file, sep='\t', index=False)
            logger.info(f"Saved crossovers to {crossovers_file}")
            
            # Save statistics
            stats_file = output_path / f"{sid}_statistics.tsv"
            stats_df = pd.DataFrame([results.statistics])
            stats_df.to_csv(stats_file, sep='\t', index=False)
            logger.info(f"Saved statistics to {stats_file}")


# ============================================================================
# Command Line Interface
# ============================================================================

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Crossover detection from single sperm sequencing data (Lu et al. 2012)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run complete pipeline
  python crossover_detection.py --input data/sperm_snps.tsv --output results/
  
  # Analyze specific sperm sample
  python crossover_detection.py --input data/sperm_snps.tsv --sperm-id S01 --output results/
  
  # Analyze specific chromosome
  python crossover_detection.py --input data/sperm_snps.tsv --chromosome chr1 --output results/
        """
    )
    
    parser.add_argument(
        '--input', '-i',
        required=True,
        help='Input TSV file with SNP genotypes'
    )
    
    parser.add_argument(
        '--phasing', '-p',
        required=False,
        help='Input TSV file with haplotype phasing information'
    )
    
    parser.add_argument(
        '--output', '-o',
        required=True,
        help='Output directory for results'
    )
    
    parser.add_argument(
        '--sperm-id',
        required=False,
        default=None,
        help='Specific sperm sample ID to analyze (default: all)'
    )
    
    parser.add_argument(
        '--chromosome',
        required=False,
        default=None,
        help='Specific chromosome to analyze (default: all)'
    )
    
    parser.add_argument(
        '--window-size',
        type=int,
        default=10,
        help='SNP window size for rough estimation (default: 10)'
    )
    
    parser.add_argument(
        '--min-fragment-size',
        type=int,
        default=300000,
        help='Minimum fragment size in bp (default: 300000)'
    )
    
    parser.add_argument(
        '--error-rate',
        type=float,
        default=0.05,
        help='SNP calling error rate (default: 0.05)'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose logging'
    )
    
    return parser.parse_args()


def main():
    """Main entry point."""
    args = parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Initialize detector
    detector = CrossoverDetector(
        window_size=args.window_size,
        min_fragment_size=args.min_fragment_size,
        error_rate=args.error_rate
    )
    
    # Load data
    detector.load_snp_data(args.input)
    
    if args.phasing:
        detector.load_phasing_data(args.phasing)
    
    # Get sperm IDs
    if args.sperm_id:
        sperm_ids = [args.sperm_id]
    else:
        sperm_ids = detector.snp_data['sperm_id'].unique() if 'sperm_id' in detector.snp_data.columns else ['default']
    
    # Run detection for each sperm
    for sperm_id in sperm_ids:
        try:
            results = detector.detect_crossovers(
                sperm_id=sperm_id,
                chromosome=args.chromosome
            )
            
            # Print summary
            print(f"\n{'='*60}")
            print(f"Sperm ID: {sperm_id}")
            print(f"{'='*60}")
            print(f"Total crossovers: {results.statistics['total_crossovers']}")
            print(f"Total haplotype blocks: {results.statistics['total_blocks']}")
            print(f"Mean resolution: {results.statistics['mean_resolution_kb']:.2f} kb")
            print(f"Mean confidence: {results.statistics['mean_confidence']:.3f}")
            
        except Exception as e:
            logger.error(f"Error processing sperm {sperm_id}: {e}")
            continue
    
    # Save results
    detector.save_results(args.output, args.sperm_id)
    print(f"\nResults saved to: {args.output}")


if __name__ == '__main__':
    main()
