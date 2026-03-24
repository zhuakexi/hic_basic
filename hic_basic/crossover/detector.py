"""
HMM-based crossover detection algorithm.

This module implements the core Hidden Markov Model for detecting
meiotic crossover events based on Lu et al. (2012) Science method.
"""

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class HaplotypeBlock:
    """
    Container for a haplotype block.

    Attributes:
        chromosome: Chromosome name (e.g., 'chr1', '1')
        start_pos: Start position (1-based)
        end_pos: End position (1-based)
        haplotype: Haplotype type ('P' for paternal, 'M' for maternal, 'U' for undetermined)
        confidence: Confidence score (0-1)
        num_snps: Number of SNPs in the block
    """
    chromosome: str
    start_pos: int
    end_pos: int
    haplotype: str
    confidence: float
    num_snps: int


@dataclass
class Crossover:
    """
    Container for a crossover event.

    Attributes:
        sperm_id: Sample ID
        chromosome: Chromosome name
        position: Crossover position (1-based)
        start_pos: Start position of the transition region
        end_pos: End position of the transition region
        from_haplotype: Haplotype before crossover ('P' or 'M')
        to_haplotype: Haplotype after crossover ('P' or 'M')
        confidence: Confidence score (0-1)
        resolution_kb: Resolution in kilobase pairs
    """
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
    """
    Container for crossover detection results.

    Attributes:
        sperm_id: Sample ID
        haplotype_blocks: List of detected haplotype blocks
        crossovers: List of detected crossover events
        statistics: Summary statistics
    """
    sperm_id: str
    haplotype_blocks: List[HaplotypeBlock] = field(default_factory=list)
    crossovers: List[Crossover] = field(default_factory=list)
    statistics: Dict = field(default_factory=dict)


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

    Attributes:
        window_size: Number of SNPs per window for rough estimation
        min_fragment_size: Minimum fragment size in bp (filter small fragments)
        error_rate: SNP calling error rate for emission probabilities
        crossover_rate: Prior crossover rate for transition probabilities
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
        self.transition_matrix: Optional[np.ndarray] = None
        self.emission_matrix: Optional[np.ndarray] = None
        self._init_transition_matrix()
        self._init_emission_matrix()

    def _init_transition_matrix(self) -> None:
        """
        Initialize transition probability matrix.

        Based on rough crossover frequency from Step 1.
        P->M and M->P transitions represent crossovers.
        """
        p_crossover = self.crossover_rate / 2
        p_stay = 1 - p_crossover
        p_undetermined = 0.01

        self.transition_matrix = np.array([
            [p_stay, p_crossover, p_undetermined],
            [p_crossover, p_stay, p_undetermined],
            [p_undetermined, p_undetermined, 1 - 2 * p_undetermined],
        ])

        self.transition_matrix = self.transition_matrix / self.transition_matrix.sum(axis=1, keepdims=True)

    def _init_emission_matrix(self) -> None:
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
            [1 - e, e, 0.01],
            [e, 1 - e, 0.01],
            [0.01, 0.01, 0.98],
        ])

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
            positions: SNP positions (sorted)
            haplotype_alleles: Array of 0 (paternal) or 1 (maternal) for each SNP

        Returns:
            List of rough crossover positions
        """
        n_snps = len(positions)
        if n_snps < self.window_size * 2:
            return []

        n_windows = n_snps - self.window_size + 1
        paternal_ratios = np.zeros(n_windows)
        window_positions = np.zeros(n_windows)

        for i in range(n_windows):
            window = haplotype_alleles[i:i + self.window_size]
            valid_snps = window[window >= 0]
            if len(valid_snps) > 0:
                paternal_ratios[i] = np.sum(valid_snps == 0) / len(valid_snps)
            else:
                paternal_ratios[i] = 0.5
            window_positions[i] = np.mean(positions[i:i + self.window_size])

        crossover_positions = []
        threshold = 0.5

        for i in range(1, len(paternal_ratios)):
            ratio_change = abs(paternal_ratios[i] - paternal_ratios[i-1])
            if ratio_change > threshold:
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

        delta = np.zeros((n_obs, n_states))
        psi = np.zeros((n_obs, n_states), dtype=int)

        delta[0, :] = 1.0 / n_states * self.emission_matrix[:, observations[0]]

        for t in range(1, n_obs):
            for j in range(n_states):
                probs = delta[t-1, :] * self.transition_matrix[:, j]
                psi[t, j] = np.argmax(probs)
                delta[t, j] = probs[psi[t, j]] * self.emission_matrix[j, observations[t]]

            delta[t, :] /= np.sum(delta[t, :]) + 1e-300

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
            if sperm_gt in ['N', '.', '-', ''] or len(sperm_gt) == 0:
                observations.append(2)
            elif sperm_gt == h1:
                observations.append(0)
            elif sperm_gt == h2:
                observations.append(1)
            else:
                observations.append(2)

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
            current_state = states[i]
            start_idx = i

            while i < n_snps and states[i] == current_state:
                i += 1

            end_idx = i - 1

            if end_idx > start_idx:
                fragment_size = positions[end_idx] - positions[start_idx]

                if fragment_size < self.min_fragment_size and current_state == self.STATE_UNDETERMINED:
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
        haplotype_alleles = self._encode_observations(sperm_genotypes, haplotype_1, haplotype_2)
        rough_crossovers = self._rough_estimation(positions, haplotype_alleles)

        if len(positions) > 0 and len(rough_crossovers) > 0:
            genome_span = positions[-1] - positions[0]
            if genome_span > 0:
                self.crossover_rate = len(rough_crossovers) / (genome_span / 1e6)
                self._init_transition_matrix()

        observations = haplotype_alleles
        states, log_prob = self._viterbi(observations)
        filtered_states = self._filter_small_fragments(positions, states)

        crossover_positions = []
        for i in range(1, len(filtered_states)):
            if filtered_states[i] != filtered_states[i-1]:
                if filtered_states[i] in [self.STATE_PATERNAL, self.STATE_MATERNAL] and \
                   filtered_states[i-1] in [self.STATE_PATERNAL, self.STATE_MATERNAL]:
                    crossover_pos = int((positions[i] + positions[i-1]) / 2)
                    crossover_positions.append(crossover_pos)

        logger.debug(f"HMM detected {len(crossover_positions)} crossovers")
        return filtered_states, crossover_positions


class CrossoverDetector:
    """
    Main class for crossover detection from single sperm sequencing data.

    This class provides the main interface for running crossover detection
    on SNP genotype data with phasing information.

    Attributes:
        window_size: SNP window size for rough estimation
        min_fragment_size: Minimum fragment size to keep (bp)
        error_rate: SNP calling error rate
        snp_data: Loaded SNP genotype DataFrame
        phasing_data: Loaded phasing information DataFrame
        results: Dictionary of results by sperm ID

    Example:
        >>> detector = CrossoverDetector()
        >>> detector.load_snp_data("snps.tsv")
        >>> detector.load_phasing_data("phasing.tsv")
        >>> results = detector.detect_crossovers("S01")
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
                      Expected columns: chromosome, position, snp_id, sperm_genotype, diploid_genotype
        """
        logger.info(f"Loading SNP data from {snp_file}")
        self.snp_data = pd.read_csv(snp_file, sep='\t')
        logger.info(f"Loaded {len(self.snp_data)} SNPs")

    def load_phasing_data(self, phasing_file: str) -> None:
        """
        Load phasing information.

        Args:
            phasing_file: Path to TSV file with haplotype phasing
                          Expected columns: chromosome, position, snp_id, haplotype_1, haplotype_2
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

        Raises:
            ValueError: If SNP data not loaded
        """
        logger.info(f"Detecting crossovers for sperm {sperm_id}")

        if self.snp_data is None:
            raise ValueError("SNP data not loaded. Call load_snp_data() first.")

        results = CrossoverResults(sperm_id=sperm_id)
        chromosomes = [chromosome] if chromosome else self.snp_data['chromosome'].unique()

        total_crossovers = 0
        total_blocks = 0

        for chrom in chromosomes:
            chrom_data = self.snp_data[self.snp_data['chromosome'] == chrom].copy()

            if len(chrom_data) == 0:
                continue

            if self.phasing_data is not None:
                chrom_data = chrom_data.merge(
                    self.phasing_data,
                    on=['chromosome', 'position', 'snp_id'],
                    how='left'
                )

            positions = chrom_data['position'].values
            sperm_genotypes = chrom_data['sperm_genotype'].tolist()

            if 'haplotype_1' in chrom_data.columns:
                haplotype_1 = chrom_data['haplotype_1'].fillna('N').tolist()
                haplotype_2 = chrom_data['haplotype_2'].fillna('N').tolist()
            else:
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

            states, crossover_positions = self.hmm_detector.detect(
                positions=positions,
                sperm_genotypes=sperm_genotypes,
                haplotype_1=haplotype_1,
                haplotype_2=haplotype_2
            )

            blocks = self._states_to_blocks(
                chromosome=chrom,
                positions=positions,
                states=states
            )
            results.haplotype_blocks.extend(blocks)
            total_blocks += len(blocks)

            for i, pos in enumerate(crossover_positions):
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

            while i < n_snps and states[i] == current_state:
                i += 1

            end_idx = i - 1
            end_pos = positions[end_idx]
            num_snps = end_idx - start_idx + 1

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
        from pathlib import Path
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        sperm_ids = [sperm_id] if sperm_id else list(self.results.keys())

        for sid in sperm_ids:
            if sid not in self.results:
                logger.warning(f"No results for sperm {sid}")
                continue

            results = self.results[sid]

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

            stats_file = output_path / f"{sid}_statistics.tsv"
            stats_df = pd.DataFrame([results.statistics])
            stats_df.to_csv(stats_file, sep='\t', index=False)
            logger.info(f"Saved statistics to {stats_file}")
