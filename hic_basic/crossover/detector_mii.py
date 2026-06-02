"""
HMM for sister-chromatid composition inference from SNP read counts.

This module implements a Hidden Markov Model for inferring the local composition
of two sister chromatids along a chromosome in secondary spermatocytes (mII cells).

Reference:
    See issues/mII_crossover/alg.md for algorithm details.

Hidden States:
    - PP: both sister chromatids carry haplotype P
    - PM: one chromatid carries P, the other carries M
    - MM: both chromatids carry haplotype M
    - U: uncertain state

Observations:
    - (ref, alt, other): read counts at each SNP position

Emission Model:
    - Multinomial distribution: (r,a,x) | s ~ Multinomial(n, θ_s)
"""

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def get_default_transition_params() -> Dict[str, float]:
    """
    Get default transition matrix parameters.
    
    These parameters represent a conservative model with:
    - Low transition rates between biological states (PP, PM, MM)
    - Rare entry into U state
    - Moderate exit from U state
    
    Returns:
        Dictionary with keys: alpha, beta, gamma, rho
    """
    return {
        'alpha': 0.01,   # PP↔PM, MM↔PM transition probability
        'beta': 0.02,    # PM→PP, PM→MM transition probability
        'gamma': 0.001,  # Biological state → U transition probability
        'rho': 0.1,      # U → Biological state transition probability
    }

@dataclass
class SisterChromatidBlock:
    """
    Container for a sister chromatid block.

    Attributes:
        chromosome: Chromosome name
        start_pos: Start position (1-based)
        end_pos: End position (1-based)
        state: State type ('PP', 'PM', 'MM', or 'U')
        confidence: Confidence score (0-1)
        num_snps: Number of SNPs in the block
    """
    chromosome: str
    start_pos: int
    end_pos: int
    state: str
    confidence: float
    num_snps: int


@dataclass
class SisterChromatidResults:
    """
    Container for sister chromatid inference results.

    Attributes:
        sample_id: Sample identifier
        chromosome: Chromosome name
        states: List of inferred states for each SNP
        positions: SNP positions
        blocks: List of contiguous blocks with same state
        statistics: Summary statistics
    """
    sample_id: str
    chromosome: str
    states: List[str] = field(default_factory=list)
    positions: List[int] = field(default_factory=list)
    blocks: List[SisterChromatidBlock] = field(default_factory=list)
    statistics: Dict = field(default_factory=dict)


@dataclass
class MiiCrossoverSimulation:
    """
    Container for simulated mII crossover read counts.

    Attributes:
        read_counts: SNP-level read counts and true hidden labels.
        crossover_positions: Per-chromatid crossover positions.
        parameters: Simulation parameters.
    """
    read_counts: pd.DataFrame
    crossover_positions: Dict[str, List[int]]
    parameters: Dict


class SisterChromatidHMM:
    """
    Hidden Markov Model for sister chromatid composition inference.

    States:
        PP (0): both chromatids are haplotype P
        PM (1): one chromatid P, one chromatid M
        MM (2): both chromatids are haplotype M
        U  (3): uncertain

    Observations:
        (ref, alt, other) read counts at each SNP

    Emission Model:
        Multinomial distribution with state-specific probability vectors

    Transition Matrix Structure:
        A = [[1-α-γ,   α,      0,      γ],
             [β/2,     1-β-γ,  β/2,    γ],
             [0,       α,      1-α-γ,  γ],
             [ρ,       ρ,      ρ,      1-3ρ]]

    Attributes:
        alpha: PP<->PM and MM<->PM transition probability
        beta: PM->PP and PM->MM transition probability
        gamma: probability of entering U state
        rho: probability of leaving U state
        emission_params: Dictionary of emission probability vectors θ_s

    Example:
        >>> hmm = SisterChromatidHMM()
        >>> states = hmm.decode(read_counts_df)
    """

    # State indices
    STATE_PP = 0
    STATE_PM = 1
    STATE_MM = 2
    STATE_U = 3

    STATE_NAMES = {0: 'PP', 1: 'PM', 2: 'MM', 3: 'U'}
    STATE_NAMES_REV = {'PP': 0, 'PM': 1, 'MM': 2, 'U': 3}

    def __init__(
        self,
        alpha: float = 0.01,
        beta: float = 0.02,
        gamma: float = 0.001,
        rho: float = 0.1,
        emission_params: Optional[Dict[str, Tuple[float, float, float]]] = None,
        initial_params: Optional[Tuple[float, float, float, float]] = None
    ):
        """
        Initialize the HMM.

        Args:
            alpha: PP<->PM and MM<->PM transition probability
            beta: PM->PP and PM->MM transition probability
            gamma: probability of entering U state
            rho: probability of leaving U state
            emission_params: Emission probability vectors for each state
                Dict mapping state name to (p_ref, p_alt, p_other)
            initial_params: Initial state distribution (π_PP, π_PM, π_MM, π_U)
        """
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.rho = rho

        # Default emission parameters
        if emission_params is None:
            self.emission_params = {
                'PP': (0.95, 0.03, 0.02),
                'PM': (0.49, 0.49, 0.02),
                'MM': (0.03, 0.95, 0.02),
                'U': (0.25, 0.25, 0.50),
            }
        else:
            self.emission_params = emission_params

        # Default initial distribution
        if initial_params is None:
            self.initial_params = (0.3, 0.4, 0.3, 0.0)
        else:
            self.initial_params = initial_params

        # Build matrices
        self.transition_matrix = self._build_transition_matrix()
        self.emission_matrix = self._build_emission_matrix()
        self.initial_vector = np.array(self.initial_params)

    def _build_transition_matrix(self) -> np.ndarray:
        """
        Build the transition probability matrix.

        Returns:
            4x4 transition matrix
        """
        alpha, beta, gamma, rho = self.alpha, self.beta, self.gamma, self.rho

        A = np.array([
            # PP       PM       MM       U
            [1 - alpha - gamma, alpha, 0, gamma],  # From PP
            [beta / 2, 1 - beta - gamma, beta / 2, gamma],  # From PM
            [0, alpha, 1 - alpha - gamma, gamma],  # From MM
            [rho, rho, rho, 1 - 3 * rho],  # From U
        ])

        # Ensure rows sum to 1
        A = A / A.sum(axis=1, keepdims=True)
        return A

    def _build_emission_matrix(self) -> np.ndarray:
        """
        Build the emission probability matrix.

        Returns:
            4x3 emission matrix (states x observation categories)
        """
        B = np.zeros((4, 3))
        for state_idx, state_name in self.STATE_NAMES.items():
            params = self.emission_params.get(state_name, (1/3, 1/3, 1/3))
            B[state_idx, :] = params

        # Ensure rows sum to 1
        B = B / B.sum(axis=1, keepdims=True)
        return B

    def _multinomial_log_emission(
        self,
        ref_count: int,
        alt_count: int,
        other_count: int,
        state_idx: int
    ) -> float:
        """
        Calculate log emission probability for multinomial distribution.

        Args:
            ref_count: Reference allele read count
            alt_count: Alternative allele read count
            other_count: Other allele read count
            state_idx: State index (0-3)

        Returns:
            Log emission probability
        """
        n = ref_count + alt_count + other_count
        if n == 0:
            return -np.inf

        p = self.emission_matrix[state_idx, :]

        # Log multinomial probability (ignoring constant factorial term)
        log_prob = (
            ref_count * np.log(p[0] + 1e-300) +
            alt_count * np.log(p[1] + 1e-300) +
            other_count * np.log(p[2] + 1e-300)
        )

        return log_prob

    def _viterbi(
        self,
        ref_counts: np.ndarray,
        alt_counts: np.ndarray,
        other_counts: np.ndarray
    ) -> Tuple[np.ndarray, float]:
        """
        Viterbi algorithm for most likely state sequence.

        Args:
            ref_counts: Array of reference allele counts
            alt_counts: Array of alternative allele counts
            other_counts: Array of other allele counts

        Returns:
            Tuple of (best_path, log_probability)
        """
        n_snps = len(ref_counts)
        n_states = 4

        # Initialize Viterbi matrices
        delta = np.full((n_snps, n_states), -np.inf)
        psi = np.zeros((n_snps, n_states), dtype=int)

        # Initialization (t=0)
        for j in range(n_states):
            log_emit = self._multinomial_log_emission(
                ref_counts[0], alt_counts[0], other_counts[0], j
            )
            delta[0, j] = np.log(self.initial_vector[j] + 1e-300) + log_emit

        # Recursion
        for t in range(1, n_snps):
            for j in range(n_states):
                # Find best previous state
                probs = delta[t-1, :] + np.log(self.transition_matrix[:, j] + 1e-300)
                psi[t, j] = np.argmax(probs)
                best_prev = probs[psi[t, j]]

                # Add emission probability
                log_emit = self._multinomial_log_emission(
                    ref_counts[t], alt_counts[t], other_counts[t], j
                )
                delta[t, j] = best_prev + log_emit

        # Termination
        best_path = np.zeros(n_snps, dtype=int)
        best_path[-1] = np.argmax(delta[-1, :])
        log_prob = np.max(delta[-1, :])

        # Traceback
        for t in range(n_snps - 2, -1, -1):
            best_path[t] = psi[t + 1, best_path[t + 1]]

        return best_path, log_prob

    def decode(
        self,
        read_counts: pd.DataFrame
    ) -> Tuple[List[str], float]:
        """
        Decode the most likely state sequence from read counts.

        Args:
            read_counts: DataFrame with columns 'ref', 'alt', 'other'
                        or 'r', 'a', 'x'

        Returns:
            Tuple of (state_names, log_probability)
        """
        # Extract counts
        if 'ref' in read_counts.columns:
            ref_counts = read_counts['ref'].values
            alt_counts = read_counts['alt'].values
            other_counts = read_counts['other'].values
        elif 'r' in read_counts.columns:
            ref_counts = read_counts['r'].values
            alt_counts = read_counts['a'].values
            other_counts = read_counts['x'].values
        else:
            raise ValueError("DataFrame must have 'ref/alt/other' or 'r/a/x' columns")

        # Run Viterbi
        best_path, log_prob = self._viterbi(ref_counts, alt_counts, other_counts)

        # Convert to state names
        state_names = [self.STATE_NAMES[s] for s in best_path]

        return state_names, log_prob

    def decode_with_blocks(
        self,
        read_counts: pd.DataFrame,
        chromosome: str,
        positions: np.ndarray,
        sample_id: str = "unknown"
    ) -> SisterChromatidResults:
        """
        Decode states and group into contiguous blocks.

        Args:
            read_counts: DataFrame with read counts
            chromosome: Chromosome name
            positions: SNP positions
            sample_id: Sample identifier

        Returns:
            SisterChromatidResults object
        """
        states, log_prob = self.decode(read_counts)

        # Group into blocks
        blocks = self._states_to_blocks(
            chromosome=chromosome,
            positions=positions,
            states=states
        )

        # Calculate statistics
        state_counts = pd.Series(states).value_counts()
        statistics = {
            'total_snps': len(states),
            'log_probability': log_prob,
            'state_counts': state_counts.to_dict(),
            'n_blocks': len(blocks),
        }

        results = SisterChromatidResults(
            sample_id=sample_id,
            chromosome=chromosome,
            states=states,
            positions=positions.tolist(),
            blocks=blocks,
            statistics=statistics
        )

        return results

    def _states_to_blocks(
        self,
        chromosome: str,
        positions: np.ndarray,
        states: List[str]
    ) -> List[SisterChromatidBlock]:
        """
        Convert state sequence to contiguous blocks.

        Args:
            chromosome: Chromosome name
            positions: SNP positions
            states: List of state names

        Returns:
            List of SisterChromatidBlock objects
        """
        blocks = []
        n_snps = len(positions)

        if n_snps == 0:
            return blocks

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
            confidence = sum(s == current_state for s in block_states) / len(block_states)

            block = SisterChromatidBlock(
                chromosome=chromosome,
                start_pos=start_pos,
                end_pos=end_pos,
                state=current_state,
                confidence=confidence,
                num_snps=num_snps
            )
            blocks.append(block)

        return blocks


def estimate_emissions_from_counts(
    read_counts: pd.DataFrame,
    rough_labels: pd.Series
) -> Dict[str, Tuple[float, float, float]]:
    """
    Estimate emission probabilities from read counts and rough labels.

    Args:
        read_counts: DataFrame with columns 'ref', 'alt', 'other'
        rough_labels: Series with rough state labels ('PP', 'PM', 'MM', 'U')

    Returns:
        Dictionary mapping state name to (p_ref, p_alt, p_other)
    """
    emission_params = {}

    for state in ['PP', 'PM', 'MM', 'U']:
        mask = rough_labels == state
        if mask.sum() == 0:
            # Use default if no SNPs assigned to this state
            if state == 'PP':
                emission_params[state] = (0.95, 0.03, 0.02)
            elif state == 'PM':
                emission_params[state] = (0.49, 0.49, 0.02)
            elif state == 'MM':
                emission_params[state] = (0.03, 0.95, 0.02)
            else:
                emission_params[state] = (0.25, 0.25, 0.50)
            continue

        state_counts = read_counts[mask]
        total_ref = state_counts['ref'].sum()
        total_alt = state_counts['alt'].sum()
        total_other = state_counts['other'].sum()
        total = total_ref + total_alt + total_other

        if total > 0:
            emission_params[state] = (
                total_ref / total,
                total_alt / total,
                total_other / total
            )
        else:
            emission_params[state] = (1/3, 1/3, 1/3)

    logger.info(f"Estimated emission parameters: {emission_params}")
    return emission_params


def estimate_transitions_rough(
    rough_labels: pd.Series,
    pseudocount: float = 1.0
) -> Tuple[float, float, float, float]:
    """
    Estimate transition probabilities from rough labels.

    Self transitions are included in the row totals so the returned parameters
    represent per-adjacent-observation probabilities, not conditional
    probabilities after a non-self jump has already occurred. The returned rho
    is the per-target U->biological-state probability used by the HMM; the
    total U exit probability is 3*rho.

    Args:
        rough_labels: Series with rough state labels
        pseudocount: Pseudocount for smoothing (deprecated, kept for API compatibility)

    Returns:
        Tuple of (alpha, beta, gamma, rho)
    """
    states = ['PP', 'PM', 'MM', 'U']
    defaults = get_default_transition_params()

    # Count the full adjacent transition table, including self transitions.
    # The HMM transition parameters are per adjacent observation; excluding
    # self transitions estimates conditional jump composition instead and can
    # severely deflate state persistence.
    transitions = {(src, dst): 0 for src in states for dst in states}
    for i in range(1, len(rough_labels)):
        prev_state = rough_labels.iloc[i - 1]
        curr_state = rough_labels.iloc[i]
        key = (prev_state, curr_state)
        if key in transitions:
            transitions[key] += 1

    total_transitions = sum(transitions.values())
    if total_transitions == 0:
        logger.info(
            "Using default transition parameters: "
            f"alpha={defaults['alpha']:.4f}, beta={defaults['beta']:.4f}, "
            f"gamma={defaults['gamma']:.4f}, rho={defaults['rho']:.4f}"
        )
        return defaults['alpha'], defaults['beta'], defaults['gamma'], defaults['rho']

    row_total = {
        state: sum(transitions[(state, dst)] for dst in states)
        for state in states
    }

    # alpha: pure state (PP/MM) -> mixed PM per-step probability.
    # Direct PP<->MM rough jumps are treated as pure-state switching evidence,
    # since the structured HMM can only express that path through PM.
    pure_total = row_total['PP'] + row_total['MM']
    alpha_num = (
        transitions[('PP', 'PM')] + transitions[('MM', 'PM')] +
        transitions[('PP', 'MM')] + transitions[('MM', 'PP')]
    )
    alpha = alpha_num / pure_total if pure_total > 0 else defaults['alpha']

    # beta: mixed PM -> pure state per-step probability.
    pm_total = row_total['PM']
    beta_num = transitions[('PM', 'PP')] + transitions[('PM', 'MM')]
    beta = beta_num / pm_total if pm_total > 0 else defaults['beta']

    # gamma: biological states -> U per-step probability.
    bio_total = row_total['PP'] + row_total['PM'] + row_total['MM']
    u_enter = transitions[('PP', 'U')] + transitions[('PM', 'U')] + transitions[('MM', 'U')]
    gamma = u_enter / bio_total if bio_total > 0 else defaults['gamma']

    # rho is the per-target U->biological probability. The total probability
    # of leaving U is 3*rho in the transition matrix.
    u_total = row_total['U']
    u_leave = transitions[('U', 'PP')] + transitions[('U', 'PM')] + transitions[('U', 'MM')]
    rho = (u_leave / u_total / 3.0) if u_total > 0 else defaults['rho']
    rho = min(rho, 0.32)  # Clamp to keep 1-3*rho positive.

    # Ensure transition matrix elements are non-negative
    # For PP and MM rows: alpha + gamma <= 1
    # For PM row: beta + gamma <= 1
    if alpha + gamma > 0.99:
        # Normalize alpha and gamma proportionally
        scale = 0.99 / (alpha + gamma)
        alpha *= scale
        gamma *= scale
    
    if beta + gamma > 0.99:
        # Normalize beta and gamma proportionally
        scale = 0.99 / (beta + gamma)
        beta *= scale
        gamma *= scale

    logger.info(f"Estimated transition parameters: alpha={alpha:.4f}, beta={beta:.4f}, gamma={gamma:.4f}, rho={rho:.4f}")
    return alpha, beta, gamma, rho


def rough_classification(
    read_counts: pd.DataFrame,
    depth_threshold: int = 2,
    other_threshold: float = 0.4,
    balance_threshold: float = 0.25
) -> pd.Series:
    """
    Perform rough classification of SNPs based on read counts.

    Args:
        read_counts: DataFrame with columns 'ref', 'alt', 'other'
        depth_threshold: Minimum total depth for confident call (default: 2)
        other_threshold: Maximum other fraction for confident call (default: 0.4)
        balance_threshold: Threshold for classifying as balanced (PM) (default: 0.25)
                          PP: ref_frac > 0.75, MM: alt_frac > 0.75, PM: |ref-alt| < 0.25.
                          Ratio-intermediate observations are labeled AMBIG, not U.

    Returns:
        Series with rough state labels. U is reserved for low-depth or high-other
        observations; AMBIG marks ratio-intermediate observations.
    """
    labels = []

    for _, row in read_counts.iterrows():
        ref = row['ref']
        alt = row['alt']
        other = row['other']
        total = ref + alt + other

        if total < depth_threshold:
            labels.append('U')
            continue

        other_frac = other / total
        if other_frac > other_threshold:
            labels.append('U')
            continue

        ref_frac = ref / total
        alt_frac = alt / total

        # Classify based on allele balance
        # PP: ref_frac > 0.75, MM: alt_frac > 0.75, PM: |ref_frac - alt_frac| < 0.25
        if ref_frac > 1 - balance_threshold:  # ref_frac > 0.75
            labels.append('PP')
        elif alt_frac > 1 - balance_threshold:  # alt_frac > 0.75
            labels.append('MM')
        elif abs(ref_frac - alt_frac) < balance_threshold and ref_frac > 0.2 and alt_frac > 0.2:
            labels.append('PM')
        else:
            labels.append('AMBIG')

    return pd.Series(labels)


def simulate_mii_crossover_counts(
    chromosome: str = "chr1",
    chromosome_length: int = 100_000_000,
    snp_spacing: int = 1_000,
    crossover_rate_per_mb: float = 0.02,
    mean_depth: float = 5.0,
    sequencing_error_rate: float = 0.005,
    depth_cv: float = 0.55,
    local_depth_window_snps: int = 100,
    random_state: Optional[Union[int, np.random.Generator]] = None,
) -> MiiCrossoverSimulation:
    """
    Simulate SNP read counts for mII sister-chromatid crossover inference.

    The simulated mII signal is the aggregate of two independent chromatid
    haplotype paths. Each chromatid starts as either P or M, switches haplotype
    at crossover events, and contributes equally to the observed read pool.
    Read counts are generated with WGS-like local depth variation, Poisson
    sampling, and random sequencing errors.

    Args:
        chromosome: Chromosome name assigned to simulated SNPs.
        chromosome_length: Chromosome length in bp.
        snp_spacing: Distance between adjacent SNPs in bp.
        crossover_rate_per_mb: Per-chromatid crossover probability per Mb.
        mean_depth: Mean total read depth per SNP.
        sequencing_error_rate: Per-read probability of observing the wrong base.
        depth_cv: Coefficient of variation for local depth multipliers.
        local_depth_window_snps: Number of SNPs sharing a local depth bias.
        random_state: Optional seed or numpy Generator for reproducibility.

    Returns:
        MiiCrossoverSimulation with read-count DataFrame, crossover positions,
        and parameter metadata.
    """
    if chromosome_length <= 0:
        raise ValueError("chromosome_length must be positive")
    if snp_spacing <= 0:
        raise ValueError("snp_spacing must be positive")
    if crossover_rate_per_mb < 0:
        raise ValueError("crossover_rate_per_mb must be non-negative")
    if mean_depth < 0:
        raise ValueError("mean_depth must be non-negative")
    if not 0 <= sequencing_error_rate < 1:
        raise ValueError("sequencing_error_rate must be in [0, 1)")
    if depth_cv < 0:
        raise ValueError("depth_cv must be non-negative")
    if local_depth_window_snps <= 0:
        raise ValueError("local_depth_window_snps must be positive")

    if isinstance(random_state, np.random.Generator):
        rng = random_state
    else:
        rng = np.random.default_rng(random_state)

    positions = np.arange(1, chromosome_length + 1, snp_spacing, dtype=np.int64)
    n_snps = len(positions)
    if n_snps == 0:
        raise ValueError("No SNPs generated; check chromosome_length and snp_spacing")

    interval_mb = snp_spacing / 1_000_000
    crossover_prob = min(crossover_rate_per_mb * interval_mb, 1.0)

    hap_paths = []
    crossover_positions: Dict[str, List[int]] = {}
    for chromatid_idx in range(2):
        start_hap = int(rng.integers(0, 2))  # 1 means P/ref, 0 means M/alt.
        switches = rng.random(max(n_snps - 1, 0)) < crossover_prob
        path = np.empty(n_snps, dtype=np.int8)
        path[0] = start_hap
        if n_snps > 1:
            path[1:] = (start_hap + np.cumsum(switches)) % 2
        hap_paths.append(path)
        crossover_positions[f"chromatid_{chromatid_idx + 1}"] = positions[1:][switches].astype(int).tolist()

    hap1, hap2 = hap_paths
    n_ref_copies = hap1 + hap2
    state_by_ref_copies = np.array(["MM", "PM", "PP"], dtype=object)
    true_states = state_by_ref_copies[n_ref_copies]

    n_windows = int(np.ceil(n_snps / local_depth_window_snps))
    if depth_cv == 0:
        window_bias = np.ones(n_windows)
    else:
        gamma_shape = 1.0 / (depth_cv ** 2)
        gamma_scale = depth_cv ** 2
        window_bias = rng.gamma(shape=gamma_shape, scale=gamma_scale, size=n_windows)
    local_bias = np.repeat(window_bias, local_depth_window_snps)[:n_snps]
    per_snp_bias = rng.lognormal(mean=-0.5 * 0.15 ** 2, sigma=0.15, size=n_snps)
    expected_depth = mean_depth * local_bias * per_snp_bias
    depth = rng.poisson(expected_depth)

    err = sequencing_error_rate
    ref_copy_fraction = n_ref_copies / 2.0
    alt_copy_fraction = 1.0 - ref_copy_fraction
    p_ref = ref_copy_fraction * (1.0 - err) + alt_copy_fraction * (err / 3.0)
    p_alt = alt_copy_fraction * (1.0 - err) + ref_copy_fraction * (err / 3.0)
    p_other = np.full(n_snps, 2.0 * err / 3.0)
    probs = np.vstack([p_ref, p_alt, p_other]).T
    probs = probs / probs.sum(axis=1, keepdims=True)

    counts = np.array([rng.multinomial(int(d), p) for d, p in zip(depth, probs)])
    read_counts = pd.DataFrame({
        "chromosome": chromosome,
        "position": positions,
        "ref": counts[:, 0],
        "alt": counts[:, 1],
        "other": counts[:, 2],
        "true_state": true_states,
        "hap1": np.where(hap1 == 1, "P", "M"),
        "hap2": np.where(hap2 == 1, "P", "M"),
        "expected_depth": expected_depth,
        "depth_bias": local_bias,
    })

    parameters = {
        "chromosome": chromosome,
        "chromosome_length": chromosome_length,
        "snp_spacing": snp_spacing,
        "n_snps": n_snps,
        "crossover_rate_per_mb": crossover_rate_per_mb,
        "mean_depth": mean_depth,
        "sequencing_error_rate": sequencing_error_rate,
        "depth_cv": depth_cv,
        "local_depth_window_snps": local_depth_window_snps,
    }

    return MiiCrossoverSimulation(
        read_counts=read_counts,
        crossover_positions=crossover_positions,
        parameters=parameters,
    )
