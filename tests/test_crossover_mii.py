"""
Tests for sister chromatid crossover detection module (mII cells).
"""

import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from hic_basic.crossover.detector_mii import (
    SisterChromatidHMM,
    SisterChromatidBlock,
    SisterChromatidResults,
    estimate_emissions_from_counts,
    estimate_transitions_rough,
    rough_classification,
)


class TestSisterChromatidHMM:
    """Tests for SisterChromatidHMM class."""

    def test_init(self):
        """Test initialization with default parameters."""
        hmm = SisterChromatidHMM()
        assert hmm.alpha == 0.01
        assert hmm.beta == 0.02
        assert hmm.gamma == 0.001
        assert hmm.rho == 0.1
        assert hmm.transition_matrix.shape == (4, 4)
        assert hmm.emission_matrix.shape == (4, 3)

    def test_transition_matrix_structure(self):
        """Test transition matrix has correct structure."""
        hmm = SisterChromatidHMM()
        A = hmm.transition_matrix

        # Check rows sum to 1
        assert np.allclose(A.sum(axis=1), 1.0)

        # Check PP->MM and MM->PP are 0 (or very small)
        assert A[0, 2] < 0.01  # PP->MM
        assert A[2, 0] < 0.01  # MM->PP

        # Check symmetry in PM transitions
        assert np.isclose(A[1, 0], A[1, 2])  # PM->PP ≈ PM->MM

    def test_emission_matrix_defaults(self):
        """Test default emission probabilities."""
        hmm = SisterChromatidHMM()
        B = hmm.emission_matrix

        # PP should have high ref probability
        assert B[0, 0] > 0.9  # PP: ref

        # MM should have high alt probability
        assert B[2, 1] > 0.9  # MM: alt

        # PM should have balanced ref and alt
        assert np.isclose(B[1, 0], B[1, 1], atol=0.05)  # PM: balanced

        # U should have high other probability
        assert B[3, 2] > 0.4  # U: other

    def test_multinomial_log_emission(self):
        """Test multinomial log emission probability."""
        hmm = SisterChromatidHMM()

        # High ref counts should favor PP state
        log_prob_pp = hmm._multinomial_log_emission(20, 1, 0, 0)  # PP state
        log_prob_mm = hmm._multinomial_log_emission(20, 1, 0, 2)  # MM state
        assert log_prob_pp > log_prob_mm

        # Balanced counts should favor PM state
        log_prob_pm = hmm._multinomial_log_emission(10, 10, 0, 1)  # PM state
        log_prob_pp_balanced = hmm._multinomial_log_emission(10, 10, 0, 0)  # PP state
        assert log_prob_pm > log_prob_pp_balanced

        # High other counts should favor U state
        log_prob_u = hmm._multinomial_log_emission(2, 1, 9, 3)  # U state
        log_prob_pp_other = hmm._multinomial_log_emission(2, 1, 9, 0)  # PP state
        assert log_prob_u > log_prob_pp_other

    def test_viterbi_simple(self):
        """Test Viterbi algorithm with simple pattern."""
        hmm = SisterChromatidHMM()

        # Create clear pattern: PP -> PM -> MM
        ref_counts = np.array([20, 20, 20, 10, 10, 1, 1, 1])
        alt_counts = np.array([1, 1, 1, 10, 10, 20, 20, 20])
        other_counts = np.array([0, 0, 0, 0, 0, 0, 0, 0])

        states, log_prob = hmm._viterbi(ref_counts, alt_counts, other_counts)

        assert len(states) == 8
        assert np.isfinite(log_prob)

        # First SNPs should be PP, last SNPs should be MM
        assert states[0] == SisterChromatidHMM.STATE_PP
        assert states[-1] == SisterChromatidHMM.STATE_MM

    def test_decode(self):
        """Test decode method with DataFrame input."""
        hmm = SisterChromatidHMM()

        read_counts = pd.DataFrame({
            'ref': [20, 20, 10, 1, 1],
            'alt': [1, 1, 10, 20, 20],
            'other': [0, 0, 0, 0, 0],
        })

        states, log_prob = hmm.decode(read_counts)

        assert len(states) == 5
        assert all(s in ['PP', 'PM', 'MM', 'U'] for s in states)

    def test_decode_with_blocks(self):
        """Test decode_with_blocks method."""
        hmm = SisterChromatidHMM()

        n_snps = 100
        positions = np.sort(np.random.choice(1000000, n_snps, replace=False))

        # Create data with clear blocks
        ref_counts = []
        alt_counts = []
        for i in range(n_snps):
            if i < 30:
                ref_counts.append(20)
                alt_counts.append(1)
            elif i < 70:
                ref_counts.append(10)
                alt_counts.append(10)
            else:
                ref_counts.append(1)
                alt_counts.append(20)

        read_counts = pd.DataFrame({
            'ref': ref_counts,
            'alt': alt_counts,
            'other': [0] * n_snps,
        })

        results = hmm.decode_with_blocks(
            read_counts=read_counts,
            chromosome='chr1',
            positions=positions,
            sample_id='test'
        )

        assert isinstance(results, SisterChromatidResults)
        assert results.sample_id == 'test'
        assert results.chromosome == 'chr1'
        assert len(results.states) == n_snps
        assert len(results.blocks) >= 1
        assert 'total_snps' in results.statistics

    def test_states_to_blocks(self):
        """Test converting states to blocks."""
        hmm = SisterChromatidHMM()

        positions = np.array([1000, 2000, 3000, 4000, 5000, 6000])
        states = ['PP', 'PP', 'PM', 'PM', 'PM', 'MM']

        blocks = hmm._states_to_blocks('chr1', positions, states)

        assert len(blocks) == 3
        assert blocks[0].state == 'PP'
        assert blocks[0].num_snps == 2
        assert blocks[1].state == 'PM'
        assert blocks[1].num_snps == 3
        assert blocks[2].state == 'MM'
        assert blocks[2].num_snps == 1


class TestRoughClassification:
    """Tests for rough_classification function."""

    def test_basic_classification(self):
        """Test basic rough classification."""
        read_counts = pd.DataFrame({
            'ref': [20, 1, 10, 2],
            'alt': [1, 20, 10, 1],
            'other': [0, 0, 0, 8],
        })

        labels = rough_classification(read_counts)

        assert len(labels) == 4
        assert labels.iloc[0] == 'PP'  # High ref
        assert labels.iloc[1] == 'MM'  # High alt
        assert labels.iloc[3] == 'U'   # High other

    def test_low_depth(self):
        """Test classification with low depth."""
        read_counts = pd.DataFrame({
            'ref': [1, 1],
            'alt': [1, 0],
            'other': [0, 1],
        })

        labels = rough_classification(read_counts, depth_threshold=5)

        # All should be U due to low depth
        assert all(l == 'U' for l in labels)

    def test_balance_threshold(self):
        """Test effect of balance threshold."""
        read_counts = pd.DataFrame({
            'ref': [15, 12, 10],
            'alt': [5, 8, 10],
            'other': [0, 0, 0],
        })

        # Strict balance threshold
        labels_strict = rough_classification(read_counts, balance_threshold=0.2)

        # Loose balance threshold
        labels_loose = rough_classification(read_counts, balance_threshold=0.4)

        # More SNPs should be classified as PM with loose threshold
        n_pm_strict = sum(l == 'PM' for l in labels_strict)
        n_pm_loose = sum(l == 'PM' for l in labels_loose)
        assert n_pm_loose >= n_pm_strict


class TestParameterEstimation:
    """Tests for parameter estimation functions."""

    def test_estimate_emissions_from_counts(self):
        """Test emission parameter estimation."""
        read_counts = pd.DataFrame({
            'ref': [20, 20, 1, 1, 10],
            'alt': [1, 1, 20, 20, 10],
            'other': [0, 0, 0, 0, 0],
        })

        rough_labels = pd.Series(['PP', 'PP', 'MM', 'MM', 'PM'])

        emission_params = estimate_emissions_from_counts(read_counts, rough_labels)

        assert 'PP' in emission_params
        assert 'MM' in emission_params
        assert 'PM' in emission_params

        # PP should have high ref probability
        assert emission_params['PP'][0] > 0.9

        # MM should have high alt probability
        assert emission_params['MM'][1] > 0.9

    def test_estimate_transitions_rough(self):
        """Test transition parameter estimation."""
        # Create pattern: PP -> PP -> PM -> PM -> MM
        rough_labels = pd.Series(['PP', 'PP', 'PM', 'PM', 'MM'])

        alpha, beta, gamma, rho = estimate_transitions_rough(rough_labels)

        # Parameters should be in valid range [0, 1]
        # Note: gamma can be 0 if no transitions to U state
        assert 0 <= alpha <= 1
        assert 0 <= beta <= 1
        assert 0 <= gamma <= 1
        assert 0 <= rho <= 0.32  # rho is clamped to max 0.32

    def test_estimate_transitions_with_pseudocount(self):
        """Test transition estimation with no transitions (uses default values)."""
        # Single state - should return default values
        rough_labels = pd.Series(['PP', 'PP', 'PP'])

        alpha, beta, gamma, rho = estimate_transitions_rough(
            rough_labels, pseudocount=1.0  # pseudocount is deprecated but accepted
        )

        n = len(rough_labels)
        # Default values: alpha=1/n, beta=1/n, gamma=0.01/n, rho=min(100/n, 0.32)
        assert abs(alpha - 1.0 / n) < 1e-10
        assert abs(beta - 1.0 / n) < 1e-10
        assert abs(gamma - 0.01 / n) < 1e-10
        assert abs(rho - min(100.0 / n, 0.32)) < 1e-10


class TestDataClasses:
    """Tests for data classes."""

    def test_sister_chromatid_block(self):
        """Test SisterChromatidBlock dataclass."""
        block = SisterChromatidBlock(
            chromosome='chr1',
            start_pos=1000,
            end_pos=50000,
            state='PP',
            confidence=0.95,
            num_snps=50
        )

        assert block.chromosome == 'chr1'
        assert block.state == 'PP'
        assert 0 <= block.confidence <= 1

    def test_sister_chromatid_results(self):
        """Test SisterChromatidResults dataclass."""
        results = SisterChromatidResults(
            sample_id='S01',
            chromosome='chr1'
        )

        assert results.sample_id == 'S01'
        assert results.chromosome == 'chr1'
        assert len(results.states) == 0
        assert isinstance(results.statistics, dict)


class TestIntegration:
    """Integration tests for the full HMM pipeline."""

    def test_full_pipeline_synthetic_data(self):
        """Test full pipeline with synthetic data."""
        # Generate synthetic data with known state pattern
        np.random.seed(42)
        n_snps = 200
        positions = np.sort(np.random.choice(1000000, n_snps, replace=False))

        # True states: PP (0-50), PM (50-100), MM (100-150), PP (150-200)
        true_states = ['PP'] * 50 + ['PM'] * 50 + ['MM'] * 50 + ['PP'] * 50

        # Generate read counts based on true states
        ref_counts = []
        alt_counts = []
        other_counts = []

        for state in true_states:
            if state == 'PP':
                ref_counts.append(np.random.randint(15, 25))
                alt_counts.append(np.random.randint(0, 3))
            elif state == 'PM':
                ref_counts.append(np.random.randint(8, 15))
                alt_counts.append(np.random.randint(8, 15))
            elif state == 'MM':
                ref_counts.append(np.random.randint(0, 3))
                alt_counts.append(np.random.randint(15, 25))
            other_counts.append(np.random.randint(0, 2))

        read_counts = pd.DataFrame({
            'ref': ref_counts,
            'alt': alt_counts,
            'other': other_counts,
        })

        # Run HMM
        hmm = SisterChromatidHMM()
        results = hmm.decode_with_blocks(
            read_counts=read_counts,
            chromosome='chr1',
            positions=positions,
            sample_id='synthetic'
        )

        # Check results
        assert len(results.states) == n_snps
        assert len(results.blocks) >= 1

        # Calculate accuracy (roughly, ignoring phase)
        inferred_states = results.states
        correct = sum(1 for i, s in enumerate(inferred_states) if s == true_states[i])
        accuracy = correct / n_snps

        # Should have reasonable accuracy (>50%)
        assert accuracy > 0.5


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
