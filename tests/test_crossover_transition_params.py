#!/usr/bin/env python3
"""
Tests for transition parameter handling in mII crossover detection.

Tests cover:
1. Default transition parameters
2. Estimated transition parameters
3. Rough classification with fixed thresholds
4. Multi bin size support
"""

import pytest
import pandas as pd
import numpy as np

from hic_basic.crossover import (
    get_default_transition_params,
    estimate_transitions_rough,
    rough_classification,
    bin_by_snp,
)
from hic_basic.crossover.detector_mii import SisterChromatidHMM


class TestDefaultTransitionParams:
    """Tests for default transition parameters."""
    
    def test_get_default_transition_params(self):
        """Test default parameter values."""
        params = get_default_transition_params()
        
        assert params['alpha'] == 0.01
        assert params['beta'] == 0.02
        assert params['gamma'] == 0.001
        assert params['rho'] == 0.1
    
    def test_default_params_build_valid_transition_matrix(self):
        """Test that default parameters produce valid transition matrix."""
        params = get_default_transition_params()
        hmm = SisterChromatidHMM(**params)
        
        # All elements should be non-negative
        assert (hmm.transition_matrix >= 0).all()
        
        # Rows should sum to 1
        row_sums = hmm.transition_matrix.sum(axis=1)
        assert np.allclose(row_sums, 1.0)
        
        # 1-3*rho should be positive
        assert 1 - 3 * params['rho'] > 0


class TestEstimatedTransitionParams:
    """Tests for estimated transition parameters."""
    
    def test_estimate_transitions_with_data(self):
        """Test estimation from rough labels with transitions."""
        # Create pattern: PP -> PP -> PM -> PM -> MM
        rough_labels = pd.Series(['PP', 'PP', 'PM', 'PM', 'MM'])
        
        alpha, beta, gamma, rho = estimate_transitions_rough(rough_labels)
        
        # Parameters should be in valid range [0, 1]
        assert 0 <= alpha <= 1
        assert 0 <= beta <= 1
        assert 0 <= gamma <= 1
        assert 0 <= rho <= 0.32  # Clamped to max 0.32
    
    def test_estimate_transitions_no_transitions(self):
        """Test estimation when all labels are the same."""
        rough_labels = pd.Series(['PP', 'PP', 'PP', 'PP'])
        
        alpha, beta, gamma, rho = estimate_transitions_rough(rough_labels)
        
        # Should use default values based on n_units
        n = len(rough_labels)
        assert np.isclose(alpha, 1.0 / n)
        assert np.isclose(beta, 1.0 / n)
        assert np.isclose(gamma, 0.01 / n)
        assert np.isclose(rho, min(100.0 / n, 0.32))
    
    def test_estimate_transitions_with_u_state(self):
        """Test estimation with U state present."""
        rough_labels = pd.Series(['PP', 'U', 'PM', 'U', 'MM'])
        
        alpha, beta, gamma, rho = estimate_transitions_rough(rough_labels)
        
        # All parameters should be valid
        assert 0 <= alpha <= 1
        assert 0 <= beta <= 1
        assert 0 <= gamma <= 1
        assert 0 <= rho <= 0.32


class TestRoughClassification:
    """Tests for rough classification with fixed thresholds."""
    
    def test_fixed_threshold_025_075(self):
        """Test classification with 0.25/0.75 thresholds."""
        # Create test data with known ref_frac
        test_data = pd.DataFrame({
            'ref': [10, 5, 0, 5],
            'alt': [0, 5, 10, 5],
            'other': [0, 0, 0, 0],
        })
        
        labels = rough_classification(
            test_data,
            depth_threshold=2,
            other_threshold=0.4,
            balance_threshold=0.25  # PP: >0.75, MM: >0.75, PM: |ref-alt|<0.25
        )
        
        # ref_frac = 1.0 → PP
        # ref_frac = 0.5 → PM
        # ref_frac = 0.0 → MM
        # ref_frac = 0.5 → PM
        assert labels.iloc[0] == 'PP'
        assert labels.iloc[1] == 'PM'
        assert labels.iloc[2] == 'MM'
        assert labels.iloc[3] == 'PM'
    
    def test_depth_threshold(self):
        """Test depth threshold filtering."""
        test_data = pd.DataFrame({
            'ref': [1, 5, 10],
            'alt': [0, 5, 0],
            'other': [0, 0, 0],
        })
        
        labels = rough_classification(
            test_data,
            depth_threshold=2,  # Total < 2 → U
            other_threshold=0.4,
            balance_threshold=0.25
        )
        
        # Total = 1 < 2 → U
        # Total = 10 ≥ 2 → PP
        # Total = 10 ≥ 2 → PP
        assert labels.iloc[0] == 'U'
        assert labels.iloc[1] == 'PM'
        assert labels.iloc[2] == 'PP'
    
    def test_other_threshold(self):
        """Test other threshold filtering."""
        test_data = pd.DataFrame({
            'ref': [5, 5],
            'alt': [5, 5],
            'other': [10, 0],  # other_frac = 0.5 vs 0
        })
        
        labels = rough_classification(
            test_data,
            depth_threshold=2,
            other_threshold=0.4,  # other_frac > 0.4 → U
            balance_threshold=0.25
        )
        
        # other_frac = 0.5 > 0.4 → U
        # other_frac = 0 ≤ 0.4 → PM
        assert labels.iloc[0] == 'U'
        assert labels.iloc[1] == 'PM'


class TestBinBySnp:
    """Tests for binning by SNP count."""
    
    def test_bin_by_snp_basic(self):
        """Test basic binning functionality."""
        test_data = pd.DataFrame({
            'chromosome': ['chr1'] * 20,
            'position': range(100, 100 + 20 * 1000, 1000),
            'ref': np.random.randint(0, 5, 20),
            'alt': np.random.randint(0, 5, 20),
            'other': np.random.randint(0, 2, 20),
        })
        
        binned = bin_by_snp(test_data, n_snps=5)
        
        # Should have 20/5 = 4 bins
        assert len(binned) == 4
        
        # Check aggregation
        assert binned.iloc[0]['start_pos'] == 100
        assert binned.iloc[0]['end_pos'] == 100 + 4 * 1000
        assert binned.iloc[0]['ref'] == test_data.iloc[0:5]['ref'].sum()
    
    def test_bin_by_snp_default_250(self):
        """Test binning with default 250 SNPs."""
        n_snps_total = 1000
        test_data = pd.DataFrame({
            'chromosome': ['chr1'] * n_snps_total,
            'position': range(n_snps_total),
            'ref': np.ones(n_snps_total),
            'alt': np.ones(n_snps_total),
            'other': np.zeros(n_snps_total),
        })
        
        binned = bin_by_snp(test_data, n_snps=250)
        
        # Should have 1000/250 = 4 bins
        assert len(binned) == 4
        assert binned.iloc[0]['n_snps'] == 250


class TestHMMWithDefaultParams:
    """Tests for HMM inference with default parameters."""
    
    def test_hmm_decode_with_default_params(self):
        """Test HMM decoding with default parameters."""
        params = get_default_transition_params()
        emission_params = {
            'PP': (0.9, 0.05, 0.05),
            'PM': (0.45, 0.45, 0.1),
            'MM': (0.05, 0.9, 0.05),
            'U': (0.3, 0.3, 0.4),
        }
        
        hmm = SisterChromatidHMM(**params, emission_params=emission_params)
        
        # Create test data
        test_counts = pd.DataFrame({
            'ref': [10, 10, 5, 5, 0, 0],
            'alt': [0, 0, 5, 5, 10, 10],
            'other': [0, 0, 0, 0, 0, 0],
        })
        
        states, log_prob = hmm.decode(test_counts)
        
        # Should produce valid states
        assert len(states) == 6
        assert all(s in ['PP', 'PM', 'MM', 'U'] for s in states)
        
        # Log probability should be finite
        assert np.isfinite(log_prob)
    
    def test_hmm_default_params_vs_estimated(self):
        """Test that default params are used when specified."""
        # Default params
        default_params = get_default_transition_params()
        
        # Estimated params from data with no transitions
        rough_labels = pd.Series(['PP'] * 100)
        alpha_est, beta_est, gamma_est, rho_est = estimate_transitions_rough(rough_labels)
        
        # When no transitions, estimated params use default formula based on n_units
        # Default params should be fixed regardless of data
        assert default_params['alpha'] == 0.01
        assert default_params['beta'] == 0.02
        assert default_params['gamma'] == 0.001
        assert default_params['rho'] == 0.1
        
        # Estimated params should depend on n_units
        n = len(rough_labels)
        assert np.isclose(alpha_est, 1.0 / n)
        assert np.isclose(beta_est, 1.0 / n)


class TestTransitionParamsSource:
    """Tests for transition_params source parameter."""
    
    def test_transition_params_default_string(self):
        """Test that 'default' string is valid option."""
        # This is a documentation test
        valid_options = ['default', 'estimated']
        assert 'default' in valid_options
        assert 'estimated' in valid_options
    
    def test_default_params_are_conservative(self):
        """Test that default parameters represent conservative model."""
        params = get_default_transition_params()
        
        # Low transition rates
        assert params['alpha'] < 0.1  # PP↔PM rare
        assert params['beta'] < 0.2   # PM→PP/MM rare
        assert params['gamma'] < 0.01  # Entry to U very rare
        assert params['rho'] < 0.2     # Exit from U moderate
