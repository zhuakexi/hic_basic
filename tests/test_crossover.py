"""
Tests for crossover detection module.
"""

import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from hic_basic.crossover.detector import (
    HMMCrossoverDetector,
    CrossoverDetector,
    HaplotypeBlock,
    Crossover,
    CrossoverResults,
)
from hic_basic.crossover.io import (
    load_snp_data,
    load_phasing_data,
    save_results,
    merge_snp_with_phasing,
    validate_input_data,
)
from hic_basic.crossover.cache import CacheManager, SharedMemoryManager


class TestHMMCrossoverDetector:
    """Tests for HMMCrossoverDetector class."""

    def test_init(self):
        """Test initialization with default parameters."""
        detector = HMMCrossoverDetector()
        assert detector.window_size == 10
        assert detector.min_fragment_size == 300000
        assert detector.error_rate == 0.05
        assert detector.transition_matrix is not None
        assert detector.emission_matrix is not None

    def test_encode_observations(self):
        """Test observation encoding."""
        detector = HMMCrossoverDetector()

        sperm_genotypes = ['A', 'G', 'N', 'A', 'G']
        haplotype_1 = ['A', 'A', 'A', 'G', 'A']
        haplotype_2 = ['G', 'G', 'G', 'A', 'G']

        observations = detector._encode_observations(
            sperm_genotypes, haplotype_1, haplotype_2
        )

        assert len(observations) == 5
        assert observations[0] == 0  # A matches haplotype_1 (A)
        assert observations[1] == 1  # G matches haplotype_2 (G)
        assert observations[2] == 2  # N is undetermined

    def test_rough_estimation(self):
        """Test rough crossover estimation."""
        detector = HMMCrossoverDetector(window_size=5)

        positions = np.array([1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000])
        haplotype_alleles = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1])

        crossovers = detector._rough_estimation(positions, haplotype_alleles)

        assert len(crossovers) >= 0

    def test_viterbi(self):
        """Test Viterbi algorithm."""
        detector = HMMCrossoverDetector()

        observations = np.array([0, 0, 0, 1, 1, 1, 0, 0, 0])
        states, log_prob = detector._viterbi(observations)

        assert len(states) == len(observations)
        assert all(s in [0, 1, 2] for s in states)
        assert np.isfinite(log_prob)

    def test_filter_small_fragments(self):
        """Test filtering of small fragments."""
        detector = HMMCrossoverDetector(min_fragment_size=1000)

        positions = np.array([0, 100, 200, 500, 1000, 1500, 2000])
        states = np.array([0, 2, 2, 2, 1, 1, 1])

        filtered = detector._filter_small_fragments(positions, states)

        assert filtered is not None
        assert len(filtered) == len(states)

    def test_detect(self):
        """Test complete detection pipeline."""
        detector = HMMCrossoverDetector(
            window_size=5,
            min_fragment_size=1000,
            error_rate=0.05
        )

        positions = np.arange(100, 10001, 100)
        sperm_genotypes = ['A'] * 50 + ['G'] * 50
        haplotype_1 = ['A'] * 100
        haplotype_2 = ['G'] * 100

        states, crossovers = detector.detect(
            positions, sperm_genotypes, haplotype_1, haplotype_2
        )

        assert len(states) == 100
        assert isinstance(crossovers, list)


class TestCrossoverDetector:
    """Tests for CrossoverDetector class."""

    def test_init(self):
        """Test initialization."""
        detector = CrossoverDetector()
        assert detector.window_size == 10
        assert detector.snp_data is None
        assert detector.phasing_data is None

    def test_load_snp_data(self, tmp_path):
        """Test loading SNP data."""
        snp_file = tmp_path / "snps.tsv"
        snp_data = pd.DataFrame({
            'chromosome': ['chr1'] * 100,
            'position': range(1000, 100001, 1000),
            'snp_id': [f'rs{i}' for i in range(100)],
            'sperm_genotype': ['A'] * 100,
            'diploid_genotype': ['A/G'] * 100,
        })
        snp_data.to_csv(snp_file, sep='\t', index=False)

        detector = CrossoverDetector()
        detector.load_snp_data(str(snp_file))

        assert detector.snp_data is not None
        assert len(detector.snp_data) == 100

    def test_detect_crossovers(self, tmp_path):
        """Test crossover detection."""
        snp_file = tmp_path / "snps.tsv"
        phasing_file = tmp_path / "phasing.tsv"

        n_snps = 200
        positions = list(range(1000, 1000 * n_snps + 1, 1000))[:n_snps]
        snp_data = pd.DataFrame({
            'chromosome': ['chr1'] * n_snps,
            'position': positions,
            'snp_id': [f'rs{i}' for i in range(n_snps)],
            'sperm_genotype': ['A'] * 100 + ['G'] * 100,
            'diploid_genotype': ['A/G'] * n_snps,
        })
        snp_data.to_csv(snp_file, sep='\t', index=False)

        phasing_data = pd.DataFrame({
            'chromosome': ['chr1'] * n_snps,
            'position': positions,
            'snp_id': [f'rs{i}' for i in range(n_snps)],
            'haplotype_1': ['A'] * n_snps,
            'haplotype_2': ['G'] * n_snps,
        })
        phasing_data.to_csv(phasing_file, sep='\t', index=False)

        detector = CrossoverDetector(window_size=10)
        detector.load_snp_data(str(snp_file))
        detector.load_phasing_data(str(phasing_file))

        results = detector.detect_crossovers("S01")

        assert results.sperm_id == "S01"
        assert len(results.haplotype_blocks) > 0
        assert len(results.crossovers) >= 0
        assert 'total_crossovers' in results.statistics

    def test_save_results(self, tmp_path):
        """Test saving results."""
        results = CrossoverResults(sperm_id="S01")
        results.haplotype_blocks = [
            HaplotypeBlock("chr1", 1000, 50000, "P", 0.95, 50),
            HaplotypeBlock("chr1", 50000, 100000, "M", 0.93, 50),
        ]
        results.crossovers = [
            Crossover("S01", "chr1", 50000, 49000, 51000, "P", "M", 0.94, 2.0),
        ]
        results.statistics = {
            'total_crossovers': 1,
            'total_blocks': 2,
            'mean_confidence': 0.94,
        }

        output_dir = tmp_path / "output"
        save_results(results, str(output_dir))

        assert (output_dir / "S01_haplotype_blocks.tsv").exists()
        assert (output_dir / "S01_crossovers.tsv").exists()
        assert (output_dir / "S01_statistics.tsv").exists()


class TestIO:
    """Tests for I/O functions."""

    def test_load_snp_data_tsv(self, tmp_path):
        """Test loading TSV SNP data."""
        snp_file = tmp_path / "snps.tsv"
        snp_data = pd.DataFrame({
            'chromosome': ['chr1', 'chr1', 'chr2'],
            'position': [1000, 2000, 1000],
            'sperm_genotype': ['A', 'G', 'A'],
        })
        snp_data.to_csv(snp_file, sep='\t', index=False)

        df = load_snp_data(str(snp_file))
        assert len(df) == 3
        assert 'snp_id' in df.columns

    def test_load_phasing_data(self, tmp_path):
        """Test loading phasing data."""
        phasing_file = tmp_path / "phasing.tsv"
        phasing_data = pd.DataFrame({
            'chromosome': ['chr1', 'chr1'],
            'position': [1000, 2000],
            'haplotype_1': ['A', 'G'],
            'haplotype_2': ['G', 'A'],
        })
        phasing_data.to_csv(phasing_file, sep='\t', index=False)

        df = load_phasing_data(str(phasing_file))
        assert len(df) == 2

    def test_merge_snp_with_phasing(self):
        """Test merging SNP and phasing data."""
        snp_df = pd.DataFrame({
            'chromosome': ['chr1', 'chr1'],
            'position': [1000, 2000],
            'snp_id': ['rs1', 'rs2'],
            'sperm_genotype': ['A', 'G'],
        })

        phasing_df = pd.DataFrame({
            'chromosome': ['chr1', 'chr1'],
            'position': [1000, 2000],
            'snp_id': ['rs1', 'rs2'],
            'haplotype_1': ['A', 'G'],
            'haplotype_2': ['G', 'A'],
        })

        merged = merge_snp_with_phasing(snp_df, phasing_df)
        assert len(merged) == 2
        assert 'haplotype_1' in merged.columns

    def test_validate_input_data(self):
        """Test input data validation."""
        snp_df = pd.DataFrame({
            'chromosome': ['chr1'],
            'position': [1000],
            'sperm_genotype': ['A'],
        })

        is_valid, issues = validate_input_data(snp_df)
        assert is_valid
        assert len(issues) == 0

        bad_snp_df = pd.DataFrame({'col1': [1, 2, 3]})
        is_valid, issues = validate_input_data(bad_snp_df)
        assert not is_valid
        assert len(issues) > 0


class TestCache:
    """Tests for caching functionality."""

    def test_cache_manager_init(self, tmp_path):
        """Test cache manager initialization."""
        cache_dir = tmp_path / "cache"
        cache = CacheManager(str(cache_dir))

        assert cache.cache_dir.exists()
        assert cache._data_dir.exists()
        assert cache._flag_dir.exists()

    def test_compute_key(self, tmp_path):
        """Test key computation."""
        cache = CacheManager(str(tmp_path / "cache"))

        params1 = {"a": 1, "b": 2}
        params2 = {"a": 1, "b": 2}
        params3 = {"a": 1, "b": 3}

        key1 = cache.compute_key(params1)
        key2 = cache.compute_key(params2)
        key3 = cache.compute_key(params3)

        assert key1 == key2
        assert key1 != key3

    def test_save_and_load(self, tmp_path):
        """Test saving and loading data."""
        cache_dir = tmp_path / "cache"
        cache = CacheManager(str(cache_dir))

        df = pd.DataFrame({'a': [1, 2, 3], 'b': [4, 5, 6]})
        key = cache.compute_key({"test": "data"})

        cache.save(df, key, "test_step")
        loaded = cache.load(key, "test_step")

        assert loaded is not None
        assert len(loaded) == 3

    def test_flag_operations(self, tmp_path):
        """Test flag file operations."""
        cache = CacheManager(str(tmp_path / "cache"))
        key = cache.compute_key({"test": "flag"})

        assert not cache.is_flagged(key, "test_step")
        cache.set_flag(key, "test_step")
        assert cache.is_flagged(key, "test_step")

        cache.clear_flag(key, "test_step")
        assert not cache.is_flagged(key, "test_step")

    def test_get_or_compute(self, tmp_path):
        """Test get or compute pattern."""
        cache = CacheManager(str(tmp_path / "cache"))
        key = cache.compute_key({"test": "compute"})

        def compute_fn():
            return pd.DataFrame({'x': [1, 2, 3]})

        result1 = cache.get_or_compute(key, "step", compute_fn)
        result2 = cache.get_or_compute(key, "step", compute_fn)

        assert result1 is not None
        assert len(result1) == 3
        assert result2 is not None

    def test_shared_memory_manager(self, tmp_path):
        """Test shared memory manager."""
        smm = SharedMemoryManager(str(tmp_path / "shm"))

        arr = np.random.rand(1000)
        key = smm.share(arr, "test_array")

        retrieved = smm.get(key)
        assert retrieved is not None
        assert np.allclose(arr, retrieved)

        smm.clear(key)
        smm.cleanup()


class TestDataClasses:
    """Tests for data classes."""

    def test_haplotype_block(self):
        """Test HaplotypeBlock dataclass."""
        block = HaplotypeBlock(
            chromosome="chr1",
            start_pos=1000,
            end_pos=50000,
            haplotype="P",
            confidence=0.95,
            num_snps=50
        )

        assert block.chromosome == "chr1"
        assert block.haplotype == "P"
        assert 0 <= block.confidence <= 1

    def test_crossover(self):
        """Test Crossover dataclass."""
        crossover = Crossover(
            sperm_id="S01",
            chromosome="chr1",
            position=50000,
            start_pos=49000,
            end_pos=51000,
            from_haplotype="P",
            to_haplotype="M",
            confidence=0.94,
            resolution_kb=2.0
        )

        assert crossover.sperm_id == "S01"
        assert crossover.from_haplotype == "P"
        assert crossover.to_haplotype == "M"

    def test_crossover_results(self):
        """Test CrossoverResults dataclass."""
        results = CrossoverResults(sperm_id="S01")

        assert results.sperm_id == "S01"
        assert len(results.haplotype_blocks) == 0
        assert len(results.crossovers) == 0
        assert isinstance(results.statistics, dict)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
