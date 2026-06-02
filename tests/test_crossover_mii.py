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
    simulate_mii_crossover_counts,
)
from hic_basic.crossover.pipeline_mii import (
    detect_sister_chromatid_states_end2end,
    estimate_mii_hmm_parameters,
    plot_mii_assignment_track,
    plot_mii_final_assignment,
)
from hic_basic.plot.track import plot_genome_phased_read_count_track


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

    def test_ambiguous_ratio_is_not_u(self):
        """Intermediate ref/alt ratios should be AMBIG, not U."""
        read_counts = pd.DataFrame({
            'ref': [70],
            'alt': [30],
            'other': [0],
        })

        labels = rough_classification(read_counts)

        assert labels.iloc[0] == 'AMBIG'

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

        assert abs(alpha - 0.0) < 1e-10
        assert abs(beta - 0.02) < 1e-10
        assert abs(gamma - 0.0) < 1e-10
        assert abs(rho - 0.1) < 1e-10


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

    def test_mii_crossover_simulation_requested_example(self):
        """Test the requested 100 Mb, 1 kb SNP, 5x mII simulation setup."""
        sim = simulate_mii_crossover_counts(random_state=7)
        df = sim.read_counts

        assert len(df) == 100_000
        assert set(["chromosome", "position", "ref", "alt", "other", "true_state"]).issubset(df.columns)
        assert df["position"].iloc[0] == 1
        assert df["position"].diff().dropna().eq(1_000).all()
        assert set(df["true_state"]).issubset({"PP", "PM", "MM"})

        observed_depth = df[["ref", "alt", "other"]].sum(axis=1)
        assert 3.5 < observed_depth.mean() < 6.5
        assert "chromatid_1" in sim.crossover_positions
        assert "chromatid_2" in sim.crossover_positions

    def test_mii_estimate_only_simulated_binned250_with_diagnostics(self):
        """Estimate parameters from simulated mII data and preserve diagnostic PNGs."""
        sim = simulate_mii_crossover_counts(random_state=11)
        out_dir = Path("issues/mII_crossover_sim_untracked")
        out_dir.mkdir(parents=True, exist_ok=True)
        prefix = out_dir / "simulated_100Mb_1kb_depth5_bin250"

        params = estimate_mii_hmm_parameters(
            sim.read_counts,
            transition_params="estimated",
            binning_method="by_snp",
            binning_snp_number=250,
            params_vis_prefix=prefix,
            verbose=1,
        )

        assert len(params["counts_df"]) == 400
        assert params["rough_labels"] is not None
        assert params["emission_source"] == "estimated_with_by_snp_u_prior"
        assert params["transition_source"] == "estimated"
        assert params["transition_matrix"].shape == (4, 4)
        assert np.allclose(params["transition_matrix"].sum(axis=1), 1.0)
        assert set(params["emission_params"]) == {"PP", "PM", "MM", "U"}
        assert "effective_depth_threshold" in params
        assert params["effective_depth_threshold"].iloc[0] == 500

        for suffix in [
            "rough_label_classification",
            "rough_decision_distributions",
            "rough_segment_lengths",
            "rough_label_filtering",
            "emission_matrix",
            "transition_matrix",
        ]:
            png = Path(f"{prefix}.{suffix}.png")
            assert png.exists()
            assert png.stat().st_size > 0

    def test_transition_estimation_skips_short_rough_segments(self):
        """Short raw rough-label segments should become gaps, not joined transitions."""
        import hic_basic.crossover.pipeline_mii as pipeline_mii

        counts = pd.DataFrame({
            "chromosome": ["chr1"] * 5,
            "position": [1, 2, 3, 4, 5],
            "ref": [10, 10, 0, 5, 5],
            "alt": [0, 0, 0, 5, 5],
            "other": [0, 0, 10, 0, 0],
        })
        labels = pd.Series(["PP", "PP", "U", "PM", "PM"])
        _, _, raw_gamma, _ = estimate_transitions_rough(labels)

        segments = pipeline_mii._rough_label_segments(counts, labels)
        filtered = pipeline_mii._filter_short_rough_segments(
            labels,
            segments,
            min_segment_bins=2,
        )
        transition_labels = pipeline_mii._labels_with_chromosome_gap_separators(
            counts,
            filtered,
        )
        filtered_alpha, _, filtered_gamma, _ = estimate_transitions_rough(transition_labels)

        assert filtered.tolist() == ["PP", "PP", pipeline_mii.GAP_LABEL, "PM", "PM"]
        assert raw_gamma > 0
        assert filtered_gamma == 0
        assert filtered_alpha == 0

    def test_by_snp_estimated_transition_uses_u_prior(self):
        """SNP-binned transition estimation should not estimate U transitions from rough U."""
        counts = pd.DataFrame({
            "chromosome": ["chr1"] * 7,
            "position": np.arange(7),
            "ref": [20, 20, 0, 0, 0, 10, 10],
            "alt": [0, 0, 0, 0, 0, 10, 10],
            "other": [0, 0, 20, 20, 20, 0, 0],
        })

        params = estimate_mii_hmm_parameters(
            counts,
            transition_params="estimated",
            binning_method="by_snp",
            binning_snp_number=1,
            transition_min_segment_bins=1,
            transition_u_gamma_factor=0.1,
            transition_u_rho_prior=0.32,
        )

        assert params["transition_scale_info"]["by_snp_u_prior_applied"] is True
        assert params["transition_scale_info"]["u_labels_treated_as_gaps_for_transition"] is True
        assert "U" not in set(params["transition_rough_labels"])
        assert "U" not in set(params["filtered_rough_labels"])
        assert params["transition_params"]["gamma"] == 0.001
        assert params["transition_params"]["rho"] == 0.32
        assert params["emission_params"]["U"] == (0.475, 0.475, 0.05)
        assert params["u_emission_prior_applied"] is True

    def test_ambiguous_labels_are_gaps_not_u_emission(self):
        """By-SNP AMBIG bins should not train U emission or transition counts."""
        import hic_basic.crossover.pipeline_mii as pipeline_mii

        counts = pd.DataFrame({
            "chromosome": ["chr1"] * 5,
            "position": np.arange(5),
            "ref": [20, 70, 10, 0, 10],
            "alt": [0, 30, 10, 20, 10],
            "other": [0, 0, 0, 0, 0],
        })

        params = estimate_mii_hmm_parameters(
            counts,
            transition_params="estimated",
            binning_method="by_snp",
            binning_snp_number=1,
            transition_min_segment_bins=1,
        )

        assert params["rough_labels"].iloc[1] == "AMBIG"
        assert params["filtered_rough_labels"].iloc[1] == pipeline_mii.GAP_LABEL
        assert params["emission_params"]["U"] == (0.475, 0.475, 0.05)

    def test_estimate_diagnostic_plot_scaling(self, monkeypatch, tmp_path):
        """Diagnostic plots should use log/pseudo-log x and robust count y limits."""
        import hic_basic.crossover.pipeline_mii as pipeline_mii

        captured = []

        def capture_fig(fig, output_path, width=1200, height=700):
            captured.append(fig)

        monkeypatch.setattr(pipeline_mii, "_write_plotly_png", capture_fig)

        normal_rows = pd.DataFrame({
            "chromosome": ["chr1"] * 200,
            "position": np.arange(200),
            "ref": [10] * 200,
            "alt": [12] * 200,
            "other": [0] * 200,
        })
        extreme_row = pd.DataFrame({
            "chromosome": ["chr1"],
            "position": [201],
            "ref": [10000],
            "alt": [0],
            "other": [0],
        })
        counts = pd.concat([normal_rows, extreme_row], ignore_index=True)
        labels = pd.Series(["PM"] * len(counts))

        pipeline_mii._plot_rough_decision_distributions_png(
            counts_df=counts,
            output_path=tmp_path / "decision.png",
            effective_depth_threshold=2,
            other_threshold=0.4,
            balance_threshold=0.25,
        )
        decision_fig = captured[-1]
        assert "log10(total depth + 1)" in decision_fig.layout.xaxis.title.text
        assert "log10(other fraction + eps)" in decision_fig.layout.xaxis2.title.text

        pipeline_mii._plot_rough_label_count_png(
            counts_df=counts,
            rough_labels=labels,
            output_path=tmp_path / "label_count.png",
        )
        label_count_fig = captured[-1]
        assert label_count_fig.layout.yaxis2.range[1] < 10000
        assert any(
            "y-axis capped at p99.5" in annotation.text
            for annotation in label_count_fig.layout.annotations
        )

        pipeline_mii._plot_rough_segment_lengths_png(
            segment_table=pipeline_mii._rough_label_segments(counts, labels),
            min_segment_bins=3,
            output_path=tmp_path / "segments.png",
        )
        segment_fig = captured[-1]
        assert "Raw rough label segment length distribution" in segment_fig.layout.title.text
        assert any(
            "transition_min_segment_bins=3" in annotation.text
            for annotation in segment_fig.layout.annotations
        )

        filtered_labels = labels.copy()
        filtered_labels.iloc[0] = pipeline_mii.GAP_LABEL
        pipeline_mii._plot_raw_filtered_label_png(
            counts_df=counts,
            raw_labels=labels,
            filtered_labels=filtered_labels,
            output_path=tmp_path / "filtering.png",
        )
        filtering_fig = captured[-1]
        assert "Raw rough labels vs transition-filtered labels" in filtering_fig.layout.title.text

        pipeline_mii._plot_assigned_label_count_png(
            counts_df=counts,
            assigned_labels=labels,
            output_path=tmp_path / "assigned.png",
        )
        assigned_fig = captured[-1]
        assert "Final assigned labels" in assigned_fig.layout.title.text
        assert assigned_fig.layout.plot_bgcolor == "white"
        assert assigned_fig.layout.xaxis.showgrid is False
        assert assigned_fig.layout.xaxis.showline is True
        assert assigned_fig.layout.yaxis.showgrid is False
        assert assigned_fig.layout.yaxis.showline is True

        public_fig = pipeline_mii.plot_mii_final_assignment(
            counts.assign(state=labels),
            region="chr1:50-120",
            output_path=tmp_path / "public_assigned.png",
        )
        ref_trace = next(trace for trace in public_fig.data if trace.name == "ref")
        assert min(ref_trace.x) >= 50
        assert max(ref_trace.x) <= 120
        assert public_fig.layout.paper_bgcolor == "white"

    def test_independent_mii_tracks_write_review_pngs(self):
        """Independent count, label, and final-assignment plots should write review PNGs."""
        output_dir = Path("tests/output/crossover_mii_tracks")
        output_dir.mkdir(parents=True, exist_ok=True)

        positions = np.arange(0, 1_200_000, 100_000)
        counts = pd.DataFrame({
            "chromosome": ["chr1"] * 6 + ["chr2"] * 6,
            "start_pos": np.tile(positions[:6], 2),
            "end_pos": np.tile(positions[:6] + 99_999, 2),
            "position": np.tile(positions[:6] + 50_000, 2),
            "ref": [22, 20, 12, 6, 3, 2, 4, 5, 11, 18, 24, 28],
            "alt": [2, 4, 11, 17, 25, 27, 24, 22, 12, 7, 4, 2],
            "other": [0, 1, 1, 2, 1, 0, 0, 2, 1, 1, 0, 0],
            "state": ["PP", "PP", "PM", "PM", "MM", "MM", "MM", "MM", "PM", "PM", "PP", "PP"],
        })

        count_png = output_dir / "plot_genome_phased_read_count_track.png"
        label_png = output_dir / "plot_mii_assignment_track.png"
        final_png = output_dir / "plot_mii_final_assignment.png"

        count_fig = plot_genome_phased_read_count_track(
            counts,
            output_path=count_png,
            title="Synthetic phased read counts",
        )
        label_fig = plot_mii_assignment_track(
            counts[["chromosome", "start_pos", "end_pos", "position", "state"]],
            output_path=label_png,
            title="Synthetic mII labels",
        )
        final_fig = plot_mii_final_assignment(
            counts,
            output_path=final_png,
            title="Synthetic final assignment",
        )

        assert {trace.name for trace in count_fig.data} == {"ref", "alt", "other"}
        assert any(str(trace.name).startswith("PP") for trace in label_fig.data)
        assert len(final_fig.data) == len(count_fig.data) + len(label_fig.data)
        for png in [count_png, label_png, final_png]:
            assert png.exists()
            assert png.stat().st_size > 0

    def test_phased_read_count_track_uses_natural_chromosome_order(self):
        """Chromosome grid should follow chr1, chr2, chr10 and align to trace x positions."""
        counts = pd.DataFrame({
            "chromosome": ["chr10", "chr1", "chr2"],
            "position": [100, 100, 100],
            "ref": [10, 1, 2],
            "alt": [0, 0, 0],
            "other": [0, 0, 0],
            "state": ["MM", "PP", "PM"],
        })

        count_fig = plot_genome_phased_read_count_track(counts)
        ref_trace = next(trace for trace in count_fig.data if trace.name == "ref")
        assert ref_trace.x.tolist() == [0, 1, 2]
        assert ref_trace.customdata[:, 0].tolist() == ["chr1", "chr2", "chr10"]
        assert [shape.x0 for shape in count_fig.layout.shapes] == [1.0, 2.0]
        chrom_annotations = [
            annotation.text
            for annotation in count_fig.layout.annotations
            if annotation.text in {"chr1", "chr2", "chr10"}
        ]
        assert chrom_annotations == ["chr1", "chr2", "chr10"]

        label_fig = plot_mii_assignment_track(
            counts[["chromosome", "position", "state"]]
        )
        label_points = []
        for trace in label_fig.data:
            label_points.extend(zip(trace.x.tolist(), trace.customdata[:, 0].tolist()))
        assert [chrom for _, chrom in sorted(label_points)] == ["chr1", "chr2", "chr10"]

    def test_transition_param_unit_mb_scales_specified_params(self):
        """Specified per-Mb transition rates should scale to bin-step probabilities."""
        sim = simulate_mii_crossover_counts(chromosome_length=1_000_000, random_state=31)
        emission_params = {
            "MM": (0.01, 0.98, 0.01),
            "PM": (0.49, 0.49, 0.02),
            "PP": (0.98, 0.01, 0.01),
            "U": (0.25, 0.25, 0.50),
        }
        transition_params = {
            "alpha": 0.001,
            "beta": 0.002,
            "gamma": 0.003,
            "rho": 0.004,
        }

        params = estimate_mii_hmm_parameters(
            sim.read_counts,
            emission_params=emission_params,
            transition_params=transition_params,
            transition_param_unit="mb",
            binning_method="by_snp",
            binning_snp_number=250,
        )

        distance_mb = params["transition_scale_info"]["median_distance_mb"]
        assert np.isclose(distance_mb, 0.25)
        assert np.isclose(params["transition_params"]["alpha"], 1 - np.exp(-0.001 * 0.25))
        assert params["transition_scale_info"]["scaling_applied"] is True

    def test_end2end_estimate_use_filters_parameter_chromosomes(self, monkeypatch, tmp_path):
        """Multi-chromosome estimate_only should estimate on estimate_use subset."""
        sim1 = simulate_mii_crossover_counts(
            chromosome="chr1",
            chromosome_length=1_000_000,
            random_state=41,
        )
        sim2 = simulate_mii_crossover_counts(
            chromosome="chr2",
            chromosome_length=1_000_000,
            random_state=42,
        )
        read_counts = pd.concat([
            sim1.read_counts[["chromosome", "position", "ref", "alt", "other"]],
            sim2.read_counts[["chromosome", "position", "ref", "alt", "other"]],
        ], ignore_index=True)

        import hic_basic.crossover.pipeline_mii as pipeline_mii

        monkeypatch.setattr(
            pipeline_mii,
            "extract_read_counts",
            lambda *args, **kwargs: read_counts.copy(),
        )

        bam_file = tmp_path / "sample.bam"
        snp_file = tmp_path / "snps.tsv"
        bam_file.touch()
        snp_file.touch()

        result = detect_sister_chromatid_states_end2end(
            bam_file=str(bam_file),
            snp_file=str(snp_file),
            output_dir=str(tmp_path / "out"),
            sample_id="sim",
            chromosomes=["chr1", "chr2"],
            estimate_use=["chr1"],
            transition_params="estimated",
            binning_method="by_snp",
            binning_snp_number=250,
            estimate_only=True,
        )

        assert set(result["counts_df"]["chromosome"]) == {"chr1"}
        assert len(result["counts_df"]) == 4

    def test_end2end_multichrom_decodes_separately_with_shared_params(self, monkeypatch, tmp_path):
        """Multi-chromosome pipeline should share estimated params but decode per chromosome."""
        sim1 = simulate_mii_crossover_counts(
            chromosome="chr1",
            chromosome_length=1_000_000,
            random_state=51,
        )
        sim2 = simulate_mii_crossover_counts(
            chromosome="chr2",
            chromosome_length=1_000_000,
            random_state=52,
        )
        read_counts = pd.concat([
            sim1.read_counts[["chromosome", "position", "ref", "alt", "other"]],
            sim2.read_counts[["chromosome", "position", "ref", "alt", "other"]],
        ], ignore_index=True)

        import hic_basic.crossover.pipeline_mii as pipeline_mii

        monkeypatch.setattr(
            pipeline_mii,
            "extract_read_counts",
            lambda *args, **kwargs: read_counts.copy(),
        )

        bam_file = tmp_path / "sample.bam"
        snp_file = tmp_path / "snps.tsv"
        bam_file.touch()
        snp_file.touch()

        results = detect_sister_chromatid_states_end2end(
            bam_file=str(bam_file),
            snp_file=str(snp_file),
            output_dir=str(tmp_path / "out"),
            sample_id="sim",
            chromosomes=["chr1", "chr2"],
            estimate_use="all",
            transition_params="estimated",
            binning_method="by_snp",
            binning_snp_number=250,
        )

        assert set(results) == {"chr1", "chr2"}
        assert results["chr1"].chromosome == "chr1"
        assert results["chr2"].chromosome == "chr2"
        assert len(results["chr1"].positions) == 1000
        assert len(results["chr2"].positions) == 1000
        assert (tmp_path / "out" / "sim_hmm_params.json").exists()

    def test_end2end_estimate_only_branch_with_simulated_counts(self, monkeypatch, tmp_path):
        """Test the public pipeline estimate_only branch without requiring a real BAM."""
        sim = simulate_mii_crossover_counts(
            chromosome_length=1_000_000,
            random_state=23,
        )
        read_counts = sim.read_counts[["chromosome", "position", "ref", "alt", "other"]]

        import hic_basic.crossover.pipeline_mii as pipeline_mii

        monkeypatch.setattr(
            pipeline_mii,
            "extract_read_counts",
            lambda *args, **kwargs: read_counts.copy(),
        )

        bam_file = tmp_path / "sample.bam"
        snp_file = tmp_path / "snps.tsv"
        bam_file.touch()
        snp_file.touch()

        result = detect_sister_chromatid_states_end2end(
            bam_file=str(bam_file),
            snp_file=str(snp_file),
            output_dir=str(tmp_path / "out"),
            sample_id="sim",
            transition_params="estimated",
            binning_method="by_snp",
            binning_snp_number=250,
            estimate_only=True,
        )

        assert result["transition_source"] == "estimated"
        assert result["emission_source"] == "estimated_with_by_snp_u_prior"
        assert len(result["counts_df"]) == 4
        assert result["rough_labels"] is not None
        assert (tmp_path / "out" / "sim_hmm_params.json").exists()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
