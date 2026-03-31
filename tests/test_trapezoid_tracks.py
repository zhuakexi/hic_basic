"""
Tests for trapezoid Hi-C with multi-track visualization.

This module tests the new Track classes and plot_trapezoid_with_tracks function.
"""
import os
import unittest
from pathlib import Path
import sys

sys.path.insert(0, "/share/home/ychi/dev/hic_basic")

from hic_basic.plot.track import (
    BigWigTrack,
    GeneTrack,
    TrapezoidHeatmapTrack,
    TrackConfig
)
from hic_basic.plot.track_composer import TrackComposer, LayoutConfig
from hic_basic.plot.hic import plot_trapezoid_with_tracks


class TestBigWigTrack(unittest.TestCase):
    """Test BigWigTrack functionality."""

    def setUp(self):
        """Set up test fixtures."""
        # Use synthetic test data if available, otherwise skip
        self.test_bw = Path(os.path.dirname(__file__)) / "data" / "test_signal.bw"
        self.outdir = Path(os.path.dirname(__file__)) / "output" / "trapezoid_tracks"
        self.outdir.mkdir(parents=True, exist_ok=True)

    def test_bigwig_track_load_data(self):
        """Test BigWigTrack data loading."""
        if not self.test_bw.exists():
            self.skipTest("Test bigWig file not found")
        
        track = BigWigTrack(str(self.test_bw), bin_size=1000)
        track.load_data("chr1:1000000-2000000")
        
        self.assertIsNotNone(track._data)
        self.assertIn("start", track._data.columns)
        self.assertIn("signal", track._data.columns)
        self.assertGreater(len(track._data), 0)

    def test_bigwig_track_render(self):
        """Test BigWigTrack rendering."""
        if not self.test_bw.exists():
            self.skipTest("Test bigWig file not found")
        
        from plotly.subplots import make_subplots
        
        track = BigWigTrack(str(self.test_bw), bin_size=1000, color="#1f77b4")
        track.load_data("chr1:1000000-2000000")
        
        fig = make_subplots(rows=1, cols=1)
        track.render(fig, row=1, col=1)
        
        # Should have added a bar trace
        self.assertGreaterEqual(len(fig.data), 1)


class TestGeneTrackInTrackFramework(unittest.TestCase):
    """Test GeneTrack class in Track framework."""

    def setUp(self):
        """Set up test fixtures."""
        self.test_gtf = Path(os.path.dirname(__file__)) / "data" / "synthetic_test.gtf"
        self.outdir = Path(os.path.dirname(__file__)) / "output" / "trapezoid_tracks"
        self.outdir.mkdir(parents=True, exist_ok=True)

    def test_gene_track_load_data(self):
        """Test GeneTrack data loading."""
        track = GeneTrack(str(self.test_gtf), mode="all_transcripts")
        track.load_data("chr1:1000-9000")
        
        self.assertIsNotNone(track._genes)
        self.assertGreater(len(track._shapes), 0)

    def test_gene_track_render(self):
        """Test GeneTrack rendering."""
        from plotly.subplots import make_subplots
        
        track = GeneTrack(str(self.test_gtf), mode="all_transcripts")
        track.load_data("chr1:1000-9000")
        
        fig = make_subplots(rows=1, cols=1)
        track.render(fig, row=1, col=1)
        
        # Should have added shapes
        self.assertGreater(len(fig.layout.shapes), 0)

    def test_gene_track_visual(self):
        """Visual test for GeneTrack."""
        from plotly.subplots import make_subplots
        
        track = GeneTrack(
            str(self.test_gtf),
            mode="all_transcripts",
            show_labels=True
        )
        track.load_data("chr1:1000-9000")
        
        fig = make_subplots(rows=1, cols=1)
        track.render(fig, row=1, col=1)
        fig.update_layout(
            title="GeneTrack Visual Test",
            width=1000,
            height=400
        )
        
        out_png = self.outdir / "gene_track_track_framework.png"
        fig.write_image(str(out_png))
        self.assertTrue(out_png.exists())


class TestTrapezoidWithTracks(unittest.TestCase):
    """Test plot_trapezoid_with_tracks function."""

    def setUp(self):
        """Set up test fixtures."""
        self.synthetic_cool = Path(os.path.dirname(__file__)) / "data" / "synthetic_test.cool"
        self.synthetic_eigs = Path(os.path.dirname(__file__)) / "data" / "synthetic_eigs1.bedgraph"
        self.test_gtf = Path(os.path.dirname(__file__)) / "data" / "synthetic_test.gtf"
        self.outdir = Path(os.path.dirname(__file__)) / "output" / "trapezoid_tracks"
        self.outdir.mkdir(parents=True, exist_ok=True)

    def test_plot_trapezoid_with_tracks_basic(self):
        """Test basic functionality of plot_trapezoid_with_tracks."""
        if not self.synthetic_cool.exists():
            self.skipTest("Test data not found")
        
        # Create simple tracks
        tracks = [
            GeneTrack(str(self.test_gtf), mode="gene_only")
        ]

        fig = plot_trapezoid_with_tracks(
            coolp=str(self.synthetic_cool),
            vec_file=str(self.synthetic_eigs),
            region="chr1",
            tracks=tracks,
            title="Test: Trapezoid + Gene Track",
            strata=20,
            balance=False  # synthetic data has no weights
        )

        self.assertIsNotNone(fig)
        # Should have multiple traces (heatmap + gene shapes)
        self.assertGreaterEqual(len(fig.data), 1)

    def test_plot_trapezoid_with_tracks_visual(self):
        """Visual test for plot_trapezoid_with_tracks."""
        if not self.synthetic_cool.exists():
            self.skipTest("Test data not found")
        
        tracks = [
            GeneTrack(str(self.test_gtf), mode="all_transcripts", show_labels=True)
        ]

        fig = plot_trapezoid_with_tracks(
            coolp=str(self.synthetic_cool),
            vec_file=str(self.synthetic_eigs),
            region="chr1",
            tracks=tracks,
            title="Visual Test: Trapezoid + Multi-Track",
            strata=30,
            track_heights=[3],
            total_height=600,
            balance=False  # synthetic data has no weights
        )

        out_png = self.outdir / "trapezoid_with_tracks.png"
        fig.write_image(str(out_png))
        self.assertTrue(out_png.exists())
        print(f"Visual test output: {out_png}")

    def test_plot_trapezoid_with_tracks_custom_heights(self):
        """Test custom track heights."""
        if not self.synthetic_cool.exists():
            self.skipTest("Test data not found")
        
        tracks = [
            GeneTrack(str(self.test_gtf), mode="gene_only"),
            GeneTrack(str(self.test_gtf), mode="gene_only")
        ]

        fig = plot_trapezoid_with_tracks(
            coolp=str(self.synthetic_cool),
            vec_file=str(self.synthetic_eigs),
            region="chr1",
            tracks=tracks,
            track_heights=[2, 3],  # Second track taller
            trapezoid_height=8,
            title="Test: Custom Heights",
            balance=False  # synthetic data has no weights
        )

        self.assertIsNotNone(fig)


class TestTrackComposerForTrapezoid(unittest.TestCase):
    """Test TrackComposer.for_trapezoid_with_tracks method."""

    def test_for_trapezoid_with_tracks_layout(self):
        """Test layout creation."""
        composer = TrackComposer.for_trapezoid_with_tracks(
            n_tracks=3,
            trapezoid_height=6,
            track_heights=[2, 2, 3]
        )
        
        # Should have 4 rows (1 trapezoid + 3 tracks)
        self.assertEqual(composer.layout_config.rows, 4)
        self.assertEqual(composer.layout_config.cols, 1)
        
        # Check row heights
        expected_heights = [6, 2, 2, 3]
        self.assertEqual(composer.layout_config.row_heights, expected_heights)

    def test_for_trapezoid_with_tracks_default_heights(self):
        """Test default track heights."""
        composer = TrackComposer.for_trapezoid_with_tracks(
            n_tracks=2
        )
        
        # Default track height is 2
        expected_heights = [6, 2, 2]
        self.assertEqual(composer.layout_config.row_heights, expected_heights)

    def test_for_trapezoid_with_tracks_mismatched_heights(self):
        """Test error on mismatched track heights."""
        with self.assertRaises(ValueError):
            TrackComposer.for_trapezoid_with_tracks(
                n_tracks=3,
                track_heights=[2, 2]  # Only 2 heights for 3 tracks
            )


if __name__ == "__main__":
    unittest.main()
