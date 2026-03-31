"""
Tests for the Track architecture.

This module tests the Track base class and its implementations:
- EigenTrack
- ISTrack
- HeatmapTrack
- TrackComposer
"""

import os
import unittest
from pathlib import Path

import pandas as pd
import numpy as np
from cooler import Cooler
from plotly.subplots import make_subplots

from hic_basic.plot.track import (
    Track,
    EigenTrack,
    ISTrack,
    HeatmapTrack,
    TrackConfig,
)
from hic_basic.plot.track_composer import (
    TrackComposer,
    LayoutConfig,
)


class TestTrackConfig(unittest.TestCase):
    """Tests for TrackConfig dataclass."""

    def test_default_values(self):
        """Test default values for TrackConfig."""
        config = TrackConfig(row=1, col=1)
        self.assertEqual(config.height_ratio, 1.0)
        self.assertIsNone(config.title)
        self.assertIsNone(config.y_range)
        self.assertEqual(config.style_kwargs, {})

    def test_custom_values(self):
        """Test custom values for TrackConfig."""
        config = TrackConfig(
            row=2,
            col=3,
            height_ratio=2.0,
            title="Test Track",
            y_range=(-1, 1)
        )
        self.assertEqual(config.row, 2)
        self.assertEqual(config.col, 3)
        self.assertEqual(config.height_ratio, 2.0)
        self.assertEqual(config.title, "Test Track")
        self.assertEqual(config.y_range, (-1, 1))


class TestEigenTrack(unittest.TestCase):
    """Tests for EigenTrack class."""

    def setUp(self):
        """Set up test fixtures."""
        self.outdir = Path(os.path.dirname(__file__)) / "output" / "track"
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.synthetic_cool = Path(os.path.dirname(__file__)) / "data" / "synthetic_test.cool"
        self.synthetic_eigs1 = Path(os.path.dirname(__file__)) / "data" / "synthetic_eigs1.bedgraph"

    def test_eigen_track_init(self):
        """Test EigenTrack initialization."""
        track = EigenTrack(str(self.synthetic_eigs1))
        self.assertEqual(track.name, "eigen")
        self.assertEqual(track.eigen_col, "E1")
        self.assertEqual(track.orientation, "v")

    def test_eigen_track_custom_params(self):
        """Test EigenTrack with custom parameters."""
        track = EigenTrack(
            str(self.synthetic_eigs1),
            eigen_col="E1",
            cmap={"A": "red", "B": "blue"},
            orientation="h"
        )
        self.assertEqual(track.eigen_col, "E1")
        self.assertEqual(track.orientation, "h")

    def test_eigen_track_load_data(self):
        """Test EigenTrack data loading."""
        track = EigenTrack(str(self.synthetic_eigs1))
        clr = Cooler(str(self.synthetic_cool))
        track.load_data(clr, "chr1")

        self.assertIsNotNone(track.data)
        self.assertIn("AB", track.data.columns)
        self.assertGreater(len(track.data), 0)

    def test_eigen_track_render(self):
        """Test EigenTrack rendering."""
        track = EigenTrack(str(self.synthetic_eigs1))
        clr = Cooler(str(self.synthetic_cool))
        track.load_data(clr, "chr1")

        fig = make_subplots(rows=1, cols=1)
        track.render(fig, row=1, col=1)

        self.assertGreaterEqual(len(fig.data), 1)

    def test_eigen_track_visual(self):
        """Visual test for EigenTrack."""
        out_png = self.outdir / "eigen_track.png"

        track = EigenTrack(str(self.synthetic_eigs1), orientation="v")
        clr = Cooler(str(self.synthetic_cool))
        track.load_data(clr, "chr1")

        fig = make_subplots(rows=1, cols=1)
        track.render(fig, row=1, col=1)
        fig.write_image(str(out_png), width=800, height=400)

        self.assertTrue(out_png.exists())


class TestISTrack(unittest.TestCase):
    """Tests for ISTrack class."""

    def setUp(self):
        """Set up test fixtures."""
        self.outdir = Path(os.path.dirname(__file__)) / "output" / "track"
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.synthetic_cool = Path(os.path.dirname(__file__)) / "data" / "synthetic_test.cool"
        self.synthetic_eigs1 = Path(os.path.dirname(__file__)) / "data" / "synthetic_eigs1.bedgraph"

    def test_is_track_init(self):
        """Test ISTrack initialization."""
        track = ISTrack(str(self.synthetic_eigs1), is_col="E1")
        self.assertEqual(track.name, "IS")
        self.assertEqual(track.is_col, "E1")
        self.assertEqual(track.line_color, "black")

    def test_is_track_load_data(self):
        """Test ISTrack data loading."""
        track = ISTrack(str(self.synthetic_eigs1), is_col="E1")
        clr = Cooler(str(self.synthetic_cool))
        track.load_data(clr, "chr1")

        self.assertIsNotNone(track.data)
        self.assertGreater(len(track.data), 0)

    def test_is_track_render(self):
        """Test ISTrack rendering."""
        track = ISTrack(str(self.synthetic_eigs1), is_col="E1")
        clr = Cooler(str(self.synthetic_cool))
        track.load_data(clr, "chr1")

        fig = make_subplots(rows=1, cols=1)
        track.render(fig, row=1, col=1)

        self.assertGreaterEqual(len(fig.data), 1)


class TestHeatmapTrack(unittest.TestCase):
    """Tests for HeatmapTrack class."""

    def setUp(self):
        """Set up test fixtures."""
        self.outdir = Path(os.path.dirname(__file__)) / "output" / "track"
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.synthetic_cool = Path(os.path.dirname(__file__)) / "data" / "synthetic_test.cool"

    def test_heatmap_track_init(self):
        """Test HeatmapTrack initialization."""
        track = HeatmapTrack(str(self.synthetic_cool), balance=False)
        self.assertEqual(track.name, "heatmap")
        self.assertFalse(track.balance)

    def test_heatmap_track_load_data(self):
        """Test HeatmapTrack data loading."""
        track = HeatmapTrack(str(self.synthetic_cool), balance=False)
        track.load_data("chr1")

        self.assertIsNotNone(track._matrix)
        self.assertIsNotNone(track._bins)
        self.assertGreater(track._matrix.shape[0], 0)

    def test_heatmap_track_render(self):
        """Test HeatmapTrack rendering."""
        track = HeatmapTrack(str(self.synthetic_cool), balance=False)
        track.load_data("chr1")

        fig = make_subplots(rows=1, cols=1)
        track.render(fig, row=1, col=1)

        self.assertGreaterEqual(len(fig.data), 1)

    def test_heatmap_track_from_matrix(self):
        """Test HeatmapTrack creation from pre-computed matrix."""
        matrix = np.random.rand(10, 10)
        bins = pd.DataFrame({
            "chrom": ["chr1"] * 10,
            "start": list(range(0, 1000000, 100000)),
            "end": list(range(100000, 1100000, 100000))
        })

        track = HeatmapTrack.from_matrix(matrix, bins)
        self.assertEqual(track._matrix.shape, (10, 10))
        self.assertEqual(len(track._bins), 10)


class TestTrackComposer(unittest.TestCase):
    """Tests for TrackComposer class."""

    def setUp(self):
        """Set up test fixtures."""
        self.outdir = Path(os.path.dirname(__file__)) / "output" / "track"
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.synthetic_cool = Path(os.path.dirname(__file__)) / "data" / "synthetic_test.cool"
        self.synthetic_eigs1 = Path(os.path.dirname(__file__)) / "data" / "synthetic_eigs1.bedgraph"
        self.synthetic_eigs2 = Path(os.path.dirname(__file__)) / "data" / "synthetic_eigs2.bedgraph"

    def test_layout_config_default(self):
        """Test LayoutConfig default values."""
        config = LayoutConfig(rows=2, cols=2)
        self.assertEqual(config.rows, 2)
        self.assertEqual(config.cols, 2)
        self.assertEqual(config.row_heights, [1.0, 1.0])
        self.assertEqual(config.col_widths, [1.0, 1.0])

    def test_layout_config_custom(self):
        """Test LayoutConfig with custom values."""
        config = LayoutConfig(
            rows=3,
            cols=1,
            row_heights=[10, 2, 2],
            vertical_spacing=0.05
        )
        self.assertEqual(config.rows, 3)
        self.assertEqual(config.row_heights, [10, 2, 2])
        self.assertEqual(config.vertical_spacing, 0.05)

    def test_composer_add_track(self):
        """Test TrackComposer add_track method."""
        composer = TrackComposer(LayoutConfig(rows=2, cols=2))
        track = EigenTrack(str(self.synthetic_eigs1))

        composer.add_track(track, row=1, col=1)
        self.assertEqual(len(composer.tracks), 1)
        self.assertEqual(composer.get_track(1, 1), track)

    def test_composer_list_tracks(self):
        """Test TrackComposer list_tracks method."""
        composer = TrackComposer(LayoutConfig(rows=2, cols=2))
        track1 = EigenTrack(str(self.synthetic_eigs1))
        track2 = ISTrack(str(self.synthetic_eigs2), is_col="E1")

        composer.add_track(track1, row=1, col=1)
        composer.add_track(track2, row=2, col=1)

        tracks = composer.list_tracks()
        self.assertEqual(len(tracks), 2)

    def test_composer_build(self):
        """Test TrackComposer build method."""
        composer = TrackComposer(LayoutConfig(rows=2, cols=1))
        heatmap = HeatmapTrack(str(self.synthetic_cool), balance=False)
        eigen = EigenTrack(str(self.synthetic_eigs1))

        composer.add_track(heatmap, row=1, col=1)
        composer.add_track(eigen, row=2, col=1)

        clr = Cooler(str(self.synthetic_cool))
        fig = composer.build(region="chr1", title="Test", clr=clr)
        self.assertIsNotNone(fig)
        self.assertGreaterEqual(len(fig.data), 2)

    def test_composer_for_tiling_compartment(self):
        """Test TrackComposer preset for tiling compartment."""
        composer = TrackComposer.for_tiling_compartment()
        self.assertEqual(composer.layout_config.rows, 2)
        self.assertEqual(composer.layout_config.cols, 2)
        self.assertEqual(composer.layout_config.row_heights, [10, 1])
        self.assertEqual(composer.layout_config.col_widths, [1, 10])

    def test_composer_for_trapezoid_2tracks(self):
        """Test TrackComposer preset for trapezoid 2 tracks."""
        composer = TrackComposer.for_trapezoid_2tracks()
        self.assertEqual(composer.layout_config.rows, 3)
        self.assertEqual(composer.layout_config.cols, 1)
        self.assertEqual(composer.layout_config.row_heights, [6, 2, 2])

    def test_composer_for_heatmap_with_tracks(self):
        """Test TrackComposer preset for heatmap with tracks."""
        composer = TrackComposer.for_heatmap_with_tracks(n_tracks=3)
        self.assertEqual(composer.layout_config.rows, 4)
        self.assertEqual(composer.layout_config.cols, 1)
        self.assertEqual(composer.layout_config.row_heights, [10, 2, 2, 2])

    def test_composer_visual(self):
        """Visual test for TrackComposer."""
        out_png = self.outdir / "track_composer.png"

        composer = TrackComposer.for_heatmap_with_tracks(n_tracks=2)
        heatmap = HeatmapTrack(str(self.synthetic_cool), balance=False)
        eigen1 = EigenTrack(str(self.synthetic_eigs1))
        eigen2 = EigenTrack(str(self.synthetic_eigs2))

        composer.add_track(heatmap, row=1, col=1)
        composer.add_track(eigen1, row=2, col=1)
        composer.add_track(eigen2, row=3, col=1)

        clr = Cooler(str(self.synthetic_cool))
        fig = composer.build(region="chr1", title="Track Composer Test", clr=clr)
        fig.write_image(str(out_png), width=800, height=1000)

        self.assertTrue(out_png.exists())


if __name__ == "__main__":
    unittest.main()
