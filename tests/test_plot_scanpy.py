"""
Unit tests for hic_basic.plot.scanpy – focusing on plot_pca's `background` parameter.
"""
import sys
sys.path.insert(0, "/share/home/ychi/dev/hic_basic")

import types
import unittest

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from hic_basic.plot.scanpy import plot_pca


def _make_adata(n_cells=20, n_pcs=4, seed=0):
    """
    Build a minimal AnnData-like namespace sufficient for plot_pca.

    Attributes
    ----------
    obsm["X_pca"] : np.ndarray  shape (n_cells, n_pcs)
    obs_names     : pd.Index
    obs           : pd.DataFrame  with columns ["cell_type", "batch"]
    """
    rng = np.random.default_rng(seed)
    adata = types.SimpleNamespace()
    adata.obsm = {"X_pca": rng.standard_normal((n_cells, n_pcs))}
    adata.obs_names = pd.Index([f"cell_{i}" for i in range(n_cells)])
    adata.obs = pd.DataFrame(
        {
            "cell_type": ["A"] * (n_cells // 2) + ["B"] * (n_cells // 2),
            "batch": (["X"] * (n_cells // 4) + ["Y"] * (n_cells // 4)) * 2,
        },
        index=adata.obs_names,
    )
    return adata


class TestPlotPcaBackground(unittest.TestCase):
    """Tests for the `background` parameter of plot_pca."""

    def setUp(self):
        self.adata = _make_adata()

    # ------------------------------------------------------------------
    # helpers
    # ------------------------------------------------------------------
    def _bg_traces(self, fig):
        """Return all traces whose legendgroup is 'background'."""
        return [t for t in fig.data if getattr(t, "legendgroup", None) == "background"]

    def _fg_traces(self, fig):
        """Return all traces that are NOT the background trace."""
        return [t for t in fig.data if getattr(t, "legendgroup", None) != "background"]

    # ------------------------------------------------------------------
    # basic smoke tests (no background)
    # ------------------------------------------------------------------
    def test_no_background_returns_figure(self):
        fig = plot_pca(self.adata, color="cell_type")
        self.assertIsInstance(fig, go.Figure)

    def test_no_background_all_samples_in_traces(self):
        fig = plot_pca(self.adata, color="cell_type")
        n_points = sum(len(t.x) for t in fig.data)
        self.assertEqual(n_points, len(self.adata.obs_names))

    # ------------------------------------------------------------------
    # background=None is explicit default
    # ------------------------------------------------------------------
    def test_background_none_equivalent_to_omitted(self):
        fig1 = plot_pca(self.adata, color="cell_type")
        fig2 = plot_pca(self.adata, color="cell_type", background=None)
        self.assertEqual(len(fig1.data), len(fig2.data))

    # ------------------------------------------------------------------
    # single plot with background (default PCs=None)
    # ------------------------------------------------------------------
    def test_background_creates_grey_trace(self):
        fig = plot_pca(self.adata, color="cell_type", background="batch == 'X'")
        bg = self._bg_traces(fig)
        self.assertTrue(len(bg) >= 1, "Expected at least one background trace")
        # marker colour must be lightgrey
        self.assertEqual(bg[0].marker.color, "lightgrey")

    def test_background_total_point_count_unchanged(self):
        """Background + foreground must cover all cells exactly once."""
        bg_query = "batch == 'X'"
        n_bg = len(self.adata.obs.query(bg_query))
        n_fg = len(self.adata.obs) - n_bg
        fig = plot_pca(self.adata, color="cell_type", background=bg_query)
        bg_points = sum(len(t.x) for t in self._bg_traces(fig))
        fg_points = sum(len(t.x) for t in self._fg_traces(fig))
        self.assertEqual(bg_points, n_bg)
        self.assertEqual(fg_points, n_fg)

    def test_background_trace_drawn_before_foreground(self):
        """Background trace must be index 0 so it renders behind foreground."""
        fig = plot_pca(self.adata, color="cell_type", background="batch == 'X'")
        self.assertEqual(
            getattr(fig.data[0], "legendgroup", None),
            "background",
            "Background trace should be the first trace in the figure",
        )

    def test_background_name_label(self):
        fig = plot_pca(self.adata, color="cell_type", background="cell_type == 'A'")
        bg = self._bg_traces(fig)
        self.assertEqual(bg[0].name, "Background")

    # ------------------------------------------------------------------
    # background with multiple subplots
    # ------------------------------------------------------------------
    def test_background_multi_subplot(self):
        fig = plot_pca(
            self.adata, color="cell_type", PCs=[0, 1, 2], background="batch == 'Y'"
        )
        bg = self._bg_traces(fig)
        # Three PC pairs → three background traces (one per subplot)
        self.assertEqual(len(bg), 3)

    def test_background_multi_subplot_only_one_legend_entry(self):
        """The Background legend entry should appear exactly once."""
        fig = plot_pca(
            self.adata, color="cell_type", PCs=[0, 1, 2], background="batch == 'Y'"
        )
        bg = self._bg_traces(fig)
        visible_entries = [t for t in bg if t.showlegend]
        self.assertEqual(len(visible_entries), 1)

    def test_background_multi_subplot_total_points(self):
        bg_query = "batch == 'Y'"
        n_bg = len(self.adata.obs.query(bg_query))
        n_fg = len(self.adata.obs) - n_bg
        n_pairs = 3  # PCs=[0,1,2] → (0,1),(0,2),(1,2)
        fig = plot_pca(
            self.adata, color="cell_type", PCs=[0, 1, 2], background=bg_query
        )
        bg_points = sum(len(t.x) for t in self._bg_traces(fig))
        fg_points = sum(len(t.x) for t in self._fg_traces(fig))
        self.assertEqual(bg_points, n_bg * n_pairs)
        self.assertEqual(fg_points, n_fg * n_pairs)

    # ------------------------------------------------------------------
    # edge cases
    # ------------------------------------------------------------------
    def test_background_selects_all_cells(self):
        """When every cell matches background, foreground traces are empty."""
        fig = plot_pca(
            self.adata, color="cell_type", background="cell_type == cell_type"
        )
        bg = self._bg_traces(fig)
        self.assertEqual(sum(len(t.x) for t in bg), len(self.adata.obs_names))
        fg = self._fg_traces(fig)
        self.assertEqual(sum(len(t.x) for t in fg), 0)

    def test_background_selects_no_cells(self):
        """When no cell matches background query, result equals no-background case."""
        fig_bg = plot_pca(
            self.adata, color="cell_type", background="batch == 'NONEXISTENT'"
        )
        fig_no_bg = plot_pca(self.adata, color="cell_type")
        # Both figures should carry the same number of traces
        self.assertEqual(len(fig_bg.data), len(fig_no_bg.data))

    def test_background_explicit_pc_pairs(self):
        """background works with explicit List[List[int]] PC specification."""
        fig = plot_pca(
            self.adata,
            color="cell_type",
            PCs=[[0, 1], [2, 3]],
            background="cell_type == 'B'",
        )
        bg = self._bg_traces(fig)
        self.assertEqual(len(bg), 2)  # one per subplot


if __name__ == "__main__":
    unittest.main()
