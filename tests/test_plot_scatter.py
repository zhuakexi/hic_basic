"""
Unit tests for hic_basic.plot.scatter module.

Tests for smooth_with_ci(), lines(), and line() functions with new parameters:
- skip_smooth: bypass LOESS and plot raw y values
- error_type: "ci" for pooled CI, "sem" for per-bin SEM
- sem_col: column with pre-computed SEM values
"""
import os
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.stats import norm

# Add project to path
sys.path.insert(0, "/share/home/ychi/dev/hic_basic")

from hic_basic.plot.scatter import smooth_with_ci, lines, line, add_error_band


class TestSmoothWithCI(unittest.TestCase):
    """Test smooth_with_ci function with new parameters."""

    def setUp(self):
        """Generate synthetic data for testing."""
        np.random.seed(42)
        # Create synthetic data: 20 points with linear trend + noise
        n_points = 20
        self.df = pd.DataFrame({
            "x": np.arange(n_points),
            "y": np.arange(n_points) * 0.1 + 1.0 + np.random.randn(n_points) * 0.05
        })
        
        # Create synthetic data with SEM column
        self.df_with_sem = pd.DataFrame({
            "x": np.arange(n_points),
            "y": np.arange(n_points) * 0.1 + 1.0 + np.random.randn(n_points) * 0.05,
            "sem": np.abs(np.random.randn(n_points)) * 0.02 + 0.01  # varying SEM
        })
        
        # Output directory for visual tests
        self.outdir = Path(os.path.dirname(__file__)) / "output" / "scatter"
        self.outdir.mkdir(parents=True, exist_ok=True)

    def test_default_behavior(self):
        """Test default behavior: LOESS smoothing + pooled CI."""
        result = smooth_with_ci(self.df, x="x", y="y")
        
        # Check output structure
        self.assertIn("x", result.columns)
        self.assertIn("y", result.columns)
        self.assertIn("y_smooth", result.columns)
        self.assertIn("y_upper", result.columns)
        self.assertIn("y_lower", result.columns)
        
        # Check that smoothing was applied (y_smooth != y in general)
        # For noisy data, smoothed should differ from raw
        self.assertFalse(np.allclose(result["y_smooth"], result["y"]))
        
        # Check that CI band exists (upper > lower)
        self.assertTrue(np.all(result["y_upper"] > result["y_lower"]))
        
        # Check CI width is constant (pooled residual)
        ci_width = result["y_upper"] - result["y_lower"]
        self.assertTrue(np.allclose(ci_width, ci_width.iloc[0]))

    def test_skip_smooth_true(self):
        """Test skip_smooth=True: raw y values without LOESS."""
        result = smooth_with_ci(self.df, x="x", y="y", skip_smooth=True)
        
        # y_smooth should equal y (no smoothing)
        self.assertTrue(np.allclose(result["y_smooth"], result["y"]))
        
        # When skip_smooth=True with error_type="ci" (default),
        # residuals = y - y_smooth = 0, so sigma = 0
        # CI band has zero width (upper == lower == y_smooth)
        self.assertTrue(np.allclose(result["y_upper"], result["y_lower"]))

    def test_error_type_sem(self):
        """Test error_type="sem": per-bin SEM ribbon."""
        result = smooth_with_ci(self.df_with_sem, x="x", y="y", 
                                error_type="sem", sem_col="sem")
        
        # Check y_sem column exists
        self.assertIn("y_sem", result.columns)
        
        # Check SEM values match input
        self.assertTrue(np.allclose(result["y_sem"], self.df_with_sem["sem"].values))
        
        # Check error bounds use SEM (varying width)
        # upper - lower = 2 * z * sem, so varies with sem
        error_width = result["y_upper"] - result["y_lower"]
        # Width should vary because SEM varies
        self.assertFalse(np.allclose(error_width, error_width.iloc[0]))

    def test_skip_smooth_and_sem_combined(self):
        """Test skip_smooth=True + error_type="sem": raw values + SEM."""
        result = smooth_with_ci(self.df_with_sem, x="x", y="y", 
                                skip_smooth=True, error_type="sem", sem_col="sem")
        
        # Raw y values
        self.assertTrue(np.allclose(result["y_smooth"], result["y"]))
        
        # SEM ribbon source values
        self.assertIn("y_sem", result.columns)
        self.assertTrue(np.allclose(result["y_sem"], self.df_with_sem["sem"].values))

    def test_error_type_sem_requires_sem_col(self):
        """Test that error_type="sem" requires sem_col parameter."""
        with self.assertRaises(ValueError) as ctx:
            smooth_with_ci(self.df, x="x", y="y", error_type="sem")
        self.assertIn("sem_col must be specified", str(ctx.exception))

    def test_error_type_sem_requires_valid_column(self):
        """Test that error_type="sem" requires sem_col in DataFrame."""
        with self.assertRaises(ValueError) as ctx:
            smooth_with_ci(self.df, x="x", y="y", error_type="sem", sem_col="nonexistent")
        self.assertIn("must exist in the DataFrame", str(ctx.exception))

    def test_invalid_error_type(self):
        """Test that invalid error_type raises error."""
        with self.assertRaises(ValueError) as ctx:
            smooth_with_ci(self.df, x="x", y="y", error_type="invalid")
        self.assertIn("error_type must be 'ci' or 'sem'", str(ctx.exception))


class TestLines(unittest.TestCase):
    """Test lines function with new parameters."""

    def setUp(self):
        """Generate synthetic data for testing."""
        np.random.seed(42)
        n_points = 20
        
        # Single line data
        self.df_single = pd.DataFrame({
            "x": np.arange(n_points),
            "y1": np.arange(n_points) * 0.1 + 1.0 + np.random.randn(n_points) * 0.05
        })
        
        # Multiple lines data
        self.df_multi = pd.DataFrame({
            "x": np.arange(n_points),
            "y1": np.arange(n_points) * 0.1 + 1.0 + np.random.randn(n_points) * 0.05,
            "y2": np.arange(n_points) * 0.1 + 0.5 + np.random.randn(n_points) * 0.05
        })
        
        # Data with SEM column
        self.df_with_sem = pd.DataFrame({
            "x": np.arange(n_points),
            "y": np.arange(n_points) * 0.1 + 1.0 + np.random.randn(n_points) * 0.05,
            "sem": np.abs(np.random.randn(n_points)) * 0.02 + 0.01
        })

        # Multiple lines with dedicated SEM columns
        self.df_multi_sem = pd.DataFrame({
            "x": np.arange(n_points),
            "y1": np.arange(n_points) * 0.1 + 1.0 + np.random.randn(n_points) * 0.05,
            "y2": np.arange(n_points) * 0.08 + 0.6 + np.random.randn(n_points) * 0.05,
            "sem1": np.abs(np.random.randn(n_points)) * 0.01 + 0.01,
            "sem2": np.abs(np.random.randn(n_points)) * 0.015 + 0.015,
            "sem_shared": np.abs(np.random.randn(n_points)) * 0.02 + 0.01,
        })
        
        # Output directory
        self.outdir = Path(os.path.dirname(__file__)) / "output" / "scatter"
        self.outdir.mkdir(parents=True, exist_ok=True)

    def test_lines_basic(self):
        """Test basic lines() call."""
        fig = lines(self.df_single, x="x", y="y1")
        
        self.assertIsInstance(fig, go.Figure)
        # Should have 1 trace (line)
        self.assertEqual(len(fig.data), 1)

    def test_lines_with_error(self):
        """Test lines() with with_error=True (default CI band)."""
        fig = lines(self.df_single, x="x", y="y1", with_error=True)
        
        self.assertIsInstance(fig, go.Figure)
        # Should have 3 traces: line + 2 for CI band (upper + lower)
        self.assertEqual(len(fig.data), 3)

    def test_lines_skip_smooth(self):
        """Test lines() with skip_smooth=True."""
        fig = lines(self.df_single, x="x", y="y1", skip_smooth=True)
        
        self.assertIsInstance(fig, go.Figure)
        # Check that y values in trace match raw data
        trace = fig.data[0]
        self.assertTrue(np.allclose(trace.y, self.df_single["y1"].values))

    def test_lines_error_type_sem(self):
        """Test lines() with error_type="sem"."""
        fig = lines(self.df_with_sem, x="x", y="y", 
                   with_error=True, error_type="sem", sem_col="sem")
        
        self.assertIsInstance(fig, go.Figure)
        # Should have 3 traces: line + 2 ribbon boundary traces
        self.assertEqual(len(fig.data), 3)

        # SEM is rendered as a ribbon, not error_y on the main trace
        main_trace = fig.data[0]
        self.assertNotIn("error_y", main_trace.to_plotly_json())
        self.assertEqual(fig.data[2].fill, "tonexty")

    def test_lines_combined_skip_and_sem(self):
        """Test lines() with skip_smooth=True + error_type="sem"."""
        fig = lines(self.df_with_sem, x="x", y="y",
                   with_error=True, skip_smooth=True, 
                   error_type="sem", sem_col="sem")
        
        self.assertIsInstance(fig, go.Figure)
        # Raw values in line with SEM ribbon
        line_trace = fig.data[0]
        self.assertTrue(np.allclose(line_trace.y, self.df_with_sem["y"].values))
        self.assertEqual(len(fig.data), 3)
        self.assertNotIn("error_y", line_trace.to_plotly_json())
        self.assertEqual(fig.data[2].fill, "tonexty")

    def test_lines_multiple_y(self):
        """Test lines() with multiple y columns."""
        fig = lines(self.df_multi, x="x", y=["y1", "y2"])
        
        self.assertIsInstance(fig, go.Figure)
        # Should have 2 traces (2 lines)
        self.assertEqual(len(fig.data), 2)

    def test_lines_multiple_y_with_shared_sem_col(self):
        """Test lines() with multiple y columns sharing one SEM column."""
        fig = lines(
            self.df_multi_sem,
            x="x",
            y=["y1", "y2"],
            with_error=True,
            skip_smooth=True,
            error_type="sem",
            sem_col="sem_shared",
        )

        self.assertIsInstance(fig, go.Figure)
        self.assertEqual(len(fig.data), 6)
        self.assertEqual(fig.data[2].fill, "tonexty")
        self.assertEqual(fig.data[5].fill, "tonexty")

    def test_lines_multiple_y_with_per_line_sem_cols(self):
        """Test lines() with one SEM column per y column."""
        fig = lines(
            self.df_multi_sem,
            x="x",
            y=["y1", "y2"],
            with_error=True,
            skip_smooth=True,
            error_type="sem",
            sem_col=["sem1", "sem2"],
        )

        self.assertIsInstance(fig, go.Figure)
        self.assertEqual(len(fig.data), 6)

        y1_width = np.array(fig.data[1].y) - np.array(fig.data[2].y)
        y2_width = np.array(fig.data[4].y) - np.array(fig.data[5].y)
        z = norm.isf((1 - 0.95) / 2)
        expected_y1_width = 2 * z * self.df_multi_sem["sem1"].values
        expected_y2_width = 2 * z * self.df_multi_sem["sem2"].values

        self.assertTrue(np.allclose(y1_width, expected_y1_width))
        self.assertTrue(np.allclose(y2_width, expected_y2_width))

    def test_lines_multiple_y_sem_col_length_mismatch(self):
        """Test lines() rejects mismatched y/sem_col list lengths."""
        with self.assertRaises(ValueError) as ctx:
            lines(
                self.df_multi_sem,
                x="x",
                y=["y1", "y2"],
                with_error=True,
                error_type="sem",
                sem_col=["sem1"],
            )
        self.assertIn("same length as y", str(ctx.exception))

    def test_lines_single_y_rejects_sem_col_list(self):
        """Test lines() rejects sem_col lists when y is a single column name."""
        with self.assertRaises(ValueError) as ctx:
            lines(
                self.df_with_sem,
                x="x",
                y="y",
                with_error=True,
                error_type="sem",
                sem_col=["sem"],
            )
        self.assertIn("single column name", str(ctx.exception))

    def test_lines_visual_output(self):
        """Visual test: generate PNG for manual inspection."""
        out_png = self.outdir / "test_lines_visual.png"
        
        fig = lines(self.df_with_sem, x="x", y="y",
                   with_error=True, skip_smooth=True,
                   error_type="sem", sem_col="sem",
                  title="Test: Raw + SEM ribbon")
        
        fig.write_image(str(out_png), width=600, height=400)
        self.assertTrue(out_png.exists())


class TestLine(unittest.TestCase):
    """Test line function with new parameters."""

    def setUp(self):
        """Generate synthetic data for testing."""
        np.random.seed(42)
        n_points = 20
        
        self.df = pd.DataFrame({
            "x": np.arange(n_points),
            "y": np.arange(n_points) * 0.1 + 1.0 + np.random.randn(n_points) * 0.05,
            "group": ["A"] * 10 + ["B"] * 10,
            "sem": np.abs(np.random.randn(n_points)) * 0.02 + 0.01
        })
        
        self.outdir = Path(os.path.dirname(__file__)) / "output" / "scatter"
        self.outdir.mkdir(parents=True, exist_ok=True)

    def test_line_basic(self):
        """Test basic line() call."""
        fig = line(self.df, x="x", y="y")
        
        self.assertIsInstance(fig, go.Figure)
        # Should have 2 traces: line + scatter points
        self.assertGreaterEqual(len(fig.data), 2)

    def test_line_with_error(self):
        """Test line() with with_error=True."""
        fig = line(self.df, x="x", y="y", with_error=True)
        
        self.assertIsInstance(fig, go.Figure)
        # Should have at least 4 traces: line + points + CI band (upper + lower)
        self.assertGreaterEqual(len(fig.data), 4)

    def test_line_skip_smooth(self):
        """Test line() with skip_smooth=True."""
        fig = line(self.df, x="x", y="y", skip_smooth=True)
        
        self.assertIsInstance(fig, go.Figure)
        # Line trace should have raw values
        line_trace = fig.data[0]
        self.assertTrue(np.allclose(line_trace.y, self.df["y"].values))

    def test_line_error_type_sem(self):
        """Test line() with error_type="sem"."""
        fig = line(self.df, x="x", y="y",
                  with_error=True, error_type="sem", sem_col="sem")
        
        self.assertIsInstance(fig, go.Figure)
        # Main line trace should stay clean; SEM is rendered as a ribbon
        main_trace = fig.data[0]
        self.assertNotIn("error_y", main_trace.to_plotly_json())
        self.assertGreaterEqual(len(fig.data), 4)
        self.assertEqual(fig.data[-1].fill, "tonexty")

    def test_line_rejects_sem_col_list(self):
        """Test line() rejects sem_col lists."""
        with self.assertRaises(ValueError) as ctx:
            line(
                self.df,
                x="x",
                y="y",
                with_error=True,
                error_type="sem",
                sem_col=["sem"],
            )
        self.assertIn("single column name", str(ctx.exception))

    def test_line_with_color(self):
        """Test line() with color parameter for point coloring."""
        fig = line(self.df, x="x", y="y", color="group")
        
        self.assertIsInstance(fig, go.Figure)
        # Should have multiple traces for different color groups

    def test_line_visual_output(self):
        """Visual test: generate PNG for manual inspection."""
        out_png = self.outdir / "test_line_visual.png"
        
        fig = line(self.df, x="x", y="y",
                  with_error=True, skip_smooth=True,
                  error_type="sem", sem_col="sem",
                  title="Test line(): Raw + SEM ribbon")
        
        fig.write_image(str(out_png), width=600, height=400)
        self.assertTrue(out_png.exists())


class TestAddErrorBand(unittest.TestCase):
    """Test add_error_band function."""

    def setUp(self):
        """Generate test data for error band."""
        np.random.seed(42)
        n_points = 10
        
        # CI band data
        self.df_ci = pd.DataFrame({
            "x": np.arange(n_points),
            "y_upper": np.arange(n_points) * 0.1 + 1.2,
            "y_lower": np.arange(n_points) * 0.1 + 0.8
        })
        
        # SEM data (needs y_upper, y_lower from pre-computed SEM bounds)
        self.df_sem = pd.DataFrame({
            "x": np.arange(n_points),
            "y_smooth": np.arange(n_points) * 0.1 + 1.0,
            "y_sem": np.abs(np.random.randn(n_points)) * 0.05 + 0.02
        })
        self.df_sem["y_upper"] = self.df_sem["y_smooth"] + self.df_sem["y_sem"]
        self.df_sem["y_lower"] = self.df_sem["y_smooth"] - self.df_sem["y_sem"]

    def test_add_error_band_ci(self):
        """Test add_error_band with CI band."""
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=self.df_ci["x"], y=np.arange(10) * 0.1 + 1.0))
        
        result = add_error_band(fig, self.df_ci, error_type="ci")
        
        # Should add 2 traces (upper + lower)
        self.assertEqual(len(result.data), 3)  # 1 original + 2 band traces

    def test_add_error_band_sem(self):
        """Test add_error_band with SEM ribbon."""
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=self.df_sem["x"], y=self.df_sem["y_smooth"]))
        
        result = add_error_band(fig, self.df_sem, error_type="sem")
        
        # Should add 2 traces (upper + lower ribbon bounds)
        self.assertEqual(len(result.data), 3)
        self.assertEqual(result.data[2].fill, "tonexty")


if __name__ == "__main__":
    unittest.main()
