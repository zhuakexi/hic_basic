"""
Functional tests for MarkTrack - marking pixels on Hi-C heatmaps.

This module tests the MarkTrack class and marks functionality for:
1. Basic mark rendering (circle, rect, line shapes)
2. Multiple marks with different styles
3. Integration with HeatmapTrack
4. Integration with plot_cool_with_tracks
5. Visual verification with synthetic data

Output PNGs are saved to tests/output/mark_track/ for visual inspection.
"""
import os
import sys
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
import cooler

sys.path.insert(0, "/share/home/ychi/dev/hic_basic")

from hic_basic.plot.track import MarkTrack, HeatmapTrack
from hic_basic.plot.hic import plot_cool_with_tracks
from plotly.subplots import make_subplots


def create_synthetic_cool_for_marks(output_path, resolution=10000, n_bins=50):
    """
    Create a synthetic cooler file specifically for testing marks.
    Features:
    - Strong loops at specific positions for visual verification
    - Clear compartment pattern
    """
    chrom = "chr1"
    bins = pd.DataFrame({
        "chrom": [chrom] * n_bins,
        "start": np.arange(0, n_bins * resolution, resolution),
        "end": np.arange(resolution, (n_bins + 1) * resolution, resolution)
    })
    
    # Create matrix with distance decay
    np.random.seed(42)
    mat = np.random.poisson(lam=5, size=(n_bins, n_bins)).astype(float)
    for i in range(n_bins):
        for j in range(n_bins):
            dist = abs(i - j)
            mat[i, j] *= np.exp(-dist / 10)
    
    # Add strong, visually distinct loops at known positions
    # Loop 1: bin 5-15 (50kb-150kb)
    mat[5, 15] = mat[15, 5] = 100
    # Loop 2: bin 20-30 (200kb-300kb)
    mat[20, 30] = mat[30, 20] = 150
    # Loop 3: bin 35-45 (350kb-450kb)
    mat[35, 45] = mat[45, 35] = 200
    
    # Create cooler
    upper_tri_indices = np.triu_indices(n_bins)
    cooler.create_cooler(
        str(output_path),
        bins=bins,
        pixels=pd.DataFrame({
            "bin1_id": upper_tri_indices[0],
            "bin2_id": upper_tri_indices[1],
            "count": mat[upper_tri_indices].astype(int)
        }),
        ordered=True,
        mode='w'
    )
    return output_path


class TestMarkTrackBasic(unittest.TestCase):
    """Test basic MarkTrack functionality."""
    
    @classmethod
    def setUpClass(cls):
        """Create synthetic test data."""
        cls.test_data_dir = Path(os.path.dirname(__file__)) / "data"
        cls.test_data_dir.mkdir(parents=True, exist_ok=True)
        
        cls.synthetic_cool = cls.test_data_dir / "synthetic_mark_test.cool"
        create_synthetic_cool_for_marks(cls.synthetic_cool, resolution=10000, n_bins=50)
        
        cls.outdir = Path(os.path.dirname(__file__)) / "output" / "mark_track"
        cls.outdir.mkdir(parents=True, exist_ok=True)
    
    def test_mark_track_circle_shape(self):
        """Test circle shape mark."""
        from cooler import Cooler
        
        clr = Cooler(str(self.synthetic_cool))
        
        mark = MarkTrack(
            marks=[
                {
                    "chrom1": "chr1", "start1": 50000, "end1": 60000,
                    "chrom2": "chr1", "start2": 150000, "end2": 160000,
                    "color": "red", "shape": "circle", "width": 2
                }
            ],
            default_shape="circle"
        )
        mark.load_data(clr=clr, region="chr1:0-500000")
        
        fig = make_subplots(rows=1, cols=1)
        mark.render(fig, row=1, col=1)
        
        self.assertGreaterEqual(len(fig.layout.shapes), 1)
        self.assertEqual(fig.layout.shapes[0].type, "circle")
    
    def test_mark_track_rect_shape(self):
        """Test rectangle shape mark."""
        from cooler import Cooler
        
        clr = Cooler(str(self.synthetic_cool))
        
        mark = MarkTrack(
            marks=[
                {
                    "chrom1": "chr1", "start1": 200000, "end1": 210000,
                    "chrom2": "chr1", "start2": 300000, "end2": 310000,
                    "color": "blue", "shape": "rect", "width": 3
                }
            ]
        )
        mark.load_data(clr=clr, region="chr1:0-500000")
        
        fig = make_subplots(rows=1, cols=1)
        mark.render(fig, row=1, col=1)
        
        self.assertEqual(len(fig.layout.shapes), 1)
        self.assertEqual(fig.layout.shapes[0].type, "rect")
    
    def test_mark_track_line_shape(self):
        """Test line shape mark."""
        from cooler import Cooler
        
        clr = Cooler(str(self.synthetic_cool))
        
        mark = MarkTrack(
            marks=[
                {
                    "chrom1": "chr1", "start1": 0, "end1": 10000,
                    "chrom2": "chr1", "start2": 0, "end2": 10000,
                    "color": "green", "shape": "line", "width": 2
                }
            ]
        )
        mark.load_data(clr=clr, region="chr1:0-500000")
        
        fig = make_subplots(rows=1, cols=1)
        mark.render(fig, row=1, col=1)
        
        self.assertEqual(len(fig.layout.shapes), 1)
        self.assertEqual(fig.layout.shapes[0].type, "line")
    
    def test_mark_track_multiple_marks(self):
        """Test multiple marks with different styles."""
        from cooler import Cooler
        
        clr = Cooler(str(self.synthetic_cool))
        
        mark = MarkTrack(
            marks=[
                {
                    "chrom1": "chr1", "start1": 50000, "end1": 60000,
                    "chrom2": "chr1", "start2": 150000, "end2": 160000,
                    "color": "red", "shape": "circle", "width": 2
                },
                {
                    "chrom1": "chr1", "start1": 200000, "end1": 210000,
                    "chrom2": "chr1", "start2": 300000, "end2": 310000,
                    "color": "blue", "shape": "rect", "width": 3
                },
                {
                    "chrom1": "chr1", "start1": 350000, "end1": 360000,
                    "chrom2": "chr1", "start2": 450000, "end2": 460000,
                    "color": "green", "shape": "circle", "width": 4,
                    "label": "Strong Loop"
                }
            ]
        )
        mark.load_data(clr=clr, region="chr1:0-500000")
        
        fig = make_subplots(rows=1, cols=1)
        mark.render(fig, row=1, col=1)
        
        self.assertEqual(len(fig.layout.shapes), 3)
        # Check annotations for labeled marks
        self.assertGreaterEqual(len(fig.layout.annotations), 1)


class TestHeatmapWithMarks(unittest.TestCase):
    """Test HeatmapTrack with marks parameter."""
    
    @classmethod
    def setUpClass(cls):
        """Use existing synthetic test data."""
        cls.test_data_dir = Path(os.path.dirname(__file__)) / "data"
        cls.synthetic_cool = cls.test_data_dir / "synthetic_mark_test.cool"
        cls.outdir = Path(os.path.dirname(__file__)) / "output" / "mark_track"
        cls.outdir.mkdir(parents=True, exist_ok=True)
    
    def test_heatmap_single_mark(self):
        """Test HeatmapTrack with single mark."""
        from cooler import Cooler
        
        clr = Cooler(str(self.synthetic_cool))
        
        heatmap = HeatmapTrack(
            cool_path=str(self.synthetic_cool),
            balance=False,
            marks=[
                {
                    "chrom1": "chr1", "start1": 50000, "end1": 60000,
                    "chrom2": "chr1", "start2": 150000, "end2": 160000,
                    "color": "red", "shape": "circle", "width": 2
                }
            ]
        )
        heatmap.load_data(region="chr1:0-500000", clr=clr)
        
        fig = make_subplots(rows=1, cols=1)
        heatmap.render(fig, row=1, col=1)
        
        # Should have heatmap + mark
        self.assertGreaterEqual(len(fig.data), 1)
        self.assertGreaterEqual(len(fig.layout.shapes), 1)
    
    def test_heatmap_multiple_marks(self):
        """Test HeatmapTrack with multiple marks."""
        from cooler import Cooler
        
        clr = Cooler(str(self.synthetic_cool))
        
        heatmap = HeatmapTrack(
            cool_path=str(self.synthetic_cool),
            balance=False,
            marks=[
                {
                    "chrom1": "chr1", "start1": 50000, "end1": 60000,
                    "chrom2": "chr1", "start2": 150000, "end2": 160000,
                    "color": "red", "shape": "circle", "width": 2
                },
                {
                    "chrom1": "chr1", "start1": 200000, "end1": 210000,
                    "chrom2": "chr1", "start2": 300000, "end2": 310000,
                    "color": "blue", "shape": "rect", "width": 3
                }
            ]
        )
        heatmap.load_data(region="chr1:0-500000", clr=clr)
        
        fig = make_subplots(rows=1, cols=1)
        heatmap.render(fig, row=1, col=1)
        
        self.assertGreaterEqual(len(fig.layout.shapes), 2)


class TestMarkTrackVisual(unittest.TestCase):
    """Visual tests for MarkTrack - generates PNGs for manual inspection."""
    
    @classmethod
    def setUpClass(cls):
        """Use existing synthetic test data."""
        cls.test_data_dir = Path(os.path.dirname(__file__)) / "data"
        cls.synthetic_cool = cls.test_data_dir / "synthetic_mark_test.cool"
        cls.outdir = Path(os.path.dirname(__file__)) / "output" / "mark_track"
        cls.outdir.mkdir(parents=True, exist_ok=True)
    
    def test_visual_circle_marks(self):
        """Visual test: circle marks on heatmap."""
        out_png = self.outdir / "visual_circle_marks.png"
        
        fig = plot_cool_with_tracks(
            coolp1=str(self.synthetic_cool),
            region="chr1:0-500000",
            title="Circle Marks Test",
            balance=False,
            cmap="YlOrRd",
            vmax=100,
            marks=[
                {
                    "chrom1": "chr1", "start1": 50000, "end1": 60000,
                    "chrom2": "chr1", "start2": 150000, "end2": 160000,
                    "color": "red", "shape": "circle", "width": 2,
                    "label": "Loop 1"
                },
                {
                    "chrom1": "chr1", "start1": 200000, "end1": 210000,
                    "chrom2": "chr1", "start2": 300000, "end2": 310000,
                    "color": "blue", "shape": "circle", "width": 2,
                    "label": "Loop 2"
                },
                {
                    "chrom1": "chr1", "start1": 350000, "end1": 360000,
                    "chrom2": "chr1", "start2": 450000, "end2": 460000,
                    "color": "green", "shape": "circle", "width": 2,
                    "label": "Loop 3"
                }
            ]
        )
        
        fig.write_image(str(out_png), width=800, height=800)
        self.assertTrue(out_png.exists())
        print(f"Output: {out_png}")
    
    def test_visual_rect_marks(self):
        """Visual test: rectangle marks on heatmap."""
        out_png = self.outdir / "visual_rect_marks.png"
        
        fig = plot_cool_with_tracks(
            coolp1=str(self.synthetic_cool),
            region="chr1:0-500000",
            title="Rectangle Marks Test",
            balance=False,
            cmap="YlOrRd",
            vmax=100,
            marks=[
                {
                    "chrom1": "chr1", "start1": 50000, "end1": 60000,
                    "chrom2": "chr1", "start2": 150000, "end2": 160000,
                    "color": "red", "shape": "rect", "width": 3
                },
                {
                    "chrom1": "chr1", "start1": 200000, "end1": 210000,
                    "chrom2": "chr1", "start2": 300000, "end2": 310000,
                    "color": "blue", "shape": "rect", "width": 3
                }
            ]
        )
        
        fig.write_image(str(out_png), width=800, height=800)
        self.assertTrue(out_png.exists())
        print(f"Output: {out_png}")
    
    def test_visual_mixed_shapes(self):
        """Visual test: mixed circle and rect marks."""
        out_png = self.outdir / "visual_mixed_shapes.png"
        
        fig = plot_cool_with_tracks(
            coolp1=str(self.synthetic_cool),
            region="chr1:0-500000",
            title="Mixed Shapes Test",
            balance=False,
            cmap="YlOrRd",
            vmax=100,
            marks=[
                {
                    "chrom1": "chr1", "start1": 50000, "end1": 60000,
                    "chrom2": "chr1", "start2": 150000, "end2": 160000,
                    "color": "red", "shape": "circle", "width": 2,
                    "label": "Circle"
                },
                {
                    "chrom1": "chr1", "start1": 200000, "end1": 210000,
                    "chrom2": "chr1", "start2": 300000, "end2": 310000,
                    "color": "blue", "shape": "rect", "width": 3,
                    "label": "Rect"
                }
            ]
        )
        
        fig.write_image(str(out_png), width=800, height=800)
        self.assertTrue(out_png.exists())
        print(f"Output: {out_png}")
    
    def test_visual_marks_with_tracks(self):
        """Visual test: marks with additional tracks."""
        out_png = self.outdir / "visual_marks_with_tracks.png"
        
        from hic_basic.plot.track import EigenTrack
        
        # Create synthetic eigenvector file
        eigs_file = self.test_data_dir / "synthetic_mark_eigs.bedgraph"
        positions = np.arange(0, 50 * 10000, 10000)
        eig = np.zeros(50)
        for i in range(50):
            if (i // 10) % 2 == 0:
                eig[i] = np.random.normal(loc=0.5, scale=0.2)
            else:
                eig[i] = np.random.normal(loc=-0.5, scale=0.2)
        
        with open(eigs_file, "w") as f:
            f.write("chrom\tstart\tend\tE1\n")
            for i in range(50):
                f.write(f"chr1\t{int(positions[i])}\t{int(positions[i] + 10000)}\t{eig[i]:.4f}\n")
        
        fig = plot_cool_with_tracks(
            coolp1=str(self.synthetic_cool),
            region="chr1:0-500000",
            title="Marks + Eigen Track",
            balance=False,
            cmap="YlOrRd",
            vmax=100,
            tracks=[
                EigenTrack(str(eigs_file), eigen_col="E1")
            ],
            track_heights=[2],
            marks=[
                {
                    "chrom1": "chr1", "start1": 50000, "end1": 60000,
                    "chrom2": "chr1", "start2": 150000, "end2": 160000,
                    "color": "red", "shape": "circle", "width": 2,
                    "label": "Loop"
                }
            ]
        )
        
        fig.write_image(str(out_png), width=800, height=900)
        self.assertTrue(out_png.exists())
        print(f"Output: {out_png}")
    
    def test_visual_tiling_with_marks(self):
        """Visual test: tiling comparison with marks."""
        # Create second cooler with different loop pattern
        cool2_path = self.test_data_dir / "synthetic_mark_test2.cool"
        create_synthetic_cool_for_marks(cool2_path, resolution=10000, n_bins=50)
        
        out_png = self.outdir / "visual_tiling_with_marks.png"
        
        fig = plot_cool_with_tracks(
            coolp1=str(self.synthetic_cool),
            coolp2=str(cool2_path),
            region="chr1:0-500000",
            title="Tiling + Marks (Lower=S1, Upper=S2)",
            balance=False,
            cmap="RdBu_r",
            adjust_brightness=True,
            marks=[
                {
                    "chrom1": "chr1", "start1": 50000, "end1": 60000,
                    "chrom2": "chr1", "start2": 150000, "end2": 160000,
                    "color": "yellow", "shape": "circle", "width": 3,
                    "label": "Mark S1"
                },
                {
                    "chrom1": "chr1", "start1": 350000, "end1": 360000,
                    "chrom2": "chr1", "start2": 450000, "end2": 460000,
                    "color": "cyan", "shape": "circle", "width": 3,
                    "label": "Mark S2"
                }
            ]
        )
        
        fig.write_image(str(out_png), width=900, height=900)
        self.assertTrue(out_png.exists())
        print(f"Output: {out_png}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
