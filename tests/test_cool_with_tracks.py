"""
Functional tests for plot_cool_with_tracks function.

This module tests the new plot_cool_with_tracks function for:
1. Single cooler heatmap with tracks
2. Two cooler tiling comparison with tracks

Output PNGs are saved to tests/output/cool_with_tracks/ for visual inspection.
"""
import os
import sys
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
import cooler

sys.path.insert(0, "/share/home/ychi/dev/hic_basic")

from hic_basic.plot.hic import plot_cool_with_tracks
from hic_basic.plot.track import BigWigTrack, GeneTrack, EigenTrack, MarkTrack, HeatmapTrack


def create_synthetic_cool_file(output_path, resolution=10000, n_bins=50, sample_id=1):
    """
    Create a synthetic cooler file for testing.
    
    Parameters
    ----------
    output_path : Path
        Output cooler file path
    resolution : int
        Bin size in base pairs
    n_bins : int
        Number of bins
    sample_id : int
        Sample identifier (1 or 2) to create visibly different matrices
    """
    # Create bins
    chrom = "chr1"
    bins = pd.DataFrame({
        "chrom": [chrom] * n_bins,
        "start": np.arange(0, n_bins * resolution, resolution),
        "end": np.arange(resolution, (n_bins + 1) * resolution, resolution)
    })
    
    # Create synthetic contact matrix with distinct structure per sample
    if sample_id == 1:
        np.random.seed(42)
        base_lam = 5
        loop_strength = 50
        compartment_ratio = 1.5
    else:
        np.random.seed(123)
        base_lam = 15  # Higher baseline for sample 2
        loop_strength = 100  # Stronger loops for sample 2
        compartment_ratio = 0.7  # Weaker compartments for sample 2
    
    mat = np.random.poisson(lam=base_lam, size=(n_bins, n_bins)).astype(float)
    
    # Add distance decay
    for i in range(n_bins):
        for j in range(n_bins):
            dist = abs(i - j)
            mat[i, j] *= np.exp(-dist / 10)
    
    # Add compartments with different patterns per sample
    for i in range(n_bins):
        for j in range(n_bins):
            if sample_id == 1:
                # Sample 1: strong checkerboard
                if (i // 10 + j // 10) % 2 == 0:
                    mat[i, j] *= compartment_ratio
                else:
                    mat[i, j] *= 0.5
            else:
                # Sample 2: different compartment pattern (inverse)
                if (i // 10 + j // 10) % 2 == 0:
                    mat[i, j] *= 0.5
                else:
                    mat[i, j] *= compartment_ratio
    
    # Add sample-specific strong loops at DIFFERENT positions
    if sample_id == 1:
        # Sample 1: loops in top-left quadrant (visible in lower triangle)
        mat[5, 15] = loop_strength
        mat[8, 12] = loop_strength * 0.8
        mat[20, 25] = loop_strength * 0.6
    else:
        # Sample 2: loops in bottom-right quadrant (visible in upper triangle)
        mat[30, 40] = loop_strength * 1.5
        mat[35, 45] = loop_strength * 1.2
        mat[25, 35] = loop_strength
    
    # Create cooler - only upper triangle pixels
    # cooler requires bin1_id <= bin2_id when ordered=True
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


def create_synthetic_bigwig_file(output_path, n_bins=50, resolution=10000):
    """
    Create a synthetic bigWig file for testing.
    
    Note: This creates a simple bedGraph-like file that can be used
    with pyBigWig for testing.
    """
    try:
        import pyBigWig
    except ImportError:
        raise unittest.SkipTest("pyBigWig not available, skipping bigWig creation")
    
    # Create signal with some peaks
    positions = np.arange(0, n_bins * resolution, resolution)
    signal = np.random.exponential(scale=10, size=n_bins)
    
    # Add some peaks
    signal[10:15] += 50
    signal[30:35] += 80
    
    # Write bigWig
    header = [("chr1", n_bins * resolution)]
    with pyBigWig.open(str(output_path), "w") as bw:
        bw.addHeader(header)
        
        # Write in chunks
        for i in range(n_bins):
            bw.addEntries(
                ["chr1"] * 1,
                [int(positions[i])],
                ends=[int(positions[i] + resolution)],
                values=[float(signal[i])]
            )
    
    return output_path


def create_synthetic_eigenvector_file(output_path, n_bins=50, resolution=10000):
    """
    Create a synthetic eigenvector bedGraph file for testing.
    """
    positions = np.arange(0, n_bins * resolution, resolution)
    
    # Create compartment pattern (alternating A/B)
    eig = np.zeros(n_bins)
    for i in range(n_bins):
        if (i // 10) % 2 == 0:
            eig[i] = np.random.normal(loc=0.5, scale=0.2)  # A compartment
        else:
            eig[i] = np.random.normal(loc=-0.5, scale=0.2)  # B compartment
    
    # Write bedGraph (no header comment, use tab separator)
    with open(output_path, "w") as f:
        f.write("chrom\tstart\tend\tE1\n")
        for i in range(n_bins):
            f.write(f"chr1\t{int(positions[i])}\t{int(positions[i] + resolution)}\t{eig[i]:.4f}\n")
    
    return output_path


class TestPlotCoolWithTracksSingle(unittest.TestCase):
    """Test plot_cool_with_tracks with single cooler."""
    
    @classmethod
    def setUpClass(cls):
        """Create synthetic test data."""
        cls.test_data_dir = Path(os.path.dirname(__file__)) / "data"
        cls.test_data_dir.mkdir(parents=True, exist_ok=True)
        
        cls.synthetic_cool = cls.test_data_dir / "synthetic_cool_with_tracks.cool"
        cls.synthetic_eigs = cls.test_data_dir / "synthetic_eigs_tracks.bedgraph"
        
        # Create synthetic cooler
        create_synthetic_cool_file(cls.synthetic_cool, resolution=10000, n_bins=50)
        
        # Create synthetic eigenvector file
        create_synthetic_eigenvector_file(cls.synthetic_eigs, n_bins=50, resolution=10000)
        
        cls.outdir = Path(os.path.dirname(__file__)) / "output" / "cool_with_tracks"
        cls.outdir.mkdir(parents=True, exist_ok=True)
    
    def test_single_cooler_no_tracks(self):
        """Test single cooler without tracks."""
        out_png = self.outdir / "single_cooler_no_tracks.png"
        
        fig = plot_cool_with_tracks(
            coolp1=str(self.synthetic_cool),
            region="chr1",
            title="Single Cooler (No Tracks)",
            balance=False,
            cmap="fall"
        )
        
        fig.write_image(str(out_png), width=700, height=600)
        self.assertTrue(out_png.exists())
        print(f"Output: {out_png}")
    
    def test_single_cooler_with_eigen_track(self):
        """Test single cooler with eigenvector track."""
        out_png = self.outdir / "single_cooler_with_eigen_track.png"
        
        tracks = [
            EigenTrack(str(self.synthetic_eigs), eigen_col="E1")
        ]
        
        fig = plot_cool_with_tracks(
            coolp1=str(self.synthetic_cool),
            tracks=tracks,
            region="chr1",
            title="Single Cooler + Eigen Track",
            balance=False,
            track_heights=[2]
        )
        
        fig.write_image(str(out_png), width=700, height=700)
        self.assertTrue(out_png.exists())
        print(f"Output: {out_png}")
    
    def test_single_cooler_with_custom_vmax(self):
        """Test single cooler with custom vmax."""
        out_png = self.outdir / "single_cooler_custom_vmax.png"
        
        fig = plot_cool_with_tracks(
            coolp1=str(self.synthetic_cool),
            region="chr1",
            title="Single Cooler (vmax=20)",
            balance=False,
            vmax=20,
            cmap="fall"
        )
        
        fig.write_image(str(out_png), width=700, height=600)
        self.assertTrue(out_png.exists())
        print(f"Output: {out_png}")
    
    def test_single_cooler_ignore_diags(self):
        """Test single cooler with ignore_diags=True."""
        out_png = self.outdir / "single_cooler_ignore_diags.png"
        
        fig = plot_cool_with_tracks(
            coolp1=str(self.synthetic_cool),
            region="chr1",
            title="Single Cooler (ignore_diags=True)",
            balance=False,
            ignore_diags=True
        )
        
        fig.write_image(str(out_png), width=700, height=600)
        self.assertTrue(out_png.exists())
        print(f"Output: {out_png}")


class TestPlotCoolWithTracksTiling(unittest.TestCase):
    """Test plot_cool_with_tracks with two coolers for tiling comparison."""
    
    @classmethod
    def setUpClass(cls):
        """Create synthetic test data for tiling comparison."""
        cls.test_data_dir = Path(os.path.dirname(__file__)) / "data"
        cls.test_data_dir.mkdir(parents=True, exist_ok=True)
        
        # Create two different synthetic coolers with VISIBLE differences
        cls.synthetic_cool1 = cls.test_data_dir / "synthetic_cool1_tiling.cool"
        cls.synthetic_cool2 = cls.test_data_dir / "synthetic_cool2_tiling.cool"
        cls.synthetic_eigs1 = cls.test_data_dir / "synthetic_eigs1_tiling.bedgraph"
        cls.synthetic_eigs2 = cls.test_data_dir / "synthetic_eigs2_tiling.bedgraph"
        cls.synthetic_bw1 = cls.test_data_dir / "synthetic_signal1.bw"
        cls.synthetic_bw2 = cls.test_data_dir / "synthetic_signal2.bw"
        cls.synthetic_gtf = cls.test_data_dir / "synthetic_genes.gtf"
        
        # Create two DIFFERENT coolers (sample_id=1 vs sample_id=2)
        # Sample 1: loops in top-left (lower triangle)
        # Sample 2: loops in bottom-right (upper triangle)
        create_synthetic_cool_file(cls.synthetic_cool1, resolution=10000, n_bins=50, sample_id=1)
        create_synthetic_cool_file(cls.synthetic_cool2, resolution=10000, n_bins=50, sample_id=2)
        
        # Create two different eigenvector files
        create_synthetic_eigenvector_file(cls.synthetic_eigs1, n_bins=50, resolution=10000)
        
        # Create second eigenvector with opposite pattern
        positions = np.arange(0, 50 * 10000, 10000)
        eig = np.zeros(50)
        for i in range(50):
            if (i // 10) % 2 == 0:
                eig[i] = np.random.normal(loc=-0.5, scale=0.2)  # Opposite: B compartment
            else:
                eig[i] = np.random.normal(loc=0.5, scale=0.2)  # Opposite: A compartment
        
        with open(cls.synthetic_eigs2, "w") as f:
            f.write("chrom\tstart\tend\tE1\n")
            for i in range(50):
                f.write(f"chr1\t{int(positions[i])}\t{int(positions[i] + 10000)}\t{eig[i]:.4f}\n")
        
        # Create synthetic bigWig files (signal tracks)
        try:
            import pyBigWig
            # Signal 1: peaks at positions 5-15 (matching sample 1 loops)
            header = [("chr1", 50 * 10000)]
            with pyBigWig.open(str(cls.synthetic_bw1), "w") as bw:
                bw.addHeader(header)
                for i in range(50):
                    signal = 10 if 5 <= i <= 15 else 2
                    bw.addEntries(
                        ["chr1"],
                        [int(positions[i])],
                        ends=[int(positions[i] + 10000)],
                        values=[float(signal)]
                    )
            
            # Signal 2: peaks at positions 30-45 (matching sample 2 loops)
            with pyBigWig.open(str(cls.synthetic_bw2), "w") as bw:
                bw.addHeader(header)
                for i in range(50):
                    signal = 20 if 30 <= i <= 45 else 2
                    bw.addEntries(
                        ["chr1"],
                        [int(positions[i])],
                        ends=[int(positions[i] + 10000)],
                        values=[float(signal)]
                    )
        except ImportError:
            cls.synthetic_bw1 = None
            cls.synthetic_bw2 = None
        
        # Create synthetic GTF file for gene track
        with open(cls.synthetic_gtf, "w") as f:
            f.write("# Gene annotations for testing\n")
            # Gene 1: 5000-25000 (near sample 1 loops)
            f.write("chr1\ttest\tgene\t5000\t25000\t.\t+\t.\tgene_id \"Gene1\"; gene_name \"GeneA\";\n")
            f.write("chr1\ttest\texon\t5000\t10000\t.\t+\t.\tgene_id \"Gene1\"; gene_name \"GeneA\";\n")
            f.write("chr1\ttest\texon\t15000\t25000\t.\t+\t.\tgene_id \"Gene1\"; gene_name \"GeneA\";\n")
            # Gene 2: 30000-45000 (near sample 2 loops)
            f.write("chr1\ttest\tgene\t30000\t45000\t.\t-\t.\tgene_id \"Gene2\"; gene_name \"GeneB\";\n")
            f.write("chr1\ttest\texon\t30000\t35000\t.\t-\t.\tgene_id \"Gene2\"; gene_name \"GeneB\";\n")
            f.write("chr1\ttest\texon\t40000\t45000\t.\t-\t.\tgene_id \"Gene2\"; gene_name \"GeneB\";\n")
        
        cls.outdir = Path(os.path.dirname(__file__)) / "output" / "cool_with_tracks"
        cls.outdir.mkdir(parents=True, exist_ok=True)
    
    def test_tiling_no_tracks(self):
        """Test tiling comparison without tracks."""
        out_png = self.outdir / "tiling_no_tracks.png"
        
        fig = plot_cool_with_tracks(
            coolp1=str(self.synthetic_cool1),
            coolp2=str(self.synthetic_cool2),
            region="chr1",
            title="Tiling Comparison (No Tracks)",
            balance=False,
            adjust_brightness=True
        )
        
        fig.write_image(str(out_png), width=700, height=600)
        self.assertTrue(out_png.exists())
        print(f"Output: {out_png}")
    
    def test_tiling_with_eigen_track(self):
        """Test tiling comparison with eigenvector track."""
        out_png = self.outdir / "tiling_with_eigen_track.png"
        
        tracks = [
            EigenTrack(str(self.synthetic_eigs1), eigen_col="E1")
        ]
        
        fig = plot_cool_with_tracks(
            coolp1=str(self.synthetic_cool1),
            coolp2=str(self.synthetic_cool2),
            tracks=tracks,
            region="chr1",
            title="Tiling + Eigen Track",
            balance=False,
            adjust_brightness=True,
            track_heights=[2]
        )
        
        fig.write_image(str(out_png), width=700, height=700)
        self.assertTrue(out_png.exists())
        print(f"Output: {out_png}")
    
    def test_tiling_with_multiple_tracks(self):
        """Test tiling comparison with multiple tracks."""
        out_png = self.outdir / "tiling_with_multiple_tracks.png"
        
        tracks = [
            EigenTrack(str(self.synthetic_eigs1), eigen_col="E1"),
            EigenTrack(str(self.synthetic_eigs2), eigen_col="E1")
        ]
        
        fig = plot_cool_with_tracks(
            coolp1=str(self.synthetic_cool1),
            coolp2=str(self.synthetic_cool2),
            tracks=tracks,
            region="chr1",
            title="Tiling + 2 Eigen Tracks",
            balance=False,
            adjust_brightness=True,
            track_heights=[2, 2]
        )
        
        fig.write_image(str(out_png), width=700, height=800)
        self.assertTrue(out_png.exists())
        print(f"Output: {out_png}")
    
    def test_tiling_with_all_track_types(self):
        """Test tiling comparison with all track types: eigen, bigwig, gene."""
        out_png = self.outdir / "tiling_with_all_tracks.png"
        
        tracks = [
            EigenTrack(str(self.synthetic_eigs1), eigen_col="E1"),
        ]
        
        # Add bigWig tracks if available
        if self.synthetic_bw1 is not None:
            tracks.append(BigWigTrack(str(self.synthetic_bw1), bin_size=10000, color="#1f77b4"))
        
        # Add gene track
        tracks.append(GeneTrack(str(self.synthetic_gtf), mode="all_transcripts"))
        
        fig = plot_cool_with_tracks(
            coolp1=str(self.synthetic_cool1),
            coolp2=str(self.synthetic_cool2),
            tracks=tracks,
            region="chr1:0-500000",  # Use full region format for BigWigTrack
            title="Tiling + Eigen + BigWig + Gene",
            balance=False,
            adjust_brightness=True,
            track_heights=[2, 2, 3]  # Gene track taller
        )
        
        fig.write_image(str(out_png), width=800, height=1000)
        self.assertTrue(out_png.exists())
        print(f"Output: {out_png}")
    
    def test_tiling_no_brightness_adjust(self):
        """Test tiling comparison without brightness adjustment."""
        out_png = self.outdir / "tiling_no_brightness_adjust.png"
        
        fig = plot_cool_with_tracks(
            coolp1=str(self.synthetic_cool1),
            coolp2=str(self.synthetic_cool2),
            region="chr1",
            title="Tiling (No Brightness Adjust)",
            balance=False,
            adjust_brightness=False
        )
        
        fig.write_image(str(out_png), width=700, height=600)
        self.assertTrue(out_png.exists())
        print(f"Output: {out_png}")
    
    def test_tiling_custom_cmap(self):
        """Test tiling comparison with custom colormap."""
        out_png = self.outdir / "tiling_custom_cmap.png"
        
        fig = plot_cool_with_tracks(
            coolp1=str(self.synthetic_cool1),
            coolp2=str(self.synthetic_cool2),
            region="chr1",
            title="Tiling (RdBu_r)",
            balance=False,
            cmap="RdBu_r",
            adjust_brightness=True
        )
        
        fig.write_image(str(out_png), width=700, height=600)
        self.assertTrue(out_png.exists())
        print(f"Output: {out_png}")


class TestPlotCoolWithTracksVisual(unittest.TestCase):
    """Comprehensive visual tests for plot_cool_with_tracks."""
    
    @classmethod
    def setUpClass(cls):
        """Use existing synthetic test data."""
        cls.test_data_dir = Path(os.path.dirname(__file__)) / "data"
        
        cls.synthetic_cool1 = cls.test_data_dir / "synthetic_cool1_tiling.cool"
        cls.synthetic_cool2 = cls.test_data_dir / "synthetic_cool2_tiling.cool"
        cls.synthetic_eigs1 = cls.test_data_dir / "synthetic_eigs1_tiling.bedgraph"
        cls.synthetic_eigs2 = cls.test_data_dir / "synthetic_eigs2_tiling.bedgraph"
        cls.synthetic_bw1 = cls.test_data_dir / "synthetic_signal1.bw"
        cls.synthetic_bw2 = cls.test_data_dir / "synthetic_signal2.bw"
        cls.synthetic_gtf = cls.test_data_dir / "synthetic_genes.gtf"
        
        cls.outdir = Path(os.path.dirname(__file__)) / "output" / "cool_with_tracks"
        cls.outdir.mkdir(parents=True, exist_ok=True)
    
    def test_comprehensive_single_sample(self):
        """Comprehensive test: single sample with multiple tracks."""
        out_png = self.outdir / "comprehensive_single.png"
        
        tracks = [
            EigenTrack(str(self.synthetic_eigs1), eigen_col="E1"),
        ]
        if self.synthetic_bw1 is not None:
            tracks.append(BigWigTrack(str(self.synthetic_bw1), bin_size=10000, color="#ff7f0e"))
        tracks.append(GeneTrack(str(self.synthetic_gtf), mode="all_transcripts"))
        
        fig = plot_cool_with_tracks(
            coolp1=str(self.synthetic_cool1),
            tracks=tracks,
            region="chr1:0-500000",  # Use full region format for BigWigTrack
            title="Comprehensive: Single + Eigen + BigWig + Gene",
            balance=False,
            heatmap_height=12,
            track_heights=[2, 2, 3],
            cmap="fall",
            vmax=30
        )
        
        fig.write_image(str(out_png), width=800, height=1000)
        self.assertTrue(out_png.exists())
        print(f"Output: {out_png}")
    
    def test_comprehensive_tiling_comparison(self):
        """Comprehensive test: tiling comparison with multiple tracks."""
        out_png = self.outdir / "comprehensive_tiling.png"
        
        tracks = [
            EigenTrack(str(self.synthetic_eigs1), eigen_col="E1"),
        ]
        if self.synthetic_bw1 is not None:
            tracks.append(BigWigTrack(str(self.synthetic_bw1), bin_size=10000, color="#1f77b4"))
        if self.synthetic_bw2 is not None:
            tracks.append(BigWigTrack(str(self.synthetic_bw2), bin_size=10000, color="#2ca02c"))
        tracks.append(GeneTrack(str(self.synthetic_gtf), mode="all_transcripts"))
        
        fig = plot_cool_with_tracks(
            coolp1=str(self.synthetic_cool1),
            coolp2=str(self.synthetic_cool2),
            tracks=tracks,
            region="chr1:0-500000",  # Use full region format for BigWigTrack
            title="Comprehensive: Tiling + Eigen + 2xBigWig + Gene",
            balance=False,
            heatmap_height=12,
            track_heights=[2, 2, 2, 3],
            cmap="RdBu_r",
            adjust_brightness=True
        )
        
        fig.write_image(str(out_png), width=900, height=1100)
        self.assertTrue(out_png.exists())
        print(f"Output: {out_png}")
    
    def test_tiling_visual_emphasis(self):
        """Visual test emphasizing tiling differences between two samples."""
        out_png = self.outdir / "tiling_visual_emphasis.png"
        
        # This test specifically highlights the tiling effect
        # Sample 1 (lower triangle): loops at 5-15, 20-25
        # Sample 2 (upper triangle): loops at 30-40, 35-45
        fig = plot_cool_with_tracks(
            coolp1=str(self.synthetic_cool1),
            coolp2=str(self.synthetic_cool2),
            region="chr1",
            title="Tiling: Lower=S1(loops@5-15,20-25) Upper=S2(loops@30-45)",
            balance=False,
            cmap="YlOrRd",  # Use a colormap that makes loops stand out
            adjust_brightness=False,  # Don't adjust to see raw differences
            vmax=100
        )
        
        fig.write_image(str(out_png), width=700, height=600)
        self.assertTrue(out_png.exists())
        print(f"Output: {out_png}")


class TestMarkTrack(unittest.TestCase):
    """Test MarkTrack functionality for marking pixels on heatmaps."""
    
    @classmethod
    def setUpClass(cls):
        """Use existing synthetic test data."""
        cls.test_data_dir = Path(os.path.dirname(__file__)) / "data"
        cls.synthetic_cool1 = cls.test_data_dir / "synthetic_cool1_tiling.cool"
        cls.synthetic_cool2 = cls.test_data_dir / "synthetic_cool2_tiling.cool"
        cls.synthetic_eigs1 = cls.test_data_dir / "synthetic_eigs1_tiling.bedgraph"
        cls.outdir = Path(os.path.dirname(__file__)) / "output" / "cool_with_tracks"
        cls.outdir.mkdir(parents=True, exist_ok=True)
    
    def test_mark_track_basic(self):
        """Test basic MarkTrack functionality."""
        from plotly.subplots import make_subplots
        from cooler import Cooler
        
        clr = Cooler(str(self.synthetic_cool1))
        
        # Create mark track with a single mark
        mark = MarkTrack(
            marks=[
                {
                    "chrom1": "chr1", "start1": 50000, "end1": 60000,
                    "chrom2": "chr1", "start2": 150000, "end2": 160000,
                    "color": "red", "shape": "circle", "width": 2
                }
            ]
        )
        mark.load_data(clr=clr, region="chr1:0-500000")
        
        # Create figure and render
        fig = make_subplots(rows=1, cols=1)
        mark.render(fig, row=1, col=1)
        
        # Should have added at least one shape
        self.assertGreaterEqual(len(fig.layout.shapes), 1)
    
    def test_mark_track_multiple_shapes(self):
        """Test MarkTrack with multiple marks and shapes."""
        from plotly.subplots import make_subplots
        from cooler import Cooler
        
        clr = Cooler(str(self.synthetic_cool1))
        
        mark = MarkTrack(
            marks=[
                {
                    "chrom1": "chr1", "start1": 50000, "end1": 60000,
                    "chrom2": "chr1", "start2": 150000, "end2": 160000,
                    "color": "red", "shape": "circle", "width": 2
                },
                {
                    "chrom1": "chr1", "start1": 200000, "end1": 210000,
                    "chrom2": "chr1", "start2": 250000, "end2": 260000,
                    "color": "blue", "shape": "rect", "width": 3
                }
            ]
        )
        mark.load_data(clr=clr, region="chr1:0-500000")
        
        fig = make_subplots(rows=1, cols=1)
        mark.render(fig, row=1, col=1)
        
        # Should have added two shapes
        self.assertEqual(len(fig.layout.shapes), 2)
    
    def test_heatmap_with_marks(self):
        """Test HeatmapTrack with marks parameter."""
        from cooler import Cooler
        
        clr = Cooler(str(self.synthetic_cool1))
        
        heatmap = HeatmapTrack(
            cool_path=str(self.synthetic_cool1),
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
        
        from plotly.subplots import make_subplots
        fig = make_subplots(rows=1, cols=1)
        heatmap.render(fig, row=1, col=1)
        
        # Should have heatmap trace + mark shapes
        self.assertGreaterEqual(len(fig.data), 1)
        self.assertGreaterEqual(len(fig.layout.shapes), 1)
    
    def test_plot_cool_with_tracks_marks(self):
        """Test plot_cool_with_tracks with marks parameter."""
        out_png = self.outdir / "with_marks.png"
        
        fig = plot_cool_with_tracks(
            coolp1=str(self.synthetic_cool1),
            region="chr1:0-500000",
            title="With Marks",
            balance=False,
            marks=[
                {
                    "chrom1": "chr1", "start1": 50000, "end1": 60000,
                    "chrom2": "chr1", "start2": 150000, "end2": 160000,
                    "color": "red", "shape": "circle", "width": 2,
                    "label": "Loop"
                },
                {
                    "chrom1": "chr1", "start1": 200000, "end1": 210000,
                    "chrom2": "chr1", "start2": 250000, "end2": 260000,
                    "color": "blue", "shape": "rect", "width": 3
                }
            ]
        )
        
        fig.write_image(str(out_png), width=700, height=700)
        self.assertTrue(out_png.exists())
        print(f"Output: {out_png}")
    
    def test_tiling_with_marks(self):
        """Test tiling comparison with marks."""
        out_png = self.outdir / "tiling_with_marks.png"
        
        fig = plot_cool_with_tracks(
            coolp1=str(self.synthetic_cool1),
            coolp2=str(self.synthetic_cool2),
            region="chr1:0-500000",
            title="Tiling with Marks",
            balance=False,
            adjust_brightness=True,
            marks=[
                {
                    "chrom1": "chr1", "start1": 50000, "end1": 60000,
                    "chrom2": "chr1", "start2": 150000, "end2": 160000,
                    "color": "red", "shape": "circle", "width": 2,
                    "label": "S1 Loop"
                },
                {
                    "chrom1": "chr1", "start1": 300000, "end1": 310000,
                    "chrom2": "chr1", "start2": 400000, "end2": 410000,
                    "color": "blue", "shape": "circle", "width": 2,
                    "label": "S2 Loop"
                }
            ]
        )
        
        fig.write_image(str(out_png), width=800, height=800)
        self.assertTrue(out_png.exists())
        print(f"Output: {out_png}")


if __name__ == "__main__":
    # Run tests
    unittest.main(verbosity=2)
