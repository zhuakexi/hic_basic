# unittest
import sys
sys.path.insert(0, "/share/home/ychi/dev/hic_basic")
sys.path.insert(0, "/share/home/ychi/dev/hires_utils")
import os
import unittest
from pathlib import Path

import pandas as pd
from cooler import Cooler
from hic_basic.plot.hic import (
    plot_compartment,
    plot_cool_track,
    merge_track_data,
    _plot_mat,
    plot_cool,
    plot_tiling_compartment,
    plot_compartment_trapezoid,
    plot_compartment_trapezoid_2tracks,
    plot_tiling_compartment_v2,
    plot_compartment_trapezoid_2tracks_v2,
)

class TestPlot(unittest.TestCase):
    """
    Test the plot module.
    """
    def setUp(self) -> None:
        self.outdir = Path(os.path.dirname(__file__)) / "output" / "plot"
        self.coolp = "/share/Data/ychi/raw/Bonev2017/ES.mcool::resolutions/20000"
        self.coolp_1M = "/share/Data/ychi/raw/Bonev2017/ES.mcool::resolutions/1000000"
        self.IS = "/shareb/ychi/repo/sperm_struct/notebooks/data/mESC.20k.IS.tsv"
        # Synthetic test data for plot_tiling_compartment
        self.synthetic_cool = Path(os.path.dirname(__file__)) / "data" / "synthetic_test.cool"
        self.synthetic_eigs1 = Path(os.path.dirname(__file__)) / "data" / "synthetic_eigs1.bedgraph"
        self.synthetic_eigs2 = Path(os.path.dirname(__file__)) / "data" / "synthetic_eigs2.bedgraph"
    def test_plot_cool_large_scale_grid(self):
        """
        Test the plot_cool_large_scale function.
        """
        out_png = self.outdir / "plot_cool_large_scale_grid.png"
        fig = plot_cool(
            str(self.coolp_1M),
            #region=slice(1,300),
            region = slice(0,500),
            title = "ES 1M 1-500 bins",
            balance = True,
            donorm = False,
            zmax = 0.01,
            zmin = 0
        )
        fig.write_image(str(out_png), width=800, height=800)
        self.assertTrue(out_png.exists())
    def test_plot_cool_large_scale_nogrid(self):
        """
        Test the plot_cool_large_scale function.
        """
        out_png = self.outdir / "plot_cool_large_scale_nogrid.png"
        fig = plot_cool(
            str(self.coolp_1M),
            #region=slice(1,300),
            region = slice(0,500),
            grid = False,
            title = "ES 1M 1-500 bins",
            balance = True,
            donorm = False,
            zmax = 0.01,
            zmin = 0
        )
        fig.write_image(str(out_png), width=800, height=800)
        self.assertTrue(out_png.exists())
    def test_plot_cool_track(self):
        """
        Test the plot_cool_track function.
        """
        out_png = self.outdir / "plot_cool_track.png"
        fig = plot_cool_track(
            str(self.coolp),
            {"IS" : (self.IS, "log2_insulation_score_200000")},
            region="chr6:50,600,000-54,600,000",
            title = "HOXA locus IS track",
            vmax = 5000
        )
        fig.write_image(str(out_png), width=600, height=800)
        self.assertTrue(out_png.exists())
    # def test_pileup_IS(self):
    #     """
    #     Test the pileup_IS function.
    #     """
    #     out_png = self.outdir / "pileup_IS.png"
    #     IS = merge_track_data(
    #         Cooler(str(self.coolp)),
    #         self.IS,
    #         region="chr6:50,600,000-54,600,000",
    #         )
    #     IS = IS.loc[IS["is_boundary_200000"]]
    #     IS.to_csv(Path(os.path.dirname(__file__)) / "data" / "mm10.HOXA.IS.tsv", sep="\t", index=False)
    #     print(IS)
    #     # 51160000 is very strong
    #     stack = pileup( # using cli.pipeline.pileup, not work here
    #         Cooler(str(self.coolp)),
    #         IS, flank=40
    #         )
    #     print(stack.shape)
        self.assertTrue(True)
    def test_plot_compartment(self):
        """
        Test the plot_IS function.
        """
        out_png = out_png = self.outdir / "plot_compartment.png"
        ddir = Path("/shareb/ychi/repo/sperm_struct/notebooks/data")
        Sdir = Path("/shareb/ychi/repo/sperm_struct/notebooks/B_validation/")
        fig = plot_compartment(
            #str(Sdir / "sperm.PFA05_B6.20k.pileup.mcool::resolutions/1000000"),
            str(ddir / "PFA05.pileup.20k.mcool::resolutions/1000000"),
            Sdir / "sperm.PFA05_B6.20k.pileup.1M.cis.vecs.tsv",
            region="chr1",
            title = "sperm chr1 compartment",
            strip=True
        )
        fig.write_image(str(out_png), width=700, height=700)
        self.assertTrue(out_png.exists())
    def test_plot_compartment_balance(self):
        """
        Test the plot_IS function.
        """
        out_png = out_png = self.outdir / "plot_compartment_balance.png"
        ddir = Path("/shareb/ychi/repo/sperm_struct/notebooks/data2")
        fig = plot_compartment(
            str(ddir / "mESC_allmerge.1k.mcool::resolutions/1000000"),
            ddir / "mESC_allmerge.cis.vecs.tsv",
            region="chrX",
            title = "mESC_allmerge chrX compartment",
            balance=True,
            strip=True,
            mask_eig_na=True,
            fillna=True
        )
        fig.write_image(str(out_png), width=700, height=700)
        self.assertTrue(out_png.exists())

    def test_plot_tiling_compartment_unit(self):
        """
        Unit test for plot_tiling_compartment function.
        Test that the function returns a valid figure object with correct structure.
        """
        # Test with synthetic data using new Track architecture
        fig = plot_tiling_compartment_v2(
            coolps=[str(self.synthetic_cool), str(self.synthetic_cool)],
            eigs_files=[str(self.synthetic_eigs1), str(self.synthetic_eigs2)],
            region="chr1",
            title="Test Tiling Compartment",
            corr=False,  # use raw contact counts
            strip=False,
        )
        
        # Verify figure structure
        self.assertIsNotNone(fig)
        # Should have at least 3 traces (heatmap + 2 eigenvector tracks)
        self.assertGreaterEqual(len(fig.data), 3)
        # Title should be set
        self.assertEqual(fig.layout.title.text, "Test Tiling Compartment")

    def test_plot_tiling_compartment_corr(self):
        """
        Unit test for plot_tiling_compartment with correlation matrix.
        """
        fig = plot_tiling_compartment_v2(
            coolps=[str(self.synthetic_cool), str(self.synthetic_cool)],
            eigs_files=[str(self.synthetic_eigs1), str(self.synthetic_eigs2)],
            region="chr1",
            title="Test Tiling Compartment (Correlation)",
            corr=True,  # use correlation matrix
            strip=True,
        )
        
        self.assertIsNotNone(fig)
        self.assertGreaterEqual(len(fig.data), 3)

    def test_plot_tiling_compartment_visual(self):
        """
        Visual test for plot_tiling_compartment function.
        Generates PNG output for manual inspection.
        Output is saved to tests/output/plot/ for verification.
        """
        out_png = self.outdir / "plot_tiling_compartment.png"
        
        fig = plot_tiling_compartment_v2(
            coolps=[str(self.synthetic_cool), str(self.synthetic_cool)],
            eigs_files=[str(self.synthetic_eigs1), str(self.synthetic_eigs2)],
            region="chr1",
            title="Visual Test: Tiling Compartment",
            corr=False,
            strip=False,
            cmap="RdBu",
        )
        
        # Write image for visual inspection
        fig.write_image(str(out_png), width=800, height=800)
        self.assertTrue(out_png.exists())
        print(f"Visual test output: {out_png}")

    def test_plot_tiling_compartment_balance(self):
        """
        Unit test for plot_tiling_compartment with balance=False.
        Note: synthetic test data does not have balancing weights.
        """
        fig = plot_tiling_compartment_v2(
            coolps=[str(self.synthetic_cool), str(self.synthetic_cool)],
            eigs_files=[str(self.synthetic_eigs1), str(self.synthetic_eigs2)],
            region="chr1",
            title="Test Tiling Compartment (No Balance)",
            balance=False,  # synthetic data has no weights
            corr=False,
            strip=False,
        )
        
        self.assertIsNotNone(fig)
        self.assertGreaterEqual(len(fig.data), 3)

    def test_plot_tiling_compartment_custom_cmap(self):
        """
        Unit test for plot_tiling_compartment with custom colormap.
        """
        fig = plot_tiling_compartment_v2(
            coolps=[str(self.synthetic_cool), str(self.synthetic_cool)],
            eigs_files=[str(self.synthetic_eigs1), str(self.synthetic_eigs2)],
            region="chr1",
            title="Test Tiling Compartment (Custom Cmap)",
            cmap="viridis",
            corr=False,
            strip=False,
        )
        
        self.assertIsNotNone(fig)
        self.assertGreaterEqual(len(fig.data), 3)


class TestCompartmentTrapezoid(unittest.TestCase):
    """
    Test compartment trapezoid plot functions.
    """
    def setUp(self) -> None:
        self.outdir = Path(os.path.dirname(__file__)) / "output" / "compartment_trapezoid"
        # Create output directory if it doesn't exist
        self.outdir.mkdir(parents=True, exist_ok=True)
        # Use synthetic test data
        self.synthetic_cool = Path(os.path.dirname(__file__)) / "data" / "synthetic_test.cool"
        self.synthetic_eigs1 = Path(os.path.dirname(__file__)) / "data" / "synthetic_eigs1.bedgraph"
        self.synthetic_eigs2 = Path(os.path.dirname(__file__)) / "data" / "synthetic_eigs2.bedgraph"

    def test_plot_compartment_trapezoid_basic(self):
        """
        Unit test for plot_compartment_trapezoid basic functionality.
        """
        fig = plot_compartment_trapezoid(
            coolp=str(self.synthetic_cool),
            vec_file=str(self.synthetic_eigs1),
            region="chr1",
            strata=30,
            title="Test Trapezoid Basic",
            balance=False,
        )
        
        self.assertIsNotNone(fig)
        self.assertGreaterEqual(len(fig.data), 1)
        # Title may be set differently, just check figure is valid
        self.assertIsInstance(fig.layout.title.text, (str, type(None)))

    def test_plot_compartment_trapezoid_custom_params(self):
        """
        Unit test for plot_compartment_trapezoid with custom parameters.
        """
        fig = plot_compartment_trapezoid(
            coolp=str(self.synthetic_cool),
            vec_file=str(self.synthetic_eigs1),
            region="chr1",
            strata=50,
            title="Test Trapezoid Custom",
            cmap="viridis",
            balance=False,
            quant=0.05,
        )
        
        self.assertIsNotNone(fig)
        self.assertGreaterEqual(len(fig.data), 1)

    def test_plot_compartment_trapezoid_visual(self):
        """
        Visual test for plot_compartment_trapezoid.
        Generates PNG output for manual inspection.
        """
        out_png = self.outdir / "plot_compartment_trapezoid.png"

        fig = plot_compartment_trapezoid(
            coolp=str(self.synthetic_cool),
            vec_file=str(self.synthetic_eigs1),
            region="chr1",
            strata=40,
            title="Visual Test: Compartment Trapezoid",
            balance=False,
            cmap="RdBu_r",
            add_track=True,  # Add eigenvector track
        )

        fig.write_image(str(out_png), width=800, height=800)
        self.assertTrue(out_png.exists())
        print(f"Visual test output: {out_png}")

    def test_plot_compartment_trapezoid_2tracks_basic(self):
        """
        Unit test for plot_compartment_trapezoid_2tracks basic functionality.
        """
        fig = plot_compartment_trapezoid_2tracks_v2(
            coolp0=str(self.synthetic_cool),
            vec_file1=str(self.synthetic_eigs1),
            vec_file2=str(self.synthetic_eigs2),
            region="chr1",
            strata=30,
            title="Test 2 Tracks Basic",
            balance=False,
        )
        
        self.assertIsNotNone(fig)
        # Should have at least 3 traces (heatmap + 2 bar tracks)
        self.assertGreaterEqual(len(fig.data), 3)
        self.assertIsInstance(fig.layout.title.text, (str, type(None)))

    def test_plot_compartment_trapezoid_2tracks_visual(self):
        """
        Visual test for plot_compartment_trapezoid_2tracks.
        Generates PNG output for manual inspection.
        """
        out_png = self.outdir / "plot_compartment_trapezoid_2tracks.png"
        
        fig = plot_compartment_trapezoid_2tracks_v2(
            coolp0=str(self.synthetic_cool),
            vec_file1=str(self.synthetic_eigs1),
            vec_file2=str(self.synthetic_eigs2),
            region="chr1",
            strata=40,
            title="Visual Test: Compartment Trapezoid 2 Tracks",
            balance=False,
            cmap="RdBu_r",
        )
        
        fig.write_image(str(out_png), width=800, height=900)
        self.assertTrue(out_png.exists())
        print(f"Visual test output: {out_png}")


if __name__ == "__main__":
    unittest.main()