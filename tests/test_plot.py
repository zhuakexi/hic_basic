# unittest
import sys
sys.path.insert(0, "/share/home/ychi/dev/hic_basic")
sys.path.insert(0, "/share/home/ychi/dev/hires_utils")
import os
import unittest
from pathlib import Path

import pandas as pd
from cooler import Cooler
from hic_basic.plot.hic import plot_compartment, plot_cool_track

class TestPlot(unittest.TestCase):
    """
    Test the plot module.
    """
    def setUp(self) -> None:
        self.outdir = Path(os.path.dirname(__file__)) / "output" / "plot"
        self.coolp = "/share/Data/ychi/raw/Bonev2017/ES.mcool::resolutions/20000"
        self.IS = "/shareb/ychi/repo/sperm_struct/notebooks/data/mESC.20k.IS.tsv"
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
    def test_plot_compartment(self):
        """
        Test the plot_IS function.
        """
        out_png = os.path.join(os.path.dirname(__file__), "output", "plot_compartment.png")
        ddir = Path("/shareb/ychi/repo/sperm_struct/notebooks/Compartment")
        fig = plot_compartment(
            str(ddir / "RS.20k.pileup.s14M.1M.cool"),
            ddir / "RS.20k.pileup.s14M.1M.cis.vecs.tsv",
            region="chr1",
            title = "Compartment",
        )
        fig.savefig(out_png, dpi=300)
if __name__ == "__main__":
    unittest.main()