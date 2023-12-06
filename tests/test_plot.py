# unittest
import sys
sys.path.insert(0, "/share/home/ychi/dev/hic_basic")
sys.path.insert(0, "/share/home/ychi/dev/hires_utils")
import os
import unittest
from pathlib import Path

import pandas as pd
from cooler import Cooler
from cooltools.cli.pileup import pileup
from hic_basic.plot.hic import plot_compartment, plot_cool_track, merge_track_data, _plot_mat

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
    def test_pileup_IS(self):
        """
        Test the pileup_IS function.
        """
        out_png = self.outdir / "pileup_IS.png"
        IS = merge_track_data(
            Cooler(str(self.coolp)),
            self.IS,
            region="chr6:50,600,000-54,600,000",
            )
        IS = IS.loc[IS["is_boundary_200000"]]
        IS.to_csv(Path(os.path.dirname(__file__)) / "data" / "mm10.HOXA.IS.tsv", sep="\t", index=False)
        print(IS)
        # 51160000 is very strong
        stack = pileup(
            Cooler(str(self.coolp)),
            IS, flank=40
            )
        print(stack.shape)
        self.assertTrue(True)
    # def test_plot_compartment(self):
    #     """
    #     Test the plot_IS function.
    #     """
    #     out_png = os.path.join(os.path.dirname(__file__), "output", "plot_compartment.png")
    #     ddir = Path("/shareb/ychi/repo/sperm_struct/notebooks/Compartment")
    #     fig = plot_compartment(
    #         str(ddir / "RS.20k.pileup.s14M.1M.cool"),
    #         ddir / "RS.20k.pileup.s14M.1M.cis.vecs.tsv",
    #         region="chr1",
    #         title = "Compartment",
    #     )
    #     fig.savefig(out_png, dpi=300)
if __name__ == "__main__":
    unittest.main()