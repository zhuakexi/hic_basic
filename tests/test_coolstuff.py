import os
import sys
sys.path.insert(0, "/share/home/ychi/dev/hic_basic")
sys.path.insert(0, "/share/home/ychi/dev/hires_utils")
sys.path.insert(0, "/share/home/ychi/dev/sperm_struct")
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from hic_basic.coolstuff import cool2mat, cli_pileup, cli_expected, reset_cool_bins
from hic_basic.plot.hic import _plot_mat

class TestCoolstuff(unittest.TestCase):
    def setUp(self) -> None:
        self.outdir = Path(os.path.dirname(__file__)) / "output" / "coolstuff"
        if not self.outdir.exists():
            self.outdir.mkdir(parents=True)
        self.coolp = "/share/Data/ychi/raw/Bonev2017/ES.mcool::resolutions/20000"
        self.coolp1M = "/share/Data/ychi/raw/Bonev2017/ES.mcool::resolutions/1000000"
        self.IS = Path(os.path.dirname(__file__)) / "data" / "mm10.HOXA.IS.tsv"
    def test_cool2mat(self):
        sample_cool = self.coolp1M
        regions = [
            "chr1",
            "chr1:1000000-2000000",
            ["chr1:1,000,000-2,000,000", "chr2"],
            slice(0,-1),
            [slice(0,1000000), slice(1000000,2000000)]
        ]

        for region in regions:
            # Call the function with the sample inputs
            result = cool2mat(sample_cool, region)
            # Check the result
            print(result.shape)
            print(result.head())
            self.assertIsInstance(result, pd.DataFrame)

    def test_cli_pileup(self):
        """
        Test the cli_pileup function.
        """
        coolp = self.coolp
        features = self.IS # with 9 rows
        outfile = self.outdir / "cli_pileup.npz"
        cwd = self.outdir
        cli_pileup(
            coolp,
            features,
            outfile,
            cwd = cwd,
            flank=400000
        )
        self.assertTrue(outfile.exists())

        # check the shape of pileup
        outpng = self.outdir / "cli_pileup_bd1.png"
        loaded = np.load(outfile)
        pileup, stack = loaded["pileup"], loaded["stack"]
        self.assertEqual(
            stack.shape[2],
            pd.read_table(features).shape[0]
        )
        print(pileup.shape)

        # print first TAD border
        TADborder1 = pd.DataFrame(stack[:, :, 0])
        #print(TADborder1)
        fig = _plot_mat(
            TADborder1,
            donorm=True,
            balancing=True
        )
        fig.write_image(str(outpng), width=800, height=800)
        self.assertTrue(outpng.exists())

        # print full pileup
        outpng = self.outdir / "cli_pileup.png"
        full_pileup = pd.DataFrame(pileup)
        fig = _plot_mat(
            full_pileup,
            donorm=True,
            balancing=True
        )
        fig.write_image(str(outpng), width=800, height=800)
        self.assertTrue(outpng.exists())
    def cli_expected(self):
        """
        Test the cli_expected function.
        """
        output = self.outdir / "cli_expected.tsv"
        cli_expected(
            self.coolp1M,
            output,
            cwd=self.outdir
        )
        self.assertTrue(output.exists())
    def test_reset_cool_bins(self):
        """
        Restrict the cool bins to the standard chromosomes.
        """
        coolp = self.coolp1M
        output = self.outdir / "reset_cool_bins.1M.cool"
        reset_cool_bins(
            str(coolp),
            str(output),
            genome="mm10", chunksize=1e6
            )
        # check existence
        self.assertTrue(output.exists())
if __name__ == "__main__":
    unittest.main()

