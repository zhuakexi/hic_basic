import os
import sys
sys.path.insert(0, "/share/home/ychi/dev/hic_basic")
sys.path.insert(0, "/share/home/ychi/dev/hires_utils")
sys.path.insert(0, "/share/home/ychi/dev/sperm_struct")
import unittest

import pandas as pd
from hic_basic.binnify import GenomeIdeograph
from hic_basic.data import dupref_annote, fetch_cent_chromlen

class TestData(unittest.TestCase):
    """
    Test the data module.
    """
    def test_dupref_annote(self):
        """
        Test the dupref_annote function.
        """
        genomes = GenomeIdeograph("hg19_dip")
        bins = genomes.bins(20e3,bed=True)
        bins = bins.set_index(["chrom","start"])
        bins.drop("end",axis=1,inplace=True)

        ref = pd.read_table(
            "/share/home/ychi/dev/dip-c/color/hg19.cpg.20k.txt",
            names = ["chrom","pos","CpG"],
            index_col = ["chrom","pos"]
            )
        # dupref_annote
        annoted = dupref_annote(bins, ref)
        print(annoted.head())
        self.assertEqual(annoted.shape[1], bins.shape[1]+ref.shape[1])
        self.assertGreater(0.1*annoted.shape[0], annoted.isna().sum().max())
    def test_fetch_cent_chromlen(self):
        """
        Test the fetch_cent_chromlen function.
        """
        # fetch_cent_chromlen
        cent_chromlen = fetch_cent_chromlen("mm10")
        print(cent_chromlen["chrY"])
        self.assertTrue("chrY" in cent_chromlen.keys())
if __name__ == "__main__":
    unittest.main()