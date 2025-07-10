import os
import unittest
from pathlib import Path

import pandas as pd
from hic_basic.binnify import GenomeIdeograph
from hic_basic.data import dupref_annote, fetch_cent_chromlen, TSS

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
        print(cent_chromlen.loc["chrY"])
        self.assertTrue("chrY" in cent_chromlen.index)
        cent_chromlen = fetch_cent_chromlen("GRCh38")
        print(cent_chromlen.loc["chrY"])
        self.assertTrue("chrY" in cent_chromlen.index)
        # test dip
        cent_chromlen = fetch_cent_chromlen("mm10_dip")
        print(cent_chromlen)
        self.assertTrue("chr1(mat)" in cent_chromlen.index)
class TestGenomeFeatures(unittest.TestCase):
    """
    Test the genome feature functions.
    """
    def setUp(self):
        self.genome = "GRCh38"
        self.gtf_path = "/share/Data/ychi/genome/GRCh38/gencode.v33.primary_assembly.annotation.gtf"
        self.chrom_sizes_path = "/share/Data/ychi/genome/GRCh38/GRCh38.chr.len"
        self.out_dir = Path(__file__).parent / "output" / "data"
        if not self.out_dir.exists():
            self.out_dir.mkdir(parents=True)
    def test_get_tss_region_from_gtf(self):
        """
        Test the get_tss_region_from_gtf function.
        """
        tss_regions = TSS.get_tss_region_from_gtf(
            self.gtf_path,
            cache_file = self.out_dir / "tss_regions.csv",
        )
        print(tss_regions.head())
    def test_TSS(self):
        """
        Test the TSS class.
        """
        tss = TSS(
            self.genome
        )
        tss.compile(
            self.gtf_path,
            force = True
        )
        data = tss.data
        print(data.head())
    def test_TSS_cache(self):
        """
        Test the TSS class with cache.
        """
        tss = TSS(
            self.genome
        )
        tss.compile(
            self.gtf_path,
            force = False
        )
        data = tss.data
        print(data.head())
if __name__ == "__main__":
    unittest.main()