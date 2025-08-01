import os
import unittest
from pathlib import Path

from hic_basic.feature import TSS, CpG
class TestGenomeFeatures(unittest.TestCase):
    """
    Test the genome feature functions.
    """
    def setUp(self):
        self.genome = "GRCh38"
        self.gtf_path = "/share/Data/ychi/genome/GRCh38/gencode.v33.primary_assembly.annotation.gtf"
        self.chrom_sizes_path = "/share/Data/ychi/genome/GRCh38/GRCh38.chr.len"
        self.out_dir = Path(__file__).parent / "output" / "feature"
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
class TestSequenceFeatures(unittest.TestCase):
    """
    Test the genome feature calculate from sequences.
    """
    def setUp(self):
        self.genome = "GRCh38"
        self.fasta = "/share/Data/ychi/genome/GRCh38/raw/GRCh38.primary_assembly.genome.fa"
        self.out_dir = Path(__file__).parent / "output" / "data"
        if not self.out_dir.exists():
            self.out_dir.mkdir(parents=True)
    def test_CpG_compile(self):
        cpg = CpG(
            self.genome,
            db_dir=self.out_dir,
            binsize=1000000,
            flavor="hickit"
        )
        cpg.compile(
            self.fasta,
            force=True
        )
        cpg.load()
        print(cpg.fn)
if __name__ == "__main__":
    unittest.main()