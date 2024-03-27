import os
import sys
sys.path.insert(0, "/share/home/ychi/dev/hic_basic")
sys.path.insert(0, "/share/home/ychi/dev/hires_utils")
import unittest
from pathlib import Path

import pandas as pd
from hic_basic.binnify import GenomeIdeograph
from hic_basic.sequence import count_CpG
class TestSequence(unittest.TestCase):
    def setUp(self):
        self.outdir = Path(os.path.dirname(__file__)) / "output" / "sequence"
        self.binsize = 20000
        self.bed_df = GenomeIdeograph("mm10").bins(
            #1000000,
            self.binsize,
            bed=True,
            order=True
        )
        self.fasta = "/share/home/ychi/data/genome/GRCm38/raw/mm10.fa"
    def test_count_CpG(self):
        result = count_CpG(self.bed_df, self.fasta)
        print(result)
        self.assertIsInstance(result, pd.DataFrame)
    def test_count_CpG_compare_color(self):
        # compare with Tan2018's result
        bed_df = self.bed_df.copy()
        bed_df = bed_df[bed_df['chrom'] == 'chr1'].head(200)
        bed_df["start"] = bed_df["start"] + self.binsize//2
        bed_df["end"] = bed_df["end"] + self.binsize//2
        result = count_CpG(bed_df, self.fasta)
        print(result)
        self.assertIsInstance(result, pd.DataFrame)
        # exactly the same
        target_df = pd.DataFrame(
            [
                ["chr1",3000000,0.00321865],
                ["chr1",3020000,0.006],
                ["chr1",3040000,0.00675],
                ["chr1",3060000,0.00825]
        ])
if __name__ == "__main__":
    unittest.main()