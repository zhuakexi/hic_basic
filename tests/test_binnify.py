# unittest
import os
import unittest
import sys
sys.path.insert(0, "/share/home/ychi/dev/hic_basic")
sys.path.insert(0, "/share/home/ychi/dev/hires_utils")
sys.path.insert(0, "/share/home/ychi/dev/sperm_struct")
from pathlib import Path

import numpy as np
import pandas as pd
from hic_basic.binnify import GenomeIdeograph

class TestGenomeIdeograph(unittest.TestCase):
    def setUp(self):
        self.genome = GenomeIdeograph('mm10')
    def test_chromosomes(self):
        print(self.genome.chromosomes.index.dtype)
    def test_breaks(self):
        binsize = 1000
        breaks = self.genome.breaks(binsize)
        self.assertIsInstance(breaks, dict)
        for chrom in breaks:
            self.assertTrue(all(i <= j for i, j in zip(breaks[chrom], breaks[chrom][1:])))  # Assert that each list is sorted
    def test_bins(self):
        binsize = 1000
        bins = self.genome.bins(binsize)
        self.assertIsInstance(bins, dict)
        for chrom in bins:
            self.assertIsInstance(bins[chrom], pd.IntervalIndex)

        bins_bed = self.genome.bins(binsize, bed=True)
        self.assertIsInstance(bins_bed, pd.DataFrame)
        self.assertTrue('chrom' in bins_bed.columns)
        self.assertTrue('start' in bins_bed.columns)
        self.assertTrue('end' in bins_bed.columns)
    def test_bins_order(self):
        binsize = 1000
        bins = self.genome.bins(binsize, bed=True, order=True)
        print(bins)
        # chr1, chr2 ... not chr1, chr10, chr11, chr12, ...
        self.assertTrue(bins.query('chrom == "chr2"').index[0] < bins.query('chrom == "chr11"').index[0])
    def test_pixel_id(self):
        binsize = 1000
        row = pd.Series({
            'chrom1': 'chr1',
            'start1': 1000,
            'end1': 2000,
            'chrom2': 'chr1',
            'start2': 3000,
            'end2': 4000
        })
        pixel_id = self.genome.pixel_id(row, binsize, intra=True)
        print(type(pixel_id))
        self.assertIsInstance(pixel_id, int)
        self.assertEqual(pixel_id, 195474)
    def test_join_pixel_id(self):
        binsize = 1000000
        bedpe = Path(os.path.dirname(__file__)) / "output" / "simpute" / "cis_distance_graph_small.parquet"
        bedpe = pd.read_parquet(bedpe)
        bedpe = self.genome.join_pixel_id(bedpe, binsize)
        print(bedpe.head(10))
        self.assertTrue('pixel_id' in bedpe.columns)
        # must give uniq id
        self.assertTrue(bedpe["pixel_id"].is_unique)
    def test_join_positions(self):
        binsize = 200000
        df = pd.DataFrame(
            {
                "pixel_id" : np.random.randint(0, 1000000, 1000),
                "meanA" : np.random.rand(1000),
                "meanB" : np.random.rand(1000),
                "diff" : np.random.rand(1000),
                "pval" : np.random.rand(1000)
            }
        )
        orig_pixel_id = df["pixel_id"].copy()
        res = self.genome.join_positions(df, binsize).drop(columns=["pixel_id"])
        taret_pixel_id = self.genome.join_pixel_id(res, binsize)["pixel_id"]
        self.assertTrue(all(orig_pixel_id == taret_pixel_id))
if __name__ == '__main__':
    unittest.main()