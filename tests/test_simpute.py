# unittest
import os
import unittest
import sys
sys.path.insert(0, "/share/home/ychi/dev/hic_basic")
sys.path.insert(0, "/share/home/ychi/dev/hires_utils")
sys.path.insert(0, "/share/home/ychi/dev/sperm_struct")
from pathlib import Path

import dask.dataframe as dd
import numpy as np
import pandas as pd
from scipy.spatial.distance import euclidean
from hic_basic.simpute import boolean_radius_neighbor, cis_proximity_graph, calculate_distance, cis_distance_graph
from hires_utils.hires_io import parse_3dg

class TestSimpute(unittest.TestCase):
    def setUp(self) -> None:
        self.outdir = Path(os.path.dirname(__file__)) / "output" / "simpute"
        self._3dg_path_small = "/shareb/ychi/repo/sperm43/3dg_c/BJ8017.clean.1m.3.3dg"
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
    def test_cis_proximity_graph(self):
        _3dg_path = self._3dg_path_small
        fo = str(self.outdir / "cis_proximity_graph.cool")
        cis_proximity_graph(_3dg_path, fo, genome="mm10", binsize=1000000)
        self.assertTrue(os.path.exists(fo))

    def test_calculate_distance(self):
        row1 = pd.Series({'x': 1, 'y': 2, 'z': 3})
        row2 = pd.Series({'x': 4, 'y': 5, 'z': 6})
        result = calculate_distance(row1, row2)
        self.assertEqual(result, np.sqrt((3**2) + (3**2) + (3**2)))

    def test_cis_distance_graph(self):
        _3dg_path = self._3dg_path_small
        fo = self.outdir / "cis_distance_graph.parquet"
        cis_distance_graph(_3dg_path, fo, max_dist=20000000, binsize=1000000)
        self.assertTrue(os.path.exists(fo))
        pixels = pd.read_parquet(fo)
        #print(pixels.head(10))
        print(pixels.shape)
        self.assertTrue(isinstance(pixels, dd.DataFrame))
        # --- check if distance is correct ---
        structure = parse_3dg(_3dg_path)
        picked = pixels.sample(1000)
        for _, row in picked.iterrows():
            true_value = euclidean(
                structure.loc[(row["chrom1"], row["start1"])],
                structure.loc[(row["chrom2"], row["start2"])]
                )
            # near equal
            self.assertAlmostEqual(row["distance"], true_value, delta=1e-5)
    def test_cis_distance_graph_noout(self):
        _3dg_path = self._3dg_path_small
        fo = None
        pixels = cis_distance_graph(_3dg_path, fo, max_dist=20000000, binsize=1000000)
        print(pixels.shape)
        self.assertTrue(isinstance(pixels, dd.DataFrame))
        # --- check if distance is correct ---
        structure = parse_3dg(_3dg_path)
        picked = pixels.sample(frac=1000/structure.shape[0]**2)
        for _, row in picked.iterrows():
            true_value = euclidean(
                structure.loc[(row["chrom1"], row["start1"])],
                structure.loc[(row["chrom2"], row["start2"])]
                )
            # near equal
            self.assertAlmostEqual(row["distance"], true_value, delta=1e-5)

if __name__ == '__main__':
    unittest.main()