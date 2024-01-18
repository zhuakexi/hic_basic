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
from hic_basic.simpute import boolean_radius_neighbor, cis_proximity_graph, calculate_distance, cis_distance_graph

class TestSimpute(unittest.TestCase):
    def setUp(self) -> None:
        self.outdir = Path(os.path.dirname(__file__)) / "output" / "simpute"
        self._3dg_path = "/shareb/ychi/repo/sperm43/3dg_c/BJ8017.clean.1m.3.3dg"
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
    def test_cis_proximity_graph(self):
        _3dg_path = self._3dg_path
        fo = str(self.outdir / "cis_proximity_graph.cool")
        cis_proximity_graph(_3dg_path, fo, genome="mm10", binsize=1000000)
        self.assertTrue(os.path.exists(fo))

    def test_calculate_distance(self):
        row1 = pd.Series({'x': 1, 'y': 2, 'z': 3})
        row2 = pd.Series({'x': 4, 'y': 5, 'z': 6})
        result = calculate_distance(row1, row2)
        self.assertEqual(result, np.sqrt((3**2) + (3**2) + (3**2)))

    def test_cis_distance_graph(self):
        _3dg_path = self._3dg_path
        fo = self.outdir / "cis_distance_graph.parquet"
        cis_distance_graph(_3dg_path, fo, binsize=1000000)
        self.assertTrue(os.path.exists(fo))

if __name__ == '__main__':
    unittest.main()