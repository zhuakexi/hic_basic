import os
import sys
sys.path.insert(0, "/share/home/ychi/dev/hic_basic")
sys.path.insert(0, "/share/home/ychi/dev/hires_utils")
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from hic_basic.impute.schicluster import solve_rwr_inverse, schicluster_imputation_for_mat, schicluster_impute

class TestSchicluster(unittest.TestCase):
    def setUp(self) -> None:
        self.outdir = Path(os.path.dirname(__file__)) / "output" / "schicluster"
        if not self.outdir.exists():
            self.outdir.mkdir(parents=True)
        self.coolp = "/share/Data/ychi/raw/Bonev2017/ES.mcool::resolutions/20000"
        self.coolp1M = "/share/Data/ychi/raw/Bonev2017/ES.mcool::resolutions/1000000"
    def test_solve_rwr_inverse(self):
        # Define a sample stochastic matrix
        stoch_matrix = np.array([[0.1, 0.2, 0.7], [0.3, 0.4, 0.3], [0.2, 0.5, 0.3]])

        # Call the function with the sample input
        result = solve_rwr_inverse(stoch_matrix)
        print(result)
        self.assertIsInstance(result, np.ndarray)

    def test_schicluster_imputation_for_mat(self):
        # Define a sample matrix
        mat = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

        # Call the function with the sample input
        result = schicluster_imputation_for_mat(mat)
        print(result)
        self.assertIsInstance(result, np.ndarray)

    def test_schicluster_impute(self):
        # Define a sample cooler file path
        coolp = self.coolp1M

        # Define a sample output prefix
        outprefix = self.outdir / "schicluster_impute"

        # Call the function with the sample inputs
        chroms = ["chr1","chr2","chr3"]
        schicluster_impute(coolp, chroms, outprefix)

        # Assert that the output file exists
        self.assertTrue(os.path.exists(f"{outprefix}.chr1.parquet"))

        result = pd.read_parquet(f"{outprefix}.chr1.parquet")
        print(result)
        self.assertIsInstance(result, pd.DataFrame)

if __name__ == '__main__':
    unittest.main()