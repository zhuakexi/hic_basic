# unittest
import os
import unittest
import sys
sys.path.insert(0, "/share/home/ychi/dev/hic_basic")
sys.path.insert(0, "/share/home/ychi/dev/hires_utils")
sys.path.insert(0, "/share/home/ychi/dev/sperm_struct")
from hic_basic.DI import normalize_band, multiple_testing_correction, calculate_stats, simpleDiff_test, simpleDiff
import dask.dataframe as dd
import numpy as np
import pandas as pd

class TestDI(unittest.TestCase):
    def setUp(self) -> None:
        self.outdir = "/share/home/ychi/dev/hic_basic/tests/output/DI"
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        self._3dg_path_small = "/shareb/ychi/repo/sperm43/3dg_c/BJ8017.clean.1m.3.3dg"
        self._3dg_path_big = "/shareb/ychi/repo/sperm40_GM/3dg_c/GMO1001.clean.20k.4.3dg"
    # def test_normalize_band_dd(self):
    #     n = 1000
    #     chrom1 = np.random.choice(['chr1', 'chr2', 'chr3'], size=n)
    #     start1 = np.random.randint(0, 100000000, size=n)
    #     chrom2 = np.random.choice(['chr1', 'chr2', 'chr3'], size=n)
    #     start2 = np.random.randint(0, 100000000, size=n)
    #     df = pd.DataFrame({
    #         'chrom1': chrom1,
    #         'start1': start1,
    #         'end1': start1 + 1000000, # 1M
    #         'chrom2': chrom2,
    #         'start2': start2,
    #         'end2': start2 + 1000000, # 1M
    #         'value': np.random.rand(n)
    #     })
    #     result = normalize_band(dd.from_pandas(df, npartitions=1))
    #     self.assertIsInstance(result, pd.DataFrame)
    # def test_normalize_band_df(self):
    #     df = pd.DataFrame({
    #         'chrom1': ['chr1', 'chr1'],
    #         'start1': [100, 200],
    #         'chrom2': ['chr1', 'chr1'],
    #         'start2': [200, 300],
    #         'value': [5, 10]
    #     })
    #     result = normalize_band(df)
    #     self.assertIsInstance(result, pd.DataFrame)
    # def test_normalize_band_imputed(self):
    #     from hic_basic.simpute import cis_distance_graph
    #     _3dg_path = self._3dg_path_small
    #     structure = cis_distance_graph(_3dg_path, fo=None, max_dist=20000000, binsize=1000000)
    #     structure = normalize_band(structure).compute()
    #     print(structure.head(10))
    #     self.assertIsInstance(structure, pd.DataFrame)

    # def test_multiple_testing_correction(self):
    #     pvalues = [0.01, 0.02, 0.03, 0.04, 0.05]
    #     result = multiple_testing_correction(pvalues, correction_type="FDR")
    #     self.assertIsInstance(result, np.ndarray)

    # def test_calculate_stats(self):
    #     row = pd.Series([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    #     m1 = 5
    #     m2 = 5
    #     result = calculate_stats(row, m1, m2)
    #     self.assertIsInstance(result, (float, type(None)))

    # def test_simpleDiff_test(self):
    #     m1, m2, n, sig_n, show_n = 100, 150, 1000, 30, 100
    #     sig_rows = np.random.choice(n, size=sig_n, replace=False)
    #     # random from 0 to 1
    #     df = pd.DataFrame(np.random.rand(n, m1+m2))
    #     #print(df.head(10))
    #     df.loc[sig_rows, :m1] += 0.2
    #     df = dd.from_pandas(df, npartitions=1)
    #     result = simpleDiff_test(df, m1, m2)        #print(result.head(10))
    #     self.assertIsInstance(result, dd.DataFrame)
    #     df = result.compute()
    #     indi = df.sort_values("pvalue").index
    #     for i in range(10, show_n, 10):
    #         print(
    #             "0~{}: ".format(i),
    #             indi[0:i].intersection(sig_rows).shape[0],
    #         )
    #     self.assertLessEqual(
    #         sig_n - indi[0:sig_n].intersection(sig_rows).shape[0],
    #         int(0.1 * sig_n)
    #     ) # wrong rate less than 10%
    # Testing simpleDiff function would require a more complex setup and is not included here
    def test_simpleDiff(self):
        test_cases = [
            ([self._3dg_path_small]*30, [self._3dg_path_small]*20, "chr1", "mm10", os.path.join(self.outdir, "test_simpleDiff.small.parquet"), 0.05, 5, 2000000, 1000000, 16),
            ([self._3dg_path_big]*30, [self._3dg_path_big]*20, "chr1(mat)", "hg19_dip", os.path.join(self.outdir, "test_simpleDiff.big.parquet"), 0.05, 5, 2000000, 20000, 8),
        ]
        for groupA, groupB, chrom, genome, fo, filt_fdr, max_3d_dist, max_dist, binsize, n_jobs in test_cases:
            with self.subTest(groupA=groupA, groupB=groupB, chrom=chrom, genome=genome, fo=fo, filt_fdr=filt_fdr, max_3d_dist=max_3d_dist, max_dist=max_dist, binsize=binsize, n_jobs=n_jobs):
                simpleDiff(
                    groupA=groupA,
                    groupB=groupB,
                    chrom=chrom,
                    genome=genome,
                    fo=fo,
                    filt_fdr=filt_fdr,
                    max_3d_dist=max_3d_dist,
                    max_dist=max_dist,
                    binsize=binsize,
                    n_jobs=n_jobs
                )
                self.assertTrue(os.path.exists(fo))
        
if __name__ == '__main__':
    unittest.main()