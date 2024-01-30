# unittest
import unittest
import shutil
import sys
sys.path.insert(0, "/share/home/ychi/dev/hic_basic")
sys.path.insert(0, "/share/home/ychi/dev/hires_utils")
sys.path.insert(0, "/share/home/ychi/dev/sperm_struct")
from pathlib import Path

import dask.dataframe as dd
import numpy as np
import pandas as pd
from dask.distributed import Client, LocalCluster
from hic_basic.DI import process_file, GenomeIdeograph 

class TestProcessFile(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # 设置 Dask 客户端用于测试
        cls.cluster = LocalCluster(processes=False, threads_per_worker=2, n_workers=1, memory_limit='2GB')
        cls.client = Client(cls.cluster)
        cls.test_data_dir = Path(__file__).parent / "output/DI"

    @classmethod
    def tearDownClass(cls):
        cls.client.close()
        cls.cluster.close()

    def setUp(self):
        # 确保测试数据目录存在
        self.test_data_dir.mkdir(parents=True, exist_ok=True)

    def tearDown(self):
        # 清理测试数据
        if self.test_data_dir.exists():
            shutil.rmtree(self.test_data_dir)

    def test_process_file_with_mock_data(self):
        # 创建 Mock 数据
        mock_data = {
            'chrom1': ['chr1', 'chr1'],
            'start1': [0, 100000],
            'end1': [100000, 200000],
            'chrom2': ['chr1', 'chr1'],
            'start2': [50000, 150000],
            'distance': [0.1, 0.2]
        }
        mock_df = pd.DataFrame(mock_data)
        mock_file = self.test_data_dir / "mock_data.parquet"
        mock_df.to_parquet(mock_file)

        # 调用 process_file 函数
        result = process_file(
            str(mock_file),
            'chr1',
            'mm10',
            500000,
            100000
        )

        # 检查返回的 Dask DataFrame
        self.assertTrue(isinstance(result, dd.Series))

    def test_process_file_with_real_data(self):
        # 使用实际数据文件
        real_file = Path(__file__).parent / "output" / "simpute" / "cis_distance_graph_small.parquet"

        # 确保文件存在
        self.assertTrue(real_file.exists())

        # 调用 process_file 函数
        result = process_file(
            str(real_file),
            'chr1',
            'mm10',
            2000000,
            1000000
        )

        # 检查返回的 Dask DataFrame
        #print(result.compute())
        print(result[result.notnull()].compute())
        print("Total values:", result.compute().shape[0])
        print("Non-null values:", result[result.notnull()].compute().shape[0])
        self.assertTrue(isinstance(result, dd.Series))
        # load all data
        dfs = []
        N=3
        for imputed_path in [real_file]*N:
            dfs.append(process_file(imputed_path, "chr1", "mm10", 2000000, 1000000))
        combined_df = dd.concat(dfs, axis=1)
        combined_df.columns = [f"sample{i}" for i in range(N)]
        print("combined_df:", combined_df.compute())
        print("combined_df not null:", combined_df[combined_df.notnull().any(axis=1)].compute())

if __name__ == '__main__':
    unittest.main()

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
    # def test_simpleDiff(self):
    #     test_cases = [
    #         #([self._3dg_path_small]*30, [self._3dg_path_small]*20, "chr1", "mm10", os.path.join(self.outdir, "test_simpleDiff.small.parquet"), 0.05, 5, 2000000, 1000000, 16),
    #         #([self._3dg_path_middle]*30, [self._3dg_path_middle]*20, "chr1(mat)", "hg19_dip", os.path.join(self.outdir, "test_simpleDiff.middle.parquet"), 0.05, 5, 2000000, 200000, 8),
    #         #([self._3dg_path_big]*30, [self._3dg_path_big]*20, "chr1(mat)", "hg19_dip", os.path.join(self.outdir, "test_simpleDiff.big.parquet"), 0.05, 5, 2000000, 20000, 8),
    #         ([self._imputed_small]*30, [self._imputed_small]*20, "chr1", "mm10", os.path.join(self.outdir, "test_simpleDiff.imputed.parquet"), 0.05, 5, 2000000, 1000000, 8),
    #     ]
    #     for groupA, groupB, chrom, genome, fo, filt_fdr, max_3d_dist, max_dist, binsize, n_jobs in test_cases:
    #         with self.subTest(groupA=groupA, groupB=groupB, chrom=chrom, genome=genome, fo=fo, filt_fdr=filt_fdr, max_3d_dist=max_3d_dist, max_dist=max_dist, binsize=binsize, n_jobs=n_jobs):
    #             simpleDiff(
    #                 groupA=groupA,
    #                 groupB=groupB,
    #                 chrom=chrom,
    #                 genome=genome,
    #                 fo=fo,
    #                 filt_fdr=filt_fdr,
    #                 max_3d_dist=max_3d_dist,
    #                 max_dist=max_dist,
    #                 binsize=binsize,
    #                 n_jobs=n_jobs
    #             )
    #             self.assertTrue(os.path.exists(fo))
        
# if __name__ == '__main__':
#     unittest.main()