# unittest
import os
import unittest
import sys
sys.path.insert(0, "/share/home/ychi/dev/hic_basic")
sys.path.insert(0, "/share/home/ychi/dev/hires_utils")
sys.path.insert(0, "/share/home/ychi/dev/sperm_struct")
from pathlib import Path

import dask.array as da
import dask.dataframe as dd
import numpy as np
import pandas as pd
from scipy.spatial.distance import euclidean
from hic_basic.impute.simpute import boolean_radius_neighbor, cis_proximity_graph, cis_distance_graph, \
    cis_distance_graph_df, parse_3dg_dask, DMimpute
from hires_utils.hires_io import parse_3dg

class TestSimpute(unittest.TestCase):
    def setUp(self) -> None:
        self.outdir = Path(os.path.dirname(__file__)) / "output" / "simpute"
        self._3dg_path_small = "/shareb/ychi/repo/sperm43/3dg_c/BJ8017.clean.1m.3.3dg"
        self._3dg_path_big = "/shareb/ychi/repo/sperm40_GM/3dg_c/GMO1001.clean.20k.4.3dg"
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
    def test_cis_proximity_graph(self):
        # 1M proximity graph
        _3dg_path = self._3dg_path_small
        fo = str(self.outdir / "cis_proximity_graph.cool")
        cis_proximity_graph(_3dg_path, fo, genome="mm10", binsize=1000000)
        self.assertTrue(os.path.exists(fo))
    def test_cis_proximity_graph1(self):
        # 20K proximity graph
        _3dg_path = "/shareb/ychi/repo/sperm_struct/formal_pipeline/3dg_c/HuS02_HuSZ147.clean.20k.4.3dg"
        fo = str(self.outdir / "cis_proximity_graph1.cool")
        cis_proximity_graph(_3dg_path, fo, genome="GRCh38", binsize=20000)
        self.assertTrue(os.path.exists(fo))
    def test_DMimpute(self):
        # 1M distance matrix
        _3dg_path = self._3dg_path_small
        fo = str(self.outdir / "DMimpute.cool")
        fo = DMimpute(_3dg_path, fo, genome="mm10", max_dist=None, binsize=1000000)
        self.assertTrue(os.path.exists(fo))
    def test_parse_3dg_dask(self):
        _3dg_path = self._3dg_path_small
        df = parse_3dg_dask(_3dg_path)
        self.assertTrue(isinstance(df, dd.DataFrame))
        df = df.compute()
        print(df.head())
    def check_distance_graph(self, _3dg_path, pixels):
        # --- check if distance is correct ---
        structure = parse_3dg(_3dg_path)
        if isinstance(pixels, dd.DataFrame):
            picked = pixels.sample(frac=1000/len(pixels))
        else:
            picked = pixels.sample(1000)
        for _, row in picked.iterrows():
            true_value = euclidean(
                structure.loc[(row["chrom1"], row["start1"])],
                structure.loc[(row["chrom2"], row["start2"])]
                )
            # near equal
            self.assertAlmostEqual(row["distance"], true_value, delta=1e-5)
    
    # def test_cis_distance_graph(self):
    #     test_cases = [
    #         (self._3dg_path_small, self.outdir/"cis_distance_graph_small.parquet", 20000000, 1000000),
    #         #(self._3dg_path_big, None, 20000000, 20000),
    #         (self._3dg_path_small, None, 20000000, 1000000),
    #     ]
    #     for _3dg_path, fo, max_dist, binsize in test_cases:
    #         with self.subTest(_3dg_path=_3dg_path, fo=fo, max_dist=max_dist, binsize=binsize):
    #             pixels = cis_distance_graph(_3dg_path, fo, max_dist=max_dist, binsize=binsize)
    #             if fo is not None:
    #                 self.assertTrue(os.path.exists(fo))
    #                 pixels = dd.read_parquet(fo)
    #                 self.check_distance_graph(_3dg_path, pixels)
    #             else:
    #                 print(type(pixels))
    #                 print(pixels.head())
    #                 self.assertTrue(isinstance(pixels, dd.DataFrame))
    #                 self.check_distance_graph(_3dg_path, pixels)
    def test_cis_distance_graph_df(self):
        test_cases = [
            (self._3dg_path_small, None, None, None, 20000000, 1000000, False),
            (self._3dg_path_small, "chr1", None, None, 20000000, 1000000, False),
            (self._3dg_path_big, "chr1(mat)", "hg19_dip", None, 20000000, 20000, False),
            (self._3dg_path_small, "chr1", None, self.outdir / "cis_distance_graph_small_df.parquet", 20000000, 1000000, False),
            (self._3dg_path_big, "chr1(mat)", None, self.outdir / "cis_distance_graph_big_df.parquet", 20000000, 20000, False),
        ]
        for _3dg_path, chrom, genome, fo, max_dist, binsize, fill in test_cases:
            with self.subTest(_3dg_path=_3dg_path, chrom=chrom, genome=genome, fo=fo, 
                              max_dist=max_dist, binsize=binsize, fill=fill):
                pixels = cis_distance_graph_df(_3dg_path, chrom, genome, fo, 
                                               max_dist=max_dist, binsize=binsize, fill=fill)
                if fill is True:
                    self.assertTrue(isinstance(pixels, da))
                else:
                    if fo is not None:
                        self.assertTrue(os.path.exists(fo))
                        pixels = pd.read_parquet(fo)
                        self.check_distance_graph(_3dg_path, pixels)
                    else:
                        self.assertTrue(isinstance(pixels, pd.DataFrame))
                        self.check_distance_graph(_3dg_path, pixels)
                if (genome is not None) and (fill is False):
                    self.assertTrue(pixels["pixel_id"].is_unique)
                else:
                    pass
if __name__ == '__main__':
    unittest.main()