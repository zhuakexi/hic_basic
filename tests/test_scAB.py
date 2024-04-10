import os
import sys
import unittest
from pathlib import Path

sys.path.insert(0, "/share/home/ychi/dev/hic_basic")
sys.path.insert(0, "/share/home/ychi/dev/hires_utils")
sys.path.insert(0, "/share/home/ychi/dev/sperm_struct")

import pandas as pd
from hic_basic.scAB_embedding import s_color2
from hires_utils.hires_io import parse_3dg

class TestScAB(unittest.TestCase):
    def setUp(self) -> None:
        self._3dg_file = "/sharec/ychi/repo/sperm55/3dg_c/BJ9205.clean.200k.4.3dg"
        self.color_file = "/shareb/ychi/repo/sperm_struct/notebooks/data2/mm10.CpG.200k.tsv.gz"
        self._3dg_file_big = "/sharec/ychi/repo/sperm55/3dg_c/BJ9205.clean.20k.4.3dg"
        self.color_file_big = "/shareb/ychi/repo/sperm_struct/notebooks/data2/mm10.CpG.20k.tsv.gz"
        self.outdir = os.path.join(os.path.dirname(__file__), "output", "scAB")
        if not Path(self.outdir).exists():
            os.makedirs(self.outdir)
    def test_s_color2(self):
        _3dg = parse_3dg(self._3dg_file)
        color_data = pd.read_table(
            self.color_file,
            names = ["chrom","start","end","CpG"]
            )[['chrom','start','CpG']]
        res = s_color2(
            _3dg,
            color_data
            )
        print(res)
        self.assertTrue(
            isinstance(res, pd.DataFrame)
        )
    def test_s_color2_big(self):
        _3dg = parse_3dg(self._3dg_file_big)
        color_data = pd.read_table(
            self.color_file_big,
            names = ["chrom","start","end","CpG"]
            )[['chrom','start','CpG']]
        res = s_color2(
            _3dg,
            color_data
            )
        print(res)
        self.assertTrue(
            isinstance(res, pd.DataFrame)
        )
if __name__ == "__main__":
    unittest.main()