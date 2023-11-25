import sys
sys.path.insert(0, "/share/home/ychi/dev/hic_basic")
sys.path.insert(0, "/share/home/ychi/dev/hires_utils")
sys.path.insert(0, "/share/home/ychi/dev/sperm_struct")
import os
import unittest
from io import StringIO
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
from hic_basic.plot.render import surface_pymol
from hic_basic.structure.measure import primary_views
from hires_utils.hires_io import parse_3dg
def write_lr_view(_3dg_primary_views, outpng):
    pfigs = dict(zip(
        _3dg_primary_views["name_of_vectors"][:-1], # exclude the last(center) vector
        [_3dg_primary_views["primary_figures"][i][0] for i in range(3)] # get the first (indi:0) figure of each view
    ))
    mat = pfigs["left-right"]
    fig = go.Figure()
    fig.add_trace(
            go.Heatmap(
                z = mat,
                showlegend=False,
            )
    )
    fig.update_layout(
        height = 600,
        width = 500,
        showlegend = False,
    )
    fig.write_image(outpng)
class TestMeasure(unittest.TestCase):
    def test_primary_views(self):
        # a structure with wandering Y fragments
        _3dg_file = "/shareb/ychi/repo/sperm27/3dg_c/PD10_XD_B6J016084.clean.20k.1.3dg"
        outpng = os.path.join(os.path.dirname(__file__), "output", "withY_primary_figure.png")
        _3dg_primary_views = primary_views(
            parse_3dg(_3dg_file),
            ngrid = 64,
            keep_main=False
            )
        write_lr_view(_3dg_primary_views, outpng)
        self.assertTrue(Path(outpng).exists())
        # try to remove wandering Y fragments
        outpng = os.path.join(os.path.dirname(__file__), "output", "withoutY_primary_figure.png")
        _3dg_primary_views = primary_views(
            parse_3dg(_3dg_file),
            ngrid = 64,
            keep_main=True
            )
        write_lr_view(_3dg_primary_views, outpng)
        self.assertTrue(Path(outpng).exists())
    def test_primary_views_alot(self):
        ddir = "/share/home/ychi/dev/sperm_struct/ds_pipeline/smk/config/"
        sample_table = pd.read_csv(
            ddir + "PFA05_strain_schicluster_sample_table.csv.gz",
            index_col=0
            ).query('cell_state == "hapmal"')
        gs = sample_table["20k_g_struct1"].dropna()[:2]
        for _3dg_file in gs:
            _3dg_primary_views = primary_views(
                parse_3dg(_3dg_file),
                ngrid = 64,
                keep_main=True
                )
            outpng = os.path.join(os.path.dirname(__file__), "output", os.path.basename(_3dg_file) + "_primary_figure.png")
            write_lr_view(_3dg_primary_views, outpng)
            self.assertTrue(Path(outpng).exists())
if __name__ == "__main__":
    unittest.main()