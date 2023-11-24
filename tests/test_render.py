# unittest
import sys
sys.path.insert(0, "/share/home/ychi/dev/hic_basic")
sys.path.insert(0, "/share/home/ychi/dev/hires_utils")
sys.path.insert(0, "/share/home/ychi/dev/sperm_struct")
import os
import unittest
from io import StringIO
from pathlib import Path

import plotly.graph_objects as go
from hic_basic.plot.render import surface_pymol
from hires_utils.hires_io import parse_3dg
from hic_basic.structure.measure import primary_views
from lib.struct import sig_primary_coords
class TestRender(unittest.TestCase):
    """
    Test the render module.
    """
    def test_surface_pymol(self):
        """
        Test the surface_pymol function.
        """
        _3dg_file = "/shareb/ychi/repo/sperm43/3dg_c/BJ8017.clean.20k.3.3dg"
        tmpdir = os.path.join(os.path.dirname(__file__), "output")

        # --- compare with primary_views ---
        outpng = os.path.join(os.path.dirname(__file__), "output", "surface_primary_figure.png")
        _3dg_primary_views = primary_views(
            parse_3dg(_3dg_file),
            ngrid = 64
            )
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
        self.assertTrue(Path(outpng).exists())

        # --- direct call ---
        outpng = os.path.join(os.path.dirname(__file__), "output", "surface_pymol.png")
        surface_pymol(
            _3dg_file,
            outpng,
            tmpdir
        )
        self.assertTrue(Path(outpng).exists())

        # --- rotate call ---
        outpng = os.path.join(os.path.dirname(__file__), "output", "surface_pymol_rotate.png")
        _3dg =  sig_primary_coords(
                primary_views(parse_3dg(_3dg_file)),
                [1,1,1],
                parse_3dg(_3dg_file)
                )
        surface_pymol(
            StringIO(_3dg.to_csv(sep="\t", index=True, header=False)),
            outpng,
            tmpdir
        )
        self.assertTrue(Path(outpng).exists())
        
if __name__ == "__main__":
    unittest.main()