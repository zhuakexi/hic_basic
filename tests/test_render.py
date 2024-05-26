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
from hic_basic.plot.render import surface_territory_pymol, surface_b_pymol, clip_b_pymol, \
    surface_centelo_pymol, clip_centelo_pymol, highlight_surface_b_pymol, clip_single_centelo_pymol, \
    centelo_relpos
from hires_utils.hires_io import parse_3dg
from hic_basic.structure.measure import primary_views
from lib.struct import sig_primary_coords
def please_pymol(_3dg_file):
    _3dg = parse_3dg(_3dg_file).reset_index()
    _3dg["chr"] = _3dg["chr"].str.replace("[()]", "_", regex=True)
    _3dg = _3dg.set_index(["chr","pos"], drop=True)
    return StringIO(_3dg.to_csv(sep="\t", index=True, header=False))
class TestRender(unittest.TestCase):
    """
    Test the render module.
    """
    # --- test helper functions ---
    def test_centelo_relpos(self):
        _3dg = parse_3dg("/shareb/ychi/repo/sperm43/3dg_c/BJ8017.clean.1m.3.3dg")
        b_factor = centelo_relpos(_3dg.index, "mm10")
        self.assertTrue(b_factor.shape[0] == _3dg.shape[0])
    # --- test basic functions ---
    def test_surface_territory_pymol(self):
        """
        Test the surface_pymol function.
        """
        print("Test_surface_territory_pymol")
        outpng = os.path.join(os.path.dirname(__file__), "output", "surface_territory_pymol.png")
        surface_territory_pymol(
            please_pymol("/shareb/ychi/repo/sperm40_GM/3dg_c/GMO1001.clean.1m.4.3dg"),
            outpng,
            tmpdir=os.path.join(os.path.dirname(__file__), "output"),
            conda=None
            )
        self.assertTrue(Path(outpng).exists())
    def test_surface_b_pymol(self):
        print("Test_surface_b_pymol")
        outpng = os.path.join(os.path.dirname(__file__), "output", "surface_cpg_pymol.png")
        # cpg
        surface_b_pymol(
            "/shareb/ychi/repo/sperm40_GM/3dg_c/GMO1001.clean.1m.4.3dg",
            "/share/home/ychi/software/dip-c/color/hg19.cpg.1m.txt",
            outpng,
            tmpdir=os.path.join(os.path.dirname(__file__), "output"),
            conda=None
            )
        self.assertTrue(Path(outpng).exists())
    def test_clip_b_pymol(self):
        print("Test_clip_b_pymol")
        outpng = os.path.join(os.path.dirname(__file__), "output", "clip_cpg_pymol.png")
        # cpg
        clip_b_pymol(
            "/shareb/ychi/repo/sperm40_GM/3dg_c/GMO1001.clean.1m.4.3dg",
            "/share/home/ychi/software/dip-c/color/hg19.cpg.1m.txt",
            outpng,
            tmpdir=os.path.join(os.path.dirname(__file__), "output"),
            conda=None
            )
        self.assertTrue(Path(outpng).exists())
    def test_clip_b_pymol_target(self):
        print("Test_clip_b_pymol_target")
        outpng = os.path.join(os.path.dirname(__file__), "output", "clip_cpg_target_pymol_target.png")
        # cpg
        clip_b_pymol(
            please_pymol("/shareb/ychi/repo/sperm40_GM/3dg_c/GMO1001.clean.1m.4.3dg"),
            "/share/home/ychi/software/dip-c/color/hg19.cpg.1m.txt",
            outpng,
            #target=["chrX_mat_", "chrY_mat_"],
            cmap = "magenta green, chain chrX_mat_ or chain chrY_mat_, 0.005, 0.02",
            tmpdir=os.path.join(os.path.dirname(__file__), "output"),
            conda=None
            )
        self.assertTrue(Path(outpng).exists())
    def test_highlight_b_pymol(self):
        print("Test_highlight_b_pymol")
        outpng = os.path.join(os.path.dirname(__file__), "output", "highlight_cpg_centelo_pymol.png")
        highlight_surface_b_pymol(
            please_pymol("/shareb/ychi/repo/sperm40_GM/3dg_c/GMO1001.clean.1m.4.3dg"),
            "/share/home/ychi/software/dip-c/color/hg19.cpg.1m.txt",
            "chrX_mat_",
            outpng,
            tmpdir=os.path.join(os.path.dirname(__file__), "output"),
            conda=None
            )
        self.assertTrue(Path(outpng).exists())
    # --- test useful functions ---
    def test_surface_centelo_pymol(self):
        print("Test_surface_centelo_pymol")
        outpng = os.path.join(os.path.dirname(__file__), "output", "surface_centelo_pymol.png")
        # centelo
        surface_centelo_pymol(
            "/shareb/ychi/repo/sperm43/3dg_c/BJ8017.clean.1m.3.3dg",
            outpng,
            "mm10",
            tmpdir=os.path.join(os.path.dirname(__file__), "output"),
            conda=None
            )
        self.assertTrue(Path(outpng).exists())
    def test_clip_centelo_pymol(self):
        print("Test_clip_centelo_pymol")
        outpng = os.path.join(os.path.dirname(__file__), "output", "clip_centelo_pymol.png")
        # centelo
        clip_centelo_pymol(
            "/shareb/ychi/repo/sperm43/3dg_c/BJ8017.clean.1m.3.3dg",
            outpng,
            "mm10",
            tmpdir=os.path.join(os.path.dirname(__file__), "output"),
            conda=None
            )
        self.assertTrue(Path(outpng).exists())
    def test_clip_centelo_pymol_dip(self):
        print("Test_clip_centelo_pymol_dip")
        outpng = os.path.join(os.path.dirname(__file__), "output", "clip_centelo_pymol_dip.png")
        # centelo
        clip_centelo_pymol(
            #"/sharec/ychi/repo/sperm56_mESC/3dg_c/ESCO1008.clean.20k.5.3dg",
            "/sharec/ychi/repo/sperm56_mESC/3dg_c/ESCO1008.clean.1m.5.3dg",
            outpng,
            "mm10",
            tmpdir=os.path.join(os.path.dirname(__file__), "output"),
            conda=None,
            dupref=True
            )
    def test_clip_single_centelo_pymol(self):
        print("Test_clip_single_centelo_pymol")
        outpng = os.path.join(os.path.dirname(__file__), "output", "clip_single_centelo_pymol.png")
        # centelo
        clip_single_centelo_pymol(
            "/shareb/ychi/repo/sperm43/3dg_c/BJ8017.clean.1m.3.3dg",
            outpng,
            ["chrX","chrY"],
            "mm10",
            clip=1,
            tmpdir=os.path.join(os.path.dirname(__file__), "output"),
            conda=None
            )
        self.assertTrue(Path(outpng).exists())
if __name__ == "__main__":
    unittest.main()