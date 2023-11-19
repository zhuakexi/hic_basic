# unittest
import sys
sys.path.insert(0, "/share/home/ychi/dev/hic_basic")
sys.path.insert(0, "/share/home/ychi/dev/hires_utils")
import os
import unittest
from pathlib import Path

from hic_basic.plot.hic import plot_compartment

def test_plot_compartment(self):
    """
    Test the plot_IS function.
    """
    out_png = os.path.join(os.path.dirname(__file__), "output", "plot_compartment.png")
    ddir = Path("/shareb/ychi/repo/sperm_struct/notebooks/Compartment")
    fig = plot_compartment(
        str(ddir / "RS.20k.pileup.s14M.1M.cool"),
        ddir / "RS.20k.pileup.s14M.1M.cis.vecs.tsv",
        region="chr1",
        title = "Compartment",
    )
    fig.savefig(out_png, dpi=300)
if __name__ == "__main__":
    unittest.main()