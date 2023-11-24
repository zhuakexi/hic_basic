# unittest
import sys
sys.path.insert(0, "/share/home/ychi/dev/hic_basic")
sys.path.insert(0, "/share/home/ychi/dev/hires_utils")
import os
import unittest
def test_plot_IS(self):
    """
    Test the plot_IS function.
    """
    out_png = os.path.join(os.path.dirname(__file__), "output", "IS.png")
    pass