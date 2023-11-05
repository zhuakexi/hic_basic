# unittest
import sys
sys.path.insert(0, "/share/home/ychi/dev/hic_basic")
sys.path.insert(0, "/share/home/ychi/dev/hires_utils")
import os
import unittest

import pandas as pd
from hic_basic.territory import intermingle
def test_intermingle(self):
    """
    Test the intermingle function.
    """
    frequency_table = intermingle(
        os.path.join(os.path.dirname(__file__), "intermingle.3dg"),
        0.21,
        table=True
    )
    answer = pd.DataFrame(
        {
            "chr1(mat)" : [2,3,4,4,4,4,3,2,4,3,2,1],
            "chr1(pat)" : [0,0,0,0,1,2,3,4,2,3,3,2]
        }
    )
    pd.testing.assert_frame_equal(frequency_table, answer)

if __name__ == "__main__":
    unittest.main()