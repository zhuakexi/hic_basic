import os
import sys
sys.path.insert(0, "/share/home/ychi/dev/hic_basic")
sys.path.insert(0, "/share/home/ychi/dev/hires_utils")
sys.path.insert(0, "/share/home/ychi/dev/sperm_struct")
import unittest
from hic_basic.data import fetch_cent_chromlen

class TestData(unittest.TestCase):
    """
    Test the data module.
    """
    def test_fetch_cent_chromlen(self):
        """
        Test the fetch_cent_chromlen function.
        """
        # fetch_cent_chromlen
        cent_chromlen = fetch_cent_chromlen("mm10")
        print(cent_chromlen["chrY"])
        self.assertTrue("chrY" in cent_chromlen.keys())
if __name__ == "__main__":
    unittest.main()