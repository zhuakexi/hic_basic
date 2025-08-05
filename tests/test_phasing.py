import unittest

import pandas as pd
from hic_basic.phasing import mt_count_chrom_phased

class TestPhasing(unittest.TestCase):
    def setUp(self):
        # This method will run before each test
        self.df = pd.DataFrame(
            index = ['F1T14154N010', 'F1T1415N173', 'F1T14154N001', 'F1T14152N092', 'B6D2F1T14092']
        )
    def test_count_chrom_phased(self):
        # Test with a sample .seg.gz file
        result = mt_count_chrom_phased(
            self.df,
            input_pattern = "/shared/ychi/repo/mouse_spg_rerun/pre_seg/{sample_name}.seg.gz",
            output_pattern = None
        )
        print(result)
if __name__ == '__main__':
    unittest.main()