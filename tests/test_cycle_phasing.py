import unittest
import pandas as pd
import os
from hic_basic.cycle_phasing import dis_counts

class TestCyclePhasing(unittest.TestCase):
    def test_dis_counts(self):
        pairs_fp = os.path.join(os.path.dirname(__file__), 'data/mm10.c123.pairs.gz')
        result = dis_counts(pairs_fp)
        self.assertIsInstance(result, pd.Series)
        self.assertEqual(result.name, pairs_fp)
        self.assertEqual(len(result), 151)
        self.assertTrue(result.notna().all())
    def test_dis_count_nan(self):
        pairs_fp = os.path.join(os.path.dirname(__file__), 'data/not_exist.pairs.gz')
        result = dis_counts(pairs_fp)
        print(result)
        self.assertIsInstance(result, pd.Series)
        self.assertEqual(result.name, pairs_fp)
        self.assertEqual(len(result), 151)
        self.assertTrue(result.isna().all())
if __name__ == '__main__':
    unittest.main()
