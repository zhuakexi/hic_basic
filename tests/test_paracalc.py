import unittest
import pandas as pd
import os
from hic_basic.wet.paracalc import gen_cdps

class TestParacalc(unittest.TestCase):
    def test_gen_cdps(self):
        data = {
            'pairs_c123': [
                os.path.join(os.path.dirname(__file__), 'data/mm10.c123.pairs.gz'),
                'data/not_exist.pairs.gz'
            ]
        }
        filesp = pd.DataFrame(data, index=['sample1', 'sample2'])
        result = gen_cdps(filesp)
        #print(result)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(result.shape, (2, 151))
        self.assertTrue(result.loc['sample1'].notna().all())
        self.assertTrue(result.loc['sample2'].isna().all())

if __name__ == '__main__':
    unittest.main()
