import unittest
import os

import numpy as np
from hic_basic.hicio import parse_3dg
from hic_basic.structure.measure import angle_between_vectors, perpendicular_unit_vector, cent2telo_vector, pm_vector

class TestStructureMeasure(unittest.TestCase):

    def setUp(self):
        self.output_dir = "output"
        os.makedirs(self.output_dir, exist_ok=True)
        self._3dg_file = self.get_data_file_path("data/sample.mm10_dip.pair_no1.3dg.gz")
        self.genome = "mm10_dip"
        self.p, self.q = False, True
        self._3dg = parse_3dg(self._3dg_file)

    def get_data_file_path(self, relative_path):
        base_dir = os.path.dirname(os.path.abspath(__file__))
        return os.path.join(base_dir, relative_path)

    def test_cent2telo_vector(self):
        res, c2t_vector = cent2telo_vector(
            self._3dg, genome=self.genome, p=self.p, q=self.q, 
            show_chroms=True)

    def test_pm_vector(self):
        # ...existing test code...
        res = pm_vector(self._3dg)
        print(res)
        # result = pm_vector(...)
        # self.assertIsNotNone(result)
        # # Save result to output directory
        # with open(os.path.join(self.output_dir, "pm_vector_result.txt"), "w") as f:
        #     f.write(str(result))

    def test_angle_between_vectors(self):
        vector_a = [1, 0, 0]
        vector_b = [0, 1, 0]
        angle = angle_between_vectors(vector_a, vector_b)
        self.assertAlmostEqual(angle, 90.0)

        vector_c = [1, 0, 0]
        vector_d = [1, 1, 0]
        angle = angle_between_vectors(vector_c, vector_d)
        self.assertAlmostEqual(angle, 45.0)

    def test_perpendicular_unit_vector(self):
        vector_a = [0, 1, 0]
        vector_b = [1, 0, 0]
        result = perpendicular_unit_vector(vector_a, vector_b)
        expected_result = [0, 0, -1]
        np.testing.assert_almost_equal(result, expected_result)

    def tearDown(self):
        # Clean up output directory if needed
        # for file in os.listdir(self.output_dir):
        #     os.remove(os.path.join(self.output_dir, file))
        pass

if __name__ == "__main__":
    unittest.main()