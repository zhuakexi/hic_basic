import sys
import unittest

sys.path.insert(0, "/share/home/ychi/dev/hic_basic")
sys.path.insert(0, "/share/home/ychi/dev/hires_utils")

from hic_basic.genome import Region

class TestRegion(unittest.TestCase):
    def test_list_of_tuples_input(self):
        """Test the Region class with a list of tuples input."""
        input_region = [('chr1', 1000000), ('chr1', 2000000)]
        region = Region(input_region)
        expected_output = [('chr1', 1000000), ('chr1', 2000000)]
        self.assertEqual(region.r, expected_output)

    def test_list_of_tuples_without_position(self):
        """Test input where one of the positions is not provided, defaults should be handled."""
        # Assuming your chromosomes function and how it fetches length is correctly mocked or implemented
        input_region = [('chr1', 1000000), ('chr1',)]
        region = Region(input_region)
        expected_output = [('chr1', 1000000), ('chr1', 248956422)]  # Assuming 'chr1' end is this
        self.assertEqual(region.r, expected_output)

    def test_incorrect_input_type(self):
        """Test with incorrect input type (string), which should not be implemented yet."""
        with self.assertRaises(NotImplementedError):
            Region("chr1:1000000-2000000")

    def test_incorrect_list_length(self):
        """Test with incorrect number of elements in the input list."""
        input_region = [('chr1', 1000000)]  # Only one tuple, should be two
        with self.assertRaises(AssertionError):
            Region(input_region)

    def test_region_with_bins(self):
        """Test with binsize provided, should generate bins indexes for the region."""
        input_region = [('chr1', 248000000), ('chr2', 2000000)]
        region = Region(input_region, binsize=100000)
        print(region.bins)
        self.assertTrue(region.bins)
# Run the tests
if __name__ == '__main__':
    unittest.main()
