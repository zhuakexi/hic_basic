import os
import unittest
from pathlib import Path
from hic_basic.plot.utils import expand_colors, hex_to_rgb, plot_color_sequence

class TestExpandColors(unittest.TestCase):
    def setUp(self):
        # Input colors in hex format
        self.input_colors = [
            '#3366CC', '#DC3912', '#FF9900', '#109618', 
            '#990099', '#0099C6', '#DD4477', '#66AA00', 
            '#B82E2E', '#316395'
        ]
        # Convert hex colors to RGB tuples
        self.rgb_colors = [hex_to_rgb(color) for color in self.input_colors]
        self.output_dir = Path(os.path.dirname(__file__)) / "output" / "test_plot_utils"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.output_file = self.output_dir / "test_expand_color_sequence.png"

    def test_expand_colors_no_interpolation(self):
        # Test with no interpolation
        expanded = expand_colors(self.rgb_colors, n_interpolate=0)
        self.assertEqual(len(expanded), len(self.rgb_colors))
        self.assertEqual(expanded, self.rgb_colors)

    def test_expand_colors_with_interpolation(self):
        # Test with 2 interpolated colors between each pair
        n_interpolate = 2
        expanded = expand_colors(self.rgb_colors, n_interpolate=n_interpolate)
        expected_length = len(self.rgb_colors) + (len(self.rgb_colors) - 1) * n_interpolate
        self.assertEqual(len(expanded), expected_length)

        # Check that the first and last colors remain unchanged
        self.assertEqual(expanded[0], self.rgb_colors[0])
        self.assertEqual(expanded[-1], self.rgb_colors[-1])

        # Visualize the expanded colors
        plot_color_sequence(expanded, output_file=self.output_file)
        self.assertTrue(self.output_file.exists())

class TestPlotColorSequence(unittest.TestCase):
    def setUp(self):
        # Input colors in hex format
        self.input_colors = [
            '#3366CC', '#DC3912', '#FF9900', '#109618', 
            '#990099', '#0099C6', '#DD4477', '#66AA00', 
            '#B82E2E', '#316395'
        ]
        self.output_dir = Path(os.path.dirname(__file__)) / "output" / "test_plot_utils"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.output_file = self.output_dir / "test_color_sequence.png"

    def test_plot_color_sequence(self):
        # Test if the function saves the file correctly
        plot_color_sequence(self.input_colors, output_file=self.output_file)
        self.assertTrue(self.output_file.exists())

if __name__ == '__main__':
    unittest.main()
