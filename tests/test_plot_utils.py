import os
import unittest
from pathlib import Path
from hic_basic.plot.utils import expand_color_sequence, hex_to_rgb, plot_color_sequence, interpolate_color_rgb_linear

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
        expanded = expand_color_sequence(self.rgb_colors, n=0)
        # self.assertEqual(len(expanded), len(self.rgb_colors))
        # self.assertEqual(expanded, self.rgb_colors)

    def test_expand_colors_with_interpolation(self):
        # Test with 2 interpolated colors between each pair
        n_interpolate = 2
        expanded = expand_color_sequence(self.rgb_colors, n=n_interpolate)
        expected_length = len(self.rgb_colors) + len(self.rgb_colors) * n_interpolate
        # Check if the length of the expanded list is as expected
        self.assertEqual(len(expanded), expected_length)

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

class TestInterpolateColor(unittest.TestCase):
    def setUp(self):
        self.color1 = '#3366CC'
        self.rgb1 = hex_to_rgb(self.color1)
        self.rgb2 = (255, 255, 255) 
        self.output_dir = Path(os.path.dirname(__file__)) / "output" / "test_plot_utils"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.output_file = self.output_dir / "test_interpolate_color_rgb_white.png"

    def test_interpolation(self):
        # Interpolation factors
        expected_results = interpolate_color_rgb_linear(self.rgb1, self.rgb2, 5, cap=0.75)
        plot_color_sequence(expected_results, output_file=self.output_file)

if __name__ == '__main__':
    unittest.main()
