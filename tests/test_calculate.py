# tests/test_decorators.py
import unittest
import pandas as pd
import tempfile
import os
from pathlib import Path
import sys

# Import the module containing the mt decorator
sys.path.append(str(Path(__file__).parent.parent))
from hic_basic.calculate import mt

# Import test functions
sys.path.append(str(Path(__file__).parent))

class TestMtDecorator(unittest.TestCase):
    """Test cases for the mt decorator functionality."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.TemporaryDirectory()
        self.test_dir = Path(self.temp_dir.name)
        
        # Create some test files
        (self.test_dir / "input1.txt").write_text("test content 1")
        (self.test_dir / "input2.txt").write_text("test content 2")
        (self.test_dir / "input3.txt").write_text("test content 3")
        
        # Create test DataFrame
        self.test_data = pd.DataFrame({
            'sample_id': ['sample1', 'sample2', 'sample3'],
            'input_col1': [str(self.test_dir / "input1.txt"), 
                          str(self.test_dir / "input2.txt"), 
                          str(self.test_dir / "input3.txt")],
            'input_col2': ['extra1', 'extra2', 'extra3']
        })
        self.test_data.set_index('sample_id', inplace=True)
    
    def tearDown(self):
        """Clean up after tests."""
        self.temp_dir.cleanup()
    
    def test_basic_functionality(self):
        """
        Test basic functionality with input_cols and output_pattern.
        
        Input: DataFrame with input columns and output pattern
        Expected: Function creates output files with correct content
        """
        # Apply decorator to the simple copy function
        def mt_simple_copy():
            pass
        decorated_func = mt("test_functions")(mt_simple_copy)
        
        # Call decorated function
        result = decorated_func(
            self.test_data,
            input_cols=["input_col1"],
            output_pattern=str(self.test_dir / "output_{sample_name}.txt"),
            force=True,
            nproc=1
        )
        
        # Check that output files were created with correct content
        for sample_name in self.test_data.index:
            output_path = self.test_dir / f"output_{sample_name}.txt"
            self.assertTrue(output_path.exists())
            
            # Check content matches input
            input_path = self.test_data.loc[sample_name, "input_col1"]
            with open(input_path, 'r') as f:
                expected_content = f.read()
            
            with open(output_path, 'r') as f:
                actual_content = f.read()
            
            self.assertEqual(actual_content, expected_content)
    
    def test_multiple_input_cols(self):
        """
        Test functionality with multiple input columns.
        
        Input: DataFrame with multiple input columns
        Expected: Function receives all input columns and creates output files
        """
        # Apply decorator to the multiple inputs function
        def mt_process_multiple_inputs():
            pass
        decorated_func = mt("test_functions")(mt_process_multiple_inputs)
        
        # Call decorated function
        result = decorated_func(
            self.test_data,
            input_cols=["input_col1", "input_col2"],
            output_pattern=str(self.test_dir / "output_{sample_name}.txt"),
            force=True,
            nproc=1
        )
        
        # Check that output files were created with correct content
        for sample_name in self.test_data.index:
            output_path = self.test_dir / f"output_{sample_name}.txt"
            self.assertTrue(output_path.exists())
            
            # Check content combines both inputs
            input1 = self.test_data.loc[sample_name, "input_col1"]
            input2 = self.test_data.loc[sample_name, "input_col2"]
            
            with open(output_path, 'r') as f:
                content = f.read()
            
            self.assertEqual(content, f"{input1}_{input2}")
    
    def test_mixed_input_sources(self):
        """
        Test functionality with mixed input sources (cols and patterns).
        
        Input: DataFrame with input columns and input patterns
        Expected: Function receives inputs in correct order and creates output files
        """
        # Apply decorator to the multiple inputs function
        def mt_process_multiple_inputs():
            pass
        decorated_func = mt("test_functions")(mt_process_multiple_inputs)
        
        # Call decorated function
        result = decorated_func(
            self.test_data,
            input_cols=["input_col1"],
            input_pattern="extra_pattern_{sample_name}",
            output_pattern=str(self.test_dir / "output_{sample_name}.txt"),
            force=True,
            nproc=1
        )
        
        # Check that output files were created with correct content
        for sample_name in self.test_data.index:
            output_path = self.test_dir / f"output_{sample_name}.txt"
            self.assertTrue(output_path.exists())
            
            # Check content combines both inputs (col first, then pattern)
            col_input = self.test_data.loc[sample_name, "input_col1"]
            pattern_input = f"extra_pattern_{sample_name}"
            
            with open(output_path, 'r') as f:
                content = f.read()
            
            self.assertEqual(content, f"{col_input}_{pattern_input}")
    
    def test_multiple_output_sources(self):
        """
        Test functionality with multiple output sources.
        
        Input: DataFrame with output columns and output patterns
        Expected: Function creates multiple output files
        """
        # Add output columns to test data
        self.test_data['output_col1'] = [str(self.test_dir / "col_out1.txt"), 
                                        str(self.test_dir / "col_out2.txt"), 
                                        str(self.test_dir / "col_out3.txt")]
        
        # Apply decorator to the simple copy function
        def mt_process_multiple_outputs():
            pass
        decorated_func = mt("test_functions")(mt_process_multiple_outputs)
        
        # Call decorated function
        result = decorated_func(
            self.test_data,
            input_cols=["input_col1"],
            output_cols=["output_col1"],
            output_pattern=str(self.test_dir / "pattern_out_{sample_name}.txt"),
            force=True,
            nproc=1
        )
        
        # Check that both column-based and pattern-based output files were created
        for sample_name in self.test_data.index:
            # Check column-based output
            col_output_path = Path(self.test_data.loc[sample_name, "output_col1"])
            self.assertTrue(col_output_path.exists())
            
            # Check pattern-based output
            pattern_output_path = self.test_dir / f"pattern_out_{sample_name}.txt"
            self.assertTrue(pattern_output_path.exists())
    
    def test_force_parameter(self):
        """
        Test the force parameter behavior.
        
        Input: Existing output files and force=False
        Expected: Function calls are skipped for existing outputs
        """
        # Create output files in advance
        (self.test_dir / "output_sample1.txt").write_text("existing content")
        
        # Apply decorator to the simple copy function
        def mt_simple_copy():
            pass
        decorated_func = mt("test_functions")(mt_simple_copy)
        
        # Call with force=False (should skip existing outputs)
        result = decorated_func(
            self.test_data,
            input_cols=["input_col1"],
            output_pattern=str(self.test_dir / "output_{sample_name}.txt"),
            force=False,  # Don't force overwrite
            nproc=1
        )
        
        # Check that the existing file was not overwritten
        with open(self.test_dir / "output_sample1.txt", 'r') as f:
            content = f.read()
        self.assertEqual(content, "existing content")  # Should not be overwritten
        
        # Check that other files were created
        self.assertTrue((self.test_dir / "output_sample2.txt").exists())
        self.assertTrue((self.test_dir / "output_sample3.txt").exists())
    
    def test_concat_parameter(self):
        """
        Test the concat parameter behavior.
        
        Input: Function returns pandas Series and concat=True
        Expected: Results are concatenated into a DataFrame
        """
        # Apply decorator to the series-returning function
        def mt_return_series():
            pass
        decorated_func = mt("test_functions")(mt_return_series)
        
        # Call with concat=True
        result = decorated_func(
            self.test_data,
            input_cols=["input_col1"],
            output_pattern=str(self.test_dir / "output_{sample_name}.txt"),
            force=True,
            nproc=1,
            concat=True
        )
        
        # Result should be a DataFrame
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(result.shape[0], 3)  # 3 rows
        self.assertEqual(result.shape[1], 3)  # 3 columns (input_file, output_file, content_length)
        
        # Check that output files were created
        for sample_name in self.test_data.index:
            output_path = self.test_dir / f"output_{sample_name}.txt"
            self.assertTrue(output_path.exists())
    
    def test_error_handling(self):
        """
        Test error handling in decorated functions.
        
        Input: Function that raises an exception
        Expected: Exception is caught and error message is printed, other tasks continue
        TODO: Improve this test to actually check stderr output, now Error info is blocked by tqdm
        """
        # Apply decorator to the failing function
        def mt_failing_function():
            pass
        decorated_func = mt("test_functions")(mt_failing_function)
        
        # Capture stdout to check error messages
        import io
        from contextlib import redirect_stderr, redirect_stdout
        
        stderr_capture = io.StringIO()
        
        with redirect_stderr(stderr_capture):
            # Should not raise an exception
            result = decorated_func(
                self.test_data,
                input_cols=["input_col1"],
                output_pattern=str(self.test_dir / "output_{sample_name}.txt"),
                force=True,
                nproc=1
            )
        
        # Check that error was printed to stderr
        stderr_output = stderr_capture.getvalue()
        # self.assertIn("ERROR", stderr_output)
        # self.assertIn("This function always fails", stderr_output)
        
        # Result should be empty (all tasks failed)
        self.assertEqual(len(result), 0)


if __name__ == '__main__':
    unittest.main()