"""
Unit test template for the hic_basic package.

This template demonstrates how to structure unit tests for the hic_basic package,
including handling temporary files with tempfile and manual output verification
in the tests/output directory.
"""
import unittest
import tempfile
import os
import shutil
from pathlib import Path

# Import the modules you want to test from hic_basic
# Example: from hic_basic.module1 import some_function
#          from hic_basic.module2 import SomeClass


class TestTemplate(unittest.TestCase):
    """
    Template class for unit testing the hic_basic package.
    
    Input:
        - Various test scenarios for functions/classes in hic_basic
        - Input data for testing
        - Expected outputs for validation
    
    Output:
        - Verification of function behavior
        - Correctness of output data
        - Proper handling of edge cases and errors
    """

    def setUp(self):
        """
        Set up test fixtures before each test method.
        
        This method is called before each individual test method.
        Two common scenarios are demonstrated here:
        
        1. Using tempfile for temporary file management:
           - Advantage: Automatically cleaned up after test
           - Good for testing functions that read/write files without needing manual inspection
           - Use tempfile.NamedTemporaryFile() or tempfile.mkstemp()
        
        2. Using tests/output directory for persistent output:
           - Advantage: Output files remain for manual inspection
           - Good for testing complex output formats that need verification
           - Files placed in tests/output/[module_name]/ for organization
        """
        
        # Get the directory of this test file to calculate relative paths
        # This ensures the test works regardless of where it's run from
        test_dir = Path(__file__).parent
        self.test_data_dir = test_dir / "data"
        self.test_output_dir = test_dir / "output" / "module1"  # Adjust module name as needed
        
        # Create output directory if it doesn't exist
        self.test_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Scenario 1: Temporary directory for files that should be auto-cleaned
        self.temp_dir = tempfile.mkdtemp(prefix="hic_basic_test_")
        self.temp_file_path = os.path.join(self.temp_dir, "temp_test_file.txt")
        
        # Scenario 2: Persistent output directory for manual inspection
        self.persistent_output_path = self.test_output_dir / "test_output.txt"
        
        # Example of loading test data from tests/data/
        self.test_input_file = self.test_data_dir / "sample_input.txt"
        if self.test_input_file.exists():
            with open(self.test_input_file, 'r') as f:
                self.sample_input_data = f.read()

    def tearDown(self):
        """
        Clean up after each test method.
        
        Remove temporary directories and files created during the test.
        The persistent output files in tests/output/ are intentionally left
        for manual inspection after tests run.
        """
        # Clean up temporary directory
        shutil.rmtree(self.temp_dir)

    def test_example_with_tempfile(self):
        """
        Example test using temporary files.
        
        Input:
            - Uses temporary file/directory created in setUp
            - Sample input data loaded from tests/data/
        
        Output:
            - Verifies function behavior with temporary files
            - Tests that temporary files are properly handled
        """
        # Example: Write sample data to temporary file
        with open(self.temp_file_path, 'w') as temp_file:
            temp_file.write("Test data for temporary file")
        
        # Example: Call function from hic_basic package
        # result = some_function(self.temp_file_path)
        
        # Add assertions here to verify behavior
        # self.assertEqual(result, expected_result)
        
        # The temporary file will be automatically cleaned up in tearDown

    def test_example_with_persistent_output(self):
        """
        Example test using persistent output files for manual inspection.
        
        Input:
            - Sample input data
            - Configuration parameters
        
        Output:
            - Result written to tests/output/[module]/ for manual inspection
            - Verification that output file was created
        """
        # Example: Process data and write output to persistent location
        # You can use the persistent output path defined in setUp
        with open(self.persistent_output_path, 'w') as output_file:
            output_file.write("This output can be manually inspected")
        
        # Verify that the output file was created
        self.assertTrue(self.persistent_output_path.exists())
        
        # You can add more specific assertions based on your module's behavior
        # For example, check file size, content, etc.
        with open(self.persistent_output_path, 'r') as output_file:
            content = output_file.read()
            self.assertIn("manually inspected", content)

    def test_another_function(self):
        """
        Template for testing another function from hic_basic.
        
        Input:
            - Specific parameters for the function being tested
            - Edge cases or error conditions to verify
        
        Output:
            - Function return value validation
            - Exception handling verification
            - Side effects (file creation, etc.) validation
        """
        # Example test implementation:
        # try:
        #     result = some_function_from_hic_basic(param1, param2)
        #     self.assertEqual(result, expected_value)
        # except SomeExpectedException as e:
        #     self.fail(f"Unexpected exception: {e}")
        
        # For functions that should raise exceptions under certain conditions:
        # with self.assertRaises(ExpectedExceptionType):
        #     problematic_function(bad_parameter)
        
        pass  # Replace with actual test logic

    def test_edge_cases(self):
        """
        Test edge cases and error conditions.
        
        Input:
            - Invalid parameters
            - Boundary values
            - Missing input files
            - Malformed input data
        
        Output:
            - Proper exception handling
            - Appropriate error messages
            - Graceful degradation or clear error states
        """
        # Example: Test with invalid input
        # with self.assertRaises(ValueError):
        #     some_function_with_validation(invalid_input)
        
        # Example: Test with empty input
        # result = some_function("")
        # self.assertEqual(result, expected_for_empty_input)
        
        pass  # Replace with actual edge case tests


if __name__ == '__main__':
    # This allows running the test directly with python -m unittest test_module.py
    unittest.main()



