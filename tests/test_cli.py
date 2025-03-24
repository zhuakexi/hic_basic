import unittest
import subprocess
import os

class TestCliRender(unittest.TestCase):
    def setUp(self):
        # create tests/output/test_cli_render
        # use __file__ to get the current file path
        # clean up the output directory before running the test
        self.input_file = os.path.join(os.path.dirname(__file__), "data/test_cli_render.tsv")
        self.output_dir = os.path.join(os.path.dirname(__file__), "output", "test_cli_render")
        os.makedirs(self.output_dir, exist_ok=True)
        # Clean up the output directory before running the test
        for file in os.listdir(self.output_dir):
            file_path = os.path.join(self.output_dir, file)
            if os.path.isfile(file_path):
                os.remove(file_path)
        
    def test_surface_territory(self):
        """
        Test executing the 'hic_basic render' command with specified arguments.
        Verify the command runs successfully and outputs the expected files.
        """
        # Command to execute (split into list for safety)
        command = [
            "python",
            "-m",
            "hic_basic",
            "render",
            "-i",
            self.input_file,
            "--colname",
            "200k_g_struct1",
            "-o",
            self.output_dir,
            "--view",
            "rand",
            "--clips",
            "surf",
            "--mode",
            "territory",
            "--force"
        ]

        # Execute the command
        result = subprocess.run(
            command,
            capture_output=True,  # Capture stdout/stderr
            text=True,            # Return outputs as strings
            check=False           # Do not auto-raise exceptions (handle manually)
        )

        # Verify command execution succeeded (return code 0)
        self.assertEqual(
            result.returncode, 
            0, 
            f"Command failed with error: {result.stderr}"
        )

        # Verify output files exist (adjust based on expected outputs)
        expected_output_files = [
            os.path.join(self.output_dir, "b1.randcsurfsNone.png")
        ]
        for file in expected_output_files:
            self.assertTrue(
                os.path.exists(file),
                f"Output file {file} not found after execution"
            )

    def tearDown(self):
        """
        Clean up generated files after test execution.
        """
        # Remove generated files (uncomment if needed)
        # os.remove("test_rotation.png")
        # os.remove("test_rotation.log")

if __name__ == "__main__":
    unittest.main()