import unittest
import os
import pandas as pd
from pathlib import Path
from tempfile import TemporaryDirectory
from hic_basic.ana import Ana, load_json, dump_json  # Replace with your actual module name

class TestAna(unittest.TestCase):
    def setUp(self):
        """
        Create a temporary directory for testing.
        """
        self.temp_dir = TemporaryDirectory()
        self.home = Path(self.temp_dir.name)

    def tearDown(self):
        """
        Clean up the temporary directory after each test.
        """
        self.temp_dir.cleanup()

    def test_init_without_clear(self):
        """
        Test initialization with clear=False retains existing files.
        """
        data_path = self.home / "data.csv.gz"
        obj_path = self.home / "obj.json"
        for path in [data_path, obj_path]:
            path.touch()

        Ana(self.home, clear=False)
        self.assertTrue(data_path.exists())
        self.assertTrue(obj_path.exists())

    def test_init_with_clear(self):
        """
        Test initialization with clear=True removes existing files.
        """
        # Create some existing files
        data_path = self.home / "data.csv.gz"
        obj_path = self.home / "obj.json"
        commit_file = self.home / "commit_0001.data.csv.gz"
        
        for path in [data_path, obj_path, commit_file]:
            path.touch()
            self.assertTrue(path.exists())

        # Initialize with clear=True
        Ana(self.home, clear=True)
        
        # Data files should be recreated
        self.assertTrue(data_path.exists())
        self.assertTrue(obj_path.exists())
        # Commit files should be removed
        self.assertFalse(commit_file.exists())

    def test_init_with_data(self):
        """
        Test initialization with data parameter.
        """
        df = pd.DataFrame({"A": [1, 2], "B": [3, 4]})
        ana = Ana(self.home, data=df)
        df.index = df.index.astype("string")
        pd.testing.assert_frame_equal(ana.data, df)

    def test_init_with_obj(self):
        """
        Test initialization with obj parameter (dict of key-value lists).
        """
        obj = {"key1": [1, 2], "key2": [3, 4]}
        ana = Ana(self.home, obj=obj)
        self.assertEqual(ana.obj, obj)

    def test_data_property(self):
        """
        Test data property correctly reads and writes DataFrame.
        """
        df = pd.DataFrame({"A": [1, 2], "B": [3, 4]})
        ana = Ana(self.home, data=df)
        df.index = df.index.astype("string")
        pd.testing.assert_frame_equal(ana.data, df)

    def test_update_dataframe(self):
        """
        Test update method with DataFrame input.
        """
        df1 = pd.DataFrame({"A": [1, 2]}, index=["x", "y"])
        df2 = pd.DataFrame({"A": [3], "B": [4]}, index=["z"])
        ana = Ana(self.home, data=df1)
        ana.update(df2)
        expected = pd.concat([df1, df2], axis=0).groupby(level=0).last()
        expected.index = expected.index.astype("string")
        pd.testing.assert_frame_equal(ana.data, expected)

    def test_update_series(self):
        """
        Test update method with Series input.
        """
        ser = pd.Series([1, 2], index=["x", "y"], name="ser")
        ana = Ana(self.home)
        ana.update(ser, key="ser")
        expected = ser.to_frame()
        expected.index = expected.index.astype("string")
        # print("ana.data", ana.data)
        # print("expected", expected)
        pd.testing.assert_frame_equal(ana.data, expected)

    def test_update_dict(self):
        """
        Test update method with dict input.
        """
        d = {"a": 1, "b": 2}
        ana = Ana(self.home, clear=True)
        ana.update(d, key="dict_key")
        expected = pd.DataFrame({"dict_key": pd.Series(d)})
        expected.index = expected.index.astype("string")
        pd.testing.assert_frame_equal(ana.data, expected)

    def test_update_list(self):
        """
        Test update method with list input.
        """
        lst = [1, 2, 3]
        ana = Ana(self.home)
        ana.update(lst, key="list_key")
        self.assertEqual(ana.obj["list_key"], lst)

    def test_unsupported_type(self):
        """
        Test update method raises TypeError for unsupported data types.
        """
        ana = Ana(self.home)
        with self.assertRaises(TypeError):
            ana.update(object())

    def test_required_key_for_dict_and_list(self):
        """
        Test dict and list inputs require key to be provided.
        """
        ana = Ana(self.home)
        with self.assertRaises(AssertionError):
            ana.update([1, 2])  # Missing key
        with self.assertRaises(AssertionError):
            ana.update({"a": 1})  # Missing key

    def test_commit(self):
        """
        Test commit method creates commit files and updates metadata.
        """
        ana = Ana(self.home, clear=True, max_commits=10)
        df = pd.DataFrame({"A": [1, 2]}, index=["x", "y"])
        ana.update(df)
        
        commit_id = ana.commit("Test commit")
        
        # Check that files were created with correct naming pattern
        meta = load_json(ana.commit_meta_path)
        commit_entry = next(c for c in meta["commits"] if c["id"] == commit_id)
        filename = commit_entry["filename"]
        
        self.assertTrue((self.home / f"{filename}.data.csv.gz").exists())
        self.assertTrue((self.home / f"{filename}.obj.json").exists())
        
        # Verify metadata was updated correctly
        self.assertIn(commit_id, [c["id"] for c in meta["commits"]])
        self.assertEqual(meta["commits"][-1]["id"], commit_id)
        self.assertIn(str(commit_id), meta["file_mapping"])

    def test_circular_file_naming(self):
        """
        Test that commit files use circular naming when max_commits is exceeded.
        """
        max_commits = 5
        ana = Ana(self.home, clear=True, max_commits=max_commits, verbose=False)
        
        # Create more commits than max_commits
        for i in range(max_commits + 3):
            df = pd.DataFrame({"value": [i]}, index=[f"row_{i}"])
            ana.update(df, description=f"Commit {i}")
        
        meta = load_json(ana.commit_meta_path)
        
        # Should have exactly max_commits commits
        self.assertEqual(len(meta["commits"]), max_commits)
        
        # Check that file names are within the circular range
        for commit in meta["commits"]:
            filename = commit["filename"]
            # Extract the number from "commit_XXXX"
            file_number = int(filename.split('_')[1])
            self.assertLess(file_number, max_commits)
            self.assertGreaterEqual(file_number, 0)
            
        # Verify file mapping is correct
        self.assertEqual(len(meta["file_mapping"]), max_commits)
        
        # Check that oldest commits were properly removed
        commit_ids = [c["id"] for c in meta["commits"]]
        # Should have the most recent max_commits commits
        expected_ids = list(range(max_commits + 3 - max_commits + 1, max_commits + 3 + 1))
        self.assertEqual(commit_ids, expected_ids)

    def test_revert_functionality(self):
        """
        Test revert to previous commit.
        """
        ana = Ana(self.home, clear=True, max_commits=10, verbose=False)
        
        # Create initial data
        df1 = pd.DataFrame({"A": [1, 2]}, index=["x", "y"])
        ana.update(df1, description="Initial data")
        commit1_data = ana.data.copy()
        
        # Make a change
        df2 = pd.DataFrame({"A": [3, 4], "B": [5, 6]}, index=["z", "w"])
        ana.update(df2, description="Second update")
        
        # Revert to first commit
        new_commit_id = ana.revert(1)  # Revert to commit ID 1
        
        # Should create a new commit
        self.assertGreater(new_commit_id, 2)
        
        # Data should match the first commit state
        pd.testing.assert_frame_equal(ana.data, commit1_data)

    def test_revert_by_description(self):
        """
        Test revert using commit description.
        """
        ana = Ana(self.home, clear=True, max_commits=10, verbose=False)
        
        # Create initial data with specific description
        df1 = pd.DataFrame({"A": [1, 2]}, index=["x", "y"])
        ana.update(df1, description="First state")
        first_state = ana.data.copy()
        
        # Make changes
        df2 = pd.DataFrame({"A": [3, 4]}, index=["z", "w"])
        ana.update(df2, description="Second state")
        
        # Revert using description
        new_commit_id = ana.revert("First state")
        
        # Should revert to first state
        pd.testing.assert_frame_equal(ana.data, first_state)
        self.assertGreater(new_commit_id, 2)

    def test_commit_description_validation(self):
        """
        Test that commit descriptions are properly validated.
        """
        ana = Ana(self.home, clear=True, verbose=False)
        
        # Create first commit
        df = pd.DataFrame({"A": [1]})
        ana.update(df, description="First commit")
        
        # Should not allow duplicate descriptions
        with self.assertRaises(ValueError):
            ana.commit("First commit")
        
        # Should not allow purely numeric descriptions
        with self.assertRaises(ValueError):
            ana.commit("123")
        with self.assertRaises(ValueError):
            ana.commit(" 456 ")

    def test_file_mapping_persistence(self):
        """
        Test that file mapping is properly persisted and loaded.
        """
        # Create initial commits
        ana1 = Ana(self.home, clear=True, max_commits=5, verbose=False)
        for i in range(3):
            df = pd.DataFrame({"value": [i]})
            ana1.update(df, description=f"Commit {i}")
        
        # Reload from same directory
        ana2 = Ana(self.home, clear=False, max_commits=5, verbose=False)
        
        # File mapping should be preserved
        meta = load_json(ana2.commit_meta_path)
        self.assertIn("file_mapping", meta)
        self.assertEqual(len(meta["file_mapping"]), 3)
        
        # Should be able to revert using the mapping
        ana2.revert(1)  # Should work via file mapping

    def test_force_parameter(self):
        """
        Test the force parameter for overwriting existing data.
        """
        ana = Ana(self.home, clear=True, verbose=False)
        
        # First update without force
        df1 = pd.DataFrame({"A": [1, 2]}, index=["x", "y"])
        ana.update(df1, key="data", force=False)
        initial_shape = ana.data.shape
        
        # Try to update same key without force - should be skipped
        df2 = pd.DataFrame({"A": [3, 4]}, index=["z", "w"])
        ana.update(df2, key="data", force=False)
        self.assertEqual(ana.data.shape, initial_shape)  # No change
        
        # Update with force - should overwrite
        ana.update(df2, key="data", force=True)
        self.assertEqual(ana.data.shape[0], 2)  # Only new data

    def test_max_commits_parameter(self):
        """
        Test that max_commits parameter is respected.
        """
        max_commits = 3
        ana = Ana(self.home, clear=True, max_commits=max_commits, verbose=False)
        
        # Create more commits than max_commits
        for i in range(max_commits + 2):
            df = pd.DataFrame({"value": [i]})
            ana.update(df)
        
        meta = load_json(ana.commit_meta_path)
        self.assertEqual(len(meta["commits"]), max_commits)
        self.assertEqual(meta["max_commits"], max_commits)

    def test_commit_filename_generation(self):
        """
        Test the internal _get_commit_filename method.
        """
        ana = Ana(self.home, clear=True, max_commits=10, verbose=False)
        
        # Test various commit IDs
        test_cases = [
            (1, "commit_0001"),
            (5, "commit_0005"), 
            (10, "commit_0000"),  # 10 % 10 = 0
            (15, "commit_0005"),  # 15 % 10 = 5
            (23, "commit_0003"),  # 23 % 10 = 3
        ]
        
        for commit_id, expected_filename in test_cases:
            filename = ana._get_commit_filename(commit_id)
            self.assertEqual(filename, expected_filename)

if __name__ == "__main__":
    unittest.main()