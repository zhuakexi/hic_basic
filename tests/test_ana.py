import unittest
import os
import hashlib
import pandas as pd
import numpy as np
from pathlib import Path
from tempfile import TemporaryDirectory
from hic_basic.ana import Ana, AnaUpdateError, load_json, dump_json  # Replace with your actual module name

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

    def _sample_index(self, df):
        df.index = pd.Index([f"sample_{i}" for i in range(len(df))], dtype=pd.StringDtype())
        return df

    def _home_file_bytes(self):
        return {
            path.name: path.read_bytes()
            for path in self.home.iterdir()
            if path.is_file()
        }

    @staticmethod
    def _combined_hash(paths):
        digest = hashlib.sha256()
        for path in paths:
            digest.update(path.read_bytes())
        return digest.hexdigest()

    def test_init_without_clear(self):
        """
        Test initialization with clear=False retains existing files.
        """
        data_path = self.home / "data.csv"
        schema_path = self.home / "data_schema.json"
        obj_path = self.home / "obj.json"
        for path in [data_path, schema_path, obj_path]:
            path.touch()

        Ana(self.home, clear=False)
        self.assertTrue(data_path.exists())
        self.assertTrue(schema_path.exists())
        self.assertTrue(obj_path.exists())

    def test_init_with_clear(self):
        """
        Test initialization with clear=True removes existing files.
        """
        # Create some existing files
        data_path = self.home / "data.csv"
        schema_path = self.home / "data_schema.json"
        obj_path = self.home / "obj.json"
        commit_file = self.home / "commit_0001.data.csv"
        
        for path in [data_path, schema_path, obj_path, commit_file]:
            path.touch()
            self.assertTrue(path.exists())

        # Initialize with clear=True
        Ana(self.home, clear=True)
        
        # Data files should be recreated
        self.assertTrue(data_path.exists())
        self.assertTrue(schema_path.exists())
        self.assertTrue(obj_path.exists())
        # Commit files should be removed
        self.assertFalse(commit_file.exists())

    def test_init_with_data(self):
        """
        Test initialization with data parameter.
        """
        df = self._sample_index(pd.DataFrame({"A": [1, 2], "B": [3, 4]}))
        ana = Ana(self.home, data=df)
        # With dtype preservation, the index type should be preserved
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
        df = self._sample_index(pd.DataFrame({"A": [1, 2], "B": [3, 4]}))
        ana = Ana(self.home, data=df)
        # With dtype preservation, the index type should be preserved
        pd.testing.assert_frame_equal(ana.data, df)

    def test_update_dataframe(self):
        """
        Test update method with DataFrame input.
        """
        df1 = pd.DataFrame({"A": [1, 2]}, index=pd.Index(["x", "y"], dtype=pd.StringDtype()))
        df2 = pd.DataFrame({"A": [3], "B": [4]}, index=pd.Index(["z"], dtype=pd.StringDtype()))
        ana = Ana(self.home, data=df1)
        ana.update(df2, intent="add_samples")
        expected = pd.concat([df1, df2], axis=0).groupby(level=0).last()
        # With dtype preservation, the index type should be preserved
        pd.testing.assert_frame_equal(ana.data, expected)

    def test_update_series(self):
        """
        Test update method with Series input.
        """
        ser = pd.Series([1, 2], index=pd.Index(["x", "y"], dtype=pd.StringDtype()), name="ser")
        ana = Ana(self.home)
        ana.update(ser, key="ser", intent="add_samples")
        expected = ser.to_frame()
        # With dtype preservation, the index type should be preserved
        pd.testing.assert_frame_equal(ana.data, expected)

    def test_update_dict(self):
        """
        Test update method with dict input.
        """
        d = {"a": 1, "b": 2}
        ana = Ana(self.home, clear=True)
        ana.update(d, key="dict_key", intent="add_samples")
        expected = pd.DataFrame({"dict_key": pd.Series(d)})
        expected.index = pd.Index(expected.index, dtype=pd.StringDtype())
        # With dtype preservation, the index type should be preserved
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
        ana.update(df, intent="add_samples")
        
        commit_id = ana.commit("Test commit")
        
        # Check that files were created with correct naming pattern
        meta = load_json(ana.commit_meta_path)
        commit_entry = next(c for c in meta["commits"] if c["id"] == commit_id)
        filename = commit_entry["filename"]
        
        self.assertTrue((self.home / f"{filename}.data.csv").exists())
        self.assertTrue((self.home / f"{filename}.data_schema.json").exists())
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
            ana.update(df, description=f"Commit {i}", intent="add_samples")
        
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
        ana.update(df1, description="Initial data", intent="add_samples")
        commit1_data = ana.data.copy()
        
        # Make a change
        df2 = pd.DataFrame({"A": [3, 4], "B": [5, 6]}, index=["z", "w"])
        ana.update(df2, description="Second update", intent="add_samples")
        
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
        ana.update(df1, description="First state", intent="add_samples")
        first_state = ana.data.copy()
        
        # Make changes
        df2 = pd.DataFrame({"A": [3, 4]}, index=["z", "w"])
        ana.update(df2, description="Second state", intent="add_samples")
        
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
        ana.update(pd.DataFrame({"A": [1]}, index=["sample_1"]), description="First commit", intent="add_samples")
        
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
            df = pd.DataFrame({"value": [i]}, index=[f"sample_{i}"])
            ana1.update(df, description=f"Commit {i}", intent="add_samples")
        
        # Reload from same directory
        ana2 = Ana(self.home, clear=False, max_commits=5, verbose=False)
        
        # File mapping should be preserved
        meta = load_json(ana2.commit_meta_path)
        self.assertIn("file_mapping", meta)
        self.assertEqual(len(meta["file_mapping"]), 3)
        
        # Should be able to revert using the mapping
        ana2.revert(1)  # Should work via file mapping

    def test_idempotent_update(self):
        """
        Test that update() is idempotent: calling with the same data
        does not create additional commits.
        """
        ana = Ana(self.home, clear=True, verbose=False)
        df = pd.DataFrame({"A": [1, 2]}, index=["x", "y"])

        # First update should create a commit
        cid1 = ana.update(df, intent="add_samples")
        self.assertEqual(cid1, 1)

        # Second update with identical data should skip commit
        cid2 = ana.update(df)
        self.assertEqual(cid2, cid1)  # Same commit ID returned

        # Verify no new commits were created
        meta = load_json(ana.commit_meta_path)
        self.assertEqual(len(meta["commits"]), 1)

    def test_force_override_idempotency(self):
        """
        Test that force=True overrides idempotency and creates a new commit.
        """
        ana = Ana(self.home, clear=True, verbose=False)
        df = pd.DataFrame({"A": [1, 2]}, index=["x", "y"])

        cid1 = ana.update(df, intent="add_samples")
        self.assertEqual(cid1, 1)

        # force=True should create a new commit even with identical data
        cid2 = ana.update(df, force=True)
        self.assertGreater(cid2, cid1)

        meta = load_json(ana.commit_meta_path)
        self.assertEqual(len(meta["commits"]), 2)

    def test_update_different_data_creates_commit(self):
        """
        Test that update() with different data always creates a new commit.
        """
        ana = Ana(self.home, clear=True, verbose=False)
        df1 = pd.DataFrame({"A": [1, 2]}, index=["x", "y"])
        df2 = pd.DataFrame({"A": [3, 4]}, index=["z", "w"])

        cid1 = ana.update(df1, intent="add_samples")
        cid2 = ana.update(df2, intent="add_samples")
        self.assertGreater(cid2, cid1)

        meta = load_json(ana.commit_meta_path)
        self.assertEqual(len(meta["commits"]), 2)

    def test_max_commits_parameter(self):
        """
        Test that max_commits parameter is respected.
        """
        max_commits = 3
        ana = Ana(self.home, clear=True, max_commits=max_commits, verbose=False)
        
        # Create more commits than max_commits
        for i in range(max_commits + 2):
            ana.update(pd.DataFrame({"value": [i]}, index=[f"sample_{i}"]), intent="add_samples")
        
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

    def test_init_from_ana_different_directory(self):
        """
        Test initializing Ana from another Ana object with different directory.
        """
        # Create first Ana object with some data
        ana1 = Ana(self.home)
        data = pd.DataFrame({"A": [1, 2], "B": [3, 4]}, index=["x", "y"])
        ana1.update(data, key=None, intent="add_samples")
        ana1.update([1, 2, 3], key="list_data")
        
        # Create a second temporary directory for the new Ana object
        with TemporaryDirectory() as temp_dir2:
            # Initialize new Ana object from the first one
            ana2 = Ana(Path(temp_dir2), from_ana=ana1)
            
            # Check that data was copied correctly
            pd.testing.assert_frame_equal(ana1.data, ana2.data)
            self.assertEqual(ana1.obj, ana2.obj)
            
            # Check that they have different home directories
            self.assertNotEqual(ana1.home, ana2.home)
            
            # Check that the new Ana object has its own commit history
            self.assertEqual(len(ana2.commit_meta["commits"]), 1)
            self.assertEqual(ana2.commit_meta["commits"][0]["description"], 
                            "Initialize from another Ana object")

    def test_init_from_ana_same_directory_should_fail(self):
        """
        Test that initializing Ana from another Ana object with same directory raises ValueError.
        """
        # Create first Ana object
        ana1 = Ana(self.home)
        
        # Try to create second Ana object with same directory - should raise ValueError
        with self.assertRaises(ValueError) as context:
            Ana(self.home, from_ana=ana1)
        
        self.assertIn("must provide a different home directory", str(context.exception))

    # ===========================================
    # NEW TESTS FOR DATA TYPE PRESERVATION
    # ===========================================

    def test_pandas_categorical_dtype_preservation(self):
        """
        Test that categorical dtypes are preserved when saving and loading.
        """
        # Create DataFrame with categorical column
        df = self._sample_index(pd.DataFrame({
            'category': pd.Categorical(['a', 'b', 'c', 'a'], categories=['a', 'b', 'c', 'd'], ordered=True),
            'value': [1, 2, 3, 4]
        }))
        
        ana = Ana(self.home, data=df)
        
        # Check that dtypes are preserved
        loaded_df = ana.data
        self.assertIsInstance(loaded_df['category'].dtype, pd.CategoricalDtype)
        self.assertEqual(list(loaded_df['category'].cat.categories), ['a', 'b', 'c', 'd'])
        self.assertTrue(loaded_df['category'].cat.ordered)
        
        # Check the actual values
        pd.testing.assert_series_equal(df['category'], loaded_df['category'])

    def test_pandas_string_dtype_preservation(self):
        """
        Test that string dtypes are preserved when saving and loading.
        """
        # Create DataFrame with string dtype column (pandas 1.0+)
        df = self._sample_index(pd.DataFrame({
            'string_col': pd.Series(['hello', 'world', None, 'test'], dtype='string'),
            'value': [1, 2, 3, 4]
        }))
        
        ana = Ana(self.home, data=df)
        
        # Check that dtypes are preserved
        loaded_df = ana.data
        self.assertIsInstance(loaded_df['string_col'].dtype, pd.StringDtype)
        
        # Check the actual values
        pd.testing.assert_series_equal(df['string_col'], loaded_df['string_col'])

    def test_nullable_integer_dtype_preservation(self):
        """
        Test that nullable integer dtypes (Int64, etc.) are preserved.
        """
        # Create DataFrame with nullable integer column
        df = self._sample_index(pd.DataFrame({
            'nullable_int': pd.array([1, 2, None, 4], dtype='Int64'),
            'regular_int': [10, 20, 30, 40]
        }))
        
        ana = Ana(self.home, data=df)
        
        # Check that dtypes are preserved
        loaded_df = ana.data
        self.assertIsInstance(loaded_df['nullable_int'].dtype, pd.Int64Dtype)
        self.assertEqual(str(loaded_df['nullable_int'].dtype), 'Int64')
        
        # Check for nullable behavior (None should be preserved)
        self.assertTrue(pd.isna(loaded_df.loc["sample_2", 'nullable_int']))
        
        # Check the actual values
        pd.testing.assert_series_equal(df['nullable_int'], loaded_df['nullable_int'])

    def test_nullable_float_dtype_preservation(self):
        """
        Test that nullable float dtypes (Float64) are preserved.
        """
        # Create DataFrame with nullable float column
        df = self._sample_index(pd.DataFrame({
            'nullable_float': pd.array([1.1, None, 3.3, 4.4], dtype='Float64'),
            'regular_float': [1.1, 2.2, 3.3, 4.4]
        }))
        
        ana = Ana(self.home, data=df)
        
        # Check that dtypes are preserved
        loaded_df = ana.data
        self.assertIsInstance(loaded_df['nullable_float'].dtype, pd.Float64Dtype)
        self.assertEqual(str(loaded_df['nullable_float'].dtype), 'Float64')
        
        # Check for nullable behavior
        self.assertTrue(pd.isna(loaded_df.loc["sample_1", 'nullable_float']))
        
        # Check the actual values
        pd.testing.assert_series_equal(df['nullable_float'], loaded_df['nullable_float'])

    def test_datetime_dtype_preservation(self):
        """
        Test that datetime dtypes are preserved.
        """
        # Create DataFrame with datetime column
        df = self._sample_index(pd.DataFrame({
            'datetime': pd.to_datetime(['2023-01-01', '2023-01-02', '2023-01-03', '2023-01-04']),
            'value': [1, 2, 3, 4]
        }))
        
        ana = Ana(self.home, data=df)
        
        # Check that dtypes are preserved
        loaded_df = ana.data
        self.assertTrue(np.issubdtype(loaded_df['datetime'].dtype, np.datetime64))
            
            # Check the actual values
        pd.testing.assert_series_equal(df['datetime'], loaded_df['datetime'])

    def test_boolean_dtype_preservation(self):
        """
        Test that boolean dtypes (including nullable boolean) are preserved.
        """
        # Create DataFrame with nullable boolean column
        df = self._sample_index(pd.DataFrame({
            'nullable_bool': pd.array([True, False, None, True], dtype='boolean'),
            'regular_bool': [True, False, True, False]
        }))
        
        ana = Ana(self.home, data=df)
        
        # Check that dtypes are preserved
        loaded_df = ana.data
        self.assertIsInstance(loaded_df['nullable_bool'].dtype, pd.BooleanDtype)
        
        # Check for nullable behavior
        self.assertTrue(pd.isna(loaded_df.loc["sample_2", 'nullable_bool']))
        
        # Check the actual values
        pd.testing.assert_series_equal(df['nullable_bool'], loaded_df['nullable_bool'])

    def test_complex_dataframe_dtype_preservation(self):
        """
        Test that a DataFrame with multiple complex dtypes is fully preserved.
        """
        # Create DataFrame with various dtypes
        df = self._sample_index(pd.DataFrame({
            'category': pd.Categorical(['a', 'b', 'c', 'a']),
            'string': pd.Series(['hello', 'world', 'test', 'data'], dtype='string'),
            'nullable_int': pd.array([1, 2, None, 4], dtype='Int64'),
            'nullable_float': pd.array([1.1, None, 3.3, 4.4], dtype='Float64'),
            'datetime': pd.to_datetime(['2023-01-01', '2023-01-02', '2023-01-03', '2023-01-04']),
            'nullable_bool': pd.array([True, False, None, True], dtype='boolean'),
            'regular_int': [10, 20, 30, 40],
            'regular_float': [1.1, 2.2, 3.3, 4.4],
            'regular_bool': [True, False, True, False]
        }))
        
        ana = Ana(self.home, data=df)
        
        # Check that all dtypes are preserved
        loaded_df = ana.data
        
        # Verify each column's dtype
        self.assertIsInstance(loaded_df['category'].dtype, pd.CategoricalDtype)
        self.assertIsInstance(loaded_df['string'].dtype, pd.StringDtype)
        self.assertEqual(str(loaded_df['nullable_int'].dtype), 'Int64')
        self.assertEqual(str(loaded_df['nullable_float'].dtype), 'Float64')
        self.assertTrue(np.issubdtype(loaded_df['datetime'].dtype, np.datetime64))
        self.assertIsInstance(loaded_df['nullable_bool'].dtype, pd.BooleanDtype)
        
        # Check the entire DataFrame
        pd.testing.assert_frame_equal(df, loaded_df)

    def test_schema_file_creation(self):
        """
        Test that schema file is created when data is saved.
        """
        df = self._sample_index(pd.DataFrame({
            'category': pd.Categorical(['a', 'b', 'c']),
            'string': pd.Series(['x', 'y', 'z'], dtype='string'),
            'nullable_int': pd.array([1, 2, 3], dtype='Int64')
        }))
        
        ana = Ana(self.home, data=df)
        
        # Check that schema file exists
        self.assertTrue(ana.schema_path.exists())
        
        # Load and inspect schema
        schema = load_json(ana.schema_path)
        
        # Check schema structure
        self.assertIn('columns', schema)
        self.assertIn('category', schema['columns'])
        self.assertIn('string', schema['columns'])
        self.assertIn('nullable_int', schema['columns'])
        
        # Check schema content for categorical
        self.assertEqual(schema['columns']['category']['type'], 'categorical')
        self.assertEqual(schema['columns']['category']['categories'], ['a', 'b', 'c'])
        self.assertEqual(schema['columns']['category']['ordered'], False)

    def test_dtype_preservation_after_update(self):
        """
        Test that dtypes are preserved after updating data.
        """
        # Initial data
        df1 = pd.DataFrame({
            'category': pd.Categorical(['a', 'b'], categories=['a', 'b', 'c']),
            'nullable_int': pd.array([1, 2], dtype='Int64')
        }, index=['x', 'y'])
        
        ana = Ana(self.home, data=df1)
        
        # Update with new data
        df2 = pd.DataFrame({
            'category': pd.Categorical(['c', 'a'], categories=['a', 'b', 'c']),
            'nullable_int': pd.array([None, 4], dtype='Int64')
        }, index=['z', 'w'])
        
        ana.update(df2, intent="add_samples")
        
        # Check that dtypes are still preserved after update
        loaded_df = ana.data
        
        self.assertIsInstance(loaded_df['category'].dtype, pd.CategoricalDtype)
        self.assertEqual(str(loaded_df['nullable_int'].dtype), 'Int64')
        
        # Check that categorical categories are preserved
        self.assertEqual(list(loaded_df['category'].cat.categories), ['a', 'b', 'c'])

    def test_dtype_preservation_after_revert(self):
        """
        Test that dtypes are preserved after revert operation.
        """
        ana = Ana(self.home, clear=True, max_commits=10, verbose=False)
        
        # First commit with complex dtypes
        df1 = pd.DataFrame({
            'category': pd.Categorical(['a', 'b'], categories=['a', 'b', 'c']),
            'string': pd.Series(['x', 'y'], dtype='string')
        }, index=['i1', 'i2'])
        
        ana.update(df1, description="First state with dtypes", intent="add_samples")
        
        # Second commit with different data
        df2 = pd.DataFrame({
            'regular': [1, 2, 3],
            'other': [4, 5, 6]
        }, index=['j1', 'j2', 'j3'])
        
        ana.update(df2, description="Second state", intent="add_samples")
        
        # Revert to first commit
        ana.revert(1)
        
        # Check that dtypes from first commit are restored
        loaded_df = ana.data
        
        self.assertIsInstance(loaded_df['category'].dtype, pd.CategoricalDtype)
        self.assertIsInstance(loaded_df['string'].dtype, pd.StringDtype)
        
        # Check categorical categories
        self.assertEqual(list(loaded_df['category'].cat.categories), ['a', 'b', 'c'])

    def test_get_schema_method(self):
        """
        Test the get_schema method returns correct schema information.
        """
        df = self._sample_index(pd.DataFrame({
            'category': pd.Categorical(['a', 'b', 'c']),
            'nullable_int': pd.array([1, 2, 3], dtype='Int64'),
            'regular': [1.1, 2.2, 3.3]
        }))
        
        ana = Ana(self.home, data=df)
        
        # Get schema using method
        schema = ana.get_schema()
        
        # Check schema structure
        self.assertIn('columns', schema)
        self.assertIn('category', schema['columns'])
        self.assertIn('nullable_int', schema['columns'])
        self.assertIn('regular', schema['columns'])
        
        # Check schema content
        self.assertEqual(schema['columns']['category']['type'], 'categorical')
        self.assertEqual(schema['columns']['nullable_int']['type'], 'nullable_int')
        self.assertEqual(schema['columns']['nullable_int']['bits'], '64')
        self.assertEqual(schema['columns']['regular']['type'], 'standard')

    def test_empty_dataframe_dtype_handling(self):
        """
        Test that empty DataFrames are handled correctly.
        """
        # Create empty DataFrame with specific dtypes
        df = pd.DataFrame({
            'category': pd.Series([], dtype='category'),
            'string': pd.Series([], dtype='string'),
            'nullable_int': pd.Series([], dtype='Int64')
        })
        
        ana = Ana(self.home, data=df)
        
        # Check that empty DataFrame loads correctly
        loaded_df = ana.data
        
        self.assertEqual(len(loaded_df), 0)
        self.assertEqual(list(loaded_df.columns), ['category', 'string', 'nullable_int'])
        
        # Check dtypes (may be object for empty series, but schema should be saved)
        if len(loaded_df) > 0:
            self.assertIsInstance(loaded_df['category'].dtype, pd.CategoricalDtype)
            self.assertIsInstance(loaded_df['string'].dtype, pd.StringDtype)
            self.assertEqual(str(loaded_df['nullable_int'].dtype), 'Int64')

    def test_sample_index_canonical_string_dtype(self):
        """
        Test that the sample axis is stored as canonical string dtype.
        """
        df = pd.DataFrame(
            {'value': [1, 2, 3, 4]},
            index=pd.CategoricalIndex(['a', 'b', 'c', 'd'], name='category_idx')
        )
        
        ana = Ana(self.home, data=df)
        
        loaded_df = ana.data

        self.assertIsInstance(loaded_df.index.dtype, pd.StringDtype)
        self.assertEqual(loaded_df.index.name, 'category_idx')
        self.assertEqual(list(loaded_df.index), ['a', 'b', 'c', 'd'])

    def test_backward_compatibility_without_schema(self):
        """
        Test that data can still be loaded when schema file is missing (backward compatibility).
        """
        # Create a CSV file without schema
        df = self._sample_index(pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]}))
        csv_path = self.home / "data.csv"
        df.to_csv(csv_path)
        
        # Delete schema file if it exists
        schema_path = self.home / "data_schema.json"
        if schema_path.exists():
            schema_path.unlink()
        
        # Create Ana instance - should load without schema
        ana = Ana(self.home, clear=False)
        
        # Should still load data (though dtypes may not be preserved)
        loaded_df = ana.data
        self.assertEqual(len(loaded_df), 3)
        self.assertEqual(list(loaded_df.columns), ['A', 'B'])

    def test_schema_file_corruption_handling(self):
        """
        Test that corrupted schema file doesn't break data loading.
        """
        # Create valid data
        df = self._sample_index(pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]}))
        ana = Ana(self.home, data=df)
        
        # Corrupt the schema file
        with open(ana.schema_path, 'w') as f:
            f.write('{invalid json')
        
        # Try to reload - should handle gracefully
        ana2 = Ana(self.home, clear=False)
        
        # Should still load data (with warning)
        loaded_df = ana2.data
        self.assertEqual(len(loaded_df), 3)

    def test_file_column_auto_detection(self):
        """
        Test that columns containing file paths are auto-detected as type 'file'.
        """
        df = pd.DataFrame({
            "cool_paths": ["/data/sample1.cool", "/data/sample2.cool"],
            "bam_paths": ["/align/s1.bam", "/align/s2.bam"],
            "values": [1.0, 2.0],
        }, index=["s1", "s2"])

        ana = Ana(self.home, data=df)
        schema = ana.get_schema()

        self.assertEqual(schema["columns"]["cool_paths"]["type"], "file")
        self.assertEqual(schema["columns"]["cool_paths"]["file_format"], "cool")
        self.assertEqual(schema["columns"]["bam_paths"]["type"], "file")
        self.assertEqual(schema["columns"]["bam_paths"]["file_format"], "bam")
        self.assertEqual(schema["columns"]["values"]["type"], "standard")

    def test_column_metadata_injection(self):
        """
        Test that column_descriptions and depends_on are stored in column_meta.json.
        """
        df = pd.DataFrame({"A": [1, 2], "B": [3, 4]}, index=["x", "y"])
        ana = Ana(self.home, data=df)

        # Second update with metadata
        df2 = pd.DataFrame({"C": [5, 6]}, index=["x", "y"])
        ana.update(
            df2,
            column_descriptions={"C": "Derived column"},
            depends_on={"C": {"columns": ["A", "B"]}},
            computed_by={"file": "/scripts/calc.ipynb"},
        )

        meta = ana.column_meta("C")
        self.assertEqual(meta["description"], "Derived column")
        self.assertEqual(meta["computed_by"]["file"], "/scripts/calc.ipynb")
        self.assertEqual(meta["depends_on"]["columns"], ["A", "B"])

    def test_describe_output(self):
        """
        Test that describe() prints column and object info.
        """
        import io
        import sys

        df = pd.DataFrame({"A": [1, 2]}, index=["x", "y"])
        ana = Ana(self.home, data=df)
        ana.update([1, 2, 3], key="my_list")

        # Capture stdout
        old_stdout = sys.stdout
        sys.stdout = captured = io.StringIO()
        try:
            ana.describe()
        finally:
            sys.stdout = old_stdout

        output = captured.getvalue()
        self.assertIn("Columns (1):", output)
        self.assertIn("A", output)
        self.assertIn("Objects (1):", output)
        self.assertIn("my_list", output)

    def test_file_format_lookup(self):
        """
        Test that file_format() returns built-in definitions.
        """
        df = pd.DataFrame({
            "cool_paths": ["/data/a.cool", "/data/b.cool"],
        }, index=["a", "b"])
        ana = Ana(self.home, data=df)

        fmt = ana.file_format("cool_paths")
        self.assertIsNotNone(fmt)
        self.assertEqual(fmt["format"], "Cooler")

        # Direct key lookup
        fmt2 = ana.file_format("bam")
        self.assertEqual(fmt2["format"], "BAM")

        # Unknown key
        self.assertIsNone(ana.file_format("nonexistent_format"))

    def test_reject_rangeindex_dataframe_as_named_result_table(self):
        """
        DataFrames with non-sample RangeIndex must not be merged into ana.data.
        """
        ana = Ana(
            self.home,
            data=pd.DataFrame({"cell_type": ["2c"]}, index=["sample_1"]),
        )
        df = pd.DataFrame({"gene_id": ["g1"], "score": [1.0]})

        with self.assertRaises(AnaUpdateError) as ctx:
            ana.update(df, key="top_genes")

        self.assertIn("refusing to merge DataFrame into ana.data", str(ctx.exception))
        self.assertIn("analysis result table", str(ctx.exception))
        self.assertEqual(list(ana.data.index), ["sample_1"])
        meta = load_json(ana.commit_meta_path)
        self.assertEqual(len(meta["commits"]), 1)

    def test_manifest_dict_auto_routes_to_object(self):
        """
        Dict keys that do not match sample IDs are treated as object manifests.
        """
        ana = Ana(
            self.home,
            data=pd.DataFrame({"cell_type": ["2c"]}, index=["sample_1"]),
        )

        cid = ana.update(
            {"top_residuals_tsv": "/tmp/top_residuals.tsv"},
            key="analysis_manifest",
            description="Register analysis manifest",
        )

        self.assertEqual(ana.obj["analysis_manifest"]["top_residuals_tsv"], "/tmp/top_residuals.tsv")
        self.assertEqual(list(ana.data.index), ["sample_1"])
        self.assertEqual(cid, 2)

    def test_dict_auto_updates_existing_sample_column(self):
        """
        Dict keys matching existing sample IDs update a sample-level column.
        """
        ana = Ana(
            self.home,
            data=pd.DataFrame({"cell_type": ["2c"]}, index=["sample_1"]),
        )

        ana.update({"sample_1": True}, key="qc_pass")

        self.assertTrue(ana.data.loc["sample_1", "qc_pass"])
        self.assertNotIn("qc_pass", ana.obj)

    def test_dict_auto_refuses_accidental_new_sample(self):
        """
        Sample-looking dict keys must not add samples in auto mode.
        """
        ana = Ana(
            self.home,
            data=pd.DataFrame({"cell_type": ["2c"]}, index=["sample_1"]),
        )

        with self.assertRaises(AnaUpdateError) as ctx:
            ana.update({"sample_2": True}, key="qc_pass")

        self.assertIn("refusing to add new sample IDs", str(ctx.exception))
        self.assertNotIn("sample_2", ana.data.index)
        self.assertNotIn("qc_pass", ana.obj)

    def test_auto_dataframe_updates_existing_samples_only(self):
        """
        Auto DataFrame routing can update columns for existing sample IDs.
        """
        ana = Ana(
            self.home,
            data=pd.DataFrame({"cell_type": ["2c", "4c"]}, index=["sample_1", "sample_2"]),
        )

        ana.update(pd.DataFrame({"qc_pass": [True]}, index=["sample_1"]))

        self.assertTrue(ana.data.loc["sample_1", "qc_pass"])
        self.assertTrue(pd.isna(ana.data.loc["sample_2", "qc_pass"]))

    def test_auto_dataframe_partial_overlap_is_rejected(self):
        """
        Partial sample-axis overlap is ambiguous and must be explicit.
        """
        ana = Ana(
            self.home,
            data=pd.DataFrame({"cell_type": ["2c", "4c"]}, index=["sample_1", "sample_2"]),
        )
        update = pd.DataFrame({"qc_pass": [True, False]}, index=["sample_1", "sample_3"])

        with self.assertRaises(AnaUpdateError) as ctx:
            ana.update(update)

        self.assertIn("refusing to merge DataFrame", str(ctx.exception))
        self.assertNotIn("sample_3", list(ana.data.index))

    def test_sample_axis_operations(self):
        """
        Explicit add/drop/rename/replace intents maintain the sample axis.
        """
        ana = Ana(
            self.home,
            data=pd.DataFrame({"cell_type": ["2c", "4c"]}, index=["sample_1", "sample_2"]),
        )
        ana.update(["sample_1", "sample_2"], key="all_samples")

        ana.update(
            pd.DataFrame({"sample_id": ["sample_3"], "cell_type": ["8c"]}),
            intent="add_samples",
            description="Add sample_3",
        )
        self.assertIn("sample_3", ana.data.index)
        ana.update(["sample_1", "sample_2", "sample_3"], key="all_samples")

        ana.update({"sample_3": "sample_03"}, intent="rename_samples", description="Rename sample_3")
        self.assertIn("sample_03", ana.data.index)
        self.assertIn("sample_03", ana.obj["all_samples"])

        ana.update(["sample_2"], intent="drop_samples", description="Drop sample_2")
        self.assertNotIn("sample_2", ana.data.index)
        self.assertNotIn("sample_2", ana.obj["all_samples"])

        replacement = pd.DataFrame({"cell_type": ["2c"]}, index=["sample_1"])
        ana.update(replacement, intent="replace_samples", description="Rebuild sample axis")
        self.assertEqual(list(ana.data.index), ["sample_1"])

        meta = load_json(ana.commit_meta_path)
        last_plan = meta["commits"][-1]["update_plan"]
        self.assertEqual(last_plan["target"], "replace_samples")
        self.assertIn("sample_03", last_plan["axis_delta"]["dropped"])

    def test_sample_axis_canonical_collision_is_rejected(self):
        """
        Labels that collide after CSV/string canonicalization are unsafe.
        """
        ana = Ana(self.home, clear=True)
        df = pd.DataFrame({"A": [1, 2]}, index=pd.Index([0, "0"], dtype=object))

        with self.assertRaises(AnaUpdateError) as ctx:
            ana.update(df, intent="add_samples")

        self.assertIn("collide after CSV/string canonicalization", str(ctx.exception))

    def test_add_samples_rejects_existing_canonical_collision_before_write(self):
        """
        Existing "0" and incoming 0 collide after CSV persistence.
        """
        ana = Ana(self.home, data=pd.DataFrame({"x": [1]}, index=["0"]))

        with self.assertRaises(AnaUpdateError) as ctx:
            ana.update(pd.DataFrame({"x": [2]}, index=pd.Index([0], dtype=object)), intent="add_samples")

        self.assertIn("collides with an existing ID", str(ctx.exception))
        expected = pd.DataFrame({"x": [1]}, index=pd.Index(["0"], dtype=pd.StringDtype()))
        pd.testing.assert_frame_equal(ana.data, expected)

    def test_commit_metadata_includes_ana_system_version(self):
        """
        Ana version metadata is persisted globally and per commit.
        """
        ana = Ana(self.home, clear=True)
        ana.update(pd.DataFrame({"A": [1]}, index=["sample_1"]), intent="add_samples")

        meta = load_json(ana.commit_meta_path)
        self.assertEqual(meta["ana_system_version"], 3)
        self.assertEqual(meta["axis"]["name"], "sample_id")
        self.assertEqual(meta["commits"][-1]["ana_system_version"], 3)
        self.assertEqual(meta["commits"][-1]["axis"]["role"], "sample_axis")
        self.assertEqual(meta["commits"][-1]["update_plan"]["target"], "add_samples")
        self.assertEqual(meta["commits"][-1]["update_plan"]["axis_delta"]["added"], ["sample_1"])

    def test_add_samples_csv_roundtrip_is_semantic_noop(self):
        """An independently loaded CSV must normalize to the locked schema."""
        base = pd.DataFrame(
            {
                "batch": pd.array(["old"], dtype="string"),
                "raw_reads": pd.array([10], dtype="Int64"),
                "score": pd.array([1.5], dtype="Float64"),
                "path": pd.array(["/old/a.cool"], dtype="string"),
                "flag": pd.array([True], dtype="boolean"),
                "all_missing": pd.array([pd.NA], dtype="Float64"),
            },
            index=["old_1"],
        )
        incoming = pd.DataFrame(
            {
                "batch": pd.array(["new", "new", "new"], dtype="string"),
                "raw_reads": pd.array([100, 0, pd.NA], dtype="Int64"),
                "score": pd.array([3.25, pd.NA, 5.75], dtype="Float64"),
                "path": pd.array(
                    ["/new/a.cool", "/new/b.cool", pd.NA], dtype="string"
                ),
                "flag": pd.array([True, pd.NA, False], dtype="boolean"),
                "all_missing": pd.array([pd.NA, pd.NA, pd.NA], dtype="Float64"),
            },
            index=["new_1", "new_2", "new_3"],
        )
        ana = Ana(self.home, data=base, verbose=False)
        description = "Add heterogeneous batch metadata"
        commit_id = ana.update(
            incoming,
            intent="add_samples",
            description=description,
            computed_by={"file": "batch.ipynb"},
        )

        external_csv = self.home / "batch.csv"
        incoming.to_csv(external_csv)
        replay = pd.read_csv(external_csv, index_col=0)
        before = self._home_file_bytes()

        replay_commit_id = ana.update(
            replay,
            intent="add_samples",
            description=description,
            computed_by={"file": "batch.ipynb"},
        )

        self.assertEqual(replay_commit_id, commit_id)
        self.assertEqual(self._home_file_bytes(), before)
        self.assertEqual(str(ana.data["raw_reads"].dtype), "Int64")
        self.assertEqual(str(ana.data["flag"].dtype), "boolean")
        self.assertEqual(len(load_json(ana.commit_meta_path)["commits"]), 2)

    def test_implicit_schema_change_is_rejected_without_writes(self):
        """Existing column schemas cannot change through dtype inference."""
        ana = Ana(
            self.home,
            data=pd.DataFrame(
                {"value": pd.array([1], dtype="Int64")},
                index=["sample_1"],
            ),
            verbose=False,
        )
        before = self._home_file_bytes()

        with self.assertRaises(AnaUpdateError) as ctx:
            ana.update(
                pd.DataFrame({"value": [1.5]}, index=["sample_1"]),
                description="Introduce fractional values",
            )

        self.assertIn("does not evolve existing schemas automatically", str(ctx.exception))
        self.assertIn("schema_changes", str(ctx.exception))
        self.assertEqual(self._home_file_bytes(), before)

    def test_explicit_schema_change_is_committed_and_auditable(self):
        """Explicit dtype migration changes the schema and records its delta."""
        ana = Ana(
            self.home,
            data=pd.DataFrame(
                {"value": pd.array([1], dtype="Int64")},
                index=["sample_1"],
            ),
            verbose=False,
        )

        commit_id = ana.update(
            pd.DataFrame({"value": [1.5]}, index=["sample_1"]),
            description="Allow fractional values",
            schema_changes={"value": "Float64"},
        )

        self.assertEqual(commit_id, 2)
        self.assertEqual(str(ana.data["value"].dtype), "Float64")
        commit = load_json(ana.commit_meta_path)["commits"][-1]
        delta = commit["update_plan"]["schema_delta"]["value"]
        self.assertEqual(delta["from"]["type"], "nullable_int")
        self.assertEqual(delta["to"]["type"], "nullable_float")

    def test_category_expansion_requires_explicit_schema_change(self):
        """New category labels cannot silently alter an existing domain."""
        dtype = pd.CategoricalDtype(["2c", "4c"], ordered=True)
        ana = Ana(
            self.home,
            data=pd.DataFrame(
                {"stage": pd.Series(["2c"], dtype=dtype, index=["sample_1"])},
                index=["sample_1"],
            ),
            verbose=False,
        )
        update = pd.DataFrame({"stage": ["8c"]}, index=["sample_2"])
        before = self._home_file_bytes()

        with self.assertRaises(AnaUpdateError):
            ana.update(update, intent="add_samples", description="Add 8c sample")
        self.assertEqual(self._home_file_bytes(), before)

        expanded = pd.CategoricalDtype(["2c", "4c", "8c"], ordered=True)
        ana.update(
            update,
            intent="add_samples",
            description="Add 8c category and sample",
            schema_changes={"stage": expanded},
        )
        self.assertEqual(list(ana.data["stage"].cat.categories), ["2c", "4c", "8c"])

    def test_commit_snapshot_hash_and_revert_are_byte_exact(self):
        """Recorded hashes identify exact snapshots, including float CSV text."""
        ana = Ana(
            self.home,
            data=pd.DataFrame(
                {"value": [-1.8507079234797329e-75]},
                index=["sample_1"],
            ),
            verbose=False,
        )
        first = load_json(ana.commit_meta_path)["commits"][-1]
        stem = first["filename"]
        snapshot_paths = [
            self.home / f"{stem}.data.csv",
            self.home / f"{stem}.data_schema.json",
            self.home / f"{stem}.obj.json",
        ]
        snapshot_bytes = [path.read_bytes() for path in snapshot_paths]

        self.assertEqual(first["state_hash"], ana._compute_state_hash())
        self.assertEqual(first["state_hash"], self._combined_hash(snapshot_paths))
        self.assertEqual(ana.data_path.read_bytes(), snapshot_paths[0].read_bytes())
        self.assertEqual(ana.update(ana.data.copy()), first["id"])

        ana.update(
            pd.DataFrame({"value": [-2.0]}, index=["sample_1"]),
            description="Change float value",
        )
        ana.revert(first["id"])
        current_bytes = [
            ana.data_path.read_bytes(),
            ana.schema_path.read_bytes(),
            ana.obj_path.read_bytes(),
        ]
        self.assertEqual(current_bytes, snapshot_bytes)

    def test_rotated_commit_snapshots_preserve_hash_invariant(self):
        """Circular snapshot slots retain exact bytes for every live commit."""
        ana = Ana(self.home, clear=True, max_commits=3, verbose=False)
        for i in range(1, 7):
            ana.update(
                pd.DataFrame({"value": [i]}, index=[f"sample_{i}"]),
                intent="add_samples",
                description=f"Add sample {i}",
            )

        commits = load_json(ana.commit_meta_path)["commits"]
        self.assertEqual([commit["id"] for commit in commits], [4, 5, 6])
        for commit in commits:
            stem = commit["filename"]
            paths = [
                self.home / f"{stem}.data.csv",
                self.home / f"{stem}.data_schema.json",
                self.home / f"{stem}.obj.json",
            ]
            self.assertEqual(commit["state_hash"], self._combined_hash(paths))

if __name__ == "__main__":
    unittest.main()
