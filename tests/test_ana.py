import unittest
import os
import pandas as pd
import numpy as np
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
        df = pd.DataFrame({"A": [1, 2], "B": [3, 4]})
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
        df = pd.DataFrame({"A": [1, 2], "B": [3, 4]})
        ana = Ana(self.home, data=df)
        # With dtype preservation, the index type should be preserved
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
        # With dtype preservation, the index type should be preserved
        pd.testing.assert_frame_equal(ana.data, expected)

    def test_update_series(self):
        """
        Test update method with Series input.
        """
        ser = pd.Series([1, 2], index=["x", "y"], name="ser")
        ana = Ana(self.home)
        ana.update(ser, key="ser")
        expected = ser.to_frame()
        # With dtype preservation, the index type should be preserved
        pd.testing.assert_frame_equal(ana.data, expected)

    def test_update_dict(self):
        """
        Test update method with dict input.
        """
        d = {"a": 1, "b": 2}
        ana = Ana(self.home, clear=True)
        ana.update(d, key="dict_key")
        expected = pd.DataFrame({"dict_key": pd.Series(d)})
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
        ana.update(df)
        
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
        # `force`` only work for dict and series input 
        dat = {
            "z" : 3,
            "w" : 4
        }
        ana.update(dat, key="A", force=False)
        self.assertEqual(ana.data.shape, initial_shape)  # No change
        
        # Update with force - should overwrite
        ana.update(dat, key="A", force=True)
        self.assertEqual(ana.data.shape[0], 4)  # Add new data

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

    def test_init_from_ana_different_directory(self):
        """
        Test initializing Ana from another Ana object with different directory.
        """
        # Create first Ana object with some data
        ana1 = Ana(self.home)
        data = pd.DataFrame({"A": [1, 2], "B": [3, 4]}, index=["x", "y"])
        ana1.update(data, key=None)
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
        df = pd.DataFrame({
            'category': pd.Categorical(['a', 'b', 'c', 'a'], categories=['a', 'b', 'c', 'd'], ordered=True),
            'value': [1, 2, 3, 4]
        })
        
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
        df = pd.DataFrame({
            'string_col': pd.Series(['hello', 'world', None, 'test'], dtype='string'),
            'value': [1, 2, 3, 4]
        })
        
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
        df = pd.DataFrame({
            'nullable_int': pd.array([1, 2, None, 4], dtype='Int64'),
            'regular_int': [10, 20, 30, 40]
        })
        
        ana = Ana(self.home, data=df)
        
        # Check that dtypes are preserved
        loaded_df = ana.data
        self.assertIsInstance(loaded_df['nullable_int'].dtype, pd.Int64Dtype)
        self.assertEqual(str(loaded_df['nullable_int'].dtype), 'Int64')
        
        # Check for nullable behavior (None should be preserved)
        self.assertTrue(pd.isna(loaded_df.loc[2, 'nullable_int']))
        
        # Check the actual values
        pd.testing.assert_series_equal(df['nullable_int'], loaded_df['nullable_int'])

    def test_nullable_float_dtype_preservation(self):
        """
        Test that nullable float dtypes (Float64) are preserved.
        """
        # Create DataFrame with nullable float column
        df = pd.DataFrame({
            'nullable_float': pd.array([1.1, None, 3.3, 4.4], dtype='Float64'),
            'regular_float': [1.1, 2.2, 3.3, 4.4]
        })
        
        ana = Ana(self.home, data=df)
        
        # Check that dtypes are preserved
        loaded_df = ana.data
        self.assertIsInstance(loaded_df['nullable_float'].dtype, pd.Float64Dtype)
        self.assertEqual(str(loaded_df['nullable_float'].dtype), 'Float64')
        
        # Check for nullable behavior
        self.assertTrue(pd.isna(loaded_df.loc[1, 'nullable_float']))
        
        # Check the actual values
        pd.testing.assert_series_equal(df['nullable_float'], loaded_df['nullable_float'])

    def test_datetime_dtype_preservation(self):
        """
        Test that datetime dtypes are preserved.
        """
        # Create DataFrame with datetime column
        df = pd.DataFrame({
            'datetime': pd.to_datetime(['2023-01-01', '2023-01-02', '2023-01-03', '2023-01-04']),
            'value': [1, 2, 3, 4]
        })
        
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
        df = pd.DataFrame({
            'nullable_bool': pd.array([True, False, None, True], dtype='boolean'),
            'regular_bool': [True, False, True, False]
        })
        
        ana = Ana(self.home, data=df)
        
        # Check that dtypes are preserved
        loaded_df = ana.data
        self.assertIsInstance(loaded_df['nullable_bool'].dtype, pd.BooleanDtype)
        
        # Check for nullable behavior
        self.assertTrue(pd.isna(loaded_df.loc[2, 'nullable_bool']))
        
        # Check the actual values
        pd.testing.assert_series_equal(df['nullable_bool'], loaded_df['nullable_bool'])

    def test_complex_dataframe_dtype_preservation(self):
        """
        Test that a DataFrame with multiple complex dtypes is fully preserved.
        """
        # Create DataFrame with various dtypes
        df = pd.DataFrame({
            'category': pd.Categorical(['a', 'b', 'c', 'a']),
            'string': pd.Series(['hello', 'world', 'test', 'data'], dtype='string'),
            'nullable_int': pd.array([1, 2, None, 4], dtype='Int64'),
            'nullable_float': pd.array([1.1, None, 3.3, 4.4], dtype='Float64'),
            'datetime': pd.to_datetime(['2023-01-01', '2023-01-02', '2023-01-03', '2023-01-04']),
            'nullable_bool': pd.array([True, False, None, True], dtype='boolean'),
            'regular_int': [10, 20, 30, 40],
            'regular_float': [1.1, 2.2, 3.3, 4.4],
            'regular_bool': [True, False, True, False]
        })
        
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
        df = pd.DataFrame({
            'category': pd.Categorical(['a', 'b', 'c']),
            'string': pd.Series(['x', 'y', 'z'], dtype='string'),
            'nullable_int': pd.array([1, 2, 3], dtype='Int64')
        })
        
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
        
        ana.update(df2)
        
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
        
        ana.update(df1, description="First state with dtypes")
        
        # Second commit with different data
        df2 = pd.DataFrame({
            'regular': [1, 2, 3],
            'other': [4, 5, 6]
        }, index=['j1', 'j2', 'j3'])
        
        ana.update(df2, description="Second state")
        
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
        df = pd.DataFrame({
            'category': pd.Categorical(['a', 'b', 'c']),
            'nullable_int': pd.array([1, 2, 3], dtype='Int64'),
            'regular': [1.1, 2.2, 3.3]
        })
        
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

    def test_index_dtype_preservation(self):
        """
        Test that index dtypes are preserved when saving and loading.
        """
        # Create DataFrame with categorical index
        df = pd.DataFrame(
            {'value': [1, 2, 3, 4]},
            index=pd.CategoricalIndex(['a', 'b', 'c', 'd'], name='category_idx')
        )
        
        ana = Ana(self.home, data=df)
        
        # Check that index dtype is preserved
        loaded_df = ana.data
        
        self.assertTrue(isinstance(loaded_df.index, pd.CategoricalIndex))
        self.assertEqual(loaded_df.index.name, 'category_idx')
        self.assertEqual(list(loaded_df.index), ['a', 'b', 'c', 'd'])

    def test_backward_compatibility_without_schema(self):
        """
        Test that data can still be loaded when schema file is missing (backward compatibility).
        """
        # Create a CSV file without schema
        df = pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]})
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
        df = pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]})
        ana = Ana(self.home, data=df)
        
        # Corrupt the schema file
        with open(ana.schema_path, 'w') as f:
            f.write('{invalid json')
        
        # Try to reload - should handle gracefully
        ana2 = Ana(self.home, clear=False)
        
        # Should still load data (with warning)
        loaded_df = ana2.data
        self.assertEqual(len(loaded_df), 3)

if __name__ == "__main__":
    unittest.main()