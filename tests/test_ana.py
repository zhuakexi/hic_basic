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
        self.tag = "test_tag"

    def tearDown(self):
        """
        Clean up the temporary directory after each test.
        """
        self.temp_dir.cleanup()

    def test_init_without_clear(self):
        """
        Test initialization with clear=False retains existing files.
        """
        data_path = self.home / f"{self.tag}.data.csv.gz"
        obj_path = self.home / f"{self.tag}.obj.json"
        for path in [data_path, obj_path]:
            path.touch()

        Ana(self.home, tag=self.tag, clear=False)
        self.assertTrue(data_path.exists())
        self.assertTrue(obj_path.exists())

    def test_init_with_data(self):
        """
        Test initialization with data parameter.
        """
        df = pd.DataFrame({"A": [1, 2], "B": [3, 4]})
        ana = Ana(self.home, data=df, tag=self.tag)
        df.index = df.index.astype("string")
        pd.testing.assert_frame_equal(ana.data, df)

    def test_init_with_obj(self):
        """
        Test initialization with obj parameter (dict of key-value lists).
        """
        obj = {"key1": [1, 2], "key2": [3, 4]}
        ana = Ana(self.home, obj=obj, tag=self.tag)
        self.assertEqual(ana.obj, obj)

    def test_data_property(self):
        """
        Test data property correctly reads and writes DataFrame.
        """
        df = pd.DataFrame({"A": [1, 2], "B": [3, 4]})
        ana = Ana(self.home, data=df, tag=self.tag)
        df.index = df.index.astype("string")
        pd.testing.assert_frame_equal(ana.data, df)

    def test_update_dataframe(self):
        """
        Test update method with DataFrame input.
        """
        df1 = pd.DataFrame({"A": [1, 2]}, index=["x", "y"])
        df2 = pd.DataFrame({"A": [3], "B": [4]}, index=["z"])
        ana = Ana(self.home, data=df1, tag=self.tag)
        ana.update(df2)
        expected = pd.concat([df1, df2], axis=0).groupby(level=0).last()
        expected.index = expected.index.astype("string")
        pd.testing.assert_frame_equal(ana.data, expected)

    def test_update_series(self):
        """
        Test update method with Series input.
        """
        ser = pd.Series([1, 2], index=["x", "y"], name="ser")
        ana = Ana(self.home, tag=self.tag)
        ana.update(ser, key="ser")
        expected = ser.to_frame().T
        expected.index = expected.index.astype("string")
        pd.testing.assert_frame_equal(ana.data, expected)

    def test_update_dict(self):
        """
        Test update method with dict input.
        """
        d = {"a": 1, "b": 2}
        ana = Ana(self.home, tag=self.tag, clear=True)
        ana.update(d, key="dict_key")
        expected = pd.DataFrame({"dict_key": pd.Series(d)})
        expected.index = expected.index.astype("string")
        pd.testing.assert_frame_equal(ana.data, expected)

    def test_update_list(self):
        """
        Test update method with list input.
        """
        lst = [1, 2, 3]
        ana = Ana(self.home, tag=self.tag)
        ana.update(lst, key="list_key")
        self.assertEqual(ana.obj["list_key"], lst)

    def test_unsupported_type(self):
        """
        Test update method raises TypeError for unsupported data types.
        """
        ana = Ana(self.home, tag=self.tag)
        with self.assertRaises(TypeError):
            ana.update(object())

    def test_required_key_for_dict_and_list(self):
        """
        Test dict and list inputs require key to be provided.
        """
        ana = Ana(self.home, tag=self.tag)
        with self.assertRaises(AssertionError):
            ana.update([1, 2])  # Missing key
        with self.assertRaises(AssertionError):
            ana.update({"a": 1})  # Missing key
if __name__ == "__main__":
    unittest.main()