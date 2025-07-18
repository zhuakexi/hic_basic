# This module is designed for bioinformatic project management
from pathlib import Path

import pandas as pd
from .hicio import load_json, dump_json
class Ana:
    """
    A class to manage analysis data and operations.
    Manage intermideate file paths by careful operation of DataFrame storage.
    Store json data with tags.
    The Ana class has no inner data object, all data states are stored in files.
    """
    def __init__(self, home, tag=None, data=None, obj=None, clear=False):
        """
        Initialize the Ana class.
        Input:
            home: str, path to the home directory where data files are stored.
            data: pd.DataFrame or dict or list, initial data to be stored.
            tag: str, optional, a tag for the data file.
            clear: bool, optional, if True, clear existing data file.
        """
        self.home = Path(home)
        self.home.mkdir(parents=True, exist_ok=True)
        # initialize the tag
        if tag is None:
            self.tag = "0"
        else:
            self.tag = tag
        # initialize data files
        data_path = self.home / f"{self.tag}.data.csv.gz"
        obj_path = self.home / f"{self.tag}.obj.json"
        if data_path.exists():
            if clear:
                data_path.unlink()
                data_path.touch()
        else:
            data_path.touch()
        if obj_path.exists():
            if clear:
                obj_path.unlink()
                Ana.create_empty_json_file(obj_path)
        else:
            Ana.create_empty_json_file(obj_path)
        # initialize the data
        if not data is None:
            self.update(data, key=None)
        # initialize the object
        if not obj is None:
            for key, value in obj.items():
                self.update(value, key=key)
    @staticmethod
    def create_empty_json_file(path):
        """
        Create an empty JSON file at the specified path.
        """
        dump_json({}, path)
    @property
    def data_path(self):
        """
        Get the current data file path.
        """
        return self.home / f"{self.tag}.data.csv.gz"
    @property
    def obj_path(self):
        """
        Get the current object file path.
        """
        return self.home / f"{self.tag}.obj.json"
    @property
    def data(self):
        """
        Get the current data from db.
        """
        try:
            res = pd.read_csv(self.data_path, compression='gzip', index_col=0)
            res.index = res.index.astype("string")
        except pd.errors.EmptyDataError:
            res = pd.DataFrame()
        return res
    @property
    def obj(self):
        """
        Get current saved list
        """
        res = load_json(self.obj_path)
        return res
    def update(self, new_data, key=None):
        """
        Update the data in db.
        """
        if isinstance(new_data, pd.DataFrame):
            new_data_df = new_data.copy()
        elif isinstance(new_data, pd.Series):
            if not key is None:
                new_data = new_data.copy().rename(key)
            new_data_df = new_data.to_frame().T
        elif isinstance(new_data, dict):
            assert key is not None, "Key must be provided for dict input."
            new_data_df = pd.DataFrame(
                {
                    key: pd.Series(new_data)
                }
            )
        elif isinstance(new_data, list):
            assert key is not None, "Key must be provided for list input."
            obj = self.obj
            obj[key] = new_data
            dump_json(obj, self.obj_path)
            print("Data updated.")
            return 
        else:
            raise TypeError(f"Unsupported data type: {type(new_data)}")
        res = pd.concat([self.data, new_data_df], axis=0, join="outer")    
        res = res.groupby(level=0).last()
        res.to_csv(self.data_path, compression='gzip')
        print("Data updated.")
        return