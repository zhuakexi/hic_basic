# This module is designed for bioinformatic project management
from pathlib import Path
import re
import shutil
import pandas as pd
from .hicio import load_json, dump_json
from datetime import datetime

class Ana:
    """
    A class to manage analysis data and operations with version control.
    """
    def __init__(self, home, data=None, obj=None, clear=False, verbose=False, max_commits=50):
        """
        Initialize the Ana class with version control.
        """
        self.home = Path(home)
        self.home.mkdir(parents=True, exist_ok=True)
        self.max_commits = max_commits
        self.verbose = verbose
        
        # Initialize commit history
        self.commit_meta_path = self.home / "commits.json"
        if self.commit_meta_path.exists() and not clear:
            self.commit_meta = load_json(self.commit_meta_path)
        else:
            self.commit_meta = {
                "commits": [],
                "last_commit_count": 0,
                "max_commits": max_commits
            }
            dump_json(self.commit_meta, self.commit_meta_path)

        # Initialize data files
        data_path = self.home / "data.csv.gz"
        obj_path = self.home / "obj.json"
        
        if clear:
            if data_path.exists():
                data_path.unlink()
            if obj_path.exists():
                obj_path.unlink()
        
        if not data_path.exists():
            data_path.touch()
        if not obj_path.exists():
            Ana.create_empty_json_file(obj_path)
        
        # Process initial data
        if data is not None:
            self.update(data, key=None)
        if obj is not None:
            for key, value in obj.items():
                if not isinstance(value, list):
                    raise TypeError(f"Value for key '{key}' in obj must be a list.")
                self.update(value, key=key)
    
    @staticmethod
    def create_empty_json_file(path):
        dump_json({}, path)
    
    @property
    def data_path(self):
        return self.home / "data.csv.gz"
    
    @property
    def obj_path(self):
        return self.home / "obj.json"
    
    @property
    def data(self):
        try:
            res = pd.read_csv(self.data_path, compression='gzip', index_col=0)
            res.index = res.index.astype("string")
        except (pd.errors.EmptyDataError, FileNotFoundError):
            res = pd.DataFrame()
        return res
    
    @property
    def obj(self):
        try:
            return load_json(self.obj_path)
        except FileNotFoundError:
            return {}
    
    def commit(self, description=None):
        """Create a new commit of the current state"""
        assert isinstance(description, (str, type(None))), "Description must be a string or None"
        # Load current metadata
        meta = load_json(self.commit_meta_path)
        
        # Determine commit description
        commit_id = meta["last_commit_count"] + 1
        if description is None:
            description = str(commit_id)
        elif any(c["description"] == description for c in meta["commits"]):
            raise ValueError(f"Commit description '{description}' already exists")
        else:
            # check if description is ambiguous (e.g., purely numeric)
            # if can be converted to int, raise error
            stripped_description = description.strip()
            if re.fullmatch(r'-?\d+', stripped_description):
                raise ValueError("Commit description cannot be purely numeric")
        meta["last_commit_count"] += 1
        
        # Create commit entry
        commit_entry = {
            "id": commit_id,
            "description": description,
            "timestamp": datetime.now().isoformat()
        }
        meta["commits"].append(commit_entry)
        
        # Save data snapshot
        commit_data_path = self.home / f"{commit_id}.data.csv.gz"
        commit_obj_path = self.home / f"{commit_id}.obj.json"
        
        shutil.copy(self.data_path, commit_data_path)
        shutil.copy(self.obj_path, commit_obj_path)
        
        # Enforce commit limit
        if len(meta["commits"]) > meta["max_commits"]:
            oldest = meta["commits"].pop(0)
            (self.home / f"{oldest['id']}.data.csv.gz").unlink(missing_ok=True)
            (self.home / f"{oldest['id']}.obj.json").unlink(missing_ok=True)
        
        # Save updated metadata
        dump_json(meta, self.commit_meta_path)
        
        if self.verbose:
            print(f"Created commit: {commit_id}")
        
        return commit_id
    
    def update(self, new_data, key=None, description=None, force=False):
        """
        Update the data in db and automatically commit changes.

        Parameters
        ----------
        new_data : pandas.DataFrame, pandas.Series, dict, or list
            New data to add to the database. Behavior varies by type:
            
            - pandas.DataFrame: Entire DataFrame will be merged with existing data.
                Index will be used for alignment. The `key` parameter is ignored.
                
            - pandas.Series: Series will be converted to DataFrame and merged.
                If `key` is provided, the Series will be renamed to this column name.
                Otherwise, the Series name will be used.
                
            - dict: Dictionary will be converted to DataFrame using keys as index.
                Requires `key` parameter to specify column name for the values.
                
            - list: Data will be stored in the object store (JSON) rather than the
                main data table. Requires `key` parameter to specify the object key.
                Note: The entire list will replace any existing value for that key.

        key : str, optional
            Specifies the column name (for dict/Series) or object key (for list).
            Required for dict, list, and unnamed Series inputs.

        description : str, optional
            Custom description for the automatic commit. If not provided, uses
            an auto-generated commit identifier.
            
        force : bool, optional, default=False
            If True, overwrite existing data for the given key.
            If False, check if data already exists for the key and skip update if found.

        Notes
        -----
        For DataFrame/Series/dict inputs:
        - Data is merged with existing data using concat+groupby.last()
        - New values overwrite existing values on matching indices
        - The result is automatically committed to version history

        For list inputs:
        - Data is stored in the separate object store (JSON format)
        - The entire list replaces any previous value for the given key
        - Does not trigger data file updates, only updates object store
        """
        # Check if data already exists for the key and handle accordingly
        if not force and key is not None:
            if isinstance(new_data, list):
                # For list data, check if key exists in obj
                if key in self.obj:
                    if self.verbose:
                        print(f"Key '{key}' already exists in obj. Use force=True to overwrite.")
                    return
            elif isinstance(new_data, (pd.Series, dict)):
                # For Series/dict data, check if column exists in data
                col_name = key if key is not None else (new_data.name if hasattr(new_data, 'name') else None)
                if col_name and col_name in self.data.columns:
                    if self.verbose:
                        print(f"Column '{col_name}' already exists in data. Use force=True to overwrite.")
                    return
        
        if isinstance(new_data, pd.DataFrame):
            new_data_df = new_data.copy()
        elif isinstance(new_data, pd.Series):
            if not key is None:
                new_data = new_data.copy().rename(key)
            new_data_df = new_data.to_frame()
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
            if self.verbose:
                print("Data updated.")
            return 
        else:
            raise TypeError(f"Unsupported data type: {type(new_data)}")
        res = pd.concat([self.data, new_data_df], axis=0, join="outer")    
        res = res.groupby(level=0).last()
        res.to_csv(self.data_path, compression='gzip')
        if self.verbose:
            print("Data updated.")
        self.commit(description=description)
        if self.verbose:
            print("Data updated and committed.")
    
    def revert(self, commit_repr):
        """
        Revert to a previous commit version.
        Input:
            commit_repr: str, commit ID or description
        """
        # Validate commit exists
        meta = load_json(self.commit_meta_path)
        # check commid id by id or description
        for commit in meta["commits"]:
            if commit["id"] == commit_repr or commit["description"] == commit_repr:
                commit_id = commit["id"]
                break
        else:
            raise ValueError(f"Commit '{commit_repr}' not found in metadata")
        
        # Restore files from commit
        commit_data_path = self.home / f"{commit_id}.data.csv.gz"
        commit_obj_path = self.home / f"{commit_id}.obj.json"
        
        shutil.copy(commit_data_path, self.data_path)
        shutil.copy(commit_obj_path, self.obj_path)
        
        # Create new commit for revert action
        new_commit_id = self.commit()
        
        if self.verbose:
            print(f"Reverted to commit {commit_id} and created new commit {new_commit_id}")
        
        return new_commit_id
    
    def list_commits(self):
        """List all available commits"""
        meta = load_json(self.commit_meta_path)
        # sort by commit ID
        meta["commits"].sort(key=lambda c: c["id"])
        if self.verbose:
            print(f"Listing {len(meta['commits'])} commits:")
        for commit in meta["commits"]:
            print(f"{commit['id']}: {commit['description']} at {commit['timestamp']}")