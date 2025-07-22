# This module is designed for bioinformatic project management
from pathlib import Path
import shutil
import pandas as pd
from .hicio import load_json, dump_json
from datetime import datetime

class Ana:
    """
    A class to manage analysis data and operations with version control.
    """
    def __init__(self, home, tag=None, data=None, obj=None, clear=False, verbose=False, max_commits=50):
        """
        Initialize the Ana class with version control.
        """
        self.home = Path(home)
        self.home.mkdir(parents=True, exist_ok=True)
        self.max_commits = max_commits
        self.verbose = verbose
        
        # Initialize tag
        self.tag = tag or "0"

        # Initialize commit history
        self.commit_meta_path = self.home / f"{self.tag}.commits.json"
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
        data_path = self.home / f"{self.tag}.data.csv.gz"
        obj_path = self.home / f"{self.tag}.obj.json"
        
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
        return self.home / f"{self.tag}.data.csv.gz"
    
    @property
    def obj_path(self):
        return self.home / f"{self.tag}.obj.json"
    
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
    
    def commit(self, commit_id=None):
        """Create a new commit of the current state"""
        # Load current metadata
        meta = load_json(self.commit_meta_path)
        
        # Determine commit ID
        if commit_id is None:
            commit_id = str(meta["last_commit_count"] + 1)
            meta["last_commit_count"] += 1
        elif any(c["id"] == commit_id for c in meta["commits"]):
            raise ValueError(f"Commit ID '{commit_id}' already exists")
        
        # Create commit entry
        commit_entry = {
            "id": commit_id,
            "timestamp": datetime.now().isoformat(),
            "tag": self.tag
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
    
    def update(self, new_data, key=None, commit_id=None):
        """
        Update the data in db and automatically commit changes.
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
        self.commit(commit_id=commit_id)
        if self.verbose:
            print("Data updated and committed.")
    
    def revert(self, commit_id):
        """Revert to a previous commit version"""
        # Validate commit exists
        meta = load_json(self.commit_meta_path)
        if not any(c["id"] == commit_id for c in meta["commits"]):
            raise ValueError(f"Commit '{commit_id}' not found")
        
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
        return [c["id"] for c in meta["commits"]]