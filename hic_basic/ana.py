"""
A single-process metadata management system for bioinformatics projects.

This module provides version-controlled storage for structured data (DataFrame) 
and unstructured objects (JSON-serializable lists) with automatic commit history.
Designed for managing analysis pipelines where data provenance and reproducibility are critical.

Key Features:
- Automatic versioning with commit history
- Support for both tabular data (DataFrame) and object storage (lists)
- Configurable commit history limits with circular commit file naming
- Revert functionality to previous states
- Thread-safe for single-process use (not suitable for multiprocessing)
"""

# This module is designed for bioinformatic project management
from pathlib import Path
import re
import shutil
import pandas as pd
from .hicio import load_json, dump_json
from datetime import datetime

class Ana:
    """
    Analysis Data Manager with automatic version control.
    
    Manages two types of data storage:
    1. Structured data: Pandas DataFrames stored in compressed CSV format
    2. Object storage: JSON-serializable lists for unstructured data
    
    The system maintains a commit history similar to Git, allowing users to track changes,
    revert to previous states, and maintain data provenance throughout analysis workflows.
    
    Workflow:
    1. Initialize with a home directory (creates storage structure)
    2. Update data using update() method (automatically creates commits)
    3. Track changes using list_commits()
    4. Revert if needed using revert()
    
    Commit Storage Strategy:
    - Uses circular file naming with modulo max_commits to prevent infinite filename growth
    - Maintains mapping between logical commit IDs and physical file names in metadata
    - Automatically reuses file slots when max_commits limit is reached
    
    Example:
        >>> ana = Ana('/path/to/project')
        >>> ana.update({'sample1': 42}, key='measurements')
        >>> ana.commit('Initial data import')
    
    TODO: fix force option
    """
    
    def __init__(self, home, data=None, obj=None, clear=False, verbose=False, max_commits=50, from_ana=None):
        """
        Initialize the analysis data manager.
        
        Parameters
        ----------
        home : str or Path
            Directory path for storing all project data and metadata.
            Will be created if it doesn't exist.
        data : pandas.DataFrame, pandas.Series, dict, or list, optional
            Initial data to load into the system. If provided, will trigger
            an initial commit after loading.
        obj : dict, optional
            Initial object data to load (key-value pairs where values are lists).
            Each key-value pair will be added to object storage.
        clear : bool, default False
            If True, clear existing data in the home directory and start fresh.
        verbose : bool, default False
            If True, print operational messages for debugging.
        max_commits : int, default 50
            Maximum number of commits to maintain in history. Older commits
            are automatically purged when limit is exceeded.
        from_ana : Ana, optional
            Another Ana object to initialize from. If provided, will copy
            the data and obj from that object but create a new commit history.
            Must specify a different home directory.
            
        Raises
        ------
        TypeError
            If obj values are not lists (required for JSON serialization consistency)
        ValueError
            If from_ana is provided but home is the same as from_ana.home
        """
        self.home = Path(home)
        
        # Handle initialization from another Ana object
        if from_ana is not None:
            if self.home.resolve() == Path(from_ana.home).resolve():
                raise ValueError("When initializing from another Ana object, "
                               "you must provide a different home directory.")
            
            # Create the new directory
            self.home.mkdir(parents=True, exist_ok=True)
            
            # Copy data and obj from the source Ana object
            from_ana.data.to_csv(self.home / "data.csv.gz", compression='gzip')
            dump_json(from_ana.obj, self.home / "obj.json")
            
            # Set other parameters
            self.max_commits = max_commits
            self.verbose = verbose
            
            # Initialize fresh commit history
            self.commit_meta_path = self.home / "commits.json"
            self.commit_meta = {
                "commits": [],
                "last_commit_count": 0,
                "max_commits": max_commits,
                "file_mapping": {}
            }
            dump_json(self.commit_meta, self.commit_meta_path)
            
            # Create initial commit with copied data
            commit_id = self.commit("Initialize from another Ana object")
            # Reload commit metadata to ensure in-memory representation is up-to-date
            self.commit_meta = load_json(self.commit_meta_path)
        else:
            # Original initialization logic
            self.home.mkdir(parents=True, exist_ok=True)
            self.max_commits = max_commits
            self.verbose = verbose
            
            # Initialize commit history - stored as JSON metadata
            self.commit_meta_path = self.home / "commits.json"
            if self.commit_meta_path.exists() and not clear:
                self.commit_meta = load_json(self.commit_meta_path)
                # Ensure max_commits is updated if changed
                self.commit_meta["max_commits"] = max_commits
            else:
                self.commit_meta = {
                    "commits": [],
                    "last_commit_count": 0,
                    "max_commits": max_commits,
                    "file_mapping": {}  # Maps commit_id -> physical_filename
                }
                dump_json(self.commit_meta, self.commit_meta_path)

            # Initialize data storage files
            data_path = self.home / "data.csv.gz"
            obj_path = self.home / "obj.json"
            
            # Clear existing data if requested
            if clear:
                if data_path.exists():
                    data_path.unlink()
                if obj_path.exists():
                    obj_path.unlink()
                # Also clear any existing commit files
                for commit_file in self.home.glob("commit_*.data.csv.gz"):
                    commit_file.unlink()
                for commit_file in self.home.glob("commit_*.obj.json"):
                    commit_file.unlink()
            
            # Create empty storage files if they don't exist
            if not data_path.exists():
                data_path.touch()
            if not obj_path.exists():
                Ana.create_empty_json_file(obj_path)
            
            # Process initial data inputs
            if data is not None:
                self.update(data, key=None)
            if obj is not None:
                for key, value in obj.items():
                    if not isinstance(value, list):
                        raise TypeError(f"Value for key '{key}' in obj must be a list.")
                    self.update(value, key=key)
    
    @staticmethod
    def create_empty_json_file(path):
        """Create an empty JSON file with empty dictionary content."""
        dump_json({}, path)
    
    def _get_commit_filename(self, commit_id):
        """
        Generate physical filename for a commit using circular naming.
        
        Parameters
        ----------
        commit_id : int
            Logical commit identifier.
            
        Returns
        -------
        str
            Base filename for commit files (without extension).
            
        Notes
        -----
        Uses modulo operation with max_commits to ensure filenames cycle
        within a fixed range, preventing infinite filename number growth.
        """
        file_slot = commit_id % self.max_commits
        return f"commit_{file_slot:04d}"
    
    @property
    def data_path(self):
        """Path to the compressed CSV file storing structured data."""
        return self.home / "data.csv.gz"
    
    @property
    def obj_path(self):
        """Path to the JSON file storing unstructured object data."""
        return self.home / "obj.json"
    
    @property
    def data(self):
        """
        Retrieve the current structured data as a pandas DataFrame.
        
        Returns
        -------
        pandas.DataFrame
            Current data state. Returns empty DataFrame if no data exists.
            
        Notes
        -----
        Index is preserved as string type to maintain consistency across reads.
        Handles empty file cases gracefully.
        """
        try:
            res = pd.read_csv(self.data_path, compression='gzip', index_col=0)
            res.index = res.index.astype("string")
        except (pd.errors.EmptyDataError, FileNotFoundError):
            res = pd.DataFrame()
        return res
    
    @property
    def obj(self):
        """
        Retrieve the current object storage data as a dictionary.
        
        Returns
        -------
        dict
            Current object storage state. Returns empty dict if no objects exist.
        """
        try:
            return load_json(self.obj_path)
        except FileNotFoundError:
            return {}
    
    def commit(self, description=None):
        """
        Create a new commit snapshot of the current data state.
        
        Parameters
        ----------
        description : str, optional
            Human-readable description of the commit. If None, uses commit ID.
            Must be unique and cannot be purely numeric.
            
        Returns
        -------
        int
            Commit ID of the newly created commit.
            
        Raises
        ------
        ValueError
            If description already exists or is purely numeric.
            
        Notes
        -----
        Each commit creates physical copies of data files to enable revert functionality.
        Commit history is automatically trimmed when exceeding max_commits limit.
        Uses circular file naming to prevent infinite filename number growth.
        """
        assert isinstance(description, (str, type(None))), "Description must be a string or None"
        # Load current metadata
        meta = load_json(self.commit_meta_path)
        
        # Generate commit ID and validate description
        commit_id = meta["last_commit_count"] + 1
        if description is None:
            description = str(commit_id)
        elif any(c["description"] == description for c in meta["commits"]):
            raise ValueError(f"Commit description '{description}' already exists")
        else:
            # Prevent ambiguous numeric descriptions that could conflict with commit IDs
            stripped_description = description.strip()
            if re.fullmatch(r'-?\d+', stripped_description):
                raise ValueError("Commit description cannot be purely numeric")
        meta["last_commit_count"] = commit_id
        
        # Get physical filename for this commit using circular naming
        commit_filename = self._get_commit_filename(commit_id)
        
        # Create commit metadata entry
        commit_entry = {
            "id": commit_id,
            "description": description,
            "timestamp": datetime.now().isoformat(),
            "filename": commit_filename  # Store mapping for revert operations
        }
        meta["commits"].append(commit_entry)
        
        # Create physical snapshots of current data state
        commit_data_path = self.home / f"{commit_filename}.data.csv.gz"
        commit_obj_path = self.home / f"{commit_filename}.obj.json"
        
        shutil.copy(self.data_path, commit_data_path)
        shutil.copy(self.obj_path, commit_obj_path)
        
        # Update file mapping
        meta["file_mapping"][str(commit_id)] = commit_filename
        
        # Enforce commit limit by removing oldest commits
        if len(meta["commits"]) > meta["max_commits"]:
            oldest = meta["commits"].pop(0)
            oldest_commit_id = oldest["id"]
            oldest_filename = oldest["filename"]
            
            # Remove from file mapping
            if str(oldest_commit_id) in meta["file_mapping"]:
                del meta["file_mapping"][str(oldest_commit_id)]
            
            # Remove physical files only if they are not being used by newer commits
            # Check if any other commit is using the same filename
            filename_in_use = any(commit["filename"] == oldest_filename 
                                for commit in meta["commits"])
            
            if not filename_in_use:
                (self.home / f"{oldest_filename}.data.csv.gz").unlink(missing_ok=True)
                (self.home / f"{oldest_filename}.obj.json").unlink(missing_ok=True)
        
        # Save updated metadata
        dump_json(meta, self.commit_meta_path)
        
        if self.verbose:
            print(f"Created commit: {commit_id} (file: {commit_filename})")
        
        return commit_id
    
    def update(self, new_data, key=None, description=None, force=False):
        """
        Update data storage and automatically commit changes.
        
        This is the primary method for adding or modifying data in the system.
        Supports multiple data types with different storage strategies.
        
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
            
        force : bool, optional, default False
            If True, overwrite existing data for the given key.
            If False, check if data already exists for the key and skip update if found.

        Returns
        -------
        None

        Raises
        ------
        AssertionError
            If key is required but not provided for dict/list inputs.
        TypeError
            If unsupported data type is provided.

        Notes
        -----
        For DataFrame/Series/dict inputs:
        - Data is merged with existing data using concat+groupby.last()
        - New values overwrite existing values on matching indices
        - The result is automatically committed to version history
        - If row count changes, prints information about the change

        For list inputs:
        - Data is stored in the separate object store (JSON format)
        - The entire list replaces any previous value for the given key
        - Object store updates do not trigger automatic commits
        """
        # Check for existing data and handle based on force flag
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
        
        # Process different data types
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
        
        # Record original data shape for comparison
        original_shape = self.data.shape
        original_rows, original_cols = original_shape
        
        # Merge new data with existing data and save
        res = pd.concat([self.data, new_data_df], axis=0, join="outer")    
        res = res.groupby(level=0).last()
        
        # Get new shape and compare
        new_shape = res.shape
        new_rows, new_cols = new_shape
        
        # Print shape change information if rows changed
        if original_rows != new_rows:
            row_change = new_rows - original_rows
            row_change_type = "added" if row_change > 0 else "removed"
            print(f"Data shape changed: {original_shape} -> {new_shape}")
            print(f"Rows {row_change_type}: {abs(row_change)}")
        
        # Also print column changes if any
        if original_cols != new_cols:
            col_change = new_cols - original_cols
            col_change_type = "added" if col_change > 0 else "removed"
            print(f"Columns {col_change_type}: {abs(col_change)}")
        
        # Save the updated data
        res.to_csv(self.data_path, compression='gzip')
        
        if self.verbose:
            print("Data updated.")
        
        # Create commit for the data change
        self.commit(description=description)
        if self.verbose:
            print("Data updated and committed.")
    
    def revert(self, commit_repr):
        """
        Revert the system state to a previous commit version.
        
        Parameters
        ----------
        commit_repr : str or int
            Commit identifier, which can be either:
            - Commit ID (integer as string or number)
            - Commit description (exact string match)
            
        Returns
        -------
        int
            Commit ID of the new commit created after revert operation.
            
        Raises
        ------
        ValueError
            If the specified commit cannot be found in metadata.
            
        Notes
        -----
        This operation:
        1. Restores data files from the specified commit snapshot
        2. Creates a new commit to record the revert action
        3. Maintains the entire commit history including the revert
        
        The revert is itself versioned, allowing undo of revert operations.
        Uses the file mapping to locate the correct physical files.
        """
        # Validate commit exists in metadata
        meta = load_json(self.commit_meta_path)
        commit_id = None
        commit_filename = None
        
        for commit in meta["commits"]:
            if commit["id"] == commit_repr or commit["description"] == commit_repr:
                commit_id = commit["id"]
                commit_filename = commit["filename"]
                break
        else:
            # Also check file mapping for string representations
            if str(commit_repr) in meta["file_mapping"]:
                commit_filename = meta["file_mapping"][str(commit_repr)]
                # Find the commit entry to get the full metadata
                for commit in meta["commits"]:
                    if str(commit["id"]) == str(commit_repr):
                        commit_id = commit["id"]
                        break
        
        if commit_id is None or commit_filename is None:
            raise ValueError(f"Commit '{commit_repr}' not found in metadata")
        
        # Restore files from commit snapshot
        commit_data_path = self.home / f"{commit_filename}.data.csv.gz"
        commit_obj_path = self.home / f"{commit_filename}.obj.json"
        
        if not commit_data_path.exists() or not commit_obj_path.exists():
            raise ValueError(f"Commit files for '{commit_repr}' not found on disk")
        
        shutil.copy(commit_data_path, self.data_path)
        shutil.copy(commit_obj_path, self.obj_path)
        
        # Create new commit to record the revert action
        new_commit_id = self.commit(f"Reverted to commit {commit_id}")
        
        if self.verbose:
            print(f"Reverted to commit {commit_id} (file: {commit_filename}) and created new commit {new_commit_id}")
        
        return new_commit_id
    
    def list_commits(self):
        """
        Display all available commits in chronological order.
        
        Returns
        -------
        None
        
        Notes
        -----
        Output format: 
        [commit_id]: [description] at [timestamp] (file: [filename])
        
        Commits are sorted by ID (chronological order). Most recent commit appears last.
        Includes physical filename information for debugging.
        """
        meta = load_json(self.commit_meta_path)
        # Sort commits by ID to ensure chronological order
        meta["commits"].sort(key=lambda c: c["id"])
        if self.verbose:
            print(f"Listing {len(meta['commits'])} commits:")
        for commit in meta["commits"]:
            print(f"{commit['id']}: {commit['description']} at {commit['timestamp']} (file: {commit['filename']})")