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
- Full pandas dtype preservation using CSV + JSON schema
- Thread-safe for single-process use (not suitable for multiprocessing)
"""

# This module is designed for bioinformatic project management
from pathlib import Path
import re
import shutil
import pandas as pd
import numpy as np
import json
from datetime import datetime
from typing import Dict, Any, Optional, Union, List
import warnings
from .hicio import dump_json, load_json

class Ana:
    """
    Analysis Data Manager with automatic version control and full dtype preservation.
    
    Manages two types of data storage:
    1. Structured data: Pandas DataFrames stored in CSV format with JSON schema for dtypes
    2. Object storage: JSON-serializable lists for unstructured data
    
    The system maintains a commit history similar to Git, allowing users to track changes,
    revert to previous states, and maintain data provenance throughout analysis workflows.
    
    Data Type Preservation Strategy:
    - Main data stored as CSV for compatibility
    - Separate JSON schema file stores precise dtype information
    - Supports all pandas dtypes including categories, strings, nullable integers, etc.
    - Automatic dtype inference and storage during save
    - Accurate dtype restoration during load
    
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
        >>> df = pd.DataFrame({'A': pd.Categorical(['a', 'b', 'c']),
        ...                    'B': pd.array([1, 2, None], dtype='Int64'),
        ...                    'C': pd.Series(['x', 'y', 'z'], dtype='string')})
        >>> ana.update(df)
        >>> ana.data  # Returns DataFrame with correct dtypes: categorical, Int64, string
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
        self.verbose = verbose
        
        # Handle initialization from another Ana object
        if from_ana is not None:
            if self.home.resolve() == Path(from_ana.home).resolve():
                raise ValueError("When initializing from another Ana object, "
                               "you must provide a different home directory.")
            
            # Create the new directory
            self.home.mkdir(parents=True, exist_ok=True)
            
            # Copy data, schema, and obj from the source Ana object
            from_ana._save_data_with_schema(from_ana.data, self.home / "data.csv")
            dump_json(from_ana.obj, self.home / "obj.json")
            
            # Set other parameters
            self.max_commits = max_commits
            
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
            data_path = self.home / "data.csv"
            schema_path = self.home / "data_schema.json"
            obj_path = self.home / "obj.json"
            
            # Clear existing data if requested
            if clear:
                if data_path.exists():
                    data_path.unlink()
                if schema_path.exists():
                    schema_path.unlink()
                if obj_path.exists():
                    obj_path.unlink()
                # Also clear any existing commit files
                for commit_file in self.home.glob("commit_*.data.csv"):
                    commit_file.unlink()
                for commit_file in self.home.glob("commit_*.data_schema.json"):
                    commit_file.unlink()
                for commit_file in self.home.glob("commit_*.obj.json"):
                    commit_file.unlink()
            
            # Create empty storage files if they don't exist
            if not data_path.exists():
                data_path.touch()
            if not schema_path.exists():
                Ana.create_empty_json_file(schema_path)
            if not obj_path.exists():
                Ana.create_empty_json_file(obj_path)
            
            # Process initial data inputs
            if data is not None:
                if isinstance(data, pd.DataFrame):
                    init_df = data.copy()
                elif isinstance(data, pd.Series):
                    init_df = data.to_frame()
                elif isinstance(data, dict):
                    init_df = pd.DataFrame(data)
                else:
                    raise TypeError("Initial data must be a pandas DataFrame, Series, or dict.")

                self._save_data_with_schema(init_df, self.data_path)
                self.commit("Initial data")

            if obj is not None:
                for key, value in obj.items():
                    if not isinstance(value, list):
                        raise TypeError(f"Value for key '{key}' in obj must be a list.")
                    self.update(value, key=key)

            # Ensure schema file exists even if initial data was not provided
            if not self.schema_path.exists() and self.data_path.exists():
                try:
                    df = pd.read_csv(self.data_path, index_col=0)
                    if not df.empty or len(df.columns) > 0:
                        self._save_data_with_schema(df, self.data_path)
                    else:
                        dump_json({}, self.schema_path)
                except Exception:
                    dump_json({}, self.schema_path)
    
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
        """Path to the CSV file storing structured data."""
        return self.home / "data.csv"
    
    @property
    def schema_path(self):
        """Path to the JSON schema file storing dtype information."""
        return self.home / "data_schema.json"
    
    @property
    def obj_path(self):
        """Path to the JSON file storing unstructured object data."""
        return self.home / "obj.json"
    
    def _dtype_to_dict(self, dtype) -> Dict[str, Any]:
        """
        Convert pandas dtype to JSON-serializable dictionary.
        """
        # 处理索引为StringDtype的情况
        if hasattr(dtype, '__class__'):
            # 处理pandas的StringDtype
            if isinstance(dtype, pd.StringDtype):
                return {'type': 'string'}
            # 处理numpy的字符串类型
            elif dtype == np.dtype('O') or dtype == object:
                # 对于object类型，我们无法区分是普通对象还是字符串
                # 默认当作标准类型处理
                return {'type': 'standard', 'dtype': str(dtype)}
        
        # 处理pandas扩展类型
        if isinstance(dtype, pd.CategoricalDtype):
            dtype_dict = {
                'type': 'categorical',
                'categories': dtype.categories.tolist() if dtype.categories is not None else [],
                'ordered': dtype.ordered
            }
        elif isinstance(dtype, pd.StringDtype):
            dtype_dict = {'type': 'string'}
        elif isinstance(dtype, (pd.Int8Dtype, pd.Int16Dtype, pd.Int32Dtype, pd.Int64Dtype)):
            dtype_name = dtype.name if hasattr(dtype, 'name') else str(dtype)
            bits = dtype_name.replace('Int', '') if 'Int' in dtype_name else '64'
            dtype_dict = {'type': 'nullable_int', 'bits': bits}
        elif isinstance(dtype, (pd.Float32Dtype, pd.Float64Dtype)):
            dtype_name = dtype.name if hasattr(dtype, 'name') else str(dtype)
            bits = dtype_name.replace('Float', '') if 'Float' in dtype_name else '64'
            dtype_dict = {'type': 'nullable_float', 'bits': bits}
        elif isinstance(dtype, pd.BooleanDtype):
            dtype_dict = {'type': 'nullable_bool'}
        elif dtype == np.bool_ or dtype == bool:
            dtype_dict = {'type': 'bool'}
        elif str(dtype).startswith('datetime64'):
            dtype_dict = {'type': 'datetime64'}
        elif str(dtype).startswith('timedelta64'):
            dtype_dict = {'type': 'timedelta64'}
        elif hasattr(pd, 'PeriodDtype') and isinstance(dtype, pd.PeriodDtype):
            dtype_dict = {'type': 'period', 'freq': str(dtype.freq)}
        elif hasattr(pd, 'IntervalDtype') and isinstance(dtype, pd.IntervalDtype):
            dtype_dict = {'type': 'interval'}
        elif hasattr(pd, 'SparseDtype') and isinstance(dtype, pd.SparseDtype):
            dtype_dict = {'type': 'sparse', 'dtype': str(dtype.subtype)}
        elif hasattr(dtype, '_metadata') and hasattr(dtype, 'name'):
            # 处理其他扩展类型
            dtype_dict = {'type': 'extension', 'name': dtype.name}
        elif isinstance(dtype, pd.RangeIndex):
            # 特殊处理RangeIndex
            dtype_dict = {'type': 'range_index'}
        elif str(dtype) == 'int64' or dtype == np.int64:
            # 处理整数类型
            dtype_dict = {'type': 'standard', 'dtype': 'int64'}
        elif str(dtype) == 'float64' or dtype == np.float64:
            # 处理浮点数类型
            dtype_dict = {'type': 'standard', 'dtype': 'float64'}
        elif str(dtype) == 'object':
            # 处理对象类型
            dtype_dict = {'type': 'standard', 'dtype': 'object'}
        else:
            # 默认标准类型
            dtype_dict = {'type': 'standard', 'dtype': str(dtype)}

        return dtype_dict
    
    def _dict_to_dtype(self, dtype_dict: Dict[str, Any]):
        """
        Convert dictionary back to pandas dtype.
        
        Parameters
        ----------
        dtype_dict : dict
            Dictionary representation of dtype.
            
        Returns
        -------
        any pandas dtype
            Restored pandas dtype.
        """
        dtype_type = dtype_dict.get('type', 'standard')
        
        if dtype_type == 'categorical':
            categories = dtype_dict.get('categories', [])
            ordered = dtype_dict.get('ordered', False)
            return pd.CategoricalDtype(categories=categories, ordered=ordered)
        
        elif dtype_type == 'string':
            return pd.StringDtype()
        
        elif dtype_type == 'nullable_int':
            bits = dtype_dict.get('bits', '64')
            if bits == '8':
                return pd.Int8Dtype()
            elif bits == '16':
                return pd.Int16Dtype()
            elif bits == '32':
                return pd.Int32Dtype()
            else:
                return pd.Int64Dtype()
        
        elif dtype_type == 'nullable_float':
            bits = dtype_dict.get('bits', '64')
            if bits == '32':
                return pd.Float32Dtype()
            else:
                return pd.Float64Dtype()
        
        elif dtype_type == 'nullable_bool':
            return pd.BooleanDtype()
        
        elif dtype_type == 'bool':
            return np.bool_
        
        elif dtype_type == 'datetime64':
            return np.dtype('datetime64[ns]')
        
        elif dtype_type == 'timedelta64':
            return np.dtype('timedelta64[ns]')
        
        elif dtype_type == 'period':
            freq = dtype_dict.get('freq', 'D')
            return pd.PeriodDtype(freq=freq)
        
        elif dtype_type == 'interval':
            return pd.IntervalDtype()
        
        elif dtype_type == 'sparse':
            subtype = dtype_dict.get('dtype', 'float64')
            return pd.SparseDtype(subtype)
        
        elif dtype_type == 'range_index':
            # RangeIndex是特殊的，没有直接的dtype，返回int64作为索引类型
            return np.dtype('int64')
        
        elif dtype_type == 'extension':
            # Try to find extension dtype by name
            name = dtype_dict.get('name', '')
            for dtype in [pd.StringDtype(), pd.Int64Dtype(), pd.Float64Dtype(), 
                        pd.BooleanDtype()]:
                if hasattr(dtype, 'name') and dtype.name == name:
                    return dtype
            # 如果没找到，尝试导入ArrowDtype
            try:
                from pandas.core.arrays.arrow import ArrowDtype
                if name.startswith('arrow['):
                    return ArrowDtype.from_string(name)
            except ImportError:
                pass
        
        # Fallback: try to parse standard dtype string
        standard_dtype = dtype_dict.get('dtype', 'object')
        try:
            return np.dtype(standard_dtype)
        except:
            warnings.warn(f"Could not parse dtype: {standard_dtype}, defaulting to object")
            return np.dtype('object')
    
    def _save_data_with_schema(self, df: pd.DataFrame, data_path: Path):
        """
        Save DataFrame to CSV along with its dtype schema.
        
        Parameters
        ----------
        df : pandas.DataFrame
            DataFrame to save.
        data_path : Path
            Path where CSV file will be saved.
        """
        # Save data as CSV
        df.to_csv(data_path, index=True)
        
        # Create and save schema
        schema = {}
        
        # Save column dtypes
        schema['columns'] = {}
        for col in df.columns:
            schema['columns'][col] = self._dtype_to_dict(df[col].dtype)
        
        # 总是保存索引信息，包括RangeIndex
        schema['index'] = self._dtype_to_dict(df.index.dtype)
        schema['index_name'] = str(df.index.name) if df.index.name else None
        
        # 对于string类型的索引，确保正确记录
        if hasattr(df.index, 'dtype') and isinstance(df.index.dtype, pd.StringDtype):
            schema['index'] = {'type': 'string'}
        
        # Save metadata
        schema['shape'] = df.shape
        schema['columns_list'] = df.columns.tolist()
        
        # Save schema to adjacent JSON file with _schema suffix
        schema_path = data_path.parent / f"{data_path.stem}_schema.json"
        dump_json(schema, schema_path)
    
    def _load_data_with_schema(self, data_path: Path) -> pd.DataFrame:
        schema_path = data_path.parent / f"{data_path.stem}_schema.json"
        if not schema_path.exists():
            fallback = data_path.with_suffix('.json')
            if fallback.exists():
                schema_path = fallback
        if not schema_path.exists():
            warnings.warn(f"No schema file found for {data_path}, loading with default dtypes")
            try:
                df = pd.read_csv(data_path, index_col=0)
            except pd.errors.EmptyDataError:
                return pd.DataFrame()
            if not isinstance(df.index, (pd.MultiIndex, pd.CategoricalIndex, pd.DatetimeIndex,
                                         pd.PeriodIndex, pd.IntervalIndex, pd.TimedeltaIndex)):
                df.index = df.index.astype("string")
            return df

        try:
            schema = load_json(schema_path)
        except Exception:
            warnings.warn(f"Failed to load schema file {schema_path}, loading with default dtypes")
            try:
                df = pd.read_csv(data_path, index_col=0)
            except pd.errors.EmptyDataError:
                return pd.DataFrame()
            if not isinstance(df.index, (pd.MultiIndex, pd.CategoricalIndex, pd.DatetimeIndex,
                                         pd.PeriodIndex, pd.IntervalIndex, pd.TimedeltaIndex)):
                df.index = df.index.astype("string")
            return df

        if not schema:
            try:
                df = pd.read_csv(data_path, index_col=0)
            except pd.errors.EmptyDataError:
                return pd.DataFrame()
            if not isinstance(df.index, (pd.MultiIndex, pd.CategoricalIndex, pd.DatetimeIndex,
                                         pd.PeriodIndex, pd.IntervalIndex, pd.TimedeltaIndex)):
                df.index = df.index.astype("string")
            return df

        # 关键修改：先尝试用converters读取索引列为字符串，防止自动类型推断
        try:
            # 获取CSV文件的总列数
            with open(data_path, 'r') as f:
                first_line = f.readline()
                num_columns = len(first_line.split(','))
            
            # 如果只有索引列，需要特殊处理
            if num_columns == 1:
                df = pd.read_csv(data_path, index_col=0, header=None, names=['index'])
                df = df[[]]  # 空DataFrame，只有索引
            else:
                # 使用converters确保索引列以字符串形式读取
                df = pd.read_csv(data_path, index_col=0, converters={0: str})
        except pd.errors.EmptyDataError:
            return pd.DataFrame()
        except Exception as e:
            # 如果上述方法失败，回退到普通读取
            try:
                df = pd.read_csv(data_path, index_col=0)
            except pd.errors.EmptyDataError:
                return pd.DataFrame()

        if 'columns' in schema:
            for col, dtype_dict in schema['columns'].items():
                if col not in df.columns:
                    continue
                try:
                    dtype_type = dtype_dict.get('type', 'standard')
                    
                    if dtype_type == 'categorical':
                        categories = dtype_dict.get('categories', [])
                        ordered = dtype_dict.get('ordered', False)
                        df[col] = pd.Categorical(df[col], categories=categories, ordered=ordered)
                    
                    elif dtype_type == 'nullable_bool':
                        def convert_to_nullable_bool(val):
                            if pd.isna(val) or str(val).strip() == '':
                                return pd.NA
                            elif str(val).lower() in ['true', '1', 't', 'y', 'yes']:
                                return True
                            elif str(val).lower() in ['false', '0', 'f', 'n', 'no']:
                                return False
                            else:
                                return pd.NA
                        
                        df[col] = df[col].apply(convert_to_nullable_bool)
                        df[col] = df[col].astype(pd.BooleanDtype())
                    
                    elif dtype_type == 'nullable_int':
                        # 先转换为float处理缺失值，再转换为nullable int
                        df[col] = pd.to_numeric(df[col], errors='coerce')
                        df[col] = df[col].astype(self._dict_to_dtype(dtype_dict))
                    
                    elif dtype_type == 'nullable_float':
                        df[col] = pd.to_numeric(df[col], errors='coerce')
                        df[col] = df[col].astype(self._dict_to_dtype(dtype_dict))
                    
                    elif dtype_type == 'string':
                        df[col] = df[col].astype(pd.StringDtype())
                    
                    elif dtype_type == 'datetime64':
                        df[col] = pd.to_datetime(df[col], errors='coerce')
                    
                    elif dtype_type == 'timedelta64':
                        df[col] = pd.to_timedelta(df[col], errors='coerce')
                    
                    else:
                        dtype = self._dict_to_dtype(dtype_dict)
                        df[col] = df[col].astype(dtype)
                        
                except Exception as e:
                    warnings.warn(f"Failed to restore dtype for column '{col}': {e}")

        # 修复索引类型恢复
        if 'index' in schema:
            try:
                dtype_type = schema['index'].get('type', 'standard')
                
                if dtype_type == 'categorical':
                    categories = schema['index'].get('categories', [])
                    ordered = schema['index'].get('ordered', False)
                    df.index = pd.CategoricalIndex(df.index, categories=categories, ordered=ordered)
                
                elif dtype_type == 'string':
                    # 确保索引是字符串类型
                    df.index = df.index.astype(str)
                    # 转换为pandas的StringDtype
                    df.index = pd.Index(df.index, dtype=pd.StringDtype())
                
                elif dtype_type == 'datetime64':
                    df.index = pd.to_datetime(df.index, errors='coerce')
                
                elif dtype_type == 'nullable_int':
                    # 索引不支持nullable类型，转为普通整数
                    df.index = pd.to_numeric(df.index, errors='coerce').astype('int64')
                
                elif dtype_type == 'nullable_float':
                    # 索引不支持nullable类型，转为普通浮点数
                    df.index = pd.to_numeric(df.index, errors='coerce').astype('float64')
                
                else:
                    # 对于其他类型，尝试转换
                    dtype = self._dict_to_dtype(schema['index'])
                    try:
                        df.index = df.index.astype(dtype)
                    except:
                        # 如果转换失败，保持原样
                        pass
                        
            except Exception as e:
                warnings.warn(f"Failed to restore index dtype: {e}")

        # Only convert to string if no schema is available and it's not a special index type
        if 'index' not in schema and not isinstance(df.index, (pd.MultiIndex, pd.CategoricalIndex,
                                                               pd.DatetimeIndex, pd.PeriodIndex,
                                                               pd.IntervalIndex, pd.TimedeltaIndex,
                                                               pd.RangeIndex)):
            # For backward compatibility, only convert object/string indices to string
            if df.index.dtype == object or df.index.dtype == 'O':
                df.index = df.index.astype("string")

        if 'index_name' in schema and schema['index_name'] is not None:
            df.index.name = schema['index_name']

        return df
    
    @property
    def data(self):
        """
        Retrieve the current structured data as a pandas DataFrame with correct dtypes.
        
        Returns
        -------
        pandas.DataFrame
            Current data state with restored dtypes. Returns empty DataFrame if no data exists.
            
        Notes
        -----
        Uses schema file to restore all pandas dtypes including categorical, string,
        nullable integers, etc. Falls back to standard CSV reading if schema is missing.
        """
        try:
            df = self._load_data_with_schema(self.data_path)
            return df
        except (pd.errors.EmptyDataError, FileNotFoundError):
            return pd.DataFrame()
    
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
        Each commit creates physical copies of data files and schema files.
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
            "filename": commit_filename
        }
        meta["commits"].append(commit_entry)
        
        # Create physical snapshots of current data state
        commit_data_path = self.home / f"{commit_filename}.data.csv"
        commit_obj_path = self.home / f"{commit_filename}.obj.json"
        
        # Save current data with schema
        try:
            current_df = self._load_data_with_schema(self.data_path)
        except (pd.errors.EmptyDataError, FileNotFoundError):
            current_df = pd.DataFrame()
        self._save_data_with_schema(current_df, commit_data_path)
        
        # Copy object file if it exists
        if self.obj_path.exists():
            shutil.copy(self.obj_path, commit_obj_path)
        else:
            dump_json({}, commit_obj_path)
        
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
            filename_in_use = any(commit["filename"] == oldest_filename 
                                for commit in meta["commits"])
            
            if not filename_in_use:
                (self.home / f"{oldest_filename}.data.csv").unlink(missing_ok=True)
                (self.home / f"{oldest_filename}.data_schema.json").unlink(missing_ok=True)
                (self.home / f"{oldest_filename}.obj.json").unlink(missing_ok=True)
        
        # Save updated metadata
        dump_json(meta, self.commit_meta_path)
        
        if self.verbose:
            print(f"Created commit: {commit_id} (file: {commit_filename})")
        
        return commit_id
    
    def update(self, new_data, key=None, description=None, force=False):
        """
        Update data storage and automatically commit changes.
        
        Parameters
        ----------
        new_data : pandas.DataFrame, pandas.Series, dict, or list
            New data to add to the database.
        key : str, optional
            Specifies the column name (for dict/Series) or object key (for list).
        description : str, optional
            Custom description for the automatic commit.
        force : bool, optional, default False
            If True, overwrite existing data for the given key.
        """
        # Check for existing data and handle based on force flag
        if not force and key is not None:
            if isinstance(new_data, list):
                if key in self.obj:
                    if self.verbose:
                        print(f"Key '{key}' already exists in obj. Use force=True to overwrite.")
                    return
            elif isinstance(new_data, (pd.Series, dict)):
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
                print("Object data updated.")
            return 
        else:
            raise TypeError(f"Unsupported data type: {type(new_data)}")
        
        # Record original data shape for comparison
        original_shape = self.data.shape
        original_rows, original_cols = original_shape
        
        # Merge new data with existing data
        res = pd.concat([self.data, new_data_df], axis=0, join="outer")
        res = res.groupby(level=0, observed=True).last()
        
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
        
        # Save the updated data with schema
        self._save_data_with_schema(res, self.data_path)
        
        if self.verbose:
            print("Data updated with schema preservation.")
        
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
            Commit identifier.
            
        Returns
        -------
        int
            Commit ID of the new commit created after revert operation.
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
                for commit in meta["commits"]:
                    if str(commit["id"]) == str(commit_repr):
                        commit_id = commit["id"]
                        break
        
        if commit_id is None or commit_filename is None:
            raise ValueError(f"Commit '{commit_repr}' not found in metadata")
        
        # Restore files from commit snapshot
        commit_data_path = self.home / f"{commit_filename}.data.csv"
        commit_obj_path = self.home / f"{commit_filename}.obj.json"
        
        if not commit_data_path.exists():
            raise ValueError(f"Commit data file for '{commit_repr}' not found on disk")
        
        # Load the commit data and save it as current (this will create schema file too)
        commit_df = self._load_data_with_schema(commit_data_path)
        self._save_data_with_schema(commit_df, self.data_path)
        
        # Restore object file
        if commit_obj_path.exists():
            shutil.copy(commit_obj_path, self.obj_path)
        
        # Create new commit to record the revert action
        new_commit_id = self.commit(f"Reverted to commit {commit_id}")
        
        if self.verbose:
            print(f"Reverted to commit {commit_id} (file: {commit_filename}) and created new commit {new_commit_id}")
        
        return new_commit_id
    
    def list_commits(self):
        """
        Display all available commits in chronological order.
        """
        meta = load_json(self.commit_meta_path)
        meta["commits"].sort(key=lambda c: c["id"])
        if self.verbose:
            print(f"Listing {len(meta['commits'])} commits:")
        for commit in meta["commits"]:
            print(f"{commit['id']}: {commit['description']} at {commit['timestamp']} (file: {commit['filename']})")

    def get_schema(self):
        """
        Get the current data schema information.
        
        Returns
        -------
        dict
            Current schema information including column dtypes.
        """
        if self.schema_path.exists():
            return load_json(self.schema_path)
        return {}