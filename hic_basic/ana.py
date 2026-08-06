# This module is designed for bioinformatic project management
import inspect
import os
import re
import hashlib
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from datetime import datetime
from typing import Dict, Any, Optional, Union, List
import warnings
from .hicio import dump_json, load_json

ANA_SYSTEM_VERSION = 3

# File format extension → key mapping, used for auto-detection of file columns
_FILE_EXT_TO_FORMAT = {
    ".cool": "cool", ".mcool": "cool",
    ".scool": "scool",
    ".pairs": "pairs", ".pairsam": "pairs",
    ".bedpe": "bedpe",
    ".bam": "bam", ".cram": "cram",
    ".bed": "bed", ".bedgraph": "bedgraph", ".bedGraph": "bedgraph",
    ".bw": "bigwig", ".bigwig": "bigwig", ".bigWig": "bigwig",
    ".vcf": "vcf",
    ".fa": "fasta", ".fasta": "fasta", ".fna": "fasta",
    ".fq": "fastq", ".fastq": "fastq",
    ".gff": "gff", ".gff3": "gff", ".gtf": "gff",
    ".tsv": "tsv", ".tab": "tsv",
    ".csv": "csv",
}


class AnaUpdateError(ValueError):
    """Raised when an Ana update would violate sample-axis safety rules."""


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
        >>> df = pd.DataFrame(
        ...     {'A': pd.Categorical(['a', 'b', 'c']),
        ...      'B': pd.array([1, 2, None], dtype='Int64'),
        ...      'C': pd.array(['x', 'y', 'z'], dtype='string')},
        ...     index=['sample_a', 'sample_b', 'sample_c'],
        ... )
        >>> ana.update(df, intent='add_samples')
        >>> ana.data  # Returns DataFrame with correct dtypes: categorical, Int64, string
    """
    
    def __init__(self, home, data=None, obj=None, clear=False, verbose=False, max_commits=50, from_ana=None, initialize=True):
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
        initialize : bool, default True
            If False, construct in read-only introspection mode: the project
            directory is never created, modified, or synced (see
            ``for_introspection()``). No ``data``/``obj``/``clear``/``from_ana``.
            
        Raises
        ------
        TypeError
            If obj values are not lists (required for JSON serialization consistency)
        ValueError
            If from_ana is provided but home is the same as from_ana.home
        """
        self.home = Path(home)
        self.verbose = verbose

        if not initialize:
            assert data is None and obj is None and not clear, (
                "initialize=False (read-only introspection) accepts no "
                "data, obj, or clear"
            )

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
                "file_mapping": {},
                **self._default_system_meta(),
            }
            dump_json(self.commit_meta, self.commit_meta_path)
            
            # Create initial commit with copied data
            commit_id = self.commit("Initialize from another Ana object")
            # Reload commit metadata to ensure in-memory representation is up-to-date
            self.commit_meta = load_json(self.commit_meta_path)
        else:
            # Original initialization logic
            self.max_commits = max_commits
            self.commit_meta_path = self.home / "commits.json"
            if not initialize:
                # Read-only introspection (Ana.for_introspection): never
                # create or modify anything in the project directory.
                self.commit_meta = (
                    load_json(self.commit_meta_path)
                    if self.commit_meta_path.exists()
                    else {}
                )
            else:
                self.home.mkdir(parents=True, exist_ok=True)
                # Initialize commit history - stored as JSON metadata
                if self.commit_meta_path.exists() and not clear:
                    self.commit_meta = load_json(self.commit_meta_path)
                    # Ensure max_commits is updated if changed
                    self.commit_meta["max_commits"] = max_commits
                else:
                    self.commit_meta = {
                        "commits": [],
                        "last_commit_count": 0,
                        "max_commits": max_commits,
                        "file_mapping": {},  # Maps commit_id -> physical_filename
                        **self._default_system_meta(),
                    }
                    dump_json(self.commit_meta, self.commit_meta_path)
                self._ensure_system_meta()

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
            if initialize:
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

                self._validate_sample_index(init_df.index, context="initial ana.data")
                self._save_data_with_schema(init_df, self.data_path)
                self.commit("Initial data")

            if obj is not None:
                for key, value in obj.items():
                    if not isinstance(value, list):
                        raise TypeError(f"Value for key '{key}' in obj must be a list.")
                    self.update(value, key=key)

            # Ensure schema file exists even if initial data was not provided
            if initialize and not self.schema_path.exists() and self.data_path.exists():
                try:
                    df = pd.read_csv(self.data_path, index_col=0)
                    if not df.empty or len(df.columns) > 0:
                        self._save_data_with_schema(df, self.data_path)
                    else:
                        dump_json({}, self.schema_path)
                except Exception:
                    dump_json({}, self.schema_path)

    @classmethod
    def for_introspection(cls, home, verbose=False):
        """
        Construct an Ana without touching the project directory.

        Read-only introspection for tools like ``ana-describe``: schema,
        column metadata, and objects are loaded if present, but nothing is
        created or rewritten (no mkdir, no commits.json sync, no empty
        storage files).

        Parameters
        ----------
        home : str or Path
            Path to the Ana project directory (not created if missing).
        verbose : bool, default False
            Print operational messages.

        Returns
        -------
        Ana
            Ana instance safe for read-only use.
        """
        return cls(home, verbose=verbose, initialize=False)

    @staticmethod
    def _default_system_meta():
        return {
            "ana_system_version": ANA_SYSTEM_VERSION,
            "axis": {
                "name": "sample_id",
                "canonical_dtype": "string",
                "role": "sample_axis",
            },
        }

    def _ensure_system_meta(self):
        """Ensure persisted Ana metadata declares the current system version."""
        meta = load_json(self.commit_meta_path)
        changed = False
        for key, value in self._default_system_meta().items():
            if key not in meta:
                meta[key] = value
                changed = True
        if meta.get("ana_system_version") != ANA_SYSTEM_VERSION:
            meta["ana_system_version"] = ANA_SYSTEM_VERSION
            changed = True
        if changed:
            dump_json(meta, self.commit_meta_path)
        self.commit_meta = meta
    
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
    
    def _detect_file_format(self, series: pd.Series) -> Optional[str]:
        """
        Detect if a Series contains file paths and return the best-matching
        file format key.

        Parameters
        ----------
        series : pd.Series
            Column to inspect.

        Returns
        -------
        str or None
            File format key (e.g., ``"cool"``, ``"bam"``), or None if not
            a recognizable file column.
        """
        # Sample up to 10 non-null values
        sample = series.dropna().astype(str).head(10)
        if len(sample) == 0:
            return None

        # Count how many values look like file paths with known extensions
        format_counts: dict[str, int] = {}
        path_like_count = 0
        for val in sample:
            val_lower = val.lower().strip()
            # Check if it looks like a file path
            if any(c in val_lower for c in ["/", "\\", ".\\"]):
                path_like_count += 1
            # Check for known extensions
            for ext, fmt in _FILE_EXT_TO_FORMAT.items():
                if val_lower.endswith(ext):
                    format_counts[fmt] = format_counts.get(fmt, 0) + 1
                    break

        # Need at least 50% path-like and at least one recognized format
        n = len(sample)
        if path_like_count < n * 0.5:
            return None
        if format_counts:
            return max(format_counts, key=format_counts.get)  # type: ignore[arg-type]
        return None

    def _save_data_with_schema(self, df: pd.DataFrame, data_path: Path, schema_overrides=None):
        """
        Save DataFrame to CSV along with its dtype schema.

        Parameters
        ----------
        df : pandas.DataFrame
            DataFrame to save.
        data_path : Path
            Path where CSV file will be saved.
        schema_overrides : dict, optional
            Exact per-column schema entries to preserve while serializing a
            proposed update. Columns not present are inferred from ``df``.
        """
        df = df.copy()
        if len(df.index) > 0:
            self._validate_sample_index(df.index, context="ana.data")
            df.index = pd.Index(self._canonical_labels(df.index), name=df.index.name, dtype=pd.StringDtype())

        # Save data as CSV
        df.to_csv(data_path, index=True)

        # Create and save schema
        schema = {}

        # Save column dtypes
        schema['columns'] = {}
        for col in df.columns:
            if schema_overrides and col in schema_overrides:
                col_schema = dict(schema_overrides[col])
            else:
                col_schema = self._dtype_to_dict(df[col].dtype)
                # Auto-detect file-type columns
                fmt = self._detect_file_format(df[col])
                if fmt is not None:
                    col_schema['type'] = 'file'
                    col_schema['file_format'] = fmt
            schema['columns'][col] = col_schema
        
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

    def _compute_state_hash(self):
        """
        Compute SHA-256 hash of the current data + obj state.

        Used for idempotency detection in update().
        """
        h = hashlib.sha256()
        # Hash data file
        if self.data_path.exists():
            h.update(self.data_path.read_bytes())
        # Hash schema file
        if self.schema_path.exists():
            h.update(self.schema_path.read_bytes())
        # Hash obj file
        if self.obj_path.exists():
            h.update(self.obj_path.read_bytes())
        return h.hexdigest()

    def _validate_commit_description(self, description, meta=None):
        """Validate and resolve the description for the next commit."""
        assert isinstance(description, (str, type(None))), "Description must be a string or None"
        meta = meta or load_json(self.commit_meta_path)
        commit_id = meta["last_commit_count"] + 1
        if description is None:
            return str(commit_id)
        if any(c["description"] == description for c in meta["commits"]):
            raise ValueError(f"Commit description '{description}' already exists")
        if re.fullmatch(r'-?\d+', description.strip()):
            raise ValueError("Commit description cannot be purely numeric")
        return description

    def _ensure_state_files(self):
        """Ensure every component covered by a state hash exists."""
        if not self.data_path.exists():
            self.data_path.touch()
        if not self.schema_path.exists():
            dump_json({}, self.schema_path)
        if not self.obj_path.exists():
            dump_json({}, self.obj_path)
    
    def commit(self, description=None, state_hash=None):
        """
        Create a new commit snapshot of the current data state.

        Parameters
        ----------
        description : str, optional
            Human-readable description of the commit. If None, uses commit ID.
            Must be unique and cannot be purely numeric.
        state_hash : str, optional
            Pre-computed SHA-256 hash of the current data state. If None,
            will be computed from the data files. A provided hash must match
            the current state. Used for idempotency detection in update().

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
        self._ensure_state_files()
        self._validate_state()
        # Load current metadata
        meta = load_json(self.commit_meta_path)
        for key, value in self._default_system_meta().items():
            meta.setdefault(key, value)

        # Generate commit ID and validate description
        commit_id = meta["last_commit_count"] + 1
        description = self._validate_commit_description(description, meta=meta)
        meta["last_commit_count"] = commit_id

        # Get physical filename for this commit using circular naming
        commit_filename = self._get_commit_filename(commit_id)

        # Compute and verify the hash before copying these exact bytes.
        current_state_hash = self._compute_state_hash()
        if state_hash is None:
            state_hash = current_state_hash
        elif state_hash != current_state_hash:
            raise ValueError("Provided state_hash does not match the current Ana state.")

        # Create commit metadata entry
        commit_entry = {
            "id": commit_id,
            "description": description,
            "timestamp": datetime.now().isoformat(),
            "filename": commit_filename,
            "state_hash": state_hash,
            "ana_system_version": meta.get("ana_system_version", ANA_SYSTEM_VERSION),
            "axis": meta.get("axis", self._default_system_meta()["axis"]),
        }
        pending_update_meta = getattr(self, "_pending_update_commit_meta", None)
        if pending_update_meta:
            commit_entry["update_plan"] = pending_update_meta
        meta["commits"].append(commit_entry)
        
        # Create physical snapshots of current data state
        commit_data_path = self.home / f"{commit_filename}.data.csv"
        commit_schema_path = self.home / f"{commit_filename}.data_schema.json"
        commit_obj_path = self.home / f"{commit_filename}.obj.json"

        # Snapshot the exact bytes covered by state_hash. Re-serializing through
        # pandas can change floating-point text without changing logical values.
        shutil.copyfile(self.data_path, commit_data_path)
        shutil.copyfile(self.schema_path, commit_schema_path)
        shutil.copyfile(self.obj_path, commit_obj_path)
        
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
    
    def update(
        self, new_data, key=None, description=None, intent="auto", force=False,
        column_descriptions=None,
        computed_by=None,
        depends_on=None,
        file_format=None,
        schema_changes=None,
    ):
        """
        Update data storage and automatically commit changes.

        This operation is **idempotent**: calling update() with the same data
        will not create additional commits. Use ``force=True`` to override
        this behavior and create a new commit regardless.

        Parameters
        ----------
        new_data : pandas.DataFrame, pandas.Series, dict, or list
            New data to add to the database.
        key : str, optional
            Specifies the column name (for dict/Series) or object key (for list).
        description : str, optional
            Custom description for the automatic commit.
        intent : str, optional, default "auto"
            Update intent. Supported values are ``"auto"``, ``"data_values"``,
            ``"object"``, ``"add_samples"``, ``"drop_samples"``,
            ``"rename_samples"``, and ``"replace_samples"``.
        force : bool, optional, default False
            If True, create a new commit even when the data state is identical
            to the last commit. Useful for recording explicit checkpoints during
            iterative development.
        column_descriptions : dict, optional
            ``{col_name: "short description"}`` for new or existing columns.
            Used by ``describe()`` and for agent discoverability.
        computed_by : dict, optional
            ``{"file": "/path/to/script.ipynb", "note": "optional context"}``.
            If not provided, auto-detected from the calling frame.
        depends_on : dict, optional
            Column/object dependencies. Format::

                {"col_name": {"columns": ["a", "b"], "objects": ["x"]}}

            Can also be a flat dict mapping col_name to a list of mixed deps
            (the method will classify them automatically).
        file_format : dict, optional
            ``{col_name: format_key}`` to override auto-detected file format
            (e.g., ``{"bam_paths": "bam"}``).
        schema_changes : dict, optional
            Explicit dtype changes, for example ``{"score": "Float64"}`` or
            ``{"stage": pd.CategoricalDtype(["2c", "4c"])}``. Existing
            column schemas never change implicitly. Incoming values that can be
            losslessly represented by the existing schema are normalized to it;
            incompatible values require an explicit entry here and a non-empty
            ``description``.

        Returns
        -------
        int
            Commit ID of the new commit, or the last commit ID if the update
            was skipped due to idempotency.

        Raises
        ------
        AnaUpdateError
            If an update violates the sample axis or cannot be represented by
            an existing column schema without an explicit ``schema_changes``
            entry. Validation occurs before persistent state is replaced.
        """
        if isinstance(intent, bool):
            warnings.warn(
                "Passing force as the fourth positional argument is deprecated; "
                "use force=True.",
                DeprecationWarning,
                stacklevel=2,
            )
            force = intent
            intent = "auto"

        if schema_changes is None:
            schema_changes = {}
        if not isinstance(schema_changes, dict):
            raise TypeError("schema_changes must be a dict mapping column names to dtypes.")
        if schema_changes and (
            not isinstance(description, str) or not description.strip()
        ):
            raise AnaUpdateError("Explicit schema_changes require a non-empty description.")

        meta = {
            "column_descriptions": column_descriptions,
            "computed_by": computed_by,
            "depends_on": depends_on,
            "file_format": file_format,
            "description": description,
            "schema_changes": schema_changes,
        }
        plan = self._build_update_plan(new_data, key=key, intent=intent, meta=meta)
        self._validate_plan(plan)
        prepared = self._prepare_plan(plan, schema_changes=schema_changes)
        try:
            if not prepared["changed"] and not force:
                return self._maybe_commit(description, force=False)

            # Validate commit metadata before persistent state is replaced.
            self._validate_commit_description(description)
            touched = self._apply_plan(plan, prepared)
            self._validate_state()

            if touched["columns"]:
                self._inject_column_meta(
                    touched["columns"],
                    column_descriptions=column_descriptions,
                    computed_by=computed_by,
                    depends_on=depends_on,
                    file_format=file_format,
                )

            if touched["objects"]:
                self._inject_column_meta(
                    touched["objects"],
                    column_descriptions=column_descriptions,
                    computed_by=computed_by,
                    depends_on=depends_on,
                )

            if self.verbose:
                self._print_update_summary(plan, touched)
                print("Data updated with schema preservation.")

            self._pending_update_commit_meta = touched.get("commit_meta")
            try:
                return self._maybe_commit(description, force)
            finally:
                self._pending_update_commit_meta = None
        finally:
            self._discard_prepared_plan(prepared)

    def _build_update_plan(self, payload, key=None, intent="auto", meta=None):
        allowed = {
            "auto",
            "data_values",
            "object",
            "add_samples",
            "drop_samples",
            "rename_samples",
            "replace_samples",
        }
        if intent not in allowed:
            raise ValueError(f"Unsupported Ana update intent: {intent!r}")

        current = self.data
        meta = meta or {}

        if intent == "drop_samples":
            return {
                "target": "drop_samples",
                "payload": list(payload),
                "key": key,
                "intent": intent,
                "meta": meta,
            }
        if intent == "rename_samples":
            if not isinstance(payload, dict):
                raise TypeError("rename_samples expects a dict mapping old sample IDs to new sample IDs.")
            return {
                "target": "rename_samples",
                "payload": dict(payload),
                "key": key,
                "intent": intent,
                "meta": meta,
            }
        if intent == "replace_samples":
            df = self._payload_to_dataframe(payload, key=key)
            df = self._promote_sample_id_column(df, "replace_samples")
            return {
                "target": "replace_samples",
                "payload": df,
                "key": key,
                "intent": intent,
                "meta": meta,
            }
        if intent == "object":
            if key is None:
                raise AssertionError("Key must be provided for object updates.")
            if isinstance(payload, pd.DataFrame):
                raise self._unsafe_dataframe_error(payload, key=key)
            return {
                "target": "object",
                "payload": self._json_ready_object(payload),
                "key": key,
                "intent": intent,
                "meta": meta,
            }

        if isinstance(payload, list):
            if key is None:
                raise AssertionError("Key must be provided for list input.")
            return {
                "target": "object",
                "payload": payload,
                "key": key,
                "intent": intent,
                "meta": meta,
            }

        df = self._payload_to_dataframe(payload, key=key)
        if intent == "add_samples":
            df = self._promote_sample_id_column(df, "add_samples")
            return {
                "target": "add_samples",
                "payload": df,
                "key": key,
                "intent": intent,
                "meta": meta,
            }

        if intent == "data_values":
            return {
                "target": "data_values",
                "payload": df,
                "key": key,
                "intent": intent,
                "meta": meta,
            }

        # Conservative auto routing.
        if isinstance(payload, pd.DataFrame):
            if self._dataframe_matches_existing_samples(df, current):
                return {
                    "target": "data_values",
                    "payload": df,
                    "key": key,
                    "intent": intent,
                    "meta": meta,
                }
            raise self._unsafe_dataframe_error(df, key=key)

        if isinstance(payload, (pd.Series, dict)):
            relation = self._index_relation(df.index, current.index)
            if relation == "subset":
                return {
                    "target": "data_values",
                    "payload": df,
                    "key": key,
                    "intent": intent,
                    "meta": meta,
                }
            if relation == "disjoint":
                if key is None:
                    raise AssertionError("Key must be provided for dict or non-sample Series input.")
                if self._looks_like_sample_ids(df.index, current.index):
                    raise AnaUpdateError(
                        "AnaUpdateError: refusing to add new sample IDs with intent='auto'.\n\n"
                        "What Ana saw:\n"
                        f"- payload type: {type(payload).__name__}\n"
                        f"- key: {key}\n"
                        f"- new sample IDs: {self._canonical_labels(df.index)[:10]}\n\n"
                        "Why this is unsafe:\n"
                        "- ana.data is sample-indexed.\n"
                        "- sample-axis changes require an explicit intent.\n\n"
                        "If this is intentional sample metadata:\n"
                        "    ana.update(values, key='column_name', intent='add_samples', description='...')\n\n"
                        "If this is a manifest or summary object:\n"
                        "    ana.update(values, key='manifest_key', intent='object', description='...')"
                    )
                return {
                    "target": "object",
                    "payload": self._json_ready_object(payload),
                    "key": key,
                    "intent": intent,
                    "meta": meta,
                }
            raise AnaUpdateError(
                "AnaUpdateError: refusing partial sample-axis overlap.\n\n"
                "What Ana saw:\n"
                f"- payload type: {type(payload).__name__}\n"
                f"- key: {key}\n"
                f"- matched sample IDs: {len(set(map(str, df.index)) & set(map(str, current.index)))}\n\n"
                "Why this is unsafe:\n"
                "- some labels match ana.data sample IDs and some do not.\n"
                "- Ana cannot infer whether this is sample metadata or an object manifest.\n\n"
                "Use one of these explicit patterns:\n"
                "    ana.update(values, key='column_name', intent='data_values', description='...')\n"
                "    ana.update(values, key='manifest_key', intent='object', description='...')\n"
                "    ana.update(values, key='column_name', intent='add_samples', description='...')"
            )

        raise TypeError(f"Unsupported data type: {type(payload)}")

    def _validate_plan(self, plan):
        target = plan["target"]
        current = self.data

        if target in {"data_values", "add_samples", "replace_samples"}:
            df = plan["payload"]
            self._validate_sample_index(df.index, context=target)

        if target == "data_values":
            missing = self._missing_labels(plan["payload"].index, current.index)
            if missing:
                raise AnaUpdateError(
                    "AnaUpdateError: refusing to add sample IDs with intent='data_values'.\n\n"
                    "What Ana saw:\n"
                    f"- new sample IDs: {missing[:10]}\n\n"
                    "Use explicit sample-axis intent instead:\n"
                    "    ana.update(df, intent='add_samples', description='...')"
                )

        if target == "add_samples":
            self._validate_add_samples_axis(plan["payload"].index, current.index)

        if target == "drop_samples":
            ids = plan["payload"]
            self._validate_sample_labels(ids, context="drop_samples")
            missing = self._missing_labels(ids, current.index)
            if missing:
                raise AnaUpdateError(f"Cannot drop missing sample IDs: {missing}")

        if target == "rename_samples":
            mapping = plan["payload"]
            old_ids = list(mapping.keys())
            new_ids = list(mapping.values())
            self._validate_sample_labels(old_ids, context="rename_samples old IDs")
            self._validate_sample_labels(new_ids, context="rename_samples new IDs")
            missing = self._missing_labels(old_ids, current.index)
            if missing:
                raise AnaUpdateError(f"Cannot rename missing sample IDs: {missing}")
            existing = set(self._canonical_labels(current.index))
            old = set(self._canonical_labels(old_ids))
            conflicts = sorted(set(self._canonical_labels(new_ids)) & (existing - old))
            if conflicts:
                raise AnaUpdateError(f"New sample IDs already exist: {conflicts}")
            proposed = [mapping.get(label, label) for label in current.index]
            self._validate_sample_labels(proposed, context="renamed sample axis")

        if target == "replace_samples":
            description = (plan.get("meta") or {}).get("description")
            if not description:
                raise AnaUpdateError("replace_samples requires a non-empty description.")

        if target == "object" and plan.get("key") is None:
            raise AssertionError("Key must be provided for object updates.")

    def _prepare_plan(self, plan, schema_changes):
        """Build and validate a proposed state without modifying Ana files."""
        target = plan["target"]
        if schema_changes and target not in {"data_values", "add_samples", "replace_samples"}:
            raise AnaUpdateError(
                f"schema_changes cannot be combined with intent={target!r}. "
                "Apply the schema change with a data update."
            )

        if target == "object":
            current_obj = self.obj
            proposed_obj = dict(current_obj)
            proposed_obj[plan["key"]] = plan["payload"]
            return {
                "kind": "object",
                "changed": proposed_obj != current_obj,
                "proposed_obj": proposed_obj,
                "touched_columns": [],
                "touched_objects": [plan["key"]],
            }

        current = self.data
        original_index = list(current.index)
        schema_overrides = self._locked_schema_overrides(current)
        target_dtypes = {
            col: self._dict_to_dtype(entry)
            for col, entry in schema_overrides.items()
        }
        touched_columns = []

        if target in {"data_values", "add_samples", "replace_samples"}:
            incoming = plan["payload"].copy()
            incoming.index = pd.Index(
                self._canonical_labels(incoming.index),
                name=incoming.index.name,
                dtype=pd.StringDtype(),
            )
            current, incoming, target_dtypes, schema_overrides = self._normalize_update_schema(
                current,
                incoming,
                target_dtypes,
                schema_overrides,
                schema_changes,
            )
            touched_columns = list(dict.fromkeys(
                incoming.columns.tolist() + list(schema_changes)
            ))
            if target == "replace_samples":
                proposed = incoming
            else:
                proposed = pd.concat([current, incoming], axis=0, join="outer")
                proposed = proposed.groupby(level=0, observed=True).last()
            proposed = self._cast_dataframe_columns(
                proposed,
                target_dtypes,
                explicit_columns=set(schema_changes),
            )
        elif target == "drop_samples":
            drop_ids = self._canonical_labels(plan["payload"])
            proposed = current.drop(index=drop_ids)
        elif target == "rename_samples":
            mapping = {
                str(old).strip(): str(new).strip()
                for old, new in plan["payload"].items()
            }
            proposed = current.rename(index=mapping)
        else:
            raise ValueError(f"Unknown Ana update target: {target}")

        proposed.index = pd.Index(
            self._canonical_labels(proposed.index),
            name=proposed.index.name,
            dtype=pd.StringDtype(),
        )

        staged = self._stage_dataframe(proposed, schema_overrides=schema_overrides)
        try:
            current_defs = self._current_schema_definitions(current)
            proposed_defs = load_json(staged["schema_path"]).get("columns", {})
            implicit_changes = {
                col: {"from": current_defs[col], "to": proposed_defs[col]}
                for col in set(current_defs) & set(proposed_defs)
                if current_defs[col] != proposed_defs[col] and col not in schema_changes
            }
            if implicit_changes:
                details = ", ".join(sorted(implicit_changes))
                raise AnaUpdateError(
                    "AnaUpdateError: update would change existing column schemas "
                    f"without schema_changes: {details}."
                )

            schema_delta = {}
            for col in sorted(set(current_defs) | set(proposed_defs)):
                old = current_defs.get(col)
                new = proposed_defs.get(col)
                if old != new:
                    if col in schema_changes:
                        reason = "explicit"
                    elif old is None:
                        reason = "added"
                    else:
                        reason = "removed"
                    schema_delta[col] = {
                        "from": old,
                        "to": new,
                        "reason": reason,
                    }

            changed = (
                not self._dataframes_equal(current, staged["data"])
                or current_defs != proposed_defs
            )
            commit_meta = {
                "target": target,
                "axis_delta": self._axis_delta(original_index, list(staged["data"].index)),
            }
            if target == "rename_samples":
                commit_meta["axis_delta"] = self._axis_delta(
                    original_index,
                    list(staged["data"].index),
                    rename_map=plan["payload"],
                )
            if schema_delta:
                commit_meta["schema_delta"] = schema_delta

            staged.update({
                "kind": "data",
                "changed": changed,
                "commit_meta": commit_meta,
                "schema_delta": schema_delta,
                "touched_columns": touched_columns,
                "touched_objects": [],
            })
            return staged
        except Exception:
            self._discard_prepared_plan(staged)
            raise

    def _locked_schema_overrides(self, current):
        persisted = self.get_schema().get("columns", {})
        return {
            col: dict(persisted.get(col, self._dtype_to_dict(current[col].dtype)))
            for col in current.columns
        }

    def _normalize_update_schema(
        self,
        current,
        incoming,
        target_dtypes,
        schema_overrides,
        schema_changes,
    ):
        known_columns = set(current.columns) | set(incoming.columns)
        unknown = sorted(set(schema_changes) - known_columns)
        if unknown:
            raise AnaUpdateError(f"schema_changes refers to unknown columns: {unknown}")

        current = current.copy()
        incoming = incoming.copy()
        for col, spec in schema_changes.items():
            dtype, schema_entry = self._dtype_and_schema_from_change(spec, column=col)
            target_dtypes[col] = dtype
            schema_overrides[col] = schema_entry

        for col, dtype in target_dtypes.items():
            explicit = col in schema_changes
            if explicit and col in current.columns:
                current[col] = self._cast_series_to_dtype(
                    current[col], dtype, column=col, explicit=True
                )
            if col in incoming.columns:
                incoming[col] = self._cast_series_to_dtype(
                    incoming[col], dtype, column=col, explicit=explicit
                )

        return current, incoming, target_dtypes, schema_overrides

    def _dtype_and_schema_from_change(self, spec, column):
        if isinstance(spec, dict):
            schema_entry = dict(spec)
            dtype = self._dict_to_dtype(schema_entry)
        else:
            try:
                dtype = pd.api.types.pandas_dtype(spec)
            except (TypeError, ValueError) as exc:
                raise AnaUpdateError(
                    f"Invalid schema_changes dtype for column {column!r}: {spec!r}"
                ) from exc
            schema_entry = self._dtype_to_dict(dtype)

        if isinstance(dtype, pd.CategoricalDtype) and dtype.categories is None:
            raise AnaUpdateError(
                f"Categorical schema change for column {column!r} must declare categories, "
                "for example pd.CategoricalDtype(['a', 'b'])."
            )
        return dtype, schema_entry

    def _cast_dataframe_columns(self, df, target_dtypes, explicit_columns):
        df = df.copy()
        for col, dtype in target_dtypes.items():
            if col in df.columns:
                df[col] = self._cast_series_to_dtype(
                    df[col],
                    dtype,
                    column=col,
                    explicit=col in explicit_columns,
                )
        return df

    def _cast_series_to_dtype(self, series, dtype, column, explicit):
        try:
            if isinstance(dtype, pd.CategoricalDtype):
                converted = pd.Series(
                    pd.Categorical(series, dtype=dtype),
                    index=series.index,
                    name=series.name,
                )
                lost = series.notna() & converted.isna()
                if lost.any():
                    unknown = series[lost].astype(str).drop_duplicates().head(5).tolist()
                    raise ValueError(f"unknown categories: {unknown}")
                return converted

            if pd.api.types.is_bool_dtype(dtype):
                converted_values = []
                for value in series:
                    if pd.isna(value):
                        converted_values.append(pd.NA)
                    elif isinstance(value, (bool, np.bool_)):
                        converted_values.append(bool(value))
                    elif isinstance(value, (int, float, np.number)) and value in (0, 1):
                        converted_values.append(bool(value))
                    elif str(value).strip().lower() in {"true", "1", "t", "y", "yes"}:
                        converted_values.append(True)
                    elif str(value).strip().lower() in {"false", "0", "f", "n", "no"}:
                        converted_values.append(False)
                    else:
                        raise ValueError(f"invalid boolean value: {value!r}")
                converted = pd.Series(
                    pd.array(converted_values, dtype=pd.BooleanDtype()),
                    index=series.index,
                    name=series.name,
                )
                if not isinstance(dtype, pd.BooleanDtype):
                    converted = converted.astype(dtype)
                return converted

            if pd.api.types.is_integer_dtype(dtype):
                numeric = pd.to_numeric(series, errors="raise")
                non_null = numeric.dropna()
                if len(non_null) and not np.equal(non_null, np.floor(non_null)).all():
                    raise ValueError("fractional values cannot be represented as integers")
                return numeric.astype(dtype)

            if pd.api.types.is_float_dtype(dtype):
                return pd.to_numeric(series, errors="raise").astype(dtype)

            if isinstance(dtype, pd.StringDtype):
                return series.astype(dtype)

            if pd.api.types.is_datetime64_any_dtype(dtype):
                return pd.to_datetime(series, errors="raise").astype(dtype)

            if pd.api.types.is_timedelta64_dtype(dtype):
                return pd.to_timedelta(series, errors="raise").astype(dtype)

            return series.astype(dtype)
        except (TypeError, ValueError, OverflowError) as exc:
            if explicit:
                action = f"requested schema {dtype}"
                hint = "Correct the values or choose another explicit dtype."
            else:
                action = f"existing schema {dtype}"
                hint = (
                    "Ana does not evolve existing schemas automatically. If this change is "
                    f"intentional, pass schema_changes={{'{column}': '<new dtype>'}} and a description."
                )
            examples = series[series.notna()].head(5).tolist()
            missing_ids = series.index[series.isna()].astype(str).tolist()[:5]
            missing_note = f" Missing sample IDs: {missing_ids}." if missing_ids else ""
            raise AnaUpdateError(
                f"AnaUpdateError: column {column!r} cannot be represented by {action}. "
                f"Example incoming values: {examples}.{missing_note} {hint}"
            ) from exc

    def _stage_dataframe(self, df, schema_overrides):
        handle = tempfile.NamedTemporaryFile(
            prefix=".ana-pending-",
            suffix=".csv",
            dir=self.home,
            delete=False,
        )
        handle.close()
        data_path = Path(handle.name)
        schema_path = data_path.parent / f"{data_path.stem}_schema.json"
        try:
            self._save_data_with_schema(
                df,
                data_path,
                schema_overrides=schema_overrides,
            )
            canonical = self._load_data_with_schema(data_path)
            return {
                "data_path": data_path,
                "schema_path": schema_path,
                "data": canonical,
            }
        except Exception:
            data_path.unlink(missing_ok=True)
            schema_path.unlink(missing_ok=True)
            raise

    def _current_schema_definitions(self, current):
        return self._locked_schema_overrides(current)

    @staticmethod
    def _dataframes_equal(left, right):
        try:
            pd.testing.assert_frame_equal(
                left,
                right,
                check_dtype=True,
                check_exact=True,
                check_categorical=True,
            )
            return True
        except AssertionError:
            return False

    def _apply_plan(self, plan, prepared):
        target = plan["target"]
        touched = {
            "columns": prepared["touched_columns"],
            "objects": prepared["touched_objects"],
        }

        if prepared["changed"] and prepared["kind"] == "object":
            self._write_json_atomic(prepared["proposed_obj"], self.obj_path)
        elif prepared["changed"] and prepared["kind"] == "data":
            os.replace(prepared["data_path"], self.data_path)
            os.replace(prepared["schema_path"], self.schema_path)

        if prepared.get("commit_meta"):
            touched["commit_meta"] = prepared["commit_meta"]

        if prepared["changed"] and target == "drop_samples":
            self._sync_sample_objects(drop_ids=plan["payload"])
        elif prepared["changed"] and target == "rename_samples":
            self._sync_sample_objects(rename_map=plan["payload"])
        elif prepared["changed"] and target == "replace_samples":
            new_ids = self._canonical_labels(prepared["data"].index)
            dropped = prepared["commit_meta"]["axis_delta"]["dropped"]
            self._sync_sample_objects(valid_ids=new_ids, drop_ids=dropped)

        return touched

    def _write_json_atomic(self, obj, path):
        handle = tempfile.NamedTemporaryFile(
            prefix=f".{Path(path).stem}-pending-",
            suffix=".json",
            dir=Path(path).parent,
            delete=False,
        )
        handle.close()
        pending = Path(handle.name)
        try:
            dump_json(obj, pending)
            os.replace(pending, path)
        finally:
            pending.unlink(missing_ok=True)

    @staticmethod
    def _discard_prepared_plan(prepared):
        for key in ("data_path", "schema_path"):
            path = prepared.get(key)
            if path is not None:
                Path(path).unlink(missing_ok=True)

    def _payload_to_dataframe(self, payload, key=None):
        if isinstance(payload, pd.DataFrame):
            return payload.copy()
        if isinstance(payload, pd.Series):
            series = payload.copy()
            if key is not None:
                series = series.rename(key)
            return series.to_frame()
        if isinstance(payload, dict):
            assert key is not None, "Key must be provided for dict input."
            return pd.DataFrame({key: pd.Series(payload)})
        raise TypeError(f"Unsupported data type: {type(payload)}")

    def _json_ready_object(self, payload):
        if isinstance(payload, pd.Series):
            return payload.to_dict()
        if isinstance(payload, dict):
            return payload
        if isinstance(payload, list):
            return payload
        return payload

    def _promote_sample_id_column(self, df, context):
        if "sample_id" in df.columns and isinstance(df.index, pd.RangeIndex):
            warnings.warn(
                f"{context}: promoting 'sample_id' column to the Ana sample axis.",
                UserWarning,
                stacklevel=3,
            )
            df = df.set_index("sample_id")
        return df

    def _validate_state(self):
        df = self.data
        if len(df.index) == 0:
            return
        self._validate_sample_index(df.index, context="ana.data")

    def _validate_sample_index(self, index, context):
        if len(index) == 0:
            return
        if isinstance(index, pd.RangeIndex):
            raise AnaUpdateError(
                f"AnaUpdateError: refusing RangeIndex as {context} sample IDs.\n\n"
                "ana.data is sample-indexed. Set a real sample_id index or use "
                "intent='object' for manifests/results."
            )
        self._validate_sample_labels(index, context=context)

    def _validate_sample_labels(self, labels, context):
        seen = {}
        for label in labels:
            if pd.isna(label):
                raise AnaUpdateError(f"{context}: sample IDs cannot be null.")
            text = str(label)
            stripped = text.strip()
            if stripped == "":
                raise AnaUpdateError(f"{context}: sample IDs cannot be empty or whitespace-only.")
            if text != stripped:
                raise AnaUpdateError(
                    f"{context}: sample ID {text!r} has leading/trailing whitespace. "
                    "Canonicalize IDs before calling ana.update()."
                )
            if stripped in seen:
                raise AnaUpdateError(
                    f"{context}: sample IDs collide after CSV/string canonicalization: "
                    f"{seen[stripped]!r} and {label!r}."
                )
            seen[stripped] = label

    def _canonical_labels(self, labels):
        return [str(label).strip() for label in labels]

    def _missing_labels(self, labels, existing):
        existing_set = set(self._canonical_labels(existing))
        return [label for label in self._canonical_labels(labels) if label not in existing_set]

    def _validate_add_samples_axis(self, incoming, existing):
        existing_by_canon = {}
        for label in existing:
            existing_by_canon[str(label).strip()] = label
        for label in incoming:
            canon = str(label).strip()
            if canon in existing_by_canon and label != existing_by_canon[canon]:
                raise AnaUpdateError(
                    "add_samples: incoming sample ID collides with an existing ID after "
                    f"CSV/string canonicalization: {label!r} vs {existing_by_canon[canon]!r}."
                )

    def _index_relation(self, incoming, existing):
        if len(incoming) == 0:
            return "subset"
        incoming_set = set(self._canonical_labels(incoming))
        existing_set = set(self._canonical_labels(existing))
        if incoming_set <= existing_set:
            return "subset"
        if incoming_set.isdisjoint(existing_set):
            return "disjoint"
        return "partial"

    def _dataframe_matches_existing_samples(self, df, current):
        if isinstance(df.index, pd.RangeIndex):
            return False
        return self._index_relation(df.index, current.index) == "subset"

    def _looks_like_sample_ids(self, labels, existing):
        labels = self._canonical_labels(labels)
        existing = self._canonical_labels(existing)
        if not labels:
            return False
        if all(re.match(r"^(sample|cell|s)[_-]?\w*\d+$", label, re.IGNORECASE) for label in labels):
            return True
        existing_prefixes = {
            re.sub(r"\d+$", "", label)
            for label in existing
            if re.search(r"\d+$", label)
        }
        if existing_prefixes and all(re.sub(r"\d+$", "", label) in existing_prefixes for label in labels):
            return True
        return False

    def _axis_delta(self, old_ids, new_ids, rename_map=None):
        old = self._canonical_labels(old_ids)
        new = self._canonical_labels(new_ids)
        old_set = set(old)
        new_set = set(new)
        delta = {
            "kept": sorted(old_set & new_set),
            "added": sorted(new_set - old_set),
            "dropped": sorted(old_set - new_set),
        }
        if rename_map:
            delta["renamed"] = {str(k).strip(): str(v).strip() for k, v in rename_map.items()}
        return delta

    def _unsafe_dataframe_error(self, df, key=None):
        index_desc = type(df.index).__name__
        if isinstance(df.index, pd.RangeIndex):
            index_desc = f"RangeIndex({df.index.start}..{df.index.stop - 1})"
        return AnaUpdateError(
            "AnaUpdateError: refusing to merge DataFrame into ana.data.\n\n"
            "What Ana saw:\n"
            "- payload type: DataFrame\n"
            f"- shape: {df.shape}\n"
            f"- index: {index_desc}\n"
            f"- key: {key}\n\n"
            "Why this is unsafe:\n"
            "- ana.data is sample-indexed.\n"
            "- DataFrame key is not a named-table registry.\n"
            "- this would add non-sample row IDs to ana.data.\n\n"
            "If this is sample-level metadata:\n"
            "    df = df.set_index('sample_id')\n"
            "    ana.update(df, intent='add_samples', description='...')\n\n"
            "If this is an analysis result table:\n"
            "    out = h.ddir / '...' / 'top_residuals.tsv'\n"
            "    df.to_csv(out, sep='\\t', index=False)\n"
            "    ana.update({'top_residuals_tsv': str(out)}, key='analysis_manifest', description='...')"
        )

    def _sync_sample_objects(self, valid_ids=None, drop_ids=None, rename_map=None):
        obj = self.obj
        if not obj:
            return

        current_ids = set(self._canonical_labels(self.data.index))
        valid = set(self._canonical_labels(valid_ids)) if valid_ids is not None else current_ids
        drop = set(self._canonical_labels(drop_ids or []))
        rename = {str(k).strip(): v for k, v in (rename_map or {}).items()}
        changed = False

        for key, value in list(obj.items()):
            if not isinstance(value, list):
                if self._object_contains_labels(value, drop or set(rename)):
                    warnings.warn(
                        f"Object '{key}' contains sample IDs affected by the axis update; "
                        "Ana did not rewrite this object automatically.",
                        UserWarning,
                        stacklevel=3,
                    )
                continue

            labels = self._canonical_labels(value)
            universe = current_ids | valid | drop | set(rename) | set(map(str, rename.values()))
            if not set(labels) <= universe:
                if self._object_contains_labels(value, drop or set(rename)):
                    warnings.warn(
                        f"Object '{key}' contains sample IDs affected by the axis update; "
                        "Ana did not rewrite it because it is not a pure sample set.",
                        UserWarning,
                        stacklevel=3,
                    )
                continue

            new_value = []
            for item in value:
                canon = str(item).strip()
                if canon in drop:
                    continue
                item = rename.get(canon, item)
                if str(item).strip() in valid:
                    new_value.append(item)
            if new_value != value:
                obj[key] = new_value
                changed = True

        if changed:
            dump_json(obj, self.obj_path)

    def _object_contains_labels(self, value, labels):
        if not labels:
            return False
        if isinstance(value, dict):
            vals = list(value.keys()) + list(value.values())
        elif isinstance(value, (list, tuple, set)):
            vals = list(value)
        else:
            vals = [value]
        return bool(set(self._canonical_labels(vals)) & set(labels))

    def _print_update_summary(self, plan, touched):
        print(f"Ana update plan: {plan['target']}")
        if touched["columns"]:
            print(f"columns updated: {touched['columns']}")
        if touched["objects"]:
            print(f"objects updated: {touched['objects']}")

    def _inject_column_meta(
        self,
        columns,
        column_descriptions=None,
        computed_by=None,
        depends_on=None,
        file_format=None,
    ):
        """
        Update column metadata for the given columns.

        Parameters
        ----------
        columns : list
            Column names (or object keys) to annotate.
        column_descriptions : dict, optional
            ``{col: "description"}``.
        computed_by : dict, optional
            ``{"file": "...", "note": "..."}``. If None, auto-detected.
        depends_on : dict, optional
            ``{col: {"columns": [...], "objects": [...]}}`` or
            ``{col: [dep1, dep2, ...]}`` (flat list).
        file_format : dict, optional
            ``{col: format_key}`` for explicit file format override.
        """
        meta = self._load_column_meta()
        commit_id = None
        commit_meta_raw = load_json(self.commit_meta_path)
        if commit_meta_raw["commits"]:
            commit_id = commit_meta_raw["commits"][-1]["id"]

        # Auto-detect computed_by if not provided
        if computed_by is None:
            computed_by = self._detect_computed_by()

        for col in columns:
            if col not in meta:
                meta[col] = {}

            entry = meta[col]

            # Description
            if column_descriptions and col in column_descriptions:
                entry["description"] = column_descriptions[col]

            # computed_by (only set if not already present, or always override)
            entry["computed_by"] = computed_by

            # depends_on
            if depends_on and col in depends_on:
                dep = depends_on[col]
                if isinstance(dep, dict):
                    entry["depends_on"] = dep
                elif isinstance(dep, list):
                    # Classify deps: columns that exist in data vs. objects
                    entry["depends_on"] = {"columns": [], "objects": []}
                    for d in dep:
                        # Heuristic: if d looks like a known column or object key
                        # We keep it simple and put everything in "columns" for now
                        entry["depends_on"]["columns"].append(d)

            # file_format
            if file_format and col in file_format:
                entry["file_format"] = file_format[col]

            # commit_id
            if commit_id is not None:
                entry.setdefault("commit_id", commit_id)

        self._save_column_meta(meta)

    def _maybe_commit(self, description, force):
        """
        Commit only if the current state differs from the last commit,
        unless ``force=True``.

        Returns
        -------
        int
            New commit ID, or the last commit ID if skipped.
        """
        meta = load_json(self.commit_meta_path)
        current_hash = self._compute_state_hash()

        # Get last commit's hash (if any)
        last_hash = None
        last_commit_id = None
        if meta["commits"]:
            last_commit = meta["commits"][-1]
            last_hash = last_commit.get("state_hash")
            last_commit_id = last_commit["id"]

        # Idempotency check
        if not force and last_hash is not None and current_hash == last_hash:
            if self.verbose:
                print(f"Data unchanged, skipping commit (last commit: {last_commit_id})")
            return last_commit_id

        # Proceed with commit
        return self.commit(description=description, state_hash=current_hash)
    
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
        commit_schema_path = self.home / f"{commit_filename}.data_schema.json"
        commit_obj_path = self.home / f"{commit_filename}.obj.json"
        
        if not commit_data_path.exists():
            raise ValueError(f"Commit data file for '{commit_repr}' not found on disk")
        if not commit_schema_path.exists():
            raise ValueError(f"Commit schema file for '{commit_repr}' not found on disk")
        if not commit_obj_path.exists():
            raise ValueError(f"Commit object file for '{commit_repr}' not found on disk")

        revert_description = f"Reverted to commit {commit_id}"
        self._validate_commit_description(revert_description, meta=meta)

        # Restore the exact bytes identified by the selected commit.
        shutil.copyfile(commit_data_path, self.data_path)
        shutil.copyfile(commit_schema_path, self.schema_path)
        shutil.copyfile(commit_obj_path, self.obj_path)
        
        # Create new commit to record the revert action
        new_commit_id = self.commit(revert_description)
        
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

    # ---- Column metadata helpers ----

    @property
    def _column_meta_path(self):
        """Path to the column metadata JSON file."""
        return self.home / "column_meta.json"

    def _load_column_meta(self) -> dict:
        """Load column metadata, returning empty dict if missing."""
        path = self._column_meta_path
        if path.exists():
            return load_json(path)
        return {}

    def _save_column_meta(self, meta: dict):
        """Save column metadata."""
        dump_json(meta, self._column_meta_path)

    def _detect_computed_by(self) -> dict:
        """
        Attempt to detect the file that called ``update()``.

        Priority:
        1. VS Code Jupyter injected ``__vsc_ipynb_file__``
        2. ``ipynbname`` package (Jupyter Lab/Notebook)
        3. ``inspect.stack()`` — first frame outside hic_basic package
        4. Fallback to ``{"file": "<unknown>"}``

        Returns
        -------
        dict
            ``{"file": "<path>"}`` or ``{"file": "<unknown>"}`` with optional warning.
        """
        # 1. VS Code Jupyter
        import builtins
        vsc_file = getattr(builtins, "__vsc_ipynb_file__", None)
        if vsc_file:
            return {"file": str(vsc_file)}

        # 2. ipynbname package
        try:
            import ipynbname
            nb_path = ipynbname.path()
            if Path(nb_path).exists():
                return {"file": str(nb_path)}
        except (ImportError, FileNotFoundError, KeyError, ValueError):
            pass

        # 3. inspect.stack() — find first frame not in hic_basic
        try:
            for frame_info in inspect.stack()[2:]:
                fname = frame_info.filename
                if "hic_basic" not in fname and fname != "<string>":
                    abs_path = os.path.abspath(fname)
                    if Path(abs_path).exists():
                        return {"file": abs_path}
                    # Frame exists but file not found
                    return {"file": fname, "note": "file not found on disk"}
        except Exception:
            pass

        # 4. Fallback
        return {"file": "<unknown>", "note": "could not auto-detect caller; please declare computed_by manually"}

    def describe(self):
        """
        Print a compact, agent-friendly summary of columns and objects.

        Output includes column names, their schema type, and short
        descriptions (if available). Also lists stored objects.

        Example output::

            Columns (3):
              bam_paths        file         Aligned BAM paths per sample
              norm_contact     standard     KR-normalized contact freq
              raw_contact      standard     Raw contact counts

            Objects (2):
              peak_files       (no description)
              sample_groups    Group assignments for 3 conditions
        """
        schema = self.get_schema()
        col_meta = self._load_column_meta()
        cols = schema.get("columns", {})
        col_names = schema.get("columns_list", [])

        print(f"Columns ({len(col_names)}):")
        for col in col_names:
            col_schema = cols.get(col, {})
            col_type = col_schema.get("type", "standard")
            # Use user-provided description if available, else fallback
            meta_entry = col_meta.get(col, {})
            desc = meta_entry.get("description", "(no description)")
            print(f"  {col:<20s} {col_type:<14s} {desc}")

        # Objects
        obj = self.obj
        if obj:
            print(f"\nObjects ({len(obj)}):")
            for key in sorted(obj):
                desc = col_meta.get(key, {}).get("description", "(no description)")
                print(f"  {key:<20s} {desc}")

    def column_meta(self, key: str) -> dict:
        """
        Get full metadata for a column or object key.

        Parameters
        ----------
        key : str
            Column name or object key.

        Returns
        -------
        dict
            Metadata dict (description, computed_by, depends_on,
            file_structure, commit_id). Empty dict if not found.
        """
        return self._load_column_meta().get(key, {})

    def file_format(self, key: str) -> dict:
        """
        Look up the file format definition for a file-type column.

        Searches instance-level overrides first, then falls back to the
        built-in ``hic_basic.data.file_formats`` registry.

        Parameters
        ----------
        key : str
            File format key (e.g., ``"cool"``, ``"bam"``) or column name
            whose ``file_format`` schema value will be used for lookup.

        Returns
        -------
        dict or None
            Format definition, or None if not found.
        """
        from .data.file_formats import FILE_FORMATS as _BUILTIN_FORMATS

        # First check if key is a column name → look up its file_format
        schema = self.get_schema()
        cols = schema.get("columns", {})
        if key in cols:
            col_schema = cols[key]
            fmt_key = col_schema.get("file_format", key)
        else:
            fmt_key = key

        # Check instance-level overrides in column_meta
        col_meta = self._load_column_meta()
        instance_formats = col_meta.get("file_formats", {})
        if fmt_key in instance_formats:
            return instance_formats[fmt_key]

        # Fall back to built-in
        return _BUILTIN_FORMATS.get(fmt_key)
