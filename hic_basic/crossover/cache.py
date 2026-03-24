"""
Intermediate file caching and flag management for fault-tolerant pipeline.

This module provides caching infrastructure for the crossover detection pipeline,
enabling:
- Intermediate result persistence between pipeline steps
- Flag files to mark completed steps
- Resume capability after interruption
- Shared memory objects for efficient parallel processing
"""

import hashlib
import json
import logging
import os
import shutil
import tempfile
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class CacheEntry:
    """
    Metadata for a cached file.

    Attributes:
        file_path: Path to the cached file
        checksum: MD5 checksum of input parameters
        created_at: Creation timestamp
        step_name: Name of the pipeline step
        params: Input parameters that generated this cache
    """
    file_path: str
    checksum: str
    created_at: str
    step_name: str
    params: Dict


class CacheManager:
    """
    Manages intermediate file caching for the crossover detection pipeline.

    Features:
        - Content-addressable storage using MD5 checksums
        - Flag files to mark completed steps
        - Automatic cleanup of stale caches
        - Memory-mapped file support for large arrays

    Attributes:
        cache_dir: Base directory for cache storage
        use_parquet: Whether to use parquet format for DataFrames
        compression: Compression type for saved files

    Example:
        >>> cache = CacheManager("/tmp/crossover_cache")
        >>> key = cache.compute_key({"sam_file": "sample.bam", "window_size": 10})
        >>> cache.save(df, key, "marked_alleles")
        >>> if cache.exists(key, "marked_alleles"):
        ...     df = cache.load(key, "marked_alleles")
    """

    def __init__(
        self,
        cache_dir: str,
        use_parquet: bool = True,
        compression: str = "snappy"
    ):
        """
        Initialize cache manager.

        Args:
            cache_dir: Base directory for cache storage
            use_parquet: Whether to use parquet format for DataFrames
            compression: Compression type for saved files
        """
        self.cache_dir = Path(cache_dir)
        self.use_parquet = use_parquet
        self.compression = compression

        self.cache_dir.mkdir(parents=True, exist_ok=True)

        self._data_dir = self.cache_dir / "data"
        self._flag_dir = self.cache_dir / "flags"
        self._meta_dir = self.cache_dir / "metadata"

        for d in [self._data_dir, self._flag_dir, self._meta_dir]:
            d.mkdir(parents=True, exist_ok=True)

        self._index_file = self._meta_dir / "cache_index.json"
        self._index: Dict[str, CacheEntry] = self._load_index()

    def _load_index(self) -> Dict[str, CacheEntry]:
        """Load cache index from disk."""
        if self._index_file.exists():
            with open(self._index_file, 'r') as f:
                data = json.load(f)
                return {k: CacheEntry(**v) for k, v in data.items()}
        return {}

    def _save_index(self) -> None:
        """Save cache index to disk."""
        with open(self._index_file, 'w') as f:
            data = {k: self._entry_to_dict(v) for k, v in self._index.items()}
            json.dump(data, f, indent=2)

    def _entry_to_dict(self, entry: CacheEntry) -> Dict:
        """Convert CacheEntry to dictionary for JSON serialization."""
        return {
            "file_path": str(entry.file_path),
            "checksum": entry.checksum,
            "created_at": entry.created_at,
            "step_name": entry.step_name,
            "params": entry.params,
        }

    def compute_key(self, params: Dict[str, Any]) -> str:
        """
        Compute a unique key based on input parameters.

        Args:
            params: Dictionary of input parameters

        Returns:
            MD5 checksum as hex string
        """
        param_str = json.dumps(params, sort_keys=True, default=str)
        return hashlib.md5(param_str.encode()).hexdigest()

    def _get_data_path(
        self,
        key: str,
        step_name: str,
        extension: Optional[str] = None
    ) -> Path:
        """Get the file path for cached data."""
        if extension is None:
            extension = "parquet" if self.use_parquet else "tsv.gz"
        return self._data_dir / f"{key}.{step_name}.{extension}"

    def _get_flag_path(self, key: str, step_name: str) -> Path:
        """Get the flag file path for a completed step."""
        return self._flag_dir / f"{key}.{step_name}.done"

    def _get_meta_path(self, key: str, step_name: str) -> Path:
        """Get the metadata file path."""
        return self._meta_dir / f"{key}.{step_name}.json"

    def exists(self, key: str, step_name: str) -> bool:
        """
        Check if cached data exists for a given key and step.

        Args:
            key: Cache key (MD5 checksum)
            step_name: Name of the pipeline step

        Returns:
            True if cache exists and is valid
        """
        flag_path = self._get_flag_path(key, step_name)
        return flag_path.exists()

    def is_flagged(self, key: str, step_name: str) -> bool:
        """
        Check if a step is marked as completed.

        Args:
            key: Cache key
            step_name: Pipeline step name

        Returns:
            True if flag file exists
        """
        return self.exists(key, step_name)

    def set_flag(self, key: str, step_name: str, params: Optional[Dict] = None) -> None:
        """
        Mark a step as completed.

        Args:
            key: Cache key
            step_name: Pipeline step name
            params: Optional parameters to save with the flag
        """
        flag_path = self._get_flag_path(key, step_name)
        flag_path.touch()

        if params is not None:
            meta_path = self._get_meta_path(key, step_name)
            meta_data = {
                "key": key,
                "step_name": step_name,
                "params": params,
                "created_at": datetime.now().isoformat(),
            }
            with open(meta_path, 'w') as f:
                json.dump(meta_data, f, indent=2)

    def clear_flag(self, key: str, step_name: str) -> None:
        """
        Clear the completion flag for a step.

        Args:
            key: Cache key
            step_name: Pipeline step name
        """
        flag_path = self._get_flag_path(key, step_name)
        if flag_path.exists():
            flag_path.unlink()

        meta_path = self._get_meta_path(key, step_name)
        if meta_path.exists():
            meta_path.unlink()

    def save(
        self,
        data: Union[pd.DataFrame, np.ndarray, Dict],
        key: str,
        step_name: str,
        params: Optional[Dict] = None
    ) -> str:
        """
        Save data to cache.

        Args:
            data: Data to cache (DataFrame, array, or dict)
            key: Cache key
            step_name: Pipeline step name
            params: Optional parameters for metadata

        Returns:
            Path to the cached file
        """
        data_path = self._get_data_path(key, step_name)

        if isinstance(data, pd.DataFrame):
            if self.use_parquet:
                data.to_parquet(data_path, index=False, compression=self.compression)
            else:
                data.to_csv(data_path, sep='\t', index=False, compression='gzip')
        elif isinstance(data, np.ndarray):
            np.save(data_path.with_suffix('.npy'), data)
            data_path = data_path.with_suffix('.npy')
        elif isinstance(data, dict):
            with open(data_path.with_suffix('.json'), 'w') as f:
                json.dump(data, f, indent=2)
            data_path = data_path.with_suffix('.json')
        else:
            raise ValueError(f"Unsupported data type: {type(data)}")

        entry = CacheEntry(
            file_path=str(data_path),
            checksum=key,
            created_at=datetime.now().isoformat(),
            step_name=step_name,
            params=params or {}
        )

        cache_key = f"{key}.{step_name}"
        self._index[cache_key] = entry
        self._save_index()

        logger.debug(f"Cached {step_name} to {data_path}")
        return str(data_path)

    def load(
        self,
        key: str,
        step_name: str
    ) -> Union[pd.DataFrame, np.ndarray, Dict, None]:
        """
        Load data from cache.

        Args:
            key: Cache key
            step_name: Pipeline step name

        Returns:
            Loaded data or None if not found
        """
        cache_key = f"{key}.{step_name}"
        if cache_key not in self._index:
            return None

        entry = self._index[cache_key]
        data_path = Path(entry.file_path)

        if not data_path.exists():
            logger.warning(f"Cached file missing: {data_path}")
            return None

        suffix = data_path.suffix.lower()

        if suffix == '.parquet':
            return pd.read_parquet(data_path)
        elif suffix == '.gz':
            return pd.read_csv(data_path, sep='\t')
        elif suffix == '.npy':
            return np.load(data_path)
        elif suffix == '.json':
            with open(data_path, 'r') as f:
                return json.load(f)
        else:
            logger.warning(f"Unknown file suffix: {suffix}")
            return None

    def get_or_compute(
        self,
        key: str,
        step_name: str,
        compute_fn,
        params: Optional[Dict] = None
    ):
        """
        Get cached data or compute it if not available.

        Args:
            key: Cache key
            step_name: Pipeline step name
            compute_fn: Function to compute data if not cached
            params: Parameters to pass to compute_fn and save as metadata

        Returns:
            Cached or computed data
        """
        if self.exists(key, step_name):
            logger.debug(f"Cache hit for {step_name}")
            data = self.load(key, step_name)
            if data is not None:
                return data

        logger.debug(f"Cache miss for {step_name}, computing...")
        data = compute_fn()
        self.save(data, key, step_name, params)
        self.set_flag(key, step_name, params)
        return data

    def clear(self, key: Optional[str] = None, step_name: Optional[str] = None) -> int:
        """
        Clear cached data.

        Args:
            key: Optional cache key to clear (clear all if None)
            step_name: Optional step name to filter (clear all steps if None)

        Returns:
            Number of entries cleared
        """
        cleared = 0
        keys_to_remove = []

        for cache_key, entry in self._index.items():
            should_remove = True

            if key is not None and not cache_key.startswith(key):
                should_remove = False

            if step_name is not None and step_name not in cache_key:
                should_remove = False

            if should_remove:
                data_path = Path(entry.file_path)
                if data_path.exists():
                    data_path.unlink()

                flag_path = self._get_flag_path(entry.checksum, entry.step_name)
                if flag_path.exists():
                    flag_path.unlink()

                meta_path = self._get_meta_path(entry.checksum, entry.step_name)
                if meta_path.exists():
                    meta_path.unlink()

                keys_to_remove.append(cache_key)
                cleared += 1

        for cache_key in keys_to_remove:
            del self._index[cache_key]

        self._save_index()
        logger.info(f"Cleared {cleared} cache entries")
        return cleared

    def get_stats(self) -> Dict:
        """
        Get cache statistics.

        Returns:
            Dictionary with cache statistics
        """
        total_size = 0
        file_count = 0

        for entry in self._index.values():
            data_path = Path(entry.file_path)
            if data_path.exists():
                total_size += data_path.stat().st_size
                file_count += 1

        return {
            "total_entries": len(self._index),
            "total_files": file_count,
            "total_size_bytes": total_size,
            "total_size_mb": total_size / (1024 * 1024),
        }


class SharedMemoryManager:
    """
    Manages shared memory objects for parallel processing.

    This class provides memory-efficient sharing of large arrays
    between processes using memory-mapped files.

    Attributes:
        temp_dir: Directory for temporary memory-mapped files
        arrays: Dictionary of shared memory arrays

    Example:
        >>> smm = SharedMemoryManager()
        >>> arr = np.random.rand(1000000)
        >>> shared_key = smm.share(arr, "snp_positions")
        >>> retrieved_arr = smm.get(shared_key)
    """

    def __init__(self, temp_dir: Optional[str] = None):
        """
        Initialize shared memory manager.

        Args:
            temp_dir: Directory for temporary files (uses system temp if None)
        """
        if temp_dir is None:
            self.temp_dir = Path(tempfile.gettempdir()) / "crossover_shm"
        else:
            self.temp_dir = Path(temp_dir)

        self.temp_dir.mkdir(parents=True, exist_ok=True)
        self.arrays: Dict[str, Dict] = {}
        self._index_file = self.temp_dir / "shm_index.json"
        self._load_index()

    def _load_index(self) -> None:
        """Load shared memory index."""
        if self._index_file.exists():
            with open(self._index_file, 'r') as f:
                self.arrays = json.load(f)

    def _save_index(self) -> None:
        """Save shared memory index."""
        with open(self._index_file, 'w') as f:
            json.dump(self.arrays, f, indent=2)

    def share(
        self,
        array: np.ndarray,
        name: str,
        sample_id: Optional[str] = None
    ) -> str:
        """
        Share an array via memory-mapped file.

        Args:
            array: NumPy array to share
            name: Name for the shared array
            sample_id: Optional sample ID for namespacing

        Returns:
            Key to retrieve the shared array
        """
        if sample_id:
            key = f"{sample_id}.{name}"
        else:
            key = name

        mmap_path = self.temp_dir / f"{key}.mmap"

        with open(mmap_path, 'wb') as f:
            f.write(array.tobytes())

        self.arrays[key] = {
            "path": str(mmap_path),
            "shape": array.shape,
            "dtype": str(array.dtype),
            "size_bytes": array.nbytes,
        }
        self._save_index()

        logger.debug(f"Shared array {key} ({array.nbytes / 1e6:.2f} MB)")
        return key

    def get(self, key: str) -> Optional[np.ndarray]:
        """
        Retrieve a shared array.

        Args:
            key: Array key

        Returns:
            NumPy array or None if not found
        """
        if key not in self.arrays:
            return None

        info = self.arrays[key]
        mmap_path = Path(info["path"])

        if not mmap_path.exists():
            logger.warning(f"Shared memory file missing: {mmap_path}")
            return None

        dtype = np.dtype(info["dtype"])
        shape = tuple(info["shape"])

        with open(mmap_path, 'rb') as f:
            array = np.frombuffer(f.read(), dtype=dtype).reshape(shape)

        return array.copy()

    def share_dataframe(
        self,
        df: pd.DataFrame,
        name: str,
        sample_id: Optional[str] = None
    ) -> str:
        """
        Share a DataFrame via parquet memory-mapping.

        Args:
            df: DataFrame to share
            name: Name for the shared DataFrame
            sample_id: Optional sample ID for namespacing

        Returns:
            Key to retrieve the shared DataFrame
        """
        if sample_id:
            key = f"{sample_id}.{name}"
        else:
            key = name

        parquet_path = self.temp_dir / f"{key}.parquet"
        df.to_parquet(parquet_path, index=False)

        self.arrays[key] = {
            "path": str(parquet_path),
            "type": "dataframe",
            "columns": list(df.columns),
            "size_bytes": parquet_path.stat().st_size,
        }
        self._save_index()

        logger.debug(f"Shared DataFrame {key} ({parquet_path.stat().st_size / 1e6:.2f} MB)")
        return key

    def get_dataframe(self, key: str) -> Optional[pd.DataFrame]:
        """
        Retrieve a shared DataFrame.

        Args:
            key: DataFrame key

        Returns:
            DataFrame or None if not found
        """
        if key not in self.arrays:
            return None

        info = self.arrays[key]
        if info.get("type") != "dataframe":
            logger.warning(f"Key {key} is not a DataFrame")
            return None

        parquet_path = Path(info["path"])
        if not parquet_path.exists():
            return None

        return pd.read_parquet(parquet_path)

    def clear(self, key: Optional[str] = None) -> int:
        """
        Clear shared memory.

        Args:
            key: Optional key to clear (clear all if None)

        Returns:
            Number of entries cleared
        """
        cleared = 0
        keys_to_remove = []

        for k, info in self.arrays.items():
            if key is None or k == key or (key and k.startswith(key)):
                path = Path(info["path"])
                if path.exists():
                    path.unlink()
                keys_to_remove.append(k)
                cleared += 1

        for k in keys_to_remove:
            del self.arrays[k]

        self._save_index()
        logger.debug(f"Cleared {cleared} shared memory entries")
        return cleared

    def cleanup(self) -> None:
        """Clean up all shared memory files."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)
        self.arrays = {}
        logger.info("Cleaned up all shared memory")
