import h5py
import numpy as np
import pandas as pd
import json
import os
import scipy.sparse as sp
from pathlib import Path
from typing import Dict, Any, Optional, Union, List
import warnings


class H5ADConverter:
    """
    A converter to handle compatibility issues between different anndata versions.
    Converts h5ad files to/from generic file formats (CSV, JSON, NPY, etc.).
    
    Attributes
    ----------
    supported_layers : list
        List of supported data layers that can be converted
    """
    
    def __init__(self):
        self.supported_layers = ['X', 'layers', 'obsm', 'varm', 'obs', 'var', 'uns']
        self._dump_registry = {}
        self._load_registry = {}
        
        # Register default handlers
        self._register_default_handlers()
    
    def _register_default_handlers(self):
        """Register default dump and load handlers for all supported layers."""
        # X matrix handler
        self.register_dump_handler('X', self._dump_X)
        self.register_load_handler('X', self._load_X)
        
        # Layers handler
        self.register_dump_handler('layers', self._dump_layers)
        self.register_load_handler('layers', self._load_layers)
        
        # obs metadata handler
        self.register_dump_handler('obs', self._dump_obs)
        self.register_load_handler('obs', self._load_obs)
        
        # var metadata handler
        self.register_dump_handler('var', self._dump_var)
        self.register_load_handler('var', self._load_var)
        
        # obsm embeddings handler
        self.register_dump_handler('obsm', self._dump_obsm)
        self.register_load_handler('obsm', self._load_obsm)
        
        # varm embeddings handler
        self.register_dump_handler('varm', self._dump_varm)
        self.register_load_handler('varm', self._load_varm)
        
        # uns unstructured data handler
        self.register_dump_handler('uns', self._dump_uns)
        self.register_load_handler('uns', self._load_uns)
    
    def register_dump_handler(self, layer_name: str, handler_func):
        """
        Register a custom dump handler for a specific layer.
        
        Parameters
        ----------
        layer_name : str
            Name of the layer to register handler for
        handler_func : callable
            Function that handles dumping this layer
        """
        self._dump_registry[layer_name] = handler_func
    
    def register_load_handler(self, layer_name: str, handler_func):
        """
        Register a custom load handler for a specific layer.
        
        Parameters
        ----------
        layer_name : str
            Name of the layer to register handler for
        handler_func : callable
            Function that handles loading this layer
        """
        self._load_registry[layer_name] = handler_func
    
    def dump(self, h5ad_file: str, output_dir: str, layers: Optional[List[str]] = None) -> Dict[str, Any]:
        """
        Convert h5ad file to generic file formats.
        
        Parameters
        ----------
        h5ad_file : str
            Path to input h5ad file
        output_dir : str
            Directory to save converted files
        layers : list, optional
            List of layers to convert. If None, converts all supported layers.
            
        Returns
        -------
        dict
            Metadata about the conversion process
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        if layers is None:
            layers = self.supported_layers
        
        conversion_meta = {
            'original_file': str(Path(h5ad_file)),  # Convert to string for JSON serialization
            'converted_layers': [],
            'output_directory': str(output_path),   # Convert to string for JSON serialization
            'conversion_timestamp': pd.Timestamp.now().isoformat()
        }
        
        try:
            with h5py.File(h5ad_file, 'r') as f:
                for layer in layers:
                    if layer in self._dump_registry:
                        try:
                            self._dump_registry[layer](f, output_path)
                            conversion_meta['converted_layers'].append(layer)
                        except Exception as e:
                            warnings.warn(f"Failed to convert layer {layer}: {str(e)}")
                    else:
                        warnings.warn(f"No handler registered for layer: {layer}")
        
        except Exception as e:
            raise ValueError(f"Failed to read h5ad file {h5ad_file}: {str(e)}")
        
        # Save conversion metadata - ensure all values are JSON serializable
        serializable_meta = self._make_json_serializable(conversion_meta)
        with open(output_path / 'conversion_metadata.json', 'w') as f:
            json.dump(serializable_meta, f, indent=2)
        
        return conversion_meta
    
    def _make_json_serializable(self, obj: Any) -> Any:
        """Recursively convert objects to JSON-serializable types."""
        if isinstance(obj, (str, int, float, bool, type(None))):
            return obj
        elif isinstance(obj, (list, tuple)):
            return [self._make_json_serializable(item) for item in obj]
        elif isinstance(obj, dict):
            return {str(key): self._make_json_serializable(value) for key, value in obj.items()}
        elif isinstance(obj, (Path, pd.Timestamp)):
            return str(obj)
        elif hasattr(obj, 'isoformat'):  # Handle datetime objects
            return obj.isoformat()
        else:
            return str(obj)
    
    def load(self, input_dir: str, layers: Optional[List[str]] = None) -> Any:
        """
        Load converted files back into an AnnData object.
        
        Parameters
        ----------
        input_dir : str
            Directory containing converted files
        layers : list, optional
            List of layers to load. If None, loads all available layers.
            
        Returns
        -------
        anndata.AnnData
            Reconstructed AnnData object
        """
        try:
            import anndata
        except ImportError:
            raise ImportError("anndata package is required for loading converted data")
        
        input_path = Path(input_dir)
        
        # Load conversion metadata
        meta_file = input_path / 'conversion_metadata.json'
        if not meta_file.exists():
            warnings.warn("Conversion metadata not found. Proceeding with available files.")
            available_layers = self._detect_available_layers(input_path)
        else:
            with open(meta_file, 'r') as f:
                conversion_meta = json.load(f)
            available_layers = conversion_meta.get('converted_layers', [])
        
        if layers is None:
            layers = available_layers
        
        adata_dict = {}
        
        for layer in layers:
            if layer in self._load_registry:
                try:
                    layer_data = self._load_registry[layer](input_path)
                    adata_dict[layer] = layer_data
                except Exception as e:
                    warnings.warn(f"Failed to load layer {layer}: {str(e)}")
            else:
                warnings.warn(f"No handler registered for loading layer: {layer}")
        
        # Construct AnnData object
        if 'X' not in adata_dict:
            raise ValueError("X matrix is required to create AnnData object")
        
        adata = anndata.AnnData(X=adata_dict['X'])
        
        # Add other layers
        for layer, data in adata_dict.items():
            if layer == 'X':
                continue
            elif layer == 'obs' and data is not None:
                adata.obs = data
            elif layer == 'var' and data is not None:
                adata.var = data
            elif layer == 'obsm' and data is not None:
                adata.obsm = data
            elif layer == 'varm' and data is not None:
                adata.varm = data
            elif layer == 'layers' and data is not None:
                adata.layers = data
            elif layer == 'uns' and data is not None:
                adata.uns = data
        
        return adata
    
    def _detect_available_layers(self, input_path: Path) -> List[str]:
        """Detect which layers are available in the input directory."""
        available = []
        for layer in self.supported_layers:
            # Check for layer-specific indicators
            if layer == 'X' and (input_path / 'X_matrix.npz').exists():
                available.append(layer)
            elif layer == 'X' and (input_path / 'X_matrix.npy').exists():
                available.append(layer)
            elif layer == 'obs' and (input_path / 'obs_metadata.csv').exists():
                available.append(layer)
            elif layer == 'var' and (input_path / 'var_metadata.csv').exists():
                available.append(layer)
            elif layer == 'layers' and (input_path / 'layers').exists():
                available.append(layer)
            elif layer == 'obsm' and (input_path / 'obsm').exists():
                available.append(layer)
            elif layer == 'varm' and (input_path / 'varm').exists():
                available.append(layer)
            elif layer == 'uns' and (input_path / 'uns_data.json').exists():
                available.append(layer)
        return available
    
    # ==================== DUMP HANDLERS ====================
    
    def _dump_X(self, h5_file: h5py.File, output_path: Path):
        """Dump X matrix data."""
        if 'X' not in h5_file:
            warnings.warn("X matrix not found in h5ad file")
            return
        
        x_data = h5_file['X']
        
        # Handle sparse matrices
        if isinstance(x_data, h5py.Group) and 'data' in x_data and 'indices' in x_data and 'indptr' in x_data:
            # CSR sparse matrix
            data = x_data['data'][:]
            indices = x_data['indices'][:]
            indptr = x_data['indptr'][:]
            shape = x_data['shape'][:]
            
            sparse_matrix = sp.csr_matrix((data, indices, indptr), shape=shape)
            sp.save_npz(output_path / 'X_matrix.npz', sparse_matrix)
            
            # Save shape info
            shape_info = {'shape': shape.tolist(), 'format': 'csr'}
            with open(output_path / 'X_matrix_info.json', 'w') as f:
                json.dump(shape_info, f)
        
        else:
            # Dense matrix
            matrix_data = x_data[:]
            np.save(output_path / 'X_matrix.npy', matrix_data)
    
    def _dump_layers(self, h5_file: h5py.File, output_path: Path):
        """Dump layers data."""
        if 'layers' not in h5_file:
            return
        
        layers_dir = output_path / 'layers'
        layers_dir.mkdir(exist_ok=True)
        
        layers_group = h5_file['layers']
        layers_info = {}
        
        for layer_name in layers_group.keys():
            layer_data = layers_group[layer_name]
            
            if isinstance(layer_data, h5py.Group) and 'data' in layer_data:
                # Sparse matrix
                data = layer_data['data'][:]
                indices = layer_data['indices'][:]
                indptr = layer_data['indptr'][:]
                shape = layer_data['shape'][:]
                
                sparse_matrix = sp.csr_matrix((data, indices, indptr), shape=shape)
                sp.save_npz(layers_dir / f'{layer_name}.npz', sparse_matrix)
                layers_info[layer_name] = {'format': 'csr', 'shape': shape.tolist()}
            
            else:
                # Dense matrix
                matrix_data = layer_data[:]
                np.save(layers_dir / f'{layer_name}.npy', matrix_data)
                layers_info[layer_name] = {'format': 'dense', 'shape': matrix_data.shape}
        
        # Save layers metadata
        with open(layers_dir / 'layers_info.json', 'w') as f:
            json.dump(layers_info, f)
    
    def _dump_obs(self, h5_file: h5py.File, output_path: Path):
        """Dump observation metadata."""
        if 'obs' not in h5_file:
            return
        
        obs_group = h5_file['obs']
        
        # Convert to pandas DataFrame
        obs_data = {}
        for col_name in obs_group.keys():
            try:
                col_data = obs_group[col_name]
                
                # Handle both Dataset and Group (categorical data)
                if isinstance(col_data, h5py.Dataset):
                    # Check data type
                    dtype = col_data.dtype
                    
                    # Handle string data
                    if h5py.check_string_dtype(dtype):
                        # String data needs special handling
                        if col_data.shape is None:
                            obs_data[str(col_name)] = [col_data[()]]
                        else:
                            # Use asstr() method to properly handle strings
                            string_data = col_data.asstr()[()]
                            obs_data[str(col_name)] = string_data
                    else:
                        # Handle numeric and other data types
                        if hasattr(col_data, 'shape') and col_data.shape is not None:
                            if len(col_data.shape) > 0 and col_data.shape[0] > 0:
                                obs_data[str(col_name)] = col_data[()]
                            else:
                                # Scalar data
                                obs_data[str(col_name)] = [col_data[()]]
                        else:
                            obs_data[str(col_name)] = [col_data[()]]
                
                elif isinstance(col_data, h5py.Group):
                    # Handle categorical data (common for string columns in AnnData)
                    # Categorical data is stored as a group with 'categories' and 'codes'
                    if 'categories' in col_data and 'codes' in col_data:
                        try:
                            # Get categories (the actual string values)
                            categories_data = col_data['categories']
                            if h5py.check_string_dtype(categories_data.dtype):
                                categories = categories_data.asstr()[()]
                            else:
                                categories = categories_data[()]
                            
                            # Get codes (indices pointing to categories)
                            codes = col_data['codes'][()]
                            
                            # Convert codes to actual string values
                            string_values = []
                            for code in codes:
                                if code >= 0 and code < len(categories):
                                    string_values.append(categories[code])
                                else:
                                    # Handle missing values (code == -1)
                                    string_values.append(None)
                            
                            obs_data[str(col_name)] = string_values
                            
                        except Exception as e:
                            warnings.warn(f"Failed to process categorical column {col_name}: {str(e)}")
                            continue
                    else:
                        # For other types of groups, try to extract what we can
                        warnings.warn(f"Column {col_name} is a Group but not a categorical. Attempting to extract data...")
                        try:
                            # Try to find any datasets in the group
                            for key in col_data.keys():
                                if isinstance(col_data[key], h5py.Dataset):
                                    dataset_data = col_data[key][()]
                                    obs_data[f"{col_name}_{key}"] = dataset_data
                        except Exception as e:
                            warnings.warn(f"Failed to extract data from group column {col_name}: {str(e)}")
                            
            except Exception as e:
                warnings.warn(f"Failed to process obs column {col_name}: {str(e)}")
                continue
        
        if obs_data:
            # Ensure all columns have consistent length
            max_len = max(len(v) for v in obs_data.values()) if obs_data else 0
            for col_name in obs_data:
                if len(obs_data[col_name]) < max_len:
                    obs_data[col_name] = list(obs_data[col_name]) + [None] * (max_len - len(obs_data[col_name]))
            
            obs_df = pd.DataFrame(obs_data)
            
            # Handle index column
            if '_index' in obs_df.columns:
                obs_df = obs_df.set_index('_index')
                obs_df.index.name = None
            elif 'index' in obs_df.columns:
                obs_df = obs_df.set_index('index')
                obs_df.index.name = None
            
            # More robust string conversion
            for col in obs_df.columns:
                if obs_df[col].dtype == object:
                    # Handle bytes type
                    if len(obs_df[col]) > 0 and any(isinstance(x, bytes) for x in obs_df[col] if x is not None):
                        obs_df[col] = obs_df[col].apply(lambda x: x.decode('utf-8') if isinstance(x, bytes) else x)
                    # Handle numpy string types
                    elif len(obs_df[col]) > 0 and any(isinstance(x, (np.str_, np.bytes_)) for x in obs_df[col] if x is not None):
                        obs_df[col] = obs_df[col].astype(str)
            
            obs_df.to_csv(output_path / 'obs_metadata.csv')
        
        # Debug output
        # print(f"Processing obs columns: {list(obs_group.keys())}")
        # for col_name in obs_group.keys():
        #     col_data = obs_group[col_name]
        #     print(f"Column: {col_name}, Type: {type(col_data)}, Shape: {getattr(col_data, 'shape', 'N/A')}, Dtype: {getattr(col_data, 'dtype', 'N/A')}")
        #     if isinstance(col_data, h5py.Group):
        #         print(f"  Group keys: {list(col_data.keys())}")
        
    def _dump_var(self, h5_file: h5py.File, output_path: Path):
        """Dump variable metadata."""
        if 'var' not in h5_file:
            return
        
        var_group = h5_file['var']
        
        # Convert to pandas DataFrame
        var_data = {}
        for col_name in var_group.keys():
            try:
                col_data = var_group[col_name]
                
                # Handle both Dataset and Group (categorical data)
                if isinstance(col_data, h5py.Dataset):
                    # Check data type
                    dtype = col_data.dtype
                    
                    # Handle string data
                    if h5py.check_string_dtype(dtype):
                        # String data needs special handling
                        if col_data.shape is None:
                            var_data[str(col_name)] = [col_data[()]]
                        else:
                            # Use asstr() method to properly handle strings
                            string_data = col_data.asstr()[()]
                            var_data[str(col_name)] = string_data
                    else:
                        # Handle numeric and other data types
                        if hasattr(col_data, 'shape') and col_data.shape is not None:
                            if len(col_data.shape) > 0 and col_data.shape[0] > 0:
                                var_data[str(col_name)] = col_data[()]
                            else:
                                # Scalar data
                                var_data[str(col_name)] = [col_data[()]]
                        else:
                            var_data[str(col_name)] = [col_data[()]]
                
                elif isinstance(col_data, h5py.Group):
                    # Handle categorical data (common for string columns in AnnData)
                    # Categorical data is stored as a group with 'categories' and 'codes'
                    if 'categories' in col_data and 'codes' in col_data:
                        try:
                            # Get categories (the actual string values)
                            categories_data = col_data['categories']
                            if h5py.check_string_dtype(categories_data.dtype):
                                categories = categories_data.asstr()[()]
                            else:
                                categories = categories_data[()]
                            
                            # Get codes (indices pointing to categories)
                            codes = col_data['codes'][()]
                            
                            # Convert codes to actual string values
                            string_values = []
                            for code in codes:
                                if code >= 0 and code < len(categories):
                                    string_values.append(categories[code])
                                else:
                                    # Handle missing values (code == -1)
                                    string_values.append(None)
                            
                            var_data[str(col_name)] = string_values
                            
                        except Exception as e:
                            warnings.warn(f"Failed to process categorical column {col_name}: {str(e)}")
                            continue
                    else:
                        # For other types of groups, try to extract what we can
                        warnings.warn(f"Column {col_name} is a Group but not a categorical. Attempting to extract data...")
                        try:
                            # Try to find any datasets in the group
                            for key in col_data.keys():
                                if isinstance(col_data[key], h5py.Dataset):
                                    dataset_data = col_data[key][()]
                                    var_data[f"{col_name}_{key}"] = dataset_data
                        except Exception as e:
                            warnings.warn(f"Failed to extract data from group column {col_name}: {str(e)}")
                            
            except Exception as e:
                warnings.warn(f"Failed to process var column {col_name}: {str(e)}")
                continue
        
        if var_data:
            # Ensure all columns have consistent length
            max_len = max(len(v) for v in var_data.values()) if var_data else 0
            for col_name in var_data:
                if len(var_data[col_name]) < max_len:
                    var_data[col_name] = list(var_data[col_name]) + [None] * (max_len - len(var_data[col_name]))
            
            var_df = pd.DataFrame(var_data)
            
            # Handle index column
            if '_index' in var_df.columns:
                var_df = var_df.set_index('_index')
                var_df.index.name = None
            elif 'index' in var_df.columns:
                var_df = var_df.set_index('index')
                var_df.index.name = None
            
            # More robust string conversion
            for col in var_df.columns:
                if var_df[col].dtype == object:
                    # Handle bytes type
                    if len(var_df[col]) > 0 and any(isinstance(x, bytes) for x in var_df[col] if x is not None):
                        var_df[col] = var_df[col].apply(lambda x: x.decode('utf-8') if isinstance(x, bytes) else x)
                    # Handle numpy string types
                    elif len(var_df[col]) > 0 and any(isinstance(x, (np.str_, np.bytes_)) for x in var_df[col] if x is not None):
                        var_df[col] = var_df[col].astype(str)
            
            var_df.to_csv(output_path / 'var_metadata.csv')
    
    def _dump_obsm(self, h5_file: h5py.File, output_path: Path):
        """Dump observation embeddings."""
        if 'obsm' not in h5_file:
            return
        
        obsm_dir = output_path / 'obsm'
        obsm_dir.mkdir(exist_ok=True)
        
        obsm_group = h5_file['obsm']
        obsm_info = {}
        
        for embedding_name in obsm_group.keys():
            embedding_data = obsm_group[embedding_name][:]
            np.save(obsm_dir / f'{embedding_name}.npy', embedding_data)
            obsm_info[embedding_name] = {'shape': embedding_data.shape}
        
        # Save obsm metadata
        with open(obsm_dir / 'obsm_info.json', 'w') as f:
            json.dump(obsm_info, f)
    
    def _dump_varm(self, h5_file: h5py.File, output_path: Path):
        """Dump variable embeddings."""
        if 'varm' not in h5_file:
            return
        
        varm_dir = output_path / 'varm'
        varm_dir.mkdir(exist_ok=True)
        
        varm_group = h5_file['varm']
        varm_info = {}
        
        for embedding_name in varm_group.keys():
            embedding_data = varm_group[embedding_name][:]
            np.save(varm_dir / f'{embedding_name}.npy', embedding_data)
            varm_info[embedding_name] = {'shape': embedding_data.shape}
        
        # Save varm metadata
        with open(varm_dir / 'varm_info.json', 'w') as f:
            json.dump(varm_info, f)
    
    def _dump_uns(self, h5_file: h5py.File, output_path: Path):
        """Dump unstructured data."""
        if 'uns' not in h5_file:
            return
        
        uns_data = self._extract_uns_data(h5_file['uns'])
        
        # Save as JSON
        with open(output_path / 'uns_data.json', 'w') as f:
            json.dump(uns_data, f, default=self._json_serializer)
    
    def _extract_uns_data(self, uns_group: h5py.Group) -> Dict[str, Any]:
        """Recursively extract unstructured data from HDF5 group."""
        result = {}
        
        for key in uns_group.keys():
            item = uns_group[key]
            
            if isinstance(item, h5py.Group):
                # Recursively process subgroups
                result[str(key)] = self._extract_uns_data(item)
            elif isinstance(item, h5py.Dataset):
                # Convert dataset to appropriate Python type
                if item.shape is None:
                    result[str(key)] = None
                elif item.shape == ():
                    # Scalar value
                    result[str(key)] = item[()]
                else:
                    # Array
                    result[str(key)] = item[:]
                
                # Convert bytes to string if necessary
                if isinstance(result[str(key)], bytes):
                    result[str(key)] = result[str(key)].decode('utf-8')
                elif isinstance(result[str(key)], np.ndarray) and result[str(key)].dtype.kind in ['S', 'U']:
                    result[str(key)] = result[str(key)].astype(str).tolist()
                elif isinstance(result[str(key)], np.ndarray):
                    result[str(key)] = result[str(key)].tolist()
        
        return result
    
    def _json_serializer(self, obj):
        """Custom JSON serializer for non-serializable objects."""
        if isinstance(obj, (np.integer, np.int64, np.int32)):
            return int(obj)
        elif isinstance(obj, (np.floating, np.float64, np.float32)):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.bool_):
            return bool(obj)
        elif isinstance(obj, (Path, pd.Timestamp)):
            return str(obj)
        else:
            return str(obj)
    
    # ==================== LOAD HANDLERS ====================
    
    def _load_X(self, input_path: Path):
        """Load X matrix data."""
        # Try sparse first
        sparse_file = input_path / 'X_matrix.npz'
        if sparse_file.exists():
            return sp.load_npz(sparse_file)
        
        # Try dense
        dense_file = input_path / 'X_matrix.npy'
        if dense_file.exists():
            return np.load(dense_file)
        
        warnings.warn("X matrix files not found")
        return None
    
    def _load_layers(self, input_path: Path):
        """Load layers data."""
        layers_dir = input_path / 'layers'
        if not layers_dir.exists():
            return {}
        
        layers_info_file = layers_dir / 'layers_info.json'
        if layers_info_file.exists():
            with open(layers_info_file, 'r') as f:
                layers_info = json.load(f)
        else:
            layers_info = {}
        
        layers_data = {}
        
        for file_path in layers_dir.glob('*.npz'):
            layer_name = file_path.stem
            if layer_name != 'layers_info':
                layers_data[layer_name] = sp.load_npz(file_path)
        
        for file_path in layers_dir.glob('*.npy'):
            layer_name = file_path.stem
            if layer_name != 'layers_info':
                layers_data[layer_name] = np.load(file_path)
        
        return layers_data
    
    def _load_obs(self, input_path: Path):
        """Load observation metadata."""
        obs_file = input_path / 'obs_metadata.csv'
        if obs_file.exists():
            df = pd.read_csv(obs_file, index_col=0)
            # Ensure index is properly set
            if df.index.name == '_index':
                df.index.name = None
            return df
        return None
    
    def _load_var(self, input_path: Path):
        """Load variable metadata."""
        var_file = input_path / 'var_metadata.csv'
        if var_file.exists():
            df = pd.read_csv(var_file, index_col=0)
            # Ensure index is properly set
            if df.index.name == '_index':
                df.index.name = None
            return df
        return None
    
    def _load_obsm(self, input_path: Path):
        """Load observation embeddings."""
        obsm_dir = input_path / 'obsm'
        if not obsm_dir.exists():
            return {}
        
        obsm_data = {}
        
        for file_path in obsm_dir.glob('*.npy'):
            embedding_name = file_path.stem
            if embedding_name != 'obsm_info':
                obsm_data[embedding_name] = np.load(file_path)
        
        return obsm_data
    
    def _load_varm(self, input_path: Path):
        """Load variable embeddings."""
        varm_dir = input_path / 'varm'
        if not varm_dir.exists():
            return {}
        
        varm_data = {}
        
        for file_path in varm_dir.glob('*.npy'):
            embedding_name = file_path.stem
            if embedding_name != 'varm_info':
                varm_data[embedding_name] = np.load(file_path)
        
        return varm_data
    
    def _load_uns(self, input_path: Path):
        """Load unstructured data."""
        uns_file = input_path / 'uns_data.json'
        if uns_file.exists():
            with open(uns_file, 'r') as f:
                return json.load(f)
        return {}


# Convenience functions
def h5ad_to_files(h5ad_file: str, output_dir: str, layers: Optional[List[str]] = None) -> Dict[str, Any]:
    """
    Convert h5ad file to generic file formats.
    
    Parameters
    ----------
    h5ad_file : str
        Path to input h5ad file
    output_dir : str
        Directory to save converted files
    layers : list, optional
        List of layers to convert. If None, converts all supported layers.
        
    Returns
    -------
    dict
        Metadata about the conversion process
    TODO:
        Fix: pd.Index type in obs with col name xxx now will be dumped as xxx_mask, xxx_values in csv
    """
    print("Note: Check data type in obs before dumping. pd.Index type in obs with col name xxx now will be dumped as xxx_mask, xxx_values in csv")
    converter = H5ADConverter()
    return converter.dump(h5ad_file, output_dir, layers)


def files_to_anndata(input_dir: str, layers: Optional[List[str]] = None) -> Any:
    """
    Load converted files back into an AnnData object.
    
    Parameters
    ----------
    input_dir : str
        Directory containing converted files
    layers : list, optional
        List of layers to load. If None, loads all available layers.
        
    Returns
    -------
    anndata.AnnData
        Reconstructed AnnData object
    """
    converter = H5ADConverter()
    return converter.load(input_dir, layers)