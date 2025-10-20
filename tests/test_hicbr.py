import unittest
import tempfile
import os
import shutil
from pathlib import Path
import numpy as np
import scipy.sparse as sp
import h5py
import json
import pandas as pd

# Import the module to test
from hic_basic.hicbr import H5ADConverter, h5ad_to_files, files_to_anndata


class TestH5ADConverter(unittest.TestCase):
    """
    Unit tests for H5ADConverter class.
    
    Input:
        - Mock h5ad files with various data structures
        - Different combinations of data layers to convert
        - Edge cases for data handling
    
    Output:
        - Verification of conversion process
        - Directory structure inspection
        - Data integrity validation
    """

    def setUp(self):
        """
        Set up test fixtures before each test method.
        Creates temporary directories and mock h5ad files for testing.
        """
        # Get the directory of this test file
        test_dir = Path(__file__).parent
        self.test_output_dir = test_dir / "output" / "h5ad_converter"
        
        # Create output directory if it doesn't exist
        self.test_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Temporary directory for auto-cleanup
        self.temp_dir = tempfile.mkdtemp(prefix="h5ad_converter_test_")
        
        # Create mock h5ad files for testing
        self.mock_h5ad_files = {}
        self._create_mock_h5ad_files()
        
        # Initialize converter
        self.converter = H5ADConverter()

    def tearDown(self):
        """
        Clean up after each test method.
        """
        shutil.rmtree(self.temp_dir)

    def _create_mock_h5ad_files(self):
        """Create mock h5ad files with different data structures for testing."""
        
        # Test case 1: Basic h5ad with X, obs, var
        basic_file = Path(self.temp_dir) / "basic_test.h5ad"
        with h5py.File(basic_file, 'w') as f:
            # Create X matrix (dense)
            x_data = np.random.rand(100, 50)
            f.create_dataset('X', data=x_data)
            
            # Create obs
            obs_group = f.create_group('obs')
            obs_group.create_dataset('cell_type', data=[f'cell_{i}' for i in range(100)])
            obs_group.create_dataset('batch', data=np.random.randint(0, 3, 100))
            
            # Create var
            var_group = f.create_group('var')
            var_group.create_dataset('gene_name', data=[f'gene_{i}' for i in range(50)])
            var_group.create_dataset('highly_variable', data=np.random.choice([True, False], 50))
        
        self.mock_h5ad_files['basic'] = str(basic_file)
        
        # Test case 2: Advanced h5ad with sparse matrix and embeddings
        advanced_file = Path(self.temp_dir) / "advanced_test.h5ad"
        with h5py.File(advanced_file, 'w') as f:
            # Create X matrix (sparse CSR)
            sparse_matrix = sp.random(100, 50, density=0.1, format='csr')
            x_group = f.create_group('X')
            x_group.create_dataset('data', data=sparse_matrix.data)
            x_group.create_dataset('indices', data=sparse_matrix.indices)
            x_group.create_dataset('indptr', data=sparse_matrix.indptr)
            x_group.create_dataset('shape', data=sparse_matrix.shape)
            
            # Create obs and var
            obs_group = f.create_group('obs')
            obs_group.create_dataset('cell_type', data=[f'type_{i%5}' for i in range(100)])
            
            var_group = f.create_group('var')
            var_group.create_dataset('gene_name', data=[f'gene_{i}' for i in range(50)])
            
            # Create obsm (PCA, UMAP)
            obsm_group = f.create_group('obsm')
            obsm_group.create_dataset('X_pca', data=np.random.rand(100, 10))
            obsm_group.create_dataset('X_umap', data=np.random.rand(100, 2))
            
            # Create layers
            layers_group = f.create_group('layers')
            layers_group.create_dataset('counts', data=np.random.poisson(10, (100, 50)))
            
            # Create uns
            uns_group = f.create_group('uns')
            uns_group.create_dataset('neighbors', data=np.array([1, 2, 3]))
            uns_subgroup = uns_group.create_group('pca')
            uns_subgroup.create_dataset('variance', data=np.random.rand(10))
        
        self.mock_h5ad_files['advanced'] = str(advanced_file)
        
        # Test case 3: Minimal h5ad with only X matrix
        minimal_file = Path(self.temp_dir) / "minimal_test.h5ad"
        with h5py.File(minimal_file, 'w') as f:
            x_data = np.random.rand(20, 30)
            f.create_dataset('X', data=x_data)
        
        self.mock_h5ad_files['minimal'] = str(minimal_file)
        
        # Test case 4: String-rich h5ad with various string encodings and types
        string_rich_file = Path(self.temp_dir) / "string_rich_test.h5ad"
        with h5py.File(string_rich_file, 'w') as f:
            # Create X matrix
            x_data = np.random.rand(50, 25)
            f.create_dataset('X', data=x_data)
            
            # Create obs with various string types
            obs_group = f.create_group('obs')
            
            # Regular ASCII strings
            obs_group.create_dataset('cell_id', data=[f'cell_{i:03d}' for i in range(50)])
            
            # Unicode strings
            obs_group.create_dataset('cell_type', data=[f'类型_{i%3}' for i in range(50)])
            
            # Mixed strings with special characters
            obs_group.create_dataset('description', data=[f'sample_{i}_with_data' for i in range(50)])
            
            # Categorical strings (repeated values)
            tissue_types = ['brain', 'liver', 'heart', 'kidney', 'lung']
            obs_group.create_dataset('tissue', data=[tissue_types[i % 5] for i in range(50)])
            
            # Numeric data for comparison
            obs_group.create_dataset('n_genes', data=np.random.randint(1000, 5000, 50))
            obs_group.create_dataset('percent_mito', data=np.random.rand(50))
            
            # Create var with string columns
            var_group = f.create_group('var')
            var_group.create_dataset('gene_id', data=[f'ENSG{i:06d}' for i in range(25)])
            var_group.create_dataset('gene_name', data=[f'Gene_{chr(65 + i%26)}{i//26}' for i in range(25)])
            var_group.create_dataset('chromosome', data=[f'chr{(i % 10) + 1}' for i in range(25)])
            
            # Boolean and numeric data in var
            var_group.create_dataset('highly_variable', data=np.random.choice([True, False], 25))
            var_group.create_dataset('mean_expression', data=np.random.rand(25))

    def _print_directory_structure(self, directory_path, indent=0):
        """
        Recursively print the directory structure for inspection.
        
        Input:
            directory_path: Path to the directory to print
            indent: Current indentation level for pretty printing
        """
        directory = Path(directory_path)
        if not directory.exists():
            print(" " * indent + f"Directory does not exist: {directory}")
            return
        
        items = sorted(directory.iterdir())
        for item in items:
            if item.is_dir():
                print(" " * indent + f"📁 {item.name}/")
                self._print_directory_structure(item, indent + 2)
            else:
                size = item.stat().st_size
                print(" " * indent + f"📄 {item.name} ({size} bytes)")

    def test_string_columns_export(self):
        """
        Test export of string columns in obs and var metadata.
        
        Input:
            - h5ad file with various string column types (ASCII, Unicode, categorical)
            - Mixed data types including strings, numbers, and booleans
        
        Output:
            - All string columns properly exported to CSV
            - Correct data types preserved
            - Unicode and special characters handled correctly
        """
        # Create a test h5ad file with rich string data
        string_test_file = Path(self.temp_dir) / "string_test.h5ad"
        with h5py.File(string_test_file, 'w') as f:
            # Create X matrix
            x_data = np.random.rand(30, 20)
            f.create_dataset('X', data=x_data)
            
            # Create obs with various string types
            obs_group = f.create_group('obs')
            
            # ASCII strings
            obs_group.create_dataset('cell_id', data=[f'cell_{i:03d}' for i in range(30)])
            
            # Unicode strings
            obs_group.create_dataset('cell_type', data=[f'类型_{i%3}' for i in range(30)])
            
            # Mixed content
            obs_group.create_dataset('sample_info', data=[f'sample_{i}_data' for i in range(30)])
            
            # Numeric data for comparison
            obs_group.create_dataset('n_counts', data=np.random.randint(1000, 10000, 30))
            
            # Create var with string columns
            var_group = f.create_group('var')
            var_group.create_dataset('gene_name', data=[f'gene_{i}' for i in range(20)])
            var_group.create_dataset('chromosome', data=[f'chr{(i % 5) + 1}' for i in range(20)])

        output_dir = Path(self.temp_dir) / "string_test_output"
        
        # Perform conversion
        conversion_meta = self.converter.dump(
            str(string_test_file),
            str(output_dir),
            layers=['X', 'obs', 'var']
        )
        
        print("\n=== String Columns Export Directory Structure ===")
        self._print_directory_structure(output_dir)
        
        # Verify conversion metadata
        self.assertIn('obs', conversion_meta['converted_layers'])
        self.assertIn('var', conversion_meta['converted_layers'])
        
        # Verify files were created
        self.assertTrue((output_dir / 'obs_metadata.csv').exists())
        self.assertTrue((output_dir / 'var_metadata.csv').exists())
        
        # Load and inspect the exported CSV files
        obs_df = pd.read_csv(output_dir / 'obs_metadata.csv', index_col=0)
        var_df = pd.read_csv(output_dir / 'var_metadata.csv', index_col=0)
        
        print(f"\nExported obs columns: {list(obs_df.columns)}")
        print(f"Exported var columns: {list(var_df.columns)}")
        print(f"Obs dtypes: {obs_df.dtypes}")
        print(f"Var dtypes: {var_df.dtypes}")
        
        # Verify all expected columns are present
        expected_obs_columns = ['cell_id', 'cell_type', 'sample_info', 'n_counts']
        expected_var_columns = ['gene_name', 'chromosome']
        
        for col in expected_obs_columns:
            self.assertIn(col, obs_df.columns, f"Expected column '{col}' not found in obs metadata")
        
        for col in expected_var_columns:
            self.assertIn(col, var_df.columns, f"Expected column '{col}' not found in var metadata")
        
        # Verify string columns contain expected data
        self.assertTrue(all(obs_df['cell_id'].str.startswith('cell_')), "cell_id column doesn't contain expected string data")
        self.assertTrue(all('类型' in str(x) for x in obs_df['cell_type']), "cell_type column doesn't contain expected unicode data")
        self.assertTrue(all(obs_df['sample_info'].str.contains('sample_')), "sample_info column doesn't contain expected string data")
        
        self.assertTrue(all(var_df['gene_name'].str.startswith('gene_')), "gene_name column doesn't contain expected string data")
        self.assertTrue(all(var_df['chromosome'].str.startswith('chr')), "chromosome column doesn't contain expected string data")
        
        # Verify data integrity - check a few sample values
        self.assertEqual(obs_df.loc[0, 'cell_id'], 'cell_000')
        self.assertEqual(var_df.loc[0, 'gene_name'], 'gene_0')
        
        # Verify numeric column is properly exported
        self.assertEqual(obs_df['n_counts'].dtype, np.dtype('int64'), "n_counts column should be integer type")

    def test_basic_conversion(self):
        """
        Test basic conversion with X, obs, and var layers.
        
        Input:
            - Basic h5ad file with dense X matrix, obs, and var
            - Output directory path
        
        Output:
            - Directory structure with converted files
            - Verification of successful conversion
        """
        output_dir = Path(self.temp_dir) / "basic_conversion_output"
        
        # Perform conversion
        conversion_meta = self.converter.dump(
            self.mock_h5ad_files['basic'],
            str(output_dir),
            layers=['X', 'obs', 'var']
        )
        
        # Print directory structure
        print("\n=== Basic Conversion Directory Structure ===")
        self._print_directory_structure(output_dir)
        
        # Verify conversion metadata
        self.assertIn('X', conversion_meta['converted_layers'])
        self.assertIn('obs', conversion_meta['converted_layers'])
        self.assertIn('var', conversion_meta['converted_layers'])
        
        # Verify files were created
        self.assertTrue((output_dir / 'X_matrix.npy').exists())
        self.assertTrue((output_dir / 'obs_metadata.csv').exists())
        self.assertTrue((output_dir / 'var_metadata.csv').exists())
        self.assertTrue((output_dir / 'conversion_metadata.json').exists())

    def test_advanced_conversion(self):
        """
        Test advanced conversion with all data layers.
        
        Input:
            - Advanced h5ad file with sparse matrix, embeddings, layers, and uns
            - All supported layers for conversion
        
        Output:
            - Complex directory structure with multiple subdirectories
            - Verification of all data layers
        """
        output_dir = Path(self.temp_dir) / "advanced_conversion_output"
        
        # Perform conversion
        conversion_meta = self.converter.dump(
            self.mock_h5ad_files['advanced'],
            str(output_dir),
            layers=['X', 'obs', 'var', 'obsm', 'varm', 'layers', 'uns']
        )
        
        # Print directory structure
        print("\n=== Advanced Conversion Directory Structure ===")
        self._print_directory_structure(output_dir)
        
        # Verify all layers were converted
        expected_layers = ['X', 'obs', 'var', 'obsm', 'varm', 'layers', 'uns']
        for layer in expected_layers:
            self.assertIn(layer, conversion_meta['converted_layers'])
        
        # Verify directory structure
        self.assertTrue((output_dir / 'X_matrix.npz').exists())  # Should be sparse
        self.assertTrue((output_dir / 'X_matrix_info.json').exists())
        self.assertTrue((output_dir / 'obs_metadata.csv').exists())
        self.assertTrue((output_dir / 'var_metadata.csv').exists())
        self.assertTrue((output_dir / 'obsm').exists())
        self.assertTrue((output_dir / 'layers').exists())
        self.assertTrue((output_dir / 'uns_data.json').exists())
        
        # Verify subdirectory contents
        obsm_dir = output_dir / 'obsm'
        self.assertTrue((obsm_dir / 'X_pca.npy').exists())
        self.assertTrue((obsm_dir / 'X_umap.npy').exists())
        self.assertTrue((obsm_dir / 'obsm_info.json').exists())
        
        layers_dir = output_dir / 'layers'
        self.assertTrue((layers_dir / 'counts.npy').exists())
        self.assertTrue((layers_dir / 'layers_info.json').exists())

    def test_round_trip_conversion(self):
        """
        Test complete round-trip conversion: h5ad -> files -> anndata.
        
        Input:
            - Advanced h5ad file
            - Complete conversion and reconstruction
        
        Output:
            - Successful reconstruction of AnnData object
            - Data integrity validation
        """
        try:
            import anndata
        except ImportError:
            self.skipTest("anndata package not available")
        
        output_dir = Path(self.temp_dir) / "round_trip_output"
        
        # Convert h5ad to files
        conversion_meta = self.converter.dump(
            self.mock_h5ad_files['advanced'],
            str(output_dir)
        )
        
        print("\n=== Round-trip Conversion Directory Structure ===")
        self._print_directory_structure(output_dir)
        
        # Load files back to AnnData
        adata_reconstructed = self.converter.load(str(output_dir))
        
        # Verify basic properties
        self.assertIsInstance(adata_reconstructed, anndata.AnnData)
        self.assertEqual(adata_reconstructed.n_obs, 100)
        self.assertEqual(adata_reconstructed.n_vars, 50)
        
        # Verify data integrity
        self.assertIn('cell_type', adata_reconstructed.obs.columns)
        self.assertIn('gene_name', adata_reconstructed.var.columns)
        self.assertIn('X_pca', adata_reconstructed.obsm)
        self.assertIn('counts', adata_reconstructed.layers)

    def test_selective_layer_conversion(self):
        """
        Test conversion with only selected layers.
        
        Input:
            - Advanced h5ad file
            - Subset of layers to convert
        
        Output:
            - Directory structure with only selected layers
            - Verification of selective conversion
        """
        output_dir = Path(self.temp_dir) / "selective_output"
        
        # Convert only specific layers
        conversion_meta = self.converter.dump(
            self.mock_h5ad_files['advanced'],
            str(output_dir),
            layers=['X', 'obs', 'obsm']
        )
        
        print("\n=== Selective Conversion Directory Structure ===")
        self._print_directory_structure(output_dir)
        
        # Verify only selected layers were converted
        self.assertEqual(set(conversion_meta['converted_layers']), {'X', 'obs', 'obsm'})
        
        # Verify only expected files exist
        self.assertTrue((output_dir / 'X_matrix.npz').exists())
        self.assertTrue((output_dir / 'obs_metadata.csv').exists())
        self.assertTrue((output_dir / 'obsm').exists())
        
        # Verify other layers were NOT converted
        self.assertFalse((output_dir / 'var_metadata.csv').exists())
        self.assertFalse((output_dir / 'layers').exists())
        self.assertFalse((output_dir / 'varm').exists())

    def test_convenience_functions(self):
        """
        Test the convenience functions h5ad_to_files and files_to_anndata.
        
        Input:
            - Basic h5ad file
            - Convenience function calls
        
        Output:
            - Successful conversion using convenience functions
            - Proper function behavior
        """
        try:
            import anndata
        except ImportError:
            self.skipTest("anndata package not available")
        
        output_dir = Path(self.temp_dir) / "convenience_output"
        
        # Test h5ad_to_files
        conversion_meta = h5ad_to_files(
            self.mock_h5ad_files['basic'],
            str(output_dir)
        )
        
        print("\n=== Convenience Functions Directory Structure ===")
        self._print_directory_structure(output_dir)
        
        # Test files_to_anndata
        adata = files_to_anndata(str(output_dir))
        
        self.assertIsInstance(adata, anndata.AnnData)
        self.assertGreater(adata.n_obs, 0)
        self.assertGreater(adata.n_vars, 0)

    def test_error_handling(self):
        """
        Test error handling for invalid inputs.
        
        Input:
            - Non-existent h5ad file
            - Invalid output directory
            - Empty conversion layers
        
        Output:
            - Appropriate exceptions raised
            - Graceful error handling
        """
        # Test non-existent file
        with self.assertRaises(ValueError):
            self.converter.dump("non_existent_file.h5ad", self.temp_dir)
        
        # Test invalid output directory (read-only)
        if os.name != 'nt':  # Skip on Windows
            read_only_dir = Path(self.temp_dir) / "readonly"
            read_only_dir.mkdir()
            read_only_dir.chmod(0o444)
            
            with self.assertRaises(Exception):  # Could be various exceptions
                self.converter.dump(self.mock_h5ad_files['basic'], str(read_only_dir))
            
            # Restore permissions for cleanup
            read_only_dir.chmod(0o755)

    def test_custom_handler_registration(self):
        """
        Test custom handler registration functionality.
        
        Input:
            - Custom dump and load handlers
            - Mock h5ad file
        
        Output:
            - Successful handler registration
            - Custom handler execution
        """
        # Define custom handlers
        def custom_dump_handler(h5_file, output_path):
            """Custom dump handler that creates a test file."""
            test_file = output_path / "custom_dump_test.txt"
            test_file.write_text("Custom dump handler executed")
        
        def custom_load_handler(input_path):
            """Custom load handler that reads the test file."""
            test_file = input_path / "custom_dump_test.txt"
            if test_file.exists():
                return test_file.read_text()
            return None
        
        # Register custom handlers
        self.converter.register_dump_handler('custom_layer', custom_dump_handler)
        self.converter.register_load_handler('custom_layer', custom_load_handler)
        
        output_dir = Path(self.temp_dir) / "custom_handler_output"
        
        # Perform conversion including custom layer
        conversion_meta = self.converter.dump(
            self.mock_h5ad_files['basic'],
            str(output_dir),
            layers=['X', 'custom_layer']
        )
        
        print("\n=== Custom Handler Directory Structure ===")
        self._print_directory_structure(output_dir)
        
        # Verify custom handler was executed
        self.assertIn('custom_layer', conversion_meta['converted_layers'])
        self.assertTrue((output_dir / "custom_dump_test.txt").exists())
        
        # Test custom load handler
        result = self.converter._load_registry['custom_layer'](output_dir)
        self.assertEqual(result, "Custom dump handler executed")


if __name__ == '__main__':
    # Run the tests
    unittest.main(verbosity=2)