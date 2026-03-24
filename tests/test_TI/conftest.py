"""
Pytest fixtures for TI module tests
"""
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData


@pytest.fixture
def random_state():
    """Default random seed for reproducibility"""
    return 42


@pytest.fixture
def n_cells():
    """Default number of cells for test data"""
    return 100


@pytest.fixture
def n_genes():
    """Default number of genes for test data"""
    return 50


@pytest.fixture
def test_adata_simple(n_cells, n_genes, random_state):
    """
    Simple AnnData object for basic testing
    
    Returns:
        AnnData: Simple test data with 3 cell types
    """
    np.random.seed(random_state)
    
    # Generate simple 2D embedding with 3 clusters
    cell_types = ["A", "B", "C"]
    cells_per_type = [n_cells // 3] * 3
    cells_per_type[0] += n_cells - sum(cells_per_type)  # Adjust for remainder
    
    all_cells = []
    all_labels = []
    
    for i, (ct, n) in enumerate(zip(cell_types, cells_per_type)):
        base_x = i * 3
        base_y = 0
        cells_x = base_x + np.random.randn(n) * 0.5
        cells_y = base_y + np.random.randn(n) * 0.5
        all_cells.append(np.column_stack([cells_x, cells_y]))
        all_labels.extend([ct] * n)
    
    X_umap = np.vstack(all_cells)
    
    return AnnData(
        X=np.random.randn(n_cells, n_genes),
        obs=pd.DataFrame({'celltype': all_labels}),
        obsm={'X_umap': X_umap}
    )


@pytest.fixture
def test_adata_spermatogenesis(n_cells, n_genes, random_state):
    """
    Generate simulated spermatogenesis trajectory test data

    Simulates cell differentiation path:
    SSC → Undiff.SPG → Diff.SPG → preL → L → Z → preP → P → D → SPC → S1 → S2 → S3 → S4

    Returns:
        AnnData: Test data with:
            - obs['celltype']: Cell type labels (strings)
            - obsm['X_umap']: 2D UMAP embedding
            - uns['predefined_tree']: Predefined tree structure
    """
    np.random.seed(random_state)
    
    # Use larger n_cells for spermatogenesis to accommodate all cell types
    n_cells_actual = max(n_cells, 520)

    # Define spermatogenesis cell types (simplified)
    cell_types = [
        "SSC",           # 0: Spermatogonial stem cells
        "Undiff.SPG",    # 1: Undifferentiated spermatogonia
        "Diff.SPG",      # 2: Differentiated spermatogonia
        "preL",          # 3: Pre-leptotene
        "L",             # 4: Leptotene
        "Z",             # 5: Zygotene
        "preP",          # 6: Pre-pachytene
        "P",             # 7: Pachytene
        "D",             # 8: Diplotene
        "SPC",           # 9: Spermatids
        "S1",            # 10: Sperm 1
        "S2",            # 11: Sperm 2
        "S3",            # 12: Sperm 3
        "S4"             # 13: Sperm 4
    ]
    
    n_cell_types = len(cell_types)
    cells_per_type = [n_cells_actual // n_cell_types] * n_cell_types
    cells_per_type[0] += n_cells_actual - sum(cells_per_type)  # Adjust for remainder

    # Generate UMAP embedding (arc structure along trajectory)
    all_cells = []
    all_labels = []

    for i, (ct, n) in enumerate(zip(cell_types, cells_per_type)):
        # Parameter position along arc
        t = i / n_cell_types * np.pi  # 0 to π
        base_x = np.cos(t) * 5
        base_y = np.sin(t) * 5

        # Add cell type variation
        cells_x = base_x + np.random.randn(n) * 0.8
        cells_y = base_y + np.random.randn(n) * 0.8

        all_cells.append(np.column_stack([cells_x, cells_y]))
        all_labels.extend([ct] * n)
    
    X_umap = np.vstack(all_cells)

    # Create AnnData object
    adata = AnnData(
        X=np.random.randn(n_cells_actual, n_genes),
        obs=pd.DataFrame({'celltype': all_labels}),
        obsm={'X_umap': X_umap}
    )

    # Define predefined tree structure
    predefined_tree = {
        "SSC": ["Undiff.SPG"],
        "Undiff.SPG": ["Diff.SPG"],
        "Diff.SPG": ["preL"],
        "preL": ["L"],
        "L": ["Z"],
        "Z": ["preP"],
        "preP": ["P"],
        "P": ["D"],
        "D": ["SPC"],
        "SPC": ["S1"],
        "S1": ["S2"],
        "S2": ["S3"],
        "S3": ["S4"]
    }

    adata.uns['predefined_tree'] = predefined_tree

    return adata


@pytest.fixture
def test_adata_bifurcation(n_cells, n_genes, random_state):
    """
    Generate test data with bifurcation (e.g., L/Z branch)

    Simulates differentiation path with branching:
    preL → L → L2 → L3 → preP
           ↘ Z ↗

    Returns:
        AnnData: Test data with bifurcation
    """
    np.random.seed(random_state)
    
    # Use larger n_cells for bifurcation test
    n_cells_actual = max(n_cells, 600)

    # Define cell types (with branch)
    cell_types = [
        "SSC",
        "Undiff.SPG",
        "Diff.SPG",
        "preL",
        "L",
        "Z",             # Branch point
        "L2",
        "L3",
        "preP",
        "P",
        "D",
        "SPC"
    ]
    
    n_cell_types = len(cell_types)
    cells_per_type = [n_cells_actual // n_cell_types] * n_cell_types
    cells_per_type[0] += n_cells_actual - sum(cells_per_type)  # Adjust for remainder

    # Generate Y-shaped trajectory
    all_cells = []
    all_labels = []

    for i, (ct, n) in enumerate(zip(cell_types, cells_per_type)):
        if ct in ["L", "Z"]:
            # Branch point: two cell types at same location
            t = 0.5
            base_x = 0
            base_y = 0
        elif ct in ["L2", "L3", "preP", "P", "D", "SPC"]:
            # L branch: upper right
            idx = cell_types.index(ct) - cell_types.index("L2")
            t = idx / 6 * np.pi / 3
            base_x = np.cos(t) * 4 + 2
            base_y = np.sin(t) * 4 + 2
        else:
            # Main trunk: lower left
            idx = cell_types.index(ct)
            t = idx / 4 * np.pi / 2
            base_x = -np.cos(t) * 5
            base_y = -np.sin(t) * 5

        cells_x = base_x + np.random.randn(n) * 0.6
        cells_y = base_y + np.random.randn(n) * 0.6

        all_cells.append(np.column_stack([cells_x, cells_y]))
        all_labels.extend([ct] * n)
    
    X_umap = np.vstack(all_cells)
    
    adata = AnnData(
        X=np.random.randn(n_cells_actual, n_genes),
        obs=pd.DataFrame({'celltype': all_labels}),
        obsm={'X_umap': X_umap}
    )
    
    # Define tree structure with branch
    predefined_tree = {
        "SSC": ["Undiff.SPG"],
        "Undiff.SPG": ["Diff.SPG"],
        "Diff.SPG": ["preL"],
        "preL": ["L"],
        "L": ["Z", "L2"],  # Branch point
        "Z": ["preP"],     # Z branch
        "L2": ["L3"],      # L branch
        "L3": ["preP"],    # Convergence point
        "preP": ["P"],
        "P": ["D"],
        "D": ["SPC"]
    }
    
    adata.uns['predefined_tree'] = predefined_tree
    
    return adata


@pytest.fixture
def test_adata_integer_labels(n_cells=100, n_genes=50, random_state=42):
    """
    AnnData object with integer cluster labels (not strings)

    Returns:
        AnnData: Test data with integer labels
    """
    np.random.seed(random_state)

    # Generate data with 3 clusters
    n_clusters = 3
    cells_per_cluster = [n_cells // n_clusters] * n_clusters
    cells_per_cluster[0] += n_cells - sum(cells_per_cluster)  # Adjust for remainder

    all_cells = []
    all_labels = []

    for i, n in enumerate(cells_per_cluster):
        base_x = i * 3
        base_y = 0
        cells_x = base_x + np.random.randn(n) * 0.5
        cells_y = base_y + np.random.randn(n) * 0.5
        all_cells.append(np.column_stack([cells_x, cells_y]))
        all_labels.extend([i] * n)

    X_umap = np.vstack(all_cells)

    return AnnData(
        X=np.random.randn(n_cells, n_genes),
        obs=pd.DataFrame({'cluster': all_labels}),
        obsm={'X_umap': X_umap}
    )
