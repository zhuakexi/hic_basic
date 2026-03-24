"""
Tests for TI.slingshot.plot module - Plotly visualization functions
"""
import numpy as np
import pytest
import os
import tempfile


class TestPlotSlingshot:
    """Tests for plot_slingshot function"""
    
    def test_plot_slingshot_basic(self, test_adata_simple):
        """Test basic plot_slingshot functionality"""
        from hic_basic.TI.slingshot import Slingshot
        from hic_basic.TI.slingshot.plot import plot_slingshot
        
        slingshot = Slingshot(
            data=test_adata_simple,
            celltype_key='celltype',
            start_node="A"
        )
        slingshot.fit(num_epochs=5)
        
        # Should return a plotly figure
        fig = plot_slingshot(slingshot, show=False)
        
        assert fig is not None
        assert hasattr(fig, 'data')
        assert hasattr(fig, 'layout')
        # Should have scatter + curves
        assert len(fig.data) >= 2
    
    def test_plot_slingshot_with_save(self, test_adata_simple):
        """Test plot_slingshot with save path"""
        from hic_basic.TI.slingshot import Slingshot
        from hic_basic.TI.slingshot.plot import plot_slingshot
        
        slingshot = Slingshot(
            data=test_adata_simple,
            celltype_key='celltype',
            start_node="A"
        )
        slingshot.fit(num_epochs=5)
        
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
            temp_path = f.name
        
        try:
            fig = plot_slingshot(slingshot, save_path=temp_path, show=False)
            
            # Check file was created
            assert os.path.exists(temp_path)
            assert os.path.getsize(temp_path) > 0
        finally:
            if os.path.exists(temp_path):
                os.remove(temp_path)
    
    def test_plot_slingshot_unfitted_raises(self, test_adata_simple):
        """Test that plot_slingshot raises error for unfitted Slingshot"""
        from hic_basic.TI.slingshot import Slingshot
        from hic_basic.TI.slingshot.plot import plot_slingshot
        
        slingshot = Slingshot(
            data=test_adata_simple,
            celltype_key='celltype',
            start_node="A"
        )
        # Not fitted
        
        with pytest.raises(ValueError, match="not been fitted"):
            plot_slingshot(slingshot, show=False)
    
    def test_plot_slingshot_custom_params(self, test_adata_simple):
        """Test plot_slingshot with custom parameters"""
        from hic_basic.TI.slingshot import Slingshot
        from hic_basic.TI.slingshot.plot import plot_slingshot

        slingshot = Slingshot(
            data=test_adata_simple,
            celltype_key='celltype',
            start_node="A"
        )
        slingshot.fit(num_epochs=5)

        fig = plot_slingshot(
            slingshot,
            title="Custom Title",
            width=1000,
            height=800,
            point_size=10,
            curve_width=6,
            colorscale="Plasma",
            show=False
        )

        # Title is set in template, check layout dimensions
        assert fig.layout.width == 1000
        assert fig.layout.height == 800


class TestPlotClusters:
    """Tests for plot_clusters function"""
    
    def test_plot_clusters_basic(self, test_adata_simple):
        """Test basic plot_clusters functionality"""
        from hic_basic.TI.slingshot import Slingshot
        from hic_basic.TI.slingshot.plot import plot_clusters
        
        slingshot = Slingshot(
            data=test_adata_simple,
            celltype_key='celltype',
            start_node="A"
        )
        
        fig = plot_clusters(slingshot, show=False)
        
        assert fig is not None
        assert len(fig.data) >= 3  # At least 3 clusters + start node
    
    def test_plot_clusters_with_save(self, test_adata_simple):
        """Test plot_clusters with save path"""
        from hic_basic.TI.slingshot import Slingshot
        from hic_basic.TI.slingshot.plot import plot_clusters
        
        slingshot = Slingshot(
            data=test_adata_simple,
            celltype_key='celltype',
            start_node="A"
        )
        
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
            temp_path = f.name
        
        try:
            fig = plot_clusters(slingshot, save_path=temp_path, show=False)
            
            assert os.path.exists(temp_path)
            assert os.path.getsize(temp_path) > 0
        finally:
            if os.path.exists(temp_path):
                os.remove(temp_path)


class TestPlotCurves:
    """Tests for plot_curves function"""
    
    def test_plot_curves_basic(self, test_adata_simple):
        """Test basic plot_curves functionality"""
        from hic_basic.TI.slingshot import Slingshot
        from hic_basic.TI.slingshot.plot import plot_curves
        
        slingshot = Slingshot(
            data=test_adata_simple,
            celltype_key='celltype',
            start_node="A"
        )
        slingshot.fit(num_epochs=5)
        
        fig = plot_curves(slingshot, show=False)
        
        assert fig is not None
        # Should have curves + start node
        assert len(fig.data) >= 2
    
    def test_plot_curves_without_cells(self, test_adata_simple):
        """Test plot_curves without showing cells"""
        from hic_basic.TI.slingshot import Slingshot
        from hic_basic.TI.slingshot.plot import plot_curves
        
        slingshot = Slingshot(
            data=test_adata_simple,
            celltype_key='celltype',
            start_node="A"
        )
        slingshot.fit(num_epochs=5)
        
        fig = plot_curves(slingshot, show_cells=False, show=False)
        
        assert fig is not None
    
    def test_plot_curves_unfitted_raises(self, test_adata_simple):
        """Test that plot_curves raises error for unfitted Slingshot"""
        from hic_basic.TI.slingshot import Slingshot
        from hic_basic.TI.slingshot.plot import plot_curves
        
        slingshot = Slingshot(
            data=test_adata_simple,
            celltype_key='celltype',
            start_node="A"
        )
        
        with pytest.raises(ValueError, match="not been fitted"):
            plot_curves(slingshot, show=False)


class TestPlotWithSpermatogenesis:
    """Tests with spermatogenesis-like data"""
    
    def test_plot_slingshot_spermatogenesis(self, test_adata_spermatogenesis):
        """Test plot_slingshot with spermatogenesis data"""
        from hic_basic.TI.slingshot import Slingshot
        from hic_basic.TI.slingshot.plot import plot_slingshot
        
        predefined_tree = test_adata_spermatogenesis.uns['predefined_tree']
        
        slingshot = Slingshot(
            data=test_adata_spermatogenesis,
            celltype_key='celltype',
            predefined_tree=predefined_tree,
            start_node="SSC"
        )
        slingshot.fit(num_epochs=5)
        
        fig = plot_slingshot(slingshot, show=False)
        
        assert fig is not None
        # Should have multiple lineages
        assert len(fig.data) >= 2
    
    def test_plot_clusters_spermatogenesis(self, test_adata_spermatogenesis):
        """Test plot_clusters with spermatogenesis data"""
        from hic_basic.TI.slingshot import Slingshot
        from hic_basic.TI.slingshot.plot import plot_clusters
        
        slingshot = Slingshot(
            data=test_adata_spermatogenesis,
            celltype_key='celltype'
        )
        
        fig = plot_clusters(slingshot, show=False)
        
        assert fig is not None
        # Should have 14 cell types + start node
        assert len(fig.data) >= 14


class TestPlotWithBifurcation:
    """Tests with bifurcation data"""
    
    def test_plot_slingshot_bifurcation(self, test_adata_bifurcation):
        """Test plot_slingshot with bifurcation data"""
        from hic_basic.TI.slingshot import Slingshot
        from hic_basic.TI.slingshot.plot import plot_slingshot
        
        predefined_tree = test_adata_bifurcation.uns['predefined_tree']
        
        slingshot = Slingshot(
            data=test_adata_bifurcation,
            celltype_key='celltype',
            predefined_tree=predefined_tree,
            start_node="SSC"
        )
        slingshot.fit(num_epochs=5)
        
        fig = plot_slingshot(slingshot, show=False)
        
        assert fig is not None
        # Should have branching lineages
        assert len(fig.data) >= 2
