"""
Tests for TI.slingshot.lineage module - Lineage class
"""
import pytest
from hic_basic.TI.slingshot.lineage import Lineage


class TestLineageInit:
    """Tests for Lineage initialization"""
    
    def test_init_with_clusters(self):
        """Test initialization with cluster list"""
        clusters = [0, 1, 2, 3]
        lineage = Lineage(clusters)
        
        assert lineage.clusters == clusters
        assert len(lineage) == 4
    
    def test_init_with_string_clusters(self):
        """Test initialization with string cluster labels"""
        clusters = ["SSC", "Undiff.SPG", "Diff.SPG", "preL"]
        lineage = Lineage(clusters)
        
        assert lineage.clusters == clusters
        assert len(lineage) == 4
    
    def test_init_with_single_cluster(self):
        """Test initialization with single cluster"""
        clusters = [0]
        lineage = Lineage(clusters)
        
        assert lineage.clusters == clusters
        assert len(lineage) == 1
    
    def test_init_with_empty_clusters(self):
        """Test initialization with empty cluster list"""
        clusters = []
        lineage = Lineage(clusters)
        
        assert lineage.clusters == clusters
        assert len(lineage) == 0


class TestLineageProperties:
    """Tests for Lineage properties"""
    
    def test_start_property(self):
        """Test start property returns first cluster"""
        clusters = [0, 1, 2, 3]
        lineage = Lineage(clusters)
        
        assert lineage.start == 0
    
    def test_start_property_strings(self):
        """Test start property with string labels"""
        clusters = ["SSC", "Undiff.SPG", "Diff.SPG"]
        lineage = Lineage(clusters)
        
        assert lineage.start == "SSC"
    
    def test_end_property(self):
        """Test end property returns last cluster"""
        clusters = [0, 1, 2, 3]
        lineage = Lineage(clusters)
        
        assert lineage.end == 3
    
    def test_end_property_strings(self):
        """Test end property with string labels"""
        clusters = ["SSC", "Undiff.SPG", "Diff.SPG"]
        lineage = Lineage(clusters)
        
        assert lineage.end == "Diff.SPG"
    
    def test_start_end_same_for_single(self):
        """Test start and end are same for single cluster"""
        clusters = [42]
        lineage = Lineage(clusters)
        
        assert lineage.start == 42
        assert lineage.end == 42


class TestLineageIteration:
    """Tests for Lineage iteration"""
    
    def test_iter_integer_clusters(self):
        """Test iteration over integer clusters"""
        clusters = [0, 1, 2, 3]
        lineage = Lineage(clusters)
        
        result = list(lineage)
        assert result == clusters
    
    def test_iter_string_clusters(self):
        """Test iteration over string clusters"""
        clusters = ["A", "B", "C"]
        lineage = Lineage(clusters)
        
        result = list(lineage)
        assert result == clusters
    
    def test_iter_empty_lineage(self):
        """Test iteration over empty lineage"""
        clusters = []
        lineage = Lineage(clusters)
        
        result = list(lineage)
        assert result == []
    
    def test_iter_using_for_loop(self):
        """Test iteration using for loop"""
        clusters = [0, 1, 2]
        lineage = Lineage(clusters)
        
        result = []
        for cluster in lineage:
            result.append(cluster)
        
        assert result == clusters


class TestLineageRepresentation:
    """Tests for Lineage string representation"""
    
    def test_repr_integer_clusters(self):
        """Test string representation with integer clusters"""
        clusters = [0, 1, 2]
        lineage = Lineage(clusters)
        
        repr_str = repr(lineage)
        assert "Lineage" in repr_str
        assert "[0, 1, 2]" in repr_str or "(0, 1, 2)" in repr_str or "0, 1, 2" in repr_str
    
    def test_repr_string_clusters(self):
        """Test string representation with string clusters"""
        clusters = ["A", "B", "C"]
        lineage = Lineage(clusters)
        
        repr_str = repr(lineage)
        assert "Lineage" in repr_str
        assert "['A', 'B', 'C']" in repr_str or "('A', 'B', 'C')" in repr_str


class TestLineageLength:
    """Tests for Lineage length"""
    
    def test_len_basic(self):
        """Test basic length calculation"""
        clusters = [0, 1, 2, 3, 4]
        lineage = Lineage(clusters)
        
        assert len(lineage) == 5
    
    def test_len_empty(self):
        """Test length of empty lineage"""
        clusters = []
        lineage = Lineage(clusters)
        
        assert len(lineage) == 0
    
    def test_len_single(self):
        """Test length of single cluster lineage"""
        clusters = [42]
        lineage = Lineage(clusters)
        
        assert len(lineage) == 1
