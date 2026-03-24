"""
Tests for TI.slingshot.utils module - Utility functions
"""
import numpy as np
import pytest
from hic_basic.TI.slingshot.utils import (
    mahalanobis,
    isint,
    isstr,
    scale_to_range,
    infer_cluster_label_indices
)


class TestMahalanobis:
    """Tests for mahalanobis function"""
    
    def test_mahalanobis_basic(self):
        """Test basic Mahalanobis distance calculation"""
        X1 = np.array([1.0, 2.0])
        X2 = np.array([2.0, 3.0])
        S1 = np.eye(2)
        S2 = np.eye(2)
        
        dist = mahalanobis(X1, X2, S1, S2)
        
        # Should return a positive scalar
        assert dist.shape == (1, 1)
        assert dist[0, 0] > 0
    
    def test_mahalanobis_same_point(self):
        """Test Mahalanobis distance for same point (should be 0)"""
        X1 = np.array([1.0, 2.0])
        X2 = np.array([1.0, 2.0])
        S1 = np.eye(2)
        S2 = np.eye(2)
        
        dist = mahalanobis(X1, X2, S1, S2)
        
        # Distance to same point should be 0
        assert np.isclose(dist[0, 0], 0)
    
    def test_mahalanobis_symmetric(self):
        """Test that Mahalanobis distance is symmetric"""
        X1 = np.array([1.0, 2.0])
        X2 = np.array([3.0, 4.0])
        S1 = np.eye(2)
        S2 = np.eye(2)
        
        dist12 = mahalanobis(X1, X2, S1, S2)
        dist21 = mahalanobis(X2, X1, S1, S2)
        
        assert np.isclose(dist12[0, 0], dist21[0, 0])
    
    def test_mahalanobis_different_covariance(self):
        """Test Mahalanobis distance with different covariance matrices"""
        X1 = np.array([1.0, 2.0])
        X2 = np.array([2.0, 3.0])
        S1 = np.array([[2.0, 0.5], [0.5, 2.0]])
        S2 = np.array([[1.0, 0.0], [0.0, 1.0]])
        
        dist = mahalanobis(X1, X2, S1, S2)
        
        assert dist.shape == (1, 1)
        assert dist[0, 0] > 0
    
    def test_mahalanobis_3d(self):
        """Test Mahalanobis distance in 3D"""
        X1 = np.array([1.0, 2.0, 3.0])
        X2 = np.array([2.0, 3.0, 4.0])
        S1 = np.eye(3)
        S2 = np.eye(3)
        
        dist = mahalanobis(X1, X2, S1, S2)
        
        assert dist.shape == (1, 1)
        assert dist[0, 0] > 0


class TestIsint:
    """Tests for isint function"""
    
    def test_isint_integer(self):
        """Test isint with integer"""
        assert isint(5) == True
        assert isint(0) == True
        assert isint(-3) == True
    
    def test_isint_float(self):
        """Test isint with float"""
        assert isint(5.0) == True
        assert isint(3.14) == False
        assert isint(0.0) == True
    
    def test_isint_string_integer(self):
        """Test isint with string representing integer"""
        assert isint("5") == False  # Strings return False per implementation
        assert isint("abc") == False
    
    def test_isint_string(self):
        """Test isint with non-numeric string"""
        assert isint("hello") == False
    
    def test_isint_numpy_integer(self):
        """Test isint with numpy integer"""
        assert isint(np.int64(5)) == True
        assert isint(np.int32(3)) == True
    
    def test_isint_numpy_float(self):
        """Test isint with numpy float"""
        assert isint(np.float64(5.0)) == True
        assert isint(np.float64(3.14)) == False


class TestIsstr:
    """Tests for isstr function"""
    
    def test_isstr_string(self):
        """Test isstr with string"""
        assert isstr("hello") == True
        assert isstr("") == True
        assert isstr("123") == True
    
    def test_isstr_not_string(self):
        """Test isstr with non-string types"""
        assert isstr(123) == False
        assert isstr(3.14) == False
        assert isstr([1, 2, 3]) == False
        assert isstr(None) == False
    
    def test_isstr_numpy_string(self):
        """Test isstr with numpy string"""
        assert isstr(np.str_("hello")) == True


class TestScaleToRange:
    """Tests for scale_to_range function"""
    
    def test_scale_to_range_basic(self):
        """Test basic scaling to [0, 1]"""
        x = np.array([1, 2, 3, 4, 5])
        scaled = scale_to_range(x)
        
        assert np.isclose(scaled.min(), 0)
        assert np.isclose(scaled.max(), 1)
        assert len(scaled) == len(x)
    
    def test_scale_to_range_custom_range(self):
        """Test scaling to custom range"""
        x = np.array([1, 2, 3, 4, 5])
        scaled = scale_to_range(x, a=-1, b=1)
        
        assert np.isclose(scaled.min(), -1)
        assert np.isclose(scaled.max(), 1)
    
    def test_scale_to_range_single_value(self):
        """Test scaling with single value (edge case)"""
        x = np.array([5.0])
        scaled = scale_to_range(x)
        
        # Single value should handle gracefully (may return 0 or NaN)
        assert len(scaled) == 1
    
    def test_scale_to_range_preserves_order(self):
        """Test that scaling preserves order"""
        x = np.array([1, 3, 2, 5, 4])
        scaled = scale_to_range(x)
        
        # Order should be preserved
        assert scaled[1] > scaled[0]  # 3 > 1
        assert scaled[3] > scaled[2]  # 5 > 2
    
    def test_scale_to_range_negative_values(self):
        """Test scaling with negative values"""
        x = np.array([-5, -3, -1, 0, 1, 3, 5])
        scaled = scale_to_range(x)
        
        assert np.isclose(scaled.min(), 0)
        assert np.isclose(scaled.max(), 1)


class TestInferClusterLabelIndices:
    """Tests for infer_cluster_label_indices function"""
    
    def test_infer_cluster_label_indices_integer(self):
        """Test with integer cluster labels"""
        labels = np.array([0, 1, 2, 0, 1, 2])
        result = infer_cluster_label_indices(labels)
        
        assert np.array_equal(result, labels)
    
    def test_infer_cluster_label_indices_string(self):
        """Test with string cluster labels"""
        labels = np.array(["A", "B", "C", "A", "B", "C"])
        result = infer_cluster_label_indices(labels)
        
        # Should convert to integer indices
        assert len(result) == len(labels)
        assert result.dtype in [np.int64, np.int32, int]
        # Same labels should map to same indices
        assert result[0] == result[3]  # Both "A"
        assert result[1] == result[4]  # Both "B"
    
    def test_infer_cluster_label_indices_mixed_strings(self):
        """Test with different string labels"""
        labels = np.array(["SSC", "Undiff", "Diff", "SSC", "Undiff"])
        result = infer_cluster_label_indices(labels)
        
        assert len(result) == 5
        assert result[0] == result[3]  # Both "SSC"
        assert result[1] == result[4]  # Both "Undiff"
    
    def test_infer_cluster_label_indices_invalid_type(self):
        """Test with invalid label type"""
        labels = [1.5, 2.5, 3.5]  # Floats (not int or str)
        
        with pytest.raises(ValueError, match="Unexpected cluster label dtype"):
            infer_cluster_label_indices(labels)
    
    def test_infer_cluster_label_indices_single_label(self):
        """Test with single cluster label"""
        labels = np.array(["A", "A", "A", "A"])
        result = infer_cluster_label_indices(labels)
        
        assert len(result) == 4
        # All should map to same index
        assert len(np.unique(result)) == 1
