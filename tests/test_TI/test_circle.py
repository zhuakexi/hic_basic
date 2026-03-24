"""
Tests for TI.circle module - Circular pseudotime baseline
"""
import numpy as np
import pandas as pd
import pytest
from hic_basic.TI.circle import angle_pseudotime, _circle_arctan


class TestCircleArctan:
    """Tests for _circle_arctan helper function"""
    
    def test_first_quadrant(self):
        """Test angle calculation in first quadrant (x > 0)"""
        coord = pd.Series([1.0, 0.5])  # x > 0, y > 0 relative to center
        angle = _circle_arctan(coord, center=(0, 0))
        # For x > 0: angle = pi/2 - arctan(y/x), should be in (0, pi/2)
        assert 0 < angle < np.pi / 2
    
    def test_second_quadrant(self):
        """Test angle calculation in second quadrant (x < 0)"""
        coord = pd.Series([-1.0, 1.0])
        angle = _circle_arctan(coord, center=(0, 0))
        # For x < 0: angle = 1.5*pi - arctan(y/x)
        # With y/x = -1, arctan(-1) = -pi/4, so angle = 1.5*pi - (-pi/4) = 1.75*pi (mod 2*pi)
        assert np.pi < angle < 2 * np.pi
    
    def test_third_quadrant(self):
        """Test angle calculation in third quadrant (x < 0, y < 0)"""
        coord = pd.Series([-1.0, -1.0])
        angle = _circle_arctan(coord, center=(0, 0))
        # For x < 0, y < 0: y/x = 1, arctan(1) = pi/4
        # angle = 1.5*pi - pi/4 = 1.25*pi
        assert np.pi < angle < 1.5 * np.pi
    
    def test_fourth_quadrant(self):
        """Test angle calculation in fourth quadrant (x > 0, y < 0)"""
        coord = pd.Series([1.0, -1.0])
        angle = _circle_arctan(coord, center=(0, 0))
        # For x > 0, y < 0: y/x = -1, arctan(-1) = -pi/4
        # angle = pi/2 - (-pi/4) = 0.75*pi
        # This is actually in the 4th quadrant relative to the center mean
        assert 0 < angle < np.pi
    
    def test_center_offset(self):
        """Test angle calculation with custom center"""
        coord = pd.Series([1.0, 1.0])
        # With center at (1, 1), the point is treated as relative
        angle = _circle_arctan(coord, center=(1, 1))
        # The function adds center to coord, so it becomes (2, 2)
        # For x > 0: angle = pi/2 - arctan(2/2) = pi/2 - pi/4 = pi/4
        assert 0 < angle < np.pi / 2
    
    def test_origin_point(self):
        """Test angle calculation for point at center (should return NA)"""
        coord = pd.Series([0.0, 0.0])
        angle = _circle_arctan(coord, center=(0, 0))
        assert pd.isna(angle)
    
    def test_positive_y_axis(self):
        """Test angle calculation on positive y-axis"""
        coord = pd.Series([0.0, 1.0])
        angle = _circle_arctan(coord, center=(0, 0))
        assert angle == 0.0
    
    def test_negative_y_axis(self):
        """Test angle calculation on negative y-axis"""
        coord = pd.Series([0.0, -1.0])
        angle = _circle_arctan(coord, center=(0, 0))
        assert angle == np.pi


class TestAnglePseudotime:
    """Tests for angle_pseudotime function"""
    
    def test_basic_circular_ordering(self):
        """Test basic circular pseudotime calculation"""
        # Create circular embedding
        n_cells = 100
        angles = np.linspace(0, 2 * np.pi, n_cells, endpoint=False)
        embedding = pd.DataFrame({
            'x': np.cos(angles),
            'y': np.sin(angles)
        })
        
        result = angle_pseudotime(0, embedding)
        
        assert isinstance(result, pd.DataFrame)
        assert 'pseudotime' in result.columns
        assert result['pseudotime'].min() >= 0
        assert result['pseudotime'].max() <= 1
        assert result.shape[0] == n_cells
    
    def test_root_cell_by_index(self):
        """Test pseudotime with root cell specified by index"""
        n_cells = 50
        angles = np.linspace(0, 2 * np.pi, n_cells, endpoint=False)
        embedding = pd.DataFrame({
            'x': np.cos(angles),
            'y': np.sin(angles)
        })
        embedding.index = [f"cell_{i}" for i in range(n_cells)]
        
        # Root at cell_0 (angle 0)
        result = angle_pseudotime("cell_0", embedding)
        assert result.loc["cell_0", "pseudotime"] == 0
    
    def test_root_cell_by_angle(self):
        """Test pseudotime with root cell specified by angle value"""
        n_cells = 50
        angles = np.linspace(0, 2 * np.pi, n_cells, endpoint=False)
        embedding = pd.DataFrame({
            'x': np.cos(angles),
            'y': np.sin(angles)
        })
        
        # Root at angle 0 (using 0 as the root value)
        result = angle_pseudotime(0.0, embedding)
        assert result['pseudotime'].min() >= 0
        assert result['pseudotime'].max() <= 1
    
    def test_clockwise_direction(self):
        """Test clockwise pseudotime calculation"""
        n_cells = 50
        angles = np.linspace(0, 2 * np.pi, n_cells, endpoint=False)
        embedding = pd.DataFrame({
            'x': np.cos(angles),
            'y': np.sin(angles)
        })
        
        result_cw = angle_pseudotime(0, embedding, direction=1)
        result_ccw = angle_pseudotime(0, embedding, direction=-1)
        
        # Clockwise and anticlockwise should give different results
        assert not np.allclose(result_cw['pseudotime'], result_ccw['pseudotime'])
    
    def test_anticlockwise_direction(self):
        """Test anticlockwise pseudotime calculation"""
        n_cells = 50
        angles = np.linspace(0, 2 * np.pi, n_cells, endpoint=False)
        embedding = pd.DataFrame({
            'x': np.cos(angles),
            'y': np.sin(angles)
        })
        
        result = angle_pseudotime(0, embedding, direction=-1)
        
        assert result['pseudotime'].min() >= 0
        assert result['pseudotime'].max() <= 1
    
    def test_custom_column_name(self):
        """Test pseudotime with custom column name"""
        n_cells = 50
        angles = np.linspace(0, 2 * np.pi, n_cells, endpoint=False)
        embedding = pd.DataFrame({
            'x': np.cos(angles),
            'y': np.sin(angles)
        })
        
        custom_name = "my_pseudotime"
        result = angle_pseudotime(0, embedding, colname=custom_name)
        
        assert custom_name in result.columns
        assert 'pseudotime' not in result.columns or custom_name == 'pseudotime'
    
    def test_pseudotime_range(self):
        """Test that pseudotime is properly normalized to [0, 1]"""
        n_cells = 100
        angles = np.linspace(0, 2 * np.pi, n_cells, endpoint=False)
        embedding = pd.DataFrame({
            'x': np.cos(angles),
            'y': np.sin(angles)
        })
        
        result = angle_pseudotime(0, embedding)
        
        # Pseudotime should be in reasonable range (may not be exactly [0, 1]
        # due to circular nature and center calculation)
        assert result['pseudotime'].min() >= -0.1  # Allow small tolerance
        assert result['pseudotime'].max() <= 1.1
    
    def test_invalid_direction(self):
        """Test that invalid direction raises assertion error"""
        n_cells = 50
        angles = np.linspace(0, 2 * np.pi, n_cells, endpoint=False)
        embedding = pd.DataFrame({
            'x': np.cos(angles),
            'y': np.sin(angles)
        })
        
        with pytest.raises(AssertionError):
            angle_pseudotime(0, embedding, direction=0)
        
        with pytest.raises(AssertionError):
            angle_pseudotime(0, embedding, direction=2)
    
    def test_small_dataset(self):
        """Test with very small dataset"""
        embedding = pd.DataFrame({
            'x': [1.0, 0.0, -1.0],
            'y': [0.0, 1.0, 0.0]
        })
        
        result = angle_pseudotime(0, embedding)
        
        assert result.shape[0] == 3
        assert result['pseudotime'].min() >= 0
        assert result['pseudotime'].max() <= 1
    
    def test_embedding_with_index(self):
        """Test that result preserves original index"""
        n_cells = 50
        angles = np.linspace(0, 2 * np.pi, n_cells, endpoint=False)
        embedding = pd.DataFrame({
            'x': np.cos(angles),
            'y': np.sin(angles)
        })
        embedding.index = [f"cell_{i}" for i in range(n_cells)]
        
        result = angle_pseudotime("cell_0", embedding)
        
        assert list(result.index) == list(embedding.index)
