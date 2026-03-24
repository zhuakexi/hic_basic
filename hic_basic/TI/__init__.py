"""
Trajectory Inference (TI) 模块

包含多种轨迹推断方法：
- circle.angle_pseudotime: 圆形伪时间 baseline
- slingshot.Slingshot: Slingshot 轨迹推断
"""

from .circle import angle_pseudotime
from .slingshot import Slingshot, Lineage

__all__ = [
    'angle_pseudotime',
    'Slingshot',
    'Lineage',
]