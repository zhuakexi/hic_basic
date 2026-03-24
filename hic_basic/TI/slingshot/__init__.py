"""
Slingshot 轨迹推断模块

基于 pyslingshot 的增强版本，支持：
- 预定义树结构 (predefined_tree)
- 字符串标签 (string_labels)
- 缺失节点自动处理
"""

from .slingshot import Slingshot
from .lineage import Lineage
from .plot import plot_slingshot, plot_clusters, plot_curves
from .utils import (
    mahalanobis,
    isint,
    isstr,
    scale_to_range,
    infer_cluster_label_indices
)

__all__ = [
    'Slingshot',
    'Lineage',
    'plot_slingshot',
    'plot_clusters',
    'plot_curves',
    'mahalanobis',
    'isint',
    'isstr',
    'scale_to_range',
    'infer_cluster_label_indices',
]

__version__ = '1.0.0'
