"""
Slingshot 工具函数
"""

import numpy as np


def mahalanobis(X1, X2, S1, S2):
    """
    计算两个向量之间的马氏距离
    
    参数:
        X1: 第一个向量 (n,)
        X2: 第二个向量 (n,)
        S1: 第一个协方差矩阵 (n, n)
        S2: 第二个协方差矩阵 (n, n)
    
    返回:
        马氏距离 (标量)
    """
    S_inv = np.linalg.inv(S1 + S2)
    diff = (X1 - X2).reshape(-1, 1)
    return np.matmul(np.matmul(diff.T, S_inv), diff)


def isint(x):
    """检查是否为整数"""
    if isinstance(x, int):
        return True
    if isinstance(x, str):
        return False
    try:
        a = float(x)
        b = int(a)
    except ValueError:
        return False
    else:
        return a == b


def isstr(x):
    """检查是否为字符串"""
    return isinstance(x, str)


def scale_to_range(x, a=0, b=1):
    """
    将数组缩放到指定范围
    
    参数:
        x: 输入数组
        a: 目标最小值
        b: 目标最大值
    
    返回:
        缩放后的数组
    """
    return ((x - x.min()) / (x.max() - x.min())) * (b - a) + a


def infer_cluster_label_indices(cluster_labels):
    """
    推断聚类标签索引
    
    参数:
        cluster_labels: 字符串标签列表或整数聚类分配列表
    
    返回:
        聚类分配的 numpy 数组
    """
    if isint(cluster_labels[0]):
        return cluster_labels
    elif isstr(cluster_labels[0]):
        # 将字符串标签转换为整数索引
        return np.array([np.where(np.unique(cluster_labels) == label)[0][0] for label in cluster_labels])
    else:
        raise ValueError("Unexpected cluster label dtype.")
