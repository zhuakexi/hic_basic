"""
Lineage 类 - 表示细胞谱系
"""


class Lineage:
    """
    表示一个细胞谱系（一系列连续的细胞类型）
    """
    
    def __init__(self, clusters):
        """
        初始化 Lineage
        
        参数:
            clusters: 聚类索引列表
        """
        self.clusters = clusters
    
    def __len__(self):
        """返回谱系中的细胞类型数量"""
        return len(self.clusters)
    
    def __repr__(self):
        """返回谱系的字符串表示"""
        return 'Lineage' + str(self.clusters)
    
    def __iter__(self):
        """迭代谱系中的细胞类型"""
        for c in self.clusters:
            yield c
    
    @property
    def start(self):
        """返回起始细胞类型"""
        return self.clusters[0]
    
    @property
    def end(self):
        """返回终止细胞类型"""
        return self.clusters[-1]
