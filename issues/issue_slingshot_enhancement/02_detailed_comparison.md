# Slingshot 实现对比：`lib/slingshot.py` vs `pyslingshot`

## 概述

| 特性 | pyslingshot (外部包) | lib/slingshot.py (本地版本) |
|------|---------------------|---------------------------|
| 代码行数 | 583 行 | 1024 行 (+75%) |
| 版本 | 0.1.5 | 项目定制版 |
| 依赖 | `pcurvepy2` | `pcurvepy2` + 完整绘图模块 |

---

## 核心差异

### 1. **预定义树结构支持** (`predefined_tree`)

**pyslingshot**: ❌ 不支持

**lib/slingshot.py**: ✅ 完整支持
```python
def __init__(
    self,
    # ...
    predefined_tree: Optional[Dict[str, List[str]]] = None,  # 新增参数
):
```

**功能说明**：
- 允许用户指定已知的细胞分化树结构
- 自动处理缺失节点（祖细胞→孙细胞直连）
- 支持字符串和整数标签

**使用场景**：
```python
# 精子发生已知谱系树
predefined_tree = {
    "SSC": ["Undiff.SPG-1", "Undiff.SPG-2"],
    "Undiff.SPG-1": ["Diff.ing.SPG"],
    "Diff.ing.SPG": ["Diff.ed.SPG"],
    # ...
}
slingshot = Slingshot(
    data=adata,
    predefined_tree=predefined_tree,
    start_node="SSC"
)
```

---

### 2. **字符串标签支持** (`string_labels`)

**pyslingshot**: ❌ 仅支持整数索引

**lib/slingshot.py**: ✅ 自动检测并支持字符串标签

**实现细节**：
```python
# lib/slingshot.py 新增属性
self.string_labels = False          # 是否使用字符串标签
self.celltype_to_int = {}           # 字符串→整数映射
self.int_to_celltype = {}           # 整数→字符串映射
```

**自动映射逻辑**：
```python
if cluster_labels.dtype == object or isinstance(cluster_labels.iloc[0], str):
    self.string_labels = True
    unique_labels = sorted(cluster_labels.unique())
    self.celltype_to_int = {label: i for i, label in enumerate(unique_labels)}
    self.int_to_celltype = {i: label for label, i in self.celltype_to_int.items()}
```

---

### 3. **增强的树构建方法**

**pyslingshot**: 仅支持 MST（最小生成树）
```python
def tree(self):
    if self._tree is None:
        self.construct_mst(self.start_node)  # 只有 MST
    return self._tree
```

**lib/slingshot.py**: 支持 MST 和预定义树
```python
@property
def tree(self):
    if self._tree is None:
        if self.predefined_tree is not None:
            self.construct_tree_from_predefined()  # 新增
        else:
            self.construct_mst(self.start_node)
    return self._tree
```

**`construct_tree_from_predefined` 方法特点**：
- 递归处理缺失节点
- 支持多级跳过（祖细胞→孙细胞→曾孙细胞）
- 详细的警告输出
- 树结构验证

---

### 4. **新增属性：`label_tree`**

**lib/slingshot.py 独有**：
```python
@property
def label_tree(self):
    """返回带字符串标签的树结构"""
    if not hasattr(self, '_tree') or self._tree is None:
        return None
    
    if self.string_labels:
        string_tree = {}
        for int_parent, int_children in self._tree.items():
            parent_str = self.int_to_celltype[int_parent]
            children_str = [self.int_to_celltype[child] for child in int_children]
            string_tree[parent_str] = children_str
        return string_tree
    return self._tree
```

---

### 5. **完整的绘图模块**

**pyslingshot**: 基础绘图（分离的 `plotter.py` 模块）

**lib/slingshot.py**: 完整内嵌绘图功能
- `SlingshotPlotter` 类（~100 行）
- `generate_colormap()`: 自定义色轮生成
- `clusters()`: 支持 `color_mode='pseudotime'`
- `curves()`: 曲线绘制
- `network()`: 网络图（使用 networkx）

---

### 6. **新增工具函数**

**lib/slingshot.py 内嵌**（pyslingshot 中在单独的 `util.py`）：
```python
def mahalanobis(X1, X2, S1, S2): ...
def isint(x): ...
def isstr(x): ...
def scale_to_range(x, a=0, b=1): ...
def infer_cluster_label_indices(cluster_labels): ...
```

---

### 7. **新增方法：`list_lineages`**

```python
def list_lineages(self, cluster_to_label=None):
    """列出谱系（支持字符串标签）"""
    if cluster_to_label is None and self.string_labels:
        cluster_to_label = self.int_to_celltype
    
    for lineage in self.lineages:
        if cluster_to_label is not None:
            lineage_labels = [cluster_to_label[l] for l in lineage]
            print(', '.join(lineage_labels))
        else:
            print(', '.join([str(l) for l in lineage]))
```

---

## 核心算法对比（相同部分）

以下核心算法在两个版本中**完全相同**：

| 方法 | 功能 |
|------|------|
| `construct_mst()` | 基于马氏距离的 MST 构建 |
| `get_lineages()` | 从 MST 提取谱系 |
| `fit()` | 主曲线拟合循环 |
| `fit_lineage_curves()` | 谱系主曲线拟合 |
| `calculate_cell_weights()` | 细胞权重计算 |
| `avg_curves()` | 平均曲线计算 |
| `shrink_curves()` | 曲线收缩 |
| `shrink_branch_curves()` | 分支曲线收缩 |
| `shrinkage_percent()` | 收缩百分比计算 |
| `avg_branch_curves()` | 分支平均曲线 |

---

## 使用建议

### 使用 `pyslingshot` 的场景：
- 简单的探索性分析
- 不需要预定义树结构
- 整数标签足够

### 使用 `lib/slingshot.py` 的场景：
- ✅ **精子发生分析**（已知谱系树）
- 需要字符串细胞类型标签
- 需要可视化输出
- 需要处理缺失细胞类型
- 需要更详细的调试信息

---

## 演化历史推测

```
pyslingshot (原始 Python 实现)
     ↓
lib/slingshot.py (项目定制版)
     ├── 添加 predefined_tree 支持
     ├── 添加 string_labels 支持
     ├── 内嵌 util 函数
     ├── 扩展绘图功能
     └── 添加调试工具
```

**关键增强**：
1. **预定义树结构**：精子发生有明确的细胞分化路径，MST 可能产生不符合生物学拓扑的结构
2. **字符串标签**：直接使用细胞类型名称（如"SSC"、"preL-1"）而非整数索引
3. **缺失节点处理**：某些细胞类型可能在数据中缺失，需要自动跳过并连接祖孙节点

---

## 代码示例对比

### pyslingshot
```python
from pyslingshot import Slingshot

slingshot = Slingshot(
    data=adata,
    celltype_key='C_scANVI',
    start_node=0,  # 必须使用整数索引
    end_nodes=[5, 6]
)
slingshot.fit()
```

### lib/slingshot.py
```python
from lib.slingshot import Slingshot

# 定义已知的精子发生谱系树
predefined_tree = {
    "SSC": ["Undiff.SPG-1", "Undiff.SPG-2"],
    "Undiff.SPG-1": ["Diff.ing.SPG"],
    "Diff.ing.SPG": ["Diff.ed.SPG"],
    "Diff.ed.SPG": ["preL-1"],
    "preL-1": ["L1", "Z"],
    "L1": ["L2"],
    "L2": ["L3"],
    "L3": ["preP"],
    "preP": ["P"],
    "P": ["D"],
    "D": ["SPC7"],
    "SPC7": ["S1"],
    "S1": ["S2"],
    "S2": ["S3"],
    "S3": ["S4"]
}

slingshot = Slingshot(
    data=adata,
    celltype_key='C_scANVI',
    predefined_tree=predefined_tree,  # 使用预定义树
    start_node="SSC"  # 可直接使用字符串
)
slingshot.fit()

# 查看谱系
slingshot.list_lineages()
# 输出：SSC, Undiff.SPG-1, Diff.ing.SPG, ...
```

---

## 总结

`lib/slingshot.py` 是 `pyslingshot` 的**超集**，主要增强：

| 增强点 | 价值 |
|--------|------|
| `predefined_tree` | 整合生物学先验知识 |
| `string_labels` | 提高代码可读性 |
| 缺失节点处理 | 增强鲁棒性 |
| 内嵌绘图 | 减少依赖 |
| 调试工具 | 便于问题排查 |

**建议**：本项目继续使用 `lib/slingshot.py`，因其专为精子发生轨迹推断定制。
