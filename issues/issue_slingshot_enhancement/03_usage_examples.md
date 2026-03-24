# PR: 添加增强的 Slingshot 轨迹推断模块

## 概述

请求将增强的 Slingshot 轨迹推断模块添加到 `hic_basic` 包中。该模块基于 `pyslingshot` 包进行了重要扩展，特别适用于具有已知分化轨迹的生物学系统（如精子发生）。

---

## 动机

### 当前问题

1. **`pyslingshot` 功能限制**：
   - 仅支持整数索引，不支持字符串细胞类型标签
   - 不支持预定义的树结构
   - 无法处理缺失的细胞类型

2. **生物学需求**：
   - 精子发生具有明确的细胞分化路径
   - 某些细胞类型可能在数据中缺失
   - 需要使用生物学名称（如"SSC"、"preL"）而非整数索引

### 解决方案

将 `lib/slingshot.py` 中的增强功能整合到 `hic_basic` 中，提供：
- 预定义树结构支持
- 字符串标签自动映射
- 缺失节点智能处理

---

## 与 pyslingshot 的差异

| 特性 | pyslingshot (0.1.5) | 本模块 (增强版) |
|------|---------------------|-----------------|
| 代码行数 | 583 行 | 1024 行 |
| 预定义树结构 | ❌ | ✅ |
| 字符串标签 | ❌ | ✅ 自动检测 |
| 缺失节点处理 | ❌ | ✅ 祖孙直连 |
| 绘图功能 | 基础 | 完整内嵌 |
| 调试工具 | 少 | 详细日志 |

### 核心增强点

#### 1. 预定义树结构 (`predefined_tree`)

```python
# 精子发生已知谱系树
predefined_tree = {
    "SSC": ["Undiff.SPG-1", "Undiff.SPG-2"],
    "Undiff.SPG-1": ["Diff.ing.SPG"],
    "Diff.ing.SPG": ["preL-1"],
    "preL-1": ["L1", "Z"],  # 分支
    # ...
}

slingshot = Slingshot(
    data=adata,
    predefined_tree=predefined_tree,  # 新增参数
    start_node="SSC"  # 可直接使用字符串
)
```

#### 2. 字符串标签支持

```python
# 自动检测字符串标签并创建映射
slingshot.string_labels = True
slingshot.celltype_to_int = {"SSC": 0, "Undiff.SPG": 1, ...}
slingshot.int_to_celltype = {0: "SSC", 1: "Undiff.SPG", ...}
```

#### 3. 缺失节点处理

```python
# 如果"Undiff.SPG"缺失，自动连接"SSC"→"Diff.ing.SPG"
# 递归处理多级缺失
```

---

## 拟议的文件结构

```
hic_basic/
└── hic_basic/
    └── slingshot/
        ├── __init__.py          # 导出 Slingshot 类
        ├── slingshot.py         # 主模块（1024 行）
        ├── lineage.py           # Lineage 类（可复用）
        ├── plotter.py           # SlingshotPlotter（可复用）
        └── utils.py             # 工具函数（可复用）
```

---

## API 设计

### 类：`Slingshot`

```python
class Slingshot:
    def __init__(
        self,
        data: Union[AnnData, np.ndarray],
        cluster_labels_onehot: Optional[np.ndarray] = None,
        celltype_key: Optional[str] = None,
        obsm_key: str = 'X_umap',
        start_node: Union[str, int] = 0,
        end_nodes: Optional[List[Union[str, int]]] = None,
        predefined_tree: Optional[Dict[str, List[str]]] = None,  # 新增
        debug_level: Optional[int] = None
    )
    
    def fit(self, num_epochs: int = 10) -> None
    def get_lineages(self) -> List[Lineage]
    def list_lineages(self, cluster_to_label: Optional[Dict] = None) -> None
    def save_params(self, filepath: str) -> None
    def load_params(self, filepath: str) -> None
    
    @property
    def tree(self) -> Dict[int, List[int]]
    @property
    def label_tree(self) -> Union[Dict[int, List[int]], Dict[str, List[str]]]
    @property
    def lineages(self) -> List[Lineage]
    @property
    def curves(self) -> List[PrincipalCurve]
    @property
    def unified_pseudotime(self) -> np.ndarray
```

---

## 使用示例

### 基础用法

```python
from hic_basic.slingshot import Slingshot

# 使用预定义树
slingshot = Slingshot(
    data=adata,
    celltype_key='C_scANVI',
    predefined_tree={
        "SSC": ["Undiff.SPG"],
        "Undiff.SPG": ["Diff.SPG"],
        # ...
    },
    start_node="SSC"
)
slingshot.fit()

# 获取伪时间
pseudotime = slingshot.unified_pseudotime
```

### 高级用法

```python
# 处理缺失细胞类型
predefined_tree = {
    "SSC": ["Undiff.SPG-1", "Undiff.SPG-2"],
    "Undiff.SPG-1": ["Diff.ing.SPG"],
    # 如果"Diff.ing.SPG"缺失，自动连接"Undiff.SPG-1"→"preL-1"
}

slingshot = Slingshot(
    data=adata,
    celltype_key='C_scANVI',
    predefined_tree=predefined_tree,
    start_node="SSC",
    debug_level='verbose'  # 输出详细日志
)
slingshot.fit()

# 查看谱系
slingshot.list_lineages()
# 输出：SSC, Undiff.SPG-1, Diff.ing.SPG, preL-1, ...

# 保存/加载参数
slingshot.save_params('slingshot_params.npy')
slingshot.load_params('slingshot_params.npy')
```

---

## 测试数据

提供了两个测试数据集：

### 1. 精子发生轨迹 (`generate_test_data_spermatogenesis`)

- 模拟完整的精子发生过程
- 14 个细胞类型
- 弧形 UMAP 结构

### 2. 分支轨迹 (`generate_test_data_bifurcation`)

- 模拟 L/Z 分支
- 12 个细胞类型
- Y 形 UMAP 结构

### 测试脚本

```bash
# 运行测试
python test_slingshot_hic_basic.py

# 输出：
# - /tmp/slingshot_test_spermatogenesis_main.png
# - /tmp/slingshot_test_spermatogenesis_network.png
```

---

## 绘图功能

### `plot_slingshot_result()`

生成 4 个子图：
1. 细胞类型分布
2. 聚类中心和起始节点
3. 主曲线和谱系
4. 伪时间分布

### `plot_lineage_network()`

使用 networkx 绘制谱系网络图（hierarchical 布局）

---

## 依赖

```python
# 核心依赖
numpy
anndata
scipy
scikit-learn
pcurvepy2  # 主曲线算法
tqdm

# 可选依赖（绘图）
matplotlib
seaborn
networkx  # 网络图
pygraphviz  # networkx 布局
```

---

## 实现计划

### Phase 1: 代码重构（1-2 天）

- [ ] 拆分 `lib/slingshot.py` 为独立模块
- [ ] 提取可复用的 `Lineage` 和 `SlingshotPlotter` 类
- [ ] 整理工具函数到 `utils.py`

### Phase 2: 单元测试（1-2 天）

- [ ] 测试字符串标签映射
- [ ] 测试预定义树构建
- [ ] 测试缺失节点处理
- [ ] 测试主曲线拟合

### Phase 3: 文档（1 天）

- [ ] API 文档
- [ ] 使用示例
- [ ] 与 pyslingshot 的对比文档

### Phase 4: 集成测试（1 天）

- [ ] 在真实数据上验证
- [ ] 性能基准测试
- [ ] 边界条件测试

---

## 验收标准

1. **功能完整性**：
   - [ ] 预定义树结构正常工作
   - [ ] 字符串标签自动映射
   - [ ] 缺失节点正确处理

2. **向后兼容**：
   - [ ] 支持 `pyslingshot` 的所有现有用法
   - [ ] 新增参数均为可选

3. **测试覆盖**：
   - [ ] 单元测试覆盖率 > 80%
   - [ ] 集成测试通过

4. **文档完整**：
   - [ ] API 文档完整
   - [ ] 使用示例清晰
   - [ ] 与 pyslingshot 的差异说明

---

## 相关文件

- `docs/hic_basic_slingshot_README.md` - 模块文档
- `docs/test_slingshot_hic_basic.py` - 测试数据和绘图函数
- `docs/slingshot_comparison.md` - 详细对比文档

---

## 维护者

- **主要作者**: [Your Name]
- **审核者**: [待指定]
- **预计合并时间**: [待定]

---

## 后续计划

1. 支持多谱系伪时间整合
2. 添加差异轨迹分析
3. 支持时间序列数据
4. 优化大规模数据性能

---

## 参考资料

- Street et al. (2018). pyslingshot: Python implementation of Slingshot.
- Street et al. (2018). Slingshot: cell lineage and pseudotime inference for single-cell transcriptomics. BMC Genomics.
- 本项目 notebooks: `3.pseudotime.ipynb`, `3.1.pseudotime_slingshot.ipynb`, `3.slingshot_single_lineage.ipynb`
