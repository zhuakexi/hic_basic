# Issue: 吸收 SPG 项目增强的 Slingshot 模块

## 📋 Issue 概述

**标题**: 整合 SPG 项目增强的 Slingshot 轨迹推断模块  
**优先级**: High  
**标签**: `enhancement`, `pseudotime`, `trajectory`, `single-cell`  
**提出者**: SPG 项目团队  
**相关项目**: spg (精子发生单细胞 Hi-C 分析)

---

## 🎯 目标

将 SPG 项目中经过实际验证的 Slingshot 增强模块整合到 `hic_basic` 中，提供：
1. **预定义树结构支持** - 整合生物学先验知识
2. **字符串标签支持** - 使用细胞类型名称而非整数索引
3. **缺失节点处理** - 自动跳过缺失细胞类型

---

## 📦 交付内容

### 核心模块文件

```
hic_basic/hic_basic/slingshot/
├── __init__.py          (596 bytes)
├── slingshot.py         (28,202 bytes) - 主 Slingshot 类
├── lineage.py           (907 bytes) - Lineage 类
├── plotter.py           (5,957 bytes) - 绘图功能
└── utils.py             (1,689 bytes) - 工具函数
```

---

## 🔍 与 pyslingshot 的差异

| 特性 | pyslingshot (0.1.5) | hic_basic.slingshot (增强版) |
|------|---------------------|---------------------------|
| 代码行数 | 583 | ~700 (模块化后) |
| 预定义树结构 | ❌ | ✅ |
| 字符串标签 | ❌ | ✅ 自动检测 |
| 缺失节点处理 | ❌ | ✅ 祖孙直连 |
| 模块结构 | 单文件 | 5 个模块 |

### 核心 API 差异

```python
# pyslingshot - 仅支持整数索引
from pyslingshot import Slingshot
slingshot = Slingshot(
    data=adata,
    celltype_key='C_scANVI',
    start_node=0  # 必须使用整数
)

# hic_basic.slingshot - 支持字符串和预定义树
from hic_basic.slingshot import Slingshot
slingshot = Slingshot(
    data=adata,
    celltype_key='C_scANVI',
    predefined_tree={"SSC": ["Undiff.SPG"], ...},  # 新增
    start_node="SSC"  # 可直接使用字符串
)
```

---

## 🧪 测试方法

### 运行测试

```bash
cd /share/home/ychi/dev/hic_basic
./run_in_env.sh python -c "
from hic_basic.slingshot import Slingshot
from test_data import generate_test_data_spermatogenesis

# 生成测试数据
adata = generate_test_data_spermatogenesis(n_cells=520)

# 初始化 Slingshot
slingshot = Slingshot(
    data=adata,
    celltype_key='celltype',
    predefined_tree=adata.uns['predefined_tree'],
    start_node='SSC'
)

# 拟合
slingshot.fit(num_epochs=10)

# 验证
assert len(slingshot.lineages) > 0
print('✅ Test passed!')
"
```

---

## 📋 验收标准

### 功能完整性

- [ ] 预定义树结构正确构建
- [ ] 字符串标签自动映射
- [ ] 缺失节点正确跳过（祖孙直连）
- [ ] MST 构建正常
- [ ] 主曲线拟合收敛
- [ ] 伪时间计算正确
- [ ] 谱系列出功能正常

### 代码质量

- [ ] 模块化结构清晰
- [ ] 文档字符串完整（中文 + 英文）
- [ ] 类型注解正确
- [ ] 错误处理健壮
- [ ] 日志输出合理

### 向后兼容性

- [ ] 支持 pyslingshot 的所有现有用法
- [ ] 新增参数均为可选
- [ ] 不破坏现有 API

### 测试覆盖

- [ ] 单元测试覆盖率 > 80%
- [ ] 集成测试通过
- [ ] 在真实数据上验证

---

## 🚀 实施步骤

### Phase 1: 代码审查 (1-2 天)

- [ ] 审查 `slingshot.py` 核心算法
- [ ] 审查 `plotter.py` 绘图功能
- [ ] 审查 `utils.py` 工具函数
- [ ] 确认代码风格符合 hic_basic 规范

### Phase 2: 代码重构 (1-2 天)

- [ ] 调整导入路径
- [ ] 统一文档字符串格式
- [ ] 添加类型注解
- [ ] 优化错误处理

### Phase 3: 单元测试 (2-3 天)

- [ ] 编写 `test_slingshot.py`
- [ ] 编写 `test_lineage.py`
- [ ] 编写 `test_plotter.py`
- [ ] 编写 `test_utils.py`
- [ ] 编写 `test_predefined_tree.py`
- [ ] 编写 `test_string_labels.py`

### Phase 4: 集成测试 (1-2 天)

- [ ] 在真实 scRNA-seq 数据上测试
- [ ] 性能基准测试
- [ ] 边界条件测试
- [ ] 内存使用测试

### Phase 5: 文档完善 (1 天)

- [ ] API 文档
- [ ] 使用示例
- [ ] 教程 notebook
- [ ] 更新 README

### Phase 6: 代码合并 (1 天)

- [ ] Code review
- [ ] 修复发现的问题
- [ ] 合并到 main 分支
- [ ] 发布新版本

---

## 📚 使用示例

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
        "Diff.SPG": ["preL"],
        # ...
    },
    start_node="SSC"
)
slingshot.fit()

# 获取伪时间
pseudotime = slingshot.unified_pseudotime
```

### 处理缺失细胞类型

```python
# 如果数据中缺少"Undiff.SPG"
predefined_tree = {
    "SSC": ["Undiff.SPG"],
    "Undiff.SPG": ["Diff.SPG"],  # 自动跳过，连接 SSC→Diff.SPG
}

slingshot = Slingshot(
    data=adata,
    celltype_key='C_scANVI',
    predefined_tree=predefined_tree,
    start_node="SSC",
    debug_level='verbose'  # 输出警告信息
)
```

### 精子发生完整示例

```python
from hic_basic.slingshot import Slingshot

SPERMATOGENESIS_TREE = {
    "SSC": ["Undiff.SPG-1", "Undiff.SPG-2"],
    "Undiff.SPG-1": ["Diff.ing.SPG"],
    "Diff.ing.SPG": ["Diff.ed.SPG"],
    "Diff.ed.SPG": ["preL-1"],
    "preL-1": ["L1"],
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
    predefined_tree=SPERMATOGENESIS_TREE,
    start_node="SSC"
)
slingshot.fit(num_epochs=10)

# 查看谱系
slingshot.list_lineages()
# 输出：SSC, Undiff.SPG-1, Diff.ing.SPG, Diff.ed.SPG, preL-1, L1, L2, L3, preP, P, D, SPC7, S1, S2, S3, S4
```

---

## 🔗 相关资源

### 文档

- `01_comparison_with_pyslingshot.md` - 详细功能对比
- `02_api_documentation.md` - API 文档
- `03_usage_examples.md` - 使用示例
- `05_implementation_plan.md` - 详细实施计划

### 测试

- `04_test_data_and_functions.py` - 测试数据和绘图函数

### 外部参考

- Street et al. (2018). Slingshot: cell lineage and pseudotime inference for single-cell transcriptomics. BMC Genomics.
- pyslingshot GitHub: https://github.com/kstreet13/pyslingshot
- SPG project notebooks: `3.pseudotime.ipynb`, `3.1.pseudotime_slingshot.ipynb`

---

## 👥 参与人员

- **提案人**: SPG 项目团队
- **主要贡献**: 基于 spg 项目的 `lib/slingshot.py`
- **审核者**: [待指定]
- **实施者**: [待指定]

---

## 📅 时间线

- **Issue 创建**: 2026-03-23
- **计划开始**: [待定]
- **计划完成**: [待定]
- **预计发布时间**: hic_basic v1.x.0

---

## 📝 备注

1. 该模块已在 SPG 项目中经过实际验证
2. 特别适用于具有已知分化轨迹的生物学系统
3. 向后兼容 pyslingshot 的所有用法
4. 建议作为 hic_basic 的标准组件

---

## ✅ CheckList

- [x] Issue 目录创建
- [x] 核心模块文件准备
- [x] 文档准备
- [x] 测试数据准备
- [ ] 代码审查
- [ ] 单元测试编写
- [ ] 集成测试
- [ ] 文档完善
- [ ] 代码合并
