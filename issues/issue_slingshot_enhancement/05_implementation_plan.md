# 实施计划

## 概述

本文档详细说明如何将 SPG 项目增强的 Slingshot 模块整合到 `hic_basic` 项目中。

---

## 阶段划分

### Phase 1: 准备阶段 (1 天)

**目标**: 完成代码审查和准备工作

**任务**:
- [ ] 审查 `slingshot.py` 核心算法
- [ ] 审查 `plotter.py` 绘图功能
- [ ] 审查 `utils.py` 工具函数
- [ ] 确认代码风格符合 hic_basic 规范
- [ ] 创建 git 分支 `feature/slingshot-enhancement`

**负责人**: [待指定]  
**输出**: 代码审查报告

---

### Phase 2: 代码重构 (2 天)

**目标**: 调整代码结构，符合 hic_basic 规范

**任务**:
- [ ] 调整导入路径（使用相对导入）
- [ ] 统一文档字符串格式（中文 + 英文）
- [ ] 添加类型注解
- [ ] 优化错误处理
- [ ] 添加日志记录
- [ ] 代码格式化（black, isort）

**文件清单**:
```
hic_basic/hic_basic/slingshot/
├── __init__.py          # 模块导出
├── slingshot.py         # 主 Slingshot 类
├── lineage.py           # Lineage 类
├── plotter.py           # SlingshotPlotter 类
└── utils.py             # 工具函数
```

**负责人**: [待指定]  
**输出**: 重构后的代码

---

### Phase 3: 单元测试 (3 天)

**目标**: 编写完整的单元测试，覆盖率 > 80%

**任务**:
- [ ] 编写 `test_slingshot.py`
  - [ ] `test_init_with_anndata()`
  - [ ] `test_init_with_numpy()`
  - [ ] `test_string_labels()`
  - [ ] `test_predefined_tree()`
- [ ] 编写 `test_lineage.py`
  - [ ] `test_lineage_init()`
  - [ ] `test_lineage_iter()`
- [ ] 编写 `test_plotter.py`
  - [ ] `test_clusters_plot()`
  - [ ] `test_curves_plot()`
- [ ] 编写 `test_utils.py`
  - [ ] `test_mahalanobis()`
  - [ ] `test_infer_cluster_label_indices()`
- [ ] 编写 `test_predefined_tree.py`
  - [ ] `test_construct_tree_from_predefined()`
  - [ ] `test_missing_node_handling()`
- [ ] 编写 `test_string_labels.py`
  - [ ] `test_string_label_mapping()`
  - [ ] `test_list_lineages_with_strings()`

**测试文件结构**:
```
hic_basic/tests/test_slingshot/
├── __init__.py
├── test_slingshot.py
├── test_lineage.py
├── test_plotter.py
├── test_utils.py
├── test_predefined_tree.py
└── test_string_labels.py
```

**负责人**: [待指定]  
**输出**: 单元测试文件，覆盖率报告

---

### Phase 4: 集成测试 (2 天)

**目标**: 在真实数据上验证功能

**任务**:
- [ ] 准备测试数据集（scRNA-seq 数据）
- [ ] 运行完整流程测试
- [ ] 性能基准测试
  - [ ] 初始化时间
  - [ ] 拟合时间
  - [ ] 内存使用
- [ ] 边界条件测试
  - [ ] 空数据
  - [ ] 单细胞类型
  - [ ] 所有细胞类型缺失
- [ ] 与 pyslingshot 结果对比

**测试脚本**:
```python
# tests/integration/test_slingshot_integration.py
def test_spermatogenesis_trajectory():
    """测试精子发生轨迹推断"""
    pass

def test_bifurcation_trajectory():
    """测试分支轨迹推断"""
    pass

def test_missing_cell_types():
    """测试缺失细胞类型处理"""
    pass
```

**负责人**: [待指定]  
**输出**: 集成测试报告，性能基准报告

---

### Phase 5: 文档完善 (1 天)

**目标**: 编写完整的用户文档

**任务**:
- [ ] API 文档（Sphinx 格式）
- [ ] 使用示例
- [ ] 教程 notebook
- [ ] 更新 README
- [ ] 更新 CHANGELOG

**文档结构**:
```
hic_basic/docs/slingshot/
├── index.rst            # 主文档
├── api.rst              # API 文档
├── tutorial.ipynb       # 教程 notebook
├── examples.md          # 使用示例
└── faq.md               # 常见问题
```

**负责人**: [待指定]  
**输出**: 完整的用户文档

---

### Phase 6: 代码审查与合并 (1 天)

**目标**: 完成代码审查，合并到主分支

**任务**:
- [ ] 创建 Pull Request
- [ ] Code review
- [ ] 修复发现的问题
- [ ] 运行 CI/CD 测试
- [ ] 合并到 main 分支
- [ ] 打标签发布

**PR 模板**:
```markdown
## Description
整合 SPG 项目增强的 Slingshot 模块

## Changes
- 添加预定义树结构支持
- 添加字符串标签支持
- 添加缺失节点处理

## Testing
- 单元测试覆盖率 > 80%
- 集成测试通过
- 性能测试通过

## Checklist
- [ ] 代码审查通过
- [ ] 测试通过
- [ ] 文档完整
```

**负责人**: [待指定]  
**输出**: 合并后的代码，新版本发布

---

## 时间线

| 阶段 | 开始日期 | 结束日期 | 天数 |
|------|----------|----------|------|
| Phase 1: 准备 | Day 1 | Day 1 | 1 |
| Phase 2: 重构 | Day 2 | Day 3 | 2 |
| Phase 3: 单元测试 | Day 4 | Day 6 | 3 |
| Phase 4: 集成测试 | Day 7 | Day 8 | 2 |
| Phase 5: 文档 | Day 9 | Day 9 | 1 |
| Phase 6: 合并 | Day 10 | Day 10 | 1 |
| **总计** | | | **10 天** |

---

## 资源需求

### 人力资源

- **主要开发者**: 1-2 人
- **代码审查者**: 1 人
- **测试人员**: 1 人（可与开发者重叠）

### 计算资源

- **测试服务器**: 用于运行集成测试
- **CI/CD**: GitHub Actions 或类似服务

### 数据资源

- **测试数据**: scRNA-seq 数据集（至少 1 个）
- **参考数据**: 已知的细胞分化轨迹

---

## 风险管理

### 技术风险

| 风险 | 概率 | 影响 | 缓解措施 |
|------|------|------|----------|
| 代码风格不一致 | 中 | 低 | 提前审查，使用 linter |
| 测试覆盖率不足 | 中 | 中 | 设置覆盖率门槛 |
| 性能下降 | 低 | 中 | 性能基准测试 |
| 向后兼容性问题 | 低 | 高 | 充分的回归测试 |

### 进度风险

| 风险 | 概率 | 影响 | 缓解措施 |
|------|------|------|----------|
| 人员不足 | 中 | 中 | 提前协调资源 |
| 测试数据准备延迟 | 低 | 低 | 使用模拟数据 |
| Code review 延迟 | 中 | 中 | 提前预约审查者 |

---

## 成功标准

### 功能标准

- [ ] 预定义树结构正常工作
- [ ] 字符串标签自动映射
- [ ] 缺失节点正确处理
- [ ] 所有测试通过

### 质量标准

- [ ] 单元测试覆盖率 > 80%
- [ ] 代码审查通过
- [ ] 文档完整

### 性能标准

- [ ] 初始化时间 < 0.2s
- [ ] 拟合时间 (500 cells) < 10s
- [ ] 内存使用 < 100MB

---

## 后续计划

### 短期（1-3 个月）

- [ ] 收集用户反馈
- [ ] 修复发现的 bug
- [ ] 优化性能

### 中期（3-6 个月）

- [ ] 支持多谱系伪时间整合
- [ ] 添加差异轨迹分析
- [ ] 支持时间序列数据

### 长期（6-12 个月）

- [ ] 优化大规模数据性能
- [ ] 支持 GPU 加速
- [ ] 集成更多轨迹推断算法

---

## 联系方式

- **项目负责人**: [待指定]
- **技术支持**: [待指定]
- **问题反馈**: GitHub Issues

---

## 附录

### A. 文件清单

```
hic_basic/
├── hic_basic/
│   └── slingshot/
│       ├── __init__.py
│       ├── slingshot.py
│       ├── lineage.py
│       ├── plotter.py
│       └── utils.py
├── tests/
│   └── test_slingshot/
│       ├── __init__.py
│       ├── test_slingshot.py
│       ├── test_lineage.py
│       ├── test_plotter.py
│       ├── test_utils.py
│       ├── test_predefined_tree.py
│       └── test_string_labels.py
└── docs/
    └── slingshot/
        ├── index.rst
        ├── api.rst
        ├── tutorial.ipynb
        ├── examples.md
        └── faq.md
```

### B. 依赖清单

```python
# requirements.txt
numpy>=1.20.0
anndata>=0.7.0
scipy>=1.5.0
scikit-learn>=0.24.0
pcurvepy2>=1.0.0
tqdm>=4.50.0
matplotlib>=3.3.0
seaborn>=0.11.0
networkx>=2.5.0  # 可选，用于网络图
pygraphviz>=1.7.0  # 可选，用于 networkx 布局
```

### C. 测试命令

```bash
# 运行单元测试
pytest tests/test_slingshot/ -v

# 运行集成测试
pytest tests/integration/test_slingshot_integration.py -v

# 生成覆盖率报告
pytest tests/test_slingshot/ --cov=hic_basic.slingshot --cov-report=html
```
