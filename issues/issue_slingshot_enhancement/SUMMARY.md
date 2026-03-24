# Issue 包总结

## 📦 这是什么？

这是一个完整的 **Issue 包**，用于向 `hic_basic` 项目提案：吸收 SPG 项目中经过实际验证的 Slingshot 增强模块。

---

## 📁 文件清单

```
hic_basic_issue/
├── 00_README.md                    # Issue 主文档（必读）
├── 01_comparison_with_pyslingshot.md    # 与 pyslingshot 对比
├── 02_detailed_comparison.md       # 详细功能对比
├── 03_usage_examples.md            # 使用示例
├── 04_test_data_and_functions.py   # 测试数据和绘图函数
└── 05_implementation_plan.md       # 详细实施计划
```

---

## ✅ 我明白了什么

### 1. 核心需求

将 SPG 项目的 `lib/slingshot.py` 整合到 `hic_basic` 中，提供三个关键增强：

| 增强点 | 说明 | 价值 |
|--------|------|------|
| **预定义树结构** | 支持 `predefined_tree` 参数 | 整合生物学先验知识 |
| **字符串标签** | 自动检测并映射字符串标签 | 提高代码可读性 |
| **缺失节点处理** | 自动跳过缺失细胞类型 | 增强鲁棒性 |

### 2. 与 pyslingshot 的关系

- `hic_basic.slingshot` 是 `pyslingshot` 的**超集**
- **向后兼容** pyslingshot 的所有用法
- 新增功能均为**可选参数**
- 核心算法**完全相同**

### 3. 代码结构

```
hic_basic/slingshot/
├── __init__.py      # 模块导出
├── slingshot.py     # 主类 (~700 行)
├── lineage.py       # Lineage 类
├── plotter.py       # 绘图功能
└── utils.py         # 工具函数
```

### 4. 测试方法

提供了完整的测试脚本：
- `generate_test_data_spermatogenesis()` - 精子发生测试数据
- `generate_test_data_bifurcation()` - 分支轨迹测试数据
- `plot_slingshot_result()` - 4 合 1 绘图函数
- `plot_lineage_network()` - networkx 网络图
- `run_test_pipeline()` - 完整测试流程

### 5. 实施路径

**6 个阶段，10 天完成**:
1. Phase 1: 准备 (1 天) - 代码审查
2. Phase 2: 重构 (2 天) - 代码调整
3. Phase 3: 单元测试 (3 天) - 覆盖率 > 80%
4. Phase 4: 集成测试 (2 天) - 真实数据验证
5. Phase 5: 文档 (1 天) - 用户文档
6. Phase 6: 合并 (1 天) - Code review + 发布

---

## 📋 我打算怎么做

### 步骤 1: 审查现有代码

- [x] 阅读 `lib/slingshot.py` 源码
- [x] 对比 `pyslingshot` 源码
- [x] 识别核心差异
- [x] 确认功能完整性

### 步骤 2: 准备 Issue 文档

- [x] 创建 Issue 目录结构
- [x] 编写 Issue 主文档 (`00_README.md`)
- [x] 编写功能对比文档 (`01_comparison.md`)
- [x] 准备测试数据 (`04_test_data_and_functions.py`)
- [x] 编写实施计划 (`05_implementation_plan.md`)

### 步骤 3: 准备代码包

- [x] 模块化拆分 `slingshot.py`
  - [x] 提取 `utils.py` - 工具函数
  - [x] 提取 `lineage.py` - Lineage 类
  - [x] 提取 `plotter.py` - 绘图功能
  - [x] 保留 `slingshot.py` - 主类
- [x] 创建 `__init__.py` - 模块导出
- [x] 添加类型注解
- [x] 编写文档字符串

### 步骤 4: 准备测试

- [x] 编写测试数据生成函数
- [x] 编写绘图函数
- [x] 编写完整测试流程
- [ ] 运行测试验证（需要 hic_basic 环境）

### 步骤 5: 提交 Issue

**下一步行动**:

1. **复制文件到 hic_basic**:
   ```bash
   # 将整个 issue 包复制到 hic_basic
   cp -r /share/home/ychi/dev/spg/docs/hic_basic_issue \
         /share/home/ychi/dev/hic_basic/issues/issue_slingshot_enhancement/
   
   # 复制代码包
   cp -r /share/home/ychi/dev/hic_basic/hic_basic/slingshot \
         /share/home/ychi/dev/hic_basic/hic_basic/slingshot_backup/
   ```

2. **在 hic_basic 创建 Issue**:
   - GitHub Issues 中创建新 Issue
   - 标题：`[Enhancement] 整合 SPG 项目增强的 Slingshot 模块`
   - 标签：`enhancement`, `pseudotime`, `trajectory`
   - 关联此 issue 包

3. **通知相关人员**:
   - @hic_basic 维护者
   - @spg 项目团队

---

## 🎯 预期结果

### 短期（1-2 周）

- Issue 被接受
- 开始代码审查
- 确定实施时间表

### 中期（1-2 个月）

- 代码合并到 hic_basic
- 发布新版本（v1.x.0）
- 文档上线

### 长期（3-6 个月）

- 在更多项目中使用
- 收集用户反馈
- 持续优化

---

## 📞 联系方式

如有问题，请联系：

- **Issue 提出者**: SPG 项目团队
- **主要贡献**: 基于 spg 项目的 `lib/slingshot.py`
- **应用场景**: 精子发生轨迹推断、HiRES 数据分析

---

## 🔗 相关资源

### 文档

- `00_README.md` - Issue 主文档
- `01_comparison_with_pyslingshot.md` - 功能对比
- `05_implementation_plan.md` - 实施计划

### 代码

- `hic_basic/hic_basic/slingshot/` - 模块代码
- `04_test_data_and_functions.py` - 测试数据

### 外部参考

- pyslingshot: https://github.com/kstreet13/pyslingshot
- Slingshot paper: Street et al. (2018) BMC Genomics
- SPG notebooks: `3.pseudotime.ipynb`, `3.1.pseudotime_slingshot.ipynb`

---

## ✨ 总结

这个 Issue 包包含：

1. ✅ **完整的 Issue 文档** - 说明需求、对比、实施计划
2. ✅ **模块化代码** - 5 个模块，结构清晰
3. ✅ **测试数据和函数** - 可直接运行验证
4. ✅ **详细对比** - 与 pyslingshot 的差异说明
5. ✅ **实施计划** - 6 个阶段，10 天完成

**目标**: 让 hic_basic 项目吸收这个增强的 Slingshot 模块，使所有用户受益。

---

*最后更新：2026-03-23*
