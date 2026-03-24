# Lu 2012 Crossover Detection - 项目文件说明

## 目录结构

```
lu_2012_crossover_detection/
├── README.md                      # 主要文档：使用说明和 API
├── ALGORITHM.md                   # 算法详解：数学原理和实现细节
├── crossover_detection.py         # 主程序：Crossover 检测实现
├── test_crossover_detection.py    # 测试脚本：验证算法正确性
├── example_data.tsv               # 示例数据：输入输出格式说明
├── requirements.txt               # 依赖：Python 包列表
└── test_output/                   # 测试输出目录
    ├── test_snps.tsv              # 测试 SNP 数据
    ├── test_phasing.tsv           # 测试定相数据
    ├── TEST001_haplotype_blocks.tsv   # 单倍型区块结果
    ├── TEST001_crossovers.tsv         # Crossover 位置结果
    └── TEST001_statistics.tsv         # 统计信息结果
```

## 文件说明

### 1. `README.md`
**用途:** 主要使用文档

**内容:**
- 论文信息
- 算法概述流程图
- 输入输出数据格式
- 使用示例（命令行和 Python API）
- 参数说明
- 性能指标
- 注意事项

**何时使用:** 首次使用时阅读，了解如何使用该工具

---

### 2. `ALGORITHM.md`
**用途:** 算法技术文档

**内容:**
- 两步法详细流程
- 数学公式推导（HMM、Viterbi 算法）
- 数据结构定义
- 关键函数实现细节
- 参数优化建议
- 结果解读指南

**何时使用:** 需要深入理解算法原理或调试问题时

---

### 3. `crossover_detection.py`
**用途:** 主程序实现

**核心类:**
- `HMMCrossoverDetector`: HMM 模型实现
  - `_rough_estimation()`: 粗略估计 crossover 位置
  - `_viterbi()`: Viterbi 算法实现
  - `_encode_observations()`: 基因型编码
  - `_filter_small_fragments()`: 小片段过滤
  - `detect()`: 完整检测流程

- `CrossoverDetector`: 完整流程封装
  - `load_snp_data()`: 加载 SNP 数据
  - `load_phasing_data()`: 加载定相数据
  - `detect_crossovers()`: 执行检测
  - `save_results()`: 保存结果

**何时使用:** 直接运行或作为 Python 模块导入

---

### 4. `test_crossover_detection.py`
**用途:** 测试和验证

**功能:**
- `generate_synthetic_data()`: 生成带已知 crossover 的合成数据
- `test_hmm_detector()`: 测试 HMM 检测器性能
- `test_full_pipeline()`: 测试完整流程
- 性能指标计算（Precision, Recall, F1）

**测试结果:**
```
HMM Detector:
  Precision: 0.933
  Recall: 0.700
  F1 Score: 0.800

Full Pipeline:
  Total crossovers detected: 12
  Total haplotype blocks: 15
  Mean confidence: 1.000
```

**何时使用:** 安装后验证安装是否正确，或修改代码后回归测试

---

### 5. `example_data.tsv`
**用途:** 数据格式示例

**内容:**
- `sperm_snps.tsv`: SNP 基因型数据格式
- `phasing_info.tsv`: 单倍型定相数据格式
- `expected_output_crossovers.tsv`: 预期输出格式
- `expected_output_statistics.tsv`: 统计信息格式
- 数据格式注释和 QC 说明

**何时使用:** 准备输入数据时参考

---

### 6. `requirements.txt`
**用途:** Python 依赖列表

**依赖:**
```
numpy>=1.20.0      # 数值计算
pandas>=1.3.0      # 数据处理
scipy>=1.7.0       # 科学计算
hmmlearn>=0.2.8    # HMM 模型（可选，本实现使用自定义 HMM）
```

**安装:**
```bash
pip install -r requirements.txt
```

---

## 快速开始

### 1. 安装依赖
```bash
cd lu_2012_crossover_detection
pip install -r requirements.txt
```

### 2. 运行测试
```bash
python test_crossover_detection.py
```

### 3. 分析自己的数据
```bash
python crossover_detection.py \
    --input your_data/snps.tsv \
    --phasing your_data/phasing.tsv \
    --output results/
```

### 4. 查看结果
```bash
# 查看 crossover 位置
cat results/your_sperm_crossovers.tsv

# 查看统计信息
cat results/your_sperm_statistics.tsv
```

---

## 与其他工具的对比

| 特性 | Lu et al. (2012) 实现 | 其他工具 |
|------|----------------------|----------|
| 输入数据 | 单精子 WGS (低覆盖度) | 家系 trio/群体数据 |
| 算法 | HMM + Viterbi | 隐马尔可夫/最大似然 |
| 分辨率 | ~30 kb | 1-100 kb 不等 |
| 单细胞级别 | ✓ | 通常否 |
| 同时检测非整倍体 | ✓ | 通常否 |

---

## 常见问题

### Q1: 需要多少 SNP 才能可靠检测？
**A:** 至少 10,000 个均匀分布的 SNP。Lu et al. (2012) 使用了约 100,000-500,000 个 SNP。

### Q2: 可以检测基因转换 (gene conversion) 吗？
**A:** 当前实现过滤了<300kb 的片段（可能是基因转换）。如需检测基因转换，需要降低 `min_fragment_size` 参数并增加测序深度。

### Q3: 适用于长读长测序数据吗？
**A:** 可以。长读长数据（如 ONT, PacBio）可以提供更好的定相信息和更高的分辨率。参考 Xie et al. (2023) NAR 的实现。

### Q4: 如何处理性染色体？
**A:** X/Y 染色体需要单独处理，因为精子只携带其中一条。PAR 区域（假常染色体区域）会重组，需要特殊处理。

---

## 版本信息

- **代码版本:** 1.0.0
- **基于论文:** Lu et al. (2012) Science 338:1627-1630
- **创建日期:** 2026-03-23
- **Python 版本:** >= 3.8

---

## 许可证

本实现仅用于学术研究和教育目的。

## 致谢

算法基于 Lu et al. (2012) 的补充方法部分实现。
