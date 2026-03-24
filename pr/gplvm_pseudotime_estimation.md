# Pseudotime 估计：GPLVM 方法详解

## 一、算法目的 (Purpose)

### 1.1 核心问题

在单细胞轨迹分析中，**Pseudotime（伪时间）** 是一个关键概念：
- 它将离散的细胞按照生物学进程（如分化、发育）排序
- 每个细胞被赋予一个连续的时间值，表示其在生物过程中的相对位置
- 这是进行轨迹插值和平滑的**基础输入**

**为什么需要 GPLVM 估计 pseudotime？**

1. **没有真实时间标签**：单细胞实验通常是终点测量（endpoint assay），无法直接获取每个细胞的真实时间
2. **高维非线性**：基因表达是高维（数千基因）且非线性的过程
3. **需要降维**：将高维基因表达压缩到 1-2 维的"时间"表示
4. **不确定性建模**：概率方法可以提供估计的置信度

### 1.2 GPLVM 的作用

GPLVM (Gaussian Process Latent Variable Model) 在此框架中用于：

```
高维基因表达空间                    低维潜在空间
(数千基因 × 数千细胞)     ---->      (pseudotime × 细胞)
                                    
                                    - 第一维：pseudotime（生物学进程）
                                    - 第二维：其他潜在效应（如 batch）
```

**关键优势**：
- **概率建模**：提供不确定性估计
- **时间先验整合**：可以 Incorporate 已知的采样时间信息
- **非线性**：捕捉复杂的发育轨迹
- **平滑性**：高斯过程保证连续的 pseudotime 映射

---

## 二、在 G2G 文章中的作用 (Role in the Paper)

### 2.1 在整体流程中的位置

```
原始 scRNA-seq 数据
       ↓
预处理（标准化、过滤、HVG 选择）
       ↓
Batch 校正（scVI）
       ↓
【GPLVM Pseudotime 估计】← 本方法
       ↓
细胞按 pseudotime 排序
       ↓
高斯过程插值（在 pseudotime 上）
       ↓
轨迹对齐（G2G 核心算法）
```

### 2.2 具体应用场景

在 G2G 论文中，GPLVM 用于三个数据集的 pseudotime 估计：

| 数据集 | 先验类型 | 说明 |
|--------|----------|------|
| **ATO（人工胸腺类器官）** | 细胞采样天数 | 有真实时间信息（day -7, 0, week 1-7） |
| **T1（胎儿 T 细胞参考）** | kNN 近似时间 | 无真实时间，用参考数据集近似 |
| **CD8+（胎儿 CD8+ T 细胞）** | kNN 近似时间 | 无真实时间，用参考数据集近似 |

### 2.3 为什么这对插值至关重要

1. **插值需要时间轴**：高斯过程插值在 pseudotime 上进行，pseudotime 质量直接影响插值结果
2. **等间距插值点依赖 pseudotime 分布**：optimal binning 基于 pseudotime 分布
3. **跨数据集对齐的基础**：两个轨迹的 pseudotime 都需要归一化到 [0,1] 才能对齐

---

## 三、算法详细步骤 (Algorithm Details)

### 3.1 GPLVM 模型定义

**数学公式**：

```
给定：
- Y: 基因表达矩阵 (n_cells × n_genes)
- X: 潜在变量 (n_cells × n_latent)，即 pseudotime

GPLVM 假设：
Y = f(X) + ε

其中：
- f(·) ~ GP(0, K)：高斯过程映射
- K: 协方差核函数（如 RBF）
- ε ~ N(0, σ²I)：观测噪声

目标：
最大化边缘似然，学习 X 和 GP 超参数
```

### 3.2 完整流程

#### 步骤 1：数据准备

```
输入：
- 表达矩阵：expr[cell, gene]（已标准化、log1p 转换）
- 可选：细胞采样时间 capture_time[cell]（如果有）
- 可选：batch 信息 batch[cell]

预处理：
1. 过滤低质量细胞和基因
2. 选择高变基因（如 2,608 个基因用于 T1 数据集）
3. 标准化到零均值、单位方差
4. （可选）使用 scVI 进行 batch 校正
```

#### 步骤 2：初始化潜在变量 X

**情况 A：有采样时间信息（如 ATO 数据）**

```python
# 直接使用采样天数作为 pseudotime 先验的均值
# 例如：capture_days = [-7, 0, 7, 14, 21, 28, 35]
# 归一化到 [0, 1]
pseudotime_prior = (capture_days - min(capture_days)) / 
                   (max(capture_days) - min(capture_days))

# 初始化：X[:, 0] = pseudotime_prior
#         X[:, 1] = 0（第二维潜在变量零初始化）
```

**情况 B：无采样时间信息（如胎儿参考数据）**

```python
# 使用 kNN 从参考数据集近似时间
# 1. 在 UMAP 空间找到每个细胞的 k=3 个最近邻（来自有时间的参考数据）
# 2. 加权平均邻居的 capture time
# 3. 按细胞类型进一步细化（使用细胞类型特异性平均时间）

# 伪代码：
for each cell in reference:
    neighbors = find_knn(cell, organoid_cells, k=3, space='UMAP')
    weights = softmax(1 / distance[cell, neighbors])
    pseudotime[cell] = weighted_average(capture_time[neighbors], weights)

# 按细胞类型细化
for each cell_type:
    avg_time = mean(pseudotime[cells of this type])
    pseudotime[cells of this type] = avg_time

# 归一化到 [0, 1]
pseudotime = (pseudotime - min(pseudotime)) / 
             (max(pseudotime) - min(pseudotime))
```

#### 步骤 3：定义 GPLVM 模型

```python
# 使用 Pyro (v.1.8.0+) 实现
import pyro
import pyro.contrib.gp as gp

# 模型参数
n_latent = 2  # 2 维潜在空间
n_inducing = 32  # 稀疏 GP 的诱导点数

# 核函数：RBF
kernel = gp.kernels.RBF(input_dim=n_latent)

# GPLVM 模型
gplvm = gp.models.SparseGPRegression(
    X=X_init,  # 初始潜在变量
    Y=Y_data,  # 基因表达
    kernel=kernel,
    X_inducing=inducing_points,  # 诱导点
    likelihood=gaussian_likelihood
)
```

#### 步骤 4：优化

```python
# 优化器：Adam
optimizer = torch.optim.Adam(
    gplvm.parameters(),
    lr=0.01
)

# 训练循环
for iteration in range(2000):
    optimizer.zero_grad()
    loss = -gplvm.log_marginal_likelihood()
    loss.backward()
    optimizer.step()
    
    # 监控收敛
    if iteration % 100 == 0:
        print(f"Iteration {iteration}: Loss = {loss.item()}")

# 典型收敛：2000 次迭代
```

#### 步骤 5：提取 Pseudotime

```python
# 优化后，X 的第一维即为 pseudotime
pseudotime_raw = gplvm.X[:, 0].detach().numpy()

# 归一化到 [0, 1]
pseudotime = (pseudotime_raw - pseudotime_raw.min()) / 
             (pseudotime_raw.max() - pseudotime_raw.min())
```

#### 步骤 6：异常值校正

```python
# 按细胞类型检查并校正异常值
for each cell_type:
    pt_values = pseudotime[cells of this type]
    Q1 = np.percentile(pt_values, 25)
    Q3 = np.percentile(pt_values, 75)
    IQR = Q3 - Q1
    
    # IQR 规则识别异常值
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR
    outliers = (pt_values < lower_bound) | (pt_values > upper_bound)
    
    # 用细胞类型平均值替换异常值
    pseudotime[outliers] = np.mean(pt_values[~outliers])
```

### 3.3 关键参数设置

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `n_latent` | 2 | 潜在空间维度（第 1 维=pseudotime，第 2 维=batch 等） |
| `n_inducing` | 32 | 稀疏 GP 诱导点数 |
| `kernel` | RBF | 协方差核函数 |
| `optimizer` | Adam | 优化器 |
| `n_iterations` | 2000 | 优化迭代次数 |
| `learning_rate` | 0.01 | 学习率 |

---

## 四、为什么选择 GPLVM 而非其他方法

### 4.1 与其他 Pseudotime 方法的比较

| 方法 | 优点 | 缺点 | GPLVM 优势 |
|------|------|------|-----------|
| **Diffusion Pseudotime (DPT)** | 简单快速 | 基于图，无法整合先验 | GPLVM 可整合时间先验 |
| **Monocle3** | 流行方法 | 确定性，无不确定性 | GPLVM 是概率模型 |
| **Slingshot** | 处理分支 | 需要预定义聚类 | GPLVM 无需聚类 |
| **PAGA** | 处理复杂拓扑 | 基于图，离散 | GPLVM 连续平滑 |

### 4.2 GPLVM 的独特优势

1. **时间先验整合**：
   - 可以直接使用采样时间作为先验均值
   - 对于无时间数据，可从参考数据集迁移

2. **概率框架**：
   - 提供 pseudotime 估计的不确定性
   - 可用于下游分析的置信度加权

3. **多潜在因子解耦**：
   - 第一维：生物学时间
   - 第二维：技术效应（batch）或其他生物变异

4. **平滑性保证**：
   - 高斯过程核函数保证连续映射
   - 适合后续的高斯过程插值

---

## 五、参考文献 (References)

### 5.1 GPLVM 核心文献

1. **Lawrence, N.D. (2005)**. Probabilistic non-linear principal component analysis with Gaussian process latent variable models. *Journal of Machine Learning Research*, 6, 1783-1816.
   - GPLVM 原始论文

2. **Titsias, M.K. & Lawrence, N.D. (2010)**. Bayesian Gaussian process latent variable model. *Proceedings of the Thirteenth International Conference on Artificial Intelligence and Statistics*, 844-851.
   - 稀疏变分 GPLVM

### 5.2 单细胞应用文献

3. **Ahmed, S., Rattray, M. & Boukouvalas, A. (2018)**. GrandPrix: scaling up the Bayesian GPLVM for single-cell data. *Bioinformatics*, 35(1), 47-54.
   - 可扩展的贝叶斯 GPLVM 用于单细胞数据
   - G2G 论文引用 [79]

4. **Lalchand, V. et al. (2022)**. Modelling technical and biological effects in scRNA-seq data with scalable GPLVMs. *Proceedings of the 17th Machine Learning in Computational Biology*, 46-60.
   - 整合技术和生物效应
   - G2G 论文引用 [81]

5. **Verma, A. & Engelhardt, B.E. (2020)**. A robust nonlinear low-dimensional manifold for single cell RNA-seq data. *BMC Bioinformatics*.
   - 鲁棒的非线性流形学习
   - G2G 论文引用 [82]

6. **Lönnberg, T. et al. (2017)**. Single-cell RNA-seq and computational analysis using temporal mixture modelling resolves Th1/Tfh fate bifurcation in malaria. *Science Immunology*, 2(9), eaal2192.
   - 时间混合模型解析命运分歧
   - G2G 论文引用 [83]

7. **Kumasaka, N. et al. (2023)**. Mapping interindividual dynamics of innate immune response at single-cell resolution. *Nature Genetics*.
   - 单细胞分辨率的先天免疫反应动态
   - G2G 论文引用 [80]

### 5.3 软件实现

8. **Bingham, E. et al. (2019)**. Pyro: deep universal probabilistic programming. *Journal of Machine Learning Research*, 20, 1-6.
   - Pyro 概率编程框架
   - G2G 使用 Pyro v.1.8.0 实现 GPLVM

### 5.4 G2G 论文中的应用

9. **Sumanaweera, D. et al. (2024)**. Gene-level alignment of single-cell trajectories. *Nature Methods*.
   - DOI: https://doi.org/10.1038/s41592-024-02378-4
   - Methods 部分 "Pseudotime estimation using GPLVM"（Page 21）
   - Extended Data Fig. 6 展示 pseudotime 推断流程

---

## 六、实际应用建议 (Practical Recommendations)

### 6.1 何时有/无时间先验

**有时间先验（如时间序列实验）**：
```
- 直接使用采样时间作为高斯先验的均值
- 初始化 X[:, 0] = 归一化的采样时间
- 收敛更快，结果更可靠
```

**无时间先验（如终点测量）**：
```
- 使用 kNN 从参考数据集近似
- 或使用其他 pseudotime 方法（如 DPT）初始化
- 需要更多迭代次数，可能需手动检查质量
```

### 6.2 质量检查

```python
# 1. 检查 GPLVM 损失曲线是否收敛
plt.plot(loss_history)
plt.xlabel('Iteration')
plt.ylabel('Negative Log Marginal Likelihood')

# 2. 检查 pseudotime 与已知标记基因的相关性
from scipy.stats import spearmanr
corr, pval = spearmanr(pseudotime, known_marker_expression)

# 3. 可视化 pseudotime 在 UMAP 上的连续性
plt.scatter(umap[:, 0], umap[:, 1], c=pseudotime, cmap='viridis')
```

### 6.3 计算资源

| 数据规模 | 推荐配置 | 预计时间 |
|----------|----------|----------|
| <10,000 细胞 | CPU, 8GB RAM | 10-30 分钟 |
| 10,000-50,000 细胞 | GPU, 16GB RAM | 30-60 分钟 |
| >50,000 细胞 | GPU, 32GB RAM, 更多诱导点 | 1-2 小时 |

---

## 七、与后续插值步骤的接口

```python
# GPLVM 输出
gplvm_output = {
    'pseudotime': pseudotime,  # shape: (n_cells,)
    'pseudotime_variance': pt_var,  # 不确定性估计
    'latent_2': latent_2,  # 第二维潜在变量
}

# 直接用于插值模块
interpolation_input = {
    'pseudotime': gplvm_output['pseudotime'],
    'expression': expr_matrix,
}

# 调用插值函数
interp_result = trajectory_interpolate(
    pseudotime=interpolation_input['pseudotime'],
    expr=interpolation_input['expression'],
    method='gp',
    kernel='RBF',
    m_interp_points=14
)
```

