# 单一轨迹高斯过程插值方法

## 背景 (Background)

在单细胞轨迹分析中，原始数据通常是离散的细胞观测点，每个细胞具有连续的 pseudotime 值。为了获得平滑的基因表达动态轨迹，需要对离散观测进行插值和平滑处理。

G2G (Genes2Genes, Nature Methods 2024) 方法提出了一套基于高斯过程 (Gaussian Process, GP) 的插值框架，originally designed for aligning two trajectories。本申请将其改造用于**单一轨迹的平滑插值**，可应用于：
- Hi-C 数据的时间序列分析
- 单细胞基因表达的 pseudotime 轨迹平滑
- 任何具有连续时间/状态维度的组学数据插值

---

## 目标 (Objectives)

实现一个通用的单一轨迹高斯过程插值模块，具备以下功能：

1. **输入**：
   - 细胞/样本的 pseudotime 值（范围 [0, 1]）
   - 基因表达矩阵（细胞 × 基因）
   - 可选：协变量（如 batch 信息）

2. **输出**：
   - 在等间距时间点上的平滑表达量估计
   - 每个时间点的表达量分布（均值 + 方差）
   - 可选：插值质量评估指标

3. **核心特性**：
   - 基于高斯过程的概率插值（提供不确定性估计）
   - 自适应确定最优插值分辨率
   - 支持大规模单细胞数据（数千细胞，数千基因）

---

## 算法详细步骤 (Algorithm Details)

### 步骤 1：数据预处理

#### 1.1 表达量标准化
```
输入：原始计数矩阵 counts[cell, gene]
输出：标准化表达矩阵 expr[cell, gene]

步骤：
1. 每细胞总计数标准化到 10,000：
   normalized[cell, gene] = counts[cell, gene] / sum(counts[cell, :]) * 10000

2. log1p 转换：
   expr[cell, gene] = log(1 + normalized[cell, gene])
```

#### 1.2 Pseudotime 归一化
```
输入：原始 pseudotime 值 pseudotime_raw[cell]
输出：归一化 pseudotime ∈ [0, 1]

步骤：
pseudotime[cell] = (pseudotime_raw[cell] - min(pseudotime_raw)) / 
                   (max(pseudotime_raw) - min(pseudotime_raw))
```

#### 1.3 高变基因选择（可选）
```
输入：expr[cell, gene]
输出：高变基因子集 HVG

步骤：
1. 计算每个基因的均值和离散度
2. 选择 top N (如 2000-4000) 高变基因
```

---

### 步骤 2：确定最优插值分辨率

#### 2.1 使用 Optimal Binning 确定插值点数 m

**原理**：根据 pseudotime 分布的内在结构，自动确定合适的插值分辨率，避免过密（冗余/噪声）或过疏（信息丢失）。

```
输入：pseudotime[cell]
输出：最优插值点数 m

步骤：
1. 使用 ContinuousOptimalBinning (optBinning 包 v.0.18.0+)
   - 输入：pseudotime 数组
   - 参数：
     * n_bins: 初始分箱数（如 20）
     * min_bin_size: 最小分箱大小（如总细胞数的 5%）
     * monotonic_trend: 可选约束

2. 拟合分箱模型：
   binning = ContinuousOptimalBinning(...)
   binning.fit(pseudotime)

3. 获取最优分箱数：
   m = binning.n_bins

4. 对多个基因/条件分别计算时，取最小值保证一致性：
   m_final = min(m_ref, m_query)
```

**默认值**：如果 optimal binning 失败，使用 m = 14 作为默认值（G2G 论文经验值）。

---

### 步骤 3：定义插值时间点

#### 3.1 等间距插值网格
```
输入：插值点数 m
输出：插值时间点 grid

步骤：
grid = np.linspace(0, 1, m)

示例：m=14 时，grid = [0, 0.077, 0.154, ..., 0.923, 1.0]
```

---

### 步骤 4：高斯过程插值

#### 4.1 高斯过程基础

**模型**：
```
f(t) ~ GP(μ(t), K(t, t'))

其中：
- μ(t): 均值函数（常数为常向量或时间的函数）
- K(t, t'): 协方差核函数，控制平滑度
```

#### 4.2 核函数选择

**推荐：RBF (Radial Basis Function) 核**
```
K(t, t') = σ² * exp(-||t - t'||² / (2 * l²))

参数：
- σ²: 信号方差（signal variance）
- l: 长度尺度（length scale），控制平滑度
  * 较大的 l → 更平滑的函数
  * 较小的 l → 更灵活的函数
```

**其他可选核函数**：
- **Matérn 核**：适用于较不平滑的轨迹
- **Change Point 核**：适用于有分支/转折的轨迹
- **Periodic 核**：适用于周期性轨迹

#### 4.3 对每个基因进行 GP 插值

```
输入：
- pseudotime[cell]: 细胞的 pseudotime
- expr[cell]: 单个基因的表达量向量
- grid: 插值时间点

输出：
- mean_interp: 插值点上的均值估计
- var_interp: 插值点上的方差估计

步骤：
1. 构建训练数据：
   X_train = pseudotime.reshape(-1, 1)  # shape: (n_cells, 1)
   y_train = expr[cell]                  # shape: (n_cells,)

2. 定义 GP 模型（以 scikit-learn 为例）：
   from sklearn.gaussian_process import GaussianProcessRegressor
   from sklearn.gaussian_process.kernels import RBF, ConstantKernel
   
   kernel = ConstantKernel(σ²=1.0) * RBF(length_scale=l=0.1)
   gp = GaussianProcessRegressor(
       kernel=kernel,
       alpha=noise_variance,  # 观测噪声
       normalize_y=True,
       n_restarts_optimizer=5
   )

3. 拟合 GP 模型：
   gp.fit(X_train, y_train)

4. 在插值点上预测：
   X_grid = grid.reshape(-1, 1)  # shape: (m, 1)
   mean_interp, var_interp = gp.predict(X_grid, return_var=True)

5. 返回结果：
   return mean_interp, var_interp
```

#### 4.4 批量处理所有基因

```
输入：expr[cell, gene] 表达矩阵
输出：interp_mean[gene, time], interp_var[gene, time]

步骤：
对每个基因 g in genes:
    mean_g, var_g = GP_interpolate(pseudotime, expr[:, g], grid)
    存储结果

优化策略：
- 使用多进程并行处理（每基因独立）
- 或使用矩阵运算批量拟合（如果内存允许）
```

---

### 步骤 5：Window-based 局部加权平滑（可选增强）

**目的**：在 GP 插值前，先对原始数据进行局部平滑，减少噪声影响。

#### 5.1 定义局部邻域

```
输入：
- pseudotime[cell]
- window_size (默认 0.1，即 10% 的 pseudotime 范围)
- 插值点 t_grid[time]

输出：每个插值时间点的邻域细胞索引
```

#### 5.2 高斯加权

```
对每个插值时间点 t_i：
1. 计算距离：
   dist[cell] = |pseudotime[cell] - t_i|

2. 高斯权重：
   weight[cell] = exp(-dist[cell]² / (2 * window_size²))

3. 加权均值和方差：
   weighted_mean[g] = Σ(weight[cell] * expr[cell, g]) / Σ(weight[cell])
   weighted_var[g] = Σ(weight[cell] * (expr[cell, g] - weighted_mean[g])²) / Σ(weight[cell])
```

#### 5.3 与 GP 插值的结合

**方案 A：预处理**
- 先用 window-based 方法计算每个时间点的加权统计量
- 将这些统计量作为 GP 的观测输入

**方案 B：后处理**
- 先进行 GP 插值
- 再用 window-based 方法对 GP 输出进行平滑

**推荐**：方案 A，保留 GP 的不确定性估计能力。

---

### 步骤 6：输出格式

#### 6.1 主要输出

```python
{
    # 插值时间点 (m,)
    "time_grid": np.array([0.0, 0.077, ..., 1.0]),
    
    # 基因 × 时间点的均值估计 (n_genes, m)
    "mean": np.array(...),
    
    # 基因 × 时间点的方差估计 (n_genes, m)
    "variance": np.array(...),
    
    # 元数据
    "metadata": {
        "n_cells": ...,
        "n_genes": ...,
        "m_interp_points": ...,
        "kernel_type": "RBF",
        "window_size": 0.1,
    }
}
```

#### 6.2 可选输出

```python
{
    # 每个基因的 GP 超参数
    "gp_hyperparameters": {
        "gene1": {"length_scale": ..., "signal_var": ..., "noise_var": ...},
        ...
    },
    
    # 插值质量评估
    "quality_metrics": {
        "gene1": {"log_marginal_likelihood": ..., "rsquared": ...},
        ...
    },
    
    # 原始 pseudotime 分布信息
    "pseudotime_distribution": {
        "min": ...,
        "max": ...,
        "median": ...,
        "n_cells_per_bin": [...]
    }
}
```

---

## 参数默认值总结

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `n_hvg` | 2000-4000 | 高变基因数量 |
| `min_bin_size` | 0.05 | OptimalBinning 最小分箱比例 |
| `window_size` | 0.1 | 局部加权平滑的窗口大小 |
| `kernel` | RBF | 高斯过程核函数类型 |
| `length_scale` | 0.1-0.3 | RBF 核长度尺度（可自动优化） |
| `alpha` | 1e-6 | 观测噪声方差 |
| `n_restarts` | 5 | GP 超参数优化重启次数 |

---

## 依赖包建议

```python
# 核心依赖
numpy >= 1.20.0
scipy >= 1.7.0
scikit-learn >= 1.0.0  # GaussianProcessRegressor

# Optimal Binning
optbinning >= 0.18.0

# 可选：更强大的 GP 实现
GPyTorch >= 1.5.0  # 适用于大规模数据
# 或
GPy >= 1.10.0

# 数据处理
pandas >= 1.3.0
anndata >= 0.8.0  # 如果使用 AnnData 格式
```

---

## 预期应用场景

1. **Hi-C 时间序列插值**：
   - 输入：不同时间点的 Hi-C 接触矩阵
   - 输出：连续时间上的染色质结构动态

2. **单细胞 pseudotime 轨迹**：
   - 输入：细胞 pseudotime + 基因表达
   - 输出：平滑的基因表达动态曲线

3. **跨条件轨迹比较的预处理**：
   - 先对每个条件单独插值
   - 再进行差异分析或对齐

---

## 参考文献

1. Sumanaweera, D. et al. Gene-level alignment of single-cell trajectories. *Nature Methods* (2024). https://doi.org/10.1038/s41592-024-02378-4

2. Rasmussen, C.E. & Williams, C.K.I. Gaussian Processes for Machine Learning. MIT Press (2006).

3. Navas-Palencia, G. Optimal binning: mathematical programming formulation. arXiv:2001.08025 (2022).
