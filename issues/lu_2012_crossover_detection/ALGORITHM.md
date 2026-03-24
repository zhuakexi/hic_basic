# Lu et al. (2012) Crossover Detection 算法详解

## 目录

1. [算法概述](#算法概述)
2. [数学原理](#数学原理)
3. [实现细节](#实现细节)
4. [参数优化](#参数优化)
5. [结果解读](#结果解读)

---

## 算法概述

### 问题定义

**输入:**
- 单个精子细胞的全基因组测序数据
- 亲本二倍体基因组的 SNP 基因型
- 亲本单倍型定相信息 (phasing information)

**输出:**
- 每个精子细胞中 crossover 事件的基因组位置
- 每个 crossover 的置信度和分辨率
- 单倍型区块 (haplotype block) 的划分

### 两步法流程

```
┌─────────────────────────────────────────────────────────────────────┐
│                        Crossover Detection                           │
│                        (Lu et al. 2012 Science)                      │
└─────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────┐
│  Step 1: Rough Estimation (粗略估计)                                │
│  ─────────────────────────────────────────────────────────────────  │
│  • 沿基因组以 10 个 SNP 为窗口滑动                                       │
│  • 计算每个窗口中父本/母本 SNP 的比例                                   │
│  • 识别比例发生显著转换的位置                                          │
│  • 输出：crossover 候选位点列表                                         │
│  • 计算：crossover 频率 (用于 HMM 转移概率)                              │
└─────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────┐
│  Step 2: Accurate Determination using HMM (HMM 精确定位)            │
│  ─────────────────────────────────────────────────────────────────  │
│  • 3 状态隐马尔可夫模型：P (父本), M (母本), U (未确定)                   │
│  • 转移概率：由 Step 1 的 crossover 频率计算                               │
│  • 发射概率：由 SNP calling 错误率计算                                  │
│  • Viterbi 算法：寻找最可能的状态序列                                    │
│  • 输出：每个 SNP 的单倍型归属                                           │
└─────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────┐
│  Filtering (过滤)                                                    │
│  ─────────────────────────────────────────────────────────────────  │
│  • 识别状态转换点 (P↔M) 作为 crossover 位置                             │
│  • 过滤小于 300kb 的冲突片段 (可能是基因转换或污染)                       │
│  • 输出：最终的 crossover 位置列表                                       │
└─────────────────────────────────────────────────────────────────────┘
```

---

## 数学原理

### 1. Rough Estimation

**SNP 比例计算:**

对于包含 k 个 SNP 的窗口 w，计算父本单倍型 SNP 的比例：

$$P_{ratio}(w) = \frac{N_{paternal}(w)}{N_{total}(w)}$$

其中：
- $N_{paternal}(w)$ = 窗口 w 中匹配父本单倍型的 SNP 数量
- $N_{total}(w)$ = 窗口 w 中总的有效 SNP 数量

**转换点检测:**

$$\text{如果 } |P_{ratio}(w_i) - P_{ratio}(w_{i+1})| > \tau$$

则位置 $x = \frac{pos(w_i) + pos(w_{i+1})}{2}$ 为 crossover 候选位点

其中阈值 $\tau$ 通常设为 0.5 (表示从主要父本切换到主要母本)

**Crossover 频率估计:**

$$\lambda = \frac{N_{crossovers}}{L_{genome} / 1Mb}$$

其中：
- $N_{crossovers}$ = 检测到的 crossover 候选数量
- $L_{genome}$ = 基因组总长度 (bp)

### 2. Hidden Markov Model

**模型定义:**

隐马尔可夫模型由五元组 $\lambda = (S, V, A, B, \pi)$ 定义：

**状态空间 (S):**
$$S = \{s_0, s_1, s_2\} = \{P, M, U\}$$
- $s_0 = P$: Paternal haplotype (父本单倍型)
- $s_1 = M$: Maternal haplotype (母本单倍型)
- $s_2 = U$: Undetermined (未确定)

**观测空间 (V):**
$$V = \{v_0, v_1, v_2\} = \{allele_P, allele_M, unknown\}$$
- $v_0$: 观察到父本等位基因
- $v_1$: 观察到母本等位基因
- $v_2$: 未确定/杂合

**状态转移概率矩阵 (A):**

$$A = [a_{ij}] = \begin{bmatrix}
1 - p_c - p_u & p_c & p_u \\
p_c & 1 - p_c - p_u & p_u \\
p_u & p_u & 1 - 2p_u
\end{bmatrix}$$

其中：
- $p_c = \lambda / 2$ (crossover 概率，来自 Step 1)
- $p_u = 0.01$ (进入/退出未确定状态的概率)

**发射概率矩阵 (B):**

$$B = [b_j(k)] = \begin{bmatrix}
1 - \epsilon & \epsilon & \delta \\
\epsilon & 1 - \epsilon & \delta \\
\delta & \delta & 1 - 2\delta
\end{bmatrix}$$

其中：
- $\epsilon$ = SNP calling 错误率 (默认 0.05)
- $\delta$ = 小常数 (默认 0.01)

**初始概率分布 (π):**

$$\pi = [1/3, 1/3, 1/3]$$ (均匀分布)

### 3. Viterbi 算法

**目标:**

寻找最可能的状态序列 $S^*$：

$$S^* = \arg\max_{S} P(S | O, \lambda)$$

其中 $O = (o_1, o_2, ..., o_T)$ 是观测序列。

**动态规划:**

定义 $\delta_t(i)$ 为时刻 t 处于状态 i 的所有路径中的最大概率：

$$\delta_t(i) = \max_{s_1,...,s_{t-1}} P(s_1,...,s_{t-1}, s_t=i | o_1,...,o_t, \lambda)$$

**递推关系:**

$$\delta_t(j) = \max_{i} [\delta_{t-1}(i) \cdot a_{ij}] \cdot b_j(o_t)$$

**回溯指针:**

$$\psi_t(j) = \arg\max_{i} [\delta_{t-1}(i) \cdot a_{ij}]$$

**算法步骤:**

```
1. 初始化 (t=1):
   δ₁(i) = πᵢ · bᵢ(o₁)
   ψ₁(i) = 0

2. 递推 (t=2,3,...,T):
   δₜ(j) = maxᵢ[δₜ₋₁(i) · aᵢⱼ] · bⱼ(oₜ)
   ψₜ(j) = argmaxᵢ[δₜ₋₁(i) · aᵢⱼ]

3. 终止:
   P* = maxᵢ[δₜ(i)]
   sₜ* = argmaxᵢ[δₜ(i)]

4. 回溯 (t=T-1,...,1):
   sₜ* = ψₜ₊₁(sₜ₊₁*)
```

### 4. Crossover 位置确定

从 Viterbi 路径 $S^* = (s_1^*, s_2^*, ..., s_T^*)$ 中识别状态转换：

$$\text{Crossover 位置} = \{t | s_t^* \neq s_{t+1}^*, s_t^*, s_{t+1}^* \in \{P, M\}\}$$

具体位置计算：

$$pos_{crossover} = \frac{pos_t + pos_{t+1}}{2}$$

---

## 实现细节

### 数据结构

```python
@dataclass
class HaplotypeBlock:
    """单倍型区块"""
    chromosome: str      # 染色体
    start_pos: int       # 起始位置
    end_pos: int         # 终止位置
    haplotype: str       # 'P', 'M', 或 'U'
    confidence: float    # 置信度 (0-1)
    num_snps: int        # 包含的 SNP 数量

@dataclass
class Crossover:
    """Crossover 事件"""
    sperm_id: str        # 精子样本 ID
    chromosome: str      # 染色体
    position: int        # Crossover 位置
    start_pos: int       # 区块起始
    end_pos: int         # 区块终止
    from_haplotype: str  # 转换前单倍型
    to_haplotype: str    # 转换后单倍型
    confidence: float    # 置信度
    resolution_kb: float # 分辨率 (kb)
```

### 关键函数

#### 1. 观测编码

```python
def _encode_observations(sperm_gt, h1, h2):
    """
    将基因型编码为 HMM 观测
    
    Args:
        sperm_gt: 精子基因型
        h1: 父本单倍型等位基因
        h2: 母本单倍型等位基因
    
    Returns:
        0: 匹配父本
        1: 匹配母本
        2: 未确定
    """
    if sperm_gt in ['N', '.', '-']:
        return 2  # 未确定
    elif sperm_gt == h1:
        return 0  # 父本
    elif sperm_gt == h2:
        return 1  # 母本
    else:
        return 2  # 不匹配/未确定
```

#### 2. 小片段过滤

```python
def _filter_small_fragments(positions, states, min_size=300000):
    """
    过滤小于最小尺寸的冲突片段
    
    原因：
    - 基因转换事件通常 < 10kb
    - DNA 污染可能导致小的冲突片段
    - 真实的 crossover 间隔通常 > 1Mb
    """
    filtered = states.copy()
    i = 0
    while i < len(positions):
        # 找到当前状态的片段边界
        current_state = states[i]
        start_idx = i
        start_pos = positions[i]
        
        while i < len(positions) and states[i] == current_state:
            i += 1
        
        end_idx = i - 1
        end_pos = positions[end_idx]
        fragment_size = end_pos - start_pos
        
        # 过滤小的未确定片段
        if fragment_size < min_size and current_state == STATE_UNDETERMINED:
            if start_idx > 0:
                filtered[start_idx:end_idx+1] = states[start_idx - 1]
            elif end_idx < len(positions) - 1:
                filtered[start_idx:end_idx+1] = states[end_idx + 1]
    
    return filtered
```

### 数值稳定性

为防止概率下溢，实现中使用归一化：

```python
# 在 Viterbi 递推中
for t in range(1, n_obs):
    for j in range(n_states):
        probs = delta[t-1, :] * transition_matrix[:, j]
        psi[t, j] = np.argmax(probs)
        delta[t, j] = probs[psi[t, j]] * emission_matrix[j, observations[t]]
    
    # 归一化防止下溢
    delta[t, :] /= np.sum(delta[t, :]) + 1e-300
```

---

## 参数优化

### 默认参数

| 参数 | 默认值 | 说明 | 优化建议 |
|------|--------|------|----------|
| `window_size` | 10 | SNP 窗口大小 | 高覆盖度数据可增大到 15-20 |
| `min_fragment_size` | 300,000 bp | 最小片段大小 | 不要小于 100kb |
| `error_rate` | 0.05 | SNP calling 错误率 | 根据测序质量调整 |
| `crossover_rate` | 0.01 | 先验 crossover 率 | 由 Step 1 自动估计 |

### 参数敏感性分析

**Window Size:**
- 太小 → 噪声敏感，假阳性 crossover
- 太大 → 分辨率降低，可能漏检邻近 crossover

**Error Rate:**
- 太低 → 过度信任观测，对测序错误敏感
- 太高 → 过度平滑，可能漏检真实 crossover

**Min Fragment Size:**
- 太小 → 保留过多噪声/污染片段
- 太大 → 可能过滤掉真实的基因转换事件

---

## 结果解读

### 输出文件格式

#### 1. Haplotype Blocks (`*_haplotype_blocks.tsv`)

```
chromosome  start_pos   end_pos     haplotype  confidence  num_snps
chr1        10023456    25678901    P          0.95        1523
chr1        25678901    58234567    M          0.93        3245
chr1        58234567    89012345    P          0.91        2876
```

**解读:**
- 连续的 P 或 M 区块表示来源于同一亲本单倍型的区域
- 区块边界即为 crossover 位置
- confidence 接近 1 表示该区块内 SNP 一致性好

#### 2. Crossovers (`*_crossovers.tsv`)

```
sperm_id  chromosome  position    from_haplotype  to_haplotype  confidence  resolution_kb
S01       chr1        25678901    P               M             0.94        30.5
S01       chr1        58234567    M               P             0.92        28.3
```

**解读:**
- position: crossover 的基因组位置
- from/to_haplotype: 单倍型转换方向
- resolution_kb: 定位分辨率 (越低越好)

#### 3. Statistics (`*_statistics.tsv`)

```
sperm_id  total_crossovers  total_blocks  chromosomes_analyzed  mean_resolution_kb  mean_confidence
S01       26                27            23                    30.5                0.92
```

**预期值 (根据 Lu et al. 2012):**
- 常染色体 crossover 总数：~26 per sperm
- 平均分辨率：~30 kb
- 平均置信度：> 0.90

### 质量控制指标

**通过 QC 的标准:**
1. 总 crossover 数在 20-35 之间
2. 平均置信度 > 0.85
3. 平均分辨率 < 50 kb
4. 每条染色体至少有 1 个 crossover

**可能的数据问题:**
- crossover 数 < 15: 可能覆盖度不足或 DNA 污染
- crossover 数 > 40: 可能样本混合或算法参数不当
- 置信度 < 0.80: SNP calling 质量可能有问题

---

## 参考文献

1. Lu S, Zong C, Fan W, et al. Probing Meiotic Recombination and Aneuploidy of Single Sperm Cells by Whole-Genome Sequencing. Science. 2012;338(6114):1627-1630.

2. Rabiner LR. A tutorial on hidden Markov models and selected applications in speech recognition. Proceedings of the IEEE. 1989;77(2):257-286.

3. Viterbi A. Error bounds for convolutional codes and an asymptotically optimum decoding algorithm. IEEE Transactions on Information Theory. 1967;13(2):260-269.
