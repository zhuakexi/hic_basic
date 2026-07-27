# Issue: 测试套件长时阻塞与 GenomeIdeograph 导入断裂

状态：Resolved

解决日期：2026-07-27

实现说明：ATAC correctness 测试现使用临时 synthetic BAM/TSS 数据，真实 BAM
测试通过 `integration`/`slow` marker 显式运行；coverage 复用单个 indexed BAM handle。
`DI` 与 `impute.simpute` 已迁移到 `hic_basic.genome`，旧 `binnify` 导入保留至
`hic_basic 0.2` 并发出弃用警告。未引入可选的 `pytest-xdist`，因为仓库其余共享输出和
机器数据测试尚未完成并行安全审计。

日期：2026-07-27

范围：默认测试分层、ATAC coverage 测试、`GenomeIdeograph` 迁移兼容性、可选测试并行

## 1. 摘要

当前 `tests/` 不能作为一次可靠、可移植的默认回归测试运行，存在两个相互独立的问题：

1. `tests/test_atac.py::TestATACFunctions::test_calculate_tss_enrichment` 实际运行
   机器特定的真实数据集成任务。在当前机器上，前 159 个测试完成后，该测试单核持续
   运行超过 13 分钟。它读取 6.1 GB BAM，并对 100 个 TSS 区间逐个启动
   `pysam.depth`。
2. `hic_basic.DI` 和 `hic_basic.impute.simpute` 仍从 `hic_basic.binnify` 导入
   `GenomeIdeograph`，但该类已于提交 `8f9fee7` 移到 `hic_basic.genome`。因此完整测试
   在 collection 阶段就会报两个 `ImportError`。

这两个问题均与 `hic_basic.visualize` 和 Playwright 渲染后端无关，但会遮蔽其他模块的
真实回归结果，并造成“完整测试一直很慢或无法运行”的反复反馈。

## 2. 问题 A：ATAC 测试不是可移植的单元测试

### 2.1 复现

```bash
./run_in_env.sh pytest -vv \
  tests/test_atac.py::TestATACFunctions::test_calculate_tss_enrichment
```

完整套件在临时排除两个导入错误后也会进入同一路径：

```bash
./run_in_env.sh pytest -q tests \
  --ignore=tests/test_DI.py \
  --ignore=tests/test_simpute.py
```

观测基线：停止前为 `157 passed, 2 skipped, 0 failed`，随后在该测试内运行约
13 分钟。进程保持约一个 CPU 核心满载，不是 pytest 死锁。

### 2.2 当前实现问题

`tests/test_atac.py:17-23` 在 `setUp` 中写死以下机器路径：

- 仓库绝对路径下的染色体长度文件；
- `/shared/ychi/repo/hsv_15/...` 下的 6.1 GB BAM；
- `/share/Data/...` 下的 GTF；
- `tests/output/atac` 固定输出目录。

`test_calculate_tss_enrichment` 的执行链为：

```text
5.2 MB GRCh38 TSS BED
  -> get_tss_regions(window_size=2000)
  -> random.sample(..., 100)
  -> 对每个区间单独调用一次 pysam.depth
  -> 100 次进程启动 + 共享存储随机读取
```

`hic_basic/atac.py:223-227` 在 Python 循环中对每个区间分别调用
`pysam.depth("-b", "-", bam_path, input=...)`。这使 `sampling=100` 仍需要 100 次
外部命令启动，并反复访问同一个大型 BAM/BAM index。

该测试完成昂贵步骤后还会执行 `normalize_coverage(..., read_length, ...)`，但
`read_length` 在 `tests/test_atac.py` 中没有定义。测试也没有任何 assertion。因此当前
状态即使长时计算完成，也不能形成有效的通过结果。

### 2.3 建议修复

1. 用 `tmp_path` 生成最小、已索引的 synthetic BAM 和小型 TSS BED，建立确定性的
   `calculate_coverage` 单元测试；固定随机种子并断言 shape、距离和 coverage 数值。
2. 将真实 6.1 GB BAM 测试迁移为显式 `integration`/`slow` 测试，默认测试不运行，且
   缺少数据时明确 skip，不允许硬编码用户目录。
3. 修复 `read_length` 未定义，并为归一化结果增加数值 assertion。
4. 优化 `calculate_coverage`：持久打开一个 indexed BAM，或一次性将全部采样 BED
   交给一次 `pysam.depth`，再按区间补零和重组。不要为每个区间启动一个命令。
5. 真实数据性能测试应记录数据大小、存储位置、区间数和耗时，不纳入普通 correctness
   suite。

### 2.4 为什么 pytest 多进程不是直接修复

`pytest-xdist` 可以并行不同测试，但不能把一个测试函数拆到多个 worker。当前正式环境
也没有注册 `pytest-xdist`。直接并发 100 个 BAM depth 任务可能使共享存储竞争更严重。

完成测试隔离和固定输出目录审计后，可以将 `pytest-xdist` 作为 `scope: test` 的正式
依赖注册并锁定，使用有界并发：

```bash
./run_in_env.sh pytest -n 4 --dist loadscope -q tests
```

不建议在共享节点默认使用 `-n auto`。浏览器、BLAS/OpenMP 和共享存储测试还需要控制
每个 worker 的内部线程数，避免嵌套并发。

## 3. 问题 B：GenomeIdeograph 导入迁移不完整

### 3.1 复现

```bash
./run_in_env.sh pytest --collect-only -q \
  tests/test_DI.py tests/test_simpute.py
```

当前错误：

```text
ImportError: cannot import name 'GenomeIdeograph' from 'hic_basic.binnify'
```

受影响的生产模块：

- `hic_basic/DI.py:8`
- `hic_basic/impute/simpute.py:11`

两者在模块 import 阶段失败，因此相关函数是否被调用并不重要；任何导入这些模块的
用户、测试或 API 文档工具都会立即失败。

### 3.2 根因

提交 `8f9fee7` 的目的明确为 “move genome ideograph to genome rather than binnify”。
`GenomeIdeograph` 的实现当前位于 `hic_basic/genome.py:212`，仓库中其他活跃模块已经
使用以下新路径：

```python
from hic_basic.genome import GenomeIdeograph
```

但 `DI.py` 和 `impute/simpute.py` 遗漏了迁移。`binnify.py` 中只剩注释里的类型名，
没有兼容 re-export。

### 3.3 建议修复

1. 将两个内部导入改为 `hic_basic.genome.GenomeIdeograph`。
2. 搜索并清理所有活跃代码中的旧导入，注释和历史类型说明可单独处理。
3. 明确是否需要保留外部兼容路径：
   - 如果 `hic_basic.binnify.GenomeIdeograph` 曾是公共 API，在 `binnify.py` 增加带弃用
     说明的 re-export，并给出移除版本；
   - 如果不是公共 API，只修复内部调用，不新增永久 shim。
4. 增加模块导入 smoke test，至少覆盖：

```python
import hic_basic.DI
import hic_basic.impute.simpute
from hic_basic.genome import GenomeIdeograph
```

5. 修复后重新生成 API reference/index，并验证 `hic_basic api ...` 查询不受影响。

## 4. 验收标准

- `./run_in_env.sh pytest --collect-only -q tests` 不再出现上述两个 `ImportError`。
- 默认测试不读取 `/shared/ychi/repo/hsv_15`、`/share/Data` 或其他机器特定绝对路径。
- ATAC coverage 单元测试使用最小临时数据，在普通 CPU 环境数秒内完成，并包含有效
  assertion。
- 真实 BAM 测试有 `integration`/`slow` 标记、明确资源前提和单独执行命令。
- `DI` 与 `impute.simpute` 均可独立导入，并至少有一个最小功能测试使用新的
  `GenomeIdeograph` 路径。
- 在审计共享输出和线程嵌套后，若引入 `pytest-xdist`，`-n 4 --dist loadscope` 与串行
  suite 结果一致。

## 5. 非目标

- 本 issue 不改变 `GenomeIdeograph` 的 binning 语义。
- 本 issue 不要求把真实 BAM 基准完全改为 CPU 计算基准。
- 本 issue 不用提高 pytest worker 数来掩盖单测试算法和 I/O 设计问题。
