# hic_basic API command plan

## 背景和目标

`hic_basic` 已经有一套 API 文档和速查系统：

- `scripts/generate_api_reference.py` 扫描 `hic_basic/`，生成 API 文档。
- `docs/API_REFERENCE.md` 是人读的完整参考。
- `docs/api_index.json` 是机器读的索引。
- `docs/search_api.py` 支持按函数名、模块名、docstring keyword 查询。
- `AGENTS.md`、`README.md`、`QUICK_REFERENCE.md` 已经把这套系统写进 agent 工作流。

当前问题不是“没有 API 文档”，而是这套速查工具要求 agent 进入 `hic_basic` repo 并访问 `docs/`。在下游项目里，coding agent 更需要一个安装后可用的命令：

```bash
hic_basic api ...
```

目标是让 agent 在任何依赖 `hic_basic` 的项目中，都能通过 shell toolcall 快速完成：

1. 已知函数名时，确认 import、signature、docstring、源码位置。
2. 未知函数名时，通过模块/关键词搜索定位候选 API。
3. 需要机器解析时，用 `--json` 获取稳定结构。

这个工具应当是现有文档系统的安装后 CLI 入口，而不是另起一套文档系统。

## 当前现实约束

实现前必须承认这些现状：

- 当前 console script 只有 `hic_basic`，由 `setup.py` 注册为 `hic_basic=hic_basic.__main__:main`。
- 当前没有 `hic-basic` 或 `hb` 命令。计划中统一使用 `hic_basic api ...` 作为 canonical command。
- 当前 CLI 框架在 `hic_basic/cli/__init__.py`，通过模块里的 `name`、`description`、`add_arguments()`、`run()` 注册子命令。
- 当前 `python -m hic_basic --help` 只有 `version`、`render`、`download`、`ana-describe`。
- 当前只有 `setup.py`，没有 `pyproject.toml`、`setup.cfg`、`MANIFEST.in`。
- 如果 CLI 需要读取包内 JSON/registry，必须显式配置 package data。
- 当前 `docs/api_index.json` 的基本字段是 `name`、`signature`、`docstring`、`return_type`、`lineno`。
- 当前生成脚本会把 JSON index 里的 docstring 截断到 200 字符，这不足以支撑完整 `api show`。
- 当前仓库真实模块包括 `coolstuff`、`hicio`、`compartment`、`compare_pixels`、`embedding`、`plot.hic`、`shuffle` 等。
- 早期讨论中的 `hic_basic.compare.scc`、`hic_basic.matrix.band_zscore`、`hic_basic.eigs.eigs_cis` 是产品形态示例，不是当前真实 API，不能作为测试或 registry 默认内容。

## 设计原则

1. 复用现有系统，不复制逻辑。
2. 自动索引是主数据源，人工 registry 只是 overlay。
3. 文本输出为 agent 速读优化，短、稳定、可复制。
4. JSON 输出为工具解析优化，字段稳定，避免依赖文本格式。
5. 已知函数名必须能直达；未知任务必须能逐步收敛。
6. 不引入 embedding/RAG；先用确定性的 name/module/docstring/metadata 搜索。
7. 不要求 Phase 1 有完美 task semantics；先解决下游项目里“能查”的问题。
8. `docs/search_api.py` 短期保留兼容，不能破坏现有文档和 AGENTS.md 工作流。

## 推荐命令面

最小可用命令：

```bash
hic_basic api list [MODULE]
hic_basic api search QUERY...
hic_basic api show SYMBOL
hic_basic api where SYMBOL
```

机器可读输出：

```bash
hic_basic api list --json
hic_basic api list coolstuff --json
hic_basic api search "convert pairs" --json
hic_basic api show cool2mat --json
hic_basic api where cool2mat --json
```

后续增强命令：

```bash
hic_basic api map
hic_basic api examples [SYMBOL|SECTION]
hic_basic api source SYMBOL
hic_basic api schema [SYMBOL|--all]
```

`schema` 不列为 Phase 1 核心，因为 `--json` 已经能满足 agent 的机器解析需求。

## 数据源设计

### 自动索引

自动索引继续由 `scripts/generate_api_reference.py` 生成，作为全包 API 覆盖的主数据源。

建议生成两份同源文件：

```text
docs/api_index.json
hic_basic/data/api_index.json
```

`docs/api_index.json` 保持现有文档系统兼容。

`hic_basic/data/api_index.json` 供安装后的 `hic_basic api ...` 使用。

两份文件必须由同一个生成函数输出，避免漂移。

### JSON schema 建议

当前 schema 可以扩展为：

```json
{
  "module": "coolstuff",
  "name": "cool2mat",
  "symbol": "hic_basic.coolstuff.cool2mat",
  "import": "from hic_basic.coolstuff import cool2mat",
  "signature": "cool2mat(coolp, region, balance=False, ...)",
  "summary": "first non-empty docstring line",
  "docstring": "full docstring",
  "return_type": null,
  "source": {
    "module_path": "hic_basic/coolstuff.py",
    "lineno": 123
  }
}
```

注意：

- `docstring` 不应再截断到 200 字符。
- 如果担心 JSON 太大，可同时提供 `summary` 和完整 `docstring`；CLI 默认短输出用 `summary`。
- `symbol` 和 `import` 应在生成时写入，减少 CLI 运行时拼接错误。
- `source.module_path` 使用相对路径；`where` 命令运行时再解析安装环境里的绝对路径。

### 人工 registry overlay

不要把 registry 作为主索引。推荐后续增加轻量 overlay：

```text
hic_basic/data/api_registry.json
```

用途：

- 定义任务导向 section，例如 `cool`、`io`、`plot`、`compartment`、`tracks`、`phasing`。
- 给高频 API 增加 `tags`、`inputs`、`outputs`、`examples`、`pitfalls`、`related`。
- 改善 `api map`、`api list <section>`、`api search` 的任务语义。

不建议 Phase 1 使用 YAML，除非确认愿意新增 `PyYAML` 依赖。JSON 或 Python dict 更符合当前 `setup.py` 的低依赖状态。

overlay 示例：

```json
{
  "sections": {
    "cool": {
      "summary": "Read, write, convert, and inspect .cool/.mcool/scool data.",
      "tags": ["cool", "mcool", "scool", "pairs", "pixels"]
    },
    "plot": {
      "summary": "Plot Hi-C heatmaps, tracks, and multi-panel visualizations.",
      "tags": ["plot", "heatmap", "track", "figure"]
    }
  },
  "symbols": {
    "hic_basic.coolstuff.cool2mat": {
      "section": "cool",
      "tags": ["cool", "matrix", "region", "extract"],
      "inputs": ["cool path", "region"],
      "outputs": ["matrix"],
      "related": ["hic_basic.coolstuff.cool2pixels"]
    }
  }
}
```

CLI 读取时以自动索引为 base，再用 registry 字段覆盖或补充。

## 命令规格

### `hic_basic api list [MODULE]`

用途：

- 无参数时列出当前索引里的模块和函数数量。
- 有参数时列出匹配模块里的函数短条目。

示例：

```bash
hic_basic api list
hic_basic api list coolstuff
hic_basic api list plot.hic
hic_basic api list coolstuff --json
```

文本输出建议：

```text
Module: coolstuff
Functions: 32

cool2mat
  Signature: cool2mat(...)
  Summary: ...
  Show: hic_basic api show hic_basic.coolstuff.cool2mat

pairs2cool
  Signature: pairs2cool(...)
  Summary: ...
  Show: hic_basic api show hic_basic.coolstuff.pairs2cool
```

无 registry 时，`list <section>` 暂不实现 section 语义，只做 module/partial module match。

有 registry 后，`list cool` 可以解释为 section；如果 section 和 module 同名冲突，优先 exact module，再提示 section 候选。

### `hic_basic api search QUERY...`

用途：

- 从函数名、模块名、summary、docstring 中找候选。
- Phase 2 后加入 registry 的 section/tags/inputs/outputs/related。

示例：

```bash
hic_basic api search cool2mat
hic_basic api search convert pairs
hic_basic api search "plot cool tracks"
hic_basic api search "compartment"
```

Phase 1 搜索权重：

```text
exact function name: 100
function name contains query token: 50
module contains query token: 30
summary contains query token: 20
docstring contains query token: 10
```

Phase 2 搜索权重增加：

```text
registry tags: 40
registry section: 30
registry inputs/outputs: 20
registry related: 10
```

文本输出建议：

```text
Search: convert pairs

1. hic_basic.coolstuff.pairs2cool
   Signature: pairs2cool(...)
   Summary: ...
   Show: hic_basic api show hic_basic.coolstuff.pairs2cool

2. hic_basic.coolstuff.pairs2scool
   Signature: pairs2scool(...)
   Summary: ...
   Show: hic_basic api show hic_basic.coolstuff.pairs2scool
```

默认限制结果数量，例如 `--limit 10`。

JSON 输出字段建议：

```json
{
  "query": "convert pairs",
  "matches": [
    {
      "rank": 1,
      "score": 80,
      "symbol": "hic_basic.coolstuff.pairs2cool",
      "module": "coolstuff",
      "name": "pairs2cool",
      "signature": "pairs2cool(...)",
      "summary": "...",
      "show_command": "hic_basic api show hic_basic.coolstuff.pairs2cool"
    }
  ]
}
```

### `hic_basic api show SYMBOL`

用途：

- 已知函数名或 symbol 时，显示写代码所需信息。
- 支持短名和全名。
- 短名唯一匹配时直接显示；多重匹配时返回候选，不猜。

示例：

```bash
hic_basic api show cool2mat
hic_basic api show hic_basic.coolstuff.cool2mat
hic_basic api show cool2mat --brief
hic_basic api show cool2mat --signature
hic_basic api show cool2mat --json
```

文本输出建议：

```text
hic_basic.coolstuff.cool2mat

Import:
  from hic_basic.coolstuff import cool2mat

Signature:
  cool2mat(...)

Summary:
  ...

Docstring:
  ...

Source:
  hic_basic/coolstuff.py:123
```

`--brief` 只输出 symbol、import、signature、summary。

`--signature` 只输出 import 和 signature，便于 agent 快速写脚本。

如果 symbol 不唯一：

```text
Ambiguous symbol: plot

Matches:
  1. hic_basic.plot.hic.plot_cool
     Show: hic_basic api show hic_basic.plot.hic.plot_cool
  2. hic_basic.plot.track.plot_track
     Show: hic_basic api show hic_basic.plot.track.plot_track

Use a fully qualified symbol.
```

### `hic_basic api where SYMBOL`

用途：

- 给 agent 源码位置，方便读实现或模仿调用。
- 不直接 dump 大段源码。

示例：

```bash
hic_basic api where cool2mat
hic_basic api where hic_basic.coolstuff.cool2mat --json
```

文本输出建议：

```text
hic_basic.coolstuff.cool2mat

Source:
  /installed/path/hic_basic/coolstuff.py:123

Import:
  from hic_basic.coolstuff import cool2mat

Open:
  sed -n '123,220p' /installed/path/hic_basic/coolstuff.py
```

实现细节：

- 优先用 `importlib.import_module()` 导入模块。
- 用 `inspect.getsourcefile()` 和 `inspect.getsourcelines()` 得到安装环境中的真实路径和行号。
- 如果导入失败，fallback 到索引里的相对路径和 lineno，并明确标注 source path could not be resolved。

### `hic_basic api map`

Phase 2 后实现。

用途：

- 展示任务导向 section。
- 引导 agent 从任务进入模块/函数。

无 registry 时不建议实现假的 map。可以先让 `hic_basic api map` 输出：

```text
API map requires hic_basic/data/api_registry.json.
Use instead:
  hic_basic api list
  hic_basic api search <keywords>
```

有 registry 后输出：

```text
hic_basic API map

Sections:
  cool          Read, write, convert, and inspect .cool/.mcool/scool data.
  io            Metadata and file I/O helpers.
  plot          Hi-C heatmaps, tracks, and visualization helpers.
  compartment   Compartment score and A/B related workflows.

Next:
  hic_basic api list cool
  hic_basic api search "plot cool tracks"
  hic_basic api show cool2mat
```

### `hic_basic api examples`

Phase 3 后实现。

示例只来自 registry overlay，不从 docstring 自动猜。

如果没有 example：

```text
No curated examples for hic_basic.coolstuff.cool2mat.
Use:
  hic_basic api show hic_basic.coolstuff.cool2mat
  hic_basic api where hic_basic.coolstuff.cool2mat
```

### `hic_basic api source`

Phase 3 后实现。

用途：

- 打印受限源码片段，例如函数定义到函数结束或最多 N 行。
- 默认仍推荐 `where`，因为 agent 可以自行 `sed`。

建议参数：

```bash
hic_basic api source cool2mat
hic_basic api source cool2mat --max-lines 120
```

## 实现阶段

### Phase 0: 整合现有文档系统

目标：让现有 API 文档系统成为 CLI 的数据源。

任务：

1. 抽出 `docs/search_api.py` 中的查询逻辑，放到包内模块，例如：

   ```text
   hic_basic/api_query.py
   ```

2. 保留 `docs/search_api.py`，改成 wrapper，继续支持旧命令：

   ```bash
   python docs/search_api.py cool2mat
   python docs/search_api.py --module coolstuff
   ```

3. 修改 `scripts/generate_api_reference.py`：

   - 继续生成 `docs/API_REFERENCE.md`。
   - 继续生成 `docs/api_index.json`。
   - 新增生成 `hic_basic/data/api_index.json`。
   - JSON index 不再截断完整 docstring。
   - 增加 `summary`、`symbol`、`import`、`source.module_path` 字段。

4. 修改打包配置：

   - 在 `setup.py` 增加 `include_package_data=True` 或 `package_data`。
   - 确保 `hic_basic/data/api_index.json` 被安装。
   - 如果新增 `hic_basic/data/api_registry.json`，也要被安装。

5. 更新文档统计：

   - 当前 README/docs 中的模块数和函数数不一致。
   - 最好由生成脚本写入或改成“不写死精确数字”。

验收：

```bash
./run_in_env.sh python scripts/generate_api_reference.py
./run_in_env.sh python docs/search_api.py cool2mat
./run_in_env.sh python - <<'PY'
from importlib.resources import files
p = files("hic_basic").joinpath("data/api_index.json")
print(p.is_file())
PY
```

### Phase 1: 安装后可用的最小 API CLI

目标：在下游项目中可以直接查 API。

任务：

1. 新增 `hic_basic/cli/api.py`。
2. 在 `hic_basic/cli/__init__.py` 注册 `api` 子命令。
3. 实现：

   ```bash
   hic_basic api list [MODULE] [--json]
   hic_basic api search QUERY... [--limit N] [--json]
   hic_basic api show SYMBOL [--brief|--signature] [--json]
   hic_basic api where SYMBOL [--json]
   ```

4. 用 `importlib.resources.files("hic_basic").joinpath("data/api_index.json")` 加载包内 `api_index.json`，避免导入较重的 `hic_basic.data` 子包。
5. 支持短名和全名 symbol lookup。
6. 支持 ambiguous symbol 候选输出。
7. `where` 使用 runtime inspect 解析实际安装路径。

验收：

```bash
./run_in_env.sh hic_basic api list
./run_in_env.sh hic_basic api list coolstuff
./run_in_env.sh hic_basic api search convert pairs
./run_in_env.sh hic_basic api show cool2mat
./run_in_env.sh hic_basic api show cool2mat --json
./run_in_env.sh hic_basic api where cool2mat
```

测试建议：

- 新增 `tests/test_api_cli.py`。
- 用 `subprocess.run([sys.executable, "-m", "hic_basic", "api", ...])` 测试，避免依赖 console script 是否安装。
- 覆盖：
  - `api list` 返回 0。
  - `api list coolstuff` 包含 `cool2mat`。
  - `api search cool2mat --json` 是合法 JSON。
  - `api show cool2mat --json` 包含 `symbol`、`signature`、`import`。
  - `api where cool2mat` 包含 `.py:` 或 JSON source 字段。
  - ambiguous 短名返回非 crash 的候选。

### Phase 2: registry overlay 和任务导航

目标：让 agent 不只按函数名查，也能从任务语义收敛。

任务：

1. 新增 `hic_basic/data/api_registry.json`。
2. 先只覆盖高频模块和函数，不追求全包完整。
3. 增强查询模型：

   - section
   - tags
   - inputs
   - outputs
   - examples
   - pitfalls
   - related

4. 实现：

   ```bash
   hic_basic api map [--json]
   hic_basic api list <SECTION>
   ```

5. `api search` 加入 registry 权重。
6. `api show` 在存在 registry metadata 时显示 inputs/outputs/tags/related。

初始 section 建议基于当前真实模块：

```text
cool          coolstuff 中的 cool/mcool/scool/pairs/pixels 转换和读取
io            hicio 和 metadata/file parsing
plot          plot.hic、plot.track、plot.utils 等可视化
compartment   compartment 和 plot.hic 中的 compartment 相关函数
compare       compare_pixels、metrics、correlation/similarity 相关函数
genome        genome、binnify、sequence 等基因组坐标/分箱
phasing       phasing、cycle_phasing、crossover 相关函数
structure     structure、territory、threeD/render 相关函数
```

注意：section 是任务导向，不必等于 Python module。

验收：

```bash
./run_in_env.sh hic_basic api map
./run_in_env.sh hic_basic api list cool
./run_in_env.sh hic_basic api search "plot cool tracks"
./run_in_env.sh hic_basic api search "pairs to scool"
```

### Phase 3: examples/source/schema 增强

目标：提高 agent 写脚本成功率，但不阻塞 MVP。

任务：

1. 实现 curated examples：

   ```bash
   hic_basic api examples cool2mat
   hic_basic api examples cool
   ```

2. 实现受限源码输出：

   ```bash
   hic_basic api source cool2mat --max-lines 120
   ```

3. 实现 schema alias：

   ```bash
   hic_basic api schema cool2mat
   hic_basic api schema --all
   ```

   `schema` 可以先等价于 `show --json`，后续再稳定扩展。

4. 增加更多高频 API 的 registry metadata 和 examples。

## 文档同步计划

引入 `hic_basic api` 后，要同步更新：

- `AGENTS.md`
- `README.md`
- `QUICK_REFERENCE.md`
- `docs/README.md`
- `docs/SUMMARY.md`

更新策略：

1. 新文档主推安装后命令：

   ```bash
   hic_basic api search cool2mat
   hic_basic api list coolstuff
   hic_basic api show cool2mat
   hic_basic api where cool2mat
   ```

2. 保留 repo 内旧命令作为 fallback：

   ```bash
   python docs/search_api.py cool2mat
   ```

3. 统一执行环境描述。根据 `AGENTS.md`，项目内推荐：

   ```bash
   ./run_in_env.sh hic_basic api search cool2mat
   ./run_in_env.sh python scripts/generate_api_reference.py
   ```

4. 解决统计数字不一致。不要在多个文档中手写不同的 “530/600+/656 functions”。推荐：

   - 只在生成文档里写精确统计。
   - README/QUICK_REFERENCE 写 “auto-generated public API index”。

5. 增加 agent workflow：

   ```text
   If you know the function:
     hic_basic api show <function>
     hic_basic api where <function>

   If you only know the task:
     hic_basic api search "<task keywords>"
     hic_basic api list <module-or-section>
     hic_basic api show <chosen symbol>

   If you need structured output:
     add --json
   ```

## 向后兼容

必须保持这些入口可用：

```bash
python docs/search_api.py cool2mat
python docs/search_api.py --module coolstuff
python docs/search_api.py --keyword convert
python docs/search_api.py --list-modules
```

`docs/search_api.py` 可以变薄，但不应删除。

如果未来决定弃用，应先在文档里标记 deprecated，并保留至少一个周期。

## 不做的事情

Phase 1 不做：

- embedding search
- web/RAG
- HTML 文档站
- 完整人工维护每个 API 的 examples/pitfalls
- `hic-basic` 或 `hb` 别名
- 基于虚构模块的测试，例如 `compare.scc`、`eigs.eigs_cis`

这些可以后续讨论，但不是当前最小闭环。

## 最终验收标准

在一个只安装了 `hic_basic` 的下游项目中，agent 应该能完成：

```bash
hic_basic api search "convert pairs scool"
hic_basic api show pairs2scool
hic_basic api where pairs2scool
```

并据此写出：

```python
from hic_basic.coolstuff import pairs2scool
```

在 repo 内，旧工作流仍然可用：

```bash
python docs/search_api.py pairs2scool
python scripts/generate_api_reference.py
```

在项目标准环境中，测试通过：

```bash
./run_in_env.sh pytest -q tests/test_api_cli.py
```
