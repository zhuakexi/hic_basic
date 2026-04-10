# 公共信息统一用 Ana 管理
- 用 `ana.update(df|series|dict|list, key=..., description="...")` 写入数据，**禁止手动修改 commit 文件**。
- 首次 update 自动创建 commit；相同数据重复调用自动跳过（幂等），`force=True` 可强制记录检查点。
- 新列通过 `column_descriptions={"col": "短描述"}`、`depends_on={"col": {"columns": [...], "objects": [...]}}`、`computed_by={"file": "script.ipynb"}` 声明元数据。
- 文件路径列自动识别类型（cool/bam/bed 等），也可 `file_format={"col": "bam"}` 显式指定。

# 分析前查询有用的 Ana 输入
- **CLI 准备/核查**：`hic_basic ana-describe <home>` 查看列名+类型+描述；`--column <key>` 打印完整元数据；`--file-format cool` 查看格式定义与 reader/writer。
- **代码中**：`ana.describe()` 打印列与对象清单；`ana.data` 取 DataFrame；`ana.obj` 取字典；`ana.column_meta("col")` 查依赖/来源；`ana.file_format("cool")` 查读取函数。
