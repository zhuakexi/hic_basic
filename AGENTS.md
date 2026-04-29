# hic_basic — Agent Notes

This repository is a Python library and CLI for Hi-C / 3D genome data processing. Use this file as the single source of repo-specific guidance for automation agents.

## Quick function discovery
- Search installed API by task/name: `hic_basic api search "convert pairs"`
- Show a known function: `hic_basic api show cool2mat`
- Locate source: `hic_basic api where cool2mat`
- Search by module: `hic_basic api list coolstuff`
- Repo-local fallback: `python docs/search_api.py cool2mat`
- Programmatic index: `docs/api_index.json`
- Full reference: `docs/API_REFERENCE.md`

## Repo layout (high signal)
- Core library: `hic_basic/`
- I/O and metadata helpers: `hic_basic/hicio.py`
- Conversions/utilities: `hic_basic/coolstuff.py`
- Plotting: `hic_basic/plot/hic.py`
- CLI entrypoint: `hic_basic/__main__.py` (console script `hic_basic` in `setup.py`)
- Tests: `tests/` (also a legacy `test/` directory exists)
- Docs tooling: `scripts/generate_api_reference.py`

## Environment and execution
- **Entrypoint script**: `./run_in_env.sh` — the only project-supported compute entrypoint
- **Run commands**: `./run_in_env.sh <command>` (e.g., `./run_in_env.sh python -m hic_basic`)
- **Install editable**: `./run_in_env.sh pip install -e .`
- **Run tests**: `./run_in_env.sh pytest -q tests`
- **Self-check**: `./run_in_env.sh --self-check` (verifies environment and container setup)
- **Under the hood**: Runs commands in a Singularity container (`light_base_1_1.sif`) with the `hic_basic_v096` mamba environment

## Documentation workflow
- Public functions are auto-indexed; keep docstrings accurate.
- Regenerate docs after adding/updating public APIs:
  - `./run_in_env.sh python scripts/generate_api_reference.py`
  - Commit both code changes and updated docs (`docs/API_REFERENCE.md`, `docs/api_index.json`).

## API command system discipline
- Treat `docs/api_index.json` and the future package-installed API index as generated artifacts from the same source; do not hand-edit either index.
- When changing public functions, keep signatures and docstrings agent-usable: include the purpose, important inputs, return value, and any non-obvious constraints.
- After adding, renaming, moving, or deleting public APIs, regenerate the API docs/index and verify lookup still works for the changed symbol.
- If modifying the API lookup command itself, preserve compatibility with the existing repo-local workflow: `python docs/search_api.py <query>`, `--module`, `--keyword`, and `--list-modules`.
- Prefer adding shared API lookup logic in importable library code and making `docs/search_api.py` plus `hic_basic api ...` thin command wrappers around it.
- Keep downstream usage independent of the repo checkout: `hic_basic api ...` must read package-installed resources, not assume the current working directory is the hic_basic repository.
- Do not add task-semantic registry entries for nonexistent or planned APIs; registry/examples must refer to real importable symbols unless explicitly marked as design examples in an issue document.
- Add or update CLI tests when changing API lookup behavior, including `--json` output and ambiguous-symbol handling where relevant.

## Agent workflow (library-first)
1. **Search existing APIs** in `hic_basic` before adding new code. Prefer extending a relevant module over creating a new one.
2. **Keep changes backwards compatible** whenever possible. If a break is unavoidable:
   - add a deprecation path or migration note
   - update tests to cover both old and new behavior if feasible
3. **Add tests with minimal data**:
   - use `tests/data` and tmp paths
   - avoid absolute paths and large files
   - when shrinking data, keep a reproducible command or note in the test
4. **Avoid project-specific logic** in this repo; put shared, general utilities only.

## Project conventions
- Prefer `read_meta()` from `hic_basic.hicio` for metadata ingestion instead of ad-hoc pandas reads.
- For plotting, `plot_cool` is the common entrypoint; watch `binsize` and `vmax` for visibility.
- Reuse conversion helpers in `coolstuff.py` for pairs/scool/cool workflows and parallelization patterns.
- Avoid hard-coded absolute paths in tests; use `tests/data` and temporary paths.

## Quick examples
```python
from hic_basic.hicio import read_meta
from hic_basic.coolstuff import pairs2scool
from hic_basic.plot.hic import plot_cool

meta = read_meta("xx.meta.csv.gz")
pairs2scool(meta["pairs_c12"].to_dict(), "xx.scool", "mm10.len.tsv", 20000)
plot_cool("xx.cool", "title", "chr19", vmax=500)
```
