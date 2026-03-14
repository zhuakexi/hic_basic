# hic_basic — Agent Notes

This repository is a Python library and CLI for Hi-C / 3D genome data processing. Use this file as the single source of repo-specific guidance for automation agents.

## Quick function discovery
- Search by name: `python docs/search_api.py cool2mat`
- Search by module: `python docs/search_api.py --module coolstuff`
- Search by keyword: `python docs/search_api.py --keyword "convert"`
- Interactive search: `python docs/search_api.py --interactive`
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
  - `python scripts/generate_api_reference.py`
  - Commit both code changes and updated docs (`docs/API_REFERENCE.md`, `docs/api_index.json`).

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
