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
- Conda environment definition: `envs/hic_basic.yaml` (name: `hic_basic`, Python 3.10)
- Install editable for local use: `pip install -e .`
- Run tests: `pytest -q tests`

## Documentation workflow
- Public functions are auto-indexed; keep docstrings accurate.
- Regenerate docs after adding/updating public APIs:
  - `python scripts/generate_api_reference.py`
  - Commit both code changes and updated docs (`docs/API_REFERENCE.md`, `docs/api_index.json`).

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
