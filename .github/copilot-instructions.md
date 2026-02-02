<!-- .github/copilot-instructions.md: Guidance for AI coding agents working on hic_basic -->
# hic_basic — Copilot instructions (concise)

This file gives focused, repository-specific guidance so an AI coding agent can be immediately productive.

0) API Documentation & Function Discovery
- **Quick function search**: Use `python docs/search_api.py <function_name>` to find functions and their signatures.
  - Example: `python docs/search_api.py cool2mat` — search by name
  - Example: `python docs/search_api.py --module coolstuff` — list all functions in a module
  - Example: `python docs/search_api.py --keyword "convert"` — search by keyword in docstrings
  - Example: `python docs/search_api.py --interactive` — interactive exploration mode
- **Complete API reference**: See `docs/API_REFERENCE.md` (527 functions across 76 modules, auto-generated)
- **Programmatic access**: Load `docs/api_index.json` for structured function metadata
- **Update documentation**: After adding/modifying public functions, run `python scripts/generate_api_reference.py` to regenerate docs, then commit both code and docs together.

1) Run commands
hic_basic relies on micromamba for dependency management. To run any test commands, you must prefix them with `micromamba run -n hic_basic`. The environment definition file is located at `envs/hic_basic.yaml`.

2) Quick summary / big picture
- hic_basic is a Python library and CLI for processing Hi-C / 3D genome data. Core responsibilities:
  - I/O and metadata helpers: `hic_basic/hicio.py` (notably `read_meta`) and `hic_basic/data.py`.
  - Conversion and file utilities: `hic_basic/coolstuff.py` (functions: `pairs2cool`, `pairs2scool`, `hic_pileup`).
  - Plotting and visualization: `hic_basic/plot/hic.py` (`plot_cool`, `plot_cool_track`, `plot_cools`).
  - CLI wrappers: `hic_basic/cli/*` (e.g. `render.py`, `download.py`).
  - Higher-level analysis modules: many modules under `hic_basic/` (e.g. `impute/`, `wet/`, `pseudotime/`).

3) Why things are organized this way
- The codebase groups domain concerns: file-format conversions and caching live in `coolstuff.py`, visualization in `plot/`, and sample/meta handling in `hicio.py`. This keeps I/O and heavy data transforms separated from plotting and analysis logic.

4) Key files to read first (order matters)
- `hic_basic/hicio.py` — metadata readers and lightweight I/O utilities (call sites across CLI and tests).
- `hic_basic/coolstuff.py` — conversion and batch utilities; examples in the README and used by `wet/` utilities.
- `hic_basic/plot/hic.py` — plotting primitives and example usages (tests under `tests/test_plot.py`).
- `hic_basic/cli/render.py` — shows how CLI commands construct and consume metadata and call plot/render flows.
- `setup.py` — package entrypoint (`hic_basic` console script) and packaging hints.

5) Typical developer workflows (commands you can rely on)
- Create a Python environment: a conda YAML is provided in `envs/` (`hic_basic_py3.10_.yaml`, `hic_basic_py3.7_.yaml`). Prefer these for reproducibility.
  - Example: `conda env create -f envs/hic_basic_py3.10_.yaml` (or adapt to pip/venv if you prefer).
- Install editable package: `pip install -e .` so `python -m hic_basic` and `import hic_basic` use workspace code.
- Run unit tests: `pytest -q tests` or `python -m pytest tests` (tests use pytest and reference test-data under `tests/data`).
- Run CLI entrypoint locally: `python -m hic_basic` (maps to `hic_basic.__main__:main`). After install, the `hic_basic` console script is available.

6) Project-specific conventions and patterns
- Metadata: `read_meta()` is the canonical helper for reading compressed CSV/TSV metadata used by many modules and CLI functions — prefer it over ad-hoc pandas reads.
- File types: code expects Cooler (.cool/.mcool), scool (multi-sample cool), and pairs files. See `coolstuff.py` for conversion utilities and typical parameters.
- Plotting: `plot_cool` is the common entrypoint; callers pass either chromosome names ("chr1") or region slices. Watch `vmax` and `binsize` — tests confirm expected defaults.
- Long-running/parallel helpers: `mt_pairs2cool` and worker patterns in `coolstuff.py` — prefer reusing these helpers for parallelizing conversion tasks.
- Tests and examples sometimes reference absolute paths in older tests; prefer using `request.fspath` fixtures or test-data under `tests/data` for new tests.

7) Integration points & external dependencies
- External file formats and tools: cooler, cooler-tools ecosystems are expected (mcool, .cool). Also common bioinformatics formats (pairs). Check `envs/*.yaml` for exact packages to install.
- Notebooks: the `notebooks/` directory contains runnable examples. Some notebooks mutate `sys.path` to include `hires_utils` — if you see this, emulate the same pattern when running notebooks locally.

8) Small extractable examples (copyable references)
- Convert pairs table to scool (README & `coolstuff.py`):
  - `from hic_basic.coolstuff import pairs2scool; pairs2scool(filesp["pairs_c12"].to_dict(), "xx.scool", "mm10.len.tsv", 20000)`
- Read metadata (used throughout):
  - `from hic_basic.hicio import read_meta; meta = read_meta("xx.meta.csv.gz")`
- Plot a cool file (used in tests):
  - `from hic_basic.plot.hic import plot_cool; fig = plot_cool("xx.cool", "title", region="chr1", vmax=500)`

9) When editing code, watch for these fragile areas
- Hard-coded file paths or absolute paths in tests — prefer using `request.fspath` fixtures or test-data under `tests/data`.
- Plotting and I/O functions may assume cooler-backed arrays; avoid introducing heavy in-memory copies in hot code paths.

10) Where to add tests / quick validation
- Add unit tests under `tests/` alongside `test_*.py`. Look at `tests/test_plot.py` for plotting test patterns. Use small synthetic data files from `tests/data/` or create temporary files.

11) Documentation workflow (IMPORTANT)
- **When adding public functions**: Write clear docstrings with parameters, returns, and examples using type hints.
- **After code changes**: Run `python scripts/generate_api_reference.py` to update API documentation (regenerates `docs/API_REFERENCE.md` and `docs/api_index.json`).
- **Commit workflow**: Always commit both code changes AND updated documentation together:
  ```bash
  # 1. Make your changes
  # 2. Regenerate docs
  python scripts/generate_api_reference.py
  # 3. Commit together
  git add hic_basic/ docs/
  git commit -m "Add feature with updated docs"
  ```
- **Docstring format**: Follow NumPy/SciPy style with Parameters, Returns, and Examples sections. See existing functions for patterns.

12)  If you need more context or can't infer intent
- Point me to the specific module/PR and I will: (a) open the top-level functions used by callers, (b) list call sites across the repo, and (c) propose minimal tests demonstrating intended behavior.

-- End of file
