# QWEN.md — hic_basic

This file provides guidance for Qwen when working in the hic_basic codebase.

## Quick function discovery
- **Search by name**: `python docs/search_api.py cool2mat`
- **Search by module**: `python docs/search_api.py --module coolstuff`
- **Search by keyword**: `python docs/search_api.py --keyword "convert"`
- **Interactive search**: `python docs/search_api.py --interactive`
- **Programmatic index**: `docs/api_index.json`
- **Full reference**: `docs/API_REFERENCE.md` (527 functions across 76 modules)

## Environment and execution
- **Entrypoint script**: `./run_in_env.sh` — the only project-supported compute entrypoint
- **Run commands**: `./run_in_env.sh <command>` (e.g., `./run_in_env.sh python -m hic_basic`)
- **Install editable**: `./run_in_env.sh pip install -e .`
- **Run tests**: `./run_in_env.sh pytest -q tests`
- **Self-check**: `./run_in_env.sh --self-check` (verifies environment and container setup)
- **Under the hood**: Runs commands in a Singularity container (`light_base_1_1.sif`) with the `hic_basic_v096` mamba environment

## Repo layout
- **Core library**: `hic_basic/`
- **I/O and metadata**: `hic_basic/hicio.py` (`read_meta`), `hic_basic/data.py`
- **Conversions**: `hic_basic/coolstuff.py` (`pairs2cool`, `pairs2scool`, `hic_pileup`)
- **Plotting**: `hic_basic/plot/hic.py` (`plot_cool`, `plot_cool_track`, `plot_cools`)
- **CLI**: `hic_basic/cli/` (e.g., `render.py`, `download.py`)
- **Analysis modules**: `hic_basic/impute/`, `hic_basic/wet/`, `hic_basic/pseudotime/`
- **Tests**: `tests/` (prefer `tests/data` for fixtures)
- **Docs tooling**: `scripts/generate_api_reference.py`

## Documentation workflow
When adding or modifying public functions:
1. Write clear docstrings with Parameters, Returns, and Examples (NumPy/SciPy style)
2. Run `python scripts/generate_api_reference.py` to regenerate docs
3. Commit both code changes AND updated documentation together:
   ```bash
   git add hic_basic/ docs/
   git commit -m "Add feature with updated docs"
   ```

## Agent workflow (library-first)
1. **Search existing APIs** in `hic_basic` before adding new code. Prefer extending a relevant module over creating a new one.
2. **Keep changes backwards compatible** whenever possible. If a break is unavoidable:
   - Add a deprecation path or migration note
   - Update tests to cover both old and new behavior if feasible
3. **Add tests with minimal data**:
   - Use `tests/data` and tmp paths
   - Avoid absolute paths and large files
4. **Avoid project-specific logic** in this repo; put shared, general utilities only.

## Project conventions
- **Metadata**: Prefer `read_meta()` from `hic_basic.hicio` over ad-hoc pandas reads
- **File types**: Works with Cooler (.cool/.mcool), scool, and pairs files
- **Plotting**: `plot_cool` is the common entrypoint; watch `binsize` and `vmax`
- **Parallel patterns**: Reuse `mt_pairs2cool` and worker patterns from `coolstuff.py`
- **Tests**: Use `request.fspath` fixtures or `tests/data`; avoid hard-coded absolute paths

## Common examples
```bash
# Run a command in the environment
./run_in_env.sh python -m hic_basic --help

# Install the package in editable mode
./run_in_env.sh pip install -e .

# Run tests
./run_in_env.sh pytest -q tests

# Self-check to verify environment
./run_in_env.sh --self-check
```

```python
# Read metadata
from hic_basic.hicio import read_meta
meta = read_meta("xx.meta.csv.gz")

# Convert pairs to scool
from hic_basic.coolstuff import pairs2scool
pairs2scool(meta["pairs_c12"].to_dict(), "xx.scool", "mm10.len.tsv", 20000)

# Plot a cool file
from hic_basic.plot.hic import plot_cool
fig = plot_cool("xx.cool", "title", region="chr1", vmax=500)
```

## Fragile areas to watch
- Hard-coded file paths in tests — use fixtures or `tests/data`
- Plotting/I/O functions assume cooler-backed arrays; avoid heavy in-memory copies
- Notebooks may mutate `sys.path` to include `hires_utils`

## External dependencies
- **Core**: cooler, cooler-tools ecosystems (.cool, .mcool)
- **Bioinformatics**: pairs format
- **Package list**: See `envs/*.yaml` for exact dependencies
