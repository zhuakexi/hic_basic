# hic_basic

A Python library for processing Hi-C and 3D genome data.

## 📚 Documentation

**Complete API Reference**: [docs/API_REFERENCE.md](docs/API_REFERENCE.md)
- 530 public functions across 77 modules
- Auto-generated from source code
- Searchable with `python docs/search_api.py`

**Agent guidance**: [AGENTS.md](AGENTS.md)

Quick examples:
```bash
# Search for a function
python docs/search_api.py cool2mat

# List module functions  
python docs/search_api.py --module coolstuff

# Interactive mode
python docs/search_api.py --interactive
```

See [docs/README.md](docs/README.md) for detailed usage.

## Installation

In jupyter notebooks, add
```python
import sys
sys.path.insert(0, "$YourGitDir/hic_basic")
sys.path.insert(0, "$YourGitDir/hires_utils")
```

## Compatibility Policy (Downstream Projects)

- Downstream projects should **pin a hic_basic commit** (or tag) in their `README.md` or config.
- We aim for **backwards-compatible** changes by default.
- If a breaking change is unavoidable:
  - provide a deprecation window or migration note
  - add or update tests to cover legacy and new behavior when feasible
  - call out the change in release notes or PR description

## Contributing Checklist (PRs)

- Search existing APIs before adding new ones (prefer extending a module over creating a new one).
- Add tests that use **minimal data** from `tests/data` and tmp paths only.
- Avoid hard-coded absolute paths or project-specific logic.
- Update docstrings and regenerate API docs when public APIs change:
  - `python scripts/generate_api_reference.py`
  - commit `docs/API_REFERENCE.md` and `docs/api_index.json`

## Quick Start Examples

### Pseudo bulk analysis
pairs --> .scool --> .cool(s)  

Assume you have a meta file `xx.meta.csv.gz` that stores all pairs file path and sample name in `pairs_c12` col.  
To generate a `xx.scool` file that stores all samples, you can use the following command:  

```python
from hic_basic.coolstuff import pairs2scool, hic_pileup
from hic_basic.hicio import read_meta

filesp = read_meta("xx.meta.csv.gz")
pairs2scool(filesp["pairs_c12"].to_dict(), "xx.scool", "mm10.len.tsv", 20e3)
hic_pileup("xx.scool", {sample: cell_type, ...}, "{}.pileup.cool")
```

you need a chromosome length file, two columns, no header.  
you also need a dict that maps sample name to cell types.

### Plot Hi-C matrix

vmax and binsize are super important here.
当只能看到一条对角线时，用大binsize和小vmax。

```python
from hic_basic.plot.hic import plot_cool

# This command plot whole genome heatmap with all intra and inter chromosomal interactions.
# Make sure you use a low resolution cooler file here like 1Mb binsize for mouse genome.
plot_cool("xx.cool", "title", [slice(0,-1), slice(0,-1)], vmax=500)

# This works for mcool file with multiple resolutions.
plot_cool("xx.mcool::resolutions/1000000", "title", [slice(0,-1), slice(0,-1)], vmax=500)

# This plot inter chromosomal interaction of chr1 and chr2.
plot_cool("xx.cool", "title", ["chr1", "chr2"], vmax=500)

# This plot intra chromosome interaction of chr19
plot_cool("xx.cool", "title", "chr19", vmax=500)
```

## API Reference

For complete function documentation, see [docs/API_REFERENCE.md](docs/API_REFERENCE.md) or use the search tool:

```bash
python docs/search_api.py <function_name>
```
# plot compartment
