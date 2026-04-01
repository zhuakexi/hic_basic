# hic_basic

A Python library for processing Hi-C and 3D genome data.

## 📚 Documentation

**Complete API Reference**: [docs/API_REFERENCE.md](docs/API_REFERENCE.md)
- 600+ public functions across 85 modules
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

### Trajectory Inference (Single-cell Analysis)

```python
from hic_basic.TI.slingshot import Slingshot
from hic_basic.TI.slingshot.plot import plot_slingshot, plot_clusters

# Initialize and fit Slingshot
slingshot = Slingshot(
    data=adata,  # AnnData with obsm['X_umap']
    celltype_key='celltype',
    start_node='SSC'
)
slingshot.fit(num_epochs=10)

# Visualize trajectory (PNG output)
plot_slingshot(slingshot, save_path="trajectory.png")
plot_clusters(slingshot, save_path="clusters.png")
```

See `TI.slingshot` and `TI.slingshot.plot` modules for details.

### Hi-C Visualization with Multi-omics Tracks

```python
from hic_basic.plot.hic import plot_cool_with_tracks
from hic_basic.plot import BigWigTrack, GeneTrack, EigenTrack

# Single sample with tracks
fig = plot_cool_with_tracks(
    coolp1="sample.cool",
    tracks=[
        EigenTrack("eigvec.bedgraph"),
        BigWigTrack("H3K27ac.bw"),
        GeneTrack("annotations.gtf")
    ],
    region="chr1:0-500000"
)
fig.write_image("hic_tracks.png")

# Two samples tiling comparison
fig = plot_cool_with_tracks(
    coolp1="sample1.cool",
    coolp2="sample2.cool",
    tracks=[EigenTrack("eigvec.bedgraph")],
    region="chr1:0-500000",
    adjust_brightness=True  # Auto-balance intensities
)
fig.write_image("tiling.png")

# Mark specific pixels (loops, TAD boundaries, etc.)
fig = plot_cool_with_tracks(
    coolp1="sample.cool",
    region="chr1:0-500000",
    marks=[
        {
            "chrom1": "chr1", "start1": 50000, "end1": 60000,
            "chrom2": "chr1", "start2": 150000, "end2": 160000,
            "color": "red", "shape": "circle", "width": 2,
            "label": "Loop"
        }
    ]
)
fig.write_image("marked.png")
```

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
