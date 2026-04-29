# Quick Reference Card - hic_basic API Documentation

## 🚀 Quick Commands

```bash
# Search installed/package API
hic_basic api search cool2mat

# Show a known function
hic_basic api show cool2mat

# List all functions in a module
hic_basic api list coolstuff

# Locate source
hic_basic api where cool2mat

# Repo-local fallback
python docs/search_api.py cool2mat

# Regenerate documentation after code changes
./run_in_env.sh python scripts/generate_api_reference.py
```

## 📚 Files

| File | Use When |
|------|----------|
| `docs/API_REFERENCE.md` | Reading in editor/browser |
| `docs/api_index.json` | Programming/scripting |
| `hic_basic api ...` | Installed-package quick lookups |
| `docs/search_api.py` | Repo-local fallback lookups |
| `docs/EXAMPLES.py` | Learning the system |

## 🔥 Most Used Modules

| Module | Functions | Purpose |
|--------|-----------|---------|
| **coolstuff** | 32 | Cooler file operations |
| **hicio** | 31 | I/O and file parsing |
| **plot.render** | 27 | Visualization rendering |
| **phasing** | 21 | Cell cycle analysis |
| **plot.hic** | 22 | Hi-C heatmap plotting |
| **plot.track** | 17 | Track-based visualization (Heatmap, Eigen, IS, Gene) |
| **plot.utils** | 20 | Plotting utilities |
| **data** | 18 | Data structures |
| **TI.slingshot** | 22 | Trajectory inference (Slingshot) |
| **TI.slingshot.plot** | 3 | Plotly visualization (PNG) |
| **TI.circle** | 1 | Pseudotime baseline (circular) |
| **plot.plotly_template** | 1 | Plotly style template |
| **interpolate** | 8 | GP interpolation & grid methods |
| **plot.gene_track_parser** | 6 | GTF/GFF3 parsing for gene annotation |

## 💻 Code Snippets

### Python: Load API Index
```python
import json
with open('docs/api_index.json') as f:
    api = json.load(f)
print(f"Total functions: {sum(len(f) for f in api.values())}")
```

### Python: Search Functions
```python
# Find functions by name
for module, funcs in api.items():
    for func in funcs:
        if 'plot' in func['name']:
            print(f"{module}.{func['name']}")
```

### Python: Get Function Details
```python
# Get specific function
coolstuff_funcs = api['coolstuff']
cool2mat = next(f for f in coolstuff_funcs if f['name'] == 'cool2mat')
print(cool2mat['signature'])
print(cool2mat['docstring'])
```

### Bash: Quick Stats
```bash
# Count total functions
python -c "import json; d=json.load(open('docs/api_index.json')); print(sum(len(f) for f in d.values()))"

# List module names
python -c "import json; print('\n'.join(json.load(open('docs/api_index.json')).keys()))"
```

## 🎯 Common Tasks

### Find a Function
```bash
hic_basic api show <function_name>
```

### Explore a Module  
```bash
hic_basic api list <module_name>
```

### Search by Capability
```bash
hic_basic api search "convert"
hic_basic api search "plot"
hic_basic api search "matrix"
```

### Slingshot Trajectory Inference
```python
from hic_basic.TI.slingshot import Slingshot
from hic_basic.TI.slingshot.plot import plot_slingshot

# Fit Slingshot
slingshot = Slingshot(data=adata, celltype_key='celltype', start_node='SSC')
slingshot.fit(num_epochs=10)

# Visualize (PNG output)
plot_slingshot(slingshot, save_path="trajectory.png")

# Plot clusters by cell type
from hic_basic.TI.slingshot.plot import plot_clusters
plot_clusters(slingshot, save_path="clusters.png")

# Plot principal curves only
from hic_basic.TI.slingshot.plot import plot_curves
plot_curves(slingshot, save_path="curves.png")
```

### Circular Pseudotime (Baseline)
```python
from hic_basic.TI.circle import angle_pseudotime

# Simple circular ordering from 2D embedding
pseudotime_df = angle_pseudotime(root=0, embedding=df_2d)
```

### GP Interpolation with Grid Methods
```python
from hic_basic.interpolate import gp_interpolate, plot_gp_heatmap

# Run GP interpolation with different grid methods
result = gp_interpolate(
    adata,
    pseudotime_col='pt',
    numPts=14,
    grid_method='linspace',      # or 'quantile', 'density', 'optimal', 'custom'
    endpoint_mode='half-interval',  # avoid sparse cells at endpoints
    inplace=False
)

# Access results
grid_data = result.uns['gp_interpolate']
mean_expr = grid_data['mean']      # genes x grid_points
variance = grid_data['variance']
time_grid = grid_data['time_grid']

# Visualize single method result (publication-ready)
fig = plot_gp_heatmap(
    adata,
    pseudotime_col='pt',
    numPts=12,
    grid_method='linspace',        # single method
    endpoint_mode='half-interval',
    scaling='z_score',             # or 'zero_one'
    figsize=(8, 6),
    show_grid=True,                # show grid point positions (default)
    save_path='gp_heatmap.png'
)
```

### Grid Methods Comparison

| Method | Use Case | Endpoint Handling |
|--------|----------|-------------------|
| `linspace` | Uniform sampling | half-interval avoids t=0, t=1 |
| `quantile` | Equal cells per bin | Automatically handles endpoints |
| `density` | Dense regions | half-interval truncation |
| `optimal` | Data-driven (needs optbinning) | Automatic |
| `custom` | User-defined grid | User responsibility |

**Note**: For comparing multiple grid methods during debugging, use the internal function `_plot_grid_methods_comparison()` (not documented, for development use only).

### endpoint_mode Options
- `"half-interval"` (default): Endpoints at [half_interval, 1-half_interval]
  - Avoids GP extrapolation at t=0 and t=1 where only 1 cell exists
- `"full"`: Legacy behavior, endpoints at [0, 1]

### Plotly Template (Custom Style)
```python
from hic_basic.plot.plotly_template import apply_template, template

# Apply hic_basic style to any plotly figure
fig = go.Figure()
# ... add traces ...
fig = apply_template(fig, title="My Plot", width=800, height=600)
fig.write_image("output.png")
```

### Track Architecture (Modular Visualization)
```python
from hic_basic.plot.track import HeatmapTrack, EigenTrack, ISTrack, GeneTrack, TrackComposer

# Create individual tracks
heatmap = HeatmapTrack(cool_path="sample.cool", region="chr1:0-10000000")
eigen = EigenTrack(eig_file="eigvec.bedgraph")
is_track = ISTrack(insulation_file="IS.csv")
gene = GeneTrack(gtf_file="annotation.gtf", region="chr1:0-10000000")

# Compose multi-track figure
composer = TrackComposer(tracks=[heatmap, eigen, gene])
fig = composer.render(title="Multi-track View")
fig.write_html("output.html")
```

### Compartment Trapezoid Plot (Rotated Heatmap)
```python
from hic_basic.plot.hic import plot_compartment_trapezoid, plot_compartment_trapezoid_2tracks

# Single trapezoid heatmap
fig = plot_compartment_trapezoid(
    coolp="sample.cool",
    vec_file="eigvec.bedgraph",
    region="chr1",
    strata=50,  # height of trapezoid
    title="chr1 compartment"
)
fig.show()

# Two-track comparison (trapezoid + 2 eigenvector tracks)
fig = plot_compartment_trapezoid_2tracks(
    coolp0="sperm.cool",
    vec_file1="sperm_eigvec.bedgraph",
    vec_file2="esc_eigvec.bedgraph",
    region="chr1",
    strata=50,
    title="Sperm vs ESC compartment"
)
fig.write_html("compartment_comparison.html")
```

### Tiling Compartment Plot (Track Architecture v2)
```python
from hic_basic.plot.hic import plot_tiling_compartment

# 2x2 tiling layout: eigen tracks + combined heatmap
fig = plot_tiling_compartment(
    coolps=["sample1.cool", "sample2.cool"],
    eigs_files=["eigs1.bedgraph", "eigs2.bedgraph"],
    region="chr1:1000000-5000000",
    title="Compartment Comparison",
    corr=True  # show correlation matrix
)
fig.write_html("tiling_compartment.html")
```

### Hi-C Heatmap with Tracks (Track Architecture v3)
```python
from hic_basic.plot.hic import plot_cool_with_tracks
from hic_basic.plot import BigWigTrack, GeneTrack, EigenTrack

# Single sample with tracks
fig = plot_cool_with_tracks(
    coolp1="sample.cool",
    tracks=[
        EigenTrack("eigvec.bedgraph", eigen_col="E1"),
        BigWigTrack("H3K27ac.bw", bin_size=10000),
        GeneTrack("annotations.gtf", mode="all_transcripts")
    ],
    region="chr1:0-500000",
    title="Hi-C + Multi-omics Tracks"
)
fig.write_image("hic_with_tracks.png", width=800, height=1000)

# Two samples tiling comparison with tracks
# Lower triangle: sample1, Upper triangle: sample2
fig = plot_cool_with_tracks(
    coolp1="sample1.cool",
    coolp2="sample2.cool",
    tracks=[
        EigenTrack("eigvec.bedgraph", eigen_col="E1"),
        BigWigTrack("signal.bw", bin_size=10000)
    ],
    region="chr1:0-500000",
    title="Sample1 vs Sample2 (Tiling)",
    adjust_brightness=True,  # Auto-balance signal intensities
    cmap="RdBu_r"
)
fig.write_image("tiling_comparison.png", width=900, height=1100)

# Square heatmap option (automatically calculates height for 1:1 aspect ratio)
fig = plot_cool_with_tracks(
    coolp1="sample.cool",
    tracks=[
        BigWigTrack("H3K27ac.bw", bin_size=10000, display_name="H3K27ac"),
        GeneTrack("annotations.gtf", mode="canonical")
    ],
    region="chr1:0-500000",
    square_heatmap=True,  # Make heatmap approximately square
    track_heights=[1.0, 4.0],  # Adjust track heights (BigWig, Gene)
    heatmap_height=12
)
fig.write_image("square_heatmap.png", width=1100, height=1600)
```

**Key Features:**
- **Tiling mode**: When `coolp2` is provided, creates composite heatmap with lower triangle from `coolp1` and upper triangle from `coolp2`
- **Track types**: Supports `EigenTrack`, `BigWigTrack` (with `display_name`), `GeneTrack`, `ISTrack`, etc.
- **Brightness adjustment**: `adjust_brightness=True` auto-scales the brighter sample for better visual comparison
- **Mark pixels**: Use `marks` parameter to highlight specific position pairs with circles or rectangles
- **Square heatmap**: `square_heatmap=True` automatically calculates height for 1:1 aspect ratio
- **Track labels**: BigWigTrack `display_name` parameter sets the track label (defaults to filename)
- **Gene track**: Automatically staggers overlapping gene labels across multiple rows

### Mark Pixels on Heatmap
```python
from hic_basic.plot import plot_cool_with_tracks, MarkTrack, HeatmapTrack

# Method 1: Using plot_cool_with_tracks marks parameter
fig = plot_cool_with_tracks(
    coolp1="sample.cool",
    region="chr1:0-500000",
    marks=[
        {
            "chrom1": "chr1", "start1": 50000, "end1": 60000,
            "chrom2": "chr1", "start2": 150000, "end2": 160000,
            "color": "red",      # Mark color
            "shape": "circle",   # "circle", "rect", or "line"
            "width": 2,          # Line width
            "label": "Loop"      # Optional text label
        }
    ]
)
fig.write_image("marked_heatmap.png", width=800, height=800)

# Method 2: Using HeatmapTrack with marks
heatmap = HeatmapTrack(
    cool_path="sample.cool",
    marks=[
        {"chrom1": "chr1", "start1": 50000, "end1": 60000,
         "chrom2": "chr1", "start2": 150000, "end2": 160000,
         "color": "blue", "shape": "rect", "width": 3}
    ]
)

# Method 3: Using standalone MarkTrack
mark = MarkTrack(
    marks=[
        {"chrom1": "chr1", "start1": 50000, "end1": 60000,
         "chrom2": "chr1", "start2": 150000, "end2": 160000,
         "color": "red", "shape": "circle"}
    ],
    default_shape="circle",
    default_color="red"
)
```

**Mark Options:**
| Option | Values | Description |
|--------|--------|-------------|
| `shape` | `circle`, `rect`, `line` | Shape type |
| `color` | Any CSS color | Mark color |
| `width` | Integer | Line width (default: 2) |
| `label` | String | Optional text label (positioned outside shape) |

### Gene Track Visualization (GTF/GFF3)
```python
from hic_basic.plot.track import plot_gene_track, plot_gene_track_with_heatmap

# Plot gene annotation track
fig = plot_gene_track(
    annotation_file="gencode.gtf",
    region="chr1:100000-200000",
    mode="all_transcripts",  # or "gene_only"
    title="Gene Annotation",
    height=400,
    width=1000
)
fig.write_html("gene_track.html")

# Combined gene track + heatmap
fig = plot_gene_track_with_heatmap(
    cool_path="sample.cool",
    annotation_file="gencode.gtf",
    region="chr1:100000-200000",
    title="Hi-C + Gene Annotation"
)
fig.show()
```

### Parse GTF/GFF3 Files
```python
from hic_basic.plot.gene_track_parser import parse_annotation_file, filter_genes_in_region

# Auto-detect format (GTF or GFF3)
genes = parse_annotation_file("annotation.gtf")

# Filter genes in specific region
region_genes = filter_genes_in_region(
    genes,
    seqid="chr1",
    start=100000,
    end=200000,
    padding=10000
)

# Access gene structure
for gene_id, gene in region_genes.items():
    print(f"{gene_id}: {gene.start}-{gene.end}")
    for tx in gene.transcripts:
        print(f"  Transcript: {tx.transcript_id}")
        print(f"    Exons: {len(tx.exons)}")
        print(f"    CDS: {len(tx.cds)}")
```

### Update Documentation
```bash
# 1. Edit code and docstrings
# 2. Regenerate
./run_in_env.sh python scripts/generate_api_reference.py
# 3. Commit both
git add hic_basic/ docs/
git commit -m "Update code and docs"
```

### QC Filtering with Flexible Thresholds
```python
from hic_basic.wet.qc import c_u_filter

# Standard filtering with both bounds
qc_cells = c_u_filter(
    meta=adata.obs,
    ccol='contacts',
    ucol='umis',
    c_threshold=[150_000, 1_000_000],  # [lower, upper]
    u_threshold=[10_000, 600_000]
)

# Only upper bound (no lower limit)
qc_cells = c_u_filter(
    meta=adata.obs,
    ccol='contacts',
    ucol='umis',
    c_threshold=[None, 1_000_000],  # Only filter high contacts
    u_threshold=[None, 600_000]     # Only filter high UMIs
)

# Only lower bound (no upper limit)
qc_cells = c_u_filter(
    meta=adata.obs,
    ccol='contacts',
    ucol='umis',
    c_threshold=[150_000, None],  # Only filter low contacts
    u_threshold=[10_000, None]    # Only filter low UMIs
)
```

## 📖 Documentation Best Practices

```python
def example_function(param1: str, param2: int = 10) -> pd.DataFrame:
    """
    Brief one-line description.
    
    Longer explanation if needed.
    
    Parameters
    ----------
    param1 : str
        Description of param1
    param2 : int, optional
        Description of param2 (default: 10)
    
    Returns
    -------
    pd.DataFrame
        Description of return value
    
    Examples
    --------
    >>> result = example_function("test", 20)
    """
    pass
```

## 🔗 Quick Links

- **Full API Reference**: [docs/API_REFERENCE.md](docs/API_REFERENCE.md)
- **System Guide**: [docs/README.md](docs/README.md)
- **Examples**: [docs/EXAMPLES.py](docs/EXAMPLES.py)
- **Implementation**: [IMPLEMENTATION_SUMMARY.md](IMPLEMENTATION_SUMMARY.md)

## 📊 Statistics

- Auto-generated public API index
- Machine-readable JSON for tools and agents
- **Auto-generated** from source

## 🎓 Learning Path

1. **Start**: Read [docs/README.md](docs/README.md)
2. **Explore**: Run `python docs/EXAMPLES.py`
3. **Search**: Try `hic_basic api search "<task keywords>"`
4. **Reference**: Browse [docs/API_REFERENCE.md](docs/API_REFERENCE.md)
5. **Use**: Import and call functions!

---

**Questions?** See [docs/SUMMARY.md](docs/SUMMARY.md) for complete details.
