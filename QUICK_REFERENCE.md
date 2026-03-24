# Quick Reference Card - hic_basic API Documentation

## 🚀 Quick Commands

```bash
# Search for a function by name
python docs/search_api.py cool2mat

# Search by keyword in docstrings
python docs/search_api.py --keyword "convert"

# List all functions in a module
python docs/search_api.py --module coolstuff

# List all available modules
python docs/search_api.py --list-modules

# Interactive search mode
python docs/search_api.py --interactive

# Regenerate documentation after code changes
python scripts/generate_api_reference.py
```

## 📚 Files

| File | Use When |
|------|----------|
| `docs/API_REFERENCE.md` | Reading in editor/browser |
| `docs/api_index.json` | Programming/scripting |
| `docs/search_api.py` | Quick lookups |
| `docs/EXAMPLES.py` | Learning the system |

## 🔥 Most Used Modules

| Module | Functions | Purpose |
|--------|-----------|---------|
| **coolstuff** | 32 | Cooler file operations |
| **hicio** | 31 | I/O and file parsing |
| **plot.render** | 27 | Visualization rendering |
| **phasing** | 21 | Cell cycle analysis |
| **plot.utils** | 20 | Plotting utilities |
| **plot.hic** | 18 | Hi-C heatmap plotting |
| **data** | 18 | Data structures |
| **TI.slingshot** | 22 | Trajectory inference (Slingshot) |
| **TI.slingshot.plot** | 3 | Plotly visualization (PNG) |
| **TI.circle** | 1 | Pseudotime baseline (circular) |
| **plot.plotly_template** | 1 | Plotly style template |
| **interpolate** | 8 | GP interpolation & grid methods |

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
python docs/search_api.py <function_name>
```

### Explore a Module  
```bash
python docs/search_api.py --module <module_name>
```

### Search by Capability
```bash
python docs/search_api.py --keyword "convert"
python docs/search_api.py --keyword "plot"
python docs/search_api.py --keyword "matrix"
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
from hic_basic.interpolate import gp_interpolate, plot_grid_methods_comparison

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

# Visualize and compare multiple grid methods
fig = plot_grid_methods_comparison(
    adata,
    pseudotime_col='pt',
    numPts=12,
    grid_methods=['linspace', 'quantile', 'density'],
    endpoint_mode='half-interval',  # recommended
    scaling='z_score',              # or 'zero_one'
    figsize=(10, 10),
    save_path='grid_comparison.png'
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

### Update Documentation
```bash
# 1. Edit code and docstrings
# 2. Regenerate
python scripts/generate_api_reference.py
# 3. Commit both
git add hic_basic/ docs/
git commit -m "Update code and docs"
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

- **85** modules documented
- **600+** public functions
- **~12,000** lines of documentation
- **Auto-generated** from source

## 🎓 Learning Path

1. **Start**: Read [docs/README.md](docs/README.md)
2. **Explore**: Run `python docs/EXAMPLES.py`
3. **Search**: Try `python docs/search_api.py --interactive`
4. **Reference**: Browse [docs/API_REFERENCE.md](docs/API_REFERENCE.md)
5. **Use**: Import and call functions!

---

**Questions?** See [docs/SUMMARY.md](docs/SUMMARY.md) for complete details.
