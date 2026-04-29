# hic_basic Documentation

This directory contains auto-generated API reference documentation for the `hic_basic` library.

## 📚 Contents

- **[API_REFERENCE.md](API_REFERENCE.md)** - Complete API reference with all public functions
- **[api_index.json](api_index.json)** - Machine-readable JSON index for programmatic access
- **[search_api.py](search_api.py)** - Repo-local compatibility wrapper for API search

## 🚀 Quick Start

### View the full API reference
```bash
# Open in your browser or editor
cat docs/API_REFERENCE.md
```

### Search for specific functions
```bash
# Search installed/package API
hic_basic api search cool2mat

# Show a known function
hic_basic api show cool2mat

# List module functions
hic_basic api list coolstuff

# Locate source
hic_basic api where cool2mat

# Repo-local fallback
python docs/search_api.py cool2mat
```

## 📖 Usage for Coding Agents

This documentation is designed to be easily consumed by AI coding assistants:

1. **Quick lookup**: Use `hic_basic api search <keywords>` or check `api_index.json`
2. **Detailed reference**: Read `API_REFERENCE.md` for full function signatures and docstrings
3. **Search**: Use `hic_basic api ...` for installed-package lookup, or `search_api.py` as a repo-local fallback

### Example for AI agents:
```python
import json

# Load the API index
with open('docs/api_index.json') as f:
    api_index = json.load(f)

# Find all functions in coolstuff module
coolstuff_funcs = api_index['coolstuff']
print(f"Found {len(coolstuff_funcs)} functions in coolstuff")

# Search for specific capability
for module, functions in api_index.items():
    for func in functions:
        if 'plot' in func['name'].lower():
            print(f"{module}.{func['name']}: {func['signature']}")
```

## 🔄 Regenerating Documentation

The documentation is auto-generated from the source code. To regenerate:

```bash
./run_in_env.sh python scripts/generate_api_reference.py
```

This will:
1. Scan all Python files in `hic_basic/`
2. Extract public functions (without `_` prefix)
3. Parse signatures and docstrings
4. Generate Markdown, repo-local JSON, and package-installed JSON outputs

## 🔑 Key Modules

Popular modules with many public functions:

- **coolstuff** (32 functions) - Cooler file conversions and utilities
- **hicio** (31 functions) - I/O operations and file parsing
- **plot.render** (27 functions) - Rendering and visualization
- **phasing** (21 functions) - Cell cycle phasing analysis
- **plot.utils** (20 functions) - Plotting utilities
- **plot.hic** (18 functions) - Hi-C heatmap plotting
- **data** (18 functions) - Data structures and genome references
- **DI** (18 functions) - Directionality Index calculations

## 💡 Tips

1. **For users**: Start with the most commonly used modules (coolstuff, hicio, plot.hic)
2. **For developers**: Keep docstrings up-to-date; they feed into this documentation
3. **For automation**: Use the JSON index for building tools and workflows

## 🤖 Integration with Tools

The JSON format enables easy integration with:
- IDE autocomplete and hints
- Documentation generators (Sphinx, MkDocs)
- API testing frameworks
- Workflow automation tools

## 📝 Contributing

When adding new public functions:
1. Add clear docstrings with parameters, returns, and examples
2. Use type hints where appropriate
3. Regenerate docs with `./run_in_env.sh python scripts/generate_api_reference.py`
4. Commit both code and updated documentation

---

*Last updated: Auto-generated from source code*
