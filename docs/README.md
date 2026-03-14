# hic_basic Documentation

This directory contains auto-generated API reference documentation for the `hic_basic` library.

## 📚 Contents

- **[API_REFERENCE.md](API_REFERENCE.md)** - Complete API reference with all public functions (530 functions across 77 modules)
- **[api_index.json](api_index.json)** - Machine-readable JSON index for programmatic access
- **[search_api.py](search_api.py)** - Command-line tool to search the API

## 🚀 Quick Start

### View the full API reference
```bash
# Open in your browser or editor
cat docs/API_REFERENCE.md
```

### Search for specific functions
```bash
# Search by function name
python docs/search_api.py cool2mat

# Search by module
python docs/search_api.py --module coolstuff

# Search by keyword in docstrings
python docs/search_api.py --keyword "convert"

# Interactive search
python docs/search_api.py --interactive
```

## 📖 Usage for Coding Agents

This documentation is designed to be easily consumed by AI coding assistants:

1. **Quick lookup**: Check `api_index.json` for a list of all available functions
2. **Detailed reference**: Read `API_REFERENCE.md` for full function signatures and docstrings
3. **Search**: Use `search_api.py` to find relevant functions programmatically

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
python scripts/generate_api_reference.py
```

This will:
1. Scan all Python files in `hic_basic/`
2. Extract public functions (without `_` prefix)
3. Parse signatures and docstrings
4. Generate both Markdown and JSON outputs

## 📊 Statistics

- **Total modules**: 77
- **Total public functions**: 530
- **Documentation lines**: ~10,000

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
3. Regenerate docs with `python scripts/generate_api_reference.py`
4. Commit both code and updated documentation

---

*Last updated: Auto-generated from source code*
