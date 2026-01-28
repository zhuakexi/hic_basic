# API Documentation System Summary

## 📋 Overview

This directory contains a **complete, auto-generated API reference** for the hic_basic library, designed to be used by:
- 👨‍💻 Human developers
- 🤖 AI coding assistants  
- 🔧 Automation tools
- 📚 Documentation generators

## 🎯 Goals Achieved

1. ✅ **Convenient reference** - Quickly check what functions exist and how to use them
2. ✅ **AI-friendly** - Structured data that coding agents can easily parse
3. ✅ **Automation-ready** - JSON format for building workflows and tools
4. ✅ **Maintainable** - Auto-generated from source code docstrings

## 📁 Files

| File | Purpose | Best For |
|------|---------|----------|
| **API_REFERENCE.md** | Human-readable reference (~10k lines) | Reading in editor/browser |
| **api_index.json** | Machine-readable index | Programmatic access |
| **search_api.py** | Command-line search tool | Quick lookups |
| **EXAMPLES.py** | Usage examples | Learning the system |
| **README.md** | This guide | Understanding the system |

## 🔢 Statistics

- **76 modules** documented
- **527 public functions** indexed
- **~10,000 lines** of documentation
- **Auto-generated** from source code

## 🚀 Quick Start

### For Humans

```bash
# Search for a function
python docs/search_api.py cool2mat

# List module functions
python docs/search_api.py --module coolstuff

# Interactive mode
python docs/search_api.py --interactive

# Read full reference
less docs/API_REFERENCE.md
```

### For AI Agents

```python
import json

# Load API index
with open('docs/api_index.json') as f:
    api = json.load(f)

# Find what you need
coolstuff_funcs = api['coolstuff']
print(f"coolstuff has {len(coolstuff_funcs)} functions")

# Search by name
for module, funcs in api.items():
    for func in funcs:
        if 'plot' in func['name']:
            print(f"{module}.{func['name']}")
```

### For Automation

```python
# Generate workflow from API
import json

with open('docs/api_index.json') as f:
    api = json.load(f)

# Build a pipeline
pipeline = [
    ('hicio', 'read_meta'),
    ('coolstuff', 'pairs2cool'),
    ('plot.hic', 'plot_cool'),
]

for module, func_name in pipeline:
    func = next(f for f in api[module] if f['name'] == func_name)
    print(f"Step: {func['signature']}")
```

## 🔄 Updating Documentation

Whenever you add or modify public functions:

```bash
python scripts/generate_api_reference.py
```

This will:
1. Scan all `*.py` files in `hic_basic/`
2. Extract public functions (no `_` prefix)
3. Parse signatures and docstrings using AST
4. Generate both Markdown and JSON outputs

**Commit both code and documentation together!**

## 🏗️ Architecture

```
hic_basic/
├── hic_basic/           # Source code (76 modules)
│   ├── coolstuff.py     # Core functions
│   ├── hicio.py
│   ├── plot/
│   └── ...
├── scripts/
│   └── generate_api_reference.py  # Generator script
└── docs/
    ├── API_REFERENCE.md           # Human-readable
    ├── api_index.json             # Machine-readable
    ├── search_api.py              # Search tool
    ├── EXAMPLES.py                # Usage examples
    └── README.md                  # This file
```

## 🎯 Design Decisions

### Why both Markdown and JSON?

- **Markdown**: Easy to read, browse, and search in editors
- **JSON**: Easy to parse programmatically for tools and AI

### Why auto-generate?

- Ensures documentation stays in sync with code
- Reduces manual maintenance burden
- Encourages better docstrings

### Why include line numbers?

- Quick navigation to source code
- Easy verification of documentation
- Helpful for debugging

### What about private functions?

- Only public functions (no `_` prefix) are documented
- Private functions are implementation details
- Keeps API documentation focused and clean

## 📊 Top Modules by Function Count

1. **coolstuff** (32) - Cooler file operations
2. **hicio** (31) - I/O and file parsing  
3. **plot.render** (27) - Rendering utilities
4. **phasing** (21) - Cell cycle phasing
5. **plot.utils** (20) - Plotting helpers

[See full list in search_api.py --list-modules]

## 💡 Tips

### For Library Users

- Start with commonly used modules: `coolstuff`, `hicio`, `plot.hic`
- Use search tool to discover relevant functions
- Check docstrings for usage examples

### For Contributors

- Write clear docstrings with parameters and returns
- Use type hints where appropriate
- Run generator after adding new functions
- Commit documentation with code changes

### For AI Assistants

- Load `api_index.json` at session start
- Use it to suggest relevant functions
- Reference line numbers for context
- Generate correct import statements

## 🔗 Integration Ideas

This documentation system can be used with:

- **Sphinx** - Generate HTML documentation
- **MkDocs** - Create documentation website
- **IDE plugins** - Power autocomplete and hints
- **CI/CD** - Validate docstrings in tests
- **Workflow tools** - Auto-generate pipelines
- **API clients** - Generate client libraries

## 🤝 Contributing

When adding new functions:

1. Add clear docstrings following existing patterns
2. Include parameters, return types, and examples
3. Use type hints in signatures
4. Run `python scripts/generate_api_reference.py`
5. Commit both code and updated docs

## 📝 Docstring Best Practices

```python
def example_function(param1: str, param2: int = 10) -> pd.DataFrame:
    """
    Brief one-line description.
    
    More detailed explanation if needed. Can span
    multiple lines.
    
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
    >>> print(result.shape)
    (100, 5)
    """
    # implementation
```

## 🐛 Troubleshooting

**Documentation not updating?**
- Re-run the generator script
- Check for syntax errors in docstrings
- Verify file permissions

**Function not appearing?**
- Ensure it doesn't start with `_`
- Check if it's actually a function (not a class/variable)
- Verify the module is being scanned

**Search not working?**
- Ensure `api_index.json` exists
- Re-generate if it's out of date
- Check Python version compatibility

## 📚 Related Resources

- Main README: `../README.md`
- Contributing guide: (TBD)
- API Reference: `API_REFERENCE.md`
- Examples: `EXAMPLES.py`

---

**Last Updated**: Auto-generated from source code

**Questions?** Open an issue or check the examples!
