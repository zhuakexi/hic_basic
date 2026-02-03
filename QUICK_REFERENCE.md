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

- **76** modules documented
- **528** public functions  
- **~10,000** lines of documentation
- **Auto-generated** from source

## 🎓 Learning Path

1. **Start**: Read [docs/README.md](docs/README.md)
2. **Explore**: Run `python docs/EXAMPLES.py`
3. **Search**: Try `python docs/search_api.py --interactive`
4. **Reference**: Browse [docs/API_REFERENCE.md](docs/API_REFERENCE.md)
5. **Use**: Import and call functions!

---

**Questions?** See [docs/SUMMARY.md](docs/SUMMARY.md) for complete details.
