"""
Quick examples for using the hic_basic API documentation.

This file demonstrates how to use the generated API reference
for both manual lookup and programmatic access.
"""

import json
from pathlib import Path

# ========================================
# Example 1: Load and explore API index
# ========================================

# Load the API index
docs_dir = Path(__file__).parent
with open(docs_dir / 'api_index.json') as f:
    api_index = json.load(f)

print(f"Total modules: {len(api_index)}")
print(f"Total functions: {sum(len(funcs) for funcs in api_index.values())}")

# ========================================
# Example 2: Find all functions in a module
# ========================================

module_name = 'coolstuff'
if module_name in api_index:
    functions = api_index[module_name]
    print(f"\n{module_name} has {len(functions)} public functions:")
    for func in functions[:5]:  # Show first 5
        print(f"  - {func['name']}: {func['signature']}")

# ========================================
# Example 3: Search for functions by name pattern
# ========================================

search_term = 'plot'
matching_functions = []

for module, functions in api_index.items():
    for func in functions:
        if search_term.lower() in func['name'].lower():
            matching_functions.append((module, func))

print(f"\n Found {len(matching_functions)} functions with '{search_term}' in name:")
for module, func in matching_functions[:10]:  # Show first 10
    print(f"  - {module}.{func['name']}")

# ========================================
# Example 4: Find functions by keyword in docstring
# ========================================

keyword = 'convert'
matching_by_doc = []

for module, functions in api_index.items():
    for func in functions:
        if keyword.lower() in func['docstring'].lower():
            matching_by_doc.append((module, func))

print(f"\nFound {len(matching_by_doc)} functions mentioning '{keyword}':")
for module, func in matching_by_doc[:5]:  # Show first 5
    print(f"  - {module}.{func['name']}")
    # Print first line of docstring
    first_line = func['docstring'].split('\n')[0]
    print(f"    {first_line[:80]}...")

# ========================================
# Example 5: Get function details
# ========================================

def get_function_details(module_name, function_name):
    """Retrieve full details of a specific function."""
    if module_name not in api_index:
        return None
    
    for func in api_index[module_name]:
        if func['name'] == function_name:
            return func
    return None

# Example: Get details for cool2mat
func_details = get_function_details('coolstuff', 'cool2mat')
if func_details:
    print(f"\n\nFunction Details: coolstuff.{func_details['name']}")
    print(f"Signature: {func_details['signature']}")
    print(f"Return Type: {func_details['return_type']}")
    print(f"Line Number: {func_details['lineno']}")
    print(f"\nDocumentation:\n{func_details['docstring'][:300]}...")

# ========================================
# Example 6: Generate import statements
# ========================================

def generate_import(module_name, function_name):
    """Generate Python import statement."""
    # Convert module path to import
    import_path = f"hic_basic.{module_name.replace('/', '.')}"
    return f"from {import_path} import {function_name}"

print("\n\nGenerate imports for commonly used functions:")
common_functions = [
    ('coolstuff', 'cool2mat'),
    ('coolstuff', 'pairs2cool'),
    ('plot.hic', 'plot_cool'),
    ('hicio', 'read_meta'),
]

for module, func_name in common_functions:
    print(generate_import(module, func_name))

# ========================================
# Example 7: Find related functions
# ========================================

def find_related_functions(module_name, limit=10):
    """Find functions in the same module (related by context)."""
    if module_name not in api_index:
        return []
    return api_index[module_name][:limit]

print("\n\nRelated functions in 'plot.hic':")
related = find_related_functions('plot.hic', limit=5)
for func in related:
    print(f"  - {func['name']}: {func['signature']}")

# ========================================
# Example 8: Statistics by module
# ========================================

print("\n\nTop 10 modules by function count:")
module_stats = [(mod, len(funcs)) for mod, funcs in api_index.items()]
module_stats.sort(key=lambda x: x[1], reverse=True)

for i, (module, count) in enumerate(module_stats[:10], 1):
    print(f"  {i:2d}. {module:30s} - {count:3d} functions")

# ========================================
# Example 9: Export subset for specific purpose
# ========================================

# Export only plotting functions
plotting_functions = {}
for module, functions in api_index.items():
    if 'plot' in module:
        plotting_functions[module] = functions

print(f"\n\nPlotting modules: {len(plotting_functions)}")
print(f"Total plotting functions: {sum(len(f) for f in plotting_functions.values())}")

# Save subset
# with open(docs_dir / 'plotting_api.json', 'w') as f:
#     json.dump(plotting_functions, f, indent=2)

print("\n✓ Examples completed!")
print("\nFor more examples, see:")
print("  - docs/README.md")
print("  - docs/search_api.py --help")
