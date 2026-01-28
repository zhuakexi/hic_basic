#!/usr/bin/env python3
"""
Generate API reference documentation for hic_basic.

This script scans all Python modules in hic_basic and extracts public functions
(those without _ prefix) along with their signatures and docstrings.
"""
import ast
import inspect
import os
import re
from pathlib import Path
from typing import Dict, List, Tuple


def extract_functions_from_file(filepath: Path) -> List[Dict]:
    """Extract public functions from a Python file using AST parsing."""
    functions = []
    
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()
        
        tree = ast.parse(content, filename=str(filepath))
        
        for node in ast.walk(tree):
            if isinstance(node, ast.FunctionDef):
                # Skip private functions (starting with _)
                if node.name.startswith('_'):
                    continue
                
                # Extract function signature
                args = []
                for arg in node.args.args:
                    arg_str = arg.arg
                    if arg.annotation:
                        try:
                            arg_str += f": {ast.unparse(arg.annotation)}"
                        except:
                            pass
                    args.append(arg_str)
                
                # Handle defaults
                defaults = node.args.defaults
                num_defaults = len(defaults)
                if num_defaults > 0:
                    for i in range(num_defaults):
                        arg_idx = len(args) - num_defaults + i
                        try:
                            default_val = ast.unparse(defaults[i])
                            args[arg_idx] += f"={default_val}"
                        except:
                            args[arg_idx] += "=..."
                
                signature = f"{node.name}({', '.join(args)})"
                
                # Extract docstring
                docstring = ast.get_docstring(node) or "No documentation available."
                
                # Extract return type if available
                return_type = None
                if node.returns:
                    try:
                        return_type = ast.unparse(node.returns)
                    except:
                        pass
                
                functions.append({
                    'name': node.name,
                    'signature': signature,
                    'docstring': docstring,
                    'return_type': return_type,
                    'lineno': node.lineno
                })
    
    except Exception as e:
        print(f"Error parsing {filepath}: {e}")
    
    return functions


def generate_markdown_docs(base_path: Path, output_path: Path):
    """Generate markdown documentation for all modules."""
    
    # Collect all Python files
    module_functions = {}
    
    for py_file in sorted(base_path.rglob('*.py')):
        # Skip __pycache__ and similar
        if '__pycache__' in str(py_file) or '.egg-info' in str(py_file):
            continue
        
        # Get relative path from base
        rel_path = py_file.relative_to(base_path)
        
        # Skip __init__.py for now (can be added later)
        if py_file.name == '__init__.py':
            continue
        
        # Extract functions
        functions = extract_functions_from_file(py_file)
        
        if functions:
            module_functions[str(rel_path)] = {
                'file': py_file,
                'functions': functions
            }
    
    # Generate markdown
    lines = []
    lines.append("# hic_basic API Reference\n")
    lines.append("*Auto-generated API documentation for public functions*\n")
    lines.append(f"Total modules with public functions: {len(module_functions)}\n")
    lines.append("---\n")
    
    # Table of contents
    lines.append("## Table of Contents\n")
    for module_path in sorted(module_functions.keys()):
        module_name = module_path.replace('/', '.').replace('.py', '')
        anchor = module_path.replace('/', '-').replace('.py', '').lower()
        func_count = len(module_functions[module_path]['functions'])
        lines.append(f"- [{module_name}](#{anchor}) ({func_count} functions)")
    lines.append("\n---\n")
    
    # Detailed documentation
    for module_path in sorted(module_functions.keys()):
        data = module_functions[module_path]
        module_name = module_path.replace('/', '.').replace('.py', '')
        
        lines.append(f"\n## {module_name}\n")
        lines.append(f"**File:** `{module_path}`\n")
        lines.append(f"**Public functions:** {len(data['functions'])}\n")
        
        for func in sorted(data['functions'], key=lambda x: x['name']):
            lines.append(f"\n### `{func['signature']}`\n")
            
            if func['return_type']:
                lines.append(f"**Returns:** `{func['return_type']}`\n")
            
            # Format docstring
            docstring = func['docstring'].strip()
            if docstring:
                lines.append(f"\n{docstring}\n")
            
            lines.append(f"\n**Source:** Line {func['lineno']} in [{module_path}]({module_path}#L{func['lineno']})\n")
            lines.append("\n---\n")
    
    # Write to file
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))
    
    print(f"Generated API reference: {output_path}")
    print(f"Total modules documented: {len(module_functions)}")
    total_functions = sum(len(data['functions']) for data in module_functions.values())
    print(f"Total public functions: {total_functions}")


def generate_json_index(base_path: Path, output_path: Path):
    """Generate JSON index for programmatic access."""
    import json
    
    module_functions = {}
    
    for py_file in sorted(base_path.rglob('*.py')):
        if '__pycache__' in str(py_file) or '.egg-info' in str(py_file):
            continue
        
        rel_path = py_file.relative_to(base_path)
        
        if py_file.name == '__init__.py':
            continue
        
        functions = extract_functions_from_file(py_file)
        
        if functions:
            module_name = str(rel_path).replace('/', '.').replace('.py', '')
            module_functions[module_name] = [
                {
                    'name': f['name'],
                    'signature': f['signature'],
                    'docstring': f['docstring'][:200] + '...' if len(f['docstring']) > 200 else f['docstring'],
                    'return_type': f['return_type'],
                    'lineno': f['lineno']
                }
                for f in functions
            ]
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(module_functions, f, indent=2)
    
    print(f"Generated JSON index: {output_path}")


if __name__ == '__main__':
    # Get the repository root
    script_dir = Path(__file__).parent
    repo_root = script_dir.parent
    hic_basic_path = repo_root / 'hic_basic'
    
    # Generate markdown documentation
    docs_dir = repo_root / 'docs'
    generate_markdown_docs(hic_basic_path, docs_dir / 'API_REFERENCE.md')
    
    # Generate JSON index
    generate_json_index(hic_basic_path, docs_dir / 'api_index.json')
    
    print("\n✓ Documentation generation complete!")
    print(f"  - Markdown: {docs_dir / 'API_REFERENCE.md'}")
    print(f"  - JSON: {docs_dir / 'api_index.json'}")
