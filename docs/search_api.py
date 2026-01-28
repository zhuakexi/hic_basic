#!/usr/bin/env python3
"""
Quick search tool for hic_basic API reference.

Usage:
    python docs/search_api.py cool2mat              # Search by function name
    python docs/search_api.py --module coolstuff    # List all functions in a module
    python docs/search_api.py --keyword convert     # Search in docstrings
    python docs/search_api.py --interactive         # Interactive mode
"""
import argparse
import json
import re
import sys
from pathlib import Path
from typing import Dict, List


def load_api_index(json_path: Path) -> Dict:
    """Load the API index JSON file."""
    with open(json_path, 'r') as f:
        return json.load(f)


def search_by_name(api_index: Dict, query: str, case_sensitive: bool = False) -> List:
    """Search for functions by name."""
    results = []
    
    if not case_sensitive:
        query = query.lower()
    
    for module, functions in api_index.items():
        for func in functions:
            func_name = func['name'] if case_sensitive else func['name'].lower()
            if query in func_name:
                results.append({
                    'module': module,
                    'function': func
                })
    
    return results


def search_by_keyword(api_index: Dict, keyword: str, case_sensitive: bool = False) -> List:
    """Search for functions by keyword in docstrings."""
    results = []
    
    if not case_sensitive:
        keyword = keyword.lower()
    
    for module, functions in api_index.items():
        for func in functions:
            docstring = func['docstring'] if case_sensitive else func['docstring'].lower()
            if keyword in docstring:
                results.append({
                    'module': module,
                    'function': func
                })
    
    return results


def list_module_functions(api_index: Dict, module: str) -> List:
    """List all functions in a specific module."""
    # Try exact match first
    if module in api_index:
        return [{'module': module, 'function': f} for f in api_index[module]]
    
    # Try partial match
    results = []
    for mod_name, functions in api_index.items():
        if module in mod_name:
            results.extend([{'module': mod_name, 'function': f} for f in functions])
    
    return results


def print_results(results: List, verbose: bool = False):
    """Print search results in a readable format."""
    if not results:
        print("No results found.")
        return
    
    print(f"\n{'='*80}")
    print(f"Found {len(results)} result(s)")
    print(f"{'='*80}\n")
    
    for i, result in enumerate(results, 1):
        module = result['module']
        func = result['function']
        
        print(f"{i}. {module}.{func['name']}")
        print(f"   Signature: {func['signature']}")
        
        if func['return_type']:
            print(f"   Returns: {func['return_type']}")
        
        if verbose:
            print(f"   Documentation:")
            # Print first 3 lines of docstring
            doc_lines = func['docstring'].split('\n')[:3]
            for line in doc_lines:
                print(f"      {line}")
            if len(func['docstring'].split('\n')) > 3:
                print(f"      ...")
        else:
            # Print just the first line
            first_line = func['docstring'].split('\n')[0]
            if first_line.strip():
                print(f"   {first_line[:100]}{'...' if len(first_line) > 100 else ''}")
        
        print(f"   Source: Line {func['lineno']} in hic_basic/{module.replace('.', '/')}.py")
        print()


def interactive_mode(api_index: Dict):
    """Interactive search mode."""
    print("\n" + "="*80)
    print("hic_basic API Interactive Search")
    print("="*80)
    print("\nCommands:")
    print("  name <query>      - Search by function name")
    print("  keyword <query>   - Search in docstrings")
    print("  module <name>     - List module functions")
    print("  list              - List all modules")
    print("  quit              - Exit")
    print()
    
    while True:
        try:
            command = input("\n> ").strip()
            
            if not command:
                continue
            
            if command.lower() in ('quit', 'exit', 'q'):
                break
            
            parts = command.split(maxsplit=1)
            cmd = parts[0].lower()
            query = parts[1] if len(parts) > 1 else ""
            
            if cmd == 'list':
                modules = sorted(api_index.keys())
                print(f"\nFound {len(modules)} modules:")
                for i, mod in enumerate(modules, 1):
                    func_count = len(api_index[mod])
                    print(f"  {i:2d}. {mod:30s} ({func_count} functions)")
            
            elif cmd == 'name' and query:
                results = search_by_name(api_index, query)
                print_results(results, verbose=True)
            
            elif cmd == 'keyword' and query:
                results = search_by_keyword(api_index, query)
                print_results(results, verbose=False)
            
            elif cmd == 'module' and query:
                results = list_module_functions(api_index, query)
                print_results(results, verbose=False)
            
            else:
                print("Unknown command. Type 'quit' to exit.")
        
        except (KeyboardInterrupt, EOFError):
            print("\nExiting...")
            break


def main():
    parser = argparse.ArgumentParser(
        description='Search hic_basic API reference',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s cool2mat                    # Find functions with 'cool2mat' in name
  %(prog)s --module coolstuff          # List all functions in coolstuff
  %(prog)s --keyword "convert pairs"   # Search for keyword in docs
  %(prog)s --interactive               # Start interactive mode
  %(prog)s --list-modules              # List all available modules
        """
    )
    
    parser.add_argument('query', nargs='?', help='Search query (function name)')
    parser.add_argument('-m', '--module', help='List functions in specific module')
    parser.add_argument('-k', '--keyword', help='Search keyword in docstrings')
    parser.add_argument('-i', '--interactive', action='store_true', help='Interactive mode')
    parser.add_argument('-v', '--verbose', action='store_true', help='Show detailed output')
    parser.add_argument('-l', '--list-modules', action='store_true', help='List all modules')
    parser.add_argument('--case-sensitive', action='store_true', help='Case-sensitive search')
    
    args = parser.parse_args()
    
    # Locate the JSON index
    script_dir = Path(__file__).parent
    json_path = script_dir / 'api_index.json'
    
    if not json_path.exists():
        print(f"Error: {json_path} not found.")
        print("Run: python scripts/generate_api_reference.py")
        sys.exit(1)
    
    # Load API index
    api_index = load_api_index(json_path)
    
    # Interactive mode
    if args.interactive:
        interactive_mode(api_index)
        return
    
    # List modules
    if args.list_modules:
        modules = sorted(api_index.keys())
        print(f"\nFound {len(modules)} modules with public functions:\n")
        for i, mod in enumerate(modules, 1):
            func_count = len(api_index[mod])
            print(f"  {i:2d}. {mod:30s} ({func_count} functions)")
        print()
        return
    
    # Module listing
    if args.module:
        results = list_module_functions(api_index, args.module)
        print_results(results, verbose=args.verbose)
        return
    
    # Keyword search
    if args.keyword:
        results = search_by_keyword(api_index, args.keyword, args.case_sensitive)
        print_results(results, verbose=args.verbose)
        return
    
    # Name search (default)
    if args.query:
        results = search_by_name(api_index, args.query, args.case_sensitive)
        print_results(results, verbose=args.verbose)
        return
    
    # No arguments provided
    parser.print_help()


if __name__ == '__main__':
    main()
