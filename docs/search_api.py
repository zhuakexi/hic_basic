#!/usr/bin/env python3
"""
Quick search tool for hic_basic API reference.

This repo-local entry point is kept for compatibility. The shared query logic
lives in hic_basic.api_query and is also used by `hic_basic api ...`.
"""
import argparse
import sys
from pathlib import Path
from typing import List


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from hic_basic.api_query import (  # noqa: E402
    list_module_functions,
    list_modules,
    load_api_index,
    search_by_keyword,
    search_by_name,
)


def print_results(results: List[dict], verbose: bool = False):
    """Print search results in the legacy readable format."""
    if not results:
        print("No results found.")
        return

    print(f"\n{'=' * 80}")
    print(f"Found {len(results)} result(s)")
    print(f"{'=' * 80}\n")

    for i, entry in enumerate(results, 1):
        print(f"{i}. {entry['module']}.{entry['name']}")
        print(f"   Signature: {entry['signature']}")

        if entry.get("return_type"):
            print(f"   Returns: {entry['return_type']}")

        if verbose:
            print("   Documentation:")
            doc_lines = entry.get("docstring", "").split("\n")[:3]
            for line in doc_lines:
                print(f"      {line}")
            if len(entry.get("docstring", "").split("\n")) > 3:
                print("      ...")
        else:
            first_line = entry.get("summary", "")
            if first_line.strip():
                suffix = "..." if len(first_line) > 100 else ""
                print(f"   {first_line[:100]}{suffix}")

        source = entry.get("source", {})
        module_path = source.get("module_path", f"hic_basic/{entry['module'].replace('.', '/')}.py")
        lineno = source.get("lineno", entry.get("lineno"))
        print(f"   Source: Line {lineno} in {module_path}")
        print()


def interactive_mode(api_index):
    """Interactive search mode."""
    print("\n" + "=" * 80)
    print("hic_basic API Interactive Search")
    print("=" * 80)
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
            if command.lower() in ("quit", "exit", "q"):
                break

            parts = command.split(maxsplit=1)
            cmd = parts[0].lower()
            query = parts[1] if len(parts) > 1 else ""

            if cmd == "list":
                modules = list_modules(api_index)
                print(f"\nFound {len(modules)} modules:")
                for i, module in enumerate(modules, 1):
                    print(f"  {i:2d}. {module['module']:30s} ({module['function_count']} functions)")
            elif cmd == "name" and query:
                print_results(search_by_name(api_index, query), verbose=True)
            elif cmd == "keyword" and query:
                print_results(search_by_keyword(api_index, query), verbose=False)
            elif cmd == "module" and query:
                print_results(list_module_functions(api_index, query), verbose=False)
            else:
                print("Unknown command. Type 'quit' to exit.")
        except (KeyboardInterrupt, EOFError):
            print("\nExiting...")
            break


def main():
    parser = argparse.ArgumentParser(
        description="Search hic_basic API reference",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s cool2mat                    # Find functions with 'cool2mat' in name
  %(prog)s --module coolstuff          # List all functions in coolstuff
  %(prog)s --keyword "convert pairs"   # Search for keyword in docs
  %(prog)s --interactive               # Start interactive mode
  %(prog)s --list-modules              # List all available modules
        """,
    )

    parser.add_argument("query", nargs="?", help="Search query (function name)")
    parser.add_argument("-m", "--module", help="List functions in specific module")
    parser.add_argument("-k", "--keyword", help="Search keyword in docstrings")
    parser.add_argument("-i", "--interactive", action="store_true", help="Interactive mode")
    parser.add_argument("-v", "--verbose", action="store_true", help="Show detailed output")
    parser.add_argument("-l", "--list-modules", action="store_true", help="List all modules")
    parser.add_argument("--case-sensitive", action="store_true", help="Case-sensitive search")

    args = parser.parse_args()
    json_path = Path(__file__).parent / "api_index.json"

    if not json_path.exists():
        print(f"Error: {json_path} not found.")
        print("Run: python scripts/generate_api_reference.py")
        sys.exit(1)

    api_index = load_api_index(json_path)

    if args.interactive:
        interactive_mode(api_index)
        return

    if args.list_modules:
        modules = list_modules(api_index)
        print(f"\nFound {len(modules)} modules with public functions:\n")
        for i, module in enumerate(modules, 1):
            print(f"  {i:2d}. {module['module']:30s} ({module['function_count']} functions)")
        print()
        return

    if args.module:
        print_results(list_module_functions(api_index, args.module), verbose=args.verbose)
        return

    if args.keyword:
        print_results(search_by_keyword(api_index, args.keyword, args.case_sensitive), verbose=args.verbose)
        return

    if args.query:
        print_results(search_by_name(api_index, args.query, args.case_sensitive), verbose=args.verbose)
        return

    parser.print_help()


if __name__ == "__main__":
    main()
