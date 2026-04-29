"""CLI entry point for querying the generated hic_basic API index."""

from __future__ import annotations

import json
from typing import Any, Dict, List

from ..api_query import (
    entry_to_json,
    list_module_functions,
    list_modules,
    load_api_index,
    resolve_source,
    resolve_symbol,
    search_entries,
)


name = "api"
description = "Search and inspect the hic_basic Python API."


def add_arguments(subparser):
    subparsers = subparser.add_subparsers(dest="api_command", help="API lookup commands")

    list_parser = subparsers.add_parser("list", help="List API modules or functions in a module")
    list_parser.add_argument("module", nargs="?", help="Module name or partial module name")
    list_parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")

    search_parser = subparsers.add_parser("search", help="Search API functions by task keywords")
    search_parser.add_argument("query", nargs="+", help="Search keywords")
    search_parser.add_argument("--limit", type=int, default=10, help="Maximum number of matches")
    search_parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")

    show_parser = subparsers.add_parser("show", help="Show details for an API symbol")
    show_parser.add_argument("symbol", help="Function name or fully qualified hic_basic symbol")
    detail = show_parser.add_mutually_exclusive_group()
    detail.add_argument("--brief", action="store_true", help="Show only symbol, import, signature, and summary")
    detail.add_argument("--signature", action="store_true", help="Show only import and signature")
    show_parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")

    where_parser = subparsers.add_parser("where", help="Show source location for an API symbol")
    where_parser.add_argument("symbol", help="Function name or fully qualified hic_basic symbol")
    where_parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")


def run(args):
    api_index = load_api_index()

    if args.api_command == "list":
        _run_list(api_index, args)
    elif args.api_command == "search":
        _run_search(api_index, args)
    elif args.api_command == "show":
        _run_show(api_index, args)
    elif args.api_command == "where":
        _run_where(api_index, args)
    else:
        raise SystemExit("Run `hic_basic api --help` for available API lookup commands.")


def _run_list(api_index: Dict[str, List[Dict[str, Any]]], args):
    if args.module:
        entries = list_module_functions(api_index, args.module)
        if args.json:
            print(json.dumps({"query": args.module, "functions": [entry_to_json(e) for e in entries]}, indent=2))
            return
        if not entries:
            print(f"No functions found for module query: {args.module}")
            return
        modules = sorted({entry["module"] for entry in entries})
        if len(modules) == 1:
            print(f"Module: {modules[0]}")
        else:
            print(f"Module query: {args.module}")
            print(f"Matched modules: {', '.join(modules)}")
        print(f"Functions: {len(entries)}\n")
        for entry in entries:
            _print_function_card(entry)
        return

    modules = list_modules(api_index)
    if args.json:
        print(json.dumps({"modules": modules}, indent=2))
        return
    print(f"Found {len(modules)} modules with public functions:\n")
    for module in modules:
        print(f"  {module['module']:30s} {module['function_count']:4d} functions")


def _run_search(api_index: Dict[str, List[Dict[str, Any]]], args):
    query = " ".join(args.query)
    matches = search_entries(api_index, query, limit=args.limit)
    if args.json:
        payload = {
            "query": query,
            "matches": [
                {"rank": i, **entry_to_json(entry), "show_command": f"hic_basic api show {entry['symbol']}"}
                for i, entry in enumerate(matches, 1)
            ],
        }
        print(json.dumps(payload, indent=2))
        return

    if not matches:
        print(f"No matches found for query: {query}")
        return

    print(f"Search: {query}\n")
    for i, entry in enumerate(matches, 1):
        print(f"{i}. {entry['symbol']}")
        print(f"   Score: {entry['score']}")
        print(f"   Signature: {entry.get('signature', '')}")
        if entry.get("summary"):
            print(f"   Summary: {entry['summary']}")
        print(f"   Show: hic_basic api show {entry['symbol']}")
        print()


def _run_show(api_index: Dict[str, List[Dict[str, Any]]], args):
    status, matches = resolve_symbol(api_index, args.symbol)
    entry = _require_single_match(args.symbol, status, matches, json_output=args.json)
    if entry is None:
        return

    if args.json:
        print(json.dumps(entry_to_json(entry), indent=2))
        return

    if args.signature:
        print(entry["import"])
        print(entry.get("signature", ""))
        return

    print(entry["symbol"])
    print()
    print("Import:")
    print(f"  {entry['import']}")
    print()
    print("Signature:")
    print(f"  {entry.get('signature', '')}")

    if entry.get("summary"):
        print()
        print("Summary:")
        print(f"  {entry['summary']}")

    if not args.brief:
        print()
        print("Docstring:")
        for line in entry.get("docstring", "").splitlines() or [""]:
            print(f"  {line}")
        source = entry.get("source", {})
        module_path = source.get("module_path", f"hic_basic/{entry['module'].replace('.', '/')}.py")
        lineno = source.get("lineno")
        print()
        print("Source:")
        print(f"  {module_path}:{lineno}")


def _run_where(api_index: Dict[str, List[Dict[str, Any]]], args):
    status, matches = resolve_symbol(api_index, args.symbol)
    entry = _require_single_match(args.symbol, status, matches, json_output=args.json)
    if entry is None:
        return

    source = resolve_source(entry)
    if args.json:
        payload = entry_to_json(entry)
        payload["source"] = source
        print(json.dumps(payload, indent=2))
        return

    print(entry["symbol"])
    print()
    print("Source:")
    path = source.get("path") or source.get("module_path")
    lineno = source.get("lineno")
    print(f"  {path}:{lineno}")
    if not source.get("resolved"):
        print("  Source path could not be resolved in the installed environment.")
    print()
    print("Import:")
    print(f"  {entry['import']}")
    if path and lineno:
        end = source.get("end_lineno") or int(lineno) + 80
        print()
        print("Open:")
        print(f"  sed -n '{lineno},{end}p' {path}")


def _require_single_match(symbol: str, status: str, matches: List[Dict[str, Any]], json_output: bool = False):
    if status == "ok":
        return matches[0]

    if json_output:
        print(
            json.dumps(
                {
                    "query": symbol,
                    "status": status,
                    "matches": [entry_to_json(entry) for entry in matches],
                },
                indent=2,
            )
        )
        return None

    if status == "missing":
        print(f"No API symbol found for: {symbol}")
        return None

    print(f"Ambiguous symbol: {symbol}\n")
    print("Matches:")
    for i, entry in enumerate(matches, 1):
        print(f"  {i}. {entry['symbol']}")
        if entry.get("summary"):
            print(f"     {entry['summary']}")
        print(f"     Show: hic_basic api show {entry['symbol']}")
    print("\nUse a fully qualified symbol.")
    return None


def _print_function_card(entry: Dict[str, Any]):
    print(entry["name"])
    print(f"  Signature: {entry.get('signature', '')}")
    if entry.get("summary"):
        print(f"  Summary: {entry['summary']}")
    print(f"  Show: hic_basic api show {entry['symbol']}")
    print()
