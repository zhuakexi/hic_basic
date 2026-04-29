"""Utilities for querying the generated hic_basic API index."""

from __future__ import annotations

import importlib
import inspect
import json
from importlib.resources import files
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple


Index = Dict[str, List[Dict[str, Any]]]


def load_api_index(json_path: Optional[Path] = None) -> Index:
    """Load an API index from a path or from package resources."""
    if json_path is None:
        json_path = files("hic_basic").joinpath("data/api_index.json")
    with json_path.open("r", encoding="utf-8") as f:
        return json.load(f)


def iter_entries(api_index: Index) -> Iterable[Dict[str, Any]]:
    """Yield normalized function entries from an API index."""
    for module, functions in api_index.items():
        for func in functions:
            yield normalize_entry(module, func)


def normalize_entry(module: str, func: Dict[str, Any]) -> Dict[str, Any]:
    """Return a function entry with stable derived fields."""
    entry = dict(func)
    entry.setdefault("module", module)
    entry.setdefault("symbol", f"hic_basic.{module}.{entry['name']}")
    entry.setdefault("import", f"from hic_basic.{module} import {entry['name']}")
    entry.setdefault("summary", first_doc_line(entry.get("docstring", "")))
    source = entry.get("source") or {}
    if not source:
        source = {
            "module_path": f"hic_basic/{module.replace('.', '/')}.py",
            "lineno": entry.get("lineno"),
        }
    entry["source"] = source
    return entry


def first_doc_line(docstring: str) -> str:
    """Extract the first non-empty line from a docstring."""
    for line in (docstring or "").splitlines():
        line = line.strip()
        if line:
            return line
    return ""


def list_modules(api_index: Index) -> List[Dict[str, Any]]:
    """Return module names and public function counts."""
    return [
        {"module": module, "function_count": len(functions)}
        for module, functions in sorted(api_index.items())
    ]


def list_module_functions(api_index: Index, module_query: str) -> List[Dict[str, Any]]:
    """List functions for exact or partial module matches."""
    if module_query in api_index:
        return [normalize_entry(module_query, f) for f in api_index[module_query]]

    results: List[Dict[str, Any]] = []
    query = module_query.lower()
    for module, functions in api_index.items():
        if query in module.lower():
            results.extend(normalize_entry(module, f) for f in functions)
    return results


def search_by_name(api_index: Index, query: str, case_sensitive: bool = False) -> List[Dict[str, Any]]:
    """Search for functions by name."""
    results = []
    needle = query if case_sensitive else query.lower()
    for entry in iter_entries(api_index):
        name = entry["name"] if case_sensitive else entry["name"].lower()
        if needle in name:
            results.append(entry)
    return results


def search_by_keyword(api_index: Index, keyword: str, case_sensitive: bool = False) -> List[Dict[str, Any]]:
    """Search for functions by keyword in docstrings."""
    results = []
    needle = keyword if case_sensitive else keyword.lower()
    for entry in iter_entries(api_index):
        docstring = entry.get("docstring", "")
        haystack = docstring if case_sensitive else docstring.lower()
        if needle in haystack:
            results.append(entry)
    return results


def search_entries(api_index: Index, query: str, limit: int = 10) -> List[Dict[str, Any]]:
    """Rank API entries using deterministic keyword matching."""
    tokens = [token.lower() for token in query.split() if token.strip()]
    if not tokens:
        return []

    matches = []
    for entry in iter_entries(api_index):
        score = score_entry(entry, tokens)
        if score > 0:
            ranked = dict(entry)
            ranked["score"] = score
            matches.append(ranked)

    matches.sort(key=lambda item: (-item["score"], item["symbol"]))
    return matches[:limit]


def score_entry(entry: Dict[str, Any], tokens: List[str]) -> int:
    """Score an entry for query tokens."""
    name = entry["name"].lower()
    module = entry["module"].lower()
    summary = entry.get("summary", "").lower()
    docstring = entry.get("docstring", "").lower()

    score = 0
    query = " ".join(tokens)
    if name == query:
        score += 100

    for token in tokens:
        if token == name:
            score += 100
        elif token in name:
            score += 50
        if token in module:
            score += 30
        if token in summary:
            score += 20
        if token in docstring:
            score += 10
    return score


def resolve_symbol(api_index: Index, symbol: str) -> Tuple[str, List[Dict[str, Any]]]:
    """Resolve a short function name or fully qualified hic_basic symbol."""
    query = symbol.strip()
    if query.startswith("hic_basic."):
        query = query[len("hic_basic.") :]

    exact_matches = []
    partial_matches = []
    for entry in iter_entries(api_index):
        module_symbol = f"{entry['module']}.{entry['name']}"
        if query == entry["name"] or query == module_symbol:
            exact_matches.append(entry)
        elif query in entry["name"] or query in module_symbol:
            partial_matches.append(entry)

    if exact_matches:
        return ("ok" if len(exact_matches) == 1 else "ambiguous", exact_matches)
    if partial_matches:
        return ("ok" if len(partial_matches) == 1 else "ambiguous", partial_matches)
    return "missing", []


def resolve_source(entry: Dict[str, Any]) -> Dict[str, Any]:
    """Resolve source location in the installed environment when possible."""
    source = dict(entry.get("source") or {})
    module_name = f"hic_basic.{entry['module']}"
    try:
        module = importlib.import_module(module_name)
        obj = getattr(module, entry["name"])
        path = inspect.getsourcefile(obj)
        lines, lineno = inspect.getsourcelines(obj)
        source.update(
            {
                "path": path,
                "lineno": lineno,
                "end_lineno": lineno + len(lines) - 1,
                "resolved": True,
            }
        )
    except Exception as exc:  # pragma: no cover - depends on optional imports
        source.setdefault("path", source.get("module_path"))
        source.setdefault("lineno", source.get("lineno"))
        source["resolved"] = False
        source["error"] = str(exc)
    return source


def entry_to_json(entry: Dict[str, Any], include_resolved_source: bool = False) -> Dict[str, Any]:
    """Return a stable JSON-serializable entry payload."""
    payload = {
        "symbol": entry["symbol"],
        "module": entry["module"],
        "name": entry["name"],
        "import": entry["import"],
        "signature": entry.get("signature"),
        "summary": entry.get("summary", ""),
        "docstring": entry.get("docstring", ""),
        "return_type": entry.get("return_type"),
        "source": entry.get("source", {}),
    }
    if "score" in entry:
        payload["score"] = entry["score"]
    if include_resolved_source:
        payload["source"] = resolve_source(entry)
    return payload
