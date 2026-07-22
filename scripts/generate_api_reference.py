#!/usr/bin/env python3
"""
Generate API reference documentation for hic_basic.

This script scans Python modules in hic_basic and extracts public functions
(those without a leading underscore) into both human-readable and machine-
readable references.
"""
import ast
import json
from pathlib import Path
from typing import Dict, List


def first_doc_line(docstring: str) -> str:
    """Return the first non-empty line from a docstring."""
    for line in (docstring or "").splitlines():
        line = line.strip()
        if line:
            return line
    return ""


def module_name_from_path(rel_path: Path) -> str:
    """Convert a path under hic_basic/ to a dotted module name."""
    return str(rel_path).replace("/", ".").replace(".py", "")


def source_module_path(rel_path: Path) -> str:
    """Return a repo-relative source path for an indexed module."""
    return f"hic_basic/{str(rel_path)}"


def extract_functions_from_file(filepath: Path, module_name: str, module_path: str) -> List[Dict]:
    """Extract public functions from a Python file using AST parsing."""
    functions = []

    try:
        with open(filepath, "r", encoding="utf-8") as f:
            content = f.read()

        tree = ast.parse(content, filename=str(filepath))

        # Only module-level functions are importable using the generated import
        # statement. Class methods and nested helpers need their owning object.
        for node in tree.body:
            if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            if node.name.startswith("_"):
                continue

            def format_arg(arg):
                arg_str = arg.arg
                if arg.annotation:
                    try:
                        arg_str += f": {ast.unparse(arg.annotation)}"
                    except Exception:
                        pass
                return arg_str

            positional_nodes = list(node.args.posonlyargs) + list(node.args.args)
            positional = [format_arg(arg) for arg in positional_nodes]
            defaults = list(node.args.defaults)
            for i, default in enumerate(defaults):
                arg_idx = len(positional) - len(defaults) + i
                try:
                    positional[arg_idx] += f"={ast.unparse(default)}"
                except Exception:
                    positional[arg_idx] += "=..."

            args = []
            posonly_count = len(node.args.posonlyargs)
            for index, arg_str in enumerate(positional):
                args.append(arg_str)
                if posonly_count and index + 1 == posonly_count:
                    args.append("/")

            if node.args.vararg:
                args.append(f"*{format_arg(node.args.vararg)}")
            elif node.args.kwonlyargs:
                args.append("*")

            for arg, default in zip(node.args.kwonlyargs, node.args.kw_defaults):
                arg_str = format_arg(arg)
                if default is not None:
                    try:
                        arg_str += f"={ast.unparse(default)}"
                    except Exception:
                        arg_str += "=..."
                args.append(arg_str)

            if node.args.kwarg:
                args.append(f"**{format_arg(node.args.kwarg)}")

            signature = f"{node.name}({', '.join(args)})"
            docstring = ast.get_docstring(node) or "No documentation available."

            return_type = None
            if node.returns:
                try:
                    return_type = ast.unparse(node.returns)
                except Exception:
                    pass

            functions.append(
                {
                    "name": node.name,
                    "signature": signature,
                    "summary": first_doc_line(docstring),
                    "docstring": docstring,
                    "return_type": return_type,
                    "lineno": node.lineno,
                    "symbol": f"hic_basic.{module_name}.{node.name}",
                    "import": f"from hic_basic.{module_name} import {node.name}",
                    "source": {
                        "module_path": module_path,
                        "lineno": node.lineno,
                    },
                }
            )

    except Exception as e:
        print(f"Error parsing {filepath}: {e}")

    return functions


def collect_module_functions(base_path: Path) -> Dict[str, Dict]:
    """Collect public functions from all importable modules under base_path."""
    module_functions = {}

    for py_file in sorted(base_path.rglob("*.py")):
        if "__pycache__" in str(py_file) or ".egg-info" in str(py_file):
            continue

        rel_path = py_file.relative_to(base_path)
        if py_file.name in {"__init__.py", "__main__.py"}:
            continue

        module_name = module_name_from_path(rel_path)
        module_path = source_module_path(rel_path)
        functions = extract_functions_from_file(py_file, module_name, module_path)

        if functions:
            module_functions[module_name] = {
                "file": py_file,
                "module_path": module_path,
                "functions": functions,
            }

    return module_functions


def generate_markdown_docs(module_functions: Dict[str, Dict], output_path: Path):
    """Generate markdown documentation for all modules."""
    lines = []
    lines.append("# hic_basic API Reference\n")
    lines.append("*Auto-generated API documentation for public functions*\n")
    lines.append(f"Total modules with public functions: {len(module_functions)}\n")
    lines.append("---\n")

    lines.append("## Table of Contents\n")
    for module_name in sorted(module_functions.keys()):
        anchor = module_name.replace(".", "-").lower()
        func_count = len(module_functions[module_name]["functions"])
        lines.append(f"- [{module_name}](#{anchor}) ({func_count} functions)")
    lines.append("\n---\n")

    for module_name in sorted(module_functions.keys()):
        data = module_functions[module_name]
        module_path = data["module_path"].replace("hic_basic/", "")

        lines.append(f"\n## {module_name}\n")
        lines.append(f"**File:** `{module_path}`\n")
        lines.append(f"**Public functions:** {len(data['functions'])}\n")

        for func in sorted(data["functions"], key=lambda x: x["name"]):
            lines.append(f"\n### `{func['signature']}`\n")

            if func["return_type"]:
                lines.append(f"**Returns:** `{func['return_type']}`\n")

            docstring = func["docstring"].strip()
            if docstring:
                lines.append(f"\n{docstring}\n")

            lineno = func["source"]["lineno"]
            lines.append(f"\n**Source:** Line {lineno} in [{module_path}]({module_path}#L{lineno})\n")
            lines.append("\n---\n")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        content = "\n".join(lines)
        f.write("\n".join(line.rstrip() for line in content.splitlines()))

    print(f"Generated API reference: {output_path}")
    print(f"Total modules documented: {len(module_functions)}")
    total_functions = sum(len(data["functions"]) for data in module_functions.values())
    print(f"Total public functions: {total_functions}")


def generate_json_index(module_functions: Dict[str, Dict], output_path: Path):
    """Generate JSON index for programmatic access."""
    index = {
        module_name: sorted(data["functions"], key=lambda x: x["name"])
        for module_name, data in sorted(module_functions.items())
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(index, f, indent=2)

    print(f"Generated JSON index: {output_path}")


if __name__ == "__main__":
    script_dir = Path(__file__).parent
    repo_root = script_dir.parent
    hic_basic_path = repo_root / "hic_basic"
    docs_dir = repo_root / "docs"
    package_data_dir = hic_basic_path / "data"

    module_functions = collect_module_functions(hic_basic_path)

    generate_markdown_docs(module_functions, docs_dir / "API_REFERENCE.md")
    generate_json_index(module_functions, docs_dir / "api_index.json")
    generate_json_index(module_functions, package_data_dir / "api_index.json")

    print("\nDocumentation generation complete.")
    print(f"  - Markdown: {docs_dir / 'API_REFERENCE.md'}")
    print(f"  - JSON: {docs_dir / 'api_index.json'}")
    print(f"  - Package JSON: {package_data_dir / 'api_index.json'}")
