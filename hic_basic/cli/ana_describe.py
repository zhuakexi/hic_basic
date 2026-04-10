"""CLI subcommand: `hic-basic ana-describe <ana_home>`

Read-only introspection of an Ana project directory.
Designed for coding agents to discover available columns and objects
without writing Python code.
"""

name = "ana-describe"
description = "Describe columns and objects in an Ana project (read-only)."


def add_arguments(subparser):
    subparser.add_argument(
        "home",
        type=str,
        help="Path to the Ana project directory.",
    )
    subparser.add_argument(
        "--column", "-c",
        type=str,
        default=None,
        help="Show full metadata for a specific column or object key.",
    )
    subparser.add_argument(
        "--file-format", "-f",
        type=str,
        default=None,
        help="Look up a file format definition (e.g., 'cool', 'bam').",
    )


def run(args):
    from ..ana import Ana

    ana = Ana(args.home, verbose=False)

    if args.column:
        meta = ana.column_meta(args.column)
        if meta:
            import json
            print(json.dumps(meta, indent=2, ensure_ascii=False))
        else:
            print(f"No metadata found for '{args.column}'.")
            print("Available columns:")
            ana.describe()
    elif args.file_format:
        fmt = ana.file_format(args.file_format)
        if fmt:
            import json
            print(json.dumps(fmt, indent=2, ensure_ascii=False))
        else:
            print(f"No file format definition found for '{args.file_format}'.")
    else:
        ana.describe()
