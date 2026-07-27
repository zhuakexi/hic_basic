"""CLI for deterministic direct rendering of 3D genome structures."""

from __future__ import annotations

import json
import argparse
from pathlib import Path

from ..visualize import BrowserBackendConfig, render_structure_png, structure_scene


name = "render3d"
description = "Render a 3dg structure or standalone viewer HTML directly to PNG."


def _viewport(value: str) -> tuple[int, int]:
    try:
        width, height = value.lower().split("x", 1)
        return int(width), int(height)
    except (TypeError, ValueError) as error:
        raise argparse.ArgumentTypeError("viewport must use WIDTHxHEIGHT, for example 1440x900") from error


def add_arguments(subparser):
    subparser.add_argument("input", type=Path, help="Input five-column .3dg or standalone viewer HTML")
    subparser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Output PNG path (default: input name with .png suffix)",
    )
    subparser.add_argument(
        "--config",
        type=Path,
        help="Versioned JSON render config; controls regions, camera, layers, style, and size",
    )
    subparser.add_argument("--color", default="chrom", help="3dg territory colors or a numeric/categorical field")
    subparser.add_argument("--cmap", default="viridis", help="Matplotlib colormap for numeric 3dg fields")
    subparser.add_argument("--vmin", type=float, help="Lower continuous-color bound")
    subparser.add_argument("--vmax", type=float, help="Upper continuous-color bound")
    subparser.add_argument("--missing-color", default="#aeb7ba", help="Color for missing field values")
    subparser.add_argument(
        "--max-bond-gap",
        type=int,
        default=1_000_000,
        help="Maximum adjacent-bin genomic gap in bp (default: 1000000)",
    )
    subparser.add_argument("--no-bonds", action="store_true", help="Do not create polymer-chain bonds")
    subparser.add_argument("--browser-executable", type=Path, help="Explicit compatible Chromium executable")
    subparser.add_argument("--playwright-path", type=Path, help="Explicit Playwright Node module directory")
    subparser.add_argument("--timeout", type=float, default=180, help="Browser timeout in seconds")
    subparser.add_argument(
        "--viewport",
        type=_viewport,
        default=(1440, 900),
        help="Headless monitor viewport (default: 1440x900)",
    )
    subparser.add_argument(
        "--hardware-rendering",
        action="store_true",
        help="Allow the browser to select hardware WebGL instead of forcing SwiftShader",
    )
    subparser.add_argument("--force", action="store_true", help="Overwrite an existing output")
    subparser.add_argument("--json", action="store_true", help="Emit the complete machine-readable render report")


def run(args):
    output = args.output or args.input.with_suffix(".png")
    if args.input.suffix.lower() in {".html", ".htm"}:
        source = args.input
    else:
        source = structure_scene(
            args.input,
            color=args.color,
            cmap=args.cmap,
            vmin=args.vmin,
            vmax=args.vmax,
            missing_color=args.missing_color,
            polymer_bonds=not args.no_bonds,
            max_bond_gap=args.max_bond_gap,
        )
    backend = BrowserBackendConfig(
        executable_path=args.browser_executable,
        playwright_path=args.playwright_path,
        timeout_seconds=args.timeout,
        software_rendering=not args.hardware_rendering,
        viewport=args.viewport,
    )
    result = render_structure_png(
        source,
        output,
        config=args.config,
        backend=backend,
        overwrite=args.force,
    )
    report = result.to_dict()
    if args.json:
        print(json.dumps(report, indent=2, sort_keys=True))
    else:
        print(
            f"Wrote {result.output} ({result.width}x{result.height}, {result.profile}, "
            f"{result.size_bytes:,} bytes)"
        )
