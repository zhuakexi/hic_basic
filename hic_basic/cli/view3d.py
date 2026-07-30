"""CLI for portable, browser-based 3D genome visualization."""

from __future__ import annotations

import json
from pathlib import Path

from ..visualize import structure_scene, structure_scene_with_feature, view_structure
from ..visualize.html import EXPORT_PROFILES


name = "view3d"
description = "Write a standalone interactive HTML viewer for a 3dg structure."


def add_arguments(subparser):
    subparser.add_argument("input", type=Path, help="Input five-column .3dg file")
    subparser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Output HTML path (default: input name with .html suffix)",
    )
    subparser.add_argument(
        "--color",
        default="chrom",
        help="Territory colors ('chrom') or a numeric/categorical input field",
    )
    subparser.add_argument(
        "--feature",
        type=Path,
        help="Headerless chrom/pos/value TSV used for per-particle and polymer-bond colors",
    )
    subparser.add_argument(
        "--feature-name",
        default="feature",
        help="Scientific feature label stored with the continuous color mapping",
    )
    subparser.add_argument(
        "--no-duplicate-haplotypes",
        action="store_true",
        help="Require feature chromosome names to match haplotype suffixes exactly",
    )
    subparser.add_argument("--cmap", default="viridis", help="Matplotlib colormap for numeric fields")
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
    subparser.add_argument("--title", help="Viewer title")
    subparser.add_argument(
        "--region",
        action="append",
        default=[],
        help="Initial chromosome/range selection; repeat for a union",
    )
    subparser.add_argument("--background", default="#151918", help="Initial monitor background color")
    subparser.add_argument("--particle-size", type=float, default=1.0, help="Initial particle size multiplier")
    subparser.add_argument("--bond-size", type=float, default=1.0, help="Initial bond width in pixels")
    subparser.add_argument(
        "--profile",
        choices=EXPORT_PROFILES,
        default="molecular-studio-v1",
        help="Default high-quality PNG profile",
    )
    subparser.add_argument("--force", action="store_true", help="Overwrite an existing output")
    subparser.add_argument("--json", action="store_true", help="Emit machine-readable result metadata")


def run(args):
    output = args.output or args.input.with_suffix(".html")
    if args.feature is not None:
        if args.color != "chrom":
            raise ValueError("--color and --feature are mutually exclusive")
        scene = structure_scene_with_feature(
            args.input,
            args.feature,
            feature_name=args.feature_name,
            duplicate_haplotypes=not args.no_duplicate_haplotypes,
            cmap=args.cmap,
            vmin=args.vmin,
            vmax=args.vmax,
            missing_color=args.missing_color,
            polymer_bonds=not args.no_bonds,
            max_bond_gap=args.max_bond_gap,
        )
    else:
        scene = structure_scene(
            args.input,
            color=args.color,
            cmap=args.cmap,
            vmin=args.vmin,
            vmax=args.vmax,
            missing_color=args.missing_color,
            polymer_bonds=not args.no_bonds,
            max_bond_gap=args.max_bond_gap,
        )
    viewer = view_structure(
        scene,
        title=args.title,
        regions=args.region,
        background=args.background,
        particle_size=args.particle_size,
        bond_size=args.bond_size,
        export_profile=args.profile,
    )
    result = viewer.write_html(output, overwrite=args.force)
    payload = {
        "output": str(result),
        "size_bytes": result.stat().st_size,
        "particles": scene.particle_count,
        "bonds": scene.bond_count,
        "chromosomes": len(scene.chromosomes),
        "resolution": scene.resolution,
        "max_bond_gap": scene.bonds.max_genomic_gap,
        "scene_sha256": scene.checksum,
    }
    if args.json:
        print(json.dumps(payload, indent=2))
    else:
        print(
            f"Wrote {result} ({scene.particle_count:,} particles, "
            f"{scene.bond_count:,} polymer bonds, {result.stat().st_size:,} bytes)"
        )
