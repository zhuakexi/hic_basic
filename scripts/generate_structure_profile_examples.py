#!/usr/bin/env python3
"""Generate the structure-render profile examples used by the documentation."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

from PIL import Image

from hic_basic.visualize import (
    CameraConfig,
    ExportStyle,
    LayerConfig,
    RasterConfig,
    SectionConfig,
    StructureRenderConfig,
    ViewConfig,
    render_structure_png,
    structure_scene,
)


PROFILES = (
    "molecular-studio-v1",
    "section-pymol-v1",
    "section-studio-v1",
    "xu2025-v2",
    "xu2025-v1",
    "publication-v1",
)


def optimize_png(path: Path) -> None:
    """Remove the opaque alpha channel and apply maximum lossless compression."""
    temporary = path.with_name(f".{path.stem}.optimized.png")
    with Image.open(path) as image:
        image.convert("RGB").save(temporary, optimize=True, compress_level=9)
    os.replace(temporary, path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("structure", type=Path, help="Whole-genome 20 kb .3dg input")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("docs/assets/structure-profiles"),
        help="PNG destination (default: docs/assets/structure-profiles)",
    )
    parser.add_argument("--size", type=int, default=800, help="Square image size in pixels")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    scene = structure_scene(args.structure, color="chrom")
    for profile in PROFILES:
        output = args.output_dir / f"{profile}.png"
        config = StructureRenderConfig(
            view=ViewConfig(
                camera=CameraConfig(mode="fit", direction=(1.08, 0.76, 1.32), fov=34),
                section=SectionConfig(enabled=False),
            ),
            layers=LayerConfig(
                show_particles=False,
                show_bonds=True,
                show_annotations=False,
            ),
            style=ExportStyle(profile=profile),
            raster=RasterConfig(width=args.size, height=args.size),
        )
        result = render_structure_png(scene, output, config=config)
        optimize_png(output)
        print(f"{profile}: {result.width}x{result.height}, {output.stat().st_size:,} bytes")


if __name__ == "__main__":
    main()
