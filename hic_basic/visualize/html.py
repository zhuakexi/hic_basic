"""Standalone HTML output for renderer-independent 3D genome scenes."""

from __future__ import annotations

import base64
import hashlib
import html
import json
import os
import tempfile
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence, Tuple, Union

import numpy as np
from matplotlib import colors as mpl_colors

from .config import (
    RENDER_CONFIG_SCHEMA_VERSION,
    EXPORT_PROFILES,
    ExportStyle,
    RenderConfigInput,
    StructureRenderConfig,
    ViewConfig,
    load_render_config,
)
from .scene import StructureScene
from .structure import ColorInput, StructureInput, structure_scene


RENDERER_PROTOCOL_VERSION = 2


@dataclass(frozen=True)
class StructureViewer:
    """Configuration for a portable browser view of a StructureScene."""

    scene: StructureScene
    title: str = "3D genome structure"
    regions: Tuple[str, ...] = ()
    background: str = "#151918"
    particle_size: float = 1.0
    bond_size: float = 1.0
    export_profile: str = "molecular-studio-v1"
    render_config: Optional[StructureRenderConfig] = None

    def __post_init__(self) -> None:
        if not isinstance(self.scene, StructureScene):
            raise TypeError("scene must be a StructureScene")
        normalized_background = mpl_colors.to_hex(self.background, keep_alpha=False)
        if not 0.5 <= float(self.particle_size) <= 3.0:
            raise ValueError("particle_size must be between 0.5 and 3.0")
        if not 1.0 <= float(self.bond_size) <= 6.0:
            raise ValueError("bond_size must be between 1.0 and 6.0")
        if self.export_profile not in EXPORT_PROFILES:
            choices = ", ".join(EXPORT_PROFILES)
            raise ValueError(f"export_profile must be one of: {choices}")
        if self.render_config is None:
            render_config = StructureRenderConfig(
                view=ViewConfig(regions=tuple(str(region) for region in self.regions)),
                style=ExportStyle(profile=self.export_profile),
            )
        else:
            render_config = load_render_config(self.render_config)
        object.__setattr__(self, "background", normalized_background)
        object.__setattr__(self, "particle_size", float(self.particle_size))
        object.__setattr__(self, "bond_size", float(self.bond_size))
        object.__setattr__(self, "regions", render_config.view.regions)
        object.__setattr__(self, "export_profile", render_config.style.profile)
        object.__setattr__(self, "render_config", render_config)

    def to_html(self) -> str:
        """Return a self-contained HTML document with no network dependencies."""
        payload = _scene_payload(self)
        bundle = _asset_text("viewer.min.js").replace("</script>", "<\\/script>")
        styles = _asset_text("viewer.css").replace("</style>", "<\\/style>")
        payload_json = json.dumps(payload, ensure_ascii=True, separators=(",", ":"))
        payload_json = payload_json.replace("</", "<\\/")
        title = html.escape(self.title, quote=False)
        return _html_document(
            title=title,
            payload_json=payload_json,
            bundle=bundle,
            styles=styles,
            background=self.background,
            particle_size=self.particle_size,
            bond_size=self.bond_size,
            export_profile=self.export_profile,
        )

    def write_html(self, output: Union[str, os.PathLike], overwrite: bool = True) -> Path:
        """Atomically write a standalone HTML artifact and return its path."""
        output_path = Path(output).expanduser()
        if output_path.exists() and not overwrite:
            raise FileExistsError(f"output already exists: {output_path}")
        if not output_path.parent.exists():
            output_path.parent.mkdir(parents=True)
        document = self.to_html()
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(
                mode="w",
                encoding="utf-8",
                prefix=f".{output_path.name}.",
                suffix=".tmp",
                dir=output_path.parent,
                delete=False,
            ) as stream:
                temporary = Path(stream.name)
                stream.write(document)
                stream.flush()
                os.fsync(stream.fileno())
            os.replace(temporary, output_path)
        finally:
            if temporary is not None and temporary.exists():
                temporary.unlink()
        return output_path


def view_structure(
    scene: StructureScene,
    *,
    title: Optional[str] = None,
    regions: Optional[Sequence[str]] = None,
    background: str = "#151918",
    particle_size: float = 1.0,
    bond_size: float = 1.0,
    export_profile: str = "molecular-studio-v1",
    render_config: RenderConfigInput = None,
) -> StructureViewer:
    """Create a configurable offline browser viewer for a structure scene.

    The returned object is lightweight and does not start a server or renderer.
    Call ``write_html`` to create one portable HTML file. ``regions`` accepts
    chromosome names or inclusive ranges such as ``chr1(mat):0-50Mb``; multiple
    entries are displayed as a union. When supplied, ``render_config`` is the
    authoritative region/profile request embedded for direct PNG rendering;
    the older ``regions`` and ``export_profile`` arguments remain compatible.
    """
    return StructureViewer(
        scene=scene,
        title=str(title or scene.dataset_name),
        regions=tuple(regions or ()),
        background=background,
        particle_size=particle_size,
        bond_size=bond_size,
        export_profile=export_profile,
        render_config=render_config,
    )


def write_structure_html(
    structure: StructureInput,
    output: Union[str, os.PathLike],
    *,
    color: ColorInput = "chrom",
    palette: Optional[Mapping[Any, Any]] = None,
    cmap: str = "viridis",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    missing_color: Any = "#aeb7ba",
    polymer_bonds: bool = True,
    max_bond_gap: int = 1_000_000,
    dataset_name: Optional[str] = None,
    metadata: Optional[Mapping[str, Any]] = None,
    title: Optional[str] = None,
    regions: Optional[Sequence[str]] = None,
    background: str = "#151918",
    particle_size: float = 1.0,
    bond_size: float = 1.0,
    export_profile: str = "molecular-studio-v1",
    render_config: RenderConfigInput = None,
    overwrite: bool = True,
) -> Path:
    """Build a 3D genome scene and write it as one portable interactive HTML.

    This convenience API performs no external process invocation. JavaScript,
    CSS, typed arrays, colors, and polymer bonds are embedded in the output, so
    the result can be downloaded and opened directly or served over Remote SSH.
    ``render_config`` embeds the same typed request accepted by
    ``render_structure_png``.

    Returns
    -------
    pathlib.Path
        Path to the atomically written HTML artifact.
    """
    scene = structure_scene(
        structure,
        color=color,
        palette=palette,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        missing_color=missing_color,
        polymer_bonds=polymer_bonds,
        max_bond_gap=max_bond_gap,
        dataset_name=dataset_name,
        metadata=metadata,
    )
    viewer = view_structure(
        scene,
        title=title,
        regions=regions,
        background=background,
        particle_size=particle_size,
        bond_size=bond_size,
        export_profile=export_profile,
        render_config=render_config,
    )
    return viewer.write_html(output, overwrite=overwrite)


def _scene_payload(viewer: StructureViewer) -> dict:
    scene = viewer.scene
    payload = {
        "schema_version": scene.schema_version,
        "renderer_protocol_version": RENDERER_PROTOCOL_VERSION,
        "render_config_schema_version": RENDER_CONFIG_SCHEMA_VERSION,
        "dataset": {
            "name": scene.dataset_name,
            "source": scene.source,
            "sha256": scene.checksum,
            "resolution": scene.resolution,
            "max_bond_gap": scene.bonds.max_genomic_gap,
            "metadata": _json_compatible(scene.metadata),
        },
        "stats": {
            "particles": scene.particle_count,
            "bonds": scene.bond_count,
            "chromosomes": len(scene.chromosomes),
            "bbox": scene.bbox,
        },
        "chromosomes": scene.chromosomes,
        "groups": [
            {
                "name": group.name,
                "code": group.code,
                "start": group.start,
                "count": group.count,
                "bond_start": group.bond_start,
                "bond_count": group.bond_count,
                "bbox": group.bbox,
                "color": group.color,
            }
            for group in scene.groups
        ],
        "viewer": {
            "base_point_size": scene.base_point_size,
            "background": viewer.background,
            "particle_size": viewer.particle_size,
            "bond_size": viewer.bond_size,
            "export_profile": viewer.export_profile,
            "initial_regions": viewer.regions,
            "render_config": viewer.render_config.to_dict(),
        },
        "arrays": {
            "positions": _encode_array(scene.particles.positions),
            "colors": _encode_array(scene.particles.colors),
            "chromosome_codes": _encode_array(scene.particles.chromosome_codes),
            "genomic_positions": _encode_array(scene.particles.genomic_positions),
            "bond_indices": _encode_array(scene.bonds.endpoints),
        },
        "generated_at_ms": 0,
    }
    json.dumps(payload, ensure_ascii=True)
    return payload


def _encode_array(array: np.ndarray) -> dict:
    contiguous = np.ascontiguousarray(array)
    raw = contiguous.tobytes()
    return {
        "data": base64.b64encode(raw).decode("ascii"),
        "dtype": _dtype_name(contiguous.dtype),
        "shape": list(contiguous.shape),
        "byte_order": "little",
        "sha256": hashlib.sha256(raw).hexdigest(),
    }


def _dtype_name(dtype: np.dtype) -> str:
    dtype = np.dtype(dtype)
    if dtype == np.dtype("float32"):
        return "float32"
    if dtype == np.dtype("uint8"):
        return "uint8"
    if dtype == np.dtype("uint16"):
        return "uint16"
    if dtype == np.dtype("uint32"):
        return "uint32"
    raise TypeError(f"unsupported viewer array dtype: {dtype}")


def _json_compatible(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _json_compatible(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_compatible(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    return value


def _asset_text(name: str) -> str:
    resource = files("hic_basic").joinpath("data", "viewer", name)
    try:
        return resource.read_text(encoding="utf-8")
    except FileNotFoundError as error:
        raise RuntimeError(f"installed hic_basic viewer resource is missing: {name}") from error


def _html_document(
    *,
    title: str,
    payload_json: str,
    bundle: str,
    styles: str,
    background: str,
    particle_size: float,
    bond_size: float,
    export_profile: str,
) -> str:
    profile_options = [
        ("molecular-studio-v1", "Molecular studio v1"),
        ("section-studio-v1", "Section studio v1"),
        ("section-pymol-v1", "Section PyMOL v1"),
        ("xu2025-v2", "Xu 2025 v2"),
        ("xu2025-v1", "Xu 2025 v1 (legacy)"),
        ("publication-v1", "Publication v1 (legacy)"),
    ]
    options = "\n".join(
        f'          <option value="{value}"{" selected" if value == export_profile else ""}>{label}</option>'
        for value, label in profile_options
    )
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width,initial-scale=1">
  <title>{title}</title>
  <style>{styles}</style>
</head>
<body data-ready="0" style="background:{background}">
  <main id="app">
    <div id="viewport" aria-label="Interactive 3D genome structure"></div>
    <header class="topbar">
      <div class="identity">
        <h1>{title}</h1>
        <div class="dataset-meta">
          <span id="dataset-name"></span>
          <span id="resolution"></span>
          <span id="particle-count"></span>
          <span id="bond-count"></span>
        </div>
      </div>
      <div class="controls" aria-label="Viewer controls">
        <label class="toggle"><input id="particles-toggle" type="checkbox" checked><span>Particles</span></label>
        <label class="toggle"><input id="bonds-toggle" type="checkbox" checked><span>Polymer bonds</span></label>
        <label class="toggle"><input id="section-toggle" type="checkbox"><span>Section</span></label>
        <label class="field field-section"><span>Thickness</span>
          <input id="section-thickness" type="range" min="0.1" max="10" step="0.1" value="2">
          <output id="section-thickness-value">2</output>
        </label>
        <label class="field field-section"><span>Offset</span>
          <input id="section-offset" type="range" min="-10" max="10" step="0.1" value="0">
          <output id="section-offset-value">0</output>
        </label>
        <button id="regions-button" class="region-button" type="button"
                aria-expanded="false" aria-controls="region-panel">
          Regions <span id="region-summary" class="control-value">All</span>
        </button>
        <label class="field field-size"><span>Particle size</span>
          <input id="point-size" type="range" min="0.5" max="3" step="0.25" value="{particle_size:g}">
        </label>
        <label class="field field-size"><span>Bond size</span>
          <input id="bond-size" type="range" min="1" max="6" step="0.5" value="{bond_size:g}">
        </label>
        <label class="field field-color" title="Monitor background"><span>Background</span>
          <input id="background-color" type="color" value="{background}">
        </label>
        <label class="field field-profile"><span>PNG style</span><select id="export-profile">
{options}
        </select></label>
        <button id="fit-view" type="button" title="Fit the active structure to view">Fit</button>
        <button id="save-png" type="button" title="Render a publication PNG">PNG</button>
      </div>
    </header>
    <section id="region-panel" class="region-panel" aria-label="Genome regions" hidden>
      <form id="region-form" class="region-form">
        <label for="region-input">Genome region</label>
        <div class="region-entry">
          <input id="region-input" type="text" placeholder="chr1(mat):0-50Mb" autocomplete="off">
          <button type="submit">Add</button>
        </div>
      </form>
      <output id="region-error" class="region-error" aria-live="polite"></output>
      <div id="region-list" class="region-list"></div>
      <button id="clear-regions" class="clear-regions" type="button">All genome</button>
    </section>
    <footer class="bottombar">
      <output id="fps-status">Loading</output>
      <output id="renderer-status">Three.js</output>
      <output id="selection-status">Preparing scene</output>
    </footer>
    <div id="tooltip" class="tooltip" role="status"></div>
    <div id="error-state" class="error-state" role="alert"></div>
  </main>
  <script id="scene-data" type="application/json">{payload_json}</script>
  <script>{bundle}</script>
</body>
</html>
"""
