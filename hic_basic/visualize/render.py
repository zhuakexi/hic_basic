"""Direct PNG rendering for 3D genome scenes and standalone viewer HTML."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path
from typing import Any, Mapping, Optional, Union

import pandas as pd
from PIL import Image

from .config import BrowserBackendConfig, RenderConfigInput, StructureRenderConfig, load_render_config
from .html import RENDERER_PROTOCOL_VERSION, StructureViewer, view_structure
from .scene import StructureScene
from .structure import StructureInput, structure_scene


RenderSource = Union[StructureInput, StructureScene, StructureViewer]


@dataclass(frozen=True)
class RenderResult:
    """Validated output and provenance returned by a direct PNG render."""

    output: Path
    width: int
    height: int
    size_bytes: int
    profile: str
    scene_sha256: str
    metadata: Mapping[str, Any]

    def to_dict(self) -> dict:
        """Return a JSON-compatible rendering report."""
        return {
            "output": str(self.output),
            "width": self.width,
            "height": self.height,
            "size_bytes": self.size_bytes,
            "profile": self.profile,
            "scene_sha256": self.scene_sha256,
            **dict(self.metadata),
        }


def render_structure_png(
    source: RenderSource,
    output: Union[str, os.PathLike],
    *,
    config: RenderConfigInput = None,
    backend: Optional[BrowserBackendConfig] = None,
    overwrite: bool = True,
) -> RenderResult:
    """Render a 3dg structure, scene, viewer, or standalone HTML to PNG.

    The function uses the repository-pinned Playwright runtime and a cached
    headless Chromium. It does not start a server, open a GUI, use the network,
    or download browser assets. ``config`` accepts a ``StructureRenderConfig``,
    a strict mapping, JSON text, or a JSON file path. A ``StructureViewer``
    supplies its embedded config when ``config`` is omitted.

    Parameters
    ----------
    source
        A 3dg path, native structure DataFrame, ``StructureScene``,
        ``StructureViewer``, or compatible standalone viewer HTML.
    output
        Destination PNG path. The image is validated before atomic replacement.
    config
        Reproducible region, camera, cross-section slab, layer, profile,
        colorbar, style, and raster settings.
    backend
        Optional browser executable/module paths, timeout, viewport, and
        software-rendering policy. These execution details are kept separate
        from the scientific render configuration.
    overwrite
        Whether an existing output may be replaced.

    Returns
    -------
    RenderResult
        Output dimensions, file size, scene checksum, and backend/render report.

    Notes
    -----
    Install the matching Chromium asset once with
    ``./run_in_env.sh --network internet playwright install --only-shell chromium``
    when the backend reports that its executable is missing.
    """
    output_path = Path(output).expanduser()
    if output_path.exists() and not overwrite:
        raise FileExistsError(f"output already exists: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    backend_config = backend or BrowserBackendConfig()
    if not isinstance(backend_config, BrowserBackendConfig):
        raise TypeError("backend must be a BrowserBackendConfig")

    supplied_config = None if config is None else load_render_config(config)
    html_source, temporary_html, effective_config = _prepare_html_source(source, supplied_config, output_path.parent)
    temporary_png = None
    request_path = None
    started_at = time.monotonic()
    try:
        with tempfile.NamedTemporaryFile(
            prefix=f".{output_path.name}.", suffix=".png", dir=output_path.parent, delete=False
        ) as stream:
            temporary_png = Path(stream.name)
        temporary_png.unlink()

        playwright_module = _playwright_module(backend_config)
        node = shutil.which("node")
        if node is None:
            raise RuntimeError("Node.js is missing from the active hic_basic environment")
        driver = files("hic_basic").joinpath("data", "viewer", "render_driver.cjs")
        launch_args = list(backend_config.extra_args)
        if backend_config.software_rendering:
            launch_args = [
                "--use-angle=swiftshader",
                "--enable-unsafe-swiftshader",
                "--disable-gpu-sandbox",
                *launch_args,
            ]
        request = {
            "protocol_version": RENDERER_PROTOCOL_VERSION,
            "html_path": str(html_source.resolve()),
            "png_path": str(temporary_png.resolve()),
            "render_config": effective_config.to_dict() if effective_config is not None else None,
            "browser": {
                **backend_config.to_dict(),
                "args": launch_args,
            },
        }
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=f".{output_path.name}.",
            suffix=".render.json",
            dir=output_path.parent,
            delete=False,
        ) as stream:
            request_path = Path(stream.name)
            json.dump(request, stream, ensure_ascii=True)

        environment = os.environ.copy()
        environment["HIC_BASIC_PLAYWRIGHT_MODULE"] = str(playwright_module)
        library_path = str(Path(sys.prefix) / "lib")
        existing_library_path = environment.get("LD_LIBRARY_PATH")
        environment["LD_LIBRARY_PATH"] = (
            f"{library_path}{os.pathsep}{existing_library_path}"
            if existing_library_path
            else library_path
        )
        completed = subprocess.run(
            [node, str(driver), str(request_path)],
            capture_output=True,
            text=True,
            check=False,
            timeout=backend_config.timeout_seconds + 15,
            env=environment,
        )
        if completed.returncode:
            _raise_backend_error(completed.stderr or completed.stdout)
        try:
            metadata = json.loads(completed.stdout.strip())
        except json.JSONDecodeError as error:
            raise RuntimeError(f"renderer returned invalid metadata: {completed.stdout!r}") from error
        _validate_png(temporary_png, metadata)
        os.replace(temporary_png, output_path)
        temporary_png = None
        metadata["total_time_ms"] = round((time.monotonic() - started_at) * 1000)
        result_metadata = _normalize_metadata(metadata)
        return RenderResult(
            output=output_path,
            width=int(metadata.pop("width")),
            height=int(metadata.pop("height")),
            size_bytes=output_path.stat().st_size,
            profile=str(metadata.pop("profile")),
            scene_sha256=str(metadata.pop("sceneSha256", "")),
            metadata=result_metadata,
        )
    except subprocess.TimeoutExpired as error:
        raise RuntimeError(
            f"structure rendering exceeded the {backend_config.timeout_seconds:g} second browser timeout"
        ) from error
    finally:
        for path in (request_path, temporary_png, temporary_html):
            if path is not None and Path(path).exists():
                Path(path).unlink()


def _prepare_html_source(
    source: RenderSource,
    config: Optional[StructureRenderConfig],
    temporary_directory: Path,
) -> tuple[Path, Optional[Path], Optional[StructureRenderConfig]]:
    if isinstance(source, StructureViewer):
        effective_config = config or source.render_config
        viewer = source if config is None else view_structure(source.scene, title=source.title, render_config=config)
    elif isinstance(source, StructureScene):
        effective_config = config or StructureRenderConfig()
        viewer = view_structure(source, render_config=effective_config)
    elif isinstance(source, pd.DataFrame):
        effective_config = config or StructureRenderConfig()
        viewer = view_structure(structure_scene(source), render_config=effective_config)
    else:
        path = Path(source).expanduser()
        if path.suffix.lower() in {".html", ".htm"}:
            if not path.is_file():
                raise FileNotFoundError(f"standalone structure HTML does not exist: {path}")
            return path, None, config
        effective_config = config or StructureRenderConfig()
        viewer = view_structure(structure_scene(path), render_config=effective_config)

    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        prefix=".hic-basic-render.",
        suffix=".html",
        dir=temporary_directory,
        delete=False,
    ) as stream:
        temporary_html = Path(stream.name)
        stream.write(viewer.to_html())
    return temporary_html, temporary_html, effective_config


def _playwright_module(backend: BrowserBackendConfig) -> Path:
    if backend.playwright_path is not None:
        candidate = Path(backend.playwright_path)
        if candidate.is_file():
            candidate = candidate.resolve()
        if candidate.name == "cli.js":
            candidate = candidate.parent
    else:
        candidate = Path(sys.prefix) / "lib" / "node_modules" / "playwright"
        if not candidate.is_dir():
            executable = shutil.which("playwright")
            if executable is not None:
                candidate = Path(executable).resolve().parent
    if not (candidate / "package.json").is_file():
        raise RuntimeError(
            "Playwright is missing from the active hic_basic environment; rebuild it from the registered lock"
        )
    if backend.executable_path is not None and not Path(backend.executable_path).is_file():
        raise FileNotFoundError(f"Chromium executable does not exist: {backend.executable_path}")
    return candidate


def _validate_png(path: Path, metadata: Mapping[str, Any]) -> None:
    if not path.is_file():
        raise RuntimeError("browser completed without creating a PNG")
    try:
        with Image.open(path) as image:
            if image.format != "PNG":
                raise RuntimeError(f"browser output is {image.format}, not PNG")
            dimensions = image.size
            image.verify()
    except Exception as error:
        raise RuntimeError(f"browser produced an invalid PNG: {error}") from error
    expected = (int(metadata["width"]), int(metadata["height"]))
    if dimensions != expected:
        raise RuntimeError(f"PNG dimensions {dimensions} do not match renderer report {expected}")


def _normalize_metadata(metadata: Mapping[str, Any]) -> dict:
    renderer = metadata.get("renderer", {})
    backend = metadata.get("backend", {})
    return {
        "render_passes": int(metadata.get("renderPasses", 0)),
        "geometry": metadata.get("geometryStats"),
        "request": metadata.get("request", {}),
        "selection": metadata.get("selection", {}),
        "annotations": metadata.get("annotations", []),
        "section": metadata.get("section"),
        "colorbar": metadata.get("colorbar"),
        "camera": metadata.get("camera", {}),
        "renderer": {
            "three_revision": renderer.get("threeRevision"),
            "software": renderer.get("software"),
            "label": renderer.get("label"),
        },
        "backend": {
            "browser_version": backend.get("browserVersion"),
            "playwright_version": backend.get("playwrightVersion"),
            "protocol_version": backend.get("protocolVersion"),
            "render_config_schema_version": backend.get("renderConfigSchemaVersion"),
        },
        "render_time_ms": int(metadata.get("renderTimeMs", 0)),
        "total_time_ms": int(metadata.get("total_time_ms", 0)),
    }


def _raise_backend_error(detail: str) -> None:
    clean_detail = detail.strip()
    if "Executable doesn't exist" in clean_detail or "browserType.launch" in clean_detail:
        raise RuntimeError(
            "the pinned Chromium asset is unavailable or could not start. Install it once with "
            "'./run_in_env.sh --network internet playwright install --only-shell chromium'. "
            f"Backend detail: {clean_detail}"
        )
    raise RuntimeError(f"headless structure renderer failed: {clean_detail}")
