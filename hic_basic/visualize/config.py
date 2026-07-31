"""Typed, serializable configuration for direct structure rendering."""

from __future__ import annotations

import json
import math
import numbers
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence, Tuple, Union

from matplotlib import colors as mpl_colors


RENDER_CONFIG_SCHEMA_VERSION = 1
EXPORT_PROFILES = (
    "molecular-studio-v1",
    "section-studio-v1",
    "section-pymol-v1",
    "xu2025-v2",
    "xu2025-v1",
    "publication-v1",
)


def _strict_mapping(value: Any, allowed: set[str], name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise TypeError(f"{name} must be a mapping")
    unknown = sorted(set(value) - allowed)
    if unknown:
        raise ValueError(f"unknown {name} fields: {', '.join(unknown)}")
    return value


def _finite_number(value: Any, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite number")
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be a finite number") from error
    if not math.isfinite(result):
        raise ValueError(f"{name} must be a finite number")
    return result


def _integer(value: Any, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, numbers.Integral):
        raise ValueError(f"{name} must be an integer")
    return int(value)


def _vector3(value: Optional[Sequence[float]], name: str) -> Optional[Tuple[float, float, float]]:
    if value is None:
        return None
    try:
        valid_length = not isinstance(value, (str, bytes)) and len(value) == 3
    except TypeError:
        valid_length = False
    if not valid_length:
        raise ValueError(f"{name} must contain exactly three numbers")
    return tuple(_finite_number(item, name) for item in value)


def _positive_width(value: Optional[float], name: str) -> Optional[float]:
    if value is None:
        return None
    result = _finite_number(value, name)
    if not 0 < result <= 64:
        raise ValueError(f"{name} must be greater than 0 and at most 64 reference pixels")
    return result


@dataclass(frozen=True)
class CameraConfig:
    """Camera state for a structure render.

    ``mode="fit"`` deterministically fits the selected structure, optionally
    viewing it from ``direction``. With ``mode="explicit"``, ``position`` and
    ``target`` are required.
    """

    mode: str = "fit"
    position: Optional[Tuple[float, float, float]] = None
    target: Optional[Tuple[float, float, float]] = None
    up: Tuple[float, float, float] = (0.0, 1.0, 0.0)
    fov: float = 34.0
    direction: Optional[Tuple[float, float, float]] = None

    def __post_init__(self) -> None:
        if self.mode not in {"fit", "explicit"}:
            raise ValueError("camera mode must be 'fit' or 'explicit'")
        position = _vector3(self.position, "camera position")
        target = _vector3(self.target, "camera target")
        up = _vector3(self.up, "camera up")
        direction = _vector3(self.direction, "camera direction")
        if up is None:
            raise ValueError("camera up must contain exactly three numbers")
        fov = _finite_number(self.fov, "camera fov")
        if not 1 <= fov <= 179:
            raise ValueError("camera fov must be between 1 and 179 degrees")
        if up == (0.0, 0.0, 0.0):
            raise ValueError("camera up vector must be non-zero")
        if direction == (0.0, 0.0, 0.0):
            raise ValueError("camera direction vector must be non-zero")
        if self.mode == "explicit":
            if direction is not None:
                raise ValueError("explicit camera does not accept direction")
            if position is None or target is None:
                raise ValueError("explicit camera requires position and target")
            if position == target:
                raise ValueError("camera position and target must differ")
            view_vector = tuple(left - right for left, right in zip(position, target))
            cross = (
                view_vector[1] * up[2] - view_vector[2] * up[1],
                view_vector[2] * up[0] - view_vector[0] * up[2],
                view_vector[0] * up[1] - view_vector[1] * up[0],
            )
            if sum(value * value for value in cross) <= 1e-20:
                raise ValueError("camera up vector must not be parallel to the viewing direction")
        else:
            if position is not None or target is not None:
                raise ValueError("fit camera does not accept position or target")
            if direction is not None:
                cross = (
                    direction[1] * up[2] - direction[2] * up[1],
                    direction[2] * up[0] - direction[0] * up[2],
                    direction[0] * up[1] - direction[1] * up[0],
                )
                if sum(value * value for value in cross) <= 1e-20:
                    raise ValueError("camera up vector must not be parallel to the viewing direction")
        object.__setattr__(self, "position", position)
        object.__setattr__(self, "target", target)
        object.__setattr__(self, "up", up)
        object.__setattr__(self, "fov", fov)
        object.__setattr__(self, "direction", direction)

    def to_dict(self) -> dict:
        return {
            "mode": self.mode,
            "position": list(self.position) if self.position is not None else None,
            "target": list(self.target) if self.target is not None else None,
            "up": list(self.up),
            "fov": self.fov,
            "direction": list(self.direction) if self.direction is not None else None,
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "CameraConfig":
        data = _strict_mapping(
            value,
            {"mode", "position", "target", "up", "fov", "direction"},
            "camera config",
        )
        return cls(**dict(data))


@dataclass(frozen=True)
class SectionConfig:
    """Camera-aligned slab used for cross-section rendering.

    ``thickness`` and ``offset`` use the coordinate units of the 3dg input.
    The slab is centered on the camera target plus ``offset`` toward the
    camera. A disabled section preserves the complete surface view.
    """

    enabled: bool = False
    thickness: float = 2.0
    offset: float = 0.0

    def __post_init__(self) -> None:
        if not isinstance(self.enabled, bool):
            raise TypeError("section enabled must be bool")
        thickness = _finite_number(self.thickness, "section thickness")
        offset = _finite_number(self.offset, "section offset")
        if thickness <= 0:
            raise ValueError("section thickness must be greater than 0")
        object.__setattr__(self, "thickness", thickness)
        object.__setattr__(self, "offset", offset)

    def to_dict(self) -> dict:
        return {
            "enabled": self.enabled,
            "thickness": self.thickness,
            "offset": self.offset,
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "SectionConfig":
        data = _strict_mapping(value, {"enabled", "thickness", "offset"}, "section config")
        return cls(**dict(data))


@dataclass(frozen=True)
class ViewConfig:
    """Genome selection, camera, and optional cross-section slab."""

    regions: Tuple[str, ...] = ()
    camera: CameraConfig = field(default_factory=CameraConfig)
    section: SectionConfig = field(default_factory=SectionConfig)

    def __post_init__(self) -> None:
        if isinstance(self.regions, (str, bytes)):
            raise TypeError("regions must be a sequence of genome-region strings")
        regions = tuple(str(region).strip() for region in self.regions)
        if any(not region for region in regions):
            raise ValueError("render regions must not contain empty values")
        if not isinstance(self.camera, CameraConfig):
            raise TypeError("camera must be a CameraConfig")
        if not isinstance(self.section, SectionConfig):
            raise TypeError("section must be a SectionConfig")
        object.__setattr__(self, "regions", regions)

    def to_dict(self) -> dict:
        return {
            "regions": list(self.regions),
            "camera": self.camera.to_dict(),
            "section": self.section.to_dict(),
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "ViewConfig":
        data = _strict_mapping(value, {"regions", "camera", "section"}, "view config")
        camera = CameraConfig.from_dict(data.get("camera", {}))
        section = SectionConfig.from_dict(data.get("section", {}))
        return cls(regions=tuple(data.get("regions", ())), camera=camera, section=section)


@dataclass(frozen=True)
class LayerConfig:
    """Visibility overrides for base geometry and structure annotations."""

    show_particles: Optional[bool] = None
    show_bonds: Optional[bool] = None
    show_annotations: Optional[bool] = None

    def __post_init__(self) -> None:
        for name in ("show_particles", "show_bonds", "show_annotations"):
            value = getattr(self, name)
            if value is not None and not isinstance(value, bool):
                raise TypeError(f"{name} must be bool or None")
        if (
            self.show_particles is False
            and self.show_bonds is False
            and self.show_annotations is not True
        ):
            raise ValueError("at least one render layer must be visible")

    def to_dict(self) -> dict:
        return {
            "show_particles": self.show_particles,
            "show_bonds": self.show_bonds,
            "show_annotations": self.show_annotations,
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "LayerConfig":
        data = _strict_mapping(
            value,
            {"show_particles", "show_bonds", "show_annotations"},
            "layer config",
        )
        return cls(**dict(data))


@dataclass(frozen=True)
class ExportStyle:
    """Versioned profile plus explicit high-quality rendering overrides."""

    profile: str = "molecular-studio-v1"
    tube_width: Optional[float] = None
    silhouette_width: Optional[float] = None
    particle_size: Optional[float] = None
    background: Optional[str] = None
    show_colorbar: Optional[bool] = None

    def __post_init__(self) -> None:
        if self.profile not in EXPORT_PROFILES:
            choices = ", ".join(EXPORT_PROFILES)
            raise ValueError(f"export profile must be one of: {choices}")
        background = None
        if self.background is not None:
            background = mpl_colors.to_hex(self.background, keep_alpha=False)
        if self.show_colorbar is not None and not isinstance(self.show_colorbar, bool):
            raise TypeError("show_colorbar must be bool or None")
        object.__setattr__(self, "tube_width", _positive_width(self.tube_width, "tube_width"))
        object.__setattr__(
            self,
            "silhouette_width",
            _positive_width(self.silhouette_width, "silhouette_width"),
        )
        object.__setattr__(
            self,
            "particle_size",
            _positive_width(self.particle_size, "particle_size"),
        )
        object.__setattr__(self, "background", background)

    def to_dict(self) -> dict:
        return {
            "profile": self.profile,
            "tube_width": self.tube_width,
            "silhouette_width": self.silhouette_width,
            "particle_size": self.particle_size,
            "background": self.background,
            "show_colorbar": self.show_colorbar,
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "ExportStyle":
        allowed = {
            "profile",
            "tube_width",
            "silhouette_width",
            "particle_size",
            "background",
            "show_colorbar",
        }
        data = _strict_mapping(value, allowed, "export style")
        return cls(**dict(data))


@dataclass(frozen=True)
class RasterConfig:
    """Output raster dimensions.

    When ``width`` and ``height`` are omitted, ``long_edge`` is combined with
    the browser viewport aspect ratio. Supplying explicit dimensions is the
    most reproducible direct-rendering mode.
    """

    long_edge: int = 2400
    width: Optional[int] = None
    height: Optional[int] = None

    def __post_init__(self) -> None:
        long_edge = _integer(self.long_edge, "long_edge")
        if not 64 <= long_edge <= 8192:
            raise ValueError("long_edge must be between 64 and 8192 pixels")
        if (self.width is None) != (self.height is None):
            raise ValueError("raster width and height must be provided together")
        width = None if self.width is None else _integer(self.width, "raster width")
        height = None if self.height is None else _integer(self.height, "raster height")
        if width is not None and (not 64 <= width <= 8192 or not 64 <= height <= 8192):
            raise ValueError("raster width and height must be between 64 and 8192 pixels")
        object.__setattr__(self, "long_edge", long_edge)
        object.__setattr__(self, "width", width)
        object.__setattr__(self, "height", height)

    def to_dict(self) -> dict:
        return {"long_edge": self.long_edge, "width": self.width, "height": self.height}

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "RasterConfig":
        data = _strict_mapping(value, {"long_edge", "width", "height"}, "raster config")
        return cls(**dict(data))


@dataclass(frozen=True)
class StructureRenderConfig:
    """Complete reproducible request shared by monitor and direct PNG export."""

    view: ViewConfig = field(default_factory=ViewConfig)
    layers: LayerConfig = field(default_factory=LayerConfig)
    style: ExportStyle = field(default_factory=ExportStyle)
    raster: RasterConfig = field(default_factory=RasterConfig)
    schema_version: int = RENDER_CONFIG_SCHEMA_VERSION

    def __post_init__(self) -> None:
        schema_version = _integer(self.schema_version, "schema_version")
        if schema_version != RENDER_CONFIG_SCHEMA_VERSION:
            raise ValueError(
                f"unsupported render config schema_version: {schema_version}; "
                f"expected {RENDER_CONFIG_SCHEMA_VERSION}"
            )
        if not isinstance(self.view, ViewConfig):
            raise TypeError("view must be a ViewConfig")
        if not isinstance(self.layers, LayerConfig):
            raise TypeError("layers must be a LayerConfig")
        if not isinstance(self.style, ExportStyle):
            raise TypeError("style must be an ExportStyle")
        if not isinstance(self.raster, RasterConfig):
            raise TypeError("raster must be a RasterConfig")
        object.__setattr__(self, "schema_version", schema_version)

    def to_dict(self) -> dict:
        return {
            "schema_version": self.schema_version,
            "view": self.view.to_dict(),
            "layers": self.layers.to_dict(),
            "style": self.style.to_dict(),
            "raster": self.raster.to_dict(),
        }

    def to_json(self, *, indent: Optional[int] = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent, sort_keys=True)

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "StructureRenderConfig":
        allowed = {"schema_version", "view", "layers", "style", "raster"}
        data = _strict_mapping(value, allowed, "structure render config")
        return cls(
            schema_version=data.get("schema_version", RENDER_CONFIG_SCHEMA_VERSION),
            view=ViewConfig.from_dict(data.get("view", {})),
            layers=LayerConfig.from_dict(data.get("layers", {})),
            style=ExportStyle.from_dict(data.get("style", {})),
            raster=RasterConfig.from_dict(data.get("raster", {})),
        )

    @classmethod
    def from_json(cls, value: str) -> "StructureRenderConfig":
        payload = json.loads(value)
        return cls.from_dict(payload)


@dataclass(frozen=True)
class BrowserBackendConfig:
    """Execution-only settings for the pinned headless browser backend."""

    executable_path: Optional[Union[str, os.PathLike]] = None
    playwright_path: Optional[Union[str, os.PathLike]] = None
    timeout_seconds: float = 180.0
    software_rendering: bool = True
    viewport: Tuple[int, int] = (1440, 900)
    extra_args: Tuple[str, ...] = ()

    def __post_init__(self) -> None:
        timeout = _finite_number(self.timeout_seconds, "browser timeout")
        if not 1 <= timeout <= 3600:
            raise ValueError("browser timeout must be between 1 and 3600 seconds")
        if not isinstance(self.software_rendering, bool):
            raise TypeError("software_rendering must be bool")
        if isinstance(self.viewport, (str, bytes)) or len(self.viewport) != 2:
            raise ValueError("browser viewport must contain width and height")
        viewport = tuple(_integer(value, "browser viewport dimension") for value in self.viewport)
        if any(value < 64 or value > 8192 for value in viewport):
            raise ValueError("browser viewport dimensions must be between 64 and 8192")
        object.__setattr__(self, "timeout_seconds", timeout)
        object.__setattr__(self, "viewport", viewport)
        if isinstance(self.extra_args, (str, bytes)):
            raise TypeError("extra_args must be a sequence of browser arguments")
        object.__setattr__(self, "extra_args", tuple(str(value) for value in self.extra_args))
        if self.executable_path is not None:
            object.__setattr__(self, "executable_path", str(Path(self.executable_path).expanduser()))
        if self.playwright_path is not None:
            object.__setattr__(self, "playwright_path", str(Path(self.playwright_path).expanduser()))

    def to_dict(self) -> dict:
        return {
            "executable_path": self.executable_path,
            "timeout_ms": round(self.timeout_seconds * 1000),
            "software_rendering": self.software_rendering,
            "viewport": list(self.viewport),
            "extra_args": list(self.extra_args),
        }


RenderConfigInput = Optional[
    Union[StructureRenderConfig, Mapping[str, Any], str, os.PathLike]
]


def load_render_config(value: RenderConfigInput = None) -> StructureRenderConfig:
    """Load and validate a structure render configuration.

    ``value`` may be an existing config, a mapping, a JSON string, a JSON file,
    or ``None`` for defaults. Unknown fields and invalid values are rejected.
    """
    if value is None:
        return StructureRenderConfig()
    if isinstance(value, StructureRenderConfig):
        return value
    if isinstance(value, Mapping):
        return StructureRenderConfig.from_dict(value)
    serialized = os.fspath(value)
    if isinstance(serialized, str) and serialized.lstrip().startswith("{"):
        return StructureRenderConfig.from_json(serialized)
    path = Path(serialized).expanduser()
    try:
        is_file = path.is_file()
    except OSError:
        is_file = False
    if is_file:
        return StructureRenderConfig.from_json(path.read_text(encoding="utf-8"))
    return StructureRenderConfig.from_json(str(serialized))
