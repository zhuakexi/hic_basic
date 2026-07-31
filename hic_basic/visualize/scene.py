"""Validated, renderer-independent scene objects for 3D genome structures."""

from __future__ import annotations

from dataclasses import dataclass, field
import math
import numbers
from typing import Any, Mapping, Optional, Tuple, Union

import numpy as np
from matplotlib import colors as mpl_colors


def _readonly_array(value: Any, dtype: np.dtype, shape_tail: Tuple[int, ...], name: str) -> np.ndarray:
    array = np.asarray(value, dtype=dtype)
    expected_ndim = len(shape_tail) + 1
    if array.ndim != expected_ndim or tuple(array.shape[1:]) != shape_tail:
        expected = "N" + "".join(f" x {size}" for size in shape_tail)
        raise ValueError(f"{name} must have shape ({expected}); got {array.shape}")
    array = np.ascontiguousarray(array)
    array.setflags(write=False)
    return array


def _readonly_unsigned(
    value: Any,
    dtype: np.dtype,
    shape_tail: Tuple[int, ...],
    name: str,
) -> np.ndarray:
    raw = np.asarray(value)
    expected_ndim = len(shape_tail) + 1
    if raw.ndim != expected_ndim or tuple(raw.shape[1:]) != shape_tail:
        expected = "N" + "".join(f" x {size}" for size in shape_tail)
        raise ValueError(f"{name} must have shape ({expected}); got {raw.shape}")
    if raw.dtype.kind not in "uif":
        raise ValueError(f"{name} must contain non-negative integers")
    if raw.dtype.kind == "f" and (not np.isfinite(raw).all() or not np.equal(raw, np.floor(raw)).all()):
        raise ValueError(f"{name} must contain finite integers")
    if raw.size and (raw.min() < 0 or raw.max() > np.iinfo(dtype).max):
        raise ValueError(f"{name} values do not fit {np.dtype(dtype).name}")
    array = np.ascontiguousarray(raw, dtype=dtype)
    array.setflags(write=False)
    return array


@dataclass(frozen=True)
class ParticleLayer:
    """Particle coordinates and genomic identity in display order."""

    positions: np.ndarray
    colors: np.ndarray
    chromosome_codes: np.ndarray
    genomic_positions: np.ndarray

    def __post_init__(self) -> None:
        positions = _readonly_array(self.positions, np.dtype("<f4"), (3,), "positions")
        colors = _readonly_unsigned(self.colors, np.dtype("u1"), (3,), "colors")
        raw_codes = np.asarray(self.chromosome_codes)
        if raw_codes.ndim != 1:
            raise ValueError("chromosome_codes must contain one value per particle")
        maximum_code = int(raw_codes.max()) if raw_codes.size else 0
        if maximum_code <= np.iinfo(np.uint8).max:
            code_dtype = np.dtype("u1")
        elif maximum_code <= np.iinfo(np.uint16).max:
            code_dtype = np.dtype("<u2")
        else:
            code_dtype = np.dtype("<u4")
        chromosome_codes = _readonly_unsigned(raw_codes, code_dtype, (), "chromosome_codes")
        genomic_positions = _readonly_unsigned(
            self.genomic_positions,
            np.dtype("<u4"),
            (),
            "genomic_positions",
        )

        count = len(positions)
        if len(colors) != count:
            raise ValueError("colors must contain one RGB value per particle")
        if chromosome_codes.ndim != 1 or len(chromosome_codes) != count:
            raise ValueError("chromosome_codes must contain one value per particle")
        if genomic_positions.ndim != 1 or len(genomic_positions) != count:
            raise ValueError("genomic_positions must contain one value per particle")
        if not np.isfinite(positions).all():
            raise ValueError("positions contain non-finite coordinates")
        object.__setattr__(self, "positions", positions)
        object.__setattr__(self, "colors", colors)
        object.__setattr__(self, "chromosome_codes", chromosome_codes)
        object.__setattr__(self, "genomic_positions", genomic_positions)

    @property
    def count(self) -> int:
        """Number of particles in the layer."""
        return len(self.positions)


@dataclass(frozen=True)
class PolymerBondLayer:
    """Polymer-chain bonds expressed as particle endpoint indices."""

    endpoints: np.ndarray
    max_genomic_gap: int

    def __post_init__(self) -> None:
        endpoints = _readonly_unsigned(self.endpoints, np.dtype("<u4"), (2,), "endpoints")
        if int(self.max_genomic_gap) < 0:
            raise ValueError("max_genomic_gap must be non-negative")
        object.__setattr__(self, "endpoints", endpoints)
        object.__setattr__(self, "max_genomic_gap", int(self.max_genomic_gap))

    @property
    def count(self) -> int:
        """Number of polymer bonds in the layer."""
        return len(self.endpoints)


def _annotation_id(value: Any, name: str) -> str:
    result = str(value).strip()
    if not result:
        raise ValueError(f"{name} must not be empty")
    return result


def _annotation_order(value: Any, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, numbers.Integral):
        raise ValueError(f"{name} must be an integer")
    return int(value)


def _annotation_width(value: Any, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be a finite number") from error
    if not math.isfinite(result) or not 0 < result <= 64:
        raise ValueError(f"{name} must be greater than 0 and at most 64 reference pixels")
    return result


@dataclass(frozen=True)
class GenomeRegionTarget:
    """Unresolved inclusive genome-region target for a structure annotation."""

    selector: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "selector", _annotation_id(self.selector, "target selector"))


@dataclass(frozen=True)
class ResolvedGenomeTarget:
    """Genome target resolved against a scene's stable particle and bond order."""

    selector: str
    chromosomes: Tuple[str, ...]
    particle_indices: np.ndarray
    bond_indices: np.ndarray
    anchor: Tuple[float, float, float]
    bbox: Tuple[Tuple[float, float, float], Tuple[float, float, float]]

    def __post_init__(self) -> None:
        selector = _annotation_id(self.selector, "target selector")
        chromosomes = tuple(_annotation_id(value, "target chromosome") for value in self.chromosomes)
        particle_indices = _readonly_unsigned(
            self.particle_indices,
            np.dtype("<u4"),
            (),
            "annotation particle_indices",
        )
        bond_indices = _readonly_unsigned(
            self.bond_indices,
            np.dtype("<u4"),
            (),
            "annotation bond_indices",
        )
        anchor = tuple(float(value) for value in self.anchor)
        bbox = tuple(tuple(float(value) for value in point) for point in self.bbox)
        if len(anchor) != 3 or not np.isfinite(anchor).all():
            raise ValueError("annotation target anchor must contain three finite coordinates")
        if len(bbox) != 2 or any(len(point) != 3 for point in bbox) or not np.isfinite(bbox).all():
            raise ValueError("annotation target bbox must contain two finite 3D coordinates")
        if particle_indices.size == 0:
            raise ValueError("annotation target must resolve to at least one particle")
        if particle_indices.size > 1 and np.any(np.diff(particle_indices.astype(np.int64)) <= 0):
            raise ValueError("annotation particle_indices must be strictly increasing")
        if bond_indices.size > 1 and np.any(np.diff(bond_indices.astype(np.int64)) <= 0):
            raise ValueError("annotation bond_indices must be strictly increasing")
        object.__setattr__(self, "selector", selector)
        object.__setattr__(self, "chromosomes", chromosomes)
        object.__setattr__(self, "particle_indices", particle_indices)
        object.__setattr__(self, "bond_indices", bond_indices)
        object.__setattr__(self, "anchor", anchor)
        object.__setattr__(self, "bbox", bbox)


AnnotationTarget = Union[GenomeRegionTarget, ResolvedGenomeTarget]


@dataclass(frozen=True)
class TubeMark:
    """Tube mark following the real polymer bonds of an annotation target.

    Widths are diameters in reference pixels at a 1400-pixel export width.
    ``material="unlit"`` preserves the requested semantic sRGB color, while
    ``material="profile"`` opts the primary tube into the active base export
    profile. An optional outline remains an unlit semantic halo.
    """

    color: str
    tube_width: float = 8.0
    outline_width: Optional[float] = None
    outline_color: Optional[str] = None
    material: str = "unlit"
    opacity: float = 1.0
    occlusion: str = "depth_tested"
    kind: str = field(default="tube", init=False)
    requested_color: str = field(default="", init=False)
    requested_outline_color: Optional[str] = field(default=None, init=False)

    def __post_init__(self) -> None:
        requested_color = str(self.color)
        requested_outline_color = None if self.outline_color is None else str(self.outline_color)
        try:
            color = mpl_colors.to_hex(self.color, keep_alpha=False)
            outline_color = (
                None
                if self.outline_color is None
                else mpl_colors.to_hex(self.outline_color, keep_alpha=False)
            )
        except ValueError as error:
            raise ValueError(f"invalid annotation color: {error}") from error
        tube_width = _annotation_width(self.tube_width, "annotation tube_width")
        outline_width = tube_width if self.outline_width is None else _annotation_width(
            self.outline_width,
            "annotation outline_width",
        )
        if outline_width < tube_width:
            raise ValueError("annotation outline_width must be at least tube_width")
        if outline_width > tube_width and outline_color is None:
            raise ValueError("annotation outline_color is required when outline_width exceeds tube_width")
        if self.material not in {"unlit", "profile"}:
            raise ValueError("annotation material must be 'unlit' or 'profile'")
        try:
            opacity = float(self.opacity)
        except (TypeError, ValueError) as error:
            raise ValueError("annotation opacity must be a finite number") from error
        if not math.isfinite(opacity) or not 0 < opacity <= 1:
            raise ValueError("annotation opacity must be greater than 0 and at most 1")
        if self.occlusion not in {"depth_tested", "always_visible"}:
            raise ValueError("annotation occlusion must be 'depth_tested' or 'always_visible'")
        object.__setattr__(self, "color", color)
        object.__setattr__(self, "requested_color", requested_color)
        object.__setattr__(self, "tube_width", tube_width)
        object.__setattr__(self, "outline_width", outline_width)
        object.__setattr__(self, "outline_color", outline_color)
        object.__setattr__(self, "requested_outline_color", requested_outline_color)
        object.__setattr__(self, "opacity", opacity)


AnnotationMark = TubeMark


@dataclass(frozen=True)
class StructureAnnotation:
    """One semantic annotation target with one or more typed visual marks."""

    id: str
    target: AnnotationTarget
    marks: Tuple[AnnotationMark, ...]
    label: Optional[str] = None
    order: int = 0

    def __post_init__(self) -> None:
        annotation_id = _annotation_id(self.id, "annotation id")
        if not isinstance(self.target, (GenomeRegionTarget, ResolvedGenomeTarget)):
            raise TypeError("annotation target must be a genome-region target")
        marks = tuple(self.marks)
        if not marks:
            raise ValueError("annotation must contain at least one mark")
        if any(not isinstance(mark, TubeMark) for mark in marks):
            raise TypeError("unsupported annotation mark type")
        kinds = [mark.kind for mark in marks]
        if len(kinds) != len(set(kinds)):
            raise ValueError("annotation cannot contain duplicate mark kinds")
        label = None if self.label is None else str(self.label).strip()
        if label == "":
            raise ValueError("annotation label must not be empty")
        object.__setattr__(self, "id", annotation_id)
        object.__setattr__(self, "marks", marks)
        object.__setattr__(self, "label", label)
        object.__setattr__(self, "order", _annotation_order(self.order, "annotation order"))


@dataclass(frozen=True)
class AnnotationLayer:
    """Named, independently visible collection of structure annotations."""

    id: str
    annotations: Tuple[StructureAnnotation, ...]
    label: Optional[str] = None
    visible: bool = True
    order: int = 0

    def __post_init__(self) -> None:
        layer_id = _annotation_id(self.id, "annotation layer id")
        annotations = tuple(self.annotations)
        if not annotations:
            raise ValueError("annotation layer must contain at least one annotation")
        if any(not isinstance(annotation, StructureAnnotation) for annotation in annotations):
            raise TypeError("annotation layer entries must be StructureAnnotation objects")
        ids = [annotation.id for annotation in annotations]
        if len(ids) != len(set(ids)):
            raise ValueError("annotation ids must be unique within a layer")
        if not isinstance(self.visible, bool):
            raise TypeError("annotation layer visible must be bool")
        label = None if self.label is None else str(self.label).strip()
        if label == "":
            raise ValueError("annotation layer label must not be empty")
        annotations = tuple(sorted(annotations, key=lambda item: (item.order, item.id)))
        object.__setattr__(self, "id", layer_id)
        object.__setattr__(self, "annotations", annotations)
        object.__setattr__(self, "label", label)
        object.__setattr__(self, "order", _annotation_order(self.order, "annotation layer order"))


@dataclass(frozen=True)
class ChromosomeGroup:
    """Contiguous particle and bond spans for one chromosome chain."""

    name: str
    code: int
    start: int
    count: int
    bond_start: int
    bond_count: int
    bbox: Tuple[Tuple[float, float, float], Tuple[float, float, float]]
    color: Tuple[int, int, int]


@dataclass(frozen=True)
class StructureScene:
    """Validated scientific scene shared by visualization backends."""

    particles: ParticleLayer
    bonds: PolymerBondLayer
    chromosomes: Tuple[str, ...]
    groups: Tuple[ChromosomeGroup, ...]
    dataset_name: str
    source: str
    checksum: str
    resolution: Optional[int]
    bbox: Tuple[Tuple[float, float, float], Tuple[float, float, float]]
    base_point_size: float
    metadata: Mapping[str, Any] = field(default_factory=dict)
    annotation_layers: Tuple[AnnotationLayer, ...] = ()
    schema_version: int = 1

    def __post_init__(self) -> None:
        particle_count = self.particles.count
        if particle_count == 0:
            raise ValueError("a StructureScene must contain at least one particle")
        if self.bonds.count and int(self.bonds.endpoints.max()) >= particle_count:
            raise ValueError("bond endpoint refers to a particle outside the scene")
        if len(self.chromosomes) != len(self.groups):
            raise ValueError("chromosomes and groups must have the same length")
        if tuple(group.name for group in self.groups) != tuple(self.chromosomes):
            raise ValueError("group order must match chromosome order")
        if any(group.code != index for index, group in enumerate(self.groups)):
            raise ValueError("chromosome group codes must be contiguous from zero")
        if self.resolution is not None and int(self.resolution) <= 0:
            raise ValueError("resolution must be positive when provided")
        if not np.isfinite(float(self.base_point_size)) or self.base_point_size <= 0:
            raise ValueError("base_point_size must be finite and positive")

        annotation_layers = tuple(self.annotation_layers)
        if any(not isinstance(layer, AnnotationLayer) for layer in annotation_layers):
            raise TypeError("annotation_layers must contain AnnotationLayer objects")
        layer_ids = [layer.id for layer in annotation_layers]
        if len(layer_ids) != len(set(layer_ids)):
            raise ValueError("annotation layer ids must be unique")
        annotation_ids = []
        bond_owners = {}
        for layer in annotation_layers:
            for annotation in layer.annotations:
                annotation_ids.append(annotation.id)
                if not isinstance(annotation.target, ResolvedGenomeTarget):
                    raise ValueError("scene annotations must contain resolved targets")
                target = annotation.target
                if target.particle_indices.size and int(target.particle_indices[-1]) >= particle_count:
                    raise ValueError(f"annotation {annotation.id!r} refers to a particle outside the scene")
                if target.bond_indices.size and int(target.bond_indices[-1]) >= self.bonds.count:
                    raise ValueError(f"annotation {annotation.id!r} refers to a bond outside the scene")
                selected_particles = set(target.particle_indices.tolist())
                for bond_index in target.bond_indices.tolist():
                    left, right = self.bonds.endpoints[bond_index]
                    if int(left) not in selected_particles or int(right) not in selected_particles:
                        raise ValueError(f"annotation {annotation.id!r} bond endpoints are outside its target")
                    owner = bond_owners.get(bond_index)
                    if owner is not None:
                        raise ValueError(
                            f"annotations {owner!r} and {annotation.id!r} share polymer bond {bond_index}"
                        )
                    bond_owners[bond_index] = annotation.id
                if any(isinstance(mark, TubeMark) for mark in annotation.marks) and not target.bond_indices.size:
                    raise ValueError(f"tube annotation {annotation.id!r} resolves to no polymer bonds")
        if len(annotation_ids) != len(set(annotation_ids)):
            raise ValueError("annotation ids must be unique across the scene")
        annotation_layers = tuple(sorted(annotation_layers, key=lambda item: (item.order, item.id)))
        object.__setattr__(self, "annotation_layers", annotation_layers)

        expected_particle_start = 0
        expected_bond_start = 0
        codes = self.particles.chromosome_codes
        genomic_positions = self.particles.genomic_positions
        for group in self.groups:
            if group.count <= 0:
                raise ValueError(f"chromosome group {group.name!r} contains no particles")
            if group.start != expected_particle_start or group.bond_start != expected_bond_start:
                raise ValueError("chromosome particle and bond spans must be contiguous")
            particle_end = group.start + group.count
            bond_end = group.bond_start + group.bond_count
            if particle_end > particle_count or bond_end > self.bonds.count:
                raise ValueError(f"chromosome group {group.name!r} span exceeds its layer")
            if not np.all(codes[group.start:particle_end] == group.code):
                raise ValueError(f"chromosome codes do not match group {group.name!r}")
            group_positions = genomic_positions[group.start:particle_end]
            if len(group_positions) > 1 and np.any(np.diff(group_positions.astype(np.int64)) < 0):
                raise ValueError(f"genomic positions are not sorted in group {group.name!r}")

            endpoints = self.bonds.endpoints[group.bond_start:bond_end]
            if len(endpoints):
                if np.any(endpoints[:, 0] < group.start) or np.any(endpoints[:, 1] >= particle_end):
                    raise ValueError(f"polymer bond crosses chromosome group {group.name!r}")
                if np.any(endpoints[:, 1] != endpoints[:, 0] + 1):
                    raise ValueError("polymer bonds must connect adjacent sorted particles")
                gaps = genomic_positions[endpoints[:, 1]] - genomic_positions[endpoints[:, 0]]
                if np.any(gaps > self.bonds.max_genomic_gap):
                    raise ValueError("polymer bond exceeds max_genomic_gap")
            expected_particle_start = particle_end
            expected_bond_start = bond_end

        if expected_particle_start != particle_count or expected_bond_start != self.bonds.count:
            raise ValueError("chromosome groups do not cover all particles and bonds")

    @property
    def particle_count(self) -> int:
        """Number of particles in the scene."""
        return self.particles.count

    @property
    def bond_count(self) -> int:
        """Number of polymer bonds in the scene."""
        return self.bonds.count
