"""Validated, renderer-independent scene objects for 3D genome structures."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping, Optional, Tuple

import numpy as np


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
