"""Adapters from native 3dg data to validated visualization scenes."""

from __future__ import annotations

import colorsys
import hashlib
import os
import re
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
import matplotlib
from matplotlib import colors as mpl_colors

from .scene import ChromosomeGroup, ParticleLayer, PolymerBondLayer, StructureScene


StructureInput = Union[str, os.PathLike, pd.DataFrame]
ColorInput = Union[str, Sequence[Any], np.ndarray]


def _chromosome_key(name: str) -> Tuple[int, int, str]:
    """Return a stable genome-aware sort key for a chromosome name."""
    match = re.match(r"^chr(\d+|X|Y|M|MT)(?:\((mat|pat)\))?$", str(name), flags=re.IGNORECASE)
    if not match:
        return (10_000, 10_000, str(name).lower())
    token, haplotype = match.groups()
    token = token.upper()
    if token.isdigit():
        chromosome = int(token)
    else:
        chromosome = {"X": 23, "Y": 24, "M": 25, "MT": 25}[token]
    haplotype_order = {None: 0, "mat": 1, "pat": 2}[haplotype.lower() if haplotype else None]
    return (chromosome, haplotype_order, str(name).lower())


def _chromosome_color(name: str) -> Tuple[int, int, int]:
    """Return a deterministic territory color for a chromosome name."""
    chromosome, haplotype, _ = _chromosome_key(name)
    if chromosome >= 10_000:
        digest = hashlib.sha256(str(name).encode("utf-8")).digest()
        hue = int.from_bytes(digest[:4], "little") / float(2**32)
    else:
        hue = ((chromosome - 1) * 0.61803398875) % 1.0
    saturation = 0.72 if haplotype in (0, 1) else 0.82
    lightness = 0.63 if haplotype in (0, 1) else 0.50
    rgb = colorsys.hls_to_rgb(hue, lightness, saturation)
    return tuple(int(round(channel * 255)) for channel in rgb)


def structure_scene(
    structure: StructureInput,
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
) -> StructureScene:
    """Build a renderer-independent scene from a 3dg path or native DataFrame.

    Parameters
    ----------
    structure
        A five-column ``.3dg`` path or a DataFrame. DataFrames may use a
        ``(chr, pos)`` MultiIndex with ``x, y, z`` columns, or flat
        ``chrom/chr, pos, x, y, z`` columns.
    color
        ``"chrom"`` for chromosome territory colors, a DataFrame column name,
        an ``N x 3`` RGB array, or a length-N sequence of color specifications.
        Numeric columns use ``cmap``; categorical columns use ``palette`` or a
        deterministic categorical palette.
    palette
        Optional category-to-color mapping. With ``color="chrom"``, keys are
        chromosome names.
    cmap, vmin, vmax, missing_color
        Continuous color normalization and missing-value controls.
    polymer_bonds
        Whether to connect adjacent genome-bin particles into polymer chains.
    max_bond_gap
        Maximum genomic separation in base pairs for an adjacent-bin bond.
        Bonds never cross chromosome or haplotype boundaries.
    dataset_name
        Display name stored in the scene. Defaults to the input file name.
    metadata
        Additional JSON-compatible dataset metadata.

    Returns
    -------
    StructureScene
        Validated positions, colors, genomic identity, and polymer bonds. The
        scene does not contain Three.js or other renderer-specific objects.
    """
    max_bond_gap = int(max_bond_gap)
    if max_bond_gap < 0:
        raise ValueError("max_bond_gap must be non-negative")

    frame, source, default_name, source_checksum = _read_structure(structure)
    input_colors = _resolve_colors(frame, color, palette, cmap, vmin, vmax, missing_color)
    frame = frame.copy()
    frame["_r"] = input_colors[:, 0]
    frame["_g"] = input_colors[:, 1]
    frame["_b"] = input_colors[:, 2]

    chromosome_names = tuple(sorted(frame["chrom"].unique().tolist(), key=_chromosome_key))
    code_dtype = _code_dtype(len(chromosome_names))
    ordered_parts = []
    groups = []
    bond_parts = []
    observed_gaps = []
    particle_offset = 0
    bond_offset = 0

    for code, chromosome in enumerate(chromosome_names):
        part = frame.loc[frame["chrom"] == chromosome].sort_values("pos", kind="stable").reset_index(drop=True)
        coordinates = part[["x", "y", "z"]].to_numpy(dtype=np.float32)
        genomic_positions = part["pos"].to_numpy(dtype=np.int64)
        gaps = np.diff(genomic_positions)
        observed_gaps.append(gaps[gaps > 0])

        if polymer_bonds and len(part) > 1:
            local_left = np.flatnonzero(gaps <= max_bond_gap).astype(np.uint32)
            local_bonds = np.column_stack((local_left, local_left + 1)).astype(np.uint32)
            local_bonds += np.uint32(particle_offset)
        else:
            local_bonds = np.empty((0, 2), dtype=np.uint32)

        group_colors = part[["_r", "_g", "_b"]].to_numpy(dtype=np.uint8)
        representative_color = tuple(np.round(np.median(group_colors, axis=0)).astype(int).tolist())
        minimum = tuple(coordinates.min(axis=0).astype(float).tolist())
        maximum = tuple(coordinates.max(axis=0).astype(float).tolist())
        groups.append(
            ChromosomeGroup(
                name=chromosome,
                code=code,
                start=particle_offset,
                count=len(part),
                bond_start=bond_offset,
                bond_count=len(local_bonds),
                bbox=(minimum, maximum),
                color=representative_color,
            )
        )
        ordered_parts.append(part)
        bond_parts.append(local_bonds)
        particle_offset += len(part)
        bond_offset += len(local_bonds)

    ordered = pd.concat(ordered_parts, ignore_index=True)
    positions = ordered[["x", "y", "z"]].to_numpy(dtype="<f4")
    colors = ordered[["_r", "_g", "_b"]].to_numpy(dtype=np.uint8)
    chromosome_to_code = {name: code for code, name in enumerate(chromosome_names)}
    chromosome_codes = ordered["chrom"].map(chromosome_to_code).to_numpy(dtype=code_dtype)
    genomic_positions = ordered["pos"].to_numpy(dtype="<u4")
    endpoints = np.concatenate(bond_parts, axis=0).astype("<u4", copy=False)

    bbox = (
        tuple(positions.min(axis=0).astype(float).tolist()),
        tuple(positions.max(axis=0).astype(float).tolist()),
    )
    extent = np.subtract(bbox[1], bbox[0])
    base_point_size = float(max(np.linalg.norm(extent) / 620, 0.07))
    resolution = _infer_resolution(observed_gaps)
    scene_checksum = _scene_checksum(
        positions,
        colors,
        chromosome_codes,
        genomic_positions,
        endpoints,
        chromosome_names=chromosome_names,
        max_bond_gap=max_bond_gap,
    )

    scene_metadata = dict(metadata or {})
    scene_metadata.setdefault("source_checksum", source_checksum)
    return StructureScene(
        particles=ParticleLayer(positions, colors, chromosome_codes, genomic_positions),
        bonds=PolymerBondLayer(endpoints, max_bond_gap),
        chromosomes=chromosome_names,
        groups=tuple(groups),
        dataset_name=str(dataset_name or default_name),
        source=source,
        checksum=scene_checksum,
        resolution=resolution,
        bbox=bbox,
        base_point_size=base_point_size,
        metadata=scene_metadata,
    )


def _read_structure(structure: StructureInput) -> Tuple[pd.DataFrame, str, str, str]:
    if isinstance(structure, pd.DataFrame):
        frame = structure.copy()
        source = "<dataframe>"
        name = "3D genome structure"
        source_checksum = ""
    else:
        path = Path(structure).expanduser()
        if not path.is_file():
            raise FileNotFoundError(f"3dg structure does not exist: {path}")
        frame = pd.read_csv(path, sep="\t", comment="#", header=None)
        if frame.shape[1] != 5:
            raise ValueError(f"3dg input must have exactly five tab-separated columns; found {frame.shape[1]}")
        frame.columns = ["chrom", "pos", "x", "y", "z"]
        source = str(path)
        name = path.name
        source_checksum = _sha256_file(path)

    frame = _normalize_frame(frame)
    if not source_checksum:
        source_checksum = _frame_checksum(frame)
    return frame, source, name, source_checksum


def _normalize_frame(frame: pd.DataFrame) -> pd.DataFrame:
    if isinstance(frame.index, pd.MultiIndex) and frame.index.nlevels >= 2:
        index_names = list(frame.index.names)
        reset = frame.reset_index()
        chrom_column = index_names[0] if index_names[0] is not None else reset.columns[0]
        pos_column = index_names[1] if index_names[1] is not None else reset.columns[1]
        frame = reset.rename(columns={chrom_column: "chrom", pos_column: "pos"})
    else:
        frame = frame.reset_index(drop=True)

    lower_columns = {str(column).lower(): column for column in frame.columns}
    chrom_column = next(
        (lower_columns[name] for name in ("chrom", "chr", "chromosome") if name in lower_columns),
        None,
    )
    pos_column = next((lower_columns[name] for name in ("pos", "position") if name in lower_columns), None)
    coordinate_columns = [lower_columns.get(axis) for axis in ("x", "y", "z")]
    if chrom_column is None or pos_column is None or any(column is None for column in coordinate_columns):
        raise ValueError("structure DataFrame needs chrom/chr, pos, x, y, and z fields")

    rename = {chrom_column: "chrom", pos_column: "pos"}
    rename.update({column: axis for column, axis in zip(coordinate_columns, ("x", "y", "z"))})
    frame = frame.rename(columns=rename)
    if frame.empty:
        raise ValueError("structure contains no particles")
    if frame[["chrom", "pos", "x", "y", "z"]].isna().any(axis=None):
        raise ValueError("structure contains missing chromosome, position, or coordinate values")

    frame["chrom"] = frame["chrom"].astype(str)
    if (frame["chrom"].str.len() == 0).any():
        raise ValueError("chromosome names must not be empty")
    numeric = frame[["pos", "x", "y", "z"]].apply(pd.to_numeric, errors="coerce")
    if numeric.isna().any(axis=None):
        raise ValueError("position and coordinate fields must be numeric")
    positions = numeric["pos"].to_numpy(dtype=np.float64)
    if not np.isfinite(numeric.to_numpy(dtype=np.float64)).all():
        raise ValueError("structure contains non-finite positions or coordinates")
    if not np.equal(positions, np.floor(positions)).all():
        raise ValueError("genomic positions must be integers")
    if positions.min() < 0 or positions.max() > np.iinfo(np.uint32).max:
        raise ValueError("genomic positions must fit an unsigned 32-bit integer")
    frame[["pos", "x", "y", "z"]] = numeric
    frame["pos"] = frame["pos"].astype(np.int64)
    return frame


def _resolve_colors(
    frame: pd.DataFrame,
    color: ColorInput,
    palette: Optional[Mapping[Any, Any]],
    cmap: str,
    vmin: Optional[float],
    vmax: Optional[float],
    missing_color: Any,
) -> np.ndarray:
    if isinstance(color, str) and color == "chrom":
        supplied_palette = palette or {}
        resolved = {
            chromosome: _color_to_uint8(supplied_palette.get(chromosome, _chromosome_color(chromosome)))
            for chromosome in frame["chrom"].unique()
        }
        return np.asarray(
            [resolved[chromosome] for chromosome in frame["chrom"]],
            dtype=np.uint8,
        )

    if isinstance(color, str):
        if color not in frame.columns:
            raise ValueError(f"color field not found in structure DataFrame: {color}")
        values = frame[color]
        if pd.api.types.is_numeric_dtype(values):
            return _continuous_colors(values, cmap, vmin, vmax, missing_color)
        return _categorical_colors(values, palette, missing_color)

    array = np.asarray(color)
    if len(array) != len(frame):
        raise ValueError(f"direct colors must have length {len(frame)}; got {len(array)}")
    if array.ndim == 1:
        return np.asarray([_color_to_uint8(value) for value in array], dtype=np.uint8)
    if array.ndim == 2 and array.shape[1] in (3, 4):
        numeric = np.asarray(array[:, :3], dtype=np.float64)
        if not np.isfinite(numeric).all():
            raise ValueError("direct RGB colors contain non-finite values")
        if numeric.min() >= 0 and numeric.max() <= 1:
            numeric = numeric * 255
        if numeric.min() < 0 or numeric.max() > 255:
            raise ValueError("direct RGB colors must be in the 0..1 or 0..255 range")
        return np.round(numeric).astype(np.uint8)
    raise ValueError("color must be 'chrom', a DataFrame column, color sequence, or N x 3 RGB array")


def _continuous_colors(
    values: pd.Series,
    cmap: str,
    vmin: Optional[float],
    vmax: Optional[float],
    missing_color: Any,
) -> np.ndarray:
    numeric = pd.to_numeric(values, errors="coerce").to_numpy(dtype=np.float64)
    finite = np.isfinite(numeric)
    if not finite.any():
        raise ValueError("continuous color field contains no finite values")
    lower = float(np.nanmin(numeric[finite]) if vmin is None else vmin)
    upper = float(np.nanmax(numeric[finite]) if vmax is None else vmax)
    if not np.isfinite(lower) or not np.isfinite(upper) or lower > upper:
        raise ValueError("vmin and vmax must be finite with vmin <= vmax")
    normalized = np.full(len(numeric), 0.5, dtype=np.float64)
    if upper > lower:
        normalized[finite] = np.clip((numeric[finite] - lower) / (upper - lower), 0, 1)
    color_map = matplotlib.colormaps.get_cmap(cmap)
    rgb = np.round(color_map(normalized)[:, :3] * 255).astype(np.uint8)
    rgb[~finite] = _color_to_uint8(missing_color)
    return rgb


def _categorical_colors(values: pd.Series, palette: Optional[Mapping[Any, Any]], missing_color: Any) -> np.ndarray:
    categories = [value for value in pd.unique(values) if not pd.isna(value)]
    resolved = {}
    for index, category in enumerate(sorted(categories, key=lambda value: str(value))):
        if palette is not None and category in palette:
            resolved[category] = _color_to_uint8(palette[category])
        else:
            hue = (index * 0.61803398875) % 1.0
            resolved[category] = _color_to_uint8(colorsys.hls_to_rgb(hue, 0.57, 0.70))
    missing = _color_to_uint8(missing_color)
    return np.asarray([missing if pd.isna(value) else resolved[value] for value in values], dtype=np.uint8)


def _color_to_uint8(value: Any) -> Tuple[int, int, int]:
    if isinstance(value, (tuple, list, np.ndarray)) and len(value) in (3, 4):
        numeric = np.asarray(value[:3], dtype=np.float64)
        if np.isfinite(numeric).all() and numeric.min() >= 0 and numeric.max() <= 255:
            if numeric.max() <= 1:
                numeric = numeric * 255
            return tuple(np.round(numeric).astype(int).tolist())
    try:
        rgb = mpl_colors.to_rgb(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"invalid color specification: {value!r}") from error
    return tuple(int(round(channel * 255)) for channel in rgb)


def _code_dtype(count: int) -> np.dtype:
    if count <= np.iinfo(np.uint8).max + 1:
        return np.dtype("u1")
    if count <= np.iinfo(np.uint16).max + 1:
        return np.dtype("<u2")
    return np.dtype("<u4")


def _infer_resolution(gaps: Sequence[np.ndarray]) -> Optional[int]:
    non_empty = [gap for gap in gaps if len(gap)]
    if not non_empty:
        return None
    values, counts = np.unique(np.concatenate(non_empty), return_counts=True)
    return int(values[np.argmax(counts)])


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _frame_checksum(frame: pd.DataFrame) -> str:
    digest = hashlib.sha256()
    for column in ("chrom", "pos", "x", "y", "z"):
        digest.update(pd.util.hash_pandas_object(frame[column], index=False).to_numpy(dtype="<u8").tobytes())
    return digest.hexdigest()


def _scene_checksum(
    *arrays: np.ndarray,
    chromosome_names: Sequence[str],
    max_bond_gap: int,
) -> str:
    digest = hashlib.sha256()
    digest.update(b"hic_basic.structure_scene.v1\0")
    digest.update(str(int(max_bond_gap)).encode("ascii"))
    for chromosome in chromosome_names:
        encoded = chromosome.encode("utf-8")
        digest.update(len(encoded).to_bytes(4, "little"))
        digest.update(encoded)
    for array in arrays:
        contiguous = np.ascontiguousarray(array)
        digest.update(str(contiguous.dtype).encode("ascii"))
        digest.update(str(contiguous.shape).encode("ascii"))
        digest.update(contiguous.tobytes())
    return digest.hexdigest()
