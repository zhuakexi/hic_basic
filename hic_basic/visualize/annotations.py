"""Renderer-independent structure annotation construction and resolution."""

from __future__ import annotations

import hashlib
import json
import re
from dataclasses import replace
from typing import Optional, Sequence, Tuple, Union

import numpy as np

from .scene import (
    AnnotationLayer,
    GenomeRegionTarget,
    ResolvedGenomeTarget,
    StructureAnnotation,
    StructureScene,
    TubeMark,
)


AnnotationInput = Union[StructureAnnotation, Sequence[StructureAnnotation]]


def genomic_segment_annotation(
    region: str,
    *,
    id: str,
    label: Optional[str] = None,
    color: str = "#ff0000",
    tube_width: float = 8.0,
    outline_width: Optional[float] = None,
    outline_color: Optional[str] = None,
    material: str = "unlit",
    opacity: float = 1.0,
    occlusion: str = "depth_tested",
    order: int = 0,
) -> StructureAnnotation:
    """Describe a polymer-segment annotation before scene resolution.

    Parameters
    ----------
    region
        Inclusive genome region such as ``chr1(pat):20000-30000``. A complete
        interval includes both endpoint particles and only real polymer bonds
        already present in the scene.
    id, label
        Stable annotation identifier and optional display label.
    color
        Requested semantic display color. It is preserved in ``unlit`` mode.
    tube_width, outline_width, outline_color
        Independent annotation diameters and optional outline color in
        reference pixels at a 1400-pixel export width.
    material
        ``"unlit"`` bypasses profile color transforms and lighting;
        ``"profile"`` opts into the active base-polymer profile.
    opacity
        Annotation opacity in the interval ``(0, 1]``.
    occlusion
        ``"depth_tested"`` for geometry-aware occlusion or the explicit
        ``"always_visible"`` overlay policy.
    order
        Stable ordering key within the annotation layer.

    Returns
    -------
    StructureAnnotation
        Unresolved typed annotation accepted by :func:`annotate_structure`.
    """
    return StructureAnnotation(
        id=id,
        label=label,
        order=order,
        target=GenomeRegionTarget(region),
        marks=(
            TubeMark(
                color=color,
                tube_width=tube_width,
                outline_width=outline_width,
                outline_color=outline_color,
                material=material,
                opacity=opacity,
                occlusion=occlusion,
            ),
        ),
    )


def annotate_structure(
    scene: StructureScene,
    annotations: AnnotationInput,
    *,
    layer: str = "annotations",
    layer_label: Optional[str] = None,
    visible: bool = True,
    layer_order: int = 0,
) -> StructureScene:
    """Resolve and add typed annotations to an immutable structure scene.

    Genome targets are resolved once in Python against the scene's stable
    chromosome, particle, and polymer-bond order. The resulting scene embeds
    the same resolved membership for interactive HTML and direct PNG output.
    Calling the function again with the same layer id appends annotations to
    that layer; annotation ids remain unique across the complete scene.

    Parameters
    ----------
    scene
        Existing validated structure scene.
    annotations
        One annotation or a sequence of annotations. Genome-region targets may
        be unresolved; already resolved targets are retained and revalidated.
    layer, layer_label, visible, layer_order
        Stable annotation-layer identity, display label, default visibility,
        and composition ordering key.

    Returns
    -------
    StructureScene
        New scene containing resolved, immutable annotation membership and an
        annotation-aware checksum.
    """
    if not isinstance(scene, StructureScene):
        raise TypeError("scene must be a StructureScene")
    if isinstance(annotations, StructureAnnotation):
        candidates = (annotations,)
    else:
        if isinstance(annotations, (str, bytes)):
            raise TypeError("annotations must be StructureAnnotation objects")
        candidates = tuple(annotations)
    if not candidates:
        raise ValueError("annotations must not be empty")
    if any(not isinstance(annotation, StructureAnnotation) for annotation in candidates):
        raise TypeError("annotations must be StructureAnnotation objects")

    resolved = tuple(_resolve_annotation(scene, annotation) for annotation in candidates)
    existing_layers = list(scene.annotation_layers)
    matching = next((item for item in existing_layers if item.id == str(layer).strip()), None)
    if matching is None:
        existing_layers.append(
            AnnotationLayer(
                id=layer,
                label=layer_label,
                visible=visible,
                order=layer_order,
                annotations=resolved,
            )
        )
    else:
        if layer_label is not None and layer_label != matching.label:
            raise ValueError(f"annotation layer {matching.id!r} already has a different label")
        if visible != matching.visible or layer_order != matching.order:
            raise ValueError(f"annotation layer {matching.id!r} already has different layer settings")
        replacement = AnnotationLayer(
            id=matching.id,
            label=matching.label,
            visible=matching.visible,
            order=matching.order,
            annotations=matching.annotations + resolved,
        )
        existing_layers[existing_layers.index(matching)] = replacement

    layers = tuple(sorted(existing_layers, key=lambda item: (item.order, item.id)))
    checksum = _annotated_scene_checksum(scene, layers)
    return replace(
        scene,
        annotation_layers=layers,
        checksum=checksum,
        schema_version=max(scene.schema_version, 2),
    )


def _resolve_annotation(scene: StructureScene, annotation: StructureAnnotation) -> StructureAnnotation:
    if isinstance(annotation.target, ResolvedGenomeTarget):
        return annotation
    if not isinstance(annotation.target, GenomeRegionTarget):
        raise TypeError("unsupported annotation target type")
    selector = annotation.target.selector
    groups, start, end = _parse_region(scene, selector)
    particles = []
    bonds = []
    genomic_positions = scene.particles.genomic_positions
    endpoints = scene.bonds.endpoints

    for group in groups:
        particle_end = group.start + group.count
        local_positions = genomic_positions[group.start:particle_end]
        local_start = int(np.searchsorted(local_positions, start, side="left"))
        local_end = int(np.searchsorted(local_positions, end, side="right"))
        absolute_start = group.start + local_start
        absolute_end = group.start + local_end
        if absolute_end <= absolute_start:
            continue
        particles.extend(range(absolute_start, absolute_end))
        for bond_index in range(group.bond_start, group.bond_start + group.bond_count):
            left, right = endpoints[bond_index]
            if absolute_start <= int(left) < absolute_end and absolute_start <= int(right) < absolute_end:
                bonds.append(bond_index)

    if not particles:
        raise ValueError(f"annotation target {selector!r} resolves to no particles")
    particle_indices = np.asarray(sorted(set(particles)), dtype="<u4")
    bond_indices = np.asarray(sorted(set(bonds)), dtype="<u4")
    coordinates = scene.particles.positions[particle_indices]
    minimum = tuple(coordinates.min(axis=0).astype(float).tolist())
    maximum = tuple(coordinates.max(axis=0).astype(float).tolist())
    anchor = tuple(coordinates.mean(axis=0, dtype=np.float64).astype(float).tolist())
    target = ResolvedGenomeTarget(
        selector=selector,
        chromosomes=tuple(group.name for group in groups),
        particle_indices=particle_indices,
        bond_indices=bond_indices,
        anchor=anchor,
        bbox=(minimum, maximum),
    )
    return replace(annotation, target=target)


def _parse_region(scene: StructureScene, selector: str) -> Tuple[tuple, int, int]:
    value = str(selector).strip()
    colon = value.find(":")
    chromosome_token = value if colon < 0 else value[:colon].strip()
    coordinate_token = "" if colon < 0 else value[colon + 1 :].strip()
    groups = _resolve_groups(scene, chromosome_token)
    start = 0
    end = np.iinfo(np.uint32).max
    if coordinate_token:
        match = re.match(r"^(.+?)\s*-\s*(.+)$", coordinate_token)
        if match is None:
            raise ValueError("annotation region must use chr:start-end for a partial region")
        start = _parse_coordinate(match.group(1))
        end = _parse_coordinate(match.group(2))
        if start > end:
            raise ValueError("annotation region start is greater than its end")
    return groups, start, end


def _resolve_groups(scene: StructureScene, raw_name: str) -> tuple:
    name = str(raw_name).strip()
    if not name:
        raise ValueError("annotation chromosome must not be empty")
    alternatives = [name]
    if not name.lower().startswith("chr"):
        alternatives.append(f"chr{name}")
    by_name = {group.name.lower(): group for group in scene.groups}
    for alternative in alternatives:
        exact = by_name.get(alternative.lower())
        if exact is not None:
            return (exact,)
    for alternative in alternatives:
        prefix = f"{alternative.lower()}("
        homologs = tuple(group for group in scene.groups if group.name.lower().startswith(prefix))
        if homologs:
            return homologs
    raise ValueError(f"annotation chromosome not found: {raw_name}")


def _parse_coordinate(value: str) -> int:
    match = re.match(r"^([0-9][0-9,_.]*(?:\.[0-9]+)?)\s*(bp|kb|mb|gb)?$", value.strip(), re.I)
    if match is None:
        raise ValueError(f"invalid annotation coordinate: {value}")
    numeric = float(match.group(1).replace(",", "").replace("_", ""))
    multiplier = {"bp": 1, "kb": 1_000, "mb": 1_000_000, "gb": 1_000_000_000}[
        (match.group(2) or "bp").lower()
    ]
    result = round(numeric * multiplier)
    if result < 0 or result > np.iinfo(np.uint32).max:
        raise ValueError(f"annotation coordinate out of range: {value}")
    return result


def _annotated_scene_checksum(scene: StructureScene, layers: Sequence[AnnotationLayer]) -> str:
    digest = hashlib.sha256()
    digest.update(b"hic_basic.structure_scene.annotations.v1\0")
    digest.update(str(scene.bonds.max_genomic_gap).encode("ascii"))
    for chromosome in scene.chromosomes:
        encoded = chromosome.encode("utf-8")
        digest.update(len(encoded).to_bytes(4, "little"))
        digest.update(encoded)
    for array in (
        scene.particles.positions,
        scene.particles.colors,
        scene.particles.chromosome_codes,
        scene.particles.genomic_positions,
        scene.bonds.endpoints,
    ):
        contiguous = np.ascontiguousarray(array)
        digest.update(str(contiguous.dtype).encode("ascii"))
        digest.update(str(contiguous.shape).encode("ascii"))
        digest.update(contiguous.tobytes())
    for layer in layers:
        layer_record = {
            "id": layer.id,
            "label": layer.label,
            "visible": layer.visible,
            "order": layer.order,
        }
        digest.update(json.dumps(layer_record, sort_keys=True, separators=(",", ":")).encode("utf-8"))
        for annotation in layer.annotations:
            target = annotation.target
            record = {
                "id": annotation.id,
                "label": annotation.label,
                "order": annotation.order,
                "selector": target.selector,
                "chromosomes": target.chromosomes,
                "anchor": target.anchor,
                "bbox": target.bbox,
                "marks": [
                    {
                        "kind": mark.kind,
                        "requested_color": mark.requested_color,
                        "color": mark.color,
                        "tube_width": mark.tube_width,
                        "outline_width": mark.outline_width,
                        "outline_color": mark.outline_color,
                        "requested_outline_color": mark.requested_outline_color,
                        "material": mark.material,
                        "opacity": mark.opacity,
                        "occlusion": mark.occlusion,
                    }
                    for mark in annotation.marks
                ],
            }
            digest.update(json.dumps(record, sort_keys=True, separators=(",", ":")).encode("utf-8"))
            digest.update(target.particle_indices.tobytes())
            digest.update(target.bond_indices.tobytes())
    return digest.hexdigest()
