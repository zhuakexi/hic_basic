import base64
import json
import re
import subprocess
import sys
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from hic_basic.visualize import (
    AnnotationLayer,
    ExportStyle,
    PolymerBondLayer,
    SectionConfig,
    StructureRenderConfig,
    ViewConfig,
    annotate_structure,
    genomic_segment_annotation,
    structure_scene,
    structure_scene_with_feature,
    view_structure,
    write_structure_html,
)


def synthetic_structure():
    return pd.DataFrame(
        [
            ("chr1(pat)", 20_000, 4.0, 0.0, 0.0),
            ("chr1(mat)", 40_000, 2.0, 0.0, 0.0),
            ("chr2", 20_000, 7.0, 0.0, 0.0),
            ("chr1(mat)", 0, 0.0, 0.0, 0.0),
            ("chr1(mat)", 2_000_000, 3.0, 0.0, 0.0),
            ("chr1(pat)", 0, 5.0, 0.0, 0.0),
            ("chr2", 0, 6.0, 0.0, 0.0),
            ("chr1(mat)", 20_000, 1.0, 0.0, 0.0),
        ],
        columns=["chrom", "pos", "x", "y", "z"],
    )


def embedded_payload(document):
    match = re.search(r'<script id="scene-data" type="application/json">(.*?)</script>', document)
    assert match
    return json.loads(match.group(1))


def decode_array(spec):
    dtypes = {
        "float32": "<f4",
        "uint8": "u1",
        "uint16": "<u2",
        "uint32": "<u4",
    }
    return np.frombuffer(base64.b64decode(spec["data"]), dtype=dtypes[spec["dtype"]]).reshape(spec["shape"])


def test_structure_scene_orders_particles_and_builds_only_polymer_bonds():
    frame = synthetic_structure()
    direct_colors = np.column_stack(
        (
            np.arange(len(frame), dtype=np.uint8),
            np.zeros(len(frame), dtype=np.uint8),
            np.zeros(len(frame), dtype=np.uint8),
        )
    )

    scene = structure_scene(frame, color=direct_colors, max_bond_gap=1_000_000)

    assert scene.chromosomes == ("chr1(mat)", "chr1(pat)", "chr2")
    assert scene.particles.genomic_positions.tolist() == [0, 20_000, 40_000, 2_000_000, 0, 20_000, 0, 20_000]
    assert scene.particles.colors[:, 0].tolist() == [3, 7, 1, 4, 5, 0, 6, 2]
    assert scene.bonds.endpoints.tolist() == [[0, 1], [1, 2], [4, 5], [6, 7]]
    assert scene.resolution == 20_000


def test_structure_scene_accepts_native_multiindex_and_field_colors():
    frame = synthetic_structure().set_index(["chrom", "pos"])
    frame["score"] = [0.0, 1.0, np.nan, 0.5, 0.5, 0.25, 0.75, 1.0]

    scene = structure_scene(frame, color="score", cmap="viridis", vmin=0, vmax=1, missing_color="#112233")

    missing_input_index = 2
    sorted_input_indices = [3, 7, 1, 4, 5, 0, 6, 2]
    missing_output_index = sorted_input_indices.index(missing_input_index)
    assert scene.particles.colors[missing_output_index].tolist() == [17, 34, 51]


def test_structure_scene_with_feature_duplicates_haploid_track_and_records_colorbar(tmp_path):
    feature_path = tmp_path / "cpg.tsv"
    pd.DataFrame(
        [
            ("1", 0, 0.005),
            ("1", 20_000, 0.008),
            ("1", 40_000, 0.011),
            ("1", 2_000_000, 0.020),
            ("2", 0, 0.014),
            ("2", 20_000, 0.017),
        ]
    ).to_csv(feature_path, sep="\t", index=False, header=False)

    scene = structure_scene_with_feature(
        synthetic_structure(),
        feature_path,
        feature_name="CpG density",
        vmin=0.005,
        vmax=0.020,
    )

    assert scene.metadata["feature"]["matched_particles"] == 8
    assert scene.metadata["feature"]["missing_particles"] == 0
    assert scene.metadata["color_mapping"]["label"] == "CpG density"
    assert scene.metadata["color_mapping"]["cmap"] == "pymol-cpg"
    assert scene.metadata["color_mapping"]["vmin"] == 0.005
    assert scene.metadata["color_mapping"]["vmax"] == 0.020
    assert scene.particles.colors[0].tolist() == scene.particles.colors[4].tolist()


def test_structure_annotations_resolve_real_bonds_and_reject_overlap():
    base_scene = structure_scene(synthetic_structure(), color="#708087", max_bond_gap=1_000_000)
    base_checksum = base_scene.checksum
    anchor_1 = genomic_segment_annotation(
        "chr1(mat):0-20kb",
        id="anchor-1",
        color="red",
        tube_width=8,
    )
    anchor_2 = genomic_segment_annotation(
        "chr1(mat):20kb-40kb",
        id="anchor-2",
        color="#00d000",
        tube_width=9,
    )
    scene = annotate_structure(
        base_scene,
        [anchor_1, anchor_2],
        layer="anchors",
        layer_label="Loop anchors",
    )
    sequential = annotate_structure(
        annotate_structure(
            base_scene,
            anchor_1,
            layer="anchors",
            layer_label="Loop anchors",
        ),
        anchor_2,
        layer="anchors",
    )

    assert scene.schema_version == 2
    assert scene.checksum != base_checksum
    assert scene.checksum == sequential.checksum
    assert isinstance(scene.annotation_layers[0], AnnotationLayer)
    first, second = scene.annotation_layers[0].annotations
    assert first.target.particle_indices.tolist() == [0, 1]
    assert first.target.bond_indices.tolist() == [0]
    assert second.target.particle_indices.tolist() == [1, 2]
    assert second.target.bond_indices.tolist() == [1]
    assert first.marks[0].requested_color == "red"
    assert first.marks[0].color == "#ff0000"

    with pytest.raises(ValueError, match="share polymer bond"):
        annotate_structure(
            scene,
            genomic_segment_annotation(
                "chr1(mat):0-20kb",
                id="overlap",
                color="#0000ff",
            ),
            layer="more-anchors",
        )

    with pytest.raises(ValueError, match="resolves to no polymer bonds"):
        annotate_structure(
            structure_scene(synthetic_structure(), max_bond_gap=20_000),
            genomic_segment_annotation(
                "chr1(mat):40kb-2Mb",
                id="gap",
                color="#ff0000",
            ),
        )


def test_viewer_writes_standalone_versioned_payload(tmp_path):
    scene = structure_scene(synthetic_structure())
    output = tmp_path / "portable.html"
    viewer = view_structure(
        scene,
        title="A <portable> structure",
        regions=["chr1:0-40kb", "chr2"],
        background="#ffffff",
        bond_size=2.5,
        export_profile="xu2025-v2",
    )

    result = viewer.write_html(output)
    document = result.read_text(encoding="utf-8")
    payload = embedded_payload(document)

    assert result == output
    assert "<title>A &lt;portable&gt; structure</title>" in document
    assert '<script src=' not in document
    assert '<link rel="stylesheet"' not in document
    assert payload["schema_version"] == 1
    assert payload["renderer_protocol_version"] == 2
    assert payload["render_config_schema_version"] == 1
    assert payload["viewer"]["initial_regions"] == ["chr1:0-40kb", "chr2"]
    assert payload["viewer"]["bond_size"] == 2.5
    assert payload["viewer"]["export_profile"] == "xu2025-v2"
    assert payload["viewer"]["render_config"]["view"]["regions"] == ["chr1:0-40kb", "chr2"]
    assert payload["arrays"]["positions"]["dtype"] == "float32"
    assert payload["arrays"]["bond_indices"]["shape"] == [4, 2]
    assert decode_array(payload["arrays"]["bond_indices"]).tolist() == scene.bonds.endpoints.tolist()
    assert "Three.js" in document


def test_viewer_payload_serializes_typed_annotation_layers():
    scene = annotate_structure(
        structure_scene(synthetic_structure(), max_bond_gap=20_000),
        genomic_segment_annotation(
            "chr1(pat):0-20kb",
            id="anchor-pat",
            label="Paternal anchor",
            color="red",
            tube_width=10,
            outline_width=13,
            outline_color="#7f0000",
        ),
        layer="anchors",
    )

    payload = embedded_payload(view_structure(scene).to_html())
    layer = payload["annotations"]["layers"][0]
    annotation = layer["annotations"][0]

    assert payload["annotations"]["schema_version"] == 1
    assert layer["id"] == "anchors"
    assert annotation["target"]["selector"] == "chr1(pat):0-20kb"
    assert annotation["target"]["particle_count"] == 2
    assert annotation["target"]["bond_count"] == 1
    assert decode_array(annotation["target"]["particle_indices"]).tolist() == [4, 5]
    assert decode_array(annotation["target"]["bond_indices"]).tolist() == [2]
    assert annotation["marks"][0]["kind"] == "tube"
    assert annotation["marks"][0]["requested_color"] == "red"
    assert annotation["marks"][0]["color"] == "#ff0000"


def test_typed_render_config_is_embedded_as_authoritative_export_request():
    scene = structure_scene(synthetic_structure())
    config = StructureRenderConfig(
        view=ViewConfig(regions=("chr2",)),
        style=ExportStyle(profile="publication-v1", tube_width=7),
    )

    viewer = view_structure(
        scene,
        regions=["chr1"],
        export_profile="xu2025-v2",
        render_config=config,
    )
    payload = embedded_payload(viewer.to_html())

    assert viewer.regions == ("chr2",)
    assert viewer.export_profile == "publication-v1"
    assert payload["viewer"]["render_config"] == config.to_dict()


def test_section_config_is_embedded_for_monitor_and_export():
    scene = structure_scene(synthetic_structure())
    config = StructureRenderConfig(
        view=ViewConfig(section=SectionConfig(enabled=True, thickness=1.5, offset=-0.25)),
        style=ExportStyle(profile="section-studio-v1", show_colorbar=True),
    )

    payload = embedded_payload(view_structure(scene, render_config=config).to_html())

    assert payload["viewer"]["render_config"]["view"]["section"] == {
        "enabled": True,
        "thickness": 1.5,
        "offset": -0.25,
    }


def test_write_structure_html_convenience_api(tmp_path):
    input_path = tmp_path / "small.3dg"
    synthetic_structure().to_csv(input_path, sep="\t", index=False, header=False)

    result = write_structure_html(input_path, tmp_path / "small.html", max_bond_gap=20_000)

    assert result.is_file()
    assert embedded_payload(result.read_text(encoding="utf-8"))["stats"]["bonds"] == 4


def test_particle_only_scene_keeps_empty_bond_array():
    scene = structure_scene(synthetic_structure(), polymer_bonds=False)
    payload = embedded_payload(view_structure(scene).to_html())

    assert scene.bond_count == 0
    assert scene.bonds.endpoints.shape == (0, 2)
    assert payload["stats"]["bonds"] == 0
    assert payload["arrays"]["bond_indices"]["shape"] == [0, 2]


def test_scene_validation_rejects_unsigned_overflow_and_cross_chromosome_bond():
    with pytest.raises(ValueError, match="non-negative|fit"):
        PolymerBondLayer(np.array([[-1, 0]]), max_genomic_gap=1_000_000)

    scene = structure_scene(synthetic_structure())
    endpoints = scene.bonds.endpoints.copy()
    endpoints[1] = [2, 4]
    with pytest.raises(ValueError, match="crosses chromosome"):
        replace(scene, bonds=PolymerBondLayer(endpoints, scene.bonds.max_genomic_gap))


def test_view3d_cli_writes_json_report(tmp_path):
    input_path = tmp_path / "small.3dg"
    output_path = tmp_path / "small.viewer.html"
    synthetic_structure().to_csv(input_path, sep="\t", index=False, header=False)

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "hic_basic",
            "view3d",
            str(input_path),
            "-o",
            str(output_path),
            "--region",
            "chr1:0-40kb",
            "--profile",
            "xu2025-v2",
            "--force",
            "--json",
        ],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    report = json.loads(result.stdout)
    assert report["output"] == str(output_path)
    assert report["particles"] == 8
    assert report["bonds"] == 4
    assert output_path.is_file()
