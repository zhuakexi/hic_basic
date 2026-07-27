import json
from pathlib import Path

import numpy as np
import pytest
from PIL import Image

from hic_basic.visualize import (
    BrowserBackendConfig,
    CameraConfig,
    ExportStyle,
    LayerConfig,
    RasterConfig,
    StructureRenderConfig,
    ViewConfig,
    load_render_config,
    render_structure_png,
    structure_scene,
    view_structure,
)


DATA = Path(__file__).parent / "data" / "representative_structure.3dg"


def test_render_config_round_trip_and_strict_validation(tmp_path):
    config = StructureRenderConfig(
        view=ViewConfig(
            regions=("chr1(mat):0-40kb", "chr2"),
            camera=CameraConfig(
                mode="explicit",
                position=(5, 4, 6),
                target=(0, 0, 0),
                fov=38,
            ),
        ),
        layers=LayerConfig(show_particles=True, show_bonds=True),
        style=ExportStyle(
            profile="xu2025-v2",
            tube_width=5.5,
            silhouette_width=6.25,
            particle_size=4,
            background="#fefefe",
        ),
        raster=RasterConfig(width=720, height=480),
    )
    path = tmp_path / "render.json"
    path.write_text(config.to_json(), encoding="utf-8")

    assert load_render_config(path) == config
    assert load_render_config(config.to_json()) == config
    assert config.to_dict()["style"]["background"] == "#fefefe"
    with pytest.raises(ValueError, match="unknown structure render config fields"):
        load_render_config({"unknown": True})
    with pytest.raises(ValueError, match="tube_width"):
        ExportStyle(tube_width=0)
    with pytest.raises(ValueError, match="must be an integer"):
        RasterConfig(width=640.5, height=480)
    with pytest.raises(ValueError, match="position and target"):
        CameraConfig(mode="explicit", position=(1, 2, 3), target=(1, 2, 3))
    with pytest.raises(ValueError, match="parallel"):
        CameraConfig(mode="explicit", position=(0, 0, 3), target=(0, 0, 0), up=(0, 0, 1))


def _render_or_skip(source, output, config):
    try:
        return render_structure_png(
            source,
            output,
            config=config,
            backend=BrowserBackendConfig(timeout_seconds=60),
        )
    except RuntimeError as error:
        if "Playwright is missing" in str(error) or "Chromium asset is unavailable" in str(error):
            pytest.skip(str(error))
        raise


def test_direct_render_applies_width_regions_layers_and_reports_backend(tmp_path):
    scene = structure_scene(DATA)
    thin_config = StructureRenderConfig(
        view=ViewConfig(regions=("chr1",)),
        layers=LayerConfig(show_particles=False, show_bonds=True),
        style=ExportStyle(profile="xu2025-v2", tube_width=2, silhouette_width=2),
        raster=RasterConfig(width=640, height=400),
    )
    thick_config = StructureRenderConfig(
        view=thin_config.view,
        layers=thin_config.layers,
        style=ExportStyle(profile="xu2025-v2", tube_width=10, silhouette_width=10),
        raster=thin_config.raster,
    )

    thin = _render_or_skip(scene, tmp_path / "thin.png", thin_config)
    thick = _render_or_skip(scene, tmp_path / "thick.png", thick_config)

    assert (thin.width, thin.height) == (640, 400)
    assert thin.scene_sha256 == scene.checksum
    assert thin.metadata["selection"] == {"regions": ["chr1"], "particles": 10, "bonds": 8}
    assert thin.metadata["request"]["layers"] == {"show_particles": False, "show_bonds": True}
    assert thin.metadata["request"]["style"]["tube_width"] == 2
    assert thin.metadata["renderer"]["software"] is True
    assert thin.metadata["backend"]["protocol_version"] == 1

    thin_pixels = np.asarray(Image.open(thin.output).convert("RGB"))
    thick_pixels = np.asarray(Image.open(thick.output).convert("RGB"))
    thin_ink = np.any(thin_pixels < 245, axis=2).sum()
    thick_ink = np.any(thick_pixels < 245, axis=2).sum()
    assert thick_ink > thin_ink * 1.8


def test_direct_render_accepts_standalone_html_and_explicit_camera(tmp_path):
    scene = structure_scene(DATA)
    html_path = view_structure(scene).write_html(tmp_path / "scene.html")
    config = StructureRenderConfig(
        view=ViewConfig(
            regions=("chr2",),
            camera=CameraConfig(
                mode="explicit",
                position=(5, 3, 7),
                target=(0, 0, 0),
                fov=42,
            ),
        ),
        layers=LayerConfig(show_particles=True, show_bonds=False),
        style=ExportStyle(profile="publication-v1", particle_size=12, background="white"),
        raster=RasterConfig(width=480, height=320),
    )

    result = _render_or_skip(html_path, tmp_path / "html-source.png", config)

    assert result.output.is_file()
    assert result.metadata["selection"]["regions"] == ["chr2"]
    assert result.metadata["selection"]["particles"] == 2
    assert result.metadata["selection"]["bonds"] == 1
    assert result.metadata["camera"]["mode"] == "explicit"
    assert result.metadata["camera"]["position"] == [5, 3, 7]
    assert result.metadata["request"]["style"]["background"] == "#ffffff"
    assert result.metadata["render_passes"] == 1


def test_render_result_is_json_serializable(tmp_path):
    config = StructureRenderConfig(raster=RasterConfig(width=320, height=240))
    result = _render_or_skip(DATA, tmp_path / "result.png", config)

    report = result.to_dict()
    assert json.loads(json.dumps(report))["output"] == str(result.output)
    assert report["size_bytes"] == result.output.stat().st_size
