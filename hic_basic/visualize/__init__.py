"""Native, browser-based visualization for 3D genome structures."""

from .config import (
    BrowserBackendConfig,
    CameraConfig,
    ExportStyle,
    LayerConfig,
    RasterConfig,
    SectionConfig,
    StructureRenderConfig,
    ViewConfig,
    load_render_config,
)
from .html import StructureViewer, view_structure, write_structure_html
from .render import RenderResult, render_structure_png
from .scene import (
    AnnotationLayer,
    ChromosomeGroup,
    GenomeRegionTarget,
    ParticleLayer,
    PolymerBondLayer,
    ResolvedGenomeTarget,
    StructureAnnotation,
    StructureScene,
    TubeMark,
)
from .annotations import annotate_structure, genomic_segment_annotation
from .structure import structure_scene, structure_scene_with_feature

__all__ = [
    "ChromosomeGroup",
    "AnnotationLayer",
    "BrowserBackendConfig",
    "CameraConfig",
    "ExportStyle",
    "LayerConfig",
    "GenomeRegionTarget",
    "ParticleLayer",
    "PolymerBondLayer",
    "RasterConfig",
    "SectionConfig",
    "RenderResult",
    "ResolvedGenomeTarget",
    "StructureAnnotation",
    "StructureScene",
    "StructureRenderConfig",
    "StructureViewer",
    "TubeMark",
    "ViewConfig",
    "load_render_config",
    "annotate_structure",
    "genomic_segment_annotation",
    "render_structure_png",
    "structure_scene",
    "structure_scene_with_feature",
    "view_structure",
    "write_structure_html",
]
