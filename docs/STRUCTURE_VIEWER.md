# 3D Genome Viewer and Direct Renderer

`hic_basic.visualize` provides an in-process alternative to the existing PyMOL
and ChimeraX workflows for reviewing `.3dg` structures. It writes one portable
HTML file containing the scene data, CSS, and a pinned Three.js 0.185.1 bundle.
Generating the file does not start a server, invoke a GUI, translate through
mmCIF, use Node.js, or require a GPU.

The same package also provides a direct PNG entry point. It accepts native
structures or a previously generated HTML, applies a versioned render request,
and returns a machine-readable provenance report without interactive steps.

## Quick Start

From the command line:

```bash
./run_in_env.sh hic_basic view3d sample.20k.3dg -o sample.20k.html
```

From Python:

```python
from hic_basic.visualize import write_structure_html

write_structure_html(
    "sample.20k.3dg",
    "sample.20k.html",
    regions=["chr1(mat):0-50Mb", "chr2(pat):10Mb-30Mb"],
    bond_size=1.5,
    background="#151918",
    export_profile="molecular-studio-v1",
)
```

For separate scene construction and output configuration:

```python
from hic_basic.visualize import structure_scene, view_structure

scene = structure_scene(
    structure_df,
    color="score",
    cmap="viridis",
    vmin=0,
    vmax=1,
    max_bond_gap=1_000_000,
)
viewer = view_structure(scene, title="Cell 12", regions=["chr1", "chr8"])
viewer.write_html("cell12.structure.html")
```

`structure_df` may have a `(chr, pos)` MultiIndex and `x, y, z` columns, which
is the native shape returned by common 3dg parsers. A flat DataFrame with
`chrom` or `chr`, `pos`, `x`, `y`, and `z` is also accepted.

## Direct PNG Rendering

The primary programmatic entry point is `render_structure_png`:

```python
from hic_basic.visualize import (
    CameraConfig,
    ExportStyle,
    LayerConfig,
    RasterConfig,
    StructureRenderConfig,
    ViewConfig,
    render_structure_png,
    structure_scene,
)

scene = structure_scene(
    structure_df,
    color="score",
    cmap="viridis",
)
config = StructureRenderConfig(
    view=ViewConfig(
        regions=("chr1(mat):0-50Mb", "chr8(pat)"),
        camera=CameraConfig(mode="fit", fov=34),
    ),
    layers=LayerConfig(show_particles=False, show_bonds=True),
    style=ExportStyle(
        profile="xu2025-v2",
        tube_width=4.5,
        silhouette_width=5.25,
        background="#ffffff",
    ),
    raster=RasterConfig(width=2400, height=1600),
)
result = render_structure_png(scene, "cell12.png", config=config)
print(result.to_dict())
```

Accepted sources are a `.3dg` path, native DataFrame, `StructureScene`,
`StructureViewer`, or a compatible standalone HTML. Use a `StructureScene`
when custom colors or data annotations matter; a bare `.3dg` path uses default
chromosome colors. The HTML path is useful for rendering an archived review
artifact without reconstructing its Python input.

The CLI uses the same API:

```bash
./run_in_env.sh hic_basic render3d sample.20k.3dg \
  -o sample.20k.png --config render.json

./run_in_env.sh hic_basic render3d sample.20k.html \
  -o sample.from-html.png --config render.json
```

`render.json` is a strict, versioned document. Unknown fields and invalid
camera or width values fail before the browser starts:

```json
{
  "schema_version": 1,
  "view": {
    "regions": ["chr1(mat):0-50Mb", "chr8(pat)"],
    "camera": {
      "mode": "explicit",
      "position": [92.0, 61.0, 118.0],
      "target": [0.0, 0.0, 0.0],
      "up": [0.0, 1.0, 0.0],
      "fov": 34.0
    }
  },
  "layers": {
    "show_particles": false,
    "show_bonds": true
  },
  "style": {
    "profile": "molecular-studio-v1",
    "tube_width": 4.5,
    "silhouette_width": 5.0,
    "particle_size": null,
    "background": "#111716"
  },
  "raster": {
    "long_edge": 2400,
    "width": 2400,
    "height": 1600
  }
}
```

Use `camera.mode="fit"` with `position` and `target` set to `null` for a
deterministic fitted view. Explicit `width` and `height` give the strongest
reproducibility. If they are `null`, `long_edge` is combined with the headless
viewport aspect ratio.

## Opening the HTML

The output is a standalone file. On a local machine, open it directly in a
normal browser. No Python process needs to remain running.

For Remote SSH, either download the HTML with the VS Code file explorer and
open it locally, or serve its directory on the remote host:

```bash
./run_in_env.sh python -m http.server 8765 --bind 127.0.0.1 --directory results
```

Forward port `8765` in VS Code and open `http://127.0.0.1:8765/sample.20k.html`
in a local browser. VS Code `Show in Preview` is not the recommended path for
large standalone viewers: its webview security and memory policies can leave
an inline WebGL document at `Loading` even when the same file works in a normal
browser.

## Polymer Bond Contract

A bond represents polymer-chain continuity between genome-bin particles. The
adapter performs a stable sort within each chromosome or haplotype and connects
only adjacent genomic positions whose gap is at most `max_bond_gap`.

- Bonds never cross chromosome or haplotype boundaries.
- Bonds are not inferred from 3D distance.
- Spatial contacts, restraints, and annotation edges are not polymer bonds.
- Setting `polymer_bonds=False` or CLI `--no-bonds` produces a particle-only
  scene.

The default `max_bond_gap` is 1,000,000 bp for compatibility with the existing
hires_utils/mmCIF conversion behavior.

## Color Inputs

`structure_scene(..., color=...)` accepts:

- `"chrom"`: deterministic chromosome territory colors;
- a numeric DataFrame column: continuous `cmap`, `vmin`, and `vmax` mapping;
- a categorical DataFrame column: deterministic categories or an explicit
  `palette` mapping;
- a length-N sequence of Matplotlib color specifications;
- an N x 3 RGB array in either 0..1 or 0..255 range.

Colors remain aligned with particles after genomic sorting and are used by both
the monitor and high-quality PNG profiles. For custom chromosome colors:

```python
scene = structure_scene(
    structure_df,
    color="chrom",
    palette={"chr1(mat)": "#007f8b", "chr1(pat)": "#d1495b"},
)
```

## Viewer and Export

The monitor provides orbit, pan, zoom, particle picking, particle and bond
visibility, particle and bond size, background color, fit-to-selection, and a
union of multiple chromosome or genomic-region selections. An unphased region
such as `chr1:0-20Mb` expands to both homologs when the scene contains
`chr1(mat)` and `chr1(pat)`.

The PNG button and `render_structure_png` call the same versioned browser
protocol and high-quality renderer. The button captures the current camera,
regions, and layer switches; direct rendering obtains them from the config.
Available stable profiles are:

| Profile | Intended use |
| --- | --- |
| `molecular-studio-v1` | Default dark, opaque PBR tube rendering |
| `xu2025-v2` | Preferred light publication rendering |
| `xu2025-v1` | Versioned Xu 2025 baseline compatibility |
| `publication-v1` | Versioned legacy light rendering |

Parameters under an existing profile name are treated as a visual compatibility
contract. New aesthetics require a new versioned profile name. `tube_width`,
`silhouette_width`, `particle_size`, and `background` are explicit overrides;
the returned report records their resolved values. Mesh profiles default to
opaque polymer tubes without embedded particle spheres. Legacy line profiles
default to both bonds and particles for old-result compatibility.

## Rendering Environment

Playwright 1.58.1 and its CPU WebGL runtime libraries are registered in
`envs/dependencies.yml`, rendered into `environment.yml`, and pinned in the
Conda lock. Create and verify environments through the repository workflow:

```bash
make env-create
make env-check
```

Playwright's matching headless Chromium is a separately cached executable. It
is never downloaded by `render_structure_png`. Install it once after creating
or bootstrapping an environment:

```bash
make renderer-bootstrap
```

The default backend forces SwiftShader, so a GPU and display server are not
required. `BrowserBackendConfig` can select an explicit compatible browser,
timeout, viewport, or hardware rendering. Backend choices are deliberately not
part of `StructureRenderConfig`: they describe execution, not the scientific
view. Browser requests to HTTP and HTTPS are blocked during direct rendering;
the standalone HTML and scene data remain offline.

## Architecture and Runtime

```text
.3dg path / native DataFrame
            |
            v
structure_scene -> StructureScene
  particles + genomic identity + polymer bond endpoints
            |                         |
            v                         v
view_structure                 render_structure_png
            |                         ^
            v                         |
standalone HTML ----------------------+
  monitor + versioned render protocol
            |
            v
Playwright + headless Chromium -> validated PNG + RenderResult
```

`StructureScene` contains only NumPy arrays and JSON-compatible metadata. It
does not expose Three.js objects. Arrays are serialized with explicit dtype,
shape, byte order, and checksum descriptors. Writes are atomic.

The 20 kb concept dataset used for acceptance contains 285,395 particles and
285,309 bonds. It was checked in Chromium with SwiftShader software WebGL, so a
physical GPU is not required. The monitor keeps the whole genome to two draw
calls at default bond width. High-quality whole-genome tube export is a separate
operation and is expected to be slower under software rendering.

## Compatibility Boundary

Existing PyMOL/ChimeraX APIs and the legacy `hic_basic render` command remain
unchanged. `hic_basic render3d` is the direct browser renderer. Renderer
protocol version 1 accepts HTML produced by this implementation; older concept
demo HTML files without the protocol fail with an actionable mismatch instead
of relying on private page functions. A server-side path tracer and arbitrary
annotation-edge layers are not part of this API version.
