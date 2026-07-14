"""
Tests and visual examples for ArcTrack.

Output PNGs are saved to tests/output/arc_track/ for visual inspection.
"""

import os
import unittest
from pathlib import Path

from plotly.subplots import make_subplots

from hic_basic.plot.hic import plot_cool_with_tracks
from hic_basic.plot.track import ArcTrack, EigenTrack, GeneTrack, plot_arc_track
from hic_basic.plot.track_composer import LayoutConfig, TrackComposer


DATA_DIR = Path(os.path.dirname(__file__)) / "data"
OUTPUT_DIR = Path(os.path.dirname(__file__)) / "output" / "arc_track"

ARCS_FILE = DATA_DIR / "synthetic_arc_loops.bedpe"
COOL_FILE = DATA_DIR / "visual_test.cool"
EIGS_FILE = DATA_DIR / "visual_test_eigs.bedgraph"
GTF_FILE = DATA_DIR / "visual_test.gtf"
REGION = "chr1:0-5000000"


class TestArcTrack(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    def test_arc_track_load_filters_region(self):
        track = ArcTrack(ARCS_FILE)
        track.load_data(REGION)

        self.assertIsNotNone(track.data)
        self.assertEqual(len(track.data), 5)
        self.assertTrue((track.data["chrom1"] == "chr1").all())
        self.assertTrue((track.data["chrom2"] == "chr1").all())

    def test_arc_track_info_overrides_score_width(self):
        track = ArcTrack(ARCS_FILE)
        track.load_data(REGION)

        fig = make_subplots(rows=1, cols=1)
        track.render(fig, row=1, col=1)

        self.assertEqual(len(fig.data), 5)
        widths = [trace.line.width for trace in fig.data]
        self.assertIn(7, widths)
        self.assertEqual(fig.data[2].line.color, "#d62728")
        self.assertEqual(fig.data[1].line.dash, "dot")
        self.assertEqual(fig.data[4].line.dash, "dashdot")

    def test_arc_track_only_visual(self):
        out_png = OUTPUT_DIR / "arc_track_only.png"

        fig = plot_arc_track(
            ARCS_FILE,
            REGION,
            title="Arc track from synthetic BEDPE",
            width=1000,
            height=260,
            output_path=out_png,
            display_name="loops",
        )

        self.assertEqual(len(fig.data), 5)
        self.assertTrue(out_png.exists())

    def test_arc_track_with_gene_visual(self):
        out_png = OUTPUT_DIR / "arc_gene_composite.png"

        composer = TrackComposer(LayoutConfig(
            rows=2,
            cols=1,
            row_heights=[2, 3],
            vertical_spacing=0.04,
            shared_xaxes=True,
        ))
        composer.add_track(ArcTrack(ARCS_FILE, display_name="loops"), row=1, col=1)
        composer.add_track(GeneTrack(GTF_FILE, mode="gene_only", show_labels=True), row=2, col=1)
        fig = composer.build(region=REGION, title="Arc + gene tracks")
        fig.update_layout(
            width=1000,
            height=420,
            paper_bgcolor="white",
            margin=dict(l=50, r=20, t=55, b=35),
        )
        fig.write_image(str(out_png), width=1000, height=420)

        self.assertTrue(out_png.exists())

    def test_arc_track_with_hic_heatmap_visual(self):
        out_png = OUTPUT_DIR / "arc_hic_heatmap_composite.png"

        tracks = [
            ArcTrack(ARCS_FILE, display_name="loops"),
            GeneTrack(GTF_FILE, mode="gene_only", show_labels=True),
            EigenTrack(EIGS_FILE, eigen_col="E1"),
        ]
        fig = plot_cool_with_tracks(
            coolp1=COOL_FILE,
            tracks=tracks,
            region=REGION,
            track_heights=[3.5, 2.8, 1.6],
            heatmap_height=8,
            square_heatmap=True,
            title="Hi-C heatmap + arc + gene + eigen tracks",
            balance=False,
            cmap="fall",
            ignore_diags=True,
        )
        fig.update_layout(width=1000, height=1050, margin=dict(l=55, r=25, t=60, b=40))
        fig.write_image(str(out_png), width=1000, height=1050)

        self.assertTrue(out_png.exists())


if __name__ == "__main__":
    unittest.main()
