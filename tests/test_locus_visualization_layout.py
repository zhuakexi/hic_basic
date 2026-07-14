"""
Regression tests for contact-plus-track locus visualization layout.

This abstracts the EPinformer locus examples: a square Hi-C heatmap,
multi-bin topic-loop anchor marks, guide positions, and an arc track must stay
on the same genomic coordinate system.
"""

import os
import unittest
from pathlib import Path

import cooler
import numpy as np
import pandas as pd

from hic_basic.plot.hic import plot_cool_with_tracks
from hic_basic.plot.track import ArcTrack


class TestLocusVisualizationLayout(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data_dir = Path(os.path.dirname(__file__)) / "data"
        cls.outdir = Path(os.path.dirname(__file__)) / "output" / "locus_layout"
        cls.data_dir.mkdir(parents=True, exist_ok=True)
        cls.outdir.mkdir(parents=True, exist_ok=True)
        cls.cool_path = cls.data_dir / "synthetic_locus_layout.cool"
        cls._create_synthetic_cool(cls.cool_path)

    @staticmethod
    def _create_synthetic_cool(output_path, resolution=5000, n_bins=56):
        chrom = "chr1"
        bins = pd.DataFrame({
            "chrom": [chrom] * n_bins,
            "start": np.arange(0, n_bins * resolution, resolution),
            "end": np.arange(resolution, (n_bins + 1) * resolution, resolution),
        })

        mat = np.zeros((n_bins, n_bins), dtype=int)
        for i in range(n_bins):
            for j in range(i, n_bins):
                dist = j - i
                value = max(1, int(120 * np.exp(-dist / 7)))
                mat[i, j] = value
                mat[j, i] = value
        mat[30:34, 38:42] += 180
        mat[38:42, 30:34] += 180

        upper = np.triu_indices(n_bins)
        cooler.create_cooler(
            str(output_path),
            bins=bins,
            pixels=pd.DataFrame({
                "bin1_id": upper[0],
                "bin2_id": upper[1],
                "count": mat[upper],
            }),
            ordered=True,
            mode="w",
        )

    def test_square_heatmap_mark_and_track_coordinate_contract(self):
        region = "chr1:0-280000"
        tss = 151957
        tss_anchor = (150000, 170000)
        distal_anchor = (190000, 210000)
        distal_center = sum(distal_anchor) / 2

        arc_df = pd.DataFrame([
            {
                "chrom1": "chr1",
                "start1": distal_anchor[0],
                "end1": distal_anchor[1],
                "chrom2": "chr1",
                "start2": tss - 250,
                "end2": tss + 250,
                "score": 1.0,
                "info": {"color": "#d62728", "width": 4},
            }
        ])
        marks = [
            {
                "chrom1": "chr1",
                "start1": tss_anchor[0],
                "end1": tss_anchor[1],
                "chrom2": "chr1",
                "start2": distal_anchor[0],
                "end2": distal_anchor[1],
                "color": "#d62728",
                "shape": "rect",
                "width": 2,
                "label": "topic loop",
            }
        ]

        fig = plot_cool_with_tracks(
            coolp1=str(self.cool_path),
            tracks=[ArcTrack(arc_df, display_name="CRE-TSS arcs")],
            region=region,
            track_heights=[2.0],
            heatmap_height=4.0,
            square_heatmap=True,
            title="locus layout regression",
            balance=False,
            marks=marks,
            colorscale="Viridis",
        )
        out_png = self.outdir / "locus_layout_regression.png"
        fig.write_image(str(out_png), width=760, height=700)

        heatmap = fig.data[0]
        self.assertEqual(float(heatmap.x[0]), 2500.0)
        self.assertEqual(float(heatmap.y[0]), 2500.0)

        mark_shape = fig.layout.shapes[0]
        self.assertEqual(mark_shape.type, "rect")
        self.assertEqual(float(mark_shape.x0), float(tss_anchor[0]))
        self.assertEqual(float(mark_shape.x1), float(tss_anchor[1]))
        self.assertEqual(float(mark_shape.y0), float(distal_anchor[0]))
        self.assertEqual(float(mark_shape.y1), float(distal_anchor[1]))
        self.assertLessEqual(float(mark_shape.x0), tss)
        self.assertGreaterEqual(float(mark_shape.x1), tss)
        self.assertLessEqual(float(mark_shape.y0), distal_center)
        self.assertGreaterEqual(float(mark_shape.y1), distal_center)

        self.assertEqual(fig.layout.yaxis.scaleanchor, "x")
        self.assertEqual(fig.layout.yaxis.scaleratio, 1)
        self.assertEqual(fig.layout.yaxis.constrain, "domain")
        self.assertEqual(fig.layout.xaxis.constrain, "domain")
        self.assertEqual(tuple(fig.layout.xaxis.range), tuple(fig.layout.xaxis2.range))
        self.assertEqual(tuple(fig.layout.xaxis.domain), tuple(fig.layout.xaxis2.domain))
        self.assertTrue(out_png.exists())


if __name__ == "__main__":
    unittest.main()
