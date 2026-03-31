"""
Functional tests for gene track visualization.

This module tests the plot_gene_track function with synthetic data
and generates PNG output for visual inspection.
"""
import os
import unittest
from pathlib import Path
import sys

sys.path.insert(0, "/share/home/ychi/dev/hic_basic")

from hic_basic.plot.track import plot_gene_track


class TestGeneTrackVisualization(unittest.TestCase):
    """Test gene track visualization functionality."""

    def setUp(self):
        """Set up test fixtures."""
        self.test_gtf = Path(os.path.dirname(__file__)) / "data" / "synthetic_test.gtf"
        # Real data from GENCODE M25
        self.real_gtf = Path(os.path.dirname(__file__)) / "data" / "gencode_m38_chr1_Xkr4.gtf"
        self.outdir = Path(os.path.dirname(__file__)) / "output" / "gene_track"
        self.outdir.mkdir(parents=True, exist_ok=True)

    def test_plot_gene_track_basic(self):
        """
        Unit test for basic gene track plotting.
        Tests that the function returns a valid figure object.
        """
        fig = plot_gene_track(
            annotation_file=str(self.test_gtf),
            region="chr1:1000-9000",
            mode="all_transcripts",
            title="Test: All Transcripts"
        )

        self.assertIsNotNone(fig)
        # Should have shapes for transcript structures
        self.assertGreater(len(fig.layout.shapes), 0)
        # Title should be set
        self.assertEqual(fig.layout.title.text, "Test: All Transcripts")

    def test_plot_gene_track_gene_only_mode(self):
        """
        Unit test for gene_only mode.
        """
        fig = plot_gene_track(
            annotation_file=str(self.test_gtf),
            region="chr1:1000-20000",
            mode="gene_only",
            title="Test: Gene Only Mode"
        )

        self.assertIsNotNone(fig)
        self.assertGreater(len(fig.layout.shapes), 0)

    def test_plot_gene_track_minus_strand(self):
        """
        Unit test for minus strand gene visualization.
        """
        fig = plot_gene_track(
            annotation_file=str(self.test_gtf),
            region="chr1:15000-20000",
            mode="all_transcripts",
            title="Test: Minus Strand"
        )

        self.assertIsNotNone(fig)
        self.assertGreater(len(fig.layout.shapes), 0)

    def test_plot_gene_track_single_exon(self):
        """
        Unit test for single-exon gene visualization.
        """
        fig = plot_gene_track(
            annotation_file=str(self.test_gtf),
            region="chr1:25000-26000",
            mode="all_transcripts",
            title="Test: Single Exon Gene"
        )

        self.assertIsNotNone(fig)
        self.assertGreater(len(fig.layout.shapes), 0)

    def test_plot_gene_track_no_genes(self):
        """
        Test behavior when no genes are found in region.
        """
        fig = plot_gene_track(
            annotation_file=str(self.test_gtf),
            region="chr1:100000-200000",
            title="Test: No Genes"
        )

        self.assertIsNotNone(fig)
        # Should have annotation about no genes found
        self.assertGreater(len(fig.layout.annotations), 0)

    def test_plot_gene_track_visual_all_transcripts(self):
        """
        Visual test for all_transcripts mode.
        Generates PNG output for manual inspection.
        """
        out_png = self.outdir / "gene_track_all_transcripts.png"

        fig = plot_gene_track(
            annotation_file=str(self.test_gtf),
            region="chr1:1000-9000",
            mode="all_transcripts",
            title="Visual Test: All Transcripts (GENE001)",
            height=400,
            width=1200,
            show_labels=True
        )

        fig.write_image(str(out_png), width=1200, height=500)
        self.assertTrue(out_png.exists())
        print(f"Visual test output: {out_png}")

    def test_plot_gene_track_visual_gene_only(self):
        """
        Visual test for gene_only mode.
        Generates PNG output for manual inspection.
        """
        out_png = self.outdir / "gene_track_gene_only.png"

        fig = plot_gene_track(
            annotation_file=str(self.test_gtf),
            region="chr1:1000-20000",
            mode="gene_only",
            title="Visual Test: Gene Only Mode",
            height=300,
            width=1200
        )

        fig.write_image(str(out_png), width=1200, height=400)
        self.assertTrue(out_png.exists())
        print(f"Visual test output: {out_png}")

    def test_plot_gene_track_visual_full_region(self):
        """
        Visual test for full region with multiple genes.
        Generates PNG output for manual inspection.
        """
        out_png = self.outdir / "gene_track_full_region.png"

        fig = plot_gene_track(
            annotation_file=str(self.test_gtf),
            region="chr1:1-30000",
            mode="all_transcripts",
            title="Visual Test: Full Region (chr1:1-30000)",
            height=500,
            width=1400
        )

        fig.write_image(str(out_png), width=1400, height=600)
        self.assertTrue(out_png.exists())
        print(f"Visual test output: {out_png}")

    def test_plot_gene_track_visual_custom_heights(self):
        """
        Visual test with custom feature heights.
        Generates PNG output for manual inspection.
        """
        out_png = self.outdir / "gene_track_custom_heights.png"

        fig = plot_gene_track(
            annotation_file=str(self.test_gtf),
            region="chr1:1000-9000",
            mode="all_transcripts",
            title="Visual Test: Custom Heights",
            height=400,
            width=1200,
            utr_height=0.25,
            exon_height=0.35,
            cds_height=0.70,
            row_spacing=0.3
        )

        fig.write_image(str(out_png), width=1200, height=500)
        self.assertTrue(out_png.exists())
        print(f"Visual test output: {out_png}")

    def test_plot_gene_track_real_data_xkr4(self):
        """
        Visual test with real GENCODE annotation data (Xkr4 gene).
        Generates PNG output for manual inspection.
        
        Data source: GENCODE M25 (Ensembl 100)
        Gene: Xkr4 (ENSMUSG00000051951.5)
        Region: chr1:3205901-3671498
        Strand: -
        Transcripts: 3 (Xkr4-201, Xkr4-202, Xkr4-203)
        """
        out_png = self.outdir / "gene_track_real_Xkr4.png"

        fig = plot_gene_track(
            annotation_file=str(self.real_gtf),
            region="chr1:3205901-3680000",
            mode="all_transcripts",
            title="Real Data: Xkr4 (chr1:-)",
            height=400,
            width=1400,
            show_labels=True,
            label_font_size=11
        )

        fig.write_image(str(out_png), width=1400, height=500)
        self.assertTrue(out_png.exists())
        print(f"Visual test output: {out_png}")

    def test_plot_gene_track_real_data_parser(self):
        """
        Unit test for parsing real GENCODE annotation data.
        Tests that the parser correctly handles real GTF format.
        """
        from hic_basic.plot.gene_track_parser import parse_annotation_file, filter_genes_in_region
        
        genes = parse_annotation_file(self.real_gtf)
        
        # Check gene was parsed
        self.assertIn("ENSMUSG00000051951.5", genes)
        gene = genes["ENSMUSG00000051951.5"]
        
        # Check gene attributes
        self.assertEqual(gene.gene_name, "Xkr4")
        self.assertEqual(gene.seqid, "chr1")
        self.assertEqual(gene.strand, "-")
        self.assertEqual(gene.biotype, "protein_coding")
        
        # Check transcripts
        self.assertEqual(len(gene.transcripts), 3)
        tx_ids = [t.transcript_id for t in gene.transcripts]
        self.assertIn("ENSMUST00000070533.4", tx_ids)  # Xkr4-201
        self.assertIn("ENSMUST00000159265.1", tx_ids)  # Xkr4-202
        self.assertIn("ENSMUST00000162897.1", tx_ids)  # Xkr4-203
        
        # Check main transcript (Xkr4-201) has CDS and UTR
        tx_201 = next(t for t in gene.transcripts if t.transcript_id == "ENSMUST00000070533.4")
        self.assertGreater(len(tx_201.exons), 0)
        self.assertGreater(len(tx_201.cds), 0)
        # UTR may be parsed as utr or utr5/utr3
        has_utr = len(tx_201.utr) > 0 or len(tx_201.utr5) > 0 or len(tx_201.utr3) > 0
        self.assertTrue(has_utr, "Xkr4-201 should have UTR regions")


if __name__ == "__main__":
    unittest.main()
