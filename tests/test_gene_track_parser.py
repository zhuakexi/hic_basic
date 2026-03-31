"""
Unit tests for gene track parser and visualization.
"""
import os
import unittest
from pathlib import Path
import sys

sys.path.insert(0, "/share/home/ychi/dev/hic_basic")

from hic_basic.plot.gene_track_parser import (
    GTFParser,
    GFF3Parser,
    parse_annotation_file,
    filter_genes_in_region,
    Gene,
    Transcript,
    Interval
)


class TestGTFParser(unittest.TestCase):
    """Test GTF parser functionality."""

    def setUp(self):
        self.test_gtf = Path(os.path.dirname(__file__)) / "data" / "synthetic_test.gtf"

    def test_parse_gtf_basic(self):
        """Test basic GTF parsing."""
        parser = GTFParser(self.test_gtf)
        genes = parser.parse()

        self.assertIsInstance(genes, dict)
        self.assertGreater(len(genes), 0)

        # Check gene structure
        self.assertIn("GENE001", genes)
        gene = genes["GENE001"]
        self.assertIsInstance(gene, Gene)
        self.assertEqual(gene.gene_id, "GENE001")
        self.assertEqual(gene.gene_name, "TEST1")
        self.assertEqual(gene.seqid, "chr1")
        self.assertEqual(gene.strand, "+")

    def test_parse_gtf_transcripts(self):
        """Test transcript parsing from GTF."""
        parser = GTFParser(self.test_gtf)
        genes = parser.parse()

        gene = genes["GENE001"]
        self.assertEqual(len(gene.transcripts), 2)

        # Check transcript structure
        tx_ids = [t.transcript_id for t in gene.transcripts]
        self.assertIn("TXN001", tx_ids)
        self.assertIn("TXN002", tx_ids)

    def test_parse_gtf_exons(self):
        """Test exon parsing from GTF."""
        parser = GTFParser(self.test_gtf)
        genes = parser.parse()

        gene = genes["GENE001"]
        tx = next(t for t in gene.transcripts if t.transcript_id == "TXN001")

        self.assertEqual(len(tx.exons), 4)
        # Check exon coordinates
        exon_coords = [(e.start, e.end) for e in tx.exons]
        self.assertIn((1000, 1500), exon_coords)
        self.assertIn((3000, 3500), exon_coords)

    def test_parse_gtf_cds(self):
        """Test CDS parsing from GTF."""
        parser = GTFParser(self.test_gtf)
        genes = parser.parse()

        gene = genes["GENE001"]
        tx = next(t for t in gene.transcripts if t.transcript_id == "TXN001")

        self.assertGreater(len(tx.cds), 0)

    def test_parse_gtf_strand(self):
        """Test strand-specific parsing."""
        parser = GTFParser(self.test_gtf)
        genes = parser.parse()

        # Check minus strand gene
        self.assertIn("GENE002", genes)
        gene = genes["GENE002"]
        self.assertEqual(gene.strand, "-")


class TestFilterGenes(unittest.TestCase):
    """Test gene filtering by region."""

    def setUp(self):
        self.test_gtf = Path(os.path.dirname(__file__)) / "data" / "synthetic_test.gtf"
        parser = GTFParser(self.test_gtf)
        self.genes = parser.parse()

    def test_filter_in_region(self):
        """Test filtering genes within a region."""
        filtered = filter_genes_in_region(self.genes, "chr1", 1000, 10000)

        self.assertIn("GENE001", filtered)
        self.assertNotIn("GENE002", filtered)
        self.assertNotIn("GENE003", filtered)

    def test_filter_with_padding(self):
        """Test filtering with padding."""
        # GENE001 ends at 9000, with padding should include nearby genes
        filtered = filter_genes_in_region(self.genes, "chr1", 1000, 8000, padding=2000)

        self.assertIn("GENE001", filtered)
        # GENE002 starts at 15000, should not be included even with padding
        self.assertNotIn("GENE002", filtered)

    def test_filter_different_chrom(self):
        """Test filtering by chromosome."""
        filtered = filter_genes_in_region(self.genes, "chr2", 1, 100000)
        self.assertEqual(len(filtered), 0)


class TestIntervalOperations(unittest.TestCase):
    """Test interval-related operations."""

    def test_interval_creation(self):
        """Test Interval dataclass."""
        interval = Interval(start=100, end=200, phase="0")
        self.assertEqual(interval.start, 100)
        self.assertEqual(interval.end, 200)
        self.assertEqual(interval.phase, "0")


if __name__ == "__main__":
    unittest.main()
