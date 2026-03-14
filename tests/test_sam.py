"""
Unit tests for hic_basic.sam module.

Tests cover:
- Interval class operations (sort, merge, index_end, find_intv, find_ovlp)
- parse_cigar function
- get_file_type function
- is_sam_content function
- entry_align_pos function
- convert_alignment_to_entry function
- count_lines_in_file function
"""
import unittest
import tempfile
import os
from pathlib import Path

import pysam
import pandas as pd

from hic_basic.sam import (
    Interval,
    parse_cigar,
    get_file_type,
    is_sam_content,
    entry_align_pos,
    convert_alignment_to_entry,
    count_lines_in_file,
    count_reads_per_chrom,
)


class TestIntervalSort(unittest.TestCase):
    """Tests for Interval.sort method."""

    def test_sort_empty_list(self):
        """Test sorting an empty list."""
        a = []
        Interval.sort(a)
        self.assertEqual(a, [])

    def test_sort_numbers(self):
        """Test sorting a list of numbers."""
        a = [5, 2, 8, 1, 9]
        Interval.sort(a)
        self.assertEqual(a, [1, 2, 5, 8, 9])

    def test_sort_intervals(self):
        """Test sorting a list of [start, end] intervals."""
        a = [[50, 60], [10, 20], [30, 40]]
        Interval.sort(a)
        self.assertEqual(a, [[10, 20], [30, 40], [50, 60]])

    def test_sort_intervals_same_start(self):
        """Test sorting intervals with same start position."""
        a = [[10, 30], [10, 20], [10, 40]]
        Interval.sort(a)
        self.assertEqual(a, [[10, 20], [10, 30], [10, 40]])


class TestIntervalMerge(unittest.TestCase):
    """Tests for Interval.merge method."""

    def test_merge_empty_list(self):
        """Test merging an empty list."""
        a = []
        Interval.merge(a)
        self.assertEqual(a, [])

    def test_merge_no_overlap(self):
        """Test merging intervals with no overlaps."""
        a = [[10, 20], [30, 40], [50, 60]]
        Interval.merge(a)
        self.assertEqual(a, [[10, 20], [30, 40], [50, 60]])

    def test_merge_overlapping(self):
        """Test merging overlapping intervals."""
        a = [[10, 30], [20, 40], [50, 60]]
        Interval.merge(a, sorted=True)
        self.assertEqual(a, [[10, 40], [50, 60]])

    def test_merge_contained(self):
        """Test merging when one interval contains another."""
        a = [[10, 50], [20, 30], [60, 70]]
        Interval.merge(a, sorted=True)
        self.assertEqual(a, [[10, 50], [60, 70]])

    def test_merge_unsorted(self):
        """Test merging unsorted intervals."""
        a = [[50, 60], [10, 30], [20, 40]]
        Interval.merge(a, sorted=False)
        self.assertEqual(a, [[10, 40], [50, 60]])


class TestIntervalIndexEnd(unittest.TestCase):
    """Tests for Interval.index_end method."""

    def test_index_end_empty(self):
        """Test indexing an empty list."""
        a = []
        Interval.index_end(a)
        self.assertEqual(a, [])

    def test_index_end_basic(self):
        """Test indexing basic non-overlapping intervals."""
        a = [[10, 20, 'a'], [30, 40, 'b'], [50, 60, 'c']]
        Interval.index_end(a, sorted=True)
        # Each interval should have an index appended
        self.assertEqual(len(a[0]), 4)  # original 3 + index
        self.assertEqual(a[0][-1], 0)
        self.assertEqual(a[1][-1], 1)
        self.assertEqual(a[2][-1], 2)

    def test_index_end_overlapping(self):
        """Test indexing overlapping intervals."""
        a = [[10, 30, 'a'], [20, 40, 'b'], [50, 60, 'c']]
        Interval.index_end(a, sorted=True)
        # First two overlap, so they should share index
        self.assertEqual(a[0][-1], 0)
        self.assertEqual(a[1][-1], 0)  # Overlaps with first
        self.assertEqual(a[2][-1], 2)


class TestIntervalFindIntv(unittest.TestCase):
    """Tests for Interval.find_intv method."""

    def test_find_intv_empty(self):
        """Test finding in an empty list."""
        self.assertEqual(Interval.find_intv([], 5), -1)

    def test_find_intv_found(self):
        """Test finding an existing value."""
        a = [[10, 20], [30, 40], [50, 60]]
        self.assertEqual(Interval.find_intv(a, 35), 1)
        self.assertEqual(Interval.find_intv(a, 15), 0)
        self.assertEqual(Interval.find_intv(a, 55), 2)

    def test_find_intv_not_found(self):
        """Test finding a non-existing value."""
        a = [[10, 20], [30, 40], [50, 60]]
        # Value 25 is between intervals
        result = Interval.find_intv(a, 25)
        self.assertEqual(result, 0)  # Returns the interval before

    def test_find_intv_numbers(self):
        """Test finding in a list of numbers."""
        a = [10, 20, 30, 40, 50]
        self.assertEqual(Interval.find_intv(a, 30), 2)
        self.assertEqual(Interval.find_intv(a, 25), 1)


class TestIntervalFindOvlp(unittest.TestCase):
    """Tests for Interval.find_ovlp method."""

    def test_find_ovlp_empty(self):
        """Test finding overlaps in an empty list."""
        self.assertEqual(Interval.find_ovlp([], 10, 20), [])

    def test_find_ovlp_invalid_range(self):
        """Test with invalid range (start >= end)."""
        a = [[10, 20], [30, 40]]
        self.assertEqual(Interval.find_ovlp(a, 20, 10), [])

    def test_find_ovlp_single(self):
        """Test finding a single overlapping interval."""
        # Need to call index_end first to add index
        a = [[10, 20, 'a'], [30, 40, 'b'], [50, 60, 'c']]
        Interval.index_end(a, sorted=True)
        overlaps = Interval.find_ovlp(a, 35, 38)
        self.assertEqual(len(overlaps), 1)
        self.assertEqual(overlaps[0][2], 'b')

    def test_find_ovlp_multiple(self):
        """Test finding multiple overlapping intervals."""
        # Need to call index_end first to add index
        a = [[10, 25, 'a'], [20, 40, 'b'], [50, 60, 'c']]
        Interval.index_end(a, sorted=True)
        overlaps = Interval.find_ovlp(a, 22, 30)
        self.assertEqual(len(overlaps), 2)

    def test_find_ovlp_none(self):
        """Test finding no overlapping intervals."""
        # Need to call index_end first to add index
        a = [[10, 20, 'a'], [30, 40, 'b'], [50, 60, 'c']]
        Interval.index_end(a, sorted=True)
        overlaps = Interval.find_ovlp(a, 42, 48)
        self.assertEqual(overlaps, [])


class TestParseCigar(unittest.TestCase):
    """Tests for parse_cigar function."""

    def test_parse_cigar_simple(self):
        """Test parsing a simple CIGAR string."""
        self.assertEqual(parse_cigar("100M"), [(100, 'M')])

    def test_parse_cigar_complex(self):
        """Test parsing a complex CIGAR string."""
        cigar = "100M5D10I"
        result = parse_cigar(cigar)
        self.assertEqual(result, [(100, 'M'), (5, 'D'), (10, 'I')])

    def test_parse_cigar_with_clipping(self):
        """Test parsing CIGAR with soft and hard clipping."""
        cigar = "10S90M10S"
        result = parse_cigar(cigar)
        self.assertEqual(result, [(10, 'S'), (90, 'M'), (10, 'S')])

    def test_parse_cigar_empty(self):
        """Test parsing an empty CIGAR string."""
        self.assertEqual(parse_cigar(""), [])

    def test_parse_cigar_all_ops(self):
        """Test parsing CIGAR with all operation types."""
        cigar = "5S100M5I10D5H"
        result = parse_cigar(cigar)
        self.assertEqual(result, [(5, 'S'), (100, 'M'), (5, 'I'), (10, 'D'), (5, 'H')])


class TestGetFileType(unittest.TestCase):
    """Tests for get_file_type function."""

    def test_get_file_type_bam(self):
        """Test detecting BAM file type."""
        bam_path = "tests/data/synthetic_pe.bam"
        if os.path.exists(bam_path):
            file_type = get_file_type(bam_path)
            self.assertEqual(file_type, 'BAM')

    def test_get_file_type_sam(self):
        """Test detecting SAM file type."""
        # Create a temporary SAM file
        sam_content = """@HD	VN:1.0	SO:unsorted
@SQ	SN:chr1	LN:1000
read00001	0	chr1	100	60	100M	*	0	0	ACGT	*
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sam', delete=False) as f:
            f.write(sam_content)
            temp_path = f.name
        try:
            file_type = get_file_type(temp_path)
            self.assertEqual(file_type, 'SAM')
        finally:
            os.unlink(temp_path)

    def test_get_file_type_unknown(self):
        """Test detecting unknown file type."""
        # Create a temporary text file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("This is not a SAM or BAM file")
            temp_path = f.name
        try:
            file_type = get_file_type(temp_path)
            self.assertEqual(file_type, 'unknown')
        finally:
            os.unlink(temp_path)


class TestIsSamContent(unittest.TestCase):
    """Tests for is_sam_content function."""

    def test_is_sam_content_header(self):
        """Test detecting SAM content with header."""
        text = "@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000\n"
        self.assertTrue(is_sam_content(text))

    def test_is_sam_content_data(self):
        """Test detecting SAM content with data line."""
        text = "read001\t0\tchr1\t100\t60\t100M\t*\t0\t0\tACGT\t*\n"
        self.assertTrue(is_sam_content(text))

    def test_is_sam_content_not_sam(self):
        """Test detecting non-SAM content."""
        text = "This is not SAM format\n"
        self.assertFalse(is_sam_content(text))

    def test_is_sam_content_empty(self):
        """Test with empty content."""
        self.assertFalse(is_sam_content(""))

    def test_is_sam_content_insufficient_fields(self):
        """Test with line that has fewer than 11 fields."""
        text = "read001\t0\tchr1\t100\n"
        self.assertFalse(is_sam_content(text))


class TestEntryAlignPos(unittest.TestCase):
    """Tests for entry_align_pos function."""

    def test_entry_align_pos_simple(self):
        """Test with simple 100M alignment."""
        # entry format: [QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, ...]
        entry = ['read1', '0', 'chr1', '100', '60', '100M', '*', '0', '0', 'ACGT', '*']
        rs, re, qs, qe = entry_align_pos(entry)
        self.assertEqual(rs, 99)  # 0-based start
        self.assertEqual(re, 199)  # 99 + 100
        self.assertEqual(qs, 0)
        self.assertEqual(qe, 100)

    def test_entry_align_pos_with_deletion(self):
        """Test with deletion in CIGAR."""
        entry = ['read1', '0', 'chr1', '100', '60', '90M10D10M', '*', '0', '0', 'ACGT', '*']
        rs, re, qs, qe = entry_align_pos(entry)
        self.assertEqual(rs, 99)
        self.assertEqual(re, 99 + 90 + 10 + 10)  # ref: 90M + 10D + 10M = 110
        self.assertEqual(qs, 0)
        self.assertEqual(qe, 90 + 10)  # query: 90M + 10M = 100

    def test_entry_align_pos_with_insertion(self):
        """Test with insertion in CIGAR."""
        entry = ['read1', '0', 'chr1', '100', '60', '50M10I50M', '*', '0', '0', 'ACGT', '*']
        rs, re, qs, qe = entry_align_pos(entry)
        self.assertEqual(rs, 99)
        self.assertEqual(re, 99 + 50 + 50)  # ref: 50M + 50M = 100
        self.assertEqual(qs, 0)
        self.assertEqual(qe, 50 + 10 + 50)  # query: 50M + 10I + 50M = 110

    def test_entry_align_pos_with_soft_clip(self):
        """Test with soft clipping."""
        entry = ['read1', '0', 'chr1', '100', '60', '10S90M', '*', '0', '0', 'ACGT', '*']
        rs, re, qs, qe = entry_align_pos(entry)
        self.assertEqual(rs, 99)
        self.assertEqual(re, 189)  # 99 + 90 (only M consumes reference)
        self.assertEqual(qs, 10)  # Skip soft-clipped bases
        self.assertEqual(qe, 100)  # 10S + 90M = 100

    def test_entry_align_pos_reverse(self):
        """Test with reverse strand alignment."""
        # FLAG 16 = reverse strand
        entry = ['read1', '16', 'chr1', '100', '60', '100M', '*', '0', '0', 'ACGT', '*']
        rs, re, qs, qe = entry_align_pos(entry)
        self.assertEqual(rs, 99)
        self.assertEqual(re, 199)


class TestConvertAlignmentToEntry(unittest.TestCase):
    """Tests for convert_alignment_to_entry function."""

    def test_convert_alignment_basic(self):
        """Test converting a basic alignment."""
        bam_path = "tests/data/synthetic_pe.bam"
        if os.path.exists(bam_path):
            with pysam.AlignmentFile(bam_path, 'rb') as bam:
                for read in bam.fetch():
                    entry = convert_alignment_to_entry(read)
                    # Check entry has at least 11 fields
                    self.assertGreaterEqual(len(entry), 11)
                    # Check QNAME matches
                    self.assertEqual(entry[0], read.query_name)
                    # Check FLAG matches
                    self.assertEqual(entry[1], str(read.flag))
                    break  # Just test the first read

    def test_convert_alignment_roundtrip(self):
        """Test that converted entry can be used to reconstruct key info."""
        bam_path = "tests/data/synthetic_pe.bam"
        if os.path.exists(bam_path):
            with pysam.AlignmentFile(bam_path, 'rb') as bam:
                for read in bam.fetch():
                    entry = convert_alignment_to_entry(read)
                    # Verify we can extract the same chromosome
                    self.assertEqual(entry[2], read.reference_name)
                    # Verify position matches (SAM is 1-based)
                    self.assertEqual(int(entry[3]), read.reference_start + 1)
                    break


class TestCountLinesInFile(unittest.TestCase):
    """Tests for count_lines_in_file function."""

    def test_count_lines_sam(self):
        """Test counting lines in a SAM file."""
        # Create a temporary SAM file
        sam_content = """@HD	VN:1.0
@SQ	SN:chr1	LN:1000
read1	0	chr1	100	60	100M	*	0	0	ACGT	*
read2	0	chr1	200	60	100M	*	0	0	ACGT	*
read3	0	chr1	300	60	100M	*	0	0	ACGT	*
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sam', delete=False) as f:
            f.write(sam_content)
            temp_path = f.name
        try:
            count = count_lines_in_file(temp_path)
            self.assertEqual(count, 5)  # 2 header + 3 data lines
        finally:
            os.unlink(temp_path)

    def test_count_lines_empty(self):
        """Test counting lines in an empty file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            temp_path = f.name
        try:
            count = count_lines_in_file(temp_path)
            self.assertEqual(count, 0)
        finally:
            os.unlink(temp_path)


class TestIntegrationWithSyntheticBam(unittest.TestCase):
    """Integration tests using the synthetic BAM file."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures."""
        cls.bam_path = "tests/data/synthetic_pe.bam"
        if not os.path.exists(cls.bam_path):
            raise unittest.SkipTest(f"Synthetic BAM file not found at {cls.bam_path}")

    def test_bam_file_exists(self):
        """Test that the synthetic BAM file exists."""
        self.assertTrue(os.path.exists(self.bam_path))

    def test_bam_file_has_index(self):
        """Test that the BAM file has an index."""
        index_path = self.bam_path + ".bai"
        self.assertTrue(os.path.exists(index_path))

    def test_get_file_type_on_synthetic_bam(self):
        """Test get_file_type correctly identifies the synthetic BAM."""
        file_type = get_file_type(self.bam_path)
        self.assertEqual(file_type, 'BAM')

    def test_read_all_alignments(self):
        """Test reading all alignments from the synthetic BAM."""
        with pysam.AlignmentFile(self.bam_path, 'rb') as bam:
            count = 0
            chrom_counts = {}
            for read in bam.fetch():
                count += 1
                chrom = bam.get_reference_name(read.reference_id)
                chrom_counts[chrom] = chrom_counts.get(chrom, 0) + 1

            # Verify total count (2000 pairs = 4000 reads)
            self.assertEqual(count, 4000)

            # Verify chromosome distribution
            self.assertIn('chrX', chrom_counts)
            self.assertIn('chrY', chrom_counts)
            # chrX should be ~3% (around 120 reads)
            chrX_pct = chrom_counts.get('chrX', 0) / count * 100
            self.assertAlmostEqual(chrX_pct, 3.0, delta=0.5)
            # chrY should be ~0.3% (around 12 reads)
            chrY_pct = chrom_counts.get('chrY', 0) / count * 100
            self.assertAlmostEqual(chrY_pct, 0.3, delta=0.2)

    def test_convert_and_parse_all_reads(self):
        """Test convert_alignment_to_entry and parse_cigar on all reads."""
        with pysam.AlignmentFile(self.bam_path, 'rb') as bam:
            for read in bam.fetch():
                # Convert to entry
                entry = convert_alignment_to_entry(read)
                # Parse CIGAR from entry
                cigar_ops = parse_cigar(entry[5])
                # Verify CIGAR operations are valid
                for length, op in cigar_ops:
                    self.assertIsInstance(length, int)
                    self.assertIn(op, ['M', 'I', 'D', 'S', 'H', 'N', 'P'])

    def test_entry_align_pos_on_synthetic_reads(self):
        """Test entry_align_pos returns valid positions for synthetic reads."""
        with pysam.AlignmentFile(self.bam_path, 'rb') as bam:
            for read in bam.fetch():
                entry = convert_alignment_to_entry(read)
                rs, re, qs, qe = entry_align_pos(entry)
                # Reference start should be non-negative
                self.assertGreaterEqual(rs, 0)
                # Reference end should be >= start
                self.assertGreaterEqual(re, rs)
                # Query positions should be non-negative
                self.assertGreaterEqual(qs, 0)
                self.assertGreaterEqual(qe, qs)


class TestCountReadsPerChrom(unittest.TestCase):
    """Tests for count_reads_per_chrom function."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures."""
        cls.bam_path = "tests/data/synthetic_pe.bam"
        if not os.path.exists(cls.bam_path):
            raise unittest.SkipTest(f"Synthetic BAM file not found at {cls.bam_path}")

    def test_count_reads_per_chrom_basic(self):
        """Test basic read counting per chromosome."""
        result = count_reads_per_chrom(self.bam_path)

        # Check return type
        self.assertIsInstance(result, pd.Series)

        # Check total count (2000 pairs = 4000 reads)
        self.assertEqual(result.sum(), 4000)

        # Check all expected chromosomes are present
        self.assertIn('chr1', result.index)
        self.assertIn('chr2', result.index)
        self.assertIn('chr3', result.index)
        self.assertIn('chrX', result.index)
        self.assertIn('chrY', result.index)

    def test_count_reads_per_chrom_values(self):
        """Test read counts match expected values."""
        result = count_reads_per_chrom(self.bam_path)

        # Verify chromosome distribution
        # chrX should be ~3% (around 120 reads)
        chrX_pct = result.get('chrX', 0) / result.sum() * 100
        self.assertAlmostEqual(chrX_pct, 3.0, delta=0.5)

        # chrY should be ~0.3% (around 12 reads)
        chrY_pct = result.get('chrY', 0) / result.sum() * 100
        self.assertAlmostEqual(chrY_pct, 0.3, delta=0.2)

    def test_count_reads_per_chrom_return_fraction(self):
        """Test returning fractions instead of counts."""
        result = count_reads_per_chrom(self.bam_path, return_fraction=True)

        # Check return type
        self.assertIsInstance(result, pd.Series)

        # Check that values are fractions (sum should be ~1.0)
        self.assertAlmostEqual(result.sum(), 1.0, places=5)

        # Check that all values are between 0 and 1
        self.assertTrue(all(0 <= v <= 1 for v in result.values))

    def test_count_reads_per_chrom_min_mapq(self):
        """Test filtering by minimum mapping quality."""
        # Without MAPQ filter
        result_all = count_reads_per_chrom(self.bam_path, min_mapq=0)

        # With high MAPQ filter (should still pass since our synthetic BAM has MAPQ=60)
        result_high = count_reads_per_chrom(self.bam_path, min_mapq=50)

        # All reads in synthetic BAM have MAPQ=60, so both should be equal
        self.assertEqual(result_all.sum(), result_high.sum())

    def test_count_reads_per_chrom_skip_unmapped(self):
        """Test skipping unmapped reads."""
        # Our synthetic BAM has all mapped reads, so this shouldn't change counts
        result = count_reads_per_chrom(self.bam_path, skip_unmapped=True)
        result_no_skip = count_reads_per_chrom(self.bam_path, skip_unmapped=False)

        # Should be the same since all reads are mapped
        self.assertEqual(result.sum(), result_no_skip.sum())

    def test_count_reads_per_chrom_sorted_index(self):
        """Test that result is sorted by chromosome name."""
        result = count_reads_per_chrom(self.bam_path)

        # Check that index is sorted
        self.assertEqual(list(result.index), sorted(result.index))

    def test_count_reads_per_chrom_with_sam_file(self):
        """Test counting reads from a SAM file."""
        # Note: Using actual tab characters, not \t escape sequences
        # Sequence must match CIGAR length (100M = 100 bases)
        seq_100 = "A" * 100
        sam_content = "@HD\tVN:1.0\tSO:unsorted\n" \
                      "@SQ\tSN:chr1\tLN:1000\n" \
                      "@SQ\tSN:chr2\tLN:1000\n" \
                      f"read1\t0\tchr1\t100\t60\t100M\t*\t0\t0\t{seq_100}\t*\n" \
                      f"read2\t0\tchr1\t200\t60\t100M\t*\t0\t0\t{seq_100}\t*\n" \
                      f"read3\t0\tchr2\t100\t60\t100M\t*\t0\t0\t{seq_100}\t*\n" \
                      f"read4\t0\tchr2\t200\t60\t100M\t*\t0\t0\t{seq_100}\t*\n" \
                      f"read5\t0\tchr2\t300\t60\t100M\t*\t0\t0\t{seq_100}\t*\n"
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sam', delete=False) as f:
            f.write(sam_content)
            temp_path = f.name
        try:
            result = count_reads_per_chrom(temp_path)
            self.assertEqual(result.get('chr1', 0), 2)
            self.assertEqual(result.get('chr2', 0), 3)
        finally:
            os.unlink(temp_path)


if __name__ == '__main__':
    unittest.main()
