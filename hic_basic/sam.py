"""
SAM/BAM file operations for DNA mapping results.
"""
import gzip
import os
import re
from typing import List, Tuple, Optional, Union

import pysam
import pandas as pd


class Interval:
    """Interval operations for genomic data processing."""

    @staticmethod
    def sort(a: List) -> None:
        """
        Sort intervals based on their coordinates.

        Args:
            a: List of intervals to sort. Each interval can be a number or [start, end] pair.
        """
        if not a:
            return
        if isinstance(a[0], (int, float)):
            # Sort by value if intervals are numbers
            a.sort()
        else:
            # Sort by start position, then end position
            a.sort(key=lambda x: (x[0], x[1]))

    @staticmethod
    def merge(a: List, sorted: bool = True) -> None:
        """
        Merge overlapping intervals.

        Args:
            a: List of intervals to merge, each interval is [start, end]
            sorted: Whether the intervals are already sorted
        """
        if not sorted:
            Interval.sort(a)

        if not a:
            return

        k = 0
        for i in range(1, len(a)):
            # If current interval overlaps with the last merged interval
            if a[k][1] >= a[i][0]:
                # Extend the end of the merged interval if needed
                a[k][1] = max(a[k][1], a[i][1])
            else:
                # Move to next position and copy current interval
                k += 1
                a[k] = a[i][:]

        # Truncate list to contain only merged intervals
        del a[k+1:]

    @staticmethod
    def index_end(a: List, sorted: bool = True) -> None:
        """
        Add index information for efficient overlap finding.

        Args:
            a: List of intervals [start, end, ...]
            sorted: Whether the intervals are already sorted
        """
        if not a:
            return
        if not sorted:
            Interval.sort(a)

        # Initialize the first interval with index 0
        a[0].append(0)

        k = 0
        k_en = a[0][1]

        for i in range(1, len(a)):
            # If current interval starts after the end of the indexed interval
            if k_en <= a[i][0]:
                # Find the next interval that overlaps with current interval
                k += 1
                while k < i:
                    if a[k][1] > a[i][0]:
                        break
                    k += 1
                k_en = a[k][1]

            # Add the index to the current interval
            a[i].append(k)

    @staticmethod
    def find_intv(a: List, x: int) -> int:
        """
        Find the interval containing a given value using binary search.

        Args:
            a: List of intervals
            x: Value to search for

        Returns:
            Index of the interval containing x, or -1 if not found
        """
        left = -1
        right = len(a)

        if not a:
            return -1

        if isinstance(a[0], (int, float)):
            # Handle case where intervals are numbers
            while right - left > 1:
                mid = left + ((right - left) >> 1)  # Equivalent to (right - left) // 2
                if a[mid] > x:
                    right = mid
                elif a[mid] < x:
                    left = mid
                else:
                    return mid
        else:
            # Handle case where intervals are [start, end] pairs
            while right - left > 1:
                mid = left + ((right - left) >> 1)
                if a[mid][0] > x:
                    right = mid
                elif a[mid][0] < x:
                    left = mid
                else:
                    return mid

        return left

    @staticmethod
    def find_ovlp(a: List, st: int, en: int) -> List:
        """
        Find intervals overlapping with a given range.

        Args:
            a: List of intervals [start, end, ...]
            st: Start of the query range
            en: End of the query range

        Returns:
            List of overlapping intervals
        """
        if not a or st >= en:
            return []

        l = Interval.find_intv(a, st)
        k = 0 if l < 0 else a[l][-1]  # Get the index from the last element
        b = []

        for i in range(k, len(a)):
            if a[i][0] >= en:
                break
            elif st < a[i][1]:  # Check if intervals overlap
                b.append(a[i])

        return b


def parse_cigar(cigar: str) -> List[Tuple[int, str]]:
    """
    Parse CIGAR string into operations and lengths.

    Args:
        cigar: CIGAR string (e.g., "100M5D10I")

    Returns:
        List of (length, operation) tuples (e.g., [(100, 'M'), (5, 'D'), (10, 'I')])
    """
    cigar_pattern = re.compile(r'(\d+)([MIDSH])')
    matches = cigar_pattern.findall(cigar)
    return [(int(length), op) for length, op in matches]


def get_file_type(file_path):
    """
    Determine the file type by reading the file content.

    Args:
        file_path (str): Path to the file to check.

    Returns:
        str: The file type ('SAM', 'BAM', 'SAMgz', or 'unknown').
    """
    try:
        # First check file extension (quick hint)
        ext = os.path.splitext(file_path)[1].lower()

        with open(file_path, 'rb') as f:
            # Read file header magic bytes
            magic = f.read(4)
            f.seek(0)

            # 1. Check if BAM file
            # BAM is gzip compressed and starts with BAM\x01 after decompression
            if magic[:2] == b'\x1f\x8b':
                try:
                    with gzip.open(file_path, 'rb') as gz:
                        decompressed_header = gz.read(4)
                        if decompressed_header == b'BAM\x01':
                            return 'BAM'
                except (gzip.BadGzipFile, OSError):
                    # Not a valid gzip file or not BAM
                    pass

            # 2. Check if gzip compressed SAM file
            if magic[:2] == b'\x1f\x8b':
                try:
                    with gzip.open(file_path, 'rb') as gz:
                        # Read enough data to detect SAM features
                        decompressed_data = gz.read(4096)

                        # Try to decode as text
                        try:
                            text = decompressed_data.decode('utf-8', errors='ignore')
                            if is_sam_content(text):
                                return 'SAMgz'
                        except UnicodeDecodeError:
                            # Not valid text
                            pass
                except (gzip.BadGzipFile, OSError):
                    # Not a valid gzip file
                    pass

            # 3. Check if uncompressed SAM file
            # Read file content for detection
            content = f.read(4096)
            try:
                text = content.decode('utf-8', errors='ignore')
                if is_sam_content(text):
                    return 'SAM'
            except UnicodeDecodeError:
                # Not text file
                pass

    except (IOError, OSError) as e:
        print(f"Error reading file: {e}")

    return 'unknown'


def is_sam_content(text, max_lines=20):
    """
    Check if text content matches SAM format characteristics.

    Args:
        text (str): Text content to check
        max_lines (int): Maximum number of lines to check

    Returns:
        bool: True if SAM format, False otherwise
    """
    lines = text.split('\n')

    # Valid SAM header prefixes
    sam_headers = ['@HD', '@SQ', '@RG', '@PG', '@CO']

    for line in lines[:max_lines]:
        line = line.strip()
        if not line:
            continue

        # Check if valid SAM header line
        if line.startswith('@'):
            if any(line.startswith(h) for h in sam_headers):
                return True
        # Check if valid SAM data line
        elif '\t' in line:
            fields = line.split('\t')
            # SAM data line has at least 11 fields
            if len(fields) >= 11:
                return True

    return False


def entry_align_pos(entry):
    """
    Calculate reference and query start/end positions based on CIGAR string.

    Args:
        entry (list): A list representing a SAM format entry where:
                      entry[1] is the FLAG field,
                      entry[3] is the reference start position,
                      entry[5] is the CIGAR string.

    Returns:
        tuple: A tuple containing (rs, re, qs, qe) where:
               rs: Reference start position (0-based)
               re: Reference end position
               qs: Query start position
               qe: Query end position
    """
    cigar_ops = parse_cigar(entry[5])
    clip = [0, 0]  # [left_clip, right_clip]
    x = 0  # Reference coordinate
    y = 0  # Query coordinate
    read_nums = 0  # Bitwise flag for read numbers
    flag = int(entry[1])

    # Determine if on reverse strand
    rev = bool(flag & 16)

    # Parse CIGAR to get ref and query positions
    cigar_ops = parse_cigar(entry[5])
    for length, op in cigar_ops:
        # e.g.,
        #   "100M" -> ref forward 100 bases, query forward 100 bases
        #   "5I"  -> ref forward 0 bases, query forward 5 bases
        if op == 'M':  # Match or mismatch
            x += length
            y += length
        elif op == 'I':  # Insertion
            y += length
        elif op == 'D':  # Deletion
            x += length
        elif op in ['S', 'H']:  # Soft or hard clipping
            # Record clipping lengths and forward query position
            clip[0 if y == 0 else 1] = length  # Left clip if y==0, else right clip
            y += length

    rs = int(entry[3]) - 1  # Reference start position (0-based), directly from SAM
    re = rs + x  # Reference end position, calculated from CIGAR

    # Calculate query start/end based on strand
    qlen = y  # Total aligned length
    if not rev:
        qs = clip[0]  # Query start
        qe = qlen - clip[1]  # Query end
    else:
        qs = clip[1]
        qe = qlen - clip[0]

    read_num = (flag >> 6) & 0x3  # Extract read number (1 or 2)
    if read_num == 3:
        raise ValueError(f"{entry[0]}: incorrect read number flags")

    read_nums |= 1 << read_num

    if read_num == 2:  # Adjust for read 2
        qs1 = qlen - qe
        qe = qlen - qs
        qs = qs1
        rev = not rev

    return rs, re, qs, qe


def convert_alignment_to_entry(alignment):
    """
    Convert a pysam.AlignedSegment object to a list in the same format as a SAM entry.

    Args:
        alignment (pysam.AlignedSegment): An alignment object from pysam.

    Returns:
        list: A list representing the alignment in the same format as a SAM entry.
    """
    # Convert alignment to a list similar to SAM format
    entry = [
        alignment.query_name,  # QNAME (0)
        str(alignment.flag),   # FLAG (1)
        alignment.reference_name if alignment.reference_name else '*',  # RNAME (2)
        str(alignment.reference_start + 1) if alignment.reference_start >= 0 else '0',  # POS (3)
        str(alignment.mapping_quality),  # MAPQ (4)
        alignment.cigarstring if alignment.cigarstring else '*',  # CIGAR (5)
        alignment.next_reference_name if alignment.next_reference_name else '=',  # RNEXT (6)
        str(alignment.next_reference_start + 1) if alignment.next_reference_start >= 0 else '0',  # PNEXT (7)
        str(alignment.template_length) if alignment.template_length else '0',  # TLEN (8)
        alignment.query_sequence if alignment.query_sequence else '*',  # SEQ (9)
        alignment.qual if alignment.qual else '*'  # QUAL (10)
    ]

    # Add optional fields if they exist
    for tag, value in alignment.tags:
        entry.append(f"{tag}:{value}")

    return entry


def count_lines_in_file(file_path):
    """
    Count the number of lines in a file efficiently.

    Args:
        file_path (str): Path to the file.

    Returns:
        int: Number of lines in the file.
    """
    with open(file_path, 'r') as f:
        return sum(1 for _ in f)


def count_reads_per_chrom(
    file_path: str,
    min_mapq: int = 0,
    skip_unmapped: bool = True,
    skip_duplicates: bool = False,
    return_fraction: bool = False
) -> pd.Series:
    """
    Count reads per chromosome from a SAM or BAM file.

    This function efficiently counts the number of reads mapped to each chromosome
    in a SAM or BAM file, with optional filtering by mapping quality and read flags.

    Args:
        file_path (str): Path to the SAM or BAM file.
        min_mapq (int): Minimum mapping quality threshold. Reads with MAPQ < min_mapq
                       will be excluded. Defaults to 0 (no filtering).
        skip_unmapped (bool): If True, skip unmapped reads. Defaults to True.
        skip_duplicates (bool): If True, skip PCR/optical duplicates. Defaults to False.
        return_fraction (bool): If True, return fractions instead of raw counts.
                               Defaults to False.

    Returns:
        pandas.Series: A Series with chromosome names as index and read counts (or fractions)
                      as values, sorted by chromosome name.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the file is not a valid SAM or BAM file.

    Examples:
        Count reads per chromosome:

        >>> counts = count_reads_per_chrom("sample.bam")
        >>> print(counts)
        chr1     1000000
        chr2      900000
        chrX       30000
        chrY        3000
        dtype: int64

        Get fractions instead of counts:

        >>> fractions = count_reads_per_chrom("sample.bam", return_fraction=True)
        >>> print(fractions)
        chr1    0.5128
        chr2    0.4615
        chrX    0.0154
        chrY    0.0015
        dtype: float64

        Filter by mapping quality:

        >>> counts = count_reads_per_chrom("sample.bam", min_mapq=30)
    """
    # Open the file using pysam (auto-detects SAM/BAM format)
    try:
        bam_file = pysam.AlignmentFile(file_path, 'rb' if file_path.endswith('.bam') else 'r')
    except Exception as e:
        raise ValueError(f"Failed to open file '{file_path}': {e}")

    chrom_counts = {}
    total_count = 0

    for read in bam_file.fetch():
        # Skip unmapped reads if requested
        if skip_unmapped and read.is_unmapped:
            continue

        # Skip duplicates if requested
        if skip_duplicates and read.is_duplicate:
            continue

        # Skip reads with low mapping quality
        if read.mapping_quality < min_mapq:
            continue

        # Get chromosome name
        try:
            chrom = bam_file.get_reference_name(read.reference_id)
        except (ValueError, IndexError):
            chrom = '*'  # Unmapped or invalid reference

        # Count the read
        chrom_counts[chrom] = chrom_counts.get(chrom, 0) + 1
        total_count += 1

    bam_file.close()

    # Create pandas Series
    result = pd.Series(chrom_counts, dtype='int64')

    # Convert to fractions if requested
    if return_fraction and total_count > 0:
        result = result / total_count
        result = result.round(6)

    # Sort by chromosome name (natural sort for chromosomes)
    if len(result) > 0:
        result = result.sort_index()

    return result
