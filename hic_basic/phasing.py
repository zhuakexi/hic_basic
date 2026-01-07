"""
Sample sex and phasing.
"""
import gzip
import os
import re
import sys
import struct
import pandas as pd
import pysam
from typing import List, Tuple

from tqdm import tqdm
import numpy as np
from .calculate import mt
from .hicio import parse_seg, parse_sam_line
from .genome import GenomeIdeograph


### --- count alleles of vcf sites from sam --- ###
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
def parse_phased_snp(phased_snp_file: str) -> dict:
    snp_dict = {}
    with open(phased_snp_file, 'r') as f:
        for line in f:
            t = line.strip().split("\t")
            if len(t) < 4 or len(t[2]) != 1 or len(t[3]) != 1:
                continue
            if t[0] not in snp_dict:
                snp_dict[t[0]] = []
            pos = int(t[1])
            snp_dict[t[0]].append([pos - 1, pos, t[2], t[3]])
    
    for c in snp_dict:
        Interval.index_end(snp_dict[c], True)
    return snp_dict

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
        str: The file type ('SAM', 'BAM', or 'unknown').
    """
    try:
        with open(file_path, 'rb') as f:
            # 先检查是否是 gzip 压缩文件（BAM 就是 gzip）
            magic = f.read(2)
            
            if magic == b'\x1f\x8b':
                # 这是一个 gzip 文件，可能是 BAM
                try:
                    # 使用 gzip 打开并读取前 4 个解压后的字节
                    with gzip.open(file_path, 'rb') as gz:
                        decompressed_header = gz.read(4)
                        if decompressed_header == b'BAM\x01':
                            return 'BAM'
                except (gzip.BadGzipFile, OSError):
                    # 不是有效的 gzip 文件，或者不是 BAM
                    pass
            
            # 如果不是 gzip，或者解压后不是 BAM，检查 SAM
            f.seek(0)
            # 读取足够的内容来检查
            content = f.read(4096)
            
            try:
                text = content.decode('utf-8', errors='ignore')
                lines = text.split('\n')
                
                # 检查 SAM 头行
                sam_headers = ['@HD', '@SQ', '@RG', '@PG', '@CO']
                
                for line in lines[:10]:  # 检查前10行
                    line = line.strip()
                    if not line:
                        continue
                    
                    # 检查是否是有效的 SAM 头行
                    if line.startswith('@'):
                        if any(line.startswith(h) for h in sam_headers):
                            return 'SAM'
                    # 检查是否是有效的 SAM 数据行
                    elif '\t' in line:
                        fields = line.split('\t')
                        # SAM 数据行至少有11个字段
                        if len(fields) >= 11:
                            return 'SAM'
                            
            except UnicodeDecodeError:
                # 不是文本文件，排除 SAM
                pass
                
    except (IOError, OSError) as e:
        print(f"Error reading file: {e}")
    
    return 'unknown'

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
        raise ValueError(f"{t[0]}: incorrect read number flags")
    
    read_nums |= 1 << read_num
    
    if read_num == 2:  # Adjust for read 2
        qs1 = qlen - qe
        qe = qlen - qs
        qs = qs1
        rev = not rev

    return rs, re, qs, qe


def check_allele(entry, rs, v, append_features=None, min_baseq=20, verbose=0):
    """
    Check if a SAM entry contains a specific allele at a given position.

    Args:
        entry (list): A list representing a SAM format entry where:
                      entry[2] is the reference sequence name,
                      entry[5] is the CIGAR string,
                      entry[9] is the sequence,
                      entry[10] is the base quality string.
        rs (int): Reference start position (0-based).
        v (list): A list containing SNP information [position, snp_id, ref_allele, alt_allele].
        append_features (list, optional): List of additional features to append to the result.
        min_baseq (int, optional): Minimum base quality threshold. Defaults to 20.
        verbose (int, optional): Verbosity level. Defaults to 0.

    Returns:
        list or None: A list containing [chrom, pos, allele, read_name, relative_pos] and additional features if the allele matches,
                      otherwise None.
    """
    cigar_ops = parse_cigar(entry[5])
    chrom = entry[2]

    x_pos = rs
    y_pos = 0
    res = None
    for length, op in cigar_ops:
        if op == "M": # Match or mismatch
            p = v[0] # SNP ref position
            if x_pos <= p < x_pos + length:
                # SNP falls within this cigar operation
                p_adj = p - x_pos + y_pos # relative position of SNP site in read
                if p_adj < 0 or p_adj >= len(entry[9]):
                    # SNP relative position out of read sequence range
                    raise ValueError("CIGAR parsing error")
                c = entry[9][p_adj] # read base at SNP relative position
                q = ord(entry[10][p_adj]) - 33 if len(entry[10]) == len(entry[9]) else min_baseq

                if append_features is not None and len(append_features) > 0:
                    # parse entry
                    entry_sr = parse_sam_line(entry)
                    append_feature_values = entry_sr[append_features].tolist()
                else:
                    append_feature_values = []

                if q >= min_baseq:
                    # Determine phase based on allele
                    if c == v[2]:
                        res = [chrom, v[1], v[2], entry[0], p_adj] # snp chrom, snp pos, ref allele, read name, relative pos
                    elif c == v[3]:
                        res = [chrom, v[1], v[3], entry[0], p_adj] # snp chrom, snp pos, alt allele, read name, relative pos
                    else:
                        res = []
                        if verbose >= 1:
                            # don't contribute to phase if base doesn't match either allele
                            print(f'WARNING: read {entry[0]} has base {c} at position {entry[2]}:{v[1]} (not {v[2]}/{v[3]})', file=sys.stderr)
            x_pos += length
            y_pos += length
        elif op == "I":
            y_pos += length
        elif op == "D":
            x_pos += length
        elif op == "S":
            y_pos += length
    if res is not None:
        res = res + append_feature_values
    return res


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

def sam_mark_alleles(sam_file, phased_snp_file, outfile=None, append_features=None, min_mapq=20, min_baseq=20, show_progress=False, verbose=0):
    """
    Mark alleles in SAM/BAM file based on phased SNP information.

    Args:
        sam_file (str): Path to the input SAM or BAM file.
        phased_snp_file (str): Path to the phased SNP file.
        outfile (str, optional): Path to the output file. If None, returns a DataFrame. Defaults to None.
        append_features (list, optional): List of additional features to append to the output. Defaults to None.
        min_mapq (int, optional): Minimum mapping quality threshold. Defaults to 20.
        min_baseq (int, optional): Minimum base quality threshold. Defaults to 20.
        show_progress (bool, optional): Whether to show a progress bar. Defaults to False.
        verbose (int, optional): Verbosity level. Defaults to 0.

    Returns:
        str or pandas.DataFrame: Path to the output file if outfile is specified, otherwise a DataFrame containing marked alleles.
    """
    # Determine file type
    file_type = get_file_type(sam_file)
    if verbose >= 1:
        print(f"Detected file type: {file_type}", file=sys.stderr)
    
    if file_type not in ['SAM', 'BAM']:
        raise ValueError(f"Unsupported file type: {file_type}")
    
    snp_dict = parse_phased_snp(phased_snp_file)
    comments = []
    marked_reads = []
    
    if file_type == 'BAM':
        # Count total alignments for progress bar
        if show_progress:
            total_alignments = 0
            with pysam.AlignmentFile(sam_file, "rb") as bam_file:
                for _ in bam_file:
                    total_alignments += 1
        with pysam.AlignmentFile(sam_file, "rb") as bam_file:
            # Process alignments with progress bar
            if show_progress:
                alignment_iter = tqdm(bam_file, total=total_alignments, desc="Processing alignments")
            else:
                alignment_iter = bam_file
            
            for alignment in alignment_iter:
                # Convert alignment to entry format
                entry = convert_alignment_to_entry(alignment)
                
                # Skip unmapped reads
                if alignment.is_unmapped:
                    continue
                
                mapq = int(entry[4])
                if mapq < min_mapq:
                    continue

                r_s, r_e, q_s, q_e = entry_align_pos(entry)

                # Phase determination
                if snp_dict and entry[2] in snp_dict:
                    # Find overlapping SNPs according to reference positions
                    hit_snp_sites = Interval.find_ovlp(snp_dict[entry[2]], r_s, r_e)
                    if hit_snp_sites:
                        for v in hit_snp_sites:
                            marked_read = check_allele(
                                entry, r_s, v,
                                append_features=append_features,
                                min_baseq=min_baseq,
                                verbose=verbose
                                )
                            if marked_read is not None:
                                marked_reads.append(marked_read)
                    else:
                        continue
                else:
                    continue
    else:  # SAM file
        # Count total lines for progress bar
        total_lines = count_lines_in_file(sam_file)
        
        with open(sam_file, 'r') as f_in:
            # Process with progress bar if verbose > 0
            lines = f_in
            if show_progress:
                lines = tqdm(lines, total=total_lines, desc="Processing SAM lines")
            
            for i, line in enumerate(lines):
                t = line.strip().split("\t")
                
                # Handle header lines
                if t[0].startswith('@'):
                    if t[0] == '@SQ':
                        sn = None
                        ln = None
                        for i in range(1, len(t)):
                            m = re.search(r'(LN|SN):(\S+)', t[i])
                            if m:
                                if m.group(1) == 'SN':
                                    sn = m.group(2)
                                elif m.group(1) == 'LN':
                                    ln = int(m.group(2))
                        
                        if sn is None or ln is None:
                            raise ValueError("missing SN or LN at an @SQ line")
                        # print(f"#chromosome: {sn} {ln}")
                        comments.append(f"#chromosome: {sn} {ln}")
                    else:
                        continue
                    continue
                else:
                    entry = t
                    mapq = int(entry[4])
                    if mapq < min_mapq:
                        continue

                    r_s, r_e, q_s, q_e = entry_align_pos(entry)

                    # Phase determination
                    if snp_dict and entry[2] in snp_dict:
                        # Find overlapping SNPs according to reference positions
                        hit_snp_sites = Interval.find_ovlp(snp_dict[entry[2]], r_s, r_e)
                        if hit_snp_sites:
                            for v in hit_snp_sites:
                                res = check_allele(
                                    entry, r_s, v,
                                    append_features=append_features,
                                    min_baseq=min_baseq,
                                    verbose=verbose
                                    )
                                if res is not None:
                                    marked_reads.append(res)
                        else:
                            continue
                    else:
                        continue
    
    res = pd.DataFrame(
        marked_reads,
        columns=['chrom', 'pos', 'allele', 'read_name', 'relative_pos'] + (append_features if append_features else [])
        )
    if outfile is not None:
        res.to_csv(outfile, sep="\t", index=False)
        return outfile
    else:
        return res

def sam_count_alleles(df):
    """
    Count reference, alternative, and other alleles for each chromosome and site.

    This function takes a DataFrame containing genetic data and calculates the count
    of reference alleles, alternative alleles, and other alleles for each unique
    chromosome and position combination.

    Args:
        df (pandas.DataFrame): A DataFrame with columns 'chrom', 'pos', 'allele',
                               'input_ref', and 'input_alt'. Each row represents
                               an allele observation at a specific genomic location.

    Returns:
        pandas.DataFrame: A DataFrame with columns 'chrom', 'pos', 'ref', 'alt', 'other'.
                          Each row represents a unique chromosome and position combination
                          with the counts of reference alleles ('ref'), alternative alleles ('alt'),
                          and other alleles ('other').
    """
    # Create a copy to avoid modifying the original dataframe
    df_copy = df.copy()
    
    # Create a new column to categorize each allele
    df_copy['allele_category'] = np.select(
        [df_copy['allele'] == df_copy['input_ref'], 
         df_copy['allele'] == df_copy['input_alt']],
        ['ref', 'alt'],
        default='other'
    )
    
    # Group by chrom and pos, and count the occurrences of each allele category
    result = df_copy.groupby(['chrom', 'pos'])['allele_category'].value_counts().unstack(fill_value=0)
    
    # Ensure all required columns exist (ref, alt, other)
    for col in ['ref', 'alt', 'other']:
        if col not in result.columns:
            result[col] = 0
    
    # Select only the required columns
    result = result[['ref', 'alt', 'other']].reset_index()
    
    # Ensure the correct column order
    result = result[['chrom', 'pos', 'ref', 'alt', 'other']]
    
    return result

### --- calculate haplotype score on segments --- ###


def hap_score_gb(df):
    df = df.droplevel(0)
    return abs(df.loc["0"] - df.loc["1"])/abs(df.loc["0"]+df.loc["1"])
def get_chrom_hap_score(dump_dir):
    """
    Get chrom_hap_score of samples.
    Input:
        dump_dir: rd/dump generated by hires_utils
    Output:
        pd.Series
    """
    res = pd.concat((pd.read_pickle(file) for file in [os.path.join(dump_dir, i) for i in os.listdir(dump_dir)]),axis=1)
    res.columns = [i.split("_")[0] for i in os.listdir(dump_dir)]

    dfs = {}
    for index, df in res.groupby(level=0):
        dfs[index] = hap_score_gb(df)
    chrom_hap_score = pd.DataFrame(dfs).mean(axis=1)
    return chrom_hap_score

def count_chrom_phased(filename:str, *args)->pd.Series:
    """
    Count number of phased (with separable SNP) segments per chromosome from a .seg.gz file.
    Input:
        filename: Path to the .seg.gz file.
    Output:
        A pandas Series with multi-index (chromosome, phasing) and counts.
    """
    #print(f"Processing file: {filename}")
    df = parse_seg(filename)
    cp_count = df.value_counts(["chrom","phase"])
    return cp_count
@mt(
    "hic_basic.phasing"
)
def mt_count_chrom_phased(filename:str)->pd.Series:
    pass

def count_phased_seg_perbin(filename, genome="mm10",binsize=1e6):
    """
    Count the number of phased segments per bin of a given genome.
    Input:
        filename: str, path to the .seg file
        genome: str, genome name (default: "mm10")
        binsize: int, size of the bins (default: 1e6)
    Output:
        pandas DataFrame with columns: "chrom", "start", "phase", "count"
    """
    genome = GenomeIdeograph(genome)
    # --- read in the .seg file --- #
    df = parse_seg(filename)

    # --- transform to pairs format to fit append_bins method --- #
    fake_pairs = pd.DataFrame(
        {
            "chr1" : df["chrom"],
            "pos1" : df["start"],
            "chr2" : df["chrom"],
            "pos2" : df["end"],
            "phase0" : df["phase"],
            "phase1" : df["phase"],
        }
    ).assign(
        readID = ".",
        strand1 = "+",
        strand2 = "+"
    )

    fake_pairs = genome.append_bins(
        fake_pairs,
        binsize=binsize,
    )

    phase_bin_counts = fake_pairs.groupby(
        ["chrom1", "start1", "phase0"]
    ).size().reset_index(name="count")
    phase_bin_counts.columns = ["chrom", "start", "phase", "count"]
    phase_bin_counts = phase_bin_counts.set_index(["chrom", "start", "phase"])["count"]
    return phase_bin_counts