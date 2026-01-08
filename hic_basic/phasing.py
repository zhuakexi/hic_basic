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
        str: The file type ('SAM', 'BAM', 'SAMgz', or 'unknown').
    """
    try:
        # 首先检查文件扩展名（快速提示）
        ext = os.path.splitext(file_path)[1].lower()
        
        with open(file_path, 'rb') as f:
            # 读取文件开头的魔数
            magic = f.read(4)
            f.seek(0)
            
            # 1. 检查是否是 BAM 文件
            # BAM 是 gzip 压缩的，并且解压后以 BAM\x01 开头
            if magic[:2] == b'\x1f\x8b':
                try:
                    with gzip.open(file_path, 'rb') as gz:
                        decompressed_header = gz.read(4)
                        if decompressed_header == b'BAM\x01':
                            return 'BAM'
                except (gzip.BadGzipFile, OSError):
                    # 不是有效的 gzip 文件或不是 BAM
                    pass
            
            # 2. 检查是否是 gzip 压缩的 SAM 文件
            if magic[:2] == b'\x1f\x8b':
                try:
                    with gzip.open(file_path, 'rb') as gz:
                        # 读取足够的数据来检测 SAM 特征
                        decompressed_data = gz.read(4096)
                        
                        # 尝试解码为文本
                        try:
                            text = decompressed_data.decode('utf-8', errors='ignore')
                            if is_sam_content(text):
                                return 'SAMgz'
                        except UnicodeDecodeError:
                            # 不是有效的文本
                            pass
                except (gzip.BadGzipFile, OSError):
                    # 不是有效的 gzip 文件
                    pass
            
            # 3. 检查是否是未压缩的 SAM 文件
            # 读取文件内容进行检测
            content = f.read(4096)
            try:
                text = content.decode('utf-8', errors='ignore')
                if is_sam_content(text):
                    return 'SAM'
            except UnicodeDecodeError:
                # 不是文本文件
                pass
                
    except (IOError, OSError) as e:
        print(f"Error reading file: {e}")
    
    return 'unknown'


def is_sam_content(text, max_lines=20):
    """
    检查文本内容是否符合 SAM 格式特征
    
    Args:
        text (str): 要检查的文本内容
        max_lines (int): 最多检查的行数
    
    Returns:
        bool: 如果是 SAM 格式返回 True，否则返回 False
    """
    lines = text.split('\n')
    
    # SAM 头行的有效前缀
    sam_headers = ['@HD', '@SQ', '@RG', '@PG', '@CO']
    
    for line in lines[:max_lines]:
        line = line.strip()
        if not line:
            continue
        
        # 检查是否是有效的 SAM 头行
        if line.startswith('@'):
            if any(line.startswith(h) for h in sam_headers):
                return True
        # 检查是否是有效的 SAM 数据行
        elif '\t' in line:
            fields = line.split('\t')
            # SAM 数据行至少有11个字段
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
                        res = [chrom, v[1], c, entry[0], p_adj] # snp chrom, snp pos, other allele, read name, relative pos
                        if verbose >= 1:
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
    elif file_type in ["SAM", "SAMgz"]:  # SAM file
        # Count total lines for progress bar
        if file_type == "SAMgz":
            open_func = gzip.open
        else:
            open_func = open

        if show_progress:
            with open_func(sam_file, 'rt') as f:
                total_lines = sum(1 for line in f)
        
        with open_func(sam_file, 'rt') as f_in:
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
    else:
        raise ValueError(f"Unsupported file type: {file_type}")

    res = pd.DataFrame(
        marked_reads,
        columns=['chrom', 'pos', 'allele', 'read_name', 'relative_pos'] + (append_features if append_features else [])
        )
    if outfile is not None:
        res.to_csv(outfile, sep="\t", index=False)
        return outfile
    else:
        return res

def sam_count_alleles(allele_df, ref_snp_file, gt_strategy="max_allele"):
    """
    Count reference, alternative, and other alleles for each chromosome and site.

    This function takes a DataFrame containing genetic data and calculates the count
    of reference alleles, alternative alleles, and other alleles for each unique
    chromosome and position combination.
    TODO: Better GT strategies. For example, is 0/1 better than 1/4?

    Args:
        allele_df (pandas.DataFrame): A DataFrame with columns 'chrom', 'pos', 'allele',
                               'input_ref', and 'input_alt'. Each row represents
                               an allele observation at a specific genomic location.
        gt_strategy (str): Strategy for generating the GT column. 
                          "max_allele": Choose the allele with the highest count (ref or alt) and convert to ATGC
                          "no_conflict": Choose the non-zero allele if one is 0 and the other > 0, otherwise NA

    Returns:
        pandas.DataFrame: A DataFrame with columns 'chrom', 'pos', 'ref', 'alt', 'other', 'gt'.
                          Each row represents a unique chromosome and position combination
                          with the counts of reference alleles ('ref'), alternative alleles ('alt'),
                          other alleles ('other'), and the genotype ('gt').
    """
    # Create a copy to avoid modifying the original dataframe
    ref_snp_df = pd.read_table(
        ref_snp_file,
        names=['chrom', 'pos', 'input_ref', 'input_alt'],
        usecols=[0,1,2,3]
        )
    df = pd.merge(
        allele_df,
        ref_snp_df,
        on=['chrom', 'pos'],
        how='left'
    )
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

    # --- work on aggregated df to add GT column --- #

    # Get the ref and alt alleles for each chrom and pos from the original df
    allele_map = df[['chrom', 'pos', 'input_ref', 'input_alt']].drop_duplicates(subset=['chrom', 'pos']).reset_index(drop=True)
    result = result.merge(allele_map, on=['chrom', 'pos'], how='left')
    
    # Add the GT column based on the specified strategy using vectorized operations
    if gt_strategy == "max_allele":
        # Vectorized operation:
        condition1 = (result['ref'] > result['alt']) & (result['ref'] > result['other'])  # ref > alt and ref > other
        condition2 = (result['alt'] > result['ref']) & (result['alt'] > result['other'])  # alt > ref and alt > other

        result['gt'] = np.select(
            [condition1, condition2],
            [result['input_ref'], result['input_alt']],
            default=pd.NA
        )
        
    elif gt_strategy == "no_conflict":
        # Vectorized operation: if one is 0 and other > 0, choose the > 0 allele; if both > 0, return NA
        condition1 = (result['ref'] > 0) & (result['alt'] == 0) & (result['other'] == 0) # ref > 0, alt = 0, other = 0
        condition2 = (result['ref'] == 0) & (result['alt'] > 0) & (result['other'] == 0)  # ref = 0, alt > 0, other = 0
        condition3 = (result['ref'] > 0) & (result['alt'] > 0)   # both > 0
        
        result['gt'] = np.select(
            [condition1, condition2, condition3],
            [result['input_ref'], result['input_alt'], pd.NA],
            default=pd.NA
        )
        
    else:
        raise ValueError(f"Invalid gt_strategy: {gt_strategy}. Valid options are 'max_allele' or 'no_conflict'.")
    
    # Drop the temporary columns used for allele mapping
    result = result.drop(['input_ref', 'input_alt'], axis=1)
    
    # Ensure the correct column order
    result = result[['chrom', 'pos', 'ref', 'alt', 'other', 'gt']]
    
    return result

def do_sam2vcf_allele_count(sam_file, ref_snp_file, output):
    """
    Given a SAM file and a reference SNP file, generate a parquet file with allele counts of each SNP.
    
    {pipeline}[Gamete SNP Phasing Pipeline][1.count alleles]:
        `do_sam2vcf_allele_count` for each gamete SAM file to get allele counts at each SNP position.
    
    Parameters:
        sam_file (str): Path to the SAM file, can be sam.gz, bam, or sam.
        ref_snp_file (str): Path to the reference SNP file in tab-delimited format with columns: chrom, pos, ref_allele, alt_allele.
    Output (str): Path to the output parquet file.
    Note: output is aligned to ref_snp_file 
    """
    assert output.endswith(".parquet"), "Output file must be a .parquet file"
    # Load reference SNPs
    ref_snps = pd.read_table(ref_snp_file)
    allele_df = sam_mark_alleles(
        sam_file,
        ref_snp_file,
        output=None,
        show_progress=False,
        verbose=0
    )
    allele_df_count = sam_count_alleles(
        allele_df,
        ref_snp_file,
        gt_strategy="no_conflict",
    )
    snp = pd.read_table(ref_snp_file, names=["chrom", "pos", "ref_allele", "alt_allele"])
    extended_allele_count = pd.merge(
        snp,
        allele_df_count,
        how="left",
        on=["chrom", "pos"]
    )[["chrom","pos","gt"]]
    extended_allele_count.to_parquet(
        output,
        index=True
        )
    return output

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