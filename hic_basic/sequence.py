import pysam
import pandas as pd
def parse_fastq_id(line):
    """Parse a FASTQ ID line into (prefix, end) tuple or None if invalid.
    
    Args:
        line (str): The FASTQ ID line (e.g., '@M06168:... 1:N:0:1')
        
    Returns:
        tuple or None: (prefix, end) if valid, else None.
            prefix: The read identifier prefix (e.g., '@M06168:...').
            end: The read end identifier (e.g., '1').
    """
    line = line.strip()
    if len(line) < 2:
        return None
    # Split into prefix and end information
    parts = line.split(maxsplit=1)
    if len(parts) != 2:
        return None
    prefix = parts[0]
    end_info = parts[1].split(':', 1)[0]  # Extract first field (e.g., '1' or '2')
    return (prefix, end_info)

def read_fastq_ids(file_path):
    """Read FASTQ file and return a dictionary of {prefix: end}.
    
    This function correctly identifies ID lines by tracking the 4-line FASTQ structure.
    
    Args:
        file_path (str): Path to the FASTQ file
        
    Returns:
        dict: {prefix: end} mapping for all valid ID lines
    """
    ids = {}
    with open(file_path, 'r') as f:
        line_counter = 0  # Track line position in the 4-line FASTQ record
        for line in f:
            line = line.strip()
            if line_counter == 0:
                # Process ID line (first line of 4-line record)
                if not line.startswith('@'):
                    continue  # Skip invalid ID lines
                parsed = parse_fastq_id(line)
                if parsed:
                    prefix, end = parsed
                    ids[prefix] = end
            line_counter = (line_counter + 1) % 4  # Cycle through 0-3
    return ids

def compare_fastq_pairs(fq1_path, fq2_path, showstats=False):
    """Compare paired-end FASTQ files and return statistics.
    
    Args:
        fq1_path (str): Path to FASTQ file 1 (R1)
        fq2_path (str): Path to FASTQ file 2 (R2)
        
    Returns:
        dict: Statistics including:
            'total_fq1' (int): Total reads in FQ1
            'total_fq2' (int): Total reads in FQ2
            'valid_pairs' (int): Valid paired reads
            'invalid_pairs' (int): Reads with mismatched ends
            'missing_in_fq1' (int): Reads in FQ2 missing from FQ1
            'missing_in_fq2' (int): Reads in FQ1 missing from FQ2
    """
    ids1 = read_fastq_ids(fq1_path)
    ids2 = read_fastq_ids(fq2_path)
    
    stats = {
        'total_fq1': len(ids1),
        'total_fq2': len(ids2),
        'valid_pairs': 0,
        'invalid_pairs': 0,
        'missing_in_fq1': 0,
        'missing_in_fq2': 0
    }
    
    # Check FQ1 entries against FQ2
    for prefix in ids1:
        end1 = ids1[prefix]
        if prefix in ids2:
            end2 = ids2[prefix]
            # Check for valid pair (1 <-> 2)
            if (end1 == '1' and end2 == '2') or (end1 == '2' and end2 == '1'):
                stats['valid_pairs'] += 1
            else:
                stats['invalid_pairs'] += 1
        else:
            stats['missing_in_fq2'] += 1
    
    # Check FQ2 entries missing in FQ1
    for prefix in ids2:
        if prefix not in ids1:
            stats['missing_in_fq1'] += 1
    if showstats:
        print_stats(stats)
    return stats

def print_stats(stats):
    """Print comparison statistics in a human-readable format."""
    print(f"Total reads in FQ1: {stats['total_fq1']}")
    print(f"Total reads in FQ2: {stats['total_fq2']}")
    print(f"Valid paired reads: {stats['valid_pairs']}")
    print(f"Invalid pairs (end mismatch): {stats['invalid_pairs']}")
    print(f"Missing in FQ1: {stats['missing_in_fq1']} (present in FQ2)")
    print(f"Missing in FQ2: {stats['missing_in_fq2']} (present in FQ1)")

def count_CpG(bed_df:pd.DataFrame, fasta:str)-> pd.DataFrame:
    """
    Calculate CpG ratio for each bin in bed_df.
    Input:
        bed_df: pd.DataFrame, bed file with columns chrom, start, end
        fasta: str, path to the fasta file
    Output:
        pd.DataFrame, bed file with additional column CpG_ratio
    Note:
        1. omit all Ns in the sequence
        2. theoretically, max CpG ratio is 0.5 here.
        3. only bins with CpG ratio >= 0 are returned
    """
    bed_df = bed_df.copy()
    with pysam.FastaFile(fasta) as fa:
        bed_df["sequence"] = bed_df.apply(
            lambda row: fa.fetch(row['chrom'], int(row['start']), int(row['end'])),
            axis=1
            )
    def count_cg(seq):
        seq = seq.upper()
        cg_count = seq.count('CG')
        # string with ACGT
        total_count = sum([seq.count(base) for base in 'ACGT'])
        return cg_count / total_count if total_count > 0 else -1
    bed_df['CpG'] = bed_df['sequence'].apply(count_cg)
    res = bed_df.drop('sequence', axis=1)
    res = res[res['CpG'] >= 0]
    return res