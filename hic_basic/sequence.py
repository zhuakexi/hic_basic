import gzip
import random
from collections import Counter

import pysam
import pandas as pd
import plotly.graph_objects as go

from hires_utils.gcount import count_fastq
from tqdm import tqdm
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
    """Read FASTQ file (supports gzipped) and return a dictionary of {prefix: end}.
    
    This function correctly identifies ID lines by tracking the 4-line FASTQ structure.
    
    Args:
        file_path (str): Path to the FASTQ file (can be gzipped)
        
    Returns:
        dict: {prefix: end} mapping for all valid ID lines
    """
    ids = {}
    open_func = gzip.open if file_path.endswith('.gz') else open
    with open_func(file_path, 'rt') as f:
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


### --- functions to decode library structure --- ###


def search_primers(fq_fp, primers, sample_n=None, sample_r=None, seed=None):
    """
    Search primers in fastq file and return distribution of primer start positions.
    Will search forward, reverse, complement and reverse complement.
    
    Input:
        fq_fp: str, fastq file path
        primers: list of str, primers to search
            If list, search all primers as one.
        sample_n: int, optional, number of reads to sample
        sample_r: float, optional, fraction of reads to sample (0 < sample_r <= 1)
        seed: int, optional, random seed for reproducibility
    
    Output:
        res: dict, primer distribution of each form
            key: form type, value: list of start positions
    """
    # Convert primers to uppercase and process all forms
    fq_fp = str(fq_fp)
    if isinstance(primers, str):
        primers = [primers]
    primers = [p.upper() for p in primers]
    all_primer_forms = []
    for primer in primers:
        forward = primer
        reverse = primer[::-1]
        complement = primer.translate(str.maketrans('ATCG', 'TAGC'))
        rev_comp = complement[::-1]
        all_primer_forms.extend([
            (forward, 'forward'),
            (reverse, 'reverse'),
            (complement, 'complement'),
            (rev_comp, 'reverse_complement')
        ])

    # Initialize counters and result structure
    res = {
        "forward": [],
        "reverse": [],
        "complement": [],
        "reverse_complement": []
    }

    # Determine file opener based on compression
    if fq_fp.endswith(".gz"):
        open_func = gzip.open
    else:
        open_func = open

    total_reads = count_fastq(fq_fp, form="s")
    print(f"Total reads: {total_reads}")

    # Determine sampling
    sampled_indices = None
    if sample_n or sample_r:
        if seed is not None:
            random.seed(seed)
        if sample_r:
            sample_n = int(total_reads * sample_r)
        sampled_indices = sorted(random.sample(range(total_reads), sample_n))
        print(f"Sampling {len(sampled_indices)} reads.")

    with open_func(fq_fp, "rt") as fh:
        if sampled_indices:
            sampled_iter = iter(sampled_indices)
            next_sample = next(sampled_iter, None)
        else:
            next_sample = None

        for idx, lines in enumerate(tqdm(zip(*[fh]*4), total=total_reads)):
            if next_sample is not None:
                if idx < next_sample:
                    continue  # Skip to the next sampled row
                elif idx == next_sample:
                    next_sample = next(sampled_iter, None)  # Move to the next sampled index

            name, seq, _, qual = lines
            seq = seq.strip().upper()

            # Search for all primer forms in the current sequence
            for form, form_type in all_primer_forms:
                start = 0
                pos = seq.find(form, start)
                if pos != -1:
                    res[form_type].append(pos)

    # Transform to frequency
    for key in res:
        res[key] = pd.Series(Counter(res[key]))
    res = pd.concat(res, axis=1)

    if seed is not None:
        random.seed(None)  # Reset random seed
    return res

def plot_search_primers(result, output_file=None):
    """
    Visualize the primer search result using bar plots with 4 subplots.
    
    Args:
        result (pd.DataFrame): Primer search result with columns for each form.
        output_file (str, optional): Path to save the plot as an HTML file. If None, the plot is shown.
    """
    fig = make_subplots(
        rows=4, cols=1,
        subplot_titles=["forward", "reverse", "complement", "reverse_complement"],
        shared_yaxes=True
    )

    for row, col, i, form in filling_l2r_plotly(4, 1, ["forward", "reverse", "complement", "reverse_complement"]):
        if form in result.columns:
            fig.add_trace(
                go.Bar(
                    x=result.index,
                    y=result[form],
                    name=form,
                    marker_color=f"rgba({i*50}, {i*50}, 255, 0.7)",
                ),
                row=row,
                col=col
            )

    fig.update_layout(
        title="Primer Search Results",
        xaxis_title="Position",
        yaxis_title="Frequency",
        barmode="group",
        template="plotly_white",
        height=600,
        width=800
    )

    if output_file:
        # output png
        fig.write_image(output_file)
    else:
        return fig