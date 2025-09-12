import gzip
import os
import random
from collections import Counter

import pysam
import pandas as pd
import plotly.graph_objects as go
from hires_utils.gcount import count_fastq
from plotly.subplots import make_subplots
from tqdm import tqdm

def count_fastq_worker(file_path, form="s"):
    res = count_fastq(file_path, form=form)
    return pd.Series(
        [res],
        name="n_reads"
    )

def parse_fastq_id(line):
    """Parse a FASTQ ID line into (prefix, end) tuple or None if invalid.

    Args:
        line (str): The FASTQ ID line (e.g., '@M06168:... 1:N:0:1')
        
    Returns:
        tuple or None: (prefix, end) if valid, else None.
            prefix: The read identifier prefix (e.g., '@M06168:...').
            end: The read end identifier (e.g., '1').

    Notes
    -----
    Example of Read ID: `@M06168:86:000000000-LMWTD:1:1101:14802:1380 1:N:0:1`

    The Read ID is structured as follows:
    `@<instrument>:<run_id>:<flowcell_id>:<lane>:<tile>:<x_pos>:<y_pos> <read>:<is_filtered>:<control_number>:<index_seq>`

    **Detailed Field Explanation:**

    1.  **Instrument (`M06168`)**
        The unique identifier of the sequencing instrument.
        Example: `M06168` (Likely a MiSeq instrument).

    2.  **Run ID (`86`)**
        The number assigned to the specific sequencing run by the operator.

    3.  **Flowcell ID (`000000000-LMWTD`)**
        The unique barcode of the flowcell used.
        Example: `000000000-LMWTD`.

    4.  **Lane (`1`)**
        The lane number on the flowcell where the cluster was located.
        Example: `1`.

    5.  **Tile (`1101`)**
        The specific tile within the lane where the cluster was imaged.

    6.  **X-coordinate (`14802`)**
        The X-position of the cluster within the tile (in pixels).

    7.  **Y-coordinate (`1380`)**
        The Y-position of the cluster within the tile (in pixels).

    8.  **Read Number (`1`)**
        Member of a pair. `1` for the first read in a pair (R1), `2` for the second (R2).

    9.  **Filter Flag (`N`)**
        Indicates if the read passed quality filtering.
        - `Y`: Read failed quality filter (should be discarded).
        - `N`: Read passed quality filter (is reliable).

    10. **Control Number (`0`)**
        Reserved field. `0` indicates a regular sample read. Non-zero values indicate control reads.

    11. **Index Sequence Number (`1`)**
        Identifies which barcode read this is in a multiplexed sequencing run.
        - `1`: This read contains the index sequence (e.g., i7 index).
        - `2`: This read contains the second index (e.g., i5 index in a dual-index experiment).

    **Example Interpretation:**

    For the Read ID
    `@M06168:86:000000000-LMWTD:1:1101:14802:1380 1:N:0:1`:

    - **Instrument:** M06168
    - **Run ID:** 86
    - **Flowcell ID:** 000000000-LMWTD
    - **Lane:** 1
    - **Tile:** 1101
    - **Cluster Coordinates (X, Y):** (14802, 1380)
    - **Read Number:** 1 (This is the first read, R1)
    - **Filter Flag:** N (This is a high-quality pass)
    - **Control Number:** 0 (This is a regular sample read)
    - **Index Sequence Number:** 1 (This read's sequence is the i7 index barcode)

    This read is an index read (i7) from a paired-end run that passed quality control.
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


def search_primers(fq_fp, primers, sample_n=None, sample_r=None, seed=None, keep_hits=0):
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
        keep_hits: keep number of hits for further analysis, default 0
    
    Output:
        if keep_hits == 0:
            result (pd.DataFrame): Primer search result with columns for each form.
        if keep_hits > 0:
            result,
            hits (dict of list): sampled hits for each primer form.
    """
    assert(keep_hits >= 0), "keep_hits must be >= 0"
    # Convert primers to uppercase and process all forms
    fq_fp = str(fq_fp)
    if isinstance(primers, str):
        primers = [primers]
    primers = [p.upper() for p in primers]
    all_primer_forms = [] # list(primers) of list(4):[(forward, 'forward'), ...] of tuples(2)
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
    hits = {
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
                    if len(hits[form_type]) < keep_hits:
                        hits[form_type].append((name.strip(), pos, seq, qual.strip()))

    # Transform to frequency
    for key in res:
        res[key] = pd.Series(Counter(res[key]))
    res = pd.concat(res, axis=1)

    if seed is not None:
        random.seed(None)  # Reset random seed
    if keep_hits > 0:
        return res, hits
    else:
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



### --- helper functions for fastq --- ###




def add_suffix_to_fastq(input_file: str, output_file: str, skip_incomplete=True) -> None:
    """
    Process a FASTQ file by appending "/1" or "/2" (MGI/DNBSEQ flavor) to read IDs if not already present (illumina_flavor).

    For Illumina reads, the prefix of the read ID is determined according to the comment part (e.g. `/1` for `1:N:0:1`).

    Parameters:
    input_file (str): Path to the input FASTQ file (can be gzip-compressed if ending with .gz).
    output_file (str): Path to the output FASTQ file (will be gzip-compressed if input is compressed).

    Returns:
    None
    """
    # Determine if the input file is compressed based on its extension
    input_compressed = input_file.endswith('.gz')
    output_compressed = output_file.endswith('.gz')
    
    # Ensure output compression matches input compression
    if input_compressed != output_compressed:
        raise ValueError("Input and output compression formats must match.")
    
    # Open input file with appropriate handler
    if input_compressed:
        in_fh = gzip.open(input_file, 'rt')
    else:
        in_fh = open(input_file, 'r')
    
    # Open output file with appropriate handler
    if output_compressed:
        out_fh = gzip.open(output_file, 'wt')
    else:
        out_fh = open(output_file, 'w')
    
    with in_fh, out_fh:
        while True:
            # Read four lines per record
            header = in_fh.readline().rstrip()
            if not header:  # End of file
                break
            sequence = in_fh.readline().rstrip()
            plus_line = in_fh.readline().rstrip()
            quality = in_fh.readline().rstrip()
            
            # Check if we have a complete record
            if not (header and sequence and plus_line and quality):
                if skip_incomplete:
                    print("Warning: Incomplete FASTQ record found and skipped.")
                    continue
                else:
                    raise ValueError("Incomplete FASTQ record")
            
            # Parse header line: @name comment
            if header.startswith('@'):
                name_comment = header[1:].split(' ', 1)
                name = name_comment[0]
                comment = name_comment[1] if len(name_comment) > 1 else None
            else:
                raise ValueError("Invalid FASTQ header line: " + header)
            
            # Check if read name ends with /1 or /2
            if not (name.endswith('/1') or name.endswith('/2')):
                if comment:
                    # Determine suffix from comment (e.g., '1:N:0:1' -> '/1')
                    end_info = comment.split(':', 1)[0]
                    if end_info == '1':
                        name += '/1'
                    elif end_info == '2':
                        name += '/2'
                    else:
                        # If end_info is not '1' or '2', default to '/1'
                        name += '/1'
                else:
                    # If no comment, default to '/1'
                    name += '/1'
            
            # Reconstruct the header line
            new_header = f"@{name}"
            if comment:
                new_header += f" {comment}"
            
            # Write the modified record to output
            out_fh.write(f"{new_header}\n{sequence}\n{plus_line}\n{quality}\n")

def extract_first_n_reads(input_file: str, output_file: str, n: int) -> None:
    """
    Extract the first N reads from a FASTQ file.

    Reads from the input FASTQ file and writes the first N reads to the output file.
    The compression format of the output matches the input (compressed if input ends with .gz).

    Parameters:
    input_file (str): Path to the input FASTQ file (can be gzip-compressed if ending with .gz).
    output_file (str): Path to the output FASTQ file.
    n (int): Number of reads to extract.

    Returns:
    None
    """
    # Determine if the input file is compressed based on its extension
    input_compressed = input_file.endswith('.gz')
    output_compressed = output_file.endswith('.gz')
    
    # Ensure output compression matches input compression
    if input_compressed != output_compressed:
        raise ValueError("Input and output compression formats must match.")
    
    # Open input file with appropriate handler
    if input_compressed:
        f_in = gzip.open(input_file, 'rt')
    else:
        f_in = open(input_file, 'r')
    
    # Open output file with appropriate handler
    if output_compressed:
        f_out = gzip.open(output_file, 'wt')
    else:
        f_out = open(output_file, 'w')
    
    # Counter for number of reads processed
    read_count = 0
    
    with f_in, f_out:
        while read_count < n:
            # Read four lines for one complete FASTQ record
            header = f_in.readline().strip()
            if not header:  # End of file
                break
                
            sequence = f_in.readline().strip()
            plus_line = f_in.readline().strip()
            quality = f_in.readline().strip()
            
            # Write the record to output
            f_out.write(f"{header}\n{sequence}\n{plus_line}\n{quality}\n")
            
            read_count += 1
    
    print(f"Extracted {read_count} reads to {output_file}")
