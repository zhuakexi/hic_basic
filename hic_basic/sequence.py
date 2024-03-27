import pysam
import pandas as pd

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