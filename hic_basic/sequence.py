import pysam
import pandas as pd

def count_CpG(bed_df:pd.DataFrame, fasta:str)-> pd.DataFrame:
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
    bed_df['CG_ratio'] = bed_df['sequence'].apply(count_cg)
    res = bed_df.drop('sequence', axis=1)
    res = res[res['CG_ratio'] >= 0]
    return res