from functools import reduce

import dask
import dask.array as da
import dask.dataframe as dd
import numpy as np
import pandas as pd
from dask.distributed import Client
from hic_basic.binnify import GenomeIdeograph
from hic_basic.simpute import cis_distance_graph, cis_distance_graph_df
from scipy.stats import mannwhitneyu

def normalize_band(df, chrom=None, max_dist=2000000, binsize=20000):
    """
    Input:
        df: bedpe like dask dataframe, with columns chrom1, start1, end1, chrom2, start2, end2, value
            don't assume value name
        chrom: if not None, assume df has only intra pixels
        max_dist: max distance in impuation, used to estimiate partitions
        binsize: resolution of the data, used to estimiate partitions
    Output:
        df: same as input, with value normalized by band
        TODO: if max_dist or binsize is None, df is partitioned by chrom1, chrom2
    """
    value_col_name = df.columns[-1]
    df['band'] = (df['start2'] - df['start1']) // binsize
    #df = df.set_index('band', divisions=list(range(0, max_dist // binsize + 2)))
    #df = df.set_index('band', npartitions=max_dist // binsize + 2)
    df = df.set_index('band',npartitions="auto")

    # zscore normalization
    if chrom is None:
        grouper = ['chrom1', 'chrom2', 'band']
    else:
        grouper = "band"
    if isinstance(df, dd.DataFrame):
        df[f'z_normalized_{value_col_name}'] = df.groupby(grouper)[value_col_name].transform(
            lambda x: (x - x.mean()) / x.std(),
            meta=(f'z_normalized_{value_col_name}', 'f8')
        )
    else:
        df[f'z_normalized_{value_col_name}'] = df.groupby(grouper)[value_col_name].transform(
            lambda x: (x - x.mean()) / x.std()
        )

    # remove helper cols
    #df = df.drop(columns=['band', value_col_name])

    return df
def multiple_testing_correction(pvalues, correction_type="FDR"):
    """
    Consistent with R - print
    correct_pvalues_for_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05,
                                          0.069, 0.07, 0.071, 0.09, 0.1])
    from https://github.com/CoBiG2/cobig_misc_scripts/blob/master/FDR.py
    """
    from numpy import array, empty
    pvalues = array(pvalues)
    sample_size = pvalues.shape[0]
    qvalues = empty(sample_size)
    if correction_type == "Bonferroni":
        # Bonferroni correction
        qvalues = sample_size * pvalues
    elif correction_type == "Bonferroni-Holm":
        # Bonferroni-Holm correction
        values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
        values.sort()
        for rank, vals in enumerate(values):
            pvalue, i = vals
            qvalues[i] = (sample_size-rank) * pvalue
    elif correction_type == "FDR":
        # Benjamini-Hochberg, AKA - FDR test
        values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
        values.sort()
        values.reverse()
        new_values = []
        for i, vals in enumerate(values):
            rank = sample_size - i
            pvalue, index = vals
            new_values.append((sample_size/rank) * pvalue)
        for i in range(0, int(sample_size)-1):
            if new_values[i] < new_values[i+1]:
                new_values[i+1] = new_values[i]
        for i, vals in enumerate(values):
            pvalue, index = vals
            qvalues[index] = new_values[i]
    return qvalues
def test_row(row, m1, m2):
    """
    Input:
        row: a row of combined_df
        m1: number of samples in groupA
        m2: number of samples in groupB
    Output:
        p-value of mannwhitneyu test
    """
    # test
    return mannwhitneyu(row[:m1], row[m1:], alternative="two-sided")[1]
@dask.delayed
def process_file(imputed_path, chrom, binsize, chunksize):
    # load data
    df = dd.read_parquet(imputed_path)
    # zscore normalization
    df = df.repartition(npartitions=100)
    value_col_name = df.columns[-1]
    df['band'] = (df['start2'] - df['start1']) // binsize
    if chrom is None:
        grouper = ['chrom1', 'chrom2', 'band']
    else:
        grouper = "band"
    df[f'z_normalized_{value_col_name}'] = df.groupby(grouper)[value_col_name].transform(
        lambda x: (x - x.mean()) / x.std(),
        meta=(f'z_normalized_{value_col_name}', 'f8')
    )
    # add pixel_id
    df = df.map_partitions(
        lambda df : ideograph.join_pixel_id(
            df,
            binsize=binsize
        )
    )
    # extend to all possible pixels
    df = df.set_index("pixel_id", nparticitions=100)
    start_id = ideograph.pixel_id(
        pd.Series(
            [chrom, 0, binsize, chrom, 0, binsize],
            index = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]
        ),
        binsize=binsize
    )
    last_bin = ideograph.chromosomes.loc[chrom, "length"] - ideograph.chromosomes.loc[chrom, "length"] % binsize
    end_id = ideograph.pixel_id(
        pd.Series(
            [chrom, last_bin, ideograph.chromosomes.loc[chrom], chrom, last_bin, ideograph.chromosomes.loc[chrom]],
            index = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]
        ),
        binsize=binsize
    )
    all_possible_pixels = da.stack(
        [
            da.arange(start_id, end_id + 1, chunks=chunksize), # index
            da.full(end_id + 1 - start_id, np.nan, chunks=chunksize) # values
        ],
        axis=1
    )
    all_possible_pixels = dd.from_dask_array(
        da.arange(start_id, end_id + 1, chunks=chunksize),
        columns=["pixel_id"]
    )
    all_possible_pixels = all_possible_pixels.set_index("pixel_id", npartitions=100)
    extended_df = dd.merge(
        all_possible_pixels,
        df,
        how="left",
        left_index=True,
        right_index=True
    )
    return extended_df["z_normalized_distance"].to_dask_array(lengths=True)
def simpleDiff(groupA, groupB, chrom="chr1", genome="mm10", fo=None, max_3d_dist=5, 
               max_dist=2000000, binsize=20000, n_jobs=16, memory_limit='32GB'):
    """
    Input:
        groupA: imputed .parquet file path for groupA
        groupB: imputed .parquet file paths for groupB
        fo: output parquet file
        filt_fdr: filter pixels with fdr > filt_fdr, if None, no filtering
        n_jobs: number of jobs
    Output:
        df: n * (6+3) dask dataframe
            (chrom1, start1, end1, chrom2, start2, end2, meanA, meanB, diff, pvalue, qvalue)
    """
    chunksize = 100000

    client = Client(
        n_workers=n_jobs,
        threads_per_worker=1,
        memory_limit=memory_limit
        )
    m1, m2 = len(groupA), len(groupB)
    ideograph = GenomeIdeograph(genome)

    # load all data
    dfs = []
    for i, imputed_path in enumerate(groupA + groupB):
        dfs.append(process_file(imputed_path, chrom, binsize, chunksize))
    combined_df = da.stack(
        dd.compute(*dfs),
        axis=1
    )
    combined_df = dd.from_dask_array(combined_df, columns=[f"sample{i}" for i in range(m1 + m2)])
    print("combined_df.npartitions:",combined_df.npartitions)

    # calculate
    meanA = combined_df.iloc[:, :m1].mean(axis=1)
    meanB = combined_df.iloc[:, m1:m1+m2].mean(axis=1)
    #concerned_indi = ((meanA < max_3d_dist) | (meanB < max_3d_dist)) & (meanA.notnull() & meanB.notnull())
    #combined_df = combined_df.loc[concerned_indi]
    # pvalues = combined_df.map_partitions(
    #     lambda df: mannwhitneyu(df.iloc[:, :m1], df.iloc[:, m1:m1+m2], alternative="two-sided").pvalue,
    #     meta=('pvalue', 'f8')
    # )
    pvalues = combined_df.apply(test_row, axis=1, args=(m1, m2), meta=('pvalue', 'f8'))
    # meanA, meanB, diff, pvalues
    result_df = dd.concat(
        [
            meanA,
            meanB,
            meanA - meanB,
            pvalues
        ],
        axis=1
    )
    result_df.columns = ["meanA", "meanB", "diff", "pvalue"]

    # 保存结果
    result_df.to_parquet(fo)
    client.close()
    return result_df