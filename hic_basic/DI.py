from functools import reduce

import dask.array as da
import dask.dataframe as dd
import numpy as np
import pandas as pd
from dask import delayed
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
def simpleDiff_test(df, m1, m2):
    """
    Input:
        df: n * (m1 + m2) dask dataframe, n is pixels to be tested, m1 and m2 are number of samples in each group
    Output:
        df: n * 1 dask dataframe, pvalue
    """
    assert(isinstance(df, dd.DataFrame))
    stats_values = df.apply(
        lambda x: mannwhitneyu(x[:m1], x[m1:m1+m2], alternative="two-sided").pvalue,
        axis=1,
        meta=('pvalue', 'f8')
    )
    return dd.from_dask_array(stats_values, columns=["pvalue"])
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
    # read imputes
    dfs = []
    for i, imputed_path in enumerate(groupA + groupB):
        #df = cis_distance_graph_df(_3dg_path, chrom=chrom, genome=genome, fo=None, max_dist=max_dist, binsize=binsize)
        df = dd.read_parquet(
            imputed_path
            )
        df = normalize_band(df, chrom=chrom, max_dist=max_dist, binsize=binsize)
        df = ideograph.join_pixel_id(df, binsize=binsize)
        df = df.set_index("pixel_id", npartiitions="auto")
            
        # transform normalized distance to dask array
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
        length = end_id - start_id + 1
        values = da.full(length, np.nan, chunks=chunksize)
        index = da.arange(start_id, end_id + 1, chunks=chunksize)
        new_df = dd.from_dask_array(values, columns=["z_normalized_distance"])
        new_df.index = dd.from_dask_array(index)
        new_df["z_normalized_distance"] = new_df["z_normalized_distance"].fillna(df["z_normalized_distance"])
        #print(new_df.npartitions)
        dfs.append(
            new_df["z_normalized_distance"].to_dask_array(lengths=True)
        )

    combined_df = da.stack(dfs, axis=1)
    combined_df = dd.from_dask_array(combined_df, columns=[f"sample{i}" for i in range(m1 + m2)])
    print("combined_df.npartitions:",combined_df.npartitions)
    meanA = combined_df.iloc[:, :m1].mean(axis=1)
    meanB = combined_df.iloc[:, m1:m1+m2].mean(axis=1)
    concerned_indi = ((meanA < max_3d_dist) | (meanB < max_3d_dist)) & (meanA.notnull() & meanB.notnull())

    # 调用 simpleDiff_test 计算 p 值
    pvalues = simpleDiff_test(combined_df[concerned_indi], m1, m2)

    # 整理结果
    # id_df = combined_df[['chrom1', 'start1', 'chrom2', 'start2']]
    # id_df['end1'] = id_df['start1'] + binsize
    # id_df['end2'] = id_df['start2'] + binsize
    # result_df = dd.concat([id_df.loc[concerned_indi], pvalues], axis=1)
    result_df = pvalues
    result_df['meanA'] = meanA
    result_df['meanB'] = meanB
    result_df['diff'] = meanA - meanB

    # 保存结果
    result_df.to_parquet(fo)
    client.close()
    return result_df