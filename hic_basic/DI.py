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

def normalize_band(df):
    """
    Input:
        df: bedpe like dask dataframe, with columns chrom1, start1, end1, chrom2, start2, end2, value
            don't assume value name
    Output:
        df: same as input, with value normalized by band
    """
    value_col_name = df.columns[-1]
    df['band'] = df['start2'] - df['start1']

    # zscore normalization
    if isinstance(df, dd.DataFrame):
        # df = df.reset_index(drop=True)
        # grouped = df.groupby(['chrom1', 'chrom2', 'band'])[value_col_name]
        # z_scores = grouped.apply(lambda x: (x - x.mean()) / x.std(), meta=('x', 'f8'))
        # df[f'z_normalized_{value_col_name}'] = z_scores
        def zscore_with_reset_index(group):
            group = group.reset_index(drop=True)
            return (group[value_col_name] - group[value_col_name].mean()) / group[value_col_name].std(ddof=1)

        z_normalized = df.groupby(['chrom1', 'chrom2', 'band']).apply(zscore_with_reset_index, meta=('z_normalized_value', 'f8')).reset_index(drop=True)
        df[f'z_normalized_{value_col_name}'] = z_normalized
    else:
        df[f'z_normalized_{value_col_name}'] = df.groupby(['chrom1', 'chrom2', 'band'])[value_col_name].transform(lambda x: (x - x.mean()) / x.std())

    # remove helper cols
    df = df.drop(columns=['band', value_col_name])

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
def calculate_stats(row, m1, m2):
    group1 = row[:m1]
    group2 = row[m1:m1+m2]
    # Mann-Whitney U 测试, 注意 mannwhitneyu 需要非空数组
    if len(group1) > 0 and len(group2) > 0:
        pvalue = mannwhitneyu(group1, group2, alternative='two-sided').pvalue
    else:
        pvalue = np.nan
    return pvalue

def simpleDiff_test(df, m1, m2):
    """
    Input:
        df: n * (m1 + m2) dask dataframe, n is pixels to be tested, m1 and m2 are number of samples in each group
    Output:
        df: n * 1 dask dataframe, pvalue
    """
    assert(isinstance(df, dd.DataFrame))
    stats_values = df.apply(lambda row: calculate_stats(row, m1, m2), axis=1, meta=('x', 'f8'))

    stats_df = stats_values.to_frame()
    stats_df.columns = ['pvalue']

    return stats_df
def simpleDiff(groupA, groupB, chrom="chr1", genome="mm10", fo=None, filt_fdr=0.05, max_3d_dist=5, 
               max_dist=2000000, binsize=20000, n_jobs=16, memory_limit='32GB'):
    """
    Input:
        groupA: _3dg_file paths for groupA
        groupB: _3dg_file paths for groupB
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
    # 读取和处理 groupA 和 groupB 的文件
    dfs = []
    for i, _3dg_path in enumerate(groupA + groupB):
        df = cis_distance_graph_df(_3dg_path, chrom=chrom, genome=genome, fo=None, max_dist=max_dist, binsize=binsize)
        df = normalize_band(df)
        df = df.set_index("pixel_id")
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
        # new_df = dd.from_array(
        #     np.full(end_id - start_id + 1, np.nan),
        #     columns=["z_normalized_distance"],
        #     chunksize=1000000
        # )
        # new_df.index = dd.from_array(np.arange(start_id, end_id + 1))
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
    combined_df = combined_df.repartition(npartitions=n_jobs)
    print(combined_df.npartitions)
    # 计算均值
    meanA = combined_df.iloc[:, :m1].mean(axis=1)
    meanB = combined_df.iloc[:, m1:m1+m2].mean(axis=1)
    concerned_indi = (meanA < max_3d_dist) | (meanB < max_3d_dist)

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
    # FDR 校正可以在最后进行

    # 保存结果
    result_df.to_parquet(fo)
    client.close()
    return result_df