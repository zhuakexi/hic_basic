import dask.dataframe as dd
import numpy as np
from dask import delayed
from hic_basic.simpute import cis_distance_graph
from scipy.stats import zscore, mannwhitneyu

def normalize_band(df):
    """
    Input:
        df: bedpe like dataframe, with columns chrom1, start1, chrom2, start2, value
    Output:
        df: same as input, with value normalized by band
    """
    df['band'] = df['start2'] - df['start1']
    df['z_normalized_value'] = df.groupby(['chrom1', 'chrom2', 'band'])['value'].transform(lambda x: zscore(x, ddof=1))
    df = df.drop(columns=['band', 'value'])
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

def simpleDiff_test(df, m1, m2, top_threds=0.2):
    """
    Input:
        df: n * (m1 + m2) dask dataframe, n is pixels to be tested, m1 and m2 are number of samples in each group
    Output:
        df: n * 1 dask dataframe, pvalue
    """
    # 计算统计数据，但不执行实际计算
    stats_values = df.apply(lambda row: calculate_stats(row, m1, m2), axis=1, meta=('x', object))

    # 将结果转换为 DataFrame
    stats_df = stats_values.to_frame()
    stats_df.columns = ['pvalue']

    return stats_df
def simpleDiff(groupA, groupB, fo, max_3d_dist=5, max_dist=2000000, binsize=20000, n_jobs=16):
    """
    Input:
        groupA: _3dg_file paths for groupA
        groupB: _3dg_file paths for groupB
        fo: output parquet file
        n_jobs: number of jobs
    Output:
        df: n * (6+3) dask dataframe
            (chrom1, start1, end1, chrom2, start2, end2, meanA, meanB, diff, pvalue, qvalue)
    """
    m1, m2 = len(groupA), len(groupB)

    # 读取和处理 groupA 和 groupB 的文件
    dfs = []
    for _3dg_path in groupA + groupB:
        df = cis_distance_graph(_3dg_path, fo=None, max_dist=max_dist, binsize=binsize, n_jobs=1)
        df = normalize_band(df)
        df.drop
        dfs.append(df)

    # 使用 reduce 进行多个 DataFrame 的合并
    from functools import reduce
    def merge_dfs(df1, df2):
        return df1.merge(df2, on=['chrom1', 'start1', 'chrom2', 'start2'], how='outer')
    combined_df = reduce(merge_dfs, dfs)

    # keep only concerned pixels
    id_df, data_df = combined_df[0].iloc[:, :4], combined_df[0].iloc[:, 4:]
    id_df["end1"] = id_df["start1"] + binsize
    id_df["end2"] = id_df["start2"] + binsize
    id_df = id_df[["chrom1", "start1", "end1", "chrom2", "start2", "end2"]]
    meanA = data_df.iloc[:, :m1].mean(axis=1)
    meanB = data_df.iloc[:, m1:].mean(axis=1)
    concerned_indi = (meanA < max_3d_dist) | (meanB < max_3d_dist)
    # 调用 simpleDiff_test 计算 p 值
    pvalues = simpleDiff_test(data_df.loc[concerned_indi], m1, m2)
    
    # tidy up
    result_df = dd.concat([id_df.loc[concerned_indi], pvalues], axis=1)
    result_df.columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'pvalue']
    result_df["meanA"] = meanA
    result_df["meanB"] = meanB
    result_df["diff"] = meanA - meanB
    result_df["qvalue"] = multiple_testing_correction(result_df["pvalue"], correction_type="FDR")

    # 保存结果
    result_df.to_parquet(fo)
    return result_df