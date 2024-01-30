import dask
import dask.array as da
import dask.dataframe as dd
import numpy as np
import pandas as pd
from hic_basic.binnify import GenomeIdeograph
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
def test_chunk(chunk, m1, m2):
    """
    Input:
        chunk: a chunk of combined_df, pd.DataFrame
        m1: number of samples in groupA
        m2: number of samples in groupB
    Output:
        p-value of mannwhitneyu test
    """
    pvalue = mannwhitneyu(chunk.iloc[:, :m1], chunk.iloc[:, m1:m1+m2], alternative="two-sided").pvalue
    return pd.Series(pvalue, index=chunk.index)
def N_partitions(max_dist=2e6, binsize=20e3, chrom_mean_length=100e6, chunksize=100000):
    """
    Pick a number of partitions suitable for the data.
    Input:
        max_dist: max distance in impuation
        binsize: resolution of the data
        chrom_mean_lengths: mean length of chromosomes
        chunksize: number of rows suitable in memory, 100e3 ~ 6MB
    Output:
        number of partitions
    """
    nbands = max_dist // binsize
    nrows = nbands * (chrom_mean_length // binsize)
    if nrows < chunksize: # 1 if less than chunksize
        npartitions = 1
    else:
        npartitions = nrows // chunksize + 1 # at least 2
    return npartitions
def process_file(imputed_path, chrom, genome, max_dist, binsize, chunksize=10000):
    ideograph = GenomeIdeograph(genome)
    npartitions = N_partitions(
        max_dist,
        binsize,
        ideograph.chromosomes.loc[chrom, "length"],
        chunksize
        )
    # load data
    df = dd.read_parquet(imputed_path)
    # zscore normalization
    df = df.repartition(npartitions=npartitions)
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
    df = df.set_index("pixel_id", npartitions=npartitions)
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
    all_possible_pixels = all_possible_pixels.set_index("pixel_id", npartitions=npartitions)
    extended_df = dd.merge(
        all_possible_pixels,
        df,
        how="left",
        left_index=True,
        right_index=True
    )
    return extended_df["z_normalized_distance"]
def simpleDiff(groupA, groupB, chrom="chr1", genome="mm10", fo=None, max_3d_dist=5, 
               max_dist=2000000, binsize=20000, chunksize=10000):
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
    m1, m2 = len(groupA), len(groupB)

    # load all data
    dfs = []
    for imputed_path in groupA + groupB:
        dfs.append(process_file(imputed_path, chrom, genome, max_dist, binsize, chunksize))
    combined_df = dd.concat(dfs, axis=1)
    combined_df.columns = [f"sample{i}" for i in range(m1 + m2)]
    #combined_df.to_parquet("combined_df.parquet")
    # combined_df = da.stack(
    #     dd.compute(*dfs),
    #     axis=1
    # )
    # combined_df = dd.from_dask_array(combined_df, columns=[f"sample{i}" for i in range(m1 + m2)])
    print("combined_df.npartitions:", combined_df.npartitions)

    # calculate
    meanA = combined_df.iloc[:, :m1].mean(axis=1)
    meanB = combined_df.iloc[:, m1:m1+m2].mean(axis=1)
    concerned_indi = ((meanA < max_3d_dist) | (meanB < max_3d_dist)) & (meanA.notnull() & meanB.notnull())
    combined_df = combined_df.loc[concerned_indi]
    pvalues = combined_df.map_partitions( # return a series
        test_chunk,
        m1,
        m2,
        meta=('pvalue', 'f8')
    )
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
    return result_df