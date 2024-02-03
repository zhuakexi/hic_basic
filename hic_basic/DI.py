from pathlib import Path

import dask.array as da
import dask.dataframe as dd
import numpy as np
import pandas as pd
from dask import delayed
from hic_basic.binnify import GenomeIdeograph
from hic_basic.plot.hic import _plot_mat
from plotly.subplots import make_subplots
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
def mannwhitneyu_row(row, m1, m2):
    group1 = row[:m1]
    group2 = row[m1:m1+m2]
    _, p_value = mannwhitneyu(group1, group2, alternative='two-sided', use_continuity=False)
    return p_value
# def test_chunk(chunk, m1, m2):
#     """
#     Input:
#         chunk: a chunk of combined_df, pd.DataFrame
#         m1: number of samples in groupA
#         m2: number of samples in groupB
#     Output:
#         p-value of mannwhitneyu test
#     """
#     pvalue = mannwhitneyu(chunk.iloc[:, :m1], chunk.iloc[:, m1:m1+m2], alternative="two-sided").pvalue
#     return pd.Series(pvalue, index=chunk.index)
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
    print("npartitions:", npartitions)
    # load data
    df = dd.read_parquet(imputed_path)
    df = df.repartition(npartitions=npartitions)
    # filter out linear-far pixels
    df = df[(df["start1"] - df["start2"]).abs() < max_dist]
    # zscore normalization
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
    # print("simpleDiff")
    # print(df.npartitions)
    # print(df.compute())
    # print(df.partitions)
    # --- add pixel_id ---
    count_l1 = (ideograph.chromosomes["length"] // binsize).to_dict()
    cumul_l2 = ((ideograph.chromosomes["length"] // binsize) ** 2).cumsum()
    offsets_l2 = pd.Series(np.insert(cumul_l2.values, 0, 0)[:-1], index=cumul_l2.index).to_dict()
    # 计算 pixel_id
    df["pixel_id"] = (df['chrom1'].map(offsets_l2).astype(int)
                    + (df['start1'] // binsize) * (df['chrom1'].map(count_l1).astype(int))
                    +(df['start2'] // binsize))
    #print(df.compute())
    # df = df.map_partitions(
    #     ideograph.join_pixel_id,
    #     binsize=binsize,
    #     filep=imputed_path
    # )
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
    print("all_possible_pixels.npartitions:", all_possible_pixels.npartitions)
    possible_pixel_partitions = N_partitions(
        ideograph.chromosomes.loc[chrom, "length"], # pixel id is n**2 built, change to bandi version if needed
        binsize,
        ideograph.chromosomes.loc[chrom, "length"],
        chunksize
        )
    all_possible_pixels = all_possible_pixels.set_index("pixel_id", npartitions=possible_pixel_partitions)
    extended_df = dd.merge(
        all_possible_pixels,
        df,
        how="left",
        left_index=True,
        right_index=True
    )
    #print(extended_df.compute())
    return extended_df["z_normalized_distance"]
def simpleDiff(groupA, groupB, chrom="chr1", genome="mm10", fo=None, max_3d_dist=5, 
               max_dist=2000000, binsize=20000, chunksize=10000, cache=None):
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
    if cache is not None:
        combined_df.to_parquet(cache)

    # calculate
    print("start calculate")
    meanA = combined_df.iloc[:, :m1].mean(axis=1)
    meanB = combined_df.iloc[:, m1:m1+m2].mean(axis=1)
    concerned_indi = ((meanA < max_3d_dist) | (meanB < max_3d_dist)) & (meanA.notnull() & meanB.notnull())
    combined_df = combined_df.loc[concerned_indi]
    # pvalues = combined_df.map_partitions( # return a series
    #     test_chunk,
    #     m1,
    #     m2,
    #     meta=('pvalue', 'f8')
    # )
    pvalues = combined_df.apply(
        mannwhitneyu_row,
        axis=1,
        args=(m1, m2),
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
def imputed_reader(imputed_path, sample_id, genome, binsize, max_dist):
    # --- load data ---
    df = pd.read_parquet(imputed_path)
    
    # --- zscore normalization ---
    value_col_name = df.columns[-1]
    df['band'] = (df['start2'] - df['start1']) // binsize
    df[f'z_normalized_{value_col_name}'] = df.groupby("band")[value_col_name].transform(
        lambda x: (x - x.mean()) / x.std()
    )
    
    # --- add pixel_id ---
    ideograph = GenomeIdeograph(genome)
    count_l1 = (ideograph.chromosomes["length"] // binsize).to_dict()
    cumul_l2 = ((ideograph.chromosomes["length"] // binsize) ** 2).cumsum()
    offsets_l2 = pd.Series(np.insert(cumul_l2.values, 0, 0)[:-1], index=cumul_l2.index).to_dict()
    df["pixel_id"] = (df['chrom1'].map(offsets_l2).astype(int)
                    + (df['start1'] // binsize) * (df['chrom1'].map(count_l1).astype(int))
                    +(df['start2'] // binsize))
    
    # --- prepare output ---
    df = df.assign(
        sample_id = sample_id
    )
    # filter out linear-far pixels
    df = df.loc[(df["start1"] - df["start2"]).abs() < max_dist]
    return df[["pixel_id", "sample_id", f'z_normalized_{value_col_name}']]
# def simpleDiff_pivot(groupA, groupB, chrom="chr1", genome="mm10", fo=None, max_3d_dist=5, 
#                max_dist=2000000, binsize=20000, chunksize=10000, cache=None):
#     """
#     Input:
#         groupA: imputed .parquet file path for groupA
#         groupB: imputed .parquet file paths for groupB
#         fo: output parquet file
#         filt_fdr: filter pixels with fdr > filt_fdr, if None, no filtering
#         n_jobs: number of jobs
#     Output:
#         df: n * (6+3) dask dataframe
#             (chrom1, start1, end1, chrom2, start2, end2, meanA, meanB, diff, pvalue, qvalue)
#     """
#     m1, m2 = len(groupA), len(groupB)

#     # load all data
#     dfs = [
#         delayed(imputed_reader)(imputed_path, sample_id, genome, binsize, max_dist)
#         for sample_id, imputed_path in enumerate(groupA + groupB)
#         ]
#     combined_df = dd.from_delayed(
#         dfs,
#         meta = { "pixel_id": "i8", "sample_id": "i8", "z_normalized_distance": "f8" }
#         )
#     combined_df["sample_id"] = combined_df["sample_id"].astype("category")
#     combined_df["sample_id"] = combined_df["sample_id"].cat.set_categories(
#         list(range(m1 + m2))
#     )
#     combined_df["sample_id"] = combined_df["sample_id"].cat.as_known()
#     combined_df = combined_df.pivot_table(
#             index="pixel_id",
#             columns="sample_id",
#             values="z_normalized_distance",
#             aggfunc="sum"
#     )
#     combined_df.columns = [f"sample{i}" for i in range(m1 + m2)]
#     print("combined_df.npartitions:", combined_df.npartitions)
 
#     if cache is not None:
#         combined_df.to_parquet(cache)

#     # calculate
#     print("start calculate")
#     meanA = combined_df.iloc[:, :m1].mean(axis=1)
#     meanB = combined_df.iloc[:, m1:m1+m2].mean(axis=1)
#     concerned_indi = ((meanA < max_3d_dist) | (meanB < max_3d_dist)) & (meanA.notnull() & meanB.notnull())
#     combined_df = combined_df.loc[concerned_indi]
#     pvalues = combined_df.apply(
#         mannwhitneyu_row,
#         axis=1,
#         args=(m1, m2),
#         meta=('pvalue', 'f8')
#     )
#     # meanA, meanB, diff, pvalues
#     result_df = dd.concat(
#         [
#             meanA,
#             meanB,
#             meanA - meanB,
#             pvalues
#         ],
#         axis=1
#     )
#     result_df.columns = ["meanA", "meanB", "diff", "pvalue"]

#     # 保存结果
#     result_df.to_parquet(fo)
#     return result_df
def read_pvalues(pval_path, chrom):
    df = pd.read_parquet(pval_path)
    df = df.assign(chrom=chrom)
    return df
def simpleDiff_postprocess(pvalue_files, chroms, genome="mm10", binsize=20000, filt_fdr=0.05, topDI=0.2, fo=None):
    """
    Input:
        pvalue_file: pvalue file path
            meanA, meanB, diff, pvalue
        filt_fdr: filter pixels with fdr > filt_fdr
            if None, no filtering
        topDI: mark abs diff > topDI as topDI
            if topDI is None, don't mark
        fo: output parquet file
    Output:
        df: n * (6+3) dask dataframe
            (chrom1, start1, end1, chrom2, start2, end2, meanA, meanB, diff, pvalue, qvalue)
    """
    # load all data
    df = dd.from_delayed([
        delayed(read_pvalues)(pval_path, chrom)
        for chrom, pval_path in zip(chroms, pvalue_files)
        ])
    df = df.reset_index()
    df = df.compute()
    # multiple testing correction
    # return multiple_testing_correction(df["pvalue"], correction_type="FDR")[:]
    df = df.assign(qvalue = multiple_testing_correction(df["pvalue"], correction_type="FDR"))
    # filter
    if filt_fdr is not None:
        df = df[df["qvalue"] < filt_fdr]
    if topDI is not None:
        df = df.assign(topDI = df["diff"].abs() > topDI)
    # transform to bedpe
    ideograph = GenomeIdeograph(genome)
    df = ideograph.join_positions(df, binsize=binsize)
    # 保存结果
    if fo is not None:
        df.to_parquet(fo)
    return df
def project_DI(genome, chrom, start, end, binsize, DI):
    ideograph = GenomeIdeograph(genome)
    # generate na mat of target region
    pos = ideograph.bins(binsize, bed=True).query(
        'chrom == @chrom and start >= @start and start <= @end'
        )["start"].rename('pos')
    mat = pd.DataFrame(
        index=pos,
        columns=pos,
        data=np.nan
    )
    # fill with DI diff
    mat = mat.fillna(DI.set_index(["start1", "start2"])["diff"].unstack())
    return mat
def pileup_impute(imputes, genome, binsize, chrom, start, end,zscore=False,symmetric=True):
    ideograph = GenomeIdeograph(genome)
    # generate na mat of target region
    pos = ideograph.bins(binsize, bed=True).query(
        'chrom == @chrom and start >= @start and start <= @end'
        )["start"].rename('pos')
    pileups = pd.DataFrame(
        index=pos,
        columns=pos,
        data=0
    )
    valid_samples = 0
    for i, imputed in enumerate(imputes):
        imputed = str(imputed).format(chrom=chrom)
        if not Path(imputed).exists():
            print("Warning: file not found", imputed)
            continue
        valid_samples += 1
        mat = pd.DataFrame(
            index=pos,
            columns=pos,
            data=np.nan
        )
        pixels = pd.read_parquet(imputed)
        if zscore:
            # zscore normalization
            value_col_name = pixels.columns[-1]
            pixels['band'] = (pixels['start2'] - pixels['start1']) // binsize
            pixels[f'z_normalized_{value_col_name}'] = pixels.groupby("band")[value_col_name].transform(
                lambda x: (x - x.mean()) / x.std()
            )
            pixels = pixels[["start1","start2","z_normalized_distance"]]
            pixels.columns = ["start1", "start2", "distance"]
        else:
            pass
        mat = mat.fillna(pixels.set_index(["start1", "start2"])["distance"].unstack())
        mat = mat.fillna(0)
        pileups += mat
        if i % 100 == 0:
            print(i,"done.")
    # calculate average
    pileups = pileups / valid_samples
    if symmetric:
        pileups = pileups + pileups.T
    return pileups

#   --- plot ---
def plot_DI(genome, chrom, start, end, binsize, DI):
    mat = project_DI(genome, chrom, start, end, binsize, DI)
    fig = _plot_mat(
        mat,
        vmax=1,
        donorm=False,
        zmin=-1,
        cmap = "RdBu_r",
    )
    fig.update_layout(
        height = 500,
        width = 500
    )
    return fig
    #fig.show(renderer="png")
def plot_compare_impute(imputesA, imputesB, genome, binsize, chrom, start, end, zscore=True, subplot_titles=("A", "B", "A-B")):
    mats = []
    for imputes in [imputesA, imputesB]:
        mats.append(
            pileup_impute(
                imputes,
                genome,
                binsize,
                chrom,
                start,
                end,
                zscore=zscore
            )
        )
    mats.append(mats[0] - mats[1])
    figure = make_subplots(
        rows=1, cols=3,
        subplot_titles=subplot_titles
    )
    for i, mat in enumerate(mats):
        figure.add_trace(
            _plot_mat(
                mat,
                donorm = False,
                cmap="RdBu_r",
                zmin = -0.5,
                zmax = 0.5
            ).data[0],
            row=1, col=i+1
        )
    figure.update_layout(
        margin = dict(l=0, r=0, t=30, b=0),
        height = 325,
        width = 1000
    )
    return figure
chrom, start,end = "chr1", 135e6, 155e6