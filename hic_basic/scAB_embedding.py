# embedding cell using spacing CpG
# see Tan, Longzhi, Wenping Ma, Honggui Wu, Yinghui Zheng, Dong Xing, Ritchie Chen, Xiang Li, Nicholas Daley, Karl Deisseroth, and X. Sunney Xie. 
# “Changes in Genome Architecture and Transcriptional Dynamics Progress Independently of Sensory Experience during Post-Natal Brain Development.” 
# Cell 184, no. 3 (February 2021): 741-758.e17. https://doi.org/10.1016/j.cell.2020.12.032.

import gzip
from concurrent import futures
from itertools import dropwhile, repeat
from pathlib import Path

import numpy as np
import pandas as pd
import toolz
from sklearn import preprocessing
from sklearn.decomposition import PCA
from sklearn.neighbors import radius_neighbors_graph
from umap import UMAP

def get_bin_locus(pos,size):
    return int( round( float(pos) / size) ) * size
def color2(pairsf,color_file,bin_size,merge_haplotypes=True,dropXY=True):
    """
    Calculate 3D CpG density of each chromosome bin
    Input:
        pairsf: 4DN's pairs file
        color_file: reference linear CpG density
        bin_size: binsize of color_file(and the result output)
        merge_haplotypes: whether treat maternal paternal chromosome differently,
            when False, chromosomes must be like "chr1(mat)"
        dropXY: only output autosome result
    """
    # read in data
    pairs = []
    with gzip.open(pairsf,"rt") as f:
        for line in dropwhile(lambda x:x.startswith("#"),f):
            pairs.append(line.strip().split("\t"))
    if isinstance(color_file,str):
        color_data = {}
        with open(color_file,"rt") as f:
            for color_file_line in f:
                hom_name, ref_locus, color = color_file_line.strip().split("\t")
                ref_locus = int(ref_locus)
                color = float(color)
                color_data[(hom_name, ref_locus)] = color
    elif isinstance(color_file, pd.DataFrame):
        color_data = color_file.set_index(["chrom","start"])["CpG"]
    else:
        raise ValueError("color_file must be str or pd.DataFrame")
    # smoothing
    smooth_color_data = {}
    ## bin locus
    for con in pairs:
        bin_locus1 = get_bin_locus(con[2],bin_size)
        bin_locus2 = get_bin_locus(con[4],bin_size)
        if merge_haplotypes == True:
            # omit pairs-file's phasing col and dip-file's phased chrom name
            hom_name1 = con[1].split("(")[0]
            hom_name2 = con[3].split("(")[0]
        else:
            # omit pairs-file's phasing col but keep dip-file's phased chrom name
            hom_name1 = con[1]
            hom_name2 = con[3]
        bin_leg1 = (hom_name1, bin_locus1)
        bin_leg2 = (hom_name2, bin_locus2)
    ## recording con-color
        if bin_leg1 == bin_leg2:
            continue
        if bin_leg1 in color_data:
            if bin_leg2 not in smooth_color_data:
                smooth_color_data[bin_leg2] = []
            smooth_color_data[bin_leg2].append(color_data[bin_leg1])
        if bin_leg2 in color_data:
            if bin_leg1 not in smooth_color_data:
                smooth_color_data[bin_leg1] = []
            smooth_color_data[bin_leg1].append(color_data[bin_leg2])
    for bin_leg in smooth_color_data:
        smooth_color_data[bin_leg] = np.mean(smooth_color_data[bin_leg])
    if dropXY == True:
            smooth_color_data = {k:v for k,v in smooth_color_data.items()
                                 if k[0] not in ["chrX","chrY","chrX(mat)","chrX(pat)","chrY(mat)","chrY(pat)"]}
    return smooth_color_data
def s_color2(_3dg:pd.DataFrame, color_file, min_dist=3, n_jobs=12)->pd.DataFrame:
    """
    Calculate intermingling metrics for each particle in 3dg file.
    Input:
        _3dg: dataframe, result of hic_basic.hicio parse_3dg
        color_file: reference linear CpG density
        min_dist: int, minimum distance to be considered as neighbor
        n_jobs: int, number of jobs to run in parallel
    Output:
        dataframe with same index, and a scAB column
    """
    # --- load data ---
    if isinstance(color_file, pd.DataFrame):
        color_data = color_file.set_index(["chrom","start"])["CpG"]
    elif isinstance(color_file,str) or isinstance(color_file,Path):
        color_data = pd.read_table(
            color_file,names=["chrom","start","CpG"],
            index_col=["chrom","start"]
            )[["CpG"]]
    _3dg = _3dg.copy()
    _3dg_e = pd.concat(
        [
            _3dg,
            color_data
        ],
        axis=1,
        join="outer"
    ).loc[_3dg.index]
    del color_data
    # --- generate radius neighbor ---
    # get sparse matrix, basically tuple of two arrays, (row, col)
    RN = radius_neighbors_graph(
        _3dg_e[["x","y","z"]].values,
        min_dist,
        mode = "distance",
        metric = "minkowski",
        p = 2,
        include_self = False,
        n_jobs=n_jobs
        ).nonzero()
    # --- count RN chromosome distribution for each particle ---
    # seq
    #for_CpG_count = ((i, _3dg_e.iat[j, 3]) for i, j in zip(*RN))
    for_CpG_count = pd.DataFrame(
        {
            "index":RN[0],
            "neighbor_CpG":_3dg_e.iloc[RN[1]]["CpG"].values
        }
    ).values.tolist()
    del RN
    #print(for_CpG_count[200:300])
    # init
    # (CpG_sum, n_neighbors)
    init_CpG_count = (0,0)
    # binop
    def CpG_stat(acc, x):
        """
        Add value to CpG sum and neighbor count.
        Input:
            acc: original stat, 2-tuple, (CpG_sum, n_neighbors)
            x: input_line, (index, neighbor CpG value)
        Ouput:
            new_acc: updated stat
        Note: if nan, will be ignored
        """
        if np.isnan(x[1]):
            return acc
        else:
            return (acc[0]+x[1], acc[1]+1)
    # start reduce
    per_bin_count = toolz.reduceby(
        lambda x: x[0], # key
        CpG_stat,
        for_CpG_count,
        init_CpG_count
    )
    del for_CpG_count
    # --- generate frequency table ---
    CpG_count_table = pd.DataFrame(
        per_bin_count,
        index=["CpG_sum","n_neighbors"]
    ).T
    #print(CpG_count_table)
    # tidy up index
    CpG_count_table = pd.concat(
        [_3dg.reset_index(), CpG_count_table],
        axis=1,
        join="outer"
    )
    del _3dg_e
    # --- calculate meaningful metrics ---
    CpG_count_table["scAB"] = CpG_count_table["CpG_sum"] / CpG_count_table["n_neighbors"]
    return CpG_count_table
def stack_dict(ares, sample_name:None, col_thresh=0.9, row_thresh=0.9):
    """
    # stack list of dict, drop bad rows and columns
    # first clean bad cols, then clean bad rows
    # Input:
    #    ares: list of dict
    #    sample_name: name list of samples
    #    col_thresh: [0,1] larger is stricter, 0 to keep all
    #        keeping cols that at leat *ratio* of samples have nonNA value
    #    row_thresh: [0,1] larger is stricter, 0 to keep all
    #        after col cleaning, keeping samples that have nonNA value for at leat *ratio* of all attrs
    """
    new_data = pd.DataFrame(ares,index=sample_name)
    # clean col
    new_data.dropna(axis=1, thresh=np.ceil(new_data.shape[0]*col_thresh), inplace=True)
    # clean row
    new_data.dropna(axis=0, thresh=np.ceil(new_data.shape[1]*row_thresh), inplace=True)
    return new_data
def fill_color(data, color_file, grt_full=True):
    """
    # Fill in NA using reference CpG value
    # This make sense because spacing CpG is highly correlated with Sequencing CpG
    # Input:
    #    grt_full: guarantee that output dataframe is NA free,
    #        this will drop all col's that can't be filled by ref color_file
    """
    if isinstance(color_file,str):
        ref_color = pd.read_table(color_file,names=["chrom","pos","CpG"])
    elif isinstance(color_file, pd.DataFrame):
        ref_color = color_file[["chrom","start","CpG"]]
        ref_color.columns = ["chrom","pos","CpG"]
    else:
        raise ValueError("color_file must be str or pd.DataFrame")
    ref_color.set_index(["chrom","pos"],inplace=True)
    ref_color.index = ref_color.index.to_flat_index()
    # filling data with ref CpG file
    filled_data = data.fillna(value=ref_color["CpG"])
    # remove columns if ref doesn't have that locus
    if grt_full:
        filled_data = filled_data.dropna(axis=1)
    return filled_data
def calc_color2(filesp, file_col, color_file, binsize, merge_haplotypes=True, dropXY=True, col_thresh=0.9, row_thresh=0.9,fill_color=True, threads=24):
    """
    Input:
        filesp: dataframe, must have file_col col
        file_col: column in filesp that stores pairs file path
        color_file: ref bed file that stores linear CpG density
        binsize: binsize of color_file
        merge_haplotypes: whether treat maternal paternal chromosome differently,
            when False, chromosomes must be like "chr1(mat)"
        dropXY: only output autosome result
        col_thresh: [0,1] larger is stricter, 0 to keep all
            keeping cols that at leat *ratio* of samples have nonNA value
        row_thresh: [0,1] larger is stricter, 0 to keep all
            after col cleaning, keeping samples that have nonNA value for at leat *ratio* of all attrs
        threads: number of cores using
    Note:
        when binsize is small, valid col and row thresh drop dramatically, it's better to use color3 at that time
    """
    print("Calculating color...")
    with futures.ProcessPoolExecutor(threads) as pool:
        res = pool.map(
            color2,
            filesp[file_col],
            repeat(color_file,filesp.shape[0]),
            repeat(binsize, filesp.shape[0]),
            repeat(merge_haplotypes, filesp.shape[0]),
            repeat(dropXY, filesp.shape[0])
        )
    ares = list(res)
    print("Stacking color...")
    color_result = stack_dict(ares, filesp.index, col_thresh, row_thresh)
    if fill_color:
        print("Filling color...")
        color_result = fill_color(color_result, color_file)
    else:
        pass
    return color_result
def do_umap(data,ndims=30):
    # do exactly what Seurat RunUMAP do
    # seens to plot better than umap-learn default settings
    umap = UMAP(
        #n_neighbors = 30,
        n_neighbors=20,
        n_components=2,
        metric = "cosine",
        n_epochs=None,
        learning_rate = 1,
        #min_dist=0.3,
        min_dist = 0.0001,
        spread=1,
        set_op_mix_ratio=1,
        local_connectivity=1,
        repulsion_strength=1,
        negative_sample_rate=5,
        a = None,
        b = None,
        random_state=42,
        angular_rp_forest=False,
        densmap=False,
        dens_lambda=2,
        dens_frac=0.3,
        dens_var_shift=0.1,
        output_dens=False
    )
    umap_res = umap.fit_transform(data[:,:ndims])
    return umap_res
def scAB_embedding(data,ndims=30):
    rank_normed = (data.rank() - 1)/data.shape[0]
    scaler = preprocessing.StandardScaler(with_mean=True,with_std=False)
    scaled = scaler.fit_transform(rank_normed.values)
    pca = PCA(random_state=42)
    pca_res = pca.fit_transform(scaled)
    # keep first 2 u by default
    umap_res = do_umap(pca_res,ndims)
    result = pd.DataFrame(umap_res, index=data.index)
    result.columns = ["u1","u2"]
    return result