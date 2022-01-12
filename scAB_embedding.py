# embedding cell using spacing CpG
# see Tan, Longzhi, Wenping Ma, Honggui Wu, Yinghui Zheng, Dong Xing, Ritchie Chen, Xiang Li, Nicholas Daley, Karl Deisseroth, and X. Sunney Xie. 
# “Changes in Genome Architecture and Transcriptional Dynamics Progress Independently of Sensory Experience during Post-Natal Brain Development.” 
# Cell 184, no. 3 (February 2021): 741-758.e17. https://doi.org/10.1016/j.cell.2020.12.032.

import gzip
from itertools import dropwhile
import numpy as np
from sklearn import preprocessing
from sklearn.decomposition import PCA
from umap import UMAP
import pandas as pd
from concurrent import futures
from itertools import repeat

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
    color_data = {}
    with open(color_file,"rt") as f:
        for color_file_line in f:
            hom_name, ref_locus, color = color_file_line.strip().split("\t")
            ref_locus = int(ref_locus)
            color = float(color)
            color_data[(hom_name, ref_locus)] = color
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
    ref_color = pd.read_table(color_file,names=["chrom","pos","CpG"])
    ref_color.set_index(["chrom","pos"],inplace=True)
    ref_color.index = ref_color.index.to_flat_index()
    # filling data with ref CpG file
    filled_data = data.fillna(value=ref_color["CpG"])
    # remove columns if ref doesn't dave that locus
    if grt_full:
        filled_data = filled_data.dropna(axis=1)
    return filled_data
def calc_color2(filesp, file_col, color_file, binsize, merge_haplotypes=True, dropXY=True, col_thresh=0.9, row_thresh=0.9,threads=24):
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
    """
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
    color_result = stack_dict(ares, filesp.index)
    color_result = fill_color(color_result, color_file)
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