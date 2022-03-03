import anndata as ad
import pandas as pd
import numpy as np
from scipy import sparse
def fix_velocyto_names(adata, annote, verbose=False):
    # ---assining normname to new col---
    annote = pd.read_csv(annote,index_col=0,dtype={"sample_name":"string"})
    new_names = np.array(adata.obs.index)
    for norm_name in annote.index:
        # target CellID has norm_name substring
        s = adata.obs.index.str.contains(norm_name)
        # set all target sample name to norm_name
        # different target sample name may have same norm_name
        if s.sum() > 0:
            new_names[s] = norm_name
    adata.obs = adata.obs.assign(normname=new_names)
    
    # ---Mark name fix res---
    # annote has these samples(with same normname) multiple times
    dup_samples = adata.obs.loc[
        adata.obs["normname"].isin(
            adata.obs.groupby("normname").size()[(adata.obs.groupby("normname").size() > 1)].index
        )
    ]
    # annote has these samples that adata doesn't
    missing_samples = list(annote.index[~annote.index.isin(adata.obs["normname"])])
    # adata has these samples that annote doesn't
    black_samples = adata.obs[~adata.obs["normname"].isin(annote.index)]
    adata.obs = adata.obs.assign(fix_name_res="good")
    adata.obs.loc[dup_samples.index, "fix_name_res"] = "Dup"
    adata.obs.loc[black_samples.index, "fix_name_res"] = "Illegal"
    
    # ---Print Summary---
    if verbose == True:
        print("---Dup samples---\n", dup_samples)
        print("---Missing samples---\n", missing_samples)
        print("---Illegal samples---\n", black_samples)
    return adata
def clean_velocyto_names(adata, using_dup:list):
    if len(set(using_dup)) < len(using_dup):
        raise ValueError("clean_velocyto_names: using_dup list must be unique.")
    for i in using_dup:
        if i not in adata.obs[adata.obs["fix_name_res"] == "Dup"].index:
            raise ValueError("Error "+ i + " is not marked as Dup.")
    I = adata.obs.index[adata.obs["fix_name_res"]=="good"]
    I.append(pd.Index(using_dup))
    adata = adata[I]
    adata.obs.index.name = "velocyto_CellID"
    adata.obs = adata.obs.set_index("normname")
    adata.obs.index.name = "sample_name"
    adata.obs = adata.obs.drop("fix_name_res",axis=1)
    return adata
def expand_df(target_df, ref_df):
    """
    expand target_df to shape(with same index, columns) of ref_df, fill with 0
    Input:
        target_df: df to be transformed; sparse df
        ref_df: axis reference, must contain all rows and cols of target_df; sparse df
    """
    assert target_df.index.isin(ref_df.index).all() and target_df.columns.isin(ref_df.columns).all(),\
        "expand_df: target_df axis aren't both subset of ref_df"
    cofactor = ref_df.loc[
        ~ref_df.index.isin(target_df.index),
        ~ref_df.columns.isin(target_df.columns)
    ]
    # for cofactor, join on both axis are OK
    target_df = pd.concat(
        [target_df, cofactor],
        join = "outer"
    )
    target_df.fillna(0, inplace=True)
    target_df = target_df.astype(pd.SparseDtype(int, fill_value=0))
    return target_df
def create_adata(expr, velo_ad, g1, g2):
    """
    Merge different matrix to single AnnData object, intersection for obs, union for features.
    Input:
        expr: RNA count matrix from normal RNA-seq pipeline
        velo_ad: .loom from velocyto pipeline
        g1: genome 1 count matrix from SNPsplit
        g2: genome 2 count matrix from SNPsplit
    Output:
        adata, multiple-layered anndata.AnnData
    """
    # ---tidy up inputs---
    # using velo_ad transpose
    velo_ad.var_names_make_unique()
    velo_ad = velo_ad.T
    expr.columns = expr.columns.astype("string")
    # ---get common samples and genes---
    # expr and velo_ad are all (var, obs) here for in this study sample number is 
    #   much less than gene number
    g_samples = g1.columns.intersection(g2.columns)
    g_genes = g1.index.intersection(g2.index)
    samples = g_samples.intersection(expr.columns.intersection(
            velo_ad.var.index
    ))
    genes = g_genes.union(expr.index.union(velo_ad.obs.index))
    print(len(genes),len(samples))
    # ---generate seperate dfs--- 
    matrix = pd.DataFrame.sparse.from_spmatrix(velo_ad.layers["matrix"], index = velo_ad.obs.index, columns = velo_ad.var.index)
    spliced = pd.DataFrame.sparse.from_spmatrix(velo_ad.layers["spliced"], index = velo_ad.obs.index, columns = velo_ad.var.index)
    unspliced = pd.DataFrame.sparse.from_spmatrix(velo_ad.layers["unspliced"], index = velo_ad.obs.index, columns = velo_ad.var.index)
    ambiguous = pd.DataFrame.sparse.from_spmatrix(velo_ad.layers["ambiguous"], index = velo_ad.obs.index, columns = velo_ad.var.index)
    # future .X template
    pixels = pd.DataFrame.sparse.from_spmatrix(sparse.csr_matrix((len(genes),len(samples)),dtype=np.int64),index = genes, columns = samples)
    # ---prepare separate dfs---
    matrix = expand_df(matrix[samples], pixels)
    spliced = expand_df(spliced[samples], pixels)
    unspliced = expand_df(unspliced[samples], pixels)
    ambiguous = expand_df(ambiguous[samples], pixels)
    expr = expand_df(expr[samples], pixels)
    g1 = expand_df(g1[samples], pixels)
    g2 = expand_df(g2[samples], pixels)
    print(pixels.shape)
    for i in [matrix, spliced, unspliced, ambiguous, expr, g1, g2]:
        print(i.shape)
    # ---create annData object---
    new_ad = ad.AnnData(
        expr.loc[genes, samples].T,
        obs = pd.DataFrame(index=samples),
        var = pd.DataFrame(index=genes),
        layers = {
            "matrix":matrix.loc[genes, samples].T,
            "spliced":spliced.loc[genes, samples].T,
            "unspliced":unspliced.loc[genes, samples].T,
            "ambiguous":ambiguous.loc[genes, samples].T,
            "g1":g1.loc[genes, samples].T,
            "g2":g2.loc[genes, samples].T
        }
    )
    return new_ad