import anndata as ad
import pandas as pd
import numpy as np
from scipy import sparse
import os
import loompy
import scvelo as scv

from ..io import matr, read_meta
from .paracalc import gen_repli_score, gen_cdps, gen_PM_interactions
from .exp_record import add_cell_type, add_group_hour

def _merge_expr(fps, mapper, samplelist=None):
    """
    Merging umi count matrices.
    Input:
        fps: matrix file paths; list
        mapper: renaming samples(for ID fixing); dict
        samplelist: samples to keep, None for All; list
    Output:
        dataframe
    """
    samplelist = pd.Index(samplelist)
    all_mat = pd.DataFrame()
    for i in fps:
        sub_mat = matr(i,"\t")
        print("sub mats:",sub_mat.shape)
        # pd concat
        all_mat = pd.concat([all_mat, sub_mat],axis=1)
        all_mat.fillna(0, inplace=True)
    # ---rename sample id in expression matrix---
    all_mat.rename(columns=mapper,inplace=True)
    # ---check all_mat---
    if (~samplelist.isin(all_mat.columns)).sum() > 0:
        print("Warning: sample has no count data:")
        for i in samplelist[~samplelist.isin(all_mat.columns)]:
            print(i)
    print("raw merged matrix",all_mat.shape)
    # ---using only samples in samplelist---
    valid_cells = samplelist[samplelist.isin(all_mat.columns)]
    all_mat = all_mat.loc[:,valid_cells]
    # ---remove full 0 lines---
    all_mat = all_mat.loc[~(all_mat.sum(axis=1)==0),:]
    print("final merged matrix",all_mat.shape)
    return all_mat
def gen_expr(qc, outfp=None, mattail="count_matrix/counts.gene.tsv.gz"):
    """
    Generate single expression matrix that contains all (has valid task directory) qc passed samples from qc file.
    qc file must contain "task_dirp" col and use sample ID as index.
    fix sample ID with "_" for umi_tools limitation. 
    Input:
        qc: filtered dataframe
        outfp: if not None, output file rather than real dataframe
        mattail: RNA mat position according to task_dirp
    Output:
        merged expression matrix
    """
    return _merge_expr(
        [os.path.join(i,mattail) for i in qc["task_dirp"].unique()],
        mapper = dict(zip(qc.index[qc.index.str.contains("_")].str.split("_",expand=True).get_level_values(1),
                          qc.index[qc.index.str.contains("_")])),
        samplelist = qc.index.values
    )
def _trim_names(orig, ref, verbose=False):
    """
    Find standard id from ref. 
    e.g. 220110_embryo_8EZD5:2021123110 -> 2021123110. Use this to fix obs.index of velocyto output.
    Using verbose mode first to check dup, wrong, missing ids.
    Input:
        orig: pd.Index
        ref: pd.Index
    Output:
        dataframe orig as index, with normname col
    """
    new_names = np.array(orig)
    for norm_name in ref:
        # target CellID has norm_name substring
        s = orig.str.contains(norm_name)
        # set all target sample name to norm_name
        # different target sample name may have same norm_name
        if s.sum() > 0:
            new_names[s] = norm_name
    res = pd.DataFrame({"normname":new_names}, index=orig)
    # annote has these samples(with same normname) multiple times
    dup_samples = orig[
        res["normname"].isin(
            res.groupby("normname").size()[(res.groupby("normname").size() > 1)].index
        )
    ]
    # annote has these samples that adata doesn't
    missing_samples = list(ref[~ref.isin(res["normname"])])
    # adata has these samples that annote doesn't
    black_samples = res[~res["normname"].isin(ref)].index
    # ---Print Summary---
    if verbose == True:
        print("---Dup samples---\n", dup_samples)
        print("---Missing samples---\n", missing_samples)
        print("---Illegal samples---\n", black_samples)
        for i in black_samples:
            print("  ",i)
            
    # ---mark fixing result---
    res = res.assign(fix_name_res="good")
    res.loc[dup_samples, "fix_name_res"] = "Dup"
    res.loc[black_samples, "fix_name_res"] = "Illegal"
    return res
def fix_velocyto_names(adata, annote, verbose=False):
    """
    Fix obs.index of velocyto output.
    e.g. 220110_embryo_8EZD5:2021123110 -> 2021123110.
    """
    # ---assining normname to new col---
    annote = pd.read_csv(annote,index_col=0,dtype={"sample_name":"string"})
    res = _trim_names(adata.obs.index, annote.index, verbose)
    return res
def fix_velocyto_names_old(adata, annote, verbose=False):
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
# --- generate compartment strength ---
def gen_fileis(sample_table, dir_path, str_pat):
    """
    Check and generate fileis(input file tables storing file path) for downstream calculating.
    Input:
        sample_table: meta, sample_names as index, try to ensure all samples have their input file.
        dir_path: where to find those input files.
        str_pat: string pattern, gen exact file name. using {} to represent sample_name.
    Output:
        pd.DataFrame
    """
    if not sample_table.index.is_unique:
        raise ValueError("[gen_fileis] Error: sample_table index is not unique")
    a = dict(zip([str_pat.format(i) for i in sample_table.index], list(sample_table.index)))
    b = set(os.listdir(dir_path))
    hitting = a.keys() & b
    missing_sample = [a[k] for k in (a.keys() - hitting)]
    extra = list(b - hitting)
    
    if len(missing_sample):
        print("[gen_fileis] Warning: can't finds files for %d samples:" % len(missing_sample))
        print(",".join(missing_sample))
    if len(extra):
        print("[gen_fileis] Warning: find extra files for:")
        print(",".join(extra[:20]))
        if len(extra) > 20:
            print("...")
    return pd.Series([os.path.join(dir_path, k) for k in hitting], index=[a[k] for k in hitting])
def gen_cs(fileis):
    """
    Generate compartment strength table from pre-computed files.
    """
    cs = []
    for file in fileis.values:
        cs.append(pd.read_table(file, header=None))
    cs = pd.concat(cs,axis=1)
    cs.columns = fileis.index
    return cs.T
# --- create anndata object from scratch ---
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
    print("create_adata: Used samples %d" % len(samples))
    print("create_adata: Used genes %d" % len(genes))
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
# --- calculate all mid-files used for generating adata object ---
def gen_cache(qc, cache_dir, velo_files, g1_files, g2_files, cs_dir, threads=32):
    if not os.path.isdir(cache_dir):
        os.mkdir(cache_dir)
    
    print("Gen velo cache...")
    outfile = os.path.join(cache_dir, "velocyto.integrate.loom")
    if not os.path.isfile(outfile):
        loompy.combine(velo_files, outfile)
        velo = scv.read(outfile)
        a = _trim_names(pd.Index(velo.obs.index), qc.index)
        if "220110_embryo_8EZD5:2021123110" in a.index:
            a.loc["220110_embryo_8EZD5:2021123110","fix_name_res"] = "good"
        velo = velo[a.query('fix_name_res == "good"').index]
        velo.obs.index = a.query('fix_name_res == "good"')["normname"]
        velo.obs.index.name = "CellID"
        velo.write_loom(outfile)
    
    print("Gen g1 cache...")
    outfile = os.path.join(cache_dir, "g1.csv.gz")
    if not os.path.isfile(outfile):
        merged_g1 = _merge_expr(g1_files,  dict(zip(qc.index[qc.index.str.contains("_")].str.split("_",expand=True).get_level_values(1),
                              qc.index[qc.index.str.contains("_")])), qc.index.values)
        merged_g1.to_csv(outfile)
    
    print("Gen g2 cache...")
    outfile = os.path.join(cache_dir, "g2.csv.gz")
    if not os.path.isfile(outfile):
        merged_g2 = _merge_expr(g2_files,  dict(zip(qc.index[qc.index.str.contains("_")].str.split("_",expand=True).get_level_values(1),
                              qc.index[qc.index.str.contains("_")])), qc.index.values)
        merged_g2.to_csv(outfile)
    
    print("Gen cdps...")
    outfile = os.path.join(cache_dir, "cdps.csv.gz")
    if not os.path.isfile(outfile):
        cdps = gen_cdps(qc, threads)
        cdps.to_csv(outfile)
    
    print("Gen repli score...")
    outdir = os.path.join(cache_dir, "repliscore_cache",)
    outfile = os.path.join(cache_dir, "rs.csv.gz")
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    if not os.path.isfile(outfile):
        rs = gen_repli_score(qc, outdir)
        rs.to_csv(outfile)
    
    print("Gen PM interactions...")
    outfile = os.path.join(cache_dir, "PM_interaction.csv.gz")
    if not os.path.isfile(outfile):
        PM_interaction = gen_PM_interactions(qc)
        PM_interaction.index.name = "sample_name"
        PM_interaction.to_csv(outfile)
    
    print("Merging expression matrix...")
    outfile = os.path.join(cache_dir, "expr.csv.gz")
    if not os.path.isfile(outfile):
        mat = gen_expr(qc)
        mat.to_csv(outfile)
    
    # assuming g1 g2cs are precomputed and store in specific filename pattern.
    print("Gen g1 500k compartment strength...")
    outfile = os.path.join(cache_dir, "g1cs_500k.csv.gz")
    if not os.path.isfile(outfile):
        g1cs_fileis = gen_fileis(qc, cs_dir,  "{}.500000.cs.g1.txt")
        g1cs_500k = gen_cs(g1cs_fileis)
        g1cs_500k.to_csv(outfile)

    print("Gen g2 500k compartment strength...")
    outfile = os.path.join(cache_dir, "g2cs_500k.csv.gz")
    if not os.path.isfile(outfile):
        g2cs_fileis = gen_fileis(qc, cs_dir,  "{}.500000.cs.g2.txt")
        g2cs_500k = gen_cs(g2cs_fileis)
        g2cs_500k.to_csv(outfile)
    print("Done")
def gen_ad_from_cache(cache_dir, annote):
    print("reading tables...")
    expr = matr(os.path.join(cache_dir,"expr.csv.gz"))
    velo_ad = ad.read_loom(os.path.join(cache_dir, "velocyto.integrate.loom"))
    g1 = matr(os.path.join(cache_dir,"g1.csv.gz"))
    g2 = matr(os.path.join(cache_dir,"g2.csv.gz"))
    print("create anndata object...")
    res_ad = create_adata(expr, velo_ad, g1, g2)
    print("add per-cell annotations...")
    res_ad.obs = res_ad.obs.assign(g1_UMIs = res_ad.layers["g1"].sum(axis=1))
    res_ad.obs = res_ad.obs.assign(g2_UMIs = res_ad.layers["g2"].sum(axis=1))
    rs = read_meta(os.path.join(cache_dir,"rs.csv.gz"))
    pm_interaction = read_meta(os.path.join(cache_dir,"PM_interaction.csv.gz"))
    cdps = read_meta(os.path.join(cache_dir,"cdps.csv.gz"))
    res_ad.obs = pd.concat([res_ad.obs,pd.concat([rs, pm_interaction,annote[["group","partition","cell_type"]]],axis=1,join="inner")],axis=1)
    res_ad.obs = add_group_hour(res_ad.obs).loc[res_ad.obs_names]
    res_ad.obs = add_cell_type(res_ad.obs).loc[res_ad.obs.names]
    print("add uns...")
    res_ad.uns["cdps"] = cdps.loc[res_ad.obs_names]
    print("add g1 500k cs")
    g1cs = matr(os.path.join(cache_dir, "g1cs_500k.csv.gz"))
    res_ad.uns["g1cs_500k"] = g1cs.loc[res_ad.obs_names]
    print("add g2 500k cs")
    g2cs = matr(os.path.join(cache_dir, "g2cs_500k.csv.gz"))
    res_ad.uns["g2cs_500k"] = g2cs.loc[res_ad.obs_names]
    return res_ad