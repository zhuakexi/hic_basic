from functools import partial
import anndata as ad
import pandas as pd
import numpy as np
from ..repli_score import _repli_score
from scipy import sparse
import os
import loompy
import scvelo as scv

from ..io import matr, read_meta, get_chrom_contact_counts
from .paracalc import gen_repli_score, gen_cdps, gen_PM_interactions
from .exp_record import add_cell_type, add_group_hour
from .utils import two_sets
from ..phasing import get_chrom_hap_score

def _merge_expr(fps, mapper, samplelist=None, outfile=None, qc=None):
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
    if outfile is not None:
        all_mat.to_csv(outfile)
    return all_mat
def gen_expr(qc, mattail="count_matrix/counts.gene.tsv.gz"):
    """
    Generate single expression matrix that contains all (has valid task directory) qc passed samples from qc file.
    qc file must contain "task_dirp" col and use sample ID as index.
    fix sample ID with "_" for umi_tools limitation. 
    Input:
        qc: filtered dataframe
        mattail: RNA mat position according to task_dirp
    Output:
        merged expression matrix
    TODO:
        Adding `outfp`; if not None, output file rather than real dataframe
    """
    if qc.index.str.contains("_").sum() > 0:
        print("Warning: sample ID contains underscore")
        return _merge_expr(
            [os.path.join(i,mattail) for i in qc["task_dirp"].unique()],
            mapper = dict(zip(qc.index[qc.index.str.contains("_")].str.split("_",expand=True).get_level_values(1),
                            qc.index[qc.index.str.contains("_")])),
            samplelist = qc.index.values
        )
    else:
        return _merge_expr(
            [os.path.join(i,mattail) for i in qc["task_dirp"].unique()],
            mapper = dict(zip(qc.index.values,qc.index.values)),
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
    target_df.index = target_df.index.astype("string")
    target_df.columns = target_df.columns.astype("string")
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
# --- functions to generate all mid-files used for generating adata object ---
# sub function for generating g1 g2 cs...
## inferring the file name from the qc file
### Input: [pos] qc, outfile
### Output: ([],{})
def _infer_expr_inputs(qc, outfile):
    # infer expression matrix file paths and other arguments
    mattail = "count_matrix/counts.gene.tsv.gz"
    if qc.index.str.contains("_").sum() > 0:
        print("Warning: sample ID contains underscore")
        return [[os.path.join(i, mattail) for i in qc["task_dirp"].unique()]], {
            "mapper" :  dict(zip(qc.index[qc.index.str.contains("_")].str.split("_",expand=True).get_level_values(1),
                            qc.index[qc.index.str.contains("_")])),
            "samplelist" : qc.index.values,
            "outfile" : outfile
            }
    else:
        return [[os.path.join(i, mattail) for i in qc["task_dirp"].unique()]], {
            "mapper" : dict(zip(qc.index.values,qc.index.values)),
            "samplelist" : qc.index.values,
            "outfile" : outfile
            }
## generate and caching
### Input: exposed to output, must have qc and outfile key arguments
### Output: single object
def _gen_velo(velo_files, qc=None, outfile=None):
    print("Gen velo cache...")
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
    return velo
def _gen_cdps(qc=None, outfile=None):
    print("Gen cdps...")
    cdps = gen_cdps(qc)
    cdps.to_csv(outfile)
    return cdps
def _gen_repliscore(outdir, qc = None, outfile=None):
    print("Gen repli score...")
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    outfile = os.path.join(outdir, "rs.csv.gz")
    if not os.path.isfile(outfile):
        rs = gen_repli_score(qc, outdir)
        rs.to_csv(outfile)
    return rs
def _generate_PM_interactions(qc=None, outfile=None):
    print("Gen PM interactions...")
    PM_interaction = gen_PM_interactions(qc)
    PM_interaction.index.name = "sample_name"
    PM_interaction.to_csv(outfile)
    return PM_interaction
def _gen_g1_cs(cs_dir, qc=None, outfile=None):
    print("Gen cs...")
    cs_fileis = gen_fileis(qc, cs_dir,  "{}.500000.cs.g1.txt")
    cs_500k = gen_cs(cs_fileis)
    cs_500k.to_csv(outfile)
    return cs_500k
def _gen_g2_cs(cs_dir, qc=None, outfile=None):
    print("Gen cs...")
    cs_fileis = gen_fileis(qc, cs_dir,  "{}.500000.cs.g2.txt")
    cs_500k = gen_cs(cs_fileis)
    cs_500k.to_csv(outfile)
    return cs_500k
def _gen_annote(cols, qc=None, outfile=None):
    print("Gen annote...")
    annote = qc[cols]
    annote.to_csv(outfile)
    return annote
def _gen_g(fps, qc = None, outfile=None):
    print("Gen phased count matrix...")
    if qc.index.str.contains("_").sum() > 0:
        print("Warning: sample ID contains underscore")
        return _merge_expr(fps, 
            mapper =  dict(zip(qc.index[qc.index.str.contains("_")].str.split("_",expand=True).get_level_values(1),
                            qc.index[qc.index.str.contains("_")])),
            samplelist =  qc.index.values,
            outfile = outfile
            )
    else:
        return _merge_expr(fps,
            mapper = dict(zip(qc.index.values,qc.index.values)),
            samplelist = qc.index.values,
            outfile = outfile)
def _gen_chrom_hap_score(dump_dirs, qc = None, outfile=None):
    print("Gen chrom hap socre...")
    chrom_hap_scores = []
    for dir in dump_dirs:
        if not os.path.isdir(dir):
            raise ValueError("{} is not a directory".format(dir))
        chrom_hap_scores.append(get_chrom_hap_score(dir))
    chrom_hap_score = pd.DataFrame(pd.concat(chrom_hap_scores, axis=0), columns = ["chrom_hap_score"])
    chrom_hap_score.to_csv(outfile)
    return chrom_hap_score
def _gen_chrom_contacts(dump_dirs, qc = None, outfile=None):
    print("Gen chrom contacts...")
    chrom_contacts = []
    for dir in dump_dirs:
        if not os.path.isdir(dir):
            raise ValueError("{} is not a directory".format(dir))
        chrom_contacts.append(get_chrom_contact_counts(dir))
    chrom_contacts = pd.concat(chrom_contacts, axis=1)
    chrom_contacts.to_csv(outfile)
    return chrom_contacts
## generate dataframes for each layer
## tidy
def _tidy_velo(velo_ad):
    velo_ad.var_names_make_unique()
    velo_ad = velo_ad.T
    return velo_ad
def _tidy_expr(expr):
    expr.columns = expr.columns.astype("string")
    return expr
## adding annotations to adata object
def _add_g1_UMIs(adata=None, qc=None):
    return adata.obs.assign(g1_UMIs = adata.layers["g1"].sum(axis=1))
def _add_g2_UMIs(adata=None, qc=None):
    return adata.obs.assign(g2_UMIs = adata.layers["g2"].sum(axis=1))
def _add_concat(adata=None, df=None, qc = None):
    if not set(adata.obs.index).issubset(df.index):
        raise ValueError("[adding annotation]: adata.obs.index is not a subset of df.index")
    return pd.concat([adata.obs, df],axis=1, join="inner")
def create_adata_layers(ws, funcs):
    """
    Create adata object form the layers
    Input:
        ws: dict of layer name and dataframe
        funcs: functions to access layers
    Output:
        adata object
    """
    print("Creating adata object from layer dataframes...")
    # --- tidy up generated dataframes ---:
    for i in ws:
        if i in funcs["tidy"]:
            ws[i] = funcs["tidy"][i](ws[i])
    # --- get consensus gene and sample names ---
    samples = list(set.intersection(*map(set, [funcs["gi"][i][0](ws[i]) for i in ws])))
    genes = list(set.union(*map(set, [funcs["gi"][i][1](ws[i]) for i in ws])))
    # --- generate seperate dfs for layers ---
    # now in new dfs namespace
    dfs = {}
    for i in ws:
        if i in ("expr", "velo", "g1", "g2"):
            if i in funcs["gdf"]:
                dfs.update(funcs["gdf"][i](ws[i]))
            else:
                dfs[i] = ws[i]
    # --- future .X template ---
    pixels = pd.DataFrame.sparse.from_spmatrix(sparse.csr_matrix((len(genes),len(samples)),dtype=np.int64),index = genes, columns = samples)
    pixels.index = pixels.index.astype("string")
    pixels.columns = pixels.columns.astype("string")
    # --- expand dfs ---
    print("Expanding dfs...")
    dfs = {i : expand_df(dfs[i][samples], pixels) for i in dfs}
    # --- create new annData object ---
    print("creating new AnnData object")
    #print(dfs["expr"].loc[genes, samples].T.index)
    #print(pd.Index(samples, dtype="string"))
    new_ad = ad.AnnData(
        X = dfs["expr"].loc[genes, samples].T,
        obs = pd.DataFrame(index=pd.Index(samples, dtype="string")),
        var = pd.DataFrame(index=pd.Index(genes, dtype="string")),
        layers = {i : dfs[i].loc[genes, samples].T for i in dfs if i not in ["expr"]}
    )
    return new_ad
def add_obsm(adata, data, obsm_key):
    """
    Adding obsm to anndata.
    Input:
        data: low-dim representation of samples. samples x feature
        obsm_key: name of the representation
    """
    new_data = pd.DataFrame(index=adata.obs_names)

    two_sets(adata.obs_names, data.index, True)
    adata.obsm[obsm_key] = pd.concat(
        [   new_data,
            data.loc[data.index.intersection(adata.obs_names)]
            ],
        axis=1,
        join="outer"
        )
    return adata
def gen_adata(qc, cache_dir, rewrite=[], **args):
    ca = lambda x: os.path.join(cache_dir, x)
    funcs = {
        # not essentail
        # must output (list of positional arguments, dict of keyword arguments)
        "infer_inputs" :
            {
                "expr" : _infer_expr_inputs
            },
        # not essentail
        # must output single object
        "cache_files" :
            {
                "expr" : ca("expr.csv.gz"),
                "velo" : ca("velocyto.integrate.loom"),
                "g1" : ca("g1.csv.gz"),
                "g2" : ca("g2.csv.gz"),
                "cdps" : ca("cdps.csv.gz"),
                "rs" : ca("rs.csv.gz"),
                "pm" : ca("pm.csv.gz"),
                "g1cs" : ca("g1cs_500k.csv.gz"),
                "g2cs" : ca("g2cs_500k.csv.gz"),
                "annote" : ca("annote.csv.gz"),
                "chrom_hap_score" : ca("chrom_hap_score.csv.gz"),
                "chrom_contacts" : ca("chrom_contacts.csv.gz"),
            },
        # not esential
        "gen" :
            {
                "expr" : partial(_merge_expr, samplelist = qc.index.values),
                "velo" : _gen_velo,
                "g1" : _gen_g,
                "g2" : _gen_g,
                # no input needed
                "cdps" : _gen_cdps,
                # no input needed
                "rs" : _gen_repliscore,
                # no input needed
                "pm" : _generate_PM_interactions,
                # 
                "g1cs" : _gen_g1_cs,
                "g2cs" : _gen_g2_cs,
                "annote" : _gen_annote,
                "chrom_hap_score" : _gen_chrom_hap_score,
                "chrom_contacts" : _gen_chrom_contacts,
            },
        # essentail for all
        "read_files" :
            {
                "expr" : matr,
                "velo" : ad.read_loom,
                "g1" : matr,
                "g2" : matr,
                "cdps" : read_meta,
                "rs" : read_meta,
                "pm" : read_meta,
                "g1cs" : matr,
                "g2cs" : matr,
                "annote" : read_meta,
                "chrom_hap_score" : read_meta,
                "chrom_contacts" : partial(pd.read_csv, index_col=[0,1]),
            },
        "tidy" :
            {
                "expr" : _tidy_expr,
                "velo" : _tidy_velo,
            },
        # essential for all
        # [get_sample_names, get_gene_names]
        "gi" :
            {
                "expr" : [lambda x : x.columns, lambda x : x.index],
                "g1" : [lambda x : x.columns, lambda x : x.index],
                "g2" : [lambda x : x.columns, lambda x : x.index],
                # velo is transposed during tidy
                "velo" : [lambda x : x.var.index, lambda x : x.obs.index]
            },
        # essential for all
        "gdf" :
            {
                "velo" : lambda x : {
                    "matrix" : pd.DataFrame.sparse.from_spmatrix(x.layers["matrix"], index = x.obs.index, columns = x.var.index),
                    "spliced" : pd.DataFrame.sparse.from_spmatrix(x.layers["spliced"], index = x.obs.index, columns = x.var.index),
                    "unspliced" : pd.DataFrame.sparse.from_spmatrix(x.layers["unspliced"], index = x.obs.index, columns = x.var.index),
                    "ambiguous" : pd.DataFrame.sparse.from_spmatrix(x.layers["ambiguous"], index = x.obs.index, columns = x.var.index)
                }
            },
        "add_annote" :
            {
                "rs" : _add_concat,
                "g1_UMIs" : _add_g1_UMIs,
                "g2_UMIs" : _add_g2_UMIs,
                "pm" : _add_concat,
                "annote" : _add_concat,
                "chrom_hap_score" : _add_concat
            }
    }
    ws = args
    # --- infer inputs ---:
    for i in ws:
        if ws[i] is None:
            if i in funcs["infer_inputs"]:
                print("Inferring inputs for {}...".format(i))
                ws[i] = funcs["infer_inputs"][i](qc, funcs["cache_files"][i])
            else:
                print("No input for {}".format(i))
                continue
    print("Generating mid-files...")
    # --- generate and cache ---:
    for i in ws:
        if i not in funcs["cache_files"] and i not in funcs["gen"]:
            print("Direct attr {}".format(i))
            continue
        if os.path.isfile(funcs["cache_files"][i]) & (i not in rewrite):
            # already generated, read from cache
            print("Reading {} from cache...".format(i))
            ws[i] = funcs["read_files"][i](funcs["cache_files"][i])
        else:
            # generate and cache
            print("Generating {}...".format(i))
            if isinstance(ws[i], tuple):
                if (len(ws[i]) == 2):
                    if isinstance(ws[i][0], list) & isinstance(ws[i][1], dict):
                        # input contains additional arguments
                        ws[i] = funcs["gen"][i](*ws[i][0], **ws[i][1])
            else:
                # input is only list of file paths
                if ws[i] is None:
                    ws[i] = funcs["gen"][i](qc = qc, outfile=funcs["cache_files"][i])
                else:
                    ws[i] = funcs["gen"][i](ws[i], qc = qc, outfile=funcs["cache_files"][i])
    # --- seperate by layer, obsm, ... ---:
    layers = {i : ws[i] for i in ws if i in ("expr", "velo", "g1", "g2")}
    obs = {i : ws[i] for i in ws if i in ("rs", "pm", "g1_UMIs", "g2_UMIs", "annote", "chrom_hap_score")}
    uns = {i : ws[i] for i in ws if i in ("cdps", "g1cs", "g2cs")}
    # --- create adata ---:
    print("Creating adata...")
    adata = create_adata_layers(layers, funcs)
    # --- adding per-cell annotations ---:
    print("Adding per-cell annotations...")
    for i in obs:
        if i in funcs["cache_files"]:
            adata.obs = funcs["add_annote"][i](
                adata = adata, 
                df = funcs["read_files"][i](funcs["cache_files"][i]), 
                qc=qc
            )
        else:
            adata.obs = funcs["add_annote"][i](
                adata = adata,
                qc = qc
            )
    # --- adding uns ---:
    print("Adding uns...")
    for i in uns:
        adata.uns[i] = uns[i].loc[adata.obs_names]
    return adata
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
    res_ad.obs = add_cell_type(res_ad.obs).loc[res_ad.obs_names]
    print("add uns...")
    res_ad.uns["cdps"] = cdps.loc[res_ad.obs_names]
    print("add g1 500k cs")
    g1cs = matr(os.path.join(cache_dir, "g1cs_500k.csv.gz"))
    res_ad.uns["g1cs_500k"] = g1cs.loc[res_ad.obs_names]
    print("add g2 500k cs")
    g2cs = matr(os.path.join(cache_dir, "g2cs_500k.csv.gz"))
    res_ad.uns["g2cs_500k"] = g2cs.loc[res_ad.obs_names]
    return res_ad