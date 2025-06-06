from functools import partial
from time import time

import anndata as ad
import pandas as pd
import numpy as np
from ..repli_score import _repli_score
from scipy import sparse
import os
import loompy
import scvelo as scv

from ..hicio import read_umi_tools, read_expr, read_meta, get_chrom_contact_counts
from .paracalc import gen_repli_score, gen_cdps, gen_PM_interactions
from .exp_record import add_cell_type, add_group_hour
from ..utils import two_sets, gen_fileis
from ..phasing import get_chrom_hap_score
from ..pipeline.rule import bubble_flow_touched
from ..pp import _merge_expr
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
def expand_df(target_df, ref_df, count_matrix=False, fill_na_value=0):
    """
    expand target_df to shape(with same index, columns) of ref_df, fill with 0
    Input:
        target_df: df to be transformed; sparse df
        ref_df: axis reference, must contain all rows and cols of target_df; sparse df
        count_matrix: if True, use sparse matrix, else use dense matrix
    """
    if count_matrix:
        print("expand_df: `count_matrix` is True, use sparse matrix. Don't do this if your data is not integer")
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
    target_df.fillna(fill_na_value, inplace=True)
    if count_matrix:
        target_df = target_df.astype(pd.SparseDtype(int, fill_value=0))
    target_df.index = target_df.index.astype("string")
    target_df.columns = target_df.columns.astype("string")
    return target_df
# --- functions to generate all mid-files used for generating adata object ---
# sub function for generating g1 g2 cs...
## inferring the file name from the qc file
### Input: [pos] qc, outfile
### Output: ([],{})
def _infer_expr_inputs(qc, outfile):
    # infer expression matrix file paths and other arguments
    if "ref" not in qc.columns:
        print("[_infer_expr_inputs] Warning: ref column not found in qc file, using blank string for ref name")
        combs = [(i, "") for i in qc["task_dirp"].unique()]
    else:
        combs = qc[["task_dirp", "ref"]].value_counts().index.tolist()
    count_matrix_pat = bubble_flow_touched()["count_matrix"][0] # this function always return a list, so [0] is needed
    if qc.index.str.contains("_").sum() > 0:
        print("Warning: sample ID contains underscore")
        return [[count_matrix_pat.format(task_dirp=task_dirp, ref=ref) for task_dirp, ref in combs]], {
            "mapper" :  dict(zip(
                            [i.split("_")[-1] for i in qc.index if "_" in i],
                            [i for i in qc.index if "_" in i]
                            )),
            "samplelist" : qc.index.tolist(),
            "outfile" : outfile
            }
    else:
        return [[count_matrix_pat.format(task_dirp=task_dirp, ref=ref) for task_dirp, ref in combs]], {
            "mapper" : dict(zip(qc.index.tolist(),qc.index.tolist())),
            "samplelist" : qc.index.tolist(),
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
    cdps = gen_cdps(qc, range_dtype="string")
    cdps.to_csv(outfile)
    return cdps
def _gen_repliscore(outdir, qc = None, outfile=None):
    print("Gen repli score...")
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    rs = gen_repli_score(qc, outdir) # gen_repli_score has its own caching
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
            mapper =  dict(zip(
                            [i.split("_")[-1] for i in qc.index if "_" in i],
                            [i for i in qc.index if "_" in i]
                            )),
            samplelist =  qc.index.tolist(),
            outfile = outfile
            )
    else:
        return _merge_expr(fps,
            mapper = dict(zip(qc.index.tolist(),qc.index.tolist())),
            samplelist = qc.index.tolist(),
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
def _add_group_hour(adata=None, qc=None):
    data = add_group_hour(qc, "collect_hour")
    return adata.obs.assign(collect_hour = data["collect_hour"])
def _add_cell_type(adata=None, qc=None):
    data = add_cell_type(qc)
    return adata.obs.assign(cell_type = data["cell_type"])
def create_adata_layers(ws, funcs, debug=False):
    """
    Create adata object form the layers
    Input:
        ws: dict of layer name and dataframe
        funcs: functions to access layers
    Output:
        adata object
    TODO: figure out what cause the "Index of obs must match index of X." error
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
    pixels = pd.DataFrame.sparse.from_spmatrix(sparse.csr_matrix((len(samples),len(genes)),dtype=np.int64),index = samples, columns = genes)
    pixels.index = pixels.index.astype("string")
    pixels.columns = pixels.columns.astype("string")
    # --- expand dfs ---
    print("Expanding dfs...")
    dfs = {layer : expand_df(df.loc[samples], pixels) for layer, df in dfs.items()}
    # --- create new annData object ---
    print("creating new AnnData object")
    if debug:
        X = dfs["expr"].loc[samples, genes]
        obs = pd.DataFrame(index=pd.Index(samples, dtype="string"))
        var = pd.DataFrame(index=pd.Index(genes, dtype="string"))
        print(X.index)
        print(obs.index)
        print("X.index.equals(obs.index):", X.index.equals(obs.index))
        print("X.columns.equals(var.index):", X.columns.equals(var.index))
        # don't pass obs and var to AnnData, will rise "Index of obs must match index of X."
        # adata = ad.AnnData(
        #     X = X,
        #     obs = obs,
        #     var = var
        # )
    # new_ad = ad.AnnData(
    #     X = dfs["expr"].loc[genes, samples].T,
    #     obs = pd.DataFrame(index=pd.Index(samples, dtype="string")),
    #     var = pd.DataFrame(index=pd.Index(genes, dtype="string")),
    #     layers = {i : dfs[i].loc[genes, samples].T for i in dfs if i not in ["expr"]}
    # )
    new_ad = ad.AnnData(
        X = dfs["expr"].loc[samples, genes],
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
def gen_adata(qc, cache_dir, rewrite=[], debug=False, **args):
    """
    Generate anndata object from the QC data.
    Input:
        qc: QC dataframe
        cache_dir: directory to store the generated adata object
        rewrite: list of dataframe name to rewrite

        --- layers ---
        expr: expression matrix
            Generate single expression matrix that contains all (has valid task directory) qc passed samples from qc file.
            qc file must contain "task_dirp" col and use sample ID as index.
            fix sample ID with "_" for umi_tools limitation. 
        g1: phase 0 expression matrix
        g2: phase 1 expression matrix
        velo: velocity matrix
        --- uns ---
        cdps: contact decay profile
        g1cs: phase 0 compartment strength
        g2cs: phase 1 compartment strength
        --- obs ---
        g1_UMIs: phase 0 UMIs
        g2_UMIs: phase 1 UMIs
        pm: phase 0 1 interaction(contacts)
        rs: repli_score
        collect_hour: collection hour parsed from group name
        cell_type: cell type parsed from group name
    TODO: for samples without results in some uns, add proper nan values
    """
    # --- set up ---
    qc.index = qc.index.astype("string")
    qc.columns = qc.columns.astype("string")

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
                "expr" : partial(_merge_expr, samplelist = qc.index.tolist()),
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
                "expr" : read_umi_tools,
                "velo" : ad.read_loom,
                "g1" : read_umi_tools,
                "g2" : read_umi_tools,
                "cdps" : read_meta,
                "rs" : read_meta,
                "pm" : read_meta,
                "g1cs" : read_umi_tools,
                "g2cs" : read_umi_tools,
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
                "expr" : [lambda x : x.index, lambda x : x.columns],
                "g1" : [lambda x : x.index, lambda x : x.columns],
                "g2" : [lambda x : x.index, lambda x : x.columns],
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
                "chrom_hap_score" : _add_concat,
                "collect_hour" : _add_group_hour,
                "cell_type" : _add_cell_type,
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
    obs = {i : ws[i] for i in ws if i in ("rs", "pm", "g1_UMIs", "g2_UMIs", "annote", "chrom_hap_score",
                                          "chrom_contacts", "collect_hour", "cell_type", "g1_UMIs", "g2_UMIs",
                                          "cell_type")}
    uns = {i : ws[i] for i in ws if i in ("cdps", "g1cs", "g2cs")}
    # --- create adata ---:
    print("Creating adata...")
    adata = create_adata_layers(layers, funcs, debug=debug)
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