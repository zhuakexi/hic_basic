import json
import os
import gzip

import pandas as pd
import numpy as np
import h5py
from scipy.sparse import coo_matrix
import anndata as ad
from io import StringIO
def divide_name(filename):
    #home-made os.path.splitext, for it can't handle "name.a.b.c" properly
    basename = os.path.basename(filename)
    parts = basename.split(".") #split return >= 1 length list
    if len(parts) == 1:
        return parts[0], ""
    else:
        return parts[0], "."+".".join(parts[1:])
def parse_pairs(filename:str)->"Cell":
    '''
    read from 4DN's standard .pairs format
    compatible with all hickit originated pairs-like format 
    '''
    # read comments
    with gzip.open(filename,"rt") as f:
        comments = []
        chromosomes = []
        lengths = []
        for line in f.readlines():
            if line[0] != "#":
                break
            if line.startswith("#chromosome") or line.startswith("#chromsize"):
                chrom, length = line.split(":")[1].strip().split()
                chromosomes.append(chrom)
                lengths.append(int(length))
            if line.startswith("#columns:"):
                columns = line.split(":")[1].strip().split()
            ## comment lines are stored in dataframe.attrs["comment"]
            comments.append(line)
    dtype_array = {"readID":"category",
            "chr1":pd.CategoricalDtype(categories=chromosomes),
            "pos1":"int",
            "chr2":pd.CategoricalDtype(categories=chromosomes),
            "pos2":"int",
            "strand1":pd.CategoricalDtype(categories=["+","-"]),
            "strand2":pd.CategoricalDtype(categories=["+","-"]),
            "phase0":pd.CategoricalDtype(categories=["1","0","."]),
            "phase1":pd.CategoricalDtype(categories=["1","0","."]),
            "phase_prob00":"float",
            "phase_prob01":"float",
            "phase_prob10":"float",
            "phase_prob11":"float"}
    dtypes = {key:value for key, value in dtype_array.items() if key in columns}
    #read table format data
    pairs = pd.read_table(
        filename, 
        header=None, 
        comment="#",
        dtype=dtypes,
        names=columns
        )
    pairs.attrs["comments"] = comments
    pairs.attrs["name"], _ = divide_name(filename) # infer real sample name
    pairs.attrs["chromosomes"] = chromosomes
    pairs.attrs["lengths"] = lengths
    #assign column names
    #sys.stderr.write("pairs_parser: %s parsed \n" % filename)
    return pairs
def parse_pairs_like(filename:str)->pd.DataFrame:
    '''
    read from 4DN's standard .pairs format
    compatible with all hickit originated pairs-like format 
    '''
    #comment lines are stored in dataframe.attrs["comment"]
    name_array = "readID chr1 pos1 chr2 pos2 strand1 strand2 phase0 phase1 phase_prob00 phase_prob01 phase_prob10 phase_prob11".split()
    dtype_array = {"readID":"category",
            "chr1":"string",
            "pos1":"int",
            "chr2":"string",
            "pos2":"int",
            "strand1":"string",
            "strand2":"string",
            "phase0":"string",
            "phase1":"string",
            "phase_prob00":"float",
            "phase_prob01":"float",
            "phase_prob10":"float",
            "phase_prob11":"float"}
    #read comment line
    with gzip.open(filename,"rt") as f:
        comments = []
        for line in f.readlines():
            if line[0] != "#":
                break
            comments.append(line)
    #infer number of columns
    line_length = len(line.strip().split("\t"))
    #pick used eles from builtin arrays
    columns = name_array[0:line_length]
    dtypes = {key:value for key, value in dtype_array.items() if key in columns}
    #read table format data
    pairs = pd.read_table(
        filename, 
        header=None, 
        comment="#",
        dtype=dtypes,
        names=columns
        )
    pairs.attrs["comments"] = comments
    pairs.attrs["name"], _ = divide_name(filename) # infer real sample name
    #assign column names
    #sys.stderr.write("pairs_parser: %s parsed \n" % filename)
    return pairs
def parse_gtf(file,ID=True):
    """
    Parsing gtf file. Read all in memory. Extract gene_id to df if set true.
    """
    header = ["seqname","source","feature","start","end","score","strand","frame","attributes"]
    dtypes = dict(zip(header, ["category","category","category","int","int","string","category","category","category","string"]))
    gtf = pd.read_table(file, comment ="#", header=None, names = header, dtype=dtypes)
    # get gene_id from attributes
    if ID:
        ID = gtf["attributes"].str.extract(r"gene_id \"(\w+)\";", expand=True)
        ID.columns = ["gene_id"]
        gtf = pd.concat([gtf, ID],axis=1, join="inner")
    return gtf
def parse_gff(file, ID=False, Name=False):
    """
    Parsing gff file. Read all in memory. Extract ID(gene), and Name if set true.
    """
    header = ["seqid","annotation_source","feature_type","start","end","score","strand","phase","attributes"]
    dtypes = dict(zip(header,["category","category","category","int","int","string","category","category","string"]))
    gff = pd.read_table(file,
                        comment="#",
                        header=None, 
                        names = header,
                        dtype = dtypes
                       )
    additions = []
    if ID:
        ID = gff["attributes"].str.extract(r"ID=gene:(\w+);",expand=True)
        ID.columns = ["ID"]
        additions.append(ID)
    if Name:
        Name = gff["attributes"].str.extract(r"Name=(\w+);",expand=True)
        Name.columns = ["Name"]
        additions.append(Name)
    gff = pd.concat([gff] + additions, axis=1, join="inner")
    return gff
# read cooler file
def read_h5_ds(group):
    # read bottom level group into dataframe
    data = {}
    for key in group.keys():
        if group[key].dtype.type == np.string_:
            data[key] = group[key][:].astype("U")
        else:
            data[key] = group[key][:]
    return pd.DataFrame(data)
def load_cool(cool, root="/"):
    """
    Load simple cooler file as sparsematrix
    
    Parameters
    ----------
    cool : str
        Path to the input .cool file.

    Returns
    -------
    mat : scipy coo_matrix
        Hi-C contact map in COO format.
    frags : pandas DataFrame
        Table of bins matching the matrix.
    chroms : pandas DataFrame
        Table of chromosome informations.
    """
    with h5py.File(cool,"r") as f:
        frags = read_h5_ds(f[root]["bins"])
        n_frags = frags.groupby("chrom", sort=False).count()["start"]
        n_frags.name = "n_frags"
        chroms = read_h5_ds(f[root]["chroms"])
        mat = read_h5_ds(f[root]["pixels"])
    frags["id"] = frags.groupby("chrom", sort=False).cumcount() + 1
    chroms["cumul_length"] = (
        chroms.length.shift(1).fillna(0).cumsum().astype(int)
    )
    chroms = pd.concat([chroms,n_frags],axis=1)    
    n = int(max(np.amax(mat.bin1_id), np.amax(mat.bin2_id))) + 1
    shape = (n, n)
    mat = coo_matrix((mat["count"], (mat.bin1_id, mat.bin2_id)), shape=shape)

    return mat, frags, chroms
# write scipy sparse matrix to 3-col file
def write_triplet(sparseM,filep,max_coo=True):
    """
    # write sparse matrix to triplet txt, assume symmetry
    # write largest possible coordinate afront by default
    """
    sparseM.maxprint = sparseM.count_nonzero()
    if max_coo == True:
        with open(filep,"wt") as f:
            f.write("%d\n" % sparseM.shape[0])
        with open(filep,"at") as f:
            for line in StringIO(str(sparseM)):
                coord, value = line.split("\t")
                c12 = " ".join(coord.strip().strip("()").split(", "))
                f.write(" ".join([c12, value]))
    else:
        with open(filep,"wt") as f:
            for line in StringIO(str(sparseM)):
                coord, value = line.split("\t")
                c12 = " ".join(coord.strip().strip("()").split(", "))
                f.write(" ".join([c12, value]))
# --- private formats ---
def read_meta(fp):
    """
    Read general metadata, take care of sample_name.
    """
    df = pd.read_csv(fp, dtype={"sample_name":"string"},index_col=0)
    df.index.name = "sample_name"
    return df
def matr(path,sep=","):
    """
    Read umi_tools long-form output matrix.
    """
    mat = pd.read_csv(path,sep=sep,index_col=0)
    mat.columns = mat.columns.astype("string")
    mat.index = mat.index.astype("string")
    return mat
def matra(file):
    """
    Read umi-tools output to anndata object.
    """
    expr = matr(file,"\t")
    expr.columns = expr.columns.astype("str")
    expr.columns.name = "sample_name"
    expr.index = expr.index.astype("str")
    expr.index.name = "gene"
    adata = ad.AnnData(expr.T)
    return adata
# --- misc ---
def dump_json(obj, filep):
    """
    Dump to json shortcut.
    """
    with open(filep,"wt") as f:
        json.dump(obj,f)
def load_json(filep):
    """
    Load from json file shortcut.
    """
    with open(filep,"rt") as f:
        return json.load(f)