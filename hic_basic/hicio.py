import gzip
import json
import os
import pickle
import re
import tarfile
from collections import namedtuple
from io import StringIO
from pathlib import Path
from itertools import product, dropwhile

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from hires_utils.hires_io import parse_3dg
from scipy.io import mmread
from scipy.sparse import coo_matrix, csr_matrix




def divide_name(filename):
    #home-made os.path.splitext, for it can't handle "name.a.b.c" properly
    basename = os.path.basename(filename)
    parts = basename.split(".") #split return >= 1 length list
    if len(parts) == 1:
        return parts[0], ""
    else:
        return parts[0], "."+".".join(parts[1:])
def parse_seg(filename:str)->pd.DataFrame:
    """
    Read from dip-c's seg format.
    Input:
        filename: path to .seg file
    Output:
        pd.DataFrame
    """
    # compatible for zipped file
    if filename.endswith(".gz"):
        opener = gzip.open
    else:
        opener = open
    segs = []
    Seg = namedtuple("Seg", "chrom start end strand X1 Q X2".split())
    dtypes = [str, int, int, str, str, int, int]
    with opener(filename, "rt") as f:
        lines = dropwhile(lambda x:x.startswith("#"), f)
        seg_str_rows = (line.strip().split()[1:] for line in lines)
        seg_strs = (seg_str.split("!") for row in seg_str_rows for seg_str in row)
        segs = [Seg(*[dtype(x) for dtype, x in zip(dtypes, seg)]) for seg in seg_strs]
    res = pd.DataFrame(segs)
    return res
def parse_pairs(filename:str)->pd.DataFrame:
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
def parse_gtf(file,ID=True,name=True,mgi_id=True, **args):
    """
    Parsing gtf file. Read all in memory. Extract gene_id to df if set true.
    Input:
        file: path to .gtf file.
        ID: whether to parse gene_id from attribute string.
        name: whether to parse gene_name from attribute string.
        **args: arguments passed to pd.read_table
    Return:
        pd.DataFrame
    """
    header = ["seqname","source","feature","start","end","score","strand","frame","attributes"]
    dtypes = dict(zip(header, ["category","category","category","int","int","string","category","category","string"]))
    gtf = pd.read_table(file, comment ="#", header=None, names = header, dtype=dtypes,**args)
    # get gene_id from attributes
    if ID:
        ID = gtf["attributes"].str.extract(r"gene_id \"([\w,.]+)\";", expand=True)
        ID.columns = ["gene_id"]
        gtf = pd.concat([gtf, ID], axis=1, join="inner")
    if name:
        name = gtf["attributes"].str.extract(r"gene_name \"(\S+)\";", expand=True)
        name.columns = ["gene_name"]
        gtf = pd.concat([gtf, name], axis=1, join="inner")
    if mgi_id:
        mgi_id = gtf["attributes"].str.extract(r"mgi_id \"(MGI:\d+)\";", expand=True)
        mgi_id.columns = ["mgi_id"]
        gtf = pd.concat([gtf, mgi_id], axis=1, join="inner")
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
    

### --- read scRNA-seq file --- ###
def read_expr(path,sep=None)->pd.DataFrame:
    """
    Read expression matrix from file.
    Don't guarantee OV(obs*var) or VO(var*obs) output. It just read in the file.
    Input:
        path: path to matrix file
        sep: separator, if None will infer from file extension
    Output:
        expression matrix
    """
    mat = None

    # in case a dataframe is passed in
    if isinstance(path, pd.DataFrame):
        mat = path
        return mat

    # formats that don't require separator inference
    # like .parquet, .pkl, .pkl.gz
    if path.endswith("parquet"):
        mat = pd.read_parquet(path)
    elif path.endswith("pkl") or path.endswith("pkl.gz"):
        mat = pd.read_pickle(path)
    else:
        pass
    
    if mat is not None:
        return mat

    # format that requires separator inference
    if sep is None:
        if path.endswith("csv.gz") or path.endswith("csv"):
            mat = pd.read_csv(path,index_col=0)
        elif path.endswith("tsv.gz") or path.endswith("tsv"):
            mat = pd.read_table(path,index_col=0)
        else:
            print("Unknown file extension, will assume tab separated")
            mat = pd.read_table(path,index_col=0)
    else:
        mat = pd.read_table(path,sep=sep,index_col=0)
    mat.columns = mat.columns.astype("string")
    mat.index = mat.index.astype("string")
    return mat
def match_10x(file_list:list, retind=False)->dict:
    """
    Check if input file names contains a 10x result.
    Input:
        file_list: list of file path
        retind: return index of matched file instead of file path
    Output:
        {"matrix":matrix_file, "genes":gene_file, "barcodes":barcodes_file}
    """
    checker = {
        "matrix" : ["matrix"],
        "genes" : ["genes","features"],
        "barcodes" : ["barcodes"],
    }
    scores = {}
    result = {}
    file_list = [str(file) for file in file_list] # also accept iterable
    for filetype, keywords in checker.items():
        # score for each file
        scores = [
            sum([re.search(keyword, str(file).lower()) is not None for keyword in keywords])
            for file in file_list]
        # find file with highest score
        best_score = max(scores)
        if best_score > 0:
            if retind: # return index
                result[filetype] = scores.index(best_score)
            else: # return file path
                result[filetype] = file_list[scores.index(best_score)]
        else:
            result[filetype] = None
            print("Warning: No file found for %s" % filetype)
    return result
def combine_10x(annoted_file:dict)->pd.DataFrame:
    """
    Extract data from files and combine them into a expression matrix.
    Input:
        annoted_file: {"matrix":matrix_file, "genes":gene, "barcodes":barcodes}
    Output:
        a cell * gene matrix
    """
    pass
def read_10x(fp, gene_column=1, cell_column=0)->pd.DataFrame:
    """
    Read 10x-flavor matrix market file.
    Input:
        fp: file path
        gene_column: use this column from gene/feature matrix
        cell_column: use this column from cell matrix
    Output:
        a cell * gene matrix
    """
    fp = str(fp) # use end to decide input type
    if fp.endswith(".tar.gz"):
        with tarfile.open(fp, "r:gz") as tar:
            members = tar.getmembers()
            filenames = [member.name for member in members]
            annoted_members = {
                filetype : members[ind]
                for filetype, ind in match_10x(filenames, retind=True).items()
            }
            #print(annoted_members)
            for filetype, member in annoted_members.items():
                with tar.extractfile(member) as gf:
                    with gzip.open(gf,"rt") as f:
                        if filetype == "matrix":
                            data = mmread(f)
                        elif filetype == "genes":
                            genes = pd.read_table(
                                f,header=None)[gene_column].rename("var")
                        elif filetype == "barcodes":
                            barcodes = pd.read_table(
                                f,header=None)[cell_column].rename("obs")
                        else:
                            raise ValueError("Unknown file type: %s" % filetype)
        result = pd.DataFrame(data.todense(), index=genes, columns=barcodes)
    else:
        # treat as directory
        print("Treating as directory...")
        file_list = Path(fp).glob("*")
        for filetype, filename in match_10x(file_list, retind=False).items():
            if filename.endswith(".gz"):
                opener = gzip.open
            else:
                opener = open
            with opener(filename) as f:
                if filetype == "matrix":
                    # mmread accept filepath and call gzip in fact
                    # here just for consistency
                    data = mmread(f)
                elif filetype == "genes":
                    genes = pd.read_table(
                        f,header=None)[gene_column].rename("var")
                elif filetype == "barcodes":
                    barcodes = pd.read_table(
                        f,header=None)[cell_column].rename("obs")
                else:
                    raise ValueError("Unknown file type: %s" % filetype)
        result = pd.DataFrame(data.todense(), index=genes, columns=barcodes)
    return result
# read cooler file
def read_h5_ds(group):
    """
    Read .h5 bottom level group into dataframe.
    Input:
        group: h5py.Group; must have same length datasets.
    Output:
        pd.DataFrame
    """
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
# read schilcuster results
def parse_hicluster_res(embed, sample_table):
    """
    Read schicluster embedding result into dataframe.
    Input:
        embed: schicluster concat-cell output hdf5 file
        sample_table: input sample file of schicluster pipeline or list of sample names
    Examples:
        from hic_basic.io import parse_hicluster_res
        from hic_basic.scAB_embedding import do_umap
        hicluster_res = parse_hicluster_res(embed, sample_table)
        umap_res = do_umap(hicluster_res.values)
    """
    if isinstance(sample_table, str):
        sample_table = read_meta(sample_table)
        samples = sample_table.index
    else:
        samples = sample_table
    f = h5py.File(embed, "r")
    data = pd.DataFrame(f["data"][:], index=samples)
    f.close()
    data.columns = ["PC%d" % i for i in range(1, data.shape[1]+1)]
    return data
def schicluster2mat(filei):
    """
    Read in schicluster impute .hdf5 file to sparse matrix.
    Input:
        filei: schicluster imputed .hdf5 file
    Output:
        sparse matrix
    """
    with h5py.File(filei, "r") as f:
        g = f["Matrix"]
        A = csr_matrix((g["data"][()], g["indices"][()], g["indptr"][()]), g.attrs["shape"])
    return A
def schiclusterDir2mat(hiclusterdir):
    """
    Read in schicluster imputed .hdf5 files (from a dir) to sparse matrix.
    Input:
        hicluster: schicluster imputed .hdf5 file dir. e.g. f"/shareb/ychi/repo/sperm22_schicluster/imputed_matrix/100000/{chrom}/"
    Output:
        sparse matrix
    """
    Ap = sum(map(schicluster2mat, map(lambda x:os.path.join(hiclusterdir,x),os.listdir(hiclusterdir))))
    #Apc = mat_coarsen(Ap.todense(), coarsen)
    return Ap
# write scipy sparse matrix to 3-col file
def write_triplet(sparseM,filep,max_coo=False,zipping=True):
    """
    Write sparse matrix to triplet txt, assume symmetry.
    For scipy doesn't give a way to write triplet format directly.
    Input:
        sparseM: scipy sparse matrix
        filesp: output file path
        max_coo: if True, write largest possible coordinate at first row
    Output:
        indi, indj, data; separated by tab
    """
    sparseM.maxprint = sparseM.count_nonzero()
    if zipping:
        opener = gzip.open
    else:
        opener = open
    if max_coo == True:
        with opener(filep,"wt") as f:
            f.write("%d\n" % sparseM.shape[0])
        with opener(filep,"at") as f:
            for line in StringIO(str(sparseM)):
                coord, value = line.split("\t")
                c12 = " ".join(coord.strip().strip("()").split(", "))
                f.write(" ".join([c12, value]))
    else:
        with opener(filep,"wt") as f:
            for line in StringIO(str(sparseM)):
                coord, value = line.split("\t")
                c12 = " ".join(coord.strip().strip("()").split(", "))
                f.write(" ".join([c12, value]))
# --- private formats ---
def read_meta(fp):
    """
    Read general metadata, take care of sample_name.
    """
    if fp.endswith(".csv") or fp.endswith(".csv.gz"):
        df = pd.read_csv(fp,index_col=0,dtype={0:"string"})
    elif fp.endswith(".tsv") or fp.endswith(".tsv.gz"):
        df = pd.read_table(fp,index_col=0,dtype={0:"string"})
    else:
        raise ValueError("Unknown file extension")
    df.index.name = "sample_name"
    df.columns = df.columns.astype("string")
    return df
def read_umi_tools(path, sep=None)->pd.DataFrame:
    """
    Read umi_tools long-form output matrix.
    A wrapper for read_expr.
    By default umi_tools file store matrix in VO format, this function will transpose it to OV.
    See read_expr for more details.
    Input:
        path: path to matrix file
        sep: separator, if None will infer from file extension
    Output:
        obs * var matrix
    """
    return read_expr(path,sep).T
def matra(file):
    """
    Read umi-tools output to anndata object.
    """
    expr = read_umi_tools(file,"\t")
    expr.columns = expr.columns.astype("str")
    expr.columns.name = "sample_name"
    expr.index = expr.index.astype("str")
    expr.index.name = "gene"
    adata = ad.AnnData(expr.T)
    return adata
def get_chrom_contact_counts(dump_dir):
    """
    Get per-chromosome-per-phase contacts counts of samples.
    Input:
        dump_dir: rd/dump generated by hires_utils
    Output:
        pd.DataFrame
    """
    res = pd.concat((pd.read_pickle(file) for file in [os.path.join(dump_dir, i) for i in os.listdir(dump_dir)]),axis=1)
    res.columns = [i.split("_")[0] for i in os.listdir(dump_dir)]
    return res
def align_to_df(primary_views):
    """
    Extract general dataframe from primary_view function output.
    Input:
        primary_views: complex data structure list of dict
    Output:
        dataframe
    """
    mapper = {
        "head-tail":0,
        "dorsal-ventral":1,
        "left-right":2,
    }
    dfs = [[],[],[]]
    samples = []  
    for sample in primary_views:
        samples.append(sample)
        for aname, figures in zip(
            primary_views[sample]["name_of_vectors"],
            primary_views[sample]["primary_figures"]
        ):
            for figure in figures: # 2 figure for each axis
                figure[figure == np.inf] = np.nan
                dfs[mapper[aname]].append(figure.reshape(1,-1))
    index = list(product(samples, ["A","B"]))
    #return dfs
    dfs = [
        pd.DataFrame(
            np.concatenate(df,axis=0),
            index=pd.MultiIndex.from_product([samples, ["A", "B"]])
        )
        for df in dfs
    ]
    return dfs
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
def load_pickle(filep):
    """
    Load from pickle file shortcut.
    """
    with open(filep,"rb") as f:
        return pickle.load(f)
def dump_pickle(obj, filep):
    """
    Dump to pickle shortcut.
    """
    with open(filep,"wb") as f:
        pickle.dump(obj,f)
def get_ref_dir():
    """
    Return a path to src's ref
    """
    ref_dir = os.path.join(os.path.dirname(__file__), "ref")
    return os.path.join(ref_dir, "")