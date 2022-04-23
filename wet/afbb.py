# things to do after bubbleflow
# add_pairs -> add_pairs_num -> add_info ... -> pick_useful
import os
import json
from concurrent import futures
import gzip
from subprocess import check_output, CalledProcessError
import pandas as pd
import re
# add pairs num
def zcount(filename:str,target:str)->int:
    # count zipped file line number
    try:
        output_bytes = check_output([
            "zgrep","-c",target,
            filename])
    except CalledProcessError:
        return 0
    return int(output_bytes.decode("utf-8").strip())
def count_pairs(filename:str)->int:
    # count number of contacts in 4DN pairs file
    if not isinstance(filename,str):
        return None
    if not os.path.isfile(filename):
        return None
    comments = 0
    with gzip.open(filename,"rt") as f:
        for line in f:
            if line[0] == "#":
                comments += 1
    all_lines = zcount(filename, "$")
    return all_lines - comments
def add_pairs_num(annote,threads,pairs_types=["pairs_c1","pairs_c12","pairs_c123","dip"]):
    """
    count contact number of pairs file in parallel
    Input:
        annote: must have all pairs_type cols
    Output:
        annote with extra cols, pairs_type + "_num" by default
    """
    filesp = annote
    ares = {}
    for pairs_type in pairs_types:
        with futures.ProcessPoolExecutor(threads) as pool:
            res = pool.map(
                count_pairs,
                filesp[pairs_type]
                          )
        ares[pairs_type+"_num"] = list(res)
    #return ares
    return filesp.assign(**ares)
# add_pairs path
def real_file(i):
    if not isinstance(i,str):
        return i
    if os.path.isfile(i):
        return i
    else:
        return None
def add_pairs(annote,task_dirp):
    """
    try add pairs file path, by default using:
        pairs_0/*sample*.pairs.gz
        pairs_c1/*sample*.c1.pairs.gz
        pairs_c12/*sample*.c12.pairs.gz
        pairs_c123/*sample*.c123.pairs.gz
        dip/*sample*.dip.pairs.gz
    assign None if file doesn't exist
    Input:
        annote: must use sample name as index
    Output:
        anntoe with extra cols
    """
    cell_annote = annote
    pwd = task_dirp
    pairs_0 = [os.path.join(os.path.join(pwd, "pairs_0"), name + ".pairs.gz") for name in cell_annote.index]
    pairs_c1 = [os.path.join(os.path.join(pwd, "pairs_c1"), name + ".c1.pairs.gz") for name in cell_annote.index]
    pairs_c12 = [os.path.join(os.path.join(pwd, "pairs_c12"), name + ".c12.pairs.gz") for name in cell_annote.index]
    pairs_c123 = [os.path.join(os.path.join(pwd, "pairs_c123"), name + ".c123.pairs.gz") for name in cell_annote.index]
    dip = [os.path.join(os.path.join(pwd, "dip"), name + ".dip.pairs.gz") for name in cell_annote.index]
    return cell_annote.assign(
        task_dirp = pwd,
        pairs_0 = [real_file(i) for i in pairs_0],
        pairs_c1 = [real_file(i) for i in pairs_c1],
        pairs_c12 = [real_file(i) for i in pairs_c12],
        pairs_c123 = [real_file(i) for i in pairs_c123],
        dip = [real_file(i) for i in dip]
    )
# add info
def update_values(base:dict, record:str)->dict:
    with open(record, 'r') as f:
        kv = json.load(f)
    for k in kv:
        try:
            base[k].update(kv[k])
        except KeyError:
            base.update(kv)
    return base
def listdir_bytime(folder):
    # list files in directory by time, old to new
    # assume zero layer in input folder(all files no folder)
    # return full path(include folder/ prefix)
    files = [os.path.join(folder, file) for file in os.listdir(folder)]
    files.sort(key= lambda fn: os.path.getmtime(fn))
    return files
def get_info(log_d):
    base = {}
    for fname in listdir_bytime(log_d):
        #print(fname)
        if fname.split(".")[-1] == "json":
            #print(fname)
            update_values(base, fname)
    #print(base)
    info = pd.DataFrame(base)
    return info
def get_info(log_d):
    base = {}
    for fname in listdir_bytime(log_d):
        #print(fname)
        if fname.split(".")[-1] == "json":
            #print(fname)
            update_values(base, fname)
    #print(base)
    info = pd.DataFrame(base)
    return info
def add_info(annote,rd):
    info = get_info(rd)
    return pd.concat([annote,info.T],axis=1)
# add RNA-seq
def add_umis(annote, umip):
    """
    Adding umi stat to annote.
    Assign 0 if sample not present in count matrix.
    """
    df = pd.read_table(umip,index_col=0,dtype={"gene":"string"})
    df.columns = df.columns.astype("string")
    umis = df.sum()
    umis.name = "umis"
    genes = (df > 0).astype(int).sum()
    genes.name = "genes"
    new_annote = pd.concat([annote,umis,genes],axis=1)
    return new_annote.fillna({"umis":0, "genes":0}, axis=0)
def add_extra(annote):
    """
    Calculate contact_per_reads, umis_per_reads, rna_ratio.
    Handling zero-division and NA values(and set to -1).
    Output:
        new df with extra cols.
    """
    new_annote = annote.assign(
        con_per_reads = -1,
        umis_per_reads = -1,
        rna_ratio = -1
    )
    valid_con_per_reads_row = (~annote["raw_contacts"].isna())&(~annote["dna_reads"].isna())&(annote["dna_reads"]!=0)
    new_annote.loc[valid_con_per_reads_row, "con_per_reads"] = annote.loc[valid_con_per_reads_row,"raw_contacts"]/annote.loc[valid_con_per_reads_row,"dna_reads"]
    valid_umis_per_reads_row = (~annote["umis"].isna())&(~annote["rna_reads"].isna())&(annote["rna_reads"]!=0)
    new_annote.loc[valid_umis_per_reads_row, "umis_per_reads"] = annote.loc[valid_umis_per_reads_row,"umis"]/annote.loc[valid_umis_per_reads_row,"rna_reads"]
    valid_rna_ratio_row = (~annote["rna_reads"].isna())&(~annote["raw_reads"].isna())&(annote["raw_reads"]!=0)
    new_annote.loc[valid_rna_ratio_row, "rna_ratio"] = annote.loc[valid_rna_ratio_row,"rna_reads"]/annote.loc[valid_rna_ratio_row,"raw_reads"]
    return new_annote
# agg
def task_stat(task_dirp,threads=32):
    """
    Get all infomation about this task
    """
    con = pd.read_csv(os.path.join(task_dirp,"contacts_info.csv"),index_col = 0, dtype={"name":"string"})
    annote = add_pairs(con,task_dirp)
    annote = add_pairs_num(annote,threads)
    annote = add_info(annote,os.path.join(task_dirp,"rd"))
    annote = add_umis(annote,os.path.join(task_dirp,"count_matrix","counts.gene.tsv.gz"))
    annote = add_extra(annote)
    return annote
# pick useful cols
def pick_useful(annote,mode="DNA_RNA"):
    if mode == "DNA_RNA":
        cols = ["raw_reads","dna_reads","raw_contacts","raw_intra","dup_rate","contacts","intra",
                "pairs_c1_num","pairs_c12_num","pairs_c123_num","phased_ratio",
               "rna_reads","umis","rna_ratio","umis_per_reads","con_per_reads"
               ]
    elif mode == "DNA":
        cols = ["raw_reads","dna_reads","raw_contacts","raw_intra","dup_rate","contacts","intra",
                "pairs_c1_num","pairs_c12_num","pairs_c123_num","phased_ratio","1m_rmsd","con_per_reads"]
    cols = annote.columns.intersection(cols)
    return annote[cols]