# things to do after bubbleflow
# add_pairs -> add_pairs_num -> add_info ... -> pick_useful

import gzip
import json
import math
import os
import re
from concurrent import futures
from genericpath import isfile
from itertools import repeat
from subprocess import check_output, CalledProcessError

import numpy as np
import pandas as pd
import pysam

formal_feature_names = {
    "con_per_reads" : "Contacts per read",
    "umis_per_reads" : "UMIs per read",
    "umis" : "UMIs",
    "20k_rmsd" : "20k RMSD",
    "pairs_c12_num" : "Contacts after clean1,2",
    "intra" : "Intra contact percent",
    "rna_ratio" : "RNA ratio",
    "raw_reads" : "Total raw reads",
    "rna_reads" : "RNA reads",
}
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
# count mapping rate
def mapping_rate(bam,threads=4):
    bamstat = json.loads(pysam.flagstat(bam, "-O","json","-@",str(threads)))
    #print(bamstat)
    return {"primary mapped":bamstat["QC-passed reads"]["primary mapped"]/bamstat["QC-passed reads"]["total"],
            "mapped":bamstat["QC-passed reads"]["mapped"]/bamstat["QC-passed reads"]["total"]}
# add_pairs path
def real_file(i):
    if not isinstance(i,str):
        return i
    if os.path.isfile(i):
        return i
    else:
        return None
def multi_re(key_array, re_array, base_array):
    """
    Search according to regular expressions. Used in log file info extraction.
    Print finding result and return finding as dict.
    Input:
        key_array: feature name of each re target.
        re_s_array: re strings.
        base_array: divide by elements to get real ratio.
    Output:
        searching result as dict.
    """
    re_array = [re.compile(i) for i in re_array]
    def multi_re_(logfile):
        value_array = []
        
        with open(logfile, "rt") as f:
            lines = f.readlines()
            for re_ in re_array:
                for line in lines:
                    re_res = re_.search(line)
                    if re_res is not None:
                        value_array.append(float(re_res.group(1)))
        res = dict(zip(key_array, np.array(value_array)/np.array(base_array)))
        for i, j in res.items():
            print("{} : {:.2f}".format(i,j))
        return res
    return multi_re_

star_stat = multi_re(["Input_reads", "Uniquely_mapped", "Unmapped_to_short", "Unmapped_mismatch"],
         ["Number of input reads\s+\W\s+(\d+)",
         "Uniquely mapped reads %\s+\W\s+(\d+\.\d+)%",
         "too short\s+\W\s+(\d+\.\d+)%",
         "too many mismatches\s+\W\s+(\d+\.\d+)%"],
         [1,100,100,100])
fc_stat = multi_re(["FeatureCount"],["assigned alignments : \d+ \((\d+\.\d+)%\)"],[100])
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
def add_info(annote,rd):
    info = get_info(rd)
    return pd.concat([annote,info.T],axis=1)
# add RNA-seq
def add_umis(annote, umip):
    """
    Adding umi stat to annote.
    Assign 0 if sample not present in count matrix.
    Input:
        umip: umi count matrix path
    """
    if not os.path.isfile(umip):
        print("add_umis: umi file %s not exist." % umip)
        annote = annote.assign(umis=0)
        annote = annote.assign(genes=0)
        return annote
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
        con_per_reads = math.nan,
        umis_per_reads = math.nan,
        rna_ratio = math.nan
    ).astype( # otherwise will be int and cause warning if fill float values
        {
            "con_per_reads" : float,
            "umis_per_reads" : float,
            "rna_ratio" : float
        }
    )
    if "raw_contacts" in annote.columns and "dna_reads" in annote.columns:
        valid_con_per_reads_row = (~annote["raw_contacts"].isna())&(~annote["dna_reads"].isna())&(annote["dna_reads"]!=0)
        new_annote.loc[valid_con_per_reads_row, "con_per_reads"] = (
            annote.loc[valid_con_per_reads_row,"raw_contacts"]/annote.loc[valid_con_per_reads_row,"dna_reads"]
            ).astype(float) # direct dividing will generate object dtype, dna_reads or raw_contacts must have been mixed by non-int values
    if "umis" in annote.columns and "rna_reads" in annote.columns:
        valid_umis_per_reads_row = (~annote["umis"].isna())&(~annote["rna_reads"].isna())&(annote["rna_reads"]!=0)
        new_annote.loc[valid_umis_per_reads_row, "umis_per_reads"] = (
            annote.loc[valid_umis_per_reads_row,"umis"]/annote.loc[valid_umis_per_reads_row,"rna_reads"]
            ).astype(float)
    if "rna_reads" in annote.columns and "raw_reads" in annote.columns:
        valid_rna_ratio_row = (~annote["rna_reads"].isna())&(~annote["raw_reads"].isna())&(annote["raw_reads"]!=0)
        new_annote.loc[valid_rna_ratio_row, "rna_ratio"] = (
            annote.loc[valid_rna_ratio_row,"rna_reads"]/annote.loc[valid_rna_ratio_row,"raw_reads"]
            ).astype(float)
    return new_annote
def add_mapping(annote, old=False, threads=16):
    """
    Adding mapping and primary mapping rate of DNA library.
    Input:
        old: sam version compatible to old pipeline. task_dirp/sam/xxx.aln.sam.gz
    """
    sthreads = 4
    if old:
        with futures.ProcessPoolExecutor(int(threads/sthreads)) as pool:
            res = pool.map(
                mapping_rate,
                annote["task_dirp"] + "/sam/" + annote.index + ".aln.sam.gz",
                repeat(sthreads, annote.shape[0])
            )
    else:
        with futures.ProcessPoolExecutor(int(threads/sthreads)) as pool:
            res = pool.map(
                mapping_rate,
                annote["task_dirp"] + "/hic_mapped/" + annote.index + ".sorted.bam",
                repeat(sthreads, annote.shape[0])
            )
    ares = pd.DataFrame(list(res), index = annote.index)
    return pd.concat([annote, ares],axis=1,join="outer")
def check_RNA(task_dirp):
    res = {}
    res.update(star_stat(os.path.join(task_dirp,"star_mapped","Log.final.out")))
    res.update(fc_stat(os.path.join(task_dirp, "count_gene","gene_assigned.summary")))
    for i, j in res.items():
        print("{} : {:.2f}".format(i,j))
# agg
def task_stat(task_dirp, ref=None, threads=32)->pd.DataFrame:
    """
    Get all infomation about this bubble_flow task in the form of a dataframe.
    The meta dataframe is constructed by adding additional columns to "contacts_info.csv" generated by bubble_flow.
    TODO:
        make decent dtypes for every column step by step.
    """
    # read base dataframe from "contacts_info.csv"
    con = pd.read_csv(os.path.join(task_dirp,"contacts_info.csv"),index_col = 0, dtype={"name":"string"})
    # add Hi-C related info
    annote = add_pairs(con,task_dirp)
    annote = add_pairs_num(annote,threads)
    annote = add_info(annote,os.path.join(task_dirp,"rd"))
    # add RNA related info
    # add:
    #   con_per_reads: raw contacts / dna_reads
    #   umis_per_reads: umis / rna_reads
    #   rna_ratio: rna_reads / raw_reads
    if ref is not None:
        annote = add_umis(annote,os.path.join(task_dirp,"count_matrix"+"_"+ref,"counts.gene.tsv.gz"))
    else:
        annote = add_umis(annote,os.path.join(task_dirp,"count_matrix","counts.gene.tsv.gz"))
    annote = add_extra(annote)
    return annote
# pick useful cols
def pick_useful(annote,features=["basic"]):
    """
    Pick useful cols from meta dataframe.
    Input:
        annote: meta dataframe
        features: list of features to pick, default words are "DNA" and "DNA_RNA"
    Output:
        new df with useful cols.
    """
    real_features = []
    for i in features:
        if i == "DNA_RNA":
            real_features.extend(["raw_reads","dna_reads","raw_contacts","raw_intra","dup_rate","contacts","intra",
                "pairs_c1_num","pairs_c12_num","pairs_c123_num","phased_ratio",
               "rna_reads","umis","rna_ratio","umis_per_reads","con_per_reads"
               ])
        elif i == "DNA":
            real_features.extend(["raw_reads","dna_reads","raw_contacts","raw_intra","dup_rate","contacts","intra",
                "pairs_c1_num","pairs_c12_num","pairs_c123_num","phased_ratio","1m_rmsd","con_per_reads"])
        elif i == "basic":
            real_features.extend(["raw_contacts","raw_intra","dup_rate","contacts","intra","phased_ratio",
                "pairs_c1_num","pairs_c12_num","pairs_c123_num"])
        elif i == "basic_struct":
            real_features.extend(["raw_contacts","raw_intra","dup_rate","contacts","intra","phased_ratio",
                "pairs_c1_num","pairs_c12_num","pairs_c123_num","1m_rmsd","200k_rmsd","50k_rmsd","20k_rmsd"])
        else:
            real_features.append(i)
    cols = annote.columns.intersection(real_features)
    return annote[cols]