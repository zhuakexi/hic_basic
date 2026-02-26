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
from pathlib import Path
from subprocess import check_output, CalledProcessError

import numpy as np
import pandas as pd
import pysam
from tqdm import tqdm

from ..hicio import read_umi_tools
from ..pipeline.rule import bubble_flow_touched

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
    """
    Count the number of lines matching a target string in a gzip-compressed file.

    Parameters
    ----------
    filename : str
        Path to the gzip-compressed file.
    target : str
        String pattern to search for (passed to ``zgrep -c``).

    Returns
    -------
    int
        Number of matching lines; 0 if none found or file is absent.
    """
    try:
        output_bytes = check_output([
            "zgrep","-c",target,
            filename])
    except CalledProcessError:
        return 0
    return int(output_bytes.decode("utf-8").strip())
def count_pairs(filename:str)->int:
    """
    Count the number of contact records in a 4DN pairs file.

    Lines beginning with ``#`` are treated as header/comment lines and
    excluded from the count.

    Parameters
    ----------
    filename : str
        Path to a gzip-compressed ``.pairs.gz`` file.

    Returns
    -------
    int or None
        Number of contact records, or ``None`` if *filename* is not a string
        or the file does not exist.
    """
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
    """
    Compute the primary-mapped and overall mapped rates for a BAM file.

    Parameters
    ----------
    bam : str or None
        Path to the sorted BAM file. If ``None``, returns zeros.
    threads : int, optional
        Number of threads for ``pysam.flagstat`` (default 4).

    Returns
    -------
    dict
        Dictionary with keys ``"primary mapped"`` and ``"mapped"``, each
        holding the corresponding fraction of QC-passed reads.
    """
    if bam is None:
        return {"primary mapped":0, "mapped":0}
    bamstat = json.loads(pysam.flagstat(bam, "-O","json","-@",str(threads)))
    #print(bamstat)
    return {"primary mapped":bamstat["QC-passed reads"]["primary mapped"]/bamstat["QC-passed reads"]["total"],
            "mapped":bamstat["QC-passed reads"]["mapped"]/bamstat["QC-passed reads"]["total"]}
# add_pairs path
def real_file(i):
    """
    Return the path if it points to an existing file, otherwise ``None``.

    Parameters
    ----------
    i : str or any
        A file path string to validate. Non-string values are returned as-is.

    Returns
    -------
    str or None
        The original path if the file exists; ``None`` otherwise.
    """
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
    """
    Update a base dictionary with values from a JSON record file.
    
    For each key in the record, if the key exists in base, update the nested dictionary.
    If the key doesn't exist, add the entire record to base.
    
    Args:
        base: Base dictionary to be updated
        record: Path to JSON file containing key-value pairs to merge
    
    Returns:
        Updated base dictionary
    """
    with open(record, 'r') as f:
        kv = json.load(f)
    for k in kv:
        try:
            base[k].update(kv[k])
        except KeyError:
            base.update(kv)
    return base

def listdir_bytime(folder):
    """
    List all files in a directory sorted by modification time (oldest first).
    
    Args:
        folder: Path to directory to scan
        
    Returns:
        List of full file paths sorted by modification time
    """
    # list files in directory by time, old to new
    # assume zero layer in input folder(all files no folder)
    # return full path(include folder/ prefix)
    files = [os.path.join(folder, file) for file in os.listdir(folder)]
    files.sort(key= lambda fn: os.path.getmtime(fn))
    return files

def get_info(log_d):
    """
    Extract and combine information from all JSON files in a directory.
    
    Processes JSON files in chronological order, merging their contents.
    Each JSON file should contain sample information as key-value pairs.
    
    Args:
        log_d: Directory containing JSON files with sample information
        
    Returns:
        DataFrame where columns are sample IDs and rows are attributes
    """
    base = {}
    for fname in listdir_bytime(log_d):
        #print(fname)
        if fname.split(".")[-1] == "json":
            #print(fname)
            update_values(base, fname)
    #print(base)
    info = pd.DataFrame(base)
    return info

def add_info(annote, rd):
    """
    Add sample information from JSON files to an existing annotation DataFrame.
    
    Args:
        annote: Existing annotation DataFrame
        rd: Directory containing JSON files with additional sample information
        
    Returns:
        Combined DataFrame with original annotations and new sample information
    """
    info = get_info(rd)
    return pd.concat([annote,info.T],axis=1)
# add RNA-seq
def add_umis(annote, umip):
    """
    Adding umi stat to annote.
    Accept a single umi count matrix path or an iterable (list/tuple) of umi count matrix paths.
    Missing files will be skipped with a warning. If no valid files are found, assign 0 for both
    "umis" and "genes" for all samples.

    Input:
        umip: umi count matrix path or iterable of paths
    """
    # Normalize input into a list of candidate paths
    paths = []
    if umip is None:
        paths = []
    elif isinstance(umip, (list, tuple, set)):
        paths = list(umip)
    elif isinstance(umip, str):
        paths = [umip]
    else:
        # try to coerce any iterable-like into a list
        try:
            paths = list(umip)
        except Exception:
            paths = [umip]

    existing = []
    for p in paths:
        if p is None:
            continue
        if os.path.isfile(p):
            existing.append(p)
        else:
            # keep behavior but allow multiple files: warn and skip
            print("add_umis: umi file %s not exist. Skipping." % p)

    # If no UMI files found, keep zeros
    if len(existing) == 0:
        annote = annote.assign(umis=0)
        annote = annote.assign(genes=0)
        return annote

    # Prepare sample index from annotation
    sample_index = pd.Index(annote.index.astype("string"))
    
    # Initialize accumulators for umis and genes
    total_umis = pd.Series(0, index=sample_index)
    total_genes = pd.Series(0, index=sample_index)  # Direct gene count without deduplication
    
    any_valid_data = False
    
    # Process each count matrix separately
    for p in existing:
        try:
            # Read the count matrix using custom function (samples as index, genes as columns)
            df_raw = read_umi_tools(p)
        except Exception as e:
            print(f"add_umis: error reading {p}: {str(e)}. Skipping.")
            continue
            
        # Skip empty matrices
        if df_raw.empty:
            print(f"add_umis: empty matrix in file {p}. Skipping.")
            continue
            
        # Transpose to have genes as rows and samples as columns
        df = df_raw.T
        
        # Ensure proper data types
        df.index = df.index.astype("string")  # genes as string index
        df.columns = df.columns.astype("string")  # samples as string columns
        
        # Find overlapping samples
        overlapping_cols = df.columns.intersection(sample_index)
        ignored_cols = df.columns.difference(sample_index)
        
        # Handle non-overlapping cases
        if overlapping_cols.empty:
            print(f"add_umis: no overlapping samples in file {p}. Skipping.")
            continue
            
        any_valid_data = True
        
        if len(ignored_cols) > 0:
            print(f"add_umis: ignoring count-matrix columns not present in annotation: {', '.join(ignored_cols)}")
        
        # Process only overlapping samples
        df = df.loc[:, overlapping_cols]
        
        # Ensure numeric values (in case read_umi_tools returned non-numeric types)
        df = df.apply(pd.to_numeric, errors='coerce').fillna(0)
        
        # Calculate UMIs for this matrix (sum per sample)
        current_umis = df.sum(axis=0)
        # Accumulate UMIs using vectorized operation
        total_umis.loc[overlapping_cols] += current_umis
        
        # Count expressed genes (>0) for each sample in this matrix
        current_genes = (df > 0).sum(axis=0)
        # Accumulate gene counts without deduplication
        total_genes.loc[overlapping_cols] += current_genes
    
    # If no valid data was found in any matrix
    if not any_valid_data:
        print("add_umis: no overlapping samples between annotation and any count matrix. Filling zeros.")
        annote = annote.assign(umis=0)
        annote = annote.assign(genes=0)
        return annote
    
    # Create result DataFrame
    new_annote = annote.copy()
    new_annote["umis"] = total_umis
    new_annote["genes"] = total_genes
    
    return new_annote.fillna({"umis": 0, "genes": 0})
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

def add_mapping(annote, threads=16):
    """
    Add DNA-library mapping and primary-mapping rate columns to annotation.

    BAM paths are resolved via ``bubble_flow_touched()["hic_mapped"]``.  Samples
    whose BAM file is missing are processed with ``mapping_rate(None)`` and
    receive zero rates.

    Parameters
    ----------
    annote : pandas.DataFrame
        Annotation DataFrame indexed by sample name; must include a
        ``task_dirp`` column.
    threads : int, optional
        Total worker threads (divided among parallel ``pysam.flagstat`` calls,
        default 16).

    Returns
    -------
    pandas.DataFrame
        Original DataFrame with two new columns appended: ``"primary mapped"``
        and ``"mapped"``.
    """
    sthreads = 4
    file_paths = []

    # 根据 old 参数生成对应的文件路径
    # for i in range(annote.shape[0]):
    #     if old:
    #         path = annote.iloc[i]["task_dirp"] + "/sam/" + str(annote.index[i]) + ".aln.sam.gz"
    #     else:
    #         path = annote.iloc[i]["task_dirp"] + "/hic_mapped/" + str(annote.index[i]) + ".sorted.bam"
    #     file_paths.append(path)
    for sample_name, row in annote.iterrows():
        file_path = bubble_flow_touched()["hic_mapped"][0].format(
            task_dirp=row["task_dirp"],
            sample_name=sample_name
        )
        if os.path.isfile(file_path):
            file_paths.append(file_path)
        else:
            print(f"Warning: File {file_path} does not exist. Skipping this sample.")
            file_paths.append(None)

    # 使用 ProcessPoolExecutor 提交任务
    with futures.ProcessPoolExecutor(max_workers=int(threads / sthreads)) as pool:
        future_to_index = {}
        for i, file_path in enumerate(file_paths):
            future = pool.submit(mapping_rate, file_path, sthreads)
            future_to_index[future] = i

        # 使用 tqdm 显示进度条
        with tqdm(total=len(future_to_index), desc="Processing Mapping") as pbar:
            results = [None] * len(future_to_index)
            for future in futures.as_completed(future_to_index):
                idx = future_to_index[future]
                try:
                    result = future.result()
                    results[idx] = result
                except Exception as e:
                    print(f"Error in task {idx}: {e}")
                    results[idx] = None  # 或根据需求处理异常
                finally:
                    pbar.update(1)

    # 构造结果 DataFrame 并合并
    ares = pd.DataFrame(results, index=annote.index)
    return pd.concat([annote, ares], axis=1, join="outer")
def check_RNA(task_dirp):
    """
    Print RNA-seq QC statistics for a bubble_flow task directory.

    Reads STAR alignment summary (``star_mapped/Log.final.out``) and
    featureCounts assignment summary (``count_gene/gene_assigned.summary``)
    from *task_dirp*, then prints and returns the combined statistics.

    Parameters
    ----------
    task_dirp : str
        Root directory of a bubble_flow task.
    """
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
    # ref can be:
    #   None -> use task_dirp/count_matrix/counts.gene.tsv.gz
    #   str  -> use task_dirp/count_matrix_<ref>/counts.gene.tsv.gz
    #   list/tuple/set -> multiple refs, each will be looked up in task_dirp/count_matrix_<ref>/counts.gene.tsv.gz
    umip_paths = []
    if ref is None:
        umip_paths = [os.path.join(task_dirp, "count_matrix", "counts.gene.tsv.gz")]
    elif isinstance(ref, (list, tuple, set)):
        for r in ref:
            umip_paths.append(os.path.join(task_dirp, "count_matrix_" + str(r), "counts.gene.tsv.gz"))
    else:
        umip_paths = [os.path.join(task_dirp, "count_matrix_" + str(ref), "counts.gene.tsv.gz")]

    # add_umis accepts a single path or an iterable of paths; it will skip missing files and combine counts
    annote = add_umis(annote, umip_paths)
    annote = add_extra(annote)
    return annote
def task_stat_rna_wgs(task_dirp: str, ref: str) -> pd.DataFrame:
    """
    Build a QC summary DataFrame for a combined RNA + WGS bubble_flow task.

    Reads sample metadata from ``<task_dirp>/rd``, filters to a fixed set of
    QC columns (mapping rates, read counts, phasing scores, etc.), appends UMI
    and gene counts from the specified count matrix, and calculates derived
    metrics via :func:`add_extra`.

    Parameters
    ----------
    task_dirp : str
        Root directory of the bubble_flow task.
    ref : str
        Reference genome tag used to locate the count matrix at
        ``<task_dirp>/count_matrix_<ref>/counts.gene.tsv.gz``.

    Returns
    -------
    pandas.DataFrame
        Per-sample QC DataFrame with columns for mapping rates, read counts,
        phasing metrics, UMIs, genes, and derived ratios.
    """
    meta = get_info(Path(task_dirp) / "rd").T
    required_cols = [
        'R1_file', 'R2_file', 'dna_reads', 'rna_reads', 'rna_c1_reads',
        'rna_c2_reads', 'primary mapped', 'mapped', 'raw_reads', 'hap1_phased',
        'hap2_phased', 'biasedX_score', 'hap_score', 'ypercent', 'xpercent',
        'cell_state'
    ]
    dtypes = {
        "R1_file": str,
        "R2_file": str,
        "dna_reads": int,
        "rna_reads": int,
        "rna_c1_reads": int,
        "rna_c2_reads": int,
        "primary mapped": float,
        "mapped": float,
        "raw_reads": int,
        "hap1_phased": float,
        "hap2_phased": float,
        "biasedX_score": float,
        "hap_score": float,
        "ypercent": float,
        "xpercent": float,
        "cell_state": str
    }
    # Check which columns exist in the metadata
    available_cols = [col for col in required_cols if col in meta.columns]
    missing_cols = [col for col in required_cols if col not in meta.columns]
    # Print warnings for missing columns
    if missing_cols:
        print(f"Warning: The following columns are missing from the metadata and will be skipped: {missing_cols}")
    meta = meta.loc[
        :,
        available_cols
        ].astype({col: dtypes[col] for col in available_cols if col in dtypes})
    cm_path = bubble_flow_touched()["count_matrix"][0].format(
        task_dirp=task_dirp, ref=ref
    )
    meta_e = add_umis(meta, cm_path)
    meta_e = add_extra(meta_e)
    return meta_e
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