# cooler api functions
import os
import os.path as op
import tempfile
import yaml
from pathlib import Path
from typing import List, Union

import cooler
import numpy as np
import pandas as pd
from cooler import Cooler
from cooler.create import sanitize_records, aggregate_records
from cytoolz import compose # fail in HPC

from .binnify import GenomeIdeograph
from .data import chromosomes
from .hicio import schicluster2mat
from .utils import binnify

# from .utils import read_chromsizes
# def gen_bins(chromsizes_file:str, binsize:int):
#     """
#     Generate "bins" for cooler api.
#     """
#     if not op.exists(chromsizes_file):
#         raise ValueError('File "{}" not found'.format(chromsizes_file))
#     try:
#         binsize = int(binsize)
#     except ValueError:
#         raise ValueError(
#             'Expected integer binsize argument (bp), got "{}"'.format(binsize)
#         )
#     chromsizes = read_chromsizes(chromsizes_file, all_names=True)
#     bins = binnify(chromsizes, binsize)
#     return bins


# --- warm it --- #


def str2slice(clr, region:str):
    """
    Convert region string to slice object.
    Input:
        clr: cooler file object
        region: region string, eg. "chr1:1,000,000-2,000,000"
    Output:
        slice object
    """
    return slice(*clr.extent(region))
def iter_pixels(clr, chunksize=1e6, join=True):
    chunksize = int(chunksize)
    tot = clr.pixels().shape[0]
    for i in range(0, tot, chunksize):
        yield clr.pixels(join=join)[i:i+chunksize]
def cool2pixels(coolp:str, balance:bool=True)->pd.DataFrame:
    """
    Fetch all pixels from a cooler file.
    Input:
        coolp: path to cooler file
        balance: whether to include balanced values
    Output:
        pixels: DataFrame of pixels, [bin1_id, bin2_id, count, balanced]
    """
    clr = cooler.Cooler(str(coolp))
    pixels = clr.pixels()[:]
    bins = clr.bins()[:]
    pixels = cooler.annotate(pixels, bins, replace=True)
    if balance:
        pixels = pixels.assign(
            balanced = pixels.eval('weight1 * weight2 * count')
        )
    return pixels
def cool2mat(cool, region:Union[str, List[str], slice, List[slice]], balance:bool=False, min_nz=None,bad_bins=None)->pd.DataFrame:
    """
    Fetch matrix from cooler file with proper index and columns.
    Input:
        cool: cooler file path
        region: genome region to fetch.
            "chr1" or "chr1:1000000-2000000":
                ucsc styl, get intra-chrom matrix
            ["chr1:1,000,000-2,000,000", "chr2"]:
                list of region string to get inter-chrom matrix
            slice(0,-1):
                slicer syntax, treat the whole genome as a big matrix,get intra-chrom matrix
                Note: can only cross chrom boundaries when using slicer syntax
                For example slice(0,-1) for whole genome.
            [slice(0,1000000), slice(1000000,2000000)]:
                inter-chrom matrix using slicer syntax
        balance: balance or not, bool
        min_nz: minimum number of non-zero pixels to keep, bin with less pixels will be set to nan
            if None, no filtering
        bad_bins: list of bad bins to be set to nan, a dataframe with columns ["chrom","start"]
    Output:
        pd.DataFrame
    TODO:
        min_nz and bad_bins support for separate index(region1) and columns(region2) filtering for inter-chrom matrix
    """
    no_fetch = False # whether to use.fetch method, if False, using slicer syntax
    clr = cooler.Cooler(str(cool))
    if isinstance(region, str): # "chr1" or "chr1:1000000-2000000"
        region = [region]
    elif isinstance(region, list):
        assert len(region) == 2, "Only support 2 regions for inter-chrom matrix"
        if isinstance(region[0], slice) and isinstance(region[1], slice): # [slice(0,1000000), slice(1000000,2000000)]
            no_fetch = True
        else: # ["chr1:1,000,000-2,000,000", "chr2"], ["chr1", slice(0,1000000)]
            region = [region if isinstance(region, slice) else str2slice(clr, region) for region in region]
            no_fetch = True
    elif isinstance(region, slice): # slice(0,-1)
        region = [region]
        no_fetch = True
    else:
        raise ValueError("[Error in cool2mat]: Region must be str or list or slice, read doc for more info.")
    if no_fetch:
        # use slicer syntax
        mat = clr.matrix(balance=balance)[region[0],region[0]] if len(region) == 1 else clr.matrix(balance=balance)[region[0], region[1]] # because can't unpack in slicer
        mat = pd.DataFrame(mat)
        if len(region) == 1:
            index = clr.bins()[region[0]]
            columns = clr.bins()[region[0]]
        else:
            # inter-chromosome region
            index = clr.bins()[region[0]] # first region as index, same as cooler matrix fetch
            columns = clr.bins()[region[1]]
    else:
        mat = clr.matrix(balance=balance).fetch(*region)
        mat = pd.DataFrame(mat)
        if len(region) == 1:
            index = clr.bins().fetch(region[0])
            columns = clr.bins().fetch(region[0])
        else:
            # inter-chromosome region
            index = clr.bins().fetch(region[0]) # first region as index, same as cooler matrix fetch
            columns = clr.bins().fetch(region[1])
    # mat.columns = columns["start"]
    # mat.columns.name = "region 2"
    # mat.index = index["start"]
    # mat.index.name = "region 1"
    if min_nz is not None:
        index_coverage = (mat > 0).sum(axis=1)
        columns_coverage = (mat > 0).sum(axis=0)
        mat.loc[
            index_coverage[index_coverage < min_nz].index,
            :
            ] = np.nan
        mat.loc[
            :,
            columns_coverage[columns_coverage < min_nz].index
            ] = np.nan
    mat.columns = pd.MultiIndex.from_frame(columns[["chrom","start"]],names=["chrom","start"])
    mat.columns.name = "region 2"
    mat.index = pd.MultiIndex.from_frame(index[["chrom","start"]],names=["chrom","start"])
    mat.index.name = "region 1"
    if bad_bins is not None:
        bad_bins_index = pd.MultiIndex.from_frame(
            bad_bins[["chrom", "start"]]
        )
        # what if intersection is empty?
        mat.loc[bad_bins_index.intersection(mat.index), :] = np.nan
        mat.loc[:, bad_bins_index.intersection(mat.columns)] = np.nan
    return mat


### --- cool it --- ###


def gen_bins(genome, binsize):
    """
    Generate "bins" for cooler api.
    """
    try:
        binsize = int(binsize)
    except ValueError:
        raise ValueError(
            'Expected integer binsize argument (bp), got "{}"'.format(binsize)
        )
    chromsizes = chromosomes(genome)["length"]
    bins = binnify(chromsizes, binsize)
    return bins

def _easy_pixels(pairs_path, bins, chunksize=15e6, zero_based=False, tril_action="reflect"):
    """
    Easy to use func generating "pixels" for cooler api.
    """
    # --- pairs chunk read ---
    # give reader
    input_field_names = [
        'chrom1', 'pos1', 'chrom2', 'pos2',
    ]
    input_field_dtypes = {
        'chrom1': str,
        'pos1': np.int64,
        'chrom2': str,
        'pos2': np.int64,
    }
    input_field_numbers = {"chrom1":1,"pos1":2,"chrom2":3,"pos2":4}
    f_in = pairs_path

    reader = pd.read_csv(
        f_in,
        sep='\t',
        usecols=[input_field_numbers[name] for name in input_field_names],
        names=input_field_names,
        dtype=input_field_dtypes,
        iterator=True,
        chunksize=chunksize,
        comment="#"
        )
    # --- "online" binnify ---
    sanitize = sanitize_records(
        bins,
        schema='pairs',
        decode_chroms=True,
        is_one_based=not zero_based,
        tril_action=tril_action,
        sort=True,
        validate=True
        )
    aggregations = {}
    aggregate = aggregate_records(agg=aggregations, count=True, sort=False)
    pipeline = compose(aggregate, sanitize)
    return map(pipeline, reader)
def pairs2cool(pairs_path, coolpath, sizef, binsize):
    """
    Generate .cool file from 4DN .pairs file
    """
    bins = gen_bins(sizef, binsize)
    cooler.create_cooler(
        coolpath,
        bins,
        _easy_pixels(pairs_path, bins)
    )
def pixels2cool(pixels:pd.DataFrame, bins:pd.DataFrame, outfile, columns:list=["count"], dtypes:dict={"count":int}):
    """
    Dump pixels to cool file.
    Input:
        pixels: pixels dataframe or iterator of pixels dataframes,
            [chrom1, start1, (end1), chrom2, start2, (end2), vcol1, vcol2, ...]
        bins: [chrom, start, end]
        outfile: output cool file path
        columns: value columns to be stored
        dtypes: data type for each value column
    Output:
        outfile
    TODO: deal when bin1_id or bin2_id in pixels
    """
    if isinstance(pixels, pd.DataFrame):
        # single dataframe
        pixels = [pixels]
    new_bins = bins.assign(
        bin_id = bins.index
    )
    # add bin1_id
    new_bins1 = new_bins.copy()
    new_bins1.columns = ["chrom1", "start1", "end1", "bin1_id"]
    pixel_chunks_id1 = (
        pd.merge(
            pixel_chunk,
            new_bins1[['chrom1', 'start1', 'end1', 'bin1_id']],
            on=["chrom1", "start1", "end1"] if "end1" in pixel_chunk.columns else ["chrom1", "start1"],
        )
        for pixel_chunk in pixels
    )
    # add bin2_id
    new_bins2 = new_bins.copy()
    new_bins2.columns = ["chrom2", "start2", "end2", "bin2_id"]
    pixel_chunks_id1_id2 = (
        pd.merge(
            pixel_chunk,
            new_bins2[['chrom2', 'start2', 'end2', 'bin2_id']],
            on=["chrom2", "start2", "end2"] if "end2" in pixel_chunk.columns else ["chrom2", "start2"],
        )
        for pixel_chunk in pixel_chunks_id1
    )
    # prepare inputs
    new_bins = new_bins[["chrom", "start", "end"]]
    pixel_chunks_tidy = (
        pixel_chunk[[
            "bin1_id", "bin2_id", "count"
        ]].dropna(
            axis=0,
            how="any",
            subset=["bin1_id", "bin2_id"]
        ) # drop pixels not present in new_bins
        for pixel_chunk in pixel_chunks_id1_id2
    )
    cooler.create_cooler(
        outfile,
        new_bins,
        pixel_chunks_tidy,
        dtypes = dtypes,
        ordered = False, # doesn't ensure order
        symmetric_upper = False, # don't assume symmetric
        triucheck = False,
        dupcheck = True,
        boundscheck = True
    )
# schicluster io
def _schicluster_pixels(filesp, genome, binsize):
    """
    Generate pixels iterator for cooler constructor.
    Input:
        filesp: All hdf5 file paths for 1 cell, chrom as index; pd.Series
        genome: used to create shifts, eg. "mm10"
        binsize: binsize of original schicluster matrix
    Output:
        iterator of dict
    """
    shifts = gen_bins(genome,binsize)["chrom"].drop_duplicates(keep="first")
    shifts = {shifts[key]:key for key in shifts.index}
    for chrom in filesp.index:
        coo = schicluster2mat(filesp[chrom]).tocoo()
        yield {
            "bin1_id" : coo.row + shifts[chrom],
            "bin2_id" : coo.col + shifts[chrom],
            "count" : coo.data
        }
def schicluster2cool(filesp, fo, genome, binsize,
    ordered=True,ensure_sorted=True, **args):
    """
    Dump schicluster imputed matrix to .cool file.
    Input:
        filesp: All hdf5 file paths for 1 cell, chrom as index; pd.Series
        fo: output cool uri
        genome: eg. "mm10"
        binsize: binsize of original schicluster matrix
        **args: other args for cooler.create_cooler
    Output:
        fo
        write to file and return file path
    Example:
    ```
        from pathlib import Path
        import pandas as pd
        from hic_basic.data import chromosomes
        from hic_basic.coolstuff import schicluster2cool
        
        
        sample = "2021110112"
        ddir = "/shareb/ychi/repo/embryo_integrate/schicluster/imputed_matrix/100000/"

        chroms = [chrom for chrom in chromosomes(genome).index if chrom not in ["chrX","chrY"]]
        filesp = pd.Series(
            data = [
                Path(ddir)/ "{chrom}".format(chrom=chrom) / "{sample}_{chrom}_pad1_std1_rp0.5_sqrtvc.hdf5".format(sample=sample, chrom=chrom) 
                for chrom in chroms
            ],
            index = chroms,
            name = sample)
        
        schicluster2cool(filesp, "test.cool","mm10",100e3)
    ```
    """
    bins = gen_bins(
        genome,
        binsize
    )
    cooler.create_cooler(
        fo,
        bins,
        _schicluster_pixels(filesp, genome, binsize),
        dtypes = {
            "count":float # hicluster impute output is float
        },
        ordered = ordered,
        ensure_sorted = ensure_sorted,
        **args
    )
    return fo
# gam io
def parse_gam(gam_file, pixel=False):
    """
    Input:
        gam_file: path to gam file
            temporarily, only support .npz file
        pixel: a cooler's joint pixel format
    Output:
        mat: pandas.DataFrame
            if not pixel, a DataFrame with 2-level-multiindex (chrom, start) index and columns
    """
    npz = np.load(gam_file)
    if not pixel:
        index = pd.MultiIndex.from_frame(
            pd.DataFrame(
                {
                    "chrom" : pd.Series(npz["windows_0"][:,0], dtype="string"),
                    "start" : pd.Series(npz["windows_0"][:,1], dtype="int")
                }
                )
        )
        columns = pd.MultiIndex.from_frame(
            pd.DataFrame(
                {
                    "chrom" : pd.Series(npz["windows_1"][:,0], dtype="string"),
                    "start" : pd.Series(npz["windows_1"][:,1], dtype="int")
                }
                )
        )
        mat = pd.DataFrame(
            npz["scores"],
            index=index,
            columns=columns,
            dtype = "float"
        )
    else:
        values = npz["windows_0"]
        pos1 = pd.DataFrame(
            {
                "chrom1" : pd.Series(values[:,0], dtype="string"),
                "start1" : pd.Series(values[:,1], dtype="int"),
                "end1" : pd.Series(values[:,2], dtype="int")
            }
        )
        values = npz["windows_1"]
        pos2 = pd.DataFrame(
            {
                "chrom2" : pd.Series(values[:,0], dtype="string"),
                "start2" : pd.Series(values[:,1], dtype="int"),
                "end2" : pd.Series(values[:,2], dtype="int")
            }
        )
        mat = pd.DataFrame(
            npz["scores"],
            dtype = "float"
            ).stack().reset_index()
        mat.columns = ["sub_id1","sub_id2","count"]
        mat = pd.merge(
            mat,
            pos1,
            left_on="sub_id1",
            right_index=True
        )
        mat = pd.merge(
            mat,
            pos2,
            left_on="sub_id2",
            right_index=True
        )
        mat = mat[
            ["chrom1","start1","end1","chrom2","start2","end2","count"]
            ]
    return mat
def gam2cool(gam_filepat, bin_file, outfile, force=False):
    """
    Transform gam result to cool format
    Input:
        gam_filepat: path pattern to gam file, use * to represent any string,
            typically chrom names
        bin_file: path to bin file created by bedtools in the pipeline
        outfile: path to output cool file
    Output:
        outfile
    """
    if Path(outfile).exists() and not force:
        return outfile
    gam_filepat = Path(gam_filepat)
    gam_files = list(gam_filepat.parent.glob(gam_filepat.name))
    if len(gam_files) == 0:
        raise ValueError("No gam files found.")
    gam_pixels = (parse_gam(gam_file, pixel=True) for gam_file in gam_files)
    bins = pd.read_table(bin_file, names=["chrom","start","end"])
    pixels2cool(gam_pixels, bins, outfile, columns=["count"], dtypes={"count":"float"})
    return outfile
# gam_filepat = "/sharec/ychi/repo/tillsperm77_gam/matrix/1000000/dprime/QC_passed.*__*.npz"
# outfile = "/sharec/ychi/repo/tillsperm77_gam/cool/QC_passed.1m.cool"
# bin_file = "/sharec/ychi/repo/tillsperm77_gam/multicov/GRCh38.1000000.bed"
def _scool_pixels(pairs_paths:dict, bins, chunksize=15e6, zero_based=False, tril_action="reflect"):
    """
    Generate pixels dict for multi-pairs-file. Using in scool generating.
    Return:
        Cell name as key and pixel table DataFrame as value
    """
    n = len(pairs_paths)
    return {pairs_path : _easy_pixels(pairs_paths[pairs_path], bins, chunksize,zero_based, tril_action) 
        for pairs_path in pairs_paths}
# scool io
def pairs2scool(pairs_paths, coolpath, sizef, binsize):
    """
    Generate .scool file from multiple 4DN .pairs file
    Input:
        pairs_path: pairs file path; dict
            {cell_name: pairs_path}
        coolpath: output .scool file path; str
        sizef: chrom size file path, 2col(chrom:str, size:int) tsv without header
        binsize: bin size; int
    """
    bins = gen_bins(sizef, binsize)
    cooler.create_scool(
        coolpath,
        bins,
        _scool_pixels(
            pairs_paths,
            bins
        )
    )


### --- rearranging your fridge --- ###


def cools2scool(cools, scool_path):
    """
    Generate scool from cooler files.
    TODO:
        don't readin all cools together.
    Input:
        cools: dict of name-cooler file, using first file's bins as default bins
        scool_path: outputfile name
    """
    # pick first cool
    bins = cooler.Cooler(cools[next(iter(cools))]).bins()[:]
    multif_pixels = {
        name : cooler.Cooler(cools[name]).pixels()[:]
        for name in cools
    }
    cooler.create_scool(
        scool_path,
        bins,
        multif_pixels
    )
def hic_pileup(scool, grouping, cache_pattern="{}.pileup.cool",mergebuf=1e6):
    """
    Generate pileup cool file for each cluster.
    To simply pileup all cells, use `grouping = pd.Series(1, index=cells)`
    Input:
        scool: cooler's scool file path
        grouping: dict, key is cluster name, value is list of cell names
        cache_pattern: path to store resulting .cool file; str with {}
    Output:
        write cool according to cache_pattern
    """
    for cluster in grouping:
        samples = grouping[cluster]
        print("Merging " + str(cluster))
        cooler.merge_coolers(
            cache_pattern.format(str(cluster)),
            [scool+"::/cells/{}".format(i) for i in samples],
            mergebuf
        )
def reset_cool_bins(coolpin, coolpout, genome="GRCh38", chunksize=1e6):
    clr = cooler.Cooler(coolpin)
    binsize = clr.binsize
    new_bins = GenomeIdeograph(
        genome).bins(
            binsize,
            bed=True,
            order=True
        )
    new_bins = new_bins.assign(
        bin_id = new_bins.index
    )
    # add bin1_id
    new_bins1 = new_bins.copy()
    new_bins1.columns = ["chrom1", "start1", "end1", "bin1_id"]
    pixel_chunks_id1 = (
        pd.merge(
            pixel_chunk,
            new_bins1[['chrom1', 'start1', 'end1', 'bin1_id']],
            on=["chrom1", "start1", "end1"],
        )
        for pixel_chunk in iter_pixels(
            clr, chunksize=chunksize, join=True)
    )
    # add bin2_id
    new_bins2 = new_bins.copy()
    new_bins2.columns = ["chrom2", "start2", "end2", "bin2_id"]
    pixel_chunks_id1_id2 = (
        pd.merge(
            pixel_chunk,
            new_bins2[['chrom2', 'start2', 'end2', 'bin2_id']],
            on=["chrom2", "start2", "end2"],
        )
        for pixel_chunk in pixel_chunks_id1
    )
    # prepare inputs
    new_bins = new_bins[["chrom", "start", "end"]]
    pixel_chunks_tidy = (
        pixel_chunk[[
            "bin1_id", "bin2_id", "count"
        ]].dropna(
            axis=0,
            how="any",
            subset=["bin1_id", "bin2_id"]
        ) # drop pixels not present in new_bins
        for pixel_chunk in pixel_chunks_id1_id2
    )
    cooler.create_cooler(
        coolpout,
        new_bins,
        pixel_chunks_tidy
    )


### --- distiller-nf helper functions --- ###


def gen_config(
        sample_table, 
        cfg, 
        assembly="mm10",
        bwa="/share/home/ychi/data/genome/GRCm38/bwa_index/mm10.*",
        chrom_size="/share/home/ychi/data/genome/GRCm38/mm10.len.forCool.tsv",
        template="/shareb/ychi/ana/distiller-nf/project.yml",
        **args
        ):
    """
    Generate distiller.nf project.yaml file from bubble_flow sample_table.csv.
    Input:
        sample_table: any csv, index as sample_name, using sample-name and task_dirp to infer DNA reads files
        cfg: output project.yml file path
        assembly: genome assembly name, default mm10
        template: template project.yml file path
    """
    # read
    with open(template, "rt") as f:
        config = yaml.load(f, Loader=yaml.Loader)
    # editing input
    raw_reads_paths = {}
    sample_table = pd.read_csv(sample_table, index_col=0)
    for i in sample_table.index:
        raw_reads_paths[i] = {"lane1": [sample_table.loc[i, "R1_file"], sample_table.loc[i, "R2_file"]]}
    config["input"]["raw_reads_paths"] = raw_reads_paths
    config["input"]["library_groups"] = None
    # editing genome
    genome = {}
    genome["assembly_name"] = assembly
    genome["bwa_index_wildcard_path"] = bwa
    genome["chrom_sizes_path"] = chrom_size
    config["input"]["genome"] = genome
    # update other key:value
    if args is not None:
        for key in args:
            config[key].update(args[key])
    # write
    with open(cfg, "wt") as f:
        yaml.dump(config, f)
    print("Config generated.")
def stat_cool(names, ddir, genome, show=True):
    """
    Statistic technical details.
    Input:
        names: list of cooler paths
        ddir: distiller-nf directory
        genome: genome name eg. "mm10"
    Output:
        pd.DataFrame
    """
    pairs_stats_pat = str(Path(ddir)/"results"/"pairs_library"/"{sample}.{genome}.dedup.stats")
    cool_pat = str(Path(ddir)/"results"/"coolers_library"/"{sample}.{genome}.mapq_30.1000.cool")
    stats = pd.concat(
        list(map(
            lambda x : pd.read_table(                                                                                                                                                                                           
                pairs_stats_pat.format(sample=x, genome=genome),
                names=["item", x],
                index_col=0
            ),
            names
        )),
        axis = 1
    )
    annote = stats.iloc[:8,:]
    filtered_contacts = list(map(
        lambda x: Cooler(cool_pat.format(sample = x, genome=genome)).info["sum"],
        names
    ))
    annote = annote.T.assign(
        filtered_contacts = filtered_contacts
    )
    annote = annote.assign(
        contact_ratio = annote["filtered_contacts"] / annote["total"],
        trans_ratio = annote["trans"]/(annote["cis"] + annote["trans"])
    )
    annote = annote.astype(
        dtype={
            "total" : "int",
            "total_unmapped" : "int",
            "total_single_sided_mapped" : "int",
            "total_mapped" : "int",
            "total_dups" : "int",
            "total_nodups" : "int",
            "cis" : "int",
            "trans" : "int",
            "filtered_contacts" : "int",
            "contact_ratio" : "float",
            "trans_ratio" : "float"
        }
    )
    annote = annote.round({"contact_ratio":3, "trans_ratio":3})
    if show:
        A = annote.astype("string").copy()
        return A.T
    else:
        return annote
# lousy helper functions to work with cooler, cooltools stuff
# you need a conda env that names `cooler`, with cooler and cooltools installed
import subprocess
def cli_pairs2cool(filei,fileo,sizef,binsize):
    """
    Generate .cool file from 4DN .pairs file
    Input:
        filei: input pairs file
        fileo: output .cool file
        sizef: chrom size file, 2col(chrom:str, size:int) tsv
        binsize: binsize as you wish
    """
    subprocess.check_output(
        "conda run -n embryo cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 %s:%d %s %s " % (sizef,binsize,filei,fileo),
        shell=True)
def cli_mergecool(incools, outcool, force=False, conda_env=None, cwd=None):
    """
    Merge cool files with same indices to get a consensus heatmap.
    Input:
        incools: list of .cool file path
        outcool: output .cool file path
        conda_env: (optional) conda environment to run in
        cwd: (optional) working directory
    """
    if not force and Path(outcool).exists():
        print(f"File '{outcool}' already exists. Skipping execution.")
        return outcool
    if conda_env is None:
        conda_run = ""
    else:
        conda_run = f"conda run -n {conda_env}"

    # 将输入文件列表转换为字符串
    incools_str = " ".join(incools)
    
    # 构造命令行命令
    cmd = f"{conda_run} cooler merge {outcool} {incools_str}"

    try:
        subprocess.check_output(
            cmd,
            shell=True,
            cwd=cwd
        )
        return True
    except subprocess.CalledProcessError as e:
        print(e.output.decode())
        return False
def cli_downsample(coolp, output, count=100e6, cis_count=None, fraction=None, threads=8, force=False, conda_env=None, cwd=None):
    if not force and Path(output).exists():
        print(f"File '{output}' already exists. Skipping execution.")
        return output
    conda = f"conda run -n {conda_env}" if conda_env else ""
    cis_count = f"--cis-count {cis_count}" if cis_count else ""
    fraction = f"--frac {fraction}" if fraction else ""
    count = f"--count {count}" if count else ""
    cmd = f"{conda} cooler random-sample {count} {cis_count} {fraction} -p {threads} {coolp} {output}"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            cwd=cwd
        )
        return output
    except subprocess.CalledProcessError as e:
        print(e.output)
        return None
def cli_balance(coolp, threads=8, force=False, name="weight", cis_only=False, trans_only=False,
    min_nnz=None, min_count=None, mad_max=None, ignore_diags=None, conda_env=None, cwd=None):
    """
    Balance a cooler matrix.
    Input:
        coolp: cooler file path
        threads: number of threads
        force: whether to overwrite existing balanced matrix
        name: name of the balanced matrix
        cis_only: balance only cis contacts
        trans_only: balance only trans contacts
        min_nnz: minimum number of non-zero contacts
        min_count: minimum number of contacts
        mad_max: mad_max filter, lower values are more stringent
            NOTE: this dominates all other filtering options
        conda_env: conda environment to run in
        cwd: working directory
    Output:
        True if successful, False otherwise
    """
    assert not (cis_only and trans_only)
    if not force:
        clr = cooler.Cooler(str(coolp))
        if name in clr.bins().columns:
            print(f"Balanced matrix '{name}' already exists. Skipping execution.")
            return True

    if conda_env is None:
        conda_run = ""
    else:
        conda_run = f"conda run -n {conda_env}"
    force = "--force" if force else ""
    cis_only = "--cis-only" if cis_only else ""
    trans_only = "--trans-only" if trans_only else ""
    min_nnz = f"--min-nnz {min_nnz}" if min_nnz is not None else ""
    min_count = f"--min-count {min_count}" if min_count is not None else ""
    mad_max = f"--mad-max {mad_max}" if mad_max is not None else ""
    ignore_diags = f"--ignore-diags {ignore_diags}" if ignore_diags is not None else ""
    cmd = f"{conda_run} cooler balance {force} {cis_only} {trans_only} {min_nnz} {min_count} {mad_max} {ignore_diags} --name {name} --nproc {threads} {coolp}"
    
    try:
        # Use Popen to execute the command
        with subprocess.Popen(
            cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            cwd=cwd,
            preexec_fn=os.setsid  # Allows sending signal to the process group
        ) as process:
            try:
                # Wait for the process to complete
                stdout, stderr = process.communicate()
            except KeyboardInterrupt:
                # Send SIGINT to the process group to ensure all child processes are killed
                os.killpg(os.getpgid(process.pid), signal.SIGINT)
                process.kill()
                print("Process was terminated due to KeyboardInterrupt")
                raise
            if process.returncode != 0:
                print(stderr.decode())
                raise subprocess.CalledProcessError(process.returncode, cmd, output=stdout, stderr=stderr)
            return True
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode())
        return False
def cli_zoomify(coolp, output, resolutions=[20000,40000,100000,500000,1000000], force=False, threads=8, conda_env=None, cwd=None):
    conda_run = f"conda run -n {conda_env}" if conda_env else ""
    if not force and Path(output).exists():
        print(f"File '{output}' already exists. Skipping execution.")
        return output
    resolutions = ",".join(map(str, resolutions))
    cmd = f"{conda_run} cooler zoomify -n {threads} -r {resolutions} -o {output} {coolp}"
    subprocess.check_output(
        cmd,
        shell=True,
        cwd = cwd
    )
    return output
def cli_compartment(coolp, phasing_track, outprefix, view=None, conda_env=None, cwd=None, force=False):
    outprefix_path = Path(outprefix)

    if not force and outprefix_path.with_suffix(".cis.vecs.tsv").exists():
        print(f"Files '{outprefix}.cis.vecs.tsv' already exist. Skipping execution.")
        return f"{outprefix}.cis.vecs.tsv", f"{outprefix}.cis.lam.txt"

    conda_run = f"conda run -n {conda_env}" if conda_env else ""
    view_option = f"--view {view}" if view else ""
    cmd = f"{conda_run} cooltools eigs-cis {view_option} --phasing-track {phasing_track} -o {outprefix} {coolp}"

    subprocess.check_output(
        cmd,
        shell=True,
        cwd=cwd
    )
    return f"{outprefix}.cis.vecs.tsv", f"{outprefix}.cis.lam.txt"

def cli_saddle(coolp, eigv, expected, outprefix, view, conda_env=None, cwd=None, force=False):
    outprefix_path = Path(outprefix)

    if not force and outprefix_path.with_suffix(".saddledump.npz").exists():
        print(f"Files '{outprefix}.saddledump.npz' already exist. Skipping execution.")
        return f"{outprefix}.saddledump.npz", f"{outprefix}.png"

    conda_run = f"conda run -n {conda_env}" if conda_env else ""
    view_option = f"--view {view}" if view else ""
    cmd = f"{conda_run} cooltools saddle --qrange 0.02 0.98 --fig png -o {outprefix} {view_option} --strength {coolp} {eigv} {expected}"

    subprocess.check_output(
        cmd,
        shell=True,
        cwd=cwd
    )
    return f"{outprefix}.saddledump.npz", f"{outprefix}.png"
def cli_IS(coolp, output, windowsizes, balanced=True, append_raw_scores=True, threads=8, conda_env=None, cwd=None, force=False):
    if not force and Path(output).exists():
        print(f"File '{output}' already exists. Skipping execution.")
        return output
    if conda_env is None:
        conda_run = ""
    else:
        conda_run = f"conda run -n {conda_env}"
    windowsizes = " ".join(map(str, windowsizes))
    threads = f"-p {threads}"
    raw_score_option = "--append-raw-scores" if append_raw_scores else ""
    if balanced:
        cmd = f"{conda_run} cooltools insulation {threads} --threshold Li {raw_score_option} -o {output} {coolp} {windowsizes}"
    else:
        cmd = f'{conda_run} cooltools insulation {threads} --threshold Li {raw_score_option} --ignore-diags 1 --clr-weight-name "" -o {output} {coolp} {windowsizes}'
    subprocess.check_output(
        cmd,
        shell=True,
        cwd = cwd
    )
    return output
def cli_pileup(coolp, feature, output, format="BED", view=None, expected=None, flank=100000, conda_env=None, cwd=None, threads=8):
    if conda_env is None:
        conda_run = ""
    else:
        conda_run = f"conda run -n {conda_env}"
    if view is None:
        view = ""
    else:
        view = f"--view {view}"
    if expected is None:
        expected = ""
    else:
        expected = f"--expected {expected}"

    with tempfile.NamedTemporaryFile(mode="w+t", delete=True, dir=cwd if cwd is not None else "/tmp") as tmpfile:
        if isinstance(feature, str) or isinstance(feature, Path):
            # copy content from feature to tmpfile
            with open(feature, "rt") as f:
                tmpfile.write(f.read())
            tmpfile.seek(0)
            #print(tmpfile.read())
            feature = tmpfile.name
        elif isinstance(feature, pd.DataFrame):
            feature.to_csv(tmpfile, sep="\t", header=False, index=False)
            feature = tmpfile.name
        else:
            raise ValueError("feature should be a path or a dataframe")
        cmd = f"{conda_run} cooltools pileup -p {threads} {view} {expected} --flank {flank} --features-format {format} -o {output} {coolp} {feature}"
        subprocess.check_output(
            cmd,
            shell=True,
            cwd = cwd
        )
    return output
def cli_expected(coolp, output, balanced=False, view=None, ignore_diags=1, conda_env=None, cwd=None, threads=8, force=False):
    output_path = Path(output)

    if not force and output_path.exists():
        print(f"File '{output}' already exists. Skipping execution.")
        return output

    conda_run = f"conda run -n {conda_env}" if conda_env else ""
    view_option = f"--view {view}" if view else ""
    if ignore_diags is None:
        ignore_diags = f"--ignore-diags 1" if balanced else ""
    else:
        ignore_diags = f"--ignore-diags {ignore_diags}"
    if not balanced:
        cmd = f'{conda_run} cooltools expected-cis -p {threads} {ignore_diags} --clr-weight-name "" -o {output} {view_option} {coolp}'
    else:
        cmd = f"{conda_run} cooltools expected-cis -p {threads} {ignore_diags} -o {output} {view_option} {coolp}"
    print(cmd)
    subprocess.run(cmd, shell=True, cwd=cwd)
    return output
# def cli_ct_callTAD(filei,fileo):
#     """
#     cooltools(ct) call TAD
#     """
#     subprocess.check_output(
#         "conda run -n embryo cooltools diamond-insulation -o %s --append-raw-scores %s 100000" %(fileo,filei),
#         shell=True
#     )