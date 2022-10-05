# cooler api functions
import yaml
import os
import os.path as op

import numpy as np
import pandas as pd
from cytoolz import compose
from cooler.create import sanitize_records, aggregate_records
import cooler

from .utils import read_chromsizes, binnify

def gen_bins(chromsizes_file:str, binsize:int):
    """
    Generate "bins" for cooler api.
    """
    if not op.exists(chromsizes_file):
        raise ValueError('File "{}" not found'.format(chromsizes_file))
    try:
        binsize = int(binsize)
    except ValueError:
        raise ValueError(
            'Expected integer binsize argument (bp), got "{}"'.format(binsize)
        )
    chromsizes = read_chromsizes(chromsizes_file, all_names=True)
    bins = binnify(chromsizes, binsize)
    return bins

def easy_pixels(pairs_path, bins, chunksize=15e6, zero_based=False, tril_action="reflect"):
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

def scool_pixels(pairs_paths:dict, bins, chunksize=15e6, zero_based=False, tril_action="reflect"):
    """
    Generate pixels dict for multi-pairs-file. Using in scool generating.
    Return:
        Cell name as key and pixel table DataFrame as value
    """
    n = len(pairs_paths)
    return {pairs_path : easy_pixels(pairs_paths[pairs_path], bins, chunksize,zero_based, tril_action) 
        for pairs_path in pairs_paths}

def pairs2cool(pairs_path, coolpath, sizef, binsize):
    """
    Generate .cool file from 4DN .pairs file
    """
    bins = gen_bins(sizef, binsize)
    cooler.create_cooler(
        coolpath,
        bins,
        easy_pixels(pairs_path, bins)
    )
def pairs2scool(pairs_paths, coolpath, sizef, binsize):
    """
    Generate .scool file from multiple 4DN .pairs file
    """
    bins = gen_bins(sizef, binsize)
    cooler.create_scool(
        coolpath,
        bins,
        scool_pixels(
            pairs_paths,
            bins
        )
    )
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
    Input:
        scool: cooler's scool file path
        grouping: pd.Series
        cache_pattern: path to store resulting .cool file; str with {}
    Output:
        write cool according to cache_pattern
    """
    for cluster in grouping.unique():
        samples = grouping[grouping==cluster].index
        print("Merging " + str(cluster))
        cooler.merge_coolers(
            cache_pattern.format(str(cluster)),
            [scool+"::/cells/{}".format(i) for i in samples],
            mergebuf
        )

# --- distiller-nf helper functions ---
def gen_config(sample_table, output, task_dirp, template="/shareb/ychi/repo/ara_sperm_bulk1/HIC/project.yml"):
    """
    Generate distiller.nf project.yaml file from bubble_flow sample_table.csv.
    Input:
        sample_table: any csv, index as sample_name, using sample-name and task_dirp to infer DNA reads files
        output: project.yml file
        task_dirp: top folder for real input files
    """
    # read
    with open(template, "rt") as f:
        config = yaml.load(f, Loader=yaml.Loader)
    # editing    
    raw_reads_paths = {}
    sample_table = pd.read_csv(sample_table, index_col=0)
    for i in sample_table.index:
        raw_reads_paths[i] = {"lane1": [os.path.join(task_dirp, "DNA", i+"_R1.fq.gz"),os.path.join(task_dirp, "DNA", i+"_R2.fq.gz")]}
    config["input"]["raw_reads_paths"] = raw_reads_paths
    # write
    with open(output, "wt") as f:
        yaml.dump(config, f)
    print("Config generated.")

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
        "conda run -n cooler cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 %s:%d %s %s " % (sizef,binsize,filei,fileo),
        shell=True)
def cli_mergecool(incools,outcool):
    """
    Merge cool files with same indices to get a consensus heatmap.
    Input:
        incools: list of .cool file path
        outcool: output .cool file path
    """
    subprocess.check_output(
        "conda run -n cooler cooler merge %s %s" %(outcool, incools),
        shell=True
    )
def cli_cl_Balance(filei,threads=8):
    """
    cooler(cl) balance
    """
    subprocess.check_output(
        "conda run -n cooler cooler balance -p %d --force %s" % (threads,filei),
        shell=True
    )
def cli_ct_callTAD(filei,fileo):
    """
    cooltools(ct) call TAD
    """
    subprocess.check_output(
        "conda run -n cooler cooltools diamond-insulation -o %s --append-raw-scores %s 100000" %(fileo,filei),
        shell=True
    )