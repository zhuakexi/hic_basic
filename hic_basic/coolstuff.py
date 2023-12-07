# cooler api functions
import yaml
import os
import os.path as op
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from cytoolz import compose # fail in HPC
from cooler.create import sanitize_records, aggregate_records
import cooler

from .hicio import schicluster2mat
from .data import chromosomes
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
def _scool_pixels(pairs_paths:dict, bins, chunksize=15e6, zero_based=False, tril_action="reflect"):
    """
    Generate pixels dict for multi-pairs-file. Using in scool generating.
    Return:
        Cell name as key and pixel table DataFrame as value
    """
    n = len(pairs_paths)
    return {pairs_path : _easy_pixels(pairs_paths[pairs_path], bins, chunksize,zero_based, tril_action) 
        for pairs_path in pairs_paths}

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

# --- distiller-nf helper functions ---
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
def cli_mergecool(incools,outcool):
    """
    Merge cool files with same indices to get a consensus heatmap.
    Input:
        incools: list of .cool file path
        outcool: output .cool file path
    """
    subprocess.check_output(
        "conda run -n embryo cooler merge %s %s" %(outcool, incools),
        shell=True
    )
def cli_balance(coolp, threads=8, force=False, conda_env=None, cwd=None):
    if conda_env is None:
        conda_run = ""
    else:
        conda_run = f"conda run -n {conda_env}"
    if force:
        force = "--force"
    else:
        force = ""
    cmd = f"{conda_run} cooler balance {force} --nproc {threads} {coolp}"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            cwd = cwd
        )
        return True
    except subprocess.CalledProcessError as e:
        print(e.output)
        return False
def cli_IS(coolp, output, windowsizes, conda_env=None, cwd=None):
    if conda_env is None:
        conda_run = ""
    else:
        conda_run = f"conda run -n {conda_env}"
    windowsizes = " ".join(map(str, windowsizes))

    cmd = f"{conda_run} cooltools insulation --threshold Li -o {output} {coolp} {windowsizes}"
    subprocess.check_output(
        cmd,
        shell=True,
        cwd = cwd
    )
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
def cli_expected(coolp, output, conda_env=None, cwd=None, threads=8):
    if conda_env is None:
        conda_run = ""
    else:
        conda_run = f"conda run -n {conda_env}"
    subprocess.run(
        f"{conda_run} cooltools expected-cis -p {threads} -o {output} {coolp}",
        shell=True,
        cwd=cwd
    )
    return output
def cli_ct_callTAD(filei,fileo):
    """
    cooltools(ct) call TAD
    """
    subprocess.check_output(
        "conda run -n embryo cooltools diamond-insulation -o %s --append-raw-scores %s 100000" %(fileo,filei),
        shell=True
    )