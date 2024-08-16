import itertools
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
# --- input --- #
def safe_input(file_pat, with_key=False, **kwargs):
    """
    Expand file_pat to see if the input files exist, if not, skip.
    Input:
        file_pat: file pattern
        with_key: return outkeys or not, bool
        kwargs: key word arguments, to expand file_pat with wildcards
    Output:
        list of input files or (list of input files, list of outkeys)
        filtered by the existence of the input files
    """
    # transform all values to str
    file_pat = str(file_pat)
    kwargs = {
        key : [str(value) for value in values]
        for key, values in kwargs.items()
    }
    # get all combinations of values
    key_combines = [dict(zip(kwargs.keys(), combination)) for combination in itertools.product(*kwargs.values())]
    # check file existence and filter
    inputs, outkeys = [], []
    for key_combine in key_combines:
        file = file_pat.format(**key_combine)
        if Path(file).exists():
            inputs.append(file)
            outkeys.append(key_combine)
        else:
            print(f"Waring: {file} not exists, skip {key_combine}")
            pass
    if with_key:
        return inputs, outkeys
    else:
        return inputs
def check_input(file_pat, **kwargs):
    # transform all values to str
    file_pat = str(file_pat)
    kwargs = {
        key : [str(value) for value in values]
        for key, values in kwargs.items()
    }
    key_combines = [dict(zip(kwargs.keys(), combination)) for combination in itertools.product(*kwargs.values())]
    good = True
    for key_combine in key_combines:
        file = file_pat.format(**key_combine)
        if Path(file).exists():
            pass
        else:
            print(f"Waring: {file} not exists")
            good = False
    return good
def check_key(file_pat, key_to_check, **kwargs):
    """
    Check one key against other keys in the file pattern.
    For exmaple, if key_to_check is "sample", pick sample
    that has all files that are expanted from file_pat using
    other keys in kwargs.
    Input:
        file_pat: file pattern
        key_to_check: key to check, one of the keys in kwargs, str
        kwargs: other keys, dict
    Output:
        list of values for key_to_check
    """
    # transform all values to str
    file_pat = str(file_pat)
    sub_dict = {
        key : [str(value) for value in values]
        for key, values in kwargs.items()
        if key != key_to_check
    }

    input_values = kwargs[key_to_check]
    good_values = []
    for input_value in input_values:
        #print(input_value)
        sub_dict.update({key_to_check: [input_value]})
        if check_input(file_pat, **sub_dict):
            good_values.append(input_value)
        else:
            print(f"Warning: {input_value} don't have all files")
    return good_values


### --- intermediate files --- ###
#TODO: add a function to generate the file pattern for the intermediate files from smk rules


def check_exsit(path):
    """
    Check if a file exists.
    Tolerate the pd.NA, None, math.nan, np.nan
    """
    if path is pd.NA:
        res = False
    elif path is None:
        res = False
    elif path is math.nan:
        res = False
    elif path is np.nan:
        res = False
    else:
        try:
            res = Path(path).exists()
        except TypeError:
            print(path)
            res = False
    return res
def bubble_flow_touched():
    """
    This function returns a dictionary that maps the file paths pattern for
    all the files that are touched by the bubble_flow pipeline.
    """
    DNA = [
        "{task_dirp}/DNA/{sample}_R1.fq.gz",
        "{task_dirp}/DNA/{sample}_R2.fq.gz"
    ]
    RNA = [
        "{task_dirp}/RNA/{sample}_R1.fq.gz",
        "{task_dirp}/RNA/{sample}_R2.fq.gz"

    ]
    hic_mapped = [
        "{task_dirp}/hic_mapped/{sample}.sorted.bam",
        "{task_dirp}/hic_mapped/{sample}.sorted.bam.bai"
    ]
    seg = [
        "{task_dirp}/seg/{sample}.seg.gz"
    ]
    pre_seg = [
        "{task_dirp}/seg/{sample}.seg.gz"
    ]
    pairs_0 = [
        "{task_dirp}/pairs_0/{sample}.pairs.gz",
        "{task_dirp}/pairs_0/{sample}.raw_pairs.gz"
    ]
    pairs_c1 = [
        "{task_dirp}/pairs_c1/{sample}.c1.pairs.gz",
    ]
    pairs_c12 = [
        "{task_dirp}/pairs_c12/{sample}.c12.pairs.gz",
    ]
    dip = [
        "{task_dirp}/dip/{sample}.dip.pairs.gz",
    ]
    info = [
        "{task_dirp}/info/{sample}.basic.info",
        "{task_dirp}/info/{sample}.reads.info",
        "{task_dirp}/info/{sample}.rna_reads.info",
        "{task_dirp}/info/{sample}.dna_reads.info",
        "{task_dirp}/info/{sample}.seg_stat.info",
        ]
    _3dg = []
    for seed in [1,2,3,4,5]:
        for reso in ["4m","1m","200k","50k","20k"]:
            _3dg.append(
                "{{task_dirp}}/3dg/{{sample}}.{reso}.{seed}.3dg".format(
                    seed = seed,
                    reso = reso
                )
            )
    _3dg_c = []
    for seed in [1,2,3,4,5]:
        for reso in ["4m","1m","200k","50k","20k"]:
            _3dg_c.append(
                "{{task_dirp}}/3dg_c/{{sample}}.clean.{reso}.{seed}.3dg".format(
                    seed = seed,
                    reso = reso
                )
            )
    mapper = {
        "DNA" : DNA,
        "RNA" : RNA,
        "hic_mapped" : hic_mapped,
        "seg" : seg,
        "pre_seg" : pre_seg,
        "info" : info,
        "_3dg" : _3dg,
        "_3dg_c" : _3dg_c,
        "pairs_0" : pairs_0,
        "pairs_c1" : pairs_c1,
        "pairs_c12" : pairs_c12,
        "dip" : dip
    }
    return mapper
def symlink_files(sample, task_dirp, target_dir, omit=None, filepats=None):
    """
    Create symbolic links from a old task directory to a new task directory.
    Use this to partially run the pipeline.
    Input:
        sample: the sample name
        task_dirp: the old task directory
        target_dir: the new task directory
        omit: the pipeline steps to omit
        filepats: key-value pairs of file patterns created by the pipeline,
            if None, use the file patterns created by the bubble_flow
            key: the pipeline step name
            value: the file pattern, must contain {sample} and {task_dirp} as placeholders
    Return:
        None
    """
    if filepats is not None:
        mapper = filepats
    else:
        mapper = bubble_flow_touched()
    if omit is None:
        subset = mapper.keys()
    else:
        subset = set(mapper.keys()) - set(omit)
    #print(subset)
    for key in subset:
        for pat in mapper[key]:
            from_ = pat.format(sample=sample, task_dirp=task_dirp)
            to_ = pat.format(sample=sample, task_dirp=target_dir)
            # only create the link when the original file exists
            if check_exsit(from_):
                # create the parent directory if not exists
                Path(to_).parent.mkdir(parents=True, exist_ok=True)
                if Path(to_).is_symlink():
                    if Path(to_).resolve() == Path(from_).resolve():
                        # do nothing when the link is already exist
                        # and the link is the same
                        pass
                    else:
                        # update the link when the link is already exist
                        # but link to wrong file
                        Path(to_).unlink()
                        Path(to_).symlink_to(from_)
                elif Path(to_).exists():
                    # this happens when the new task is already run
                    print(f"File {to_} already exists, skip")
                    pass
                else:
                    Path(to_).symlink_to(from_)
            else:
                pass
def unlink_files(sample, target_dir, subset=None):
    mapper = bubble_flow_touched()
    if subset is None:
        subset = mapper.keys()
    for key in subset:
        for pat in mapper[key]:
            to_ = Path(pat.format(sample=sample, task_dirp=target_dir))
            if to_.is_symlink():
                to_.unlink()
            else:
                pass


# --- helpers --- #
def json_reader(file):
    with open(file, "r") as f:
        data = json.load(f)
    return data