from pathlib import Path

import pandas as pd

from ..hicio import read_meta, load_json
def _parse_bowtie2_log(log):
    """
    Parse the bowtie2 log file to get the overall alignment rate.
    """
    if not log.exists():
        return pd.NA
    with open(log) as f:
        for line in f:
            if "overall alignment rate" in line:
                return float(line.split()[0][:-1])
def _pp_json(filepat, samples, key_dict, dtypes = None):
    """
    Parse json files and return a dataframe.
    """
    data = {
        sample : load_json(str(filepat).format(sample = sample))
        for sample in samples
        }
    df = pd.DataFrame(data).T
    df = df.rename(
        columns = {
            value : key
            for key, value in key_dict.items()
        }
    )
    key_used = list(key_dict.keys())
    df = df[key_used].astype(dtypes)
    return df
def add_gam_mapping(sample_table_fp, ana_home)->pd.DataFrame:
    """
    Adding mapping and primary mapping rate to meta table.
    Input:
        sample_table_fp: path to the meta table.
        ana_home: analysis home directory, check ana_home/logs/mapping.log
            for the mapping rate.
    Output:
        sample_table: the original meta table with mapping rate added.
    """
    if isinstance(sample_table_fp, str):
        sample_table = read_meta(sample_table_fp)
    elif isinstance(sample_table_fp, pd.DataFrame):
        sample_table = sample_table_fp
    else:
        raise ValueError("sample_table_fp should be a path or a DataFrame")
    ddir = Path(ana_home)
    meta = sample_table.assign(
        mapping_rate = sample_table.apply(
            lambda row: _parse_bowtie2_log(
                ddir / "logs" / f"{row.name}.mapping.log"
            ),
            axis = 1
        )
    )
    dup_info = _pp_json(
        ddir / "dedup" / "{sample}.json",
        meta.index,
        dict(zip(
            ["read_num", "dup_read_num"],
            ["READ", "DUPLICATE TOTAL"]
            )),
        dtypes = {
            "read_num" : int,
            "dup_read_num" : int
        }
    )
    meta = pd.concat(
        [meta, dup_info],
        axis = 1
    )
    meta = meta.assign(
        unique_mapped = meta["read_num"] - meta["dup_read_num"],
    )
    return meta