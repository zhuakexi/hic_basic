import pandas as pd
import numpy as np
from bisect import bisect_left, bisect_right
import gzip
from itertools import dropwhile
def _make_repli_dict(repli_chipf)->dict:
    # Parse and make per-chromosome sorted repli_score reference
    # Input:
    #    repli_chipf: file path
    # Output:
    #    {"chr1":{
    #        "starts" : [...],
    #        "ends" : [...],
    #        "values" : [...]
    #          }
    #    ...
    #    }
    repli_ref = pd.read_table(
        repli_chipf,
        header=None,
        names=["chrom","start","end","rep_value"]
    )
    repli_dict = {}
    for chrom, dat in repli_ref.groupby("chrom"):
        dat = dat.sort_values("start")
        repli_dict[chrom] = {
            "starts":dat["start"].values,
            "ends":dat["end"].values,
            "values":dat["rep_value"].values
        }
    return repli_dict
def _expand_mean(site,expand,starts,ends,values):
    lm = bisect_left(starts,site-expand)
    rm = bisect_right(ends,site+expand)
    return values[lm:rm].mean() if values[lm:rm].size > 0 else np.nan
def _repli_score(pairsf,repli_dict,expand=10_000,c1=1,p1=2,c2=3,p2=4):
    # Calculate replication score of Hi-C pairs file
    # Input:
    #    pairsf: file path
    #    repli_dict: reference segment-timing annote
    #    expand: length of expansion from query leg position
    #    c1,p1,c2,p2, key col index of pairs file
    # Output:
    #    repli_score: raw_repli_score defined by Nagano2017 
    #    annote_ratio: ratio of legs with valid reference timing
    with gzip.open(pairsf,"rt") as f:
        pos_rs = []
        for line in dropwhile(lambda x:x.startswith("#"),f):
            eles = line.strip().split("\t")
            pos_rs.append(_expand_mean(
                int(eles[p1]),
                expand,
                repli_dict[eles[c1]]["starts"],
                repli_dict[eles[c1]]["ends"],
                repli_dict[eles[c1]]["values"]
            ))
            pos_rs.append(_expand_mean(
                int(eles[p2]),
                expand,
                repli_dict[eles[c2]]["starts"],
                repli_dict[eles[c2]]["ends"],
                repli_dict[eles[c2]]["values"]
            ))
    pos_rs = pd.Series(pos_rs)
    repli_score = pos_rs[pos_rs > 0].count() / pos_rs.shape[0]
    annote_ratio = 1 - pos_rs[pos_rs.isna()].shape[0] / pos_rs.shape[0]
    return {"repli_score":repli_score,"annote_ratio":annote_ratio}