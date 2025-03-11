import pandas as pd
import os
import gzip
from concurrent import futures
from .hicio import divide_name, parse_pairs

def window_count(distances:pd.Series, win_num)->pd.Series:
    """
    Count distribution of distance array with Nagano2017's window.
    If distances is None, return nan Series.
    Input:
        distances: distance array
        win_num: number of windows
    Output:
        window_count: distribution of distances in Nagano
    """
    breaks = [0] + [1000*2**(0.125*i) for i in range(0, win_num)]
    if distances is None:
        return pd.Series([float('nan')] * (win_num + 1), index=breaks)
    window_count = []
    for index, break_ in enumerate(breaks):
        if index == win_num:
            count = len(distances[distances >= break_])
        else:
            count = len(distances[(distances >= break_) & (distances < breaks[index + 1])])
        window_count.append(count)
    window_count = pd.Series(window_count, 
        index = breaks)
    # normalized by all intra contacts
    return window_count/len(distances)
def dis_counts(pairs_fp:str):
    """
    Get cell's intra contact's distribution in Nagano2017's window
    Using customized .pairs parser. Work for 11 column table only
    """
    if not os.path.exists(str(pairs_fp)):
        counts = window_count(None, 150)
    else:
        contacts = parse_pairs(pairs_fp)

        # get contact distance array
        intra = contacts.loc[contacts["chr1"] == contacts["chr2"]]
        distances = abs(intra["pos1"] - intra["pos2"])
        # count according to Peter's window range
        counts = window_count(distances, 150)
    counts.name = pairs_fp
    return counts