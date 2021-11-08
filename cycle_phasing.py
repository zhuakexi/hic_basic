import pandas as pd
import os
import gzip
from concurrent import futures
from .io import divide_name, parse_pairs

def window_count(distances:pd.DataFrame, win_num)->pd.Series:
    # count distribution of distance array
    ## breaks from Nagano2017:
    breaks = [0] + [1000*2**(0.125*i) for i in range(0, win_num)] 
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
def dis_counts(cell_name:str):
    # work for 11 column table only
    # get cell's intra contact's distribution in Peter's window
    # using customized .pairs parser
    contacts = parse_pairs(cell_name)

    # get contact distance array
    intra = contacts.query("chr1 == chr2")
    distances = abs(intra["pos1"] - intra["pos2"])
    # count according to Peter's window range
    counts = window_count(distances, 150)
    counts.name = cell_name
    #return counts
    return counts 