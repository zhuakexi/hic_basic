import pandas as pd
import os
import gzip
from concurrent import futures
from hires_utils.hires_utils.hires_io import divide_name, parse_pairs

def G1_attrs_pairs(cell_name:str) -> pd.Series:
    # get cell's %near and farAvgDist
    
    ## get cell's intra contact's distance array 
    contacts = pd.read_table(cell_name, header=None, comment="#")
    contacts.columns = "readID chr1 pos1 chr2 pos2 strand1 strand2 phas0 phase1".split()
    intra = contacts.query("chr1 == chr2")
    distances = abs(intra["pos1"] - intra["pos2"])
    counts = window_count(distances, 150)
    
    all_ = counts.reindex(range(38,151)).sum() # all VALID bins, 1-37 are discarded, bin 38-150
    near = counts.reindex(range(38,90)).sum() # bin 38-89
    mitotic = counts.reindex(range(90,130)).sum() # bin 90-109
    #far = counts.reindex(range(98, 151)).sum() # bin >= 98
    
    #farAvgDist = distances[distances < 4096000.0].mean() #mean contact distance considering bins >= 98
    #result = pd.Series({"near_p":near/all_*100, "farAvgDist":farAvgDist})
    result = pd.Series({"near_p":near/all_*100, "mitotic":mitotic/all_*100})
    result.name = cell_name
    
    return result
def contact_describe(cell_name:str) -> pd.Series:
    # get cell's basic statistics, defined in Nagano2017
    contacts = parse_pairs(cell_name)
    intra = contacts.query(' chr1 == chr2 ')
    distances = abs(intra["pos1"] - intra["pos2"])
    
    all_ = len(distances[23_000 < distances])
    short = len(distances[(23_000 < distances) & (distances < 2_000_000)])
    mitotic = len(distances[(2_000_000 < distances) & (distances < 12_000_000)])
    farAvg = distances[(4_500_000 < distances) & (distances < 225_000_000)]
    
    mitotic_r = mitotic/all_
    short_r = short/all_
    
    # assign to different stages on Peter's cirtera
    if mitotic_r >= 0.3 and short_r <= 0.5:
        group = "Post-M"
    elif short_r > 0.5 and short_r + 1.8*mitotic_r > 1.0:
        group = "Pre-M"
    elif short_r <= 0.63:
        group = "G1"
    elif 0.63 < short_r <= 0.785:
        group = "early-S"
    elif short_r > 0.785:
        group = "late-S/G2"
    else:
        group = "blank"
    
    return pd.Series({"short%":short_r, "mitotic%":mitotic_r, "farAvg":farAvg.mean(),"group":group })

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
def Nagano2017_cdp(cell_name:str):
    # work for 11 column table only
    # get cell's intra contact's distribution in Peter's window
    # using customized .pairs parser
    contacts = parse_pairs(cell_name)

    # get contact distance array
    intra = contacts.query("chr1 == chr2")
    distances = abs(intra["pos1"] - intra["pos2"])
    # count according to Peter's window
    counts = window_count(distances, 150)
    counts.name = cell_name
    #return counts
    return counts.reindex(range(38,151)) # only show 38-150
def dis_counts(cell_name:str):
    # work for 11 column table only
    # get cell's intra contact's distribution in Peter's window
    # using customized .pairs parser
    contacts = parse_pairs(cell_name)

    # get contact distance array
    intra = contacts.query("chr1 == chr2")
    distances = abs(intra["pos1"] - intra["pos2"])
    # count according to Peter's window
    counts = window_count(distances, 150)
    counts.name = cell_name
    #return counts
    return counts # only show 38-150