import pandas as pd
import os
import gzip
from concurrent import futures
from basic import divide_name

def window_count(distances:pd.DataFrame, win_num)->pd.Series:
    # count distribution of distance array
    windows = pd.Series([1000*2**(0.125*i) for i in range(0, win_num)]) # Peter's window
    window_count = []
    for index, point in enumerate(windows):
        if index == 0:
            count = len(distances[distances < point])
        elif index == win_num - 1:
            count = len(distances[distances >= point])
        else:
            count = len(distances[(distances >= point) & (distances < windows[index + 1])])
        window_count.append(count)
    window_count = pd.Series(window_count)
    window_count.index = range(1,win_num+1)
    # normalized by all intra contacts
    return window_count/len(distances)
def window_count_ratio_normed(distances:pd.DataFrame, win_num)->pd.Series:
    # count distribution of distance array
    windows = pd.Series([1000*2**(0.125*i) for i in range(0, win_num)]) # Peter's window
    window_count = []
    for index, point in enumerate(windows):
        if index == 0:
            count = len(distances[distances < point])
        elif index == win_num - 1:
            count = len(distances[distances >= point])
        else:
            count = len(distances[(distances >= point) & (distances < windows[index + 1])])
        window_count.append(count)
    window_count = pd.Series(window_count)
    window_count.index = range(1,win_num+1)
    # normalized by all intra contacts
    return window_count/len(distances)
def peter_dis_counts(cell_name:str, parser:"func") -> pd.Series:    
    # get cell's intra contact's distribution in Peter's window
    # using customized .pairs parser
    contacts = parser(cell_name)

    # get contact distance array
    intra = contacts.query("chr1 == chr2")
    distances = abs(intra["pos1"] - intra["pos2"])
    # count according to Peter's window
    counts = window_count(distances, 150)
    counts.name = cell_name
    #return counts
    return counts.reindex(range(38,151)) # only show 38-150
# walk around for jupyter's bug on 2nd layer function
def hap_dis_counts(cell_name:str):
    # work for 9 column table only
    # get cell's intra contact's distribution in Peter's window
    # using customized .pairs parser
    contacts = pd.read_table(cell_name, header=None, comment="#")
    contacts.columns = "readID chr1 pos1 chr2 pos2 strand1 strand2".split()

    # get contact distance array
    intra = contacts.query("chr1 == chr2")
    distances = abs(intra["pos1"] - intra["pos2"])
    # count according to Peter's window
    counts = window_count(distances, 150)
    counts.name = cell_name
    #counts.reindex(range(38,151))
    return counts.reindex(range(38,151)) # only show 38-150
def pairs_dis_counts(cell_name:str):
    # work for 11 column table only
    # get cell's intra contact's distribution in Peter's window
    # using customized .pairs parser
    contacts = pd.read_table(cell_name, header=None, comment="#")
    contacts.columns = "readID chr1 pos1 chr2 pos2 strand1 strand2 phase0 phase1".split()

    # get contact distance array
    intra = contacts.query("chr1 == chr2")
    distances = abs(intra["pos1"] - intra["pos2"])
    # count according to Peter's window
    counts = window_count(distances, 150)
    counts.name = cell_name
    #return counts
    return counts.reindex(range(38,151)) # only show 38-150

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
def collect(meta:pd.DataFrame, name_col:str, file_col:str, func:"callable", threads)->pd.DataFrame:
    # get cons for cells in dir
    file_names = meta[file_col]
    with futures.ProcessPoolExecutor(threads) as pool:
        res = pool.map(func, file_names)
    result = pd.concat(res, axis=1)
    result.columns = meta[name_col]
    return result
def parse_pairs(filename:str)->"Cell":
    '''
    read from 4DN's standard .pairs format
    compatible with all hickit originated pairs-like format 
    '''
    #comment lines are stored in dataframe.attrs["comment"]
    name_array = "readID chr1 pos1 chr2 pos2 strand1 strand2 phase0 phase1 phase_prob00 phase_prob01 phase_prob10 phase_prob11".split()
    #read comment line
    with gzip.open(filename,"rt") as f:
        comments = []
        for line in f.readlines():
            if line[0] != "#":
                break
            comments.append(line)
    #read table format data
    pairs = pd.read_table(filename, header=None, comment="#",low_memory=False)
    pairs.attrs["comments"] = comments
    pairs.attrs["name"], _ = divide_name(filename) # get real sample name
    #assign column names
    pairs.columns = name_array[0:pairs.shape[1]]
    #sys.stderr.write("pairs_parser: %s parsed \n" % filename)
    return pairs
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
    return counts.reindex(range(38,151)) # only show 38-150