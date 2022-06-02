from concurrent import futures
from itertools import dropwhile, repeat
import gzip
import json
import os
import time

import pandas as pd

from ..cycle_phasing import dis_counts
from ..repli_score import _repli_score, _make_repli_dict
# --- generate cooler files
def gen_cool(qc, cached, threads=16, sizef = "/share/home/ychi/data/genome/GRCm38/mm10.chrom.sizes", binsize=40_000, cool_addr="cool_paths.json"):
    """
    Generate cooler files from qc file.
    Using mid-files of compartment strength pipeline is better.
    """
    # --- check input ---
    valid_samples, fileis = check_input(qc, "pairs_c123")
    # --- prepare io ---
    if not os.path.isdir(os.path.join(cached, "cool", str(binsize))):
        os.makedirs(os.path.join(cached, "cool", str(binsize)))
    caching_files = [os.path.join(cached, "cool", str(binsize), "{name}.{res}.cool".format(name = i, res=binsize)) 
                     for i in valid_samples]
    # ---calc---
    # valid_samples, fileis, caching_files
    print("Calculating...")
    startt = time.time()
    with futures.ProcessPoolExecutor(threads) as pool:
        res = pool.map(
            pairs2cool,
            fileis,
            caching_files,
            repeat(sizef, len(fileis)),
            repeat(binsize, len(fileis)),
        )
    ares = list(res)
    print("Finished. Using %.2f min." % ((time.time() - startt)/60) )
    with open(cool_addr,"wt") as f:
        json.dump(dict(zip(valid_samples, caching_files)), f)
    return dict(zip(valid_samples, caching_files))
# --- generate contact decay profiles ---
def gen_cdps(filesp, threads = 32):
    # --- check pairs file ---
    valid_samples, valid_files = check_input(filesp, "pairs_c123")
    # --- calculate --- 
    with futures.ProcessPoolExecutor(threads) as pool:
        res = pool.map(
            dis_counts,
            valid_files
        )
    # --- collect data ---
    ares = list(res)
    cdps = pd.DataFrame(ares,index=valid_samples)
    cdps.columns = cdps.columns.astype("int")
    return cdps

# --- count Maternal-Paternal interaction ---
def count_MP(filename:str)->int:
    with gzip.open(filename,"rt") as f:
        inter = 0
        MP = 0
        MMPP = 0
        ALL = 0
        for line in dropwhile(lambda x:x.startswith("#"),f):
            attrs = line.strip().split()
            if attrs[1].split("(")[0] in ["chrX","chrY"] or attrs[3].split("(")[0] in ["chrX","chrY"]:
                continue
            ALL += 1
            if attrs[1] != attrs[3]:
                inter += 1
                if attrs[1].split("(")[1] == attrs[3].split("(")[1]:
                    MMPP +=1
                else:
                    MP +=1
    return ALL, inter, MP, MMPP
def gen_PM_interactions(filesp, threads=32):
    # ---check input---
    valid_samples, valid_files = check_input(filesp, "dip")
    # ---calculate---
    with futures.ProcessPoolExecutor(threads) as pool:
        res = pool.map(
            count_MP,
            valid_files
        )
    # ---collect info---
    ares = list(res)
    pm_interaction_autosome = pd.DataFrame(ares)
    pm_interaction_autosome.index = valid_samples
    pm_interaction_autosome.columns = ["All_contacts","All_inter_contacts","MP","MMPP"]
    return pm_interaction_autosome

# -- calculate repli score ---
def safe_repli_score(pairsf,ref,outf):
    if os.path.isfile(outf):
        # already has output
        return 1
    out = _repli_score(pairsf,ref)
    with open(outf,"wt") as f:
        json.dump(out,f)
    return 0
def check_input(filesp, col):
    valid_samples = []
    valid_files = []
    for sample, row in filesp.iterrows():
        file = row[col]
        if not os.path.isfile(file):
            print("Warning: input .pairs file missing.")
            print("Nofile" + file)
        else:
            valid_samples.append(sample)
            valid_files.append(file)
    return valid_samples, valid_files
def gen_repli_score(qc, cached, threads=32, ref="/share/home/ychi/data/genome/GRCm38/mm10_repli_chip.wig"):
    filesp = qc
    # ---check input---
    valid_samples, fileis = check_input(filesp, "pairs_c123")
    # --prepare io---
    repli_ref = _make_repli_dict(ref)
    caching_files = [os.path.join(cached, i+".rs.json") for i in valid_samples]
    # ---calc---
    # valid_samples, fileis, caching_files
    print("Calculating...")
    startt = time.time()
    with futures.ProcessPoolExecutor(threads) as pool:
        res = pool.map(
            safe_repli_score,
            fileis,
            repeat(repli_ref,len(fileis)),
            caching_files
        )
    # ---collect to single df---
    dats = []
    for file in caching_files:
        with open(file,"rt") as f:
            dats.append(json.load(f))
    rs = pd.DataFrame(dats)
    rs.index = valid_samples
    print("Finished. Using %.2f min." % ((time.time() - startt)/60) )
    return rs