import re
import pandas as pd

def merge_meta(*dfs):
    """
    Merging metadata object.
    Input:
        DataFrame object seperate with , .
    """
    new_meta = pd.concat(dfs,axis=0,join="outer")
    for df in dfs:
        print("adding %d samples..." % df.shape[0])
    if new_meta.shape[1] > dfs[0].shape[1]:
        print("Warning: expanding cols.")
    else:
        print("Cols good.")
    print("Resulting in %d samples." % new_meta.shape[0])
    return new_meta
def rmsd_filter(meta, reso="20k", rmsd_thres=1, samples=True):
    rmsd_cols = meta.columns[
        meta.columns.str.contains(f"{reso}_\d+_\d+_rmsd", regex=True)
        ]
    dat = meta[rmsd_cols]
    # value of the third smallest groups
    working_rmsd = dat[dat.rank(
        axis=1
    ) == 3].max(axis=1)
    if samples:
        return meta[working_rmsd < rmsd_thres].index
    else: # give boolean mask
        return working_rmsd < rmsd_thres
def last_mouse(df):
    """
    Find last batch number.
    """
    mouse_re = re.compile(r'm([\d,a-z]+)_')
    mouse = []
    nonnum = []
    for group in df["group"].unique():
        mres = mouse_re.search(group)
        if mres is None:
            print(group)
            raise ValueError("Illegal batch code.")
        else:
            try:
                num = int(mres.group(1))
            except ValueError:
                nonnum.append(mres.group(1))
                continue
            mouse.append(num)
    print("Non-num :",nonnum)
    return sorted(mouse)