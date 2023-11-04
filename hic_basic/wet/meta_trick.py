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