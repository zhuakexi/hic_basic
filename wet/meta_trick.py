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