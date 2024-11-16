from time import time

import pandas as pd

from .hicio import read_expr
#from sklearn import preprocessing
def _merge_expr(fps, mapper={}, samplelist=[], outfile=None, qc=None, T=False,dtype="int32"):
    """
    Merging umi count matrices.
    Input:
        fps: matrix file paths; list
        mapper: renaming samples(for ID fixing); dict
            blank for no change
        samplelist: samples to keep, [] for All; list
        T: True for transpose
        dtype: data type for matrix
            Warning: will do it silently
    Output:
        dataframe
    """
    all_mat = pd.DataFrame()
    T_run_start = time()
    for i, fp in enumerate(fps):
        #sub_mat = matr(i,"\t")
        sub_mat = read_expr(fp)
        sub_mat = sub_mat.astype(dtype)
        print("sub mats (%d/%d):" % (i+1, len(fps)) ,sub_mat.shape)
        # pd concat
        if T:
            # genes * samples -> samples * genes
            all_mat = pd.concat([all_mat, sub_mat.T],axis=0,join="outer")
            all_mat.fillna(0, inplace=True)
            # purge all zero columns(genes)
            all_mat = all_mat.loc[:,~(all_mat.sum(axis=0)==0)]
        else:
            all_mat = pd.concat([all_mat, sub_mat],axis=1,join="outer")
            all_mat.fillna(0, inplace=True)
            all_mat = all_mat.loc[~(all_mat.sum(axis=1)==0),:]
        print("all mats:",all_mat.shape)
    # ---rename sample id in expression matrix---
    if len(mapper) > 0: 
        if T:
            all_mat.rename(index=mapper,inplace=True)
        else:
            all_mat.rename(columns=mapper,inplace=True)

    # --- check and filtering samples ---
    if len(samplelist) > 0:
        samplelist = pd.Index(samplelist)
        # check all_mat
        if (~samplelist.isin(all_mat.columns)).sum() > 0:
            print("Warning: sample has no count data:")
            for i in samplelist[~samplelist.isin(all_mat.columns)]:
                print(i)
        print("raw merged matrix",all_mat.shape)
        # using only samples in samplelist
        valid_cells = samplelist[samplelist.isin(all_mat.columns)]
        all_mat = all_mat.loc[:,valid_cells]
        print("final merged matrix",all_mat.shape)
    T_run = (time() - T_run_start) / 60
    print("Took %.2f mins to merge" % T_run)
    if outfile is not None:
        Tio_start = time()
        all_mat.to_parquet(outfile)
        Tio = (time() - Tio_start) / 60
        print("IO took %.2f mins" % Tio)
    return all_mat
def standard_scaler(df, axis=0, with_std=False):
    """
    Return scaled dataframe.
    Input:
        df: pandas dataframe
        axis: 0 or 1. 0 to scale each col(col as feature), 1 to scale each row(row as feature).
        with_std: whether doing std scaling. No consensus but 
            usually don't scale it in omic-biology.
    Output:
        scaled dataframe. 
    """
    scaler = preprocessing.StandardScaler(with_std=with_std)
    if axis == 0:
        scaled = scaler.fit_transform(df.values)
        scaled = pd.DataFrame(scaled, index=df.index, columns=df.columns)
    elif axis == 1:
        scaled = scaler.fit_transform(df.values.T)
        scaled = pd.DataFrame(scaled.T, index=df.index, columns=df.columns)
    else:
        raise ValueError("axis should be 0 or 1.")
    return scaled