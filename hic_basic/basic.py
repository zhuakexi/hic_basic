import pandas as pd
import gzip
import os
from concurrent import futures

# All kinds of simple calculations on adata
def add_contact_describe(adata, inplace=True):
    """
    Get cell's contact decay statistics, defined in Nagano2017.
    Input:
        adata: AnnData object, must have uns["cdps"]
    Output:
        adata will extra obs cols
    """
    bounds = adata.uns["cdps"].columns.astype(int)
    short = adata.uns["cdps"].loc[:, (23_000 < bounds) & (bounds < 2_000_000)].sum(axis=1)
    mitotic = adata.uns["cdps"].loc[:, (2_000_000 < bounds) & (bounds < 12_000_000)].sum(axis=1)
    data = adata.obs.assign(short_r = short)
    data = data.assign(mitotic_r = mitotic)
    if inplace:
        adata.obs = data
        return None
    else:
        return data
def add_pmUMIs(adata, inplace=True):
    """
    Add cell's g1 UMI ratio.
    Input:
        adata: AnnData object, must have ["g1_UMIs","g2_UMIs"]
    Output:
        adata will extra obs cols
    """
    data = adata.obs.assign(pmUMIs = adata.obs["g1_UMIs"]/(adata.obs["g1_UMIs"] + adata.obs["g2_UMIs"]))
    if inplace:
        adata.obs = data
        return None
    else:
        return data