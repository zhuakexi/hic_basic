import concurrent
import sys
from .mdso.spectral_ordering_ import SpectralBaseline, SpectralOrdering
from .cycle_phasing import dis_counts
import numpy as np
import pandas as pd
from sklearn import preprocessing
from sklearn.decomposition import PCA
from sklearn.metrics.pairwise import euclidean_distances, rbf_kernel
from concurrent import futures
def euclid_kernel(cdps,n_components):
    # affinity based on PCA-euclid-distance
    # Input:
    #    n_components: number of PCA components used to calculate distance
    # Output:
    #    similarity matrix
    scaler = preprocessing.StandardScaler()
    scaled = scaler.fit_transform(cdps)
    pca = PCA()
    pca_res = pca.fit_transform(scaled)
    dm = euclidean_distances(pca_res[:,:n_components])
    sm = np.exp(-dm * 1/n_components)
    return sm
def euclid(cdps,n_components):
    scaler = preprocessing.StandardScaler()
    scaled = scaler.fit_transform(cdps)
    pca = PCA()
    pca_res = pca.fit_transform(scaled)
    dm = euclidean_distances(pca_res[:,:n_components])
    return dm
def c_rbf_kernel(cdps,n_components):
    scaler = preprocessing.StandardScaler()
    scaled = scaler.fit_transform(cdps)
    return rbf_kernel(scaled)
def _mdso(cdps:pd.DataFrame, annote:pd.DataFrame, n_components:int=6):
    """
    Input:
        cdps: contact decay profile, index must be sample name
        annote: annotaions, index must be sample name
    Output:
        list of sample_name as order
    """
    # calculate similarity matrix
    sm =pd.DataFrame(
        euclid_kernel(cdps.values,n_components),
        index = cdps.index,
        columns = cdps.index
        )
    so = SpectralOrdering()
    so_res = so.fit_transform(sm.values)
    return list(annote.iloc[so_res].index)
def spectral_ordering(filesp, threads=24):
    """
    Input:
        filesp : DataFrame, must have pairs col
    Output:
        DataFrame with additional col "order_index"
    """
    with futures.ProcessPoolExecutor(threads) as pool:
        res = pool.map(dis_counts,filesp["pairs"])
    ares = list(res)
    cdps = pd.DataFrame(ares)
    cdps.columns = cdps.columns.astype("string")
    cdps.index = filesp.index
    order = _mdso(cdps,filesp)
    order_index = pd.Series(list(range(len(order))),index=order,name="order_index")
    new_filesp = pd.concat([filesp,order_index],axis=1)
    return new_filesp