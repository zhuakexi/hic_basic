import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
import anndata as ad
def _kmeans(mat:np.ndarray, names:list=[], n_clusters:int=12, random_state:int=0, n_components:int=6) -> pd.Series:
    """
    Input:
        mat: m * n matrix
        names: name of samples; m str list
    Return:
        "kmeans_cluster" Series with names as index.
    """
    mat_pca = pca_rep(mat, n_components)
    kmeans_clusters = KMeans(n_clusters = n_clusters, random_state = random_state).fit_predict(mat_pca)
    if len(names) != mat.shape[0]:
        print("Wrong sample names.")
        return kmeans_clusters
    else:
        return pd.Series(kmeans_clusters, index=names, name="kmeans_clusters",dtype="category")
def do_kmeans(adata:ad.AnnData, layer:str, **k) -> ad.AnnData:
    """
    Input:
        adata
        layer: which layer to operate kmeans on
    Return:
        adata with new obs col.
    """
    pass