# distance metrics
from string import capwords
import numpy as np
from sklearn import preprocessing
from sklearn.decomposition import PCA
from sklearn.metrics.pairwise import euclidean_distances, rbf_kernel

def pca_euclid(mat:np.ndarray, n_components:int) -> np.ndarray:
    """
    Euclid distances in PCA-space.
    Input:
        mat: m( sample ) * n( feature ) matrix
        n_components: number of PCA components used to calculate distance.
    Output:
        m * m distance matrix
    """
    scaler = preprocessing.StandardScaler()
    scaled = scaler.fit_transform(mat)
    pca = PCA()
    pca_res = pca.fit_transform(scaled)
    dm = euclidean_distances(pca_res[:,:n_components])
    return dm
def pca_euclid_kernel(mat:np.ndarray, n_components:int) -> np.ndarray:
    """    
    Affinity based on PCA-euclid-distance
    Input:
        mat: m( sample ) * n( feature ) matrix
        n_components: number of PCA components used to calculate distance.
    Output:
        m * m similarity matrix
    """
    dm = pca_euclid(mat, n_components)
    sm = np.exp(-dm * 1/n_components)
    return sm
def c_rbf_kernel(mat):
    scaler = preprocessing.StandardScaler()
    scaled = scaler.fit_transform(mat)
    return rbf_kernel(scaled)