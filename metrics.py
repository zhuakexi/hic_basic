# distance/similarity metrics
import numpy as np
from scipy.stats import entropy
from sklearn import preprocessing
from sklearn.metrics.pairwise import euclidean_distances, rbf_kernel
from sklearn.decomposition import PCA

def _pca_rep(mat, n_components, with_std=False): # keep this to avoid circular import in embedding.py
    """
    Return pca representation of data.
    Input:
        mat: m * n matrix
        n_components: number of components to keep
        with_std: whether doing std scaling. No consensus but 
            usually don't scale it in omic-biology.
    Output:
        m * n_components matrix
    """
    scaler = preprocessing.StandardScaler(with_std=with_std)
    scaled = scaler.fit_transform(mat)
    pca = PCA(n_components)
    pca_res = pca.fit_transform(scaled)
    return pca_res
def distant_pca_euclid(mat:np.ndarray, n_components:int, with_std=False) -> np.ndarray:
    """
    Euclid distances in PCA-space.
    Input:
        mat: m( sample ) * n( feature ) matrix
        n_components: number of PCA components used to calculate distance.
    Output:
        m * m distance matrix
    """
    dm = euclidean_distances(
        _pca_rep(
            mat, n_components, with_std
            )
        )
    return dm
def kernel_pca_euclid(mat:np.ndarray, n_components:int) -> np.ndarray:
    """    
    Affinity based on PCA-euclid-distance
    Input:
        mat: m( sample ) * n( feature ) matrix
        n_components: number of PCA components used to calculate distance.
    Output:
        m * m similarity matrix
    """
    dm = distant_pca_euclid(mat, n_components)
    sm = np.exp(-dm * 1/n_components)
    return sm
def kernel_c_rbf(mat):
    scaler = preprocessing.StandardScaler()
    scaled = scaler.fit_transform(mat)
    return rbf_kernel(scaled)
def pairwise_DKL(df):
    """
    compute pairwise Kullback-Leibler Divergence.
    resulting D'(i,j) = D'(j,i) = D_KL(i,j) + D_KL(j,i).
    Input:
        df: sample x feature matrix
    Output:
        matrix: symdist
    """
    y = df.values.astype(np.float) # omit non-digest and trans
    y[y==0] = 1 # add a pseudocount
    y = y/y.sum(axis=1)[:,np.newaxis]

    celln = y.shape[0]
    dists = np.zeros((celln, celln))
    for i in range(celln):
        for j in range(i+1, celln):
            dists[i,j] = entropy(y[i,:], qk=y[j,:], base=None)
            dists[j,i] = entropy(y[j,:], qk=y[i,:], base=None)

    symdist = dists+dists.T
    return symdist