from sklearn import preprocessing
from sklearn.decomposition import PCA

def pca_rep(mat, n_components, with_std=False):
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
