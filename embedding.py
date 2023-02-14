import time
import numpy as np
import pandas as pd
import h5py
from scipy.sparse import csr_matrix, vstack
import sklearn.neighbors as NN
from sklearn import preprocessing
from sklearn.decomposition import PCA, TruncatedSVD

from .metrics import pairwise_DKL

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

def band_svd(cell_list, res, dist=10000000, dim=50):
    """
    Doing svd on cell matrix band(leg distance < dist).
    TODO: using api compatible to cools.
    """
    celllist = np.loadtxt(cell_list, dtype=np.str)

    with h5py.File(celllist[0], 'r') as f:
            ngene = f['Matrix'].attrs['shape'][0]
    idx = np.triu_indices(ngene, k=1)
    idxfilter = np.array([(yy - xx) < (dist / res + 1) for xx, yy in zip(idx[0], idx[1])])
    idx = (idx[0][idxfilter], idx[1][idxfilter])

    start_time = time.time()
    # matrix = np.zeros((len(celllist), np.sum(idxfilter)))
    matrix = []
    for i, cell in enumerate(celllist):
        with h5py.File(cell, 'r') as f:
            g = f['Matrix']
            A = csr_matrix((g['data'][()], g['indices'][()],
                        g['indptr'][()]), g.attrs['shape'])
        # matrix[i] = A[idx]
        matrix.append(csr_matrix(A[idx]))
        if i%100 ==0:
            print(i, 'cells loaded', time.time() - start_time, 'seconds')

    matrix = vstack(matrix)

    scalefactor = 100000
    matrix.data = matrix.data * scalefactor
    svd = TruncatedSVD(n_components=dim, algorithm='arpack')
    matrix_reduce = svd.fit_transform(matrix)
    matrix_reduce = np.concatenate((svd.singular_values_[None, :], matrix_reduce), axis=0)
    return matrix_reduce
def band_seg_svd(cell_list, res=100e3, segL=10e6, dist=10e6, dim=50):
    """
    Doing svd on cell matrix band(leg distance < dist) segments.
    Input:
        cell_list: file list of cell matrix.
        res: resolution of matrix(binsize in bp).
        segL: length of segment(in bp).
        dist: distance of band, only pixels within the near-diagonal band are considered.
        dim: number of components to keep.
    Return:
        list of svd results, same length as nseg.
    """
    celllist = np.loadtxt(cell_list, dtype=np.str)
    assert segL % res == 0, 'segment length should be integer times of resolution'
    ngene = int(segL//res) # nbins for each segment
    dim = min(dim, ngene)
    with h5py.File(celllist[0], 'r') as f:
            tgene = int(f['Matrix'].attrs['shape'][0]) # total bin number
    nseg = int(tgene // ngene) # number of segments, omit the last segment
    idx = np.triu_indices(ngene, k=1)
    idxfilter = np.array([(yy - xx) < (dist / res + 1) for xx, yy in zip(idx[0], idx[1])])
    idx = (idx[0][idxfilter], idx[1][idxfilter])
    start_time = time.time()
    # loading all matrix
    matrix = []
    for i, cell in enumerate(celllist):
        with h5py.File(cell, 'r') as f:
            g = f['Matrix']
            A = csr_matrix((g['data'][()], g['indices'][()],
                        g['indptr'][()]), g.attrs['shape'])
        matrix.append(A)
        if i%100 ==0:
            print(i, 'cells loaded', time.time() - start_time, 'seconds')
            #pass
    # embedding for each segment
    embeddings = []
    for seg in range(nseg):
        print(f"doing segment {seg} ...")
        lidx = (idx[0] + seg*ngene, idx[1] + seg*ngene) # local index
        #print("Segnum:", seg, "nseg", nseg, "ngene", ngene, "tgene", tgene ,"idx:", idx, "shapeM", matrix[0].shape, sep="\n")
        mat = vstack(
            [csr_matrix(m[lidx]) for m in matrix]
            )
        scalefactor = 100000
        mat.data = mat.data * scalefactor
        svd = TruncatedSVD(n_components=dim, algorithm='arpack')
        mat_reduce = svd.fit_transform(mat)
        #matrix_reduce = np.concatenate((svd.singular_values_[None, :], matrix_reduce), axis=0)
        embeddings.append(mat_reduce)
        print(f"segment {seg} done in", time.time() - start_time, 'seconds')
    return embeddings
def schic_spectral_embedding(df):
    """
    Non-linear dimensionality reduction of cdps curve
    Input:
        df: sample x feature matrix
    Output:
        2nd 3rd smallest eigenvectors of the graph Laplacian
    """
    symdist = pairwise_DKL(df)

    # construct NN graph

    k_nn = 7
    nn = NN.NearestNeighbors(n_neighbors=k_nn, metric='precomputed', n_jobs=-1)
    nn.fit(symdist)
    dist2neighs, neighs =  nn.kneighbors()

    adj = np.zeros((symdist.shape[0], symdist.shape[0]))
    for i in range(symdist.shape[0]):
        for ki in range(k_nn):
            adj[i,neighs[i,ki]] = 1
    adj = np.maximum(adj, adj.T)

    # compute spectral embedding

    lap = np.diag(adj.sum(axis=1))-adj
    vals, vecs = np.linalg.eig(lap)
    tmp = np.argsort(vals)
    vals = vals[tmp]
    vecs = vecs[:,tmp]
    
    return pd.DataFrame(vecs[:,1:3],columns=["ev1","ev2"],index=df.index)