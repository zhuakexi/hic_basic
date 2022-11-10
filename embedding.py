import time
import numpy as np
import h5py
from scipy.sparse import load_npz, save_npz, csr_matrix, vstack
from sklearn import preprocessing
from sklearn.decomposition import PCA, TruncatedSVD

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
def band_seg_svd(cell_list, res, nseg=10, dist=10000000, dim=50):
    """
    Doing svd on cell matrix band(leg distance < dist) segments.
    Return:
        list of svd results, same length as nseg.
    """
    celllist = np.loadtxt(cell_list, dtype=np.str)
    with h5py.File(celllist[0], 'r') as f:
            ngene = int(f['Matrix'].attrs['shape'][0] / nseg)
    embeddings = []
    for seg in range(nseg):
        print(f"doing segment {seg} ...")
        idx = np.triu_indices(ngene, k=1)
        idxfilter = np.array([(yy - xx) < (dist / res + 1) for xx, yy in zip(idx[0], idx[1])])
        idx = (idx[0][idxfilter] + seg*ngene, idx[1][idxfilter] + seg*ngene)

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
        #matrix_reduce = np.concatenate((svd.singular_values_[None, :], matrix_reduce), axis=0)
        embeddings.append(matrix_reduce)
    return embeddings