import numpy as np
import scipy as sp
from scipy import signal, ndimage
# schicluster
def solve_rwr_inverse(stoch_matrix, alpha = 0.05):
    m = stoch_matrix*(1-alpha)
    m = m.transpose()
    #y = sp.sparse.spdiags([1] * m.shape[0], 0, m.shape[0], m.shape[0], format = "csc")
    y = np.eye(m.shape[0])
    A = y - m

    s = None
    #A = A.todense()
    #y = y.todense()
    s = sp.linalg.solve(A, y)

    s *= alpha
    s += s.transpose()
    
    if y is not None:
        del y
    if A is not None:
        del A
    if m is not None:
        del m
    return s

def schicluster_imputation_for_mat(mat,alpha=0.05,kernel_size=3,sigma=2,if_convolve=True):
    gauss_kernel_1d = signal.gaussian(kernel_size, std=sigma)
    gauss_kernel_2d = np.outer(gauss_kernel_1d, gauss_kernel_1d)

    if if_convolve:
        # add if since snapHi-C did not convolve the matrix
        mat = ndimage.convolve(mat, gauss_kernel_2d, mode='constant', cval=0.0)

    np.fill_diagonal(mat[1:,:-1], mat[1:,:-1] + 1)
    np.fill_diagonal(mat[:-1,1:], mat[:-1,1:] + 1)
    # mat to stochastic matrix
    mat = mat / np.nansum(mat, axis = 0)

    mat = solve_rwr_inverse(mat,alpha)
    return mat    

def normalize_matrix(matrix):
    """
    z-score normalization for band
    """
    from scipy.stats import zscore
    normalized_matrix = np.zeros_like(matrix)
    for i in range(-matrix.shape[0] + 1, matrix.shape[1]):
        band = matrix.diagonal(i)
        normalized_band = zscore(band)
        
        if i >= 0:
            np.fill_diagonal(normalized_matrix[i:], normalized_band)
        else:
            np.fill_diagonal(normalized_matrix[:, -i:], normalized_band)
    
    return normalized_matrix