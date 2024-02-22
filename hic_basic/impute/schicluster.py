from typing import Optional

import numpy as np
import pandas as pd
import scipy as sp
from cooler import Cooler
from scipy import signal, ndimage
from ..coolstuff import cool2mat
# --- schicluster --- #
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
# --- imputation main --- #
def schicluster_impute(
    coolp, region, fo,
    max_dist:Optional[int] = 20000000,
    alpha:float = 0.05, kernel_size:int = 3, sigma:int = 2, if_convolve:bool = True
):
    """
    Impute Hi-C matrix using schicluster method from cooler file.
    Input:
        coolp: cooler file
        region: region to impute, normally a chromosome,
            otherwise maybe incompatible with downstream analysis.
            See cool2mat's region argument
        fo: output file, str
        max_dist: max linear distance, int
        alpha: alpha, float
        kernel_size: kernel size, int
        sigma: sigma, float
        if_convolve: if convolve, bool
    Output:
        write imputed matrix to fo
    """
    coolp = str(coolp)
    mat = cool2mat(coolp, region, balance = False)
    imputed = schicluster_imputation_for_mat(
        mat.copy().values,
        alpha, kernel_size, sigma, if_convolve
        )
    result = pd.DataFrame(
        imputed,
        index = mat.index,
        columns = mat.columns
        )
    result = result.stack().reset_index()
    result.columns = ["start1", "start2", "imputed"]
    result = result.assign(chrom = region)
    if max_dist is not None:
        result = result.query("start2 - start1 <= @max_dist and start2 - start1 >= 0")
    result = result[["chrom", "start1", "start2", "imputed"]]
    result.to_parquet(fo)