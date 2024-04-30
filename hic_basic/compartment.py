import numpy as np
import pandas as pd
def distance_law_from_mat(matrix, indices=None, log_bins=True, base=1.1):
    """Compute distance law as a function of the genomic coordinate aka P(s).
    Bin length increases exponentially with distance if log_bins is True. Works
    on dense and sparse matrices. Less precise than the one from the pairs.
    Parameters
    ----------
    matrix : numpy.array or scipy.sparse.coo_matrix
        Hi-C contact map of the chromosome on which the distance law is
        calculated.
    indices : None or numpy array
        List of indices on which to compute the distance law. For example
        compartments or expressed genes.
    log_bins : bool
        Whether the distance law should be computed on exponentially larger
        bins.
    Returns
    -------
    numpy array of floats :
        The start index of each bin.
    numpy array of floats :
        The distance law computed per bin on the diagonal
    """

    n = min(matrix.shape)
    included_bins = np.zeros(n, dtype=bool)
    if indices is None:
        included_bins[:] = True
    else:
        included_bins[indices] = True
    D = np.array(
        [
            np.nanmean(matrix.diagonal(j)[included_bins[: n - j]])
            for j in range(n)
        ]
    )
    if not log_bins:
        return np.array(range(len(D))), D
    else:
        n_bins = int(np.log(n) / np.log(base) + 1)
        logbin = np.unique(
            np.logspace(0, n_bins - 1, num=n_bins, base=base, dtype=np.int)
        )
        logbin = np.insert(logbin, 0, 0)
        logbin[-1] = min(n, logbin[-1])
        if n < logbin.shape[0]:
            print("Not enough bins. Increase logarithm base.")
            return np.array(range(len(D))), D
        logD = np.array(
            [
                np.average(D[logbin[i - 1] : logbin[i]])
                for i in range(1, len(logbin))
            ]
        )
        return logbin[:-1], logD
def extrude_full_zero(mat:pd.DataFrame)->pd.DataFrame:
    # purge all-zero lines and columns
    # using DataFrame to keep row-col index
    
    ## boolean indexing
    return mat.loc[ mat.sum(axis=1) != 0, 
                  mat.sum(axis=0) != 0] # using pandas api
def upper_to_symmetry(mat:pd.DataFrame)->pd.DataFrame:
    # filling lower triangle of the upper tirangle matrix
    # input:
    # upper tirangle matrix in DataFrame
    # output:
    # symmetry matrix in DataFrame
    to_fill = np.triu(mat.values, 1).T # using numpy api
    return mat + to_fill
def normalize_dense(M, norm="SCN", order=1, iterations=40):
    """Apply one of the many normalization types to input dense
    matrix. Will also apply any callable norms such as a user-made
    or a lambda function.
    NOTE: Legacy function for dense maps

    Parameters
    ----------
    M : 2D numpy array of floats
    norm : str
        Normalization procedure to use. Can be one of "SCN",
        "mirnylib", "frag" or "global". Can also be a user-
        defined function.
    order : int
        Defines the type of vector norm to use. See numpy.linalg.norm
        for details.
    iterations : int
        Iterations parameter when using an iterative normalization
        procedure.

    Returns
    -------
    2D numpy array of floats :
        Normalized dense matrix.
    """

    s = np.array(M, np.float64)
    floatorder = np.float64(order)

    if norm == "SCN":
        for _ in range(0, iterations):

            sumrows = s.sum(axis=1)
            maskrows = (sumrows != 0)[:, None] * (sumrows != 0)[None, :]
            sums_row = sumrows[:, None] * np.ones(sumrows.shape)[None, :]
            s[maskrows] = 1.0 * s[maskrows] / sums_row[maskrows]

            sumcols = s.sum(axis=0)
            maskcols = (sumcols != 0)[:, None] * (sumcols != 0)[None, :]
            sums_col = sumcols[None, :] * np.ones(sumcols.shape)[:, None]
            s[maskcols] = 1.0 * s[maskcols] / sums_col[maskcols]

    elif norm == "mirnylib":
        try:
            from mirnylib import numutils as ntls

            s = ntls.iterativeCorrection(s, iterations)[0]
        except ImportError as e:
            print(str(e))
            print("I can't find mirnylib.")
            print(
                "Please install it from "
                "https://bitbucket.org/mirnylab/mirnylib"
            )
            print("I will use default norm as fallback.")
            return normalize_dense(M, order=order, iterations=iterations)

    elif norm == "frag":
        for _ in range(1, iterations):
            s_norm_x = np.linalg.norm(s, ord=floatorder, axis=0)
            s_norm_y = np.linalg.norm(s, ord=floatorder, axis=1)
            s_norm = np.tensordot(s_norm_x, s_norm_y, axes=0)
            s[s_norm != 0] = 1.0 * s[s_norm != 0] / s_norm[s_norm != 0]

    elif callable(norm):
        s = norm(M)

    else:
        raise Exception(
            'Unknown norm, please specify one of ("mirnylib", "SCN", "frag")'
        )

    return (s + s.T) / 2
from scipy.linalg import eig 
def compartments(M, normalize=True, matrixonly=False):
    """A/B compartment analysis

    Perform a PCA-based A/B compartment analysis on a normalized, single
    chromosome contact map. The results are two vectors whose values (negative
    or positive) should presumably correlate with the presence of 'active'
    vs. 'inert' chromatin.

    Parameters
    ----------
    M : array_like
        The input, normalized contact map. Must be a single chromosome.
    normalize : bool
        Whether to normalize the matrix beforehand.
    matrixonly : bool
        Whether to return only the matrix used for the PCA (useful in compartment plotting).

    Returns
    -------
    PC1 : numpy.ndarray
        A vector representing the first component.
    PC2 : numpy.ndarray
        A vector representing the second component.
    """

    n = M.shape[0]
    if type(M) is not np.ndarray:
        M = np.array(M)

    if M.shape[0] != M.shape[1]:
        raise ValueError("Matrix is not square.")

    if normalize:
        N = normalize_dense(M)
    else:
        N = np.copy(M)
    # Computation of genomic distance law matrice:

    dist_mat = np.zeros((n, n))
    _, dist_vals = distance_law_from_mat(N, log_bins=False)
    for i in range(n):
        for j in range(n):
            dist_mat[i, j] = dist_vals[abs(j - i)]
    N /= dist_mat
    #N = N/dist_mat
    # Computation of the correlation matrice:"""
    N = pd.DataFrame(N).corr().values # use pandas to treat NaNs
    N[np.isnan(N)] = 0.0

    if matrixonly:
        return N

    # Computation of eigen vectors:
    (eig_val, eig_vec) = eig(N)
    PC1 = eig_vec[:, 0]
    PC2 = eig_vec[:, 1]
    return PC1, PC2
def call_compartment(df:pd.DataFrame,norm:bool=True)->tuple:
    # read in contact map, retur PC1 and PC2
    n0_df = extrude_full_zero(df)
    return compartments(
        upper_to_symmetry(
            n0_df
        ),
        normalize=norm
        )
def AB_block_ends(orig_data:pd.DataFrame, min_winSize=3, binsize=1000000):
    """
    Generate AB bed-like reference from cooltools vecs dataframe.
    Input:
        orig_data: chrom, start, end, E1
    Output:
        chrom, start, end, AB, mean_E1; one row per block
    """
    data = orig_data.copy()
    data["sign"] = np.sign(data["E1"])
    data["block"] = (data["sign"].diff(1) != 0).astype('int').cumsum()
    data["block_size"] = data.groupby("block")["block"].transform("size")
    data_valid = data[data["block_size"] >= min_winSize].copy()
    data_valid['c_start'] = data_valid.groupby('block')['start'].transform(lambda x : x.min())
    data_valid["c_end"] = data_valid.groupby('block')["end"].transform(lambda x : x.max())
    data_valid["mean_E1"] = data_valid.groupby('block')["E1"].transform(lambda x : x.mean())
    ends = data_valid.drop_duplicates(subset=["block", "c_start", "c_end"])[["chrom","mean_E1","c_start","c_end"]]
    ends = ends.assign(AB=pd.NA)
    ends.loc[ends["mean_E1"]>0, "AB"] = "A"
    ends.loc[ends["mean_E1"]<0, "AB"] = "B"
    ends = ends.rename(columns = {"c_start":"start","c_end":"end"})
    ends = ends.dropna(subset=["AB"])[["chrom","start","end","AB","mean_E1"]]
    return ends