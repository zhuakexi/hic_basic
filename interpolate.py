import numpy as np

def _gaussian_interpolate(expr, traj, winSz=0.1, numPts=200):
    """
    Interpolation of the expression data along one trajectory.
    Migrate from cellAlign package.
    Input:
        expr: raw single cell gene expression data (genes on rows, cells on columns)
        traj: a vector of pseudo-time scores of the data-points whose length equals to the number of samples
        winSz: window size of the interpolation
        numPts: number of desired interpolated points
    Output:
        ValNewPts: interpolated expression data
        trajValNewPts: trajectory values at the new points
    TODO: adding interpolated error
    """
    def gaussian_weights(x, y):
        return np.exp(-((x - y) ** 2) / winSz ** 2)
    # remove NAs from the trajectory (if exist)
    # and make sure they have same sample index
    if np.isnan(traj).any():
        traj = traj[~np.isnan(traj)]
        expr = expr[:, ~np.isnan(traj)]
    else:
        expr = expr.loc[:, traj.index]
    #generate interpolated values:
    trajValNewPts = np.linspace(min(traj.values), max(traj.values), numPts)
    weights = gaussian_weights(traj.values.reshape(traj.shape[0], 1), trajValNewPts)
    weights = weights/weights.sum(axis=0)
    ValNewPts = np.dot(expr, weights)
    return ValNewPts, trajValNewPts
    #return ValNewPts, error, trajValNewPts