import numpy as np
import pandas as pd
import anndata as ad

from .metrics import euclidean_distances
def _gaussian_interpolate(expr, traj, winSz=0.1, numPts=200):
    """
    Interpolation of the expression data along one trajectory.
    Migrate from cellAlign package.
    Input:
        expr: raw single cell gene expression data (genes on rows, cells on columns); pd.DataFrame
        traj: a vector of pseudo-time scores of the data-points whose length equals to the number of samples; pd.Series
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
    # add interpolated points' name
    # for ValNewPts, gene names as row names, Ip_N as col names
    # for trajValNewPts, Ip_N as index
    interpolated_point_names = ["Ip_" + str(i) for i in range(numPts)]
    ValNewPts = pd.DataFrame(ValNewPts, index=expr.index, columns=interpolated_point_names)
    trajValNewPts = pd.Series(trajValNewPts, index=interpolated_point_names)
    return ValNewPts, trajValNewPts
    #return ValNewPts, error, trajValNewPts
def _correct_pseudotime(ps, dm):
    """
    Correct pseudotime to make it linearly correlated with a distance metrix.
    Inspired by Alpert et al. (2018)
    Input:
        ps: pseudotime, pd.Series
        dm: distance matrix
    Output:
        ps_cr: corrected pseudotime, pd.Series
    """
    ps = ps.sort_values()
    i1 = ps.index[0] # i1 is start point
    new_ps = {}
    for i in ps.index:
        earlier = ps[ps < ps[i]]
        latter = ps[ps > ps[i]]
        earlier_metric = dm.loc[earlier.index, i1] + dm.loc[earlier.index, i] # D_i_i1 = D_x_i1 + D_x_i
        latter_metric = dm.loc[latter.index, i1] - dm.loc[latter.index, i] # D_i_i1 = D_x_i1 - D_x_i
        full_metric = pd.concat([earlier_metric, latter_metric]).values # measure point i n-1 times
        new_ps[i] = full_metric.mean() # using mean
    return pd.Series(new_ps)
def gaussian_interpolate(adata, pseudotime_col, obs_using, winSz=0.1, numPts=200, inplace=True):
    """
    Interpolation of all concerned traits in adata.
    Input:
        adata: AnnData object, must have adata.obs[pseudotime_col]
        pseudotime_col: pseudotime column name
        obs_using: list of obs column names to be interpolated
        winSz: window size of the interpolation
        numPts: number of desired interpolated points
    Output:
        new adata object numPts * adata.shape[1]
    """
    # get pseudotime
    ps = adata.obs[pseudotime_col]
    # interpolate expression
    new_expr, new_ps = _gaussian_interpolate(adata.to_df().T, ps, winSz, numPts)
    cdps, _ = _gaussian_interpolate(adata.uns["cdps"].T, ps, winSz, numPts)
    new_obs, _ = _gaussian_interpolate(adata.obs[obs_using].T, ps, winSz, numPts)
    new_ps.name = pseudotime_col + "_ip"
    new_obs = pd.concat([new_obs.T, new_ps], axis=1)
    new_adata = ad.AnnData(
        new_expr.T.loc[new_obs.index],
        obs=new_obs,
        var=adata.var,
        uns={"cdps": cdps.T})
    # add new columns to adata
    return new_adata
def correct_pseudotime(adata, pseudotime_col, n_pcs=15, inplace=True):
    """
    Correct pseudotime of adata to make it linearly correlated with a distance metrix.
    Inspired by Alpert et al. (2018)
    Input:
        adata: AnnData object, must have adata.obs[pseudotime_col]
        pseudotime_col: pseudotime column name
    Output:
        new adata object
    """
    if not inplace:
        adata = adata.copy()
    # deal with missing values
    if adata.obs[pseudotime_col].isna().any():
        sub = adata[~adata.obs[pseudotime_col].isna()]
    else:
        sub = adata
    #ps, expr = sub.obs[pseudotime_col], sub.to_df().T
    #dm = pd.DataFrame(euclidean_distances(expr)/expr.shape[0]**0.5, index=expr.index, columns=expr.index)
    # using pca-euclidean distances
    dm = pd.DataFrame(
        euclidean_distances(sub.obsm["X_pca"][:, :n_pcs]),
        index = sub.obs_names,
        columns = sub.obs_names
    )
    new_ps = _correct_pseudotime(sub.obs[pseudotime_col], dm)
    adata.obs[pseudotime_col + "_cr"] = new_ps # add new column inplace, adata.obs = xx will not work
    if not inplace:
        return adata