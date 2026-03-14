import numpy as np
import pandas as pd
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, RBF
import anndata as ad
from scipy import sparse

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
def _get_expression_frame(adata, layer="X"):
    """
    Extract expression matrix from AnnData as a dense dataframe.
    """
    if layer == "X":
        matrix = adata.X
    else:
        if layer not in adata.layers:
            raise KeyError(f"Layer '{layer}' not found in adata.layers.")
        matrix = adata.layers[layer]
    if sparse.issparse(matrix):
        matrix = matrix.toarray()
    else:
        matrix = np.asarray(matrix)
    return pd.DataFrame(matrix, index=adata.obs_names, columns=adata.var_names)
def _prepare_gp_inputs(adata, pseudotime_col, layer="X", gene_mask_key=None):
    """
    Validate and extract GP interpolation inputs from AnnData.
    """
    if pseudotime_col not in adata.obs:
        raise KeyError(f"Column '{pseudotime_col}' not found in adata.obs.")

    try:
        pseudotime = pd.to_numeric(adata.obs[pseudotime_col], errors="raise")
    except Exception as exc:
        raise ValueError(f"Column '{pseudotime_col}' must be numeric.") from exc

    valid_cells = pseudotime.notna()
    num_cells_dropped = int((~valid_cells).sum())
    pseudotime = pseudotime.loc[valid_cells]
    if pseudotime.shape[0] < 2:
        raise ValueError("Need at least two cells with non-missing pseudotime values.")

    pt_min = pseudotime.min()
    pt_max = pseudotime.max()
    if pt_min == pt_max:
        raise ValueError("Pseudotime must contain at least two unique values.")
    pseudotime = (pseudotime - pt_min) / (pt_max - pt_min)
    if pseudotime.nunique() < 2:
        raise ValueError("Need at least two unique normalized pseudotime values.")

    expr = _get_expression_frame(adata, layer=layer).loc[pseudotime.index]

    if gene_mask_key is None:
        selected_genes = adata.var_names
    else:
        if gene_mask_key not in adata.var:
            raise KeyError(f"Column '{gene_mask_key}' not found in adata.var.")
        gene_mask = adata.var[gene_mask_key]
        if gene_mask.isna().any():
            raise ValueError(f"Column '{gene_mask_key}' contains missing values.")
        if not pd.api.types.is_bool_dtype(gene_mask):
            raise ValueError(f"Column '{gene_mask_key}' must be boolean.")
        selected_genes = adata.var_names[gene_mask.to_numpy()]

    if len(selected_genes) == 0:
        raise ValueError("No genes selected for GP interpolation.")

    expr = expr.loc[:, selected_genes]
    if not np.isfinite(expr.to_numpy()).all():
        raise ValueError("Expression matrix contains non-finite values.")

    return expr, pseudotime, num_cells_dropped
def _fit_gp_gene(y_train, x_train, x_grid, alpha=1e-6, length_scale=0.1,
                 length_scale_bounds=(1e-2, 1.0), n_restarts_optimizer=5,
                 normalize_y=True):
    """
    Fit one gene with a Gaussian Process and return predictions and metadata.
    """
    kernel = ConstantKernel(1.0, constant_value_bounds=(1e-3, 1e3)) * RBF(
        length_scale=length_scale,
        length_scale_bounds=length_scale_bounds
    )
    gp = GaussianProcessRegressor(
        kernel=kernel,
        alpha=alpha,
        normalize_y=normalize_y,
        n_restarts_optimizer=n_restarts_optimizer
    )
    gp.fit(x_train, y_train)
    mean_pred, std_pred = gp.predict(x_grid, return_std=True)
    fitted_kernel = gp.kernel_
    hyperparameters = {
        "length_scale": fitted_kernel.k2.length_scale,
        "signal_variance": fitted_kernel.k1.constant_value,
        "alpha": alpha,
        "log_marginal_likelihood": gp.log_marginal_likelihood_value_,
        "kernel_repr": str(fitted_kernel),
    }
    return mean_pred, std_pred ** 2, hyperparameters
def gp_interpolate(adata, pseudotime_col, layer="X", gene_mask_key=None,
                   numPts=14, alpha=1e-6, length_scale=0.1,
                   length_scale_bounds=(1e-2, 1.0), n_restarts_optimizer=5,
                   normalize_y=True, result_key="gp_interpolate", inplace=True):
    """
    Interpolate single-trajectory expression dynamics with Gaussian Processes.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object containing expression values and pseudotime annotations.
    pseudotime_col : str
        Column in ``adata.obs`` containing numeric pseudotime values.
    layer : str, optional
        Expression source. Use ``"X"`` for ``adata.X`` or the name of an
        entry in ``adata.layers``.
    gene_mask_key : str, optional
        Boolean column in ``adata.var`` used to select genes for fitting.
    numPts : int, optional
        Number of evenly spaced interpolation points on normalized pseudotime.
    alpha : float, optional
        Fixed observation noise variance passed to the GP regressor.
    length_scale : float, optional
        Initial RBF length scale.
    length_scale_bounds : tuple, optional
        Bounds for RBF length-scale optimization.
    n_restarts_optimizer : int, optional
        Number of optimizer restarts for hyperparameter fitting.
    normalize_y : bool, optional
        Whether to normalize the target values before fitting.
    result_key : str, optional
        Key used to store outputs in ``adata.uns``.
    inplace : bool, optional
        If ``True``, write results to ``adata`` directly. Otherwise fit on a
        copy and return the copy.

    Returns
    -------
    anndata.AnnData
        AnnData with GP interpolation results stored in ``adata.uns[result_key]``.
    """
    if numPts < 2:
        raise ValueError("numPts must be at least 2.")
    if alpha <= 0:
        raise ValueError("alpha must be positive.")
    if length_scale <= 0:
        raise ValueError("length_scale must be positive.")
    if n_restarts_optimizer < 0:
        raise ValueError("n_restarts_optimizer must be non-negative.")
    if len(length_scale_bounds) != 2 or length_scale_bounds[0] <= 0 or length_scale_bounds[0] >= length_scale_bounds[1]:
        raise ValueError("length_scale_bounds must be a positive (min, max) tuple with min < max.")

    target = adata if inplace else adata.copy()
    expr, pseudotime, num_cells_dropped = _prepare_gp_inputs(
        target,
        pseudotime_col=pseudotime_col,
        layer=layer,
        gene_mask_key=gene_mask_key
    )

    grid_values = np.linspace(0.0, 1.0, numPts)
    grid_index = [f"Gp_{i}" for i in range(numPts)]
    x_train = pseudotime.to_numpy(dtype=float).reshape(-1, 1)
    x_grid = grid_values.reshape(-1, 1)

    mean_rows = []
    variance_rows = []
    hyper_rows = []
    errors = []

    for gene_name in expr.columns:
        try:
            mean_pred, var_pred, hyperparameters = _fit_gp_gene(
                expr[gene_name].to_numpy(dtype=float),
                x_train,
                x_grid,
                alpha=alpha,
                length_scale=length_scale,
                length_scale_bounds=length_scale_bounds,
                n_restarts_optimizer=n_restarts_optimizer,
                normalize_y=normalize_y
            )
        except Exception as exc:
            errors.append(f"{gene_name}: {exc}")
            continue
        mean_rows.append(pd.Series(mean_pred, index=grid_index, name=gene_name))
        variance_rows.append(pd.Series(var_pred, index=grid_index, name=gene_name))
        hyper_rows.append(pd.Series(hyperparameters, name=gene_name))

    if errors:
        shown_errors = "; ".join(errors[:5])
        if len(errors) > 5:
            shown_errors += f"; ... ({len(errors) - 5} more)"
        raise RuntimeError(
            f"GP interpolation failed for {len(errors)} gene(s): {shown_errors}"
        )

    target.uns[result_key] = {
        "time_grid": pd.Series(grid_values, index=grid_index, name="normalized_pseudotime"),
        "mean": pd.DataFrame(mean_rows),
        "variance": pd.DataFrame(variance_rows),
        "hyperparameters": pd.DataFrame(hyper_rows),
        "metadata": {
            "pseudotime_col": pseudotime_col,
            "layer": layer,
            "gene_mask_key": gene_mask_key,
            "num_cells_used": int(pseudotime.shape[0]),
            "num_cells_dropped": num_cells_dropped,
            "num_genes_fit": int(expr.shape[1]),
            "numPts": int(numPts),
            "normalized_pseudotime": True,
            "input_preprocessed_assumed": True,
        },
    }
    return target
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
