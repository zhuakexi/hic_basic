import numpy as np
import pandas as pd
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, RBF
import anndata as ad
from scipy import sparse
from scipy.stats import gaussian_kde

from .metrics import euclidean_distances


def _plot_grid_comparison(
    adata,
    pseudotime_col,
    layer="X",
    gene_mask_key=None,
    numPts=14,
    grid_methods=None,
    endpoint_mode="half-interval",
    scaling="z_score",
    figsize=(10, 12),
    dpi=100,
):
    """
    Visualize and compare multiple grid methods side-by-side.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object with expression data and pseudotime.
    pseudotime_col : str
        Column name in adata.obs containing pseudotime values.
    layer : str, optional
        Expression layer to use. Default is "X".
    gene_mask_key : str, optional
        Boolean column in adata.var for gene selection.
    numPts : int, optional
        Number of grid points. Default is 14.
    grid_methods : list, optional
        List of grid methods to compare. Default is ["linspace", "quantile", "density"].
    endpoint_mode : str, optional
        Endpoint handling mode. Default is "half-interval".
    scaling : str, optional
        Color scaling: "z_score" or "zero_one". Default is "z_score".
    figsize : tuple, optional
        Figure size (width, height) in inches. Default is (10, 12).
    dpi : int, optional
        Figure DPI. Default is 100.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The matplotlib figure object.
    axes : np.ndarray
        Array of matplotlib axes objects.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    if grid_methods is None:
        grid_methods = ["linspace", "quantile", "density"]

    # Get expression and pseudotime
    # Note: AnnData is cells x genes, we need genes x cells for plotting
    if layer == "X":
        expr = adata.X.T  # Transpose to genes x cells
    else:
        expr = adata.layers[layer].T
    if sparse.issparse(expr):
        expr = expr.toarray()

    pseudotime = adata.obs[pseudotime_col].values
    cell_order = np.argsort(pseudotime)
    expr_sorted = expr[:, cell_order]  # genes x sorted_cells
    pseudotime_sorted = pseudotime[cell_order]

    # Apply gene mask if provided
    if gene_mask_key is not None:
        gene_mask = adata.var[gene_mask_key].values
        expr_sorted = expr_sorted[gene_mask, :]

    # Normalize expression for visualization
    if scaling == "z_score":
        expr_mean = expr_sorted.mean(axis=1, keepdims=True)
        expr_std = expr_sorted.std(axis=1, keepdims=True) + 1e-6
        expr_norm = (expr_sorted - expr_mean) / expr_std
        cmap = "RdBu_r"
        vmin, vmax = -3, 3
    else:  # zero_one
        expr_min = expr_sorted.min(axis=1, keepdims=True)
        expr_max = expr_sorted.max(axis=1, keepdims=True)
        expr_norm = (expr_sorted - expr_min) / (expr_max - expr_min + 1e-6)
        cmap = "inferno"
        vmin, vmax = 0, 1

    n_methods = len(grid_methods)
    n_genes = expr_norm.shape[0]
    n_cells = expr_norm.shape[1]

    # Create figure
    fig, axes = plt.subplots(
        n_methods, 2,
        figsize=figsize,
        dpi=dpi,
        sharey=True,
        gridspec_kw={"width_ratios": [3, 2], "hspace": 0.15, "wspace": 0.1}
    )

    if n_methods == 1:
        axes = axes.reshape(1, -1)

    cell_ranks = np.arange(n_cells) / (n_cells - 1)

    for row_idx, method in enumerate(grid_methods):
        # Generate grid
        grid, _ = _build_time_grid(
            pd.Series(pseudotime),
            numPts=numPts,
            method=method,
            endpoint_mode=endpoint_mode
        )

        # Left panel: Raw data with grid lines
        ax = axes[row_idx, 0]
        extent = [0, 1, 0, n_genes]
        im = ax.imshow(
            expr_norm,
            aspect="auto",
            origin="lower",
            extent=extent,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )

        # Add grid lines at uniform positions based on grid point index
        # For visualization, grid points should be evenly spaced in rank space
        n_grid = len(grid)
        for i in range(n_grid):
            # Uniform spacing in rank space: i / (n_grid - 1)
            rank_pos = i / (n_grid - 1) if n_grid > 1 else 0.5
            ax.axvline(x=rank_pos, color="black", linestyle="--", linewidth=0.8, alpha=0.7)

        ax.set_title(f"{method}: Raw (grid lines)", fontsize=10, fontweight="bold")
        ax.set_xlabel("Cell rank (sorted by pseudotime)")
        if row_idx == 0:
            ax.set_ylabel("Gene index")

        # Right panel: GP-smoothed (schematic - show grid points)
        ax = axes[row_idx, 1]
        grid_indices = np.arange(len(grid))

        # Create a schematic representation of smoothed data
        # (In practice, you would run gp_interpolate first)
        smooth_schematic = np.random.randn(n_genes, len(grid)) * 0.5
        for i in range(n_genes):
            # Create smooth transition based on raw data
            raw_gene = expr_sorted[i, :]
            smooth_schematic[i, :] = np.interp(
                grid,
                cell_ranks,
                np.convolve(raw_gene, np.ones(5) / 5, mode="same")
            )

        # Normalize smooth data
        if scaling == "z_score":
            smooth_mean = smooth_schematic.mean(axis=1, keepdims=True)
            smooth_std = smooth_schematic.std(axis=1, keepdims=True) + 1e-6
            smooth_norm = (smooth_schematic - smooth_mean) / smooth_std
        else:
            smooth_min = smooth_schematic.min(axis=1, keepdims=True)
            smooth_max = smooth_schematic.max(axis=1, keepdims=True)
            smooth_norm = (smooth_schematic - smooth_min) / (smooth_max - smooth_min + 1e-6)

        extent = [-0.5, len(grid) - 0.5, 0, n_genes]
        im = ax.imshow(
            smooth_norm,
            aspect="auto",
            origin="lower",
            extent=extent,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )

        ax.set_title(f"{method}: GP-smoothed ({len(grid)} pts)", fontsize=10, fontweight="bold")
        ax.set_xlabel("Grid point index")
        ax.set_xticks(grid_indices)
        ax.set_xticklabels(grid_indices.astype(str))

        if row_idx == 0:
            ax.set_ylabel("Gene index")

    # Add colorbar
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax) if scaling == "z_score" else None
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, cax=cbar_ax, label="Z-score" if scaling == "z_score" else "Scaled [0, 1]")

    fig.suptitle(
        f"GP Interpolation: Grid Method Comparison (endpoint_mode='{endpoint_mode}')",
        y=0.995,
        fontsize=12,
        fontweight="bold"
    )

    return fig, axes


def plot_grid_methods_comparison(
    adata,
    pseudotime_col,
    layer="X",
    gene_mask_key=None,
    numPts=14,
    grid_methods=None,
    endpoint_mode="half-interval",
    scaling="z_score",
    figsize=(10, 12),
    dpi=100,
    save_path=None,
):
    """
    Plot comparison of different grid methods for GP interpolation.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object with expression data and pseudotime.
    pseudotime_col : str
        Column name in adata.obs containing pseudotime values.
    layer : str, optional
        Expression layer to use. Default is "X".
    gene_mask_key : str, optional
        Boolean column in adata.var for gene selection.
    numPts : int, optional
        Number of grid points. Default is 14.
    grid_methods : list, optional
        List of grid methods to compare. Default is ["linspace", "quantile", "density"].
    endpoint_mode : str, optional
        Endpoint handling mode. Default is "half-interval".
    scaling : str, optional
        Color scaling: "z_score" or "zero_one". Default is "z_score".
    figsize : tuple, optional
        Figure size (width, height) in inches. Default is (10, 12).
    dpi : int, optional
        Figure DPI. Default is 100.
    save_path : str, optional
        Path to save the figure. If None, figure is not saved.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The matplotlib figure object.
    """
    fig, _ = _plot_grid_comparison(
        adata,
        pseudotime_col,
        layer=layer,
        gene_mask_key=gene_mask_key,
        numPts=numPts,
        grid_methods=grid_methods,
        endpoint_mode=endpoint_mode,
        scaling=scaling,
        figsize=figsize,
        dpi=dpi,
    )

    if save_path is not None:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")

    return fig


def _build_time_grid(pseudotime, numPts=14, method="linspace", grid_custom=None,
                     endpoint_mode="half-interval"):
    """
    Build interpolation time grid using various strategies.

    Parameters
    ----------
    pseudotime : pd.Series
        Normalized pseudotime values (range [0, 1]).
    numPts : int, optional
        Number of interpolation points. Default is 14.
    method : str, optional
        Grid construction method. Options:
        - "linspace": evenly spaced points (default)
        - "quantile": quantile-based spacing (equal cells per bin)
        - "density": density-adaptive spacing (more points in dense regions)
        - "optimal": use optimal binning to determine numPts automatically
        - "custom": use user-provided grid (requires grid_custom)
    grid_custom : array-like, optional
        Custom grid values (only used when method="custom").
    endpoint_mode : str, optional
        How to handle endpoints. Options:
        - "half-interval" (default): endpoints at [half_interval, 1-half_interval]
          to avoid sparse cells at t=0 and t=1
        - "full": endpoints at [0, 1] (legacy behavior)

    Returns
    -------
    grid_values : np.ndarray
        Interpolation time grid.
    grid_index : list of str
        Index labels for the grid points.
    """
    if method == "custom":
        if grid_custom is None:
            raise ValueError("grid_custom must be provided when method='custom'.")
        grid_values = np.asarray(grid_custom)
        if not np.all((grid_values >= 0) & (grid_values <= 1)):
            raise ValueError("Custom grid values must be in range [0, 1].")
        if not np.all(np.diff(grid_values) > 0):
            raise ValueError("Custom grid values must be strictly increasing.")
        numPts = len(grid_values)
        grid_index = [f"Gp_{i}" for i in range(numPts)]
        # Warn about endpoint sparsity
        if grid_values[0] < 0.05 or grid_values[-1] > 0.95:
            import warnings
            warnings.warn(
                f"Custom grid endpoints ({grid_values[0]:.3f}, {grid_values[-1]:.3f}) "
                f"may have sparse cells. Consider using endpoint_mode='half-interval' "
                f"or adjusting grid to avoid t=0 and t=1."
            )
        return grid_values, grid_index

    if numPts < 2:
        raise ValueError("numPts must be at least 2.")

    # Calculate half-interval for endpoint truncation
    interval = 1.0 / (numPts - 1)
    half_interval = interval / 2
    
    # Determine endpoint range based on endpoint_mode
    if endpoint_mode == "half-interval":
        t_min = half_interval
        t_max = 1.0 - half_interval
    elif endpoint_mode == "full":
        t_min = 0.0
        t_max = 1.0
    else:
        raise ValueError(
            f"Unknown endpoint_mode='{endpoint_mode}'. "
            f"Options: 'half-interval', 'full'."
        )

    if method == "linspace":
        grid_values = np.linspace(t_min, t_max, numPts)

    elif method == "quantile":
        # Find quantiles corresponding to t_min and t_max
        pt_sorted = np.sort(pseudotime.values)
        n = len(pt_sorted)
        # Convert t_min/t_max to quantile indices
        q_min_idx = int(np.searchsorted(pt_sorted, t_min * pt_sorted.max()))
        q_max_idx = int(np.searchsorted(pt_sorted, t_max * pt_sorted.max()))
        q_min = max(0.0, q_min_idx / (n - 1))
        q_max = min(1.0, q_max_idx / (n - 1))
        quantiles = np.linspace(q_min, q_max, numPts)
        grid_values = np.quantile(pseudotime, quantiles)

    elif method == "density":
        # Density-adaptive grid: more points in dense regions
        # Use KDE to estimate density, then sample based on CDF
        pt_sorted = np.sort(pseudotime.values)
        kde = gaussian_kde(pt_sorted)
        # Evaluate density on fine grid
        x_fine = np.linspace(0, 1, 1000)
        density = kde(x_fine)
        # Compute CDF of density
        cdf = np.cumsum(density)
        cdf = cdf / cdf[-1]
        # Find CDF values corresponding to t_min and t_max
        cdf_min_idx = np.searchsorted(x_fine, t_min)
        cdf_max_idx = np.searchsorted(x_fine, t_max)
        cdf_min = cdf[cdf_min_idx] if cdf_min_idx < len(cdf) else 0
        cdf_max = cdf[cdf_max_idx] if cdf_max_idx < len(cdf) else 1
        # Inverse CDF sampling for quantiles
        quantiles = np.linspace(cdf_min, cdf_max, numPts)
        grid_values = np.interp(quantiles, cdf, x_fine)

    elif method == "optimal":
        # Use optimal binning to determine the number of points
        try:
            from optbinning import ContinuousOptimalBinning
        except ImportError:
            raise ImportError(
                "optbinning is required for method='optimal'. "
                "Install with: pip install optbinning"
            ) from None

        binning = ContinuousOptimalBinning(
            n_bins=20,
            min_bin_size=0.05,
            monotonic_trend=None
        )
        binning.fit(pseudotime.values)
        m = binning.n_bins
        # Use the bin boundaries as grid
        grid_values = binning.splits
        # Apply endpoint truncation
        grid_values = grid_values[(grid_values >= t_min) & (grid_values <= t_max)]
        # Ensure we have enough points
        if len(grid_values) < numPts:
            # Fall back to linspace with truncation
            grid_values = np.linspace(t_min, t_max, numPts)
        numPts = len(grid_values)

    else:
        raise ValueError(
            f"Unknown grid_method='{method}'. "
            f"Options: 'linspace', 'quantile', 'density', 'optimal', 'custom'."
        )

    grid_index = [f"Gp_{i}" for i in range(numPts)]
    return grid_values, grid_index


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
                   numPts=14, grid_method="linspace", grid_custom=None,
                   endpoint_mode="half-interval",
                   alpha=1e-6, length_scale=0.1,
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
        Number of interpolation points on normalized pseudotime. Default is 14.
        Ignored when ``grid_method='optimal'`` (auto-determined) or
        ``grid_method='custom'`` (uses ``grid_custom``).
    grid_method : str, optional
        Strategy for constructing the interpolation time grid. Options:
        - "linspace": evenly spaced points (default)
        - "quantile": quantile-based spacing (equal cells per bin)
        - "density": density-adaptive spacing (more points in dense regions)
        - "optimal": use optimal binning to determine numPts automatically
        - "custom": use user-provided grid (requires ``grid_custom``)
    grid_custom : array-like, optional
        Custom grid values in range [0, 1], strictly increasing.
        Only used when ``grid_method='custom'``.
    endpoint_mode : str, optional
        How to handle endpoints. Default is "half-interval". Options:
        - "half-interval": endpoints at [half_interval, 1-half_interval]
          to avoid sparse cells at t=0 and t=1 (recommended)
        - "full": endpoints at [0, 1] (legacy behavior)
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
    if alpha <= 0:
        raise ValueError("alpha must be positive.")
    if length_scale <= 0:
        raise ValueError("length_scale must be positive.")
    if n_restarts_optimizer < 0:
        raise ValueError("n_restarts_optimizer must be non-negative.")
    if len(length_scale_bounds) != 2 or length_scale_bounds[0] <= 0 or length_scale_bounds[0] >= length_scale_bounds[1]:
        raise ValueError("length_scale_bounds must be a positive (min, max) tuple with min < max.")
    if endpoint_mode not in ("half-interval", "full"):
        raise ValueError(
            f"Unknown endpoint_mode='{endpoint_mode}'. "
            f"Options: 'half-interval', 'full'."
        )

    target = adata if inplace else adata.copy()
    expr, pseudotime, num_cells_dropped = _prepare_gp_inputs(
        target,
        pseudotime_col=pseudotime_col,
        layer=layer,
        gene_mask_key=gene_mask_key
    )

    grid_values, grid_index = _build_time_grid(
        pseudotime, numPts=numPts, method=grid_method, grid_custom=grid_custom,
        endpoint_mode=endpoint_mode
    )
    numPts = len(grid_values)
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
            "grid_method": grid_method,
            "endpoint_mode": endpoint_mode,
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
