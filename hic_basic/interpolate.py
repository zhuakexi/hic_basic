import numpy as np
import pandas as pd
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, RBF
import anndata as ad
from scipy import sparse

from .metrics import euclidean_distances


# =============================================================================
# Feature trajectory analysis utilities (context-agnostic: RNA/Loop/etc.)
# =============================================================================

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
    Internal function: Visualize and compare multiple grid methods side-by-side.

    This is a debugging utility, not intended for production use.
    Use plot_gp_heatmap() for single-method visualization.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    if grid_methods is None:
        grid_methods = ["linspace", "quantile", "density"]

    # Get expression and pseudotime
    if layer == "X":
        expr = adata.X.T
    else:
        expr = adata.layers[layer].T
    if sparse.issparse(expr):
        expr = expr.toarray()

    pseudotime = adata.obs[pseudotime_col].values
    cell_order = np.argsort(pseudotime)
    expr_sorted = expr[:, cell_order]
    pseudotime_sorted = pseudotime[cell_order]

    if gene_mask_key is not None:
        gene_mask = adata.var[gene_mask_key].values
        expr_sorted = expr_sorted[gene_mask, :]

    # Normalize
    if scaling == "z_score":
        expr_mean = expr_sorted.mean(axis=1, keepdims=True)
        expr_std = expr_sorted.std(axis=1, keepdims=True) + 1e-6
        expr_norm = (expr_sorted - expr_mean) / expr_std
        cmap, vmin, vmax = "RdBu_r", -3, 3
    else:
        expr_min = expr_sorted.min(axis=1, keepdims=True)
        expr_max = expr_sorted.max(axis=1, keepdims=True)
        expr_norm = (expr_sorted - expr_min) / (expr_max - expr_min + 1e-6)
        cmap, vmin, vmax = "inferno", 0, 1

    n_methods = len(grid_methods)
    n_genes = expr_norm.shape[0]
    n_cells = expr_norm.shape[1]

    fig, axes = plt.subplots(
        n_methods, 2, figsize=figsize, dpi=dpi, sharey=True,
        gridspec_kw={"width_ratios": [3, 2], "hspace": 0.15, "wspace": 0.1}
    )
    if n_methods == 1:
        axes = axes.reshape(1, -1)

    cell_ranks = np.arange(n_cells) / (n_cells - 1)
    pt_sorted = np.sort(pseudotime)
    n_pt = len(pt_sorted)

    for row_idx, method in enumerate(grid_methods):
        grid, _ = _build_time_grid(pd.Series(pseudotime), numPts=numPts, method=method, endpoint_mode=endpoint_mode)

        # Left: Raw
        ax = axes[row_idx, 0]
        ax.imshow(expr_norm, aspect="auto", origin="lower", extent=[0, 1, 0, n_genes], cmap=cmap, vmin=vmin, vmax=vmax)
        for grid_t in grid:
            rank_pos = np.searchsorted(pt_sorted, grid_t) / (n_pt - 1)
            ax.axvline(x=rank_pos, color="black", linestyle="--", linewidth=0.8, alpha=0.7)
        ax.set_title(f"{method}: Raw", fontsize=10, fontweight="bold")
        ax.set_xlabel("Cell rank")
        if row_idx == 0:
            ax.set_ylabel("Gene index")

        # Right: Smoothed schematic
        ax = axes[row_idx, 1]
        smooth = np.zeros((n_genes, len(grid)))
        for i in range(n_genes):
            smooth[i, :] = np.interp(grid, cell_ranks, np.convolve(expr_sorted[i, :], np.ones(5)/5, mode="same"))
        if scaling == "z_score":
            smooth = (smooth - smooth.mean(axis=1, keepdims=True)) / (smooth.std(axis=1, keepdims=True) + 1e-6)
        else:
            smooth = (smooth - smooth.min(axis=1, keepdims=True)) / (smooth.max(axis=1, keepdims=True) + 1e-6)
        ax.imshow(smooth, aspect="auto", origin="lower", extent=[-0.5, len(grid)-0.5, 0, n_genes], cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_title(f"{method}: GP-smoothed", fontsize=10, fontweight="bold")
        ax.set_xlabel("Grid point index")
        ax.set_xticks(range(len(grid)))
        if row_idx == 0:
            ax.set_ylabel("Gene index")

    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax) if scaling == "z_score" else None
    fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), cax=cbar_ax, label="Z-score" if scaling == "z_score" else "Scaled")
    fig.suptitle(f"GP Interpolation: Grid Comparison (endpoint_mode='{endpoint_mode}')", y=0.995, fontsize=12, fontweight="bold")
    return fig, axes


def _plot_grid_methods_comparison(
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
    """Internal: Plot comparison of different grid methods."""
    fig, _ = _plot_grid_comparison(
        adata, pseudotime_col, layer=layer, gene_mask_key=gene_mask_key,
        numPts=numPts, grid_methods=grid_methods, endpoint_mode=endpoint_mode,
        scaling=scaling, figsize=figsize, dpi=dpi,
    )
    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
    return fig


def compute_continuous_peak_time(trajectory, grid_values):
    """
    Compute continuous peak time for each feature using parabolic interpolation.

    For each feature, fits a parabola through the peak grid point and its two
    neighbors, then computes the vertex position as the continuous peak time.

    Parameters
    ----------
    trajectory : np.ndarray
        Trajectory matrix of shape (n_features, n_grid_points).
        Each row is the GP-smoothed trajectory for one feature.
    grid_values : np.ndarray
        Grid point values of shape (n_grid_points,). The actual time/pseudotime
        values at each grid point.

    Returns
    -------
    peak_time : np.ndarray
        Continuous peak time for each feature, shape (n_features,).
        If a feature has no clear peak (monotonic or flat), returns the grid
        position of the maximum value.
    peak_value : np.ndarray
        The interpolated peak value at the continuous peak time.
    """
    trajectory = np.asarray(trajectory)
    grid_values = np.asarray(grid_values)
    
    n_features, n_grid = trajectory.shape
    
    if n_grid < 3:
        # Not enough points for parabolic interpolation
        peak_idx = np.argmax(trajectory, axis=1)
        return grid_values[peak_idx], trajectory[np.arange(n_features), peak_idx]
    
    peak_time = np.zeros(n_features)
    peak_value = np.zeros(n_features)
    
    for i in range(n_features):
        y = trajectory[i, :]
        peak_idx = np.argmax(y)
        
        # Handle boundary cases
        if peak_idx == 0:
            peak_idx = 1
        elif peak_idx == n_grid - 1:
            peak_idx = n_grid - 2
        
        # Get three points for parabolic fit: (x_{i-1}, y_{i-1}), (x_i, y_i), (x_{i+1}, y_{i+1})
        x = grid_values[peak_idx - 1 : peak_idx + 2]
        y_local = y[peak_idx - 1 : peak_idx + 2]
        
        # Fit parabola: y = ax^2 + bx + c
        # Using Lagrange interpolation formula for vertex:
        # x_vertex = -b/(2a)
        # where a, b are coefficients from quadratic fit
        
        # Alternative: use the closed-form solution for 3-point parabola vertex
        # Given points (x0,y0), (x1,y1), (x2,y2), the vertex x-coordinate is:
        # x_v = (x1^2 - x0^2)(y2 - y1) - (x2^2 - x1^2)(y1 - y0)
        #       ---------------------------------------------------
        #       2 * [(x1 - x0)(y2 - y1) - (x2 - x1)(y1 - y0)]
        
        x0, x1, x2 = x
        y0, y1, y2 = y_local
        
        denom = (x1 - x0) * (y2 - y1) - (x2 - x1) * (y1 - y0)
        
        if np.abs(denom) < 1e-10:
            # Parabola is degenerate (points are collinear or flat)
            # Fall back to discrete peak
            peak_time[i] = grid_values[peak_idx]
            peak_value[i] = y[peak_idx]
        else:
            x_vertex = ((x1**2 - x0**2) * (y2 - y1) - (x2**2 - x1**2) * (y1 - y0)) / (2 * denom)
            
            # Clamp to the interpolation range
            x_vertex = np.clip(x_vertex, x0, x2)
            
            # Compute interpolated value at vertex: y = ax^2 + bx + c
            # Using Lagrange interpolation:
            a_coef = (y2 - y1) / (x2 - x1) - (y1 - y0) / (x1 - x0)
            a_coef = a_coef / (x2 - x0)
            b_coef = (y1 - y0) / (x1 - x0) - a_coef * (x0 + x1)
            c_coef = y0 - a_coef * x0**2 - b_coef * x0
            
            y_vertex = a_coef * x_vertex**2 + b_coef * x_vertex + c_coef
            
            peak_time[i] = x_vertex
            peak_value[i] = y_vertex
    
    return peak_time, peak_value


def compute_center_of_mass(trajectory, grid_values):
    """
    Compute center of mass (weighted average time position) for each feature.

    The center of mass is computed as:
        com = sum(grid_values * |trajectory|) / sum(|trajectory|)

    Using absolute values to handle both positive and negative trajectories.

    Parameters
    ----------
    trajectory : np.ndarray
        Trajectory matrix of shape (n_features, n_grid_points).
    grid_values : np.ndarray
        Grid point values of shape (n_grid_points,).

    Returns
    -------
    com : np.ndarray
        Center of mass for each feature, shape (n_features,).
    """
    trajectory = np.asarray(trajectory)
    grid_values = np.asarray(grid_values)
    
    # Use absolute values for weighting (handles both up/down regulation)
    weights = np.abs(trajectory)
    
    # Avoid division by zero
    weight_sum = weights.sum(axis=1, keepdims=True)
    weight_sum = np.where(weight_sum < 1e-10, 1e-10, weight_sum)
    
    com = (grid_values * weights).sum(axis=1) / weight_sum.squeeze()
    
    return com


def sort_features_by_trajectory(trajectory, grid_values, 
                                 sort_keys=("peak_time", "com", "-peak_value")):
    """
    Sort features based on trajectory dynamics characteristics.

    Parameters
    ----------
    trajectory : np.ndarray
        Trajectory matrix of shape (n_features, n_grid_points).
    grid_values : np.ndarray
        Grid point values of shape (n_grid_points,).
    sort_keys : tuple of str
        List of sort keys in priority order. Supported keys:
        - "peak_time": continuous peak time (parabolic interpolation)
        - "com": center of mass
        - "peak_value": peak value at continuous peak time
        - "-peak_value": negative peak value (for descending sort)
        - "max_grid_idx": discrete grid index of maximum (fallback)

    Returns
    -------
    sorted_indices : np.ndarray
        Indices that would sort the features.
    sort_values : pd.DataFrame
        DataFrame with columns for each computed metric.
    """
    trajectory = np.asarray(trajectory)
    n_features = trajectory.shape[0]
    
    # Compute all possible metrics
    peak_time, peak_value = compute_continuous_peak_time(trajectory, grid_values)
    com = compute_center_of_mass(trajectory, grid_values)
    max_grid_idx = np.argmax(trajectory, axis=1)
    
    # Build sort values DataFrame
    sort_values = pd.DataFrame({
        "peak_time": peak_time,
        "com": com,
        "peak_value": peak_value,
        "-peak_value": -peak_value,
        "max_grid_idx": max_grid_idx,
    })
    
    # Build sort keys list, handling negative signs
    keys_to_sort = []
    for key in sort_keys:
        if key in sort_values.columns:
            keys_to_sort.append(key)
        else:
            raise ValueError(f"Unknown sort key: '{key}'. "
                           f"Supported: {list(sort_values.columns)}")
    
    # Compute sorted indices using lexsort (last key is primary)
    if len(keys_to_sort) == 0:
        sorted_indices = np.arange(n_features)
    else:
        # lexsort sorts by last column first, so reverse the order
        sort_arrays = [sort_values[key].values for key in reversed(keys_to_sort)]
        sorted_indices = np.lexsort(sort_arrays)
    
    return sorted_indices, sort_values


# =============================================================================
# Visualization functions
# =============================================================================

def plot_trajectory_heatmap(
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
    Internal function: Visualize and compare multiple grid methods side-by-side.

    This is a debugging utility, not intended for production use.
    Use plot_gp_heatmap() for single-method visualization.

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

        # Add grid lines at positions corresponding to each grid point's rank
        # For each grid point t_i, find the fraction of cells with pseudotime < t_i
        # This correctly shows where each method places its grid points in rank space
        pseudotime_sorted = np.sort(pseudotime)
        n_cells_sorted = len(pseudotime_sorted)
        for grid_t in grid:
            # Find rank position: fraction of cells with pseudotime < grid_t
            rank_pos = np.searchsorted(pseudotime_sorted, grid_t) / (n_cells_sorted - 1)
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


def _plot_grid_methods_comparison(
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
    Internal function: Plot comparison of different grid methods for GP interpolation.

    This is a debugging utility, not intended for production use.
    Use plot_gp_heatmap() for single-method visualization.

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


def plot_trajectory_heatmap(
    trajectory_raw,
    pseudotime,
    trajectory_gp=None,
    grid_values=None,
    grid_cell_positions=None,
    feature_order=None,
    feature_labels=None,
    scaling="z_score",
    cmap="RdBu_r",
    vmin=None,
    vmax=None,
    figsize=(8, 6),
    dpi=100,
    save_path=None,
    show_grid=True,
    grid_line_color="black",
    grid_line_style="--",
    grid_line_alpha=0.7,
    show_colorbar=True,
    left_title="Raw (cell space)",
    right_title="GP-smoothed (grid space)",
    suptitle=None,
):
    """
    Plot trajectory heatmap for features (genes, loops, etc.) with GP-smoothed trajectories.

    This is a context-agnostic visualization function that works with any feature type
    (RNA genes, chromatin loops, etc.) as long as the data is provided as numpy arrays.

    Parameters
    ----------
    trajectory_raw : np.ndarray or pd.DataFrame
        Raw trajectory matrix of shape (n_features, n_cells). Each row is a feature's
        trajectory across cells. Can also be a DataFrame with feature names as index.
    pseudotime : np.ndarray or pd.Series
        Pseudotime values for cells, shape (n_cells,). Used to sort cells and compute
        grid line positions.
    trajectory_gp : np.ndarray, optional
        GP-smoothed trajectory matrix of shape (n_features, n_grid_points).
        If None, only the raw trajectory panel is shown.
    grid_values : np.ndarray, optional
        Grid point values in pseudotime space, shape (n_grid_points,).
        Used for right panel x-axis labels and to compute grid_cell_positions if needed.
    grid_cell_positions : np.ndarray, optional
        Grid point positions in cell rank space [0, 1], shape (n_grid_points,).
        If None but grid_values is provided, computed from pseudotime.
    feature_order : np.ndarray, optional
        Indices to reorder features (rows). If None, features are shown in original order.
        Typically obtained from sort_features_by_trajectory().
    feature_labels : list or np.ndarray, optional
        Labels for features (rows). If None, uses integer indices.
    scaling : str, optional
        Color scaling method. Options:
        - "z_score": Z-score normalization per feature
        - "zero_one": Min-max normalization to [0, 1]
    cmap : str or matplotlib.colors.Colormap, optional
        Colormap for the heatmap. Default is "RdBu_r" (diverging).
    vmin, vmax : float, optional
        Color scale limits. If None, computed from data based on scaling.
    figsize : tuple, optional
        Figure size (width, height) in inches. Default is (8, 6).
    dpi : int, optional
        Figure DPI. Default is 100.
    save_path : str, optional
        Path to save the figure. If None, figure is not saved.
    show_grid : bool, optional
        Show vertical lines at grid point positions. Default is True.
    grid_line_color : str, optional
        Color of grid lines. Default is "black".
    grid_line_style : str, optional
        Line style of grid lines. Default is "--".
    grid_line_alpha : float, optional
        Alpha transparency of grid lines. Default is 0.7.
    show_colorbar : bool, optional
        Show colorbar. Default is True.
    left_title : str, optional
        Title for left panel (raw data). Default is "Raw (cell space)".
    right_title : str, optional
        Title for right panel (GP-smoothed). Default is "GP-smoothed (grid space)".
    suptitle : str, optional
        Figure super title. If None, auto-generated from parameters.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The matplotlib figure object.
    axes : np.ndarray
        Array of matplotlib axes objects (1 or 2 elements depending on trajectory_gp).

    Examples
    --------
    >>> # RNA/gene trajectory visualization
    >>> gp_result = gp_interpolate(adata, pseudotime_col="pt")
    >>> gp_mean = gp_result.uns["gp_interpolate"]["mean"].values
    >>> grid_vals = gp_result.uns["gp_interpolate"]["time_grid"].values
    >>> sorted_idx, _ = sort_features_by_trajectory(gp_mean, grid_vals)
    >>> fig, axes = plot_trajectory_heatmap(
    ...     trajectory_raw=adata.X.T,
    ...     pseudotime=adata.obs["pt"],
    ...     trajectory_gp=gp_mean,
    ...     grid_values=grid_vals,
    ...     feature_order=sorted_idx
    ... )

    >>> # Loop trajectory visualization (same interface)
    >>> gp_result = gp_interpolate(loop_adata, pseudotime_col="pt")
    >>> gp_mean = gp_result.uns["gp_interpolate"]["mean"].values
    >>> sorted_idx, _ = sort_features_by_trajectory(gp_mean, grid_vals)
    >>> fig, axes = plot_trajectory_heatmap(
    ...     trajectory_raw=loop_adata.X.T,
    ...     pseudotime=loop_adata.obs["pt"],
    ...     trajectory_gp=gp_mean,
    ...     grid_values=grid_vals,
    ...     feature_order=sorted_idx
    ... )
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    # Convert inputs to numpy arrays
    if isinstance(trajectory_raw, pd.DataFrame):
        if feature_labels is None:
            feature_labels = trajectory_raw.index.values
        trajectory_raw = trajectory_raw.values
    elif sparse.issparse(trajectory_raw):
        # Handle sparse matrix: convert to dense BEFORE np.asarray
        trajectory_raw = trajectory_raw.toarray()
    else:
        trajectory_raw = np.asarray(trajectory_raw)
    
    pseudotime = np.asarray(pseudotime)
    n_features, n_cells = trajectory_raw.shape
    
    # Apply feature ordering
    if feature_order is not None:
        trajectory_raw = trajectory_raw[feature_order, :]
        if feature_labels is not None:
            feature_labels = np.asarray(feature_labels)[feature_order]
    
    # Sort cells by pseudotime
    cell_order = np.argsort(pseudotime)
    trajectory_sorted = trajectory_raw[:, cell_order]
    pseudotime_sorted = np.sort(pseudotime)
    
    # Normalize expression
    if scaling == "z_score":
        traj_mean = trajectory_sorted.mean(axis=1, keepdims=True)
        traj_std = trajectory_sorted.std(axis=1, keepdims=True) + 1e-6
        trajectory_norm = (trajectory_sorted - traj_mean) / traj_std
        if vmin is None:
            vmin = -3
        if vmax is None:
            vmax = 3
    else:  # zero_one
        traj_min = trajectory_sorted.min(axis=1, keepdims=True)
        traj_max = trajectory_sorted.max(axis=1, keepdims=True)
        trajectory_norm = (trajectory_sorted - traj_min) / (traj_max - traj_min + 1e-6)
        if vmin is None:
            vmin = 0
        if vmax is None:
            vmax = 1

    # Compute grid cell positions if needed
    if grid_cell_positions is None and grid_values is not None:
        n = len(pseudotime_sorted)
        grid_cell_positions = np.searchsorted(pseudotime_sorted, grid_values) / (n - 1)

    # Determine number of panels
    show_gp = trajectory_gp is not None and grid_values is not None
    n_panels = 2 if show_gp else 1
    
    # Create figure
    if n_panels == 2:
        fig, axes = plt.subplots(
            1, 2,
            figsize=figsize,
            dpi=dpi,
            sharey=True,
            gridspec_kw={"width_ratios": [3, 2], "hspace": 0.1, "wspace": 0.1}
        )
        axes = np.asarray(axes)
    else:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        axes = np.asarray([ax])
    
    # Left panel: Raw data with grid lines
    ax = axes[0]
    extent = [0, 1, 0, n_features]
    im = ax.imshow(
        trajectory_norm,
        aspect="auto",
        origin="lower",
        extent=extent,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
    )
    
    # Add grid lines
    if show_grid and grid_cell_positions is not None:
        for pos in grid_cell_positions:
            ax.axvline(x=pos, color=grid_line_color, linestyle=grid_line_style,
                      linewidth=0.8, alpha=grid_line_alpha)
    
    ax.set_title(left_title, fontsize=11, fontweight="bold")
    ax.set_xlabel("Cell rank (sorted by pseudotime)")
    ax.set_ylabel("Feature index" if feature_labels is None else "Features")

    # Right panel: GP-smoothed
    if show_gp:
        trajectory_gp = np.asarray(trajectory_gp)
        
        # Apply feature ordering to GP trajectory
        if feature_order is not None:
            trajectory_gp = trajectory_gp[feature_order, :]
        
        n_grid = trajectory_gp.shape[1]
        grid_indices = np.arange(n_grid)
        
        # Normalize GP trajectory
        if scaling == "z_score":
            gp_mean = trajectory_gp.mean(axis=1, keepdims=True)
            gp_std = trajectory_gp.std(axis=1, keepdims=True) + 1e-6
            trajectory_gp_norm = (trajectory_gp - gp_mean) / gp_std
        else:  # zero_one
            gp_min = trajectory_gp.min(axis=1, keepdims=True)
            gp_max = trajectory_gp.max(axis=1, keepdims=True)
            trajectory_gp_norm = (trajectory_gp - gp_min) / (gp_max - gp_min + 1e-6)
        
        ax = axes[1]
        extent = [-0.5, n_grid - 0.5, 0, n_features]
        im = ax.imshow(
            trajectory_gp_norm,
            aspect="auto",
            origin="lower",
            extent=extent,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )
        
        ax.set_title(right_title, fontsize=11, fontweight="bold")
        ax.set_xlabel("Grid point index")
        ax.set_xticks(grid_indices)
        ax.set_xticklabels(grid_indices.astype(str))
        
        # Add grid lines at integer positions (grid point boundaries)
        if show_grid:
            for i in range(n_grid + 1):
                ax.axvline(x=i - 0.5, color=grid_line_color, linestyle=":",
                          linewidth=0.6, alpha=0.3)

    # Add colorbar
    if show_colorbar:
        cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
        norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax) if scaling == "z_score" else None
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, cax=cbar_ax, label="Z-score" if scaling == "z_score" else "Scaled [0, 1]")

    # Add super title
    if suptitle is None:
        if show_gp and grid_values is not None:
            suptitle = f"Trajectory Heatmap ({len(grid_values)} grid points)"
        else:
            suptitle = "Trajectory Heatmap (Raw)"
    
    if suptitle:
        fig.suptitle(suptitle, y=0.995, fontsize=12, fontweight="bold")

    if save_path is not None:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")

    return fig, axes


def plot_gp_heatmap(
    adata,
    pseudotime_col,
    layer="X",
    gene_mask_key=None,
    numPts=14,
    grid_method="linspace",
    endpoint_mode="half-interval",
    scaling="z_score",
    figsize=(8, 6),
    dpi=100,
    save_path=None,
    show_grid=True,
    show_colorbar=True,
    feature_order=None,
):
    """
    Plot a single heatmap for GP interpolation results.

    This is a convenience wrapper around plot_trajectory_heatmap() for AnnData-based
    RNA/gene expression data. For more control or for non-RNA data (e.g., loops),
    use plot_trajectory_heatmap() directly.

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
    grid_method : str, optional
        Grid construction method. Default is "linspace".
        Options: "linspace", "quantile", "density", "optimal", "custom".
    endpoint_mode : str, optional
        Endpoint handling mode. Default is "half-interval".
        Options: "half-interval", "full".
    scaling : str, optional
        Color scaling. Default is "z_score".
        Options: "z_score", "zero_one".
    figsize : tuple, optional
        Figure size (width, height) in inches. Default is (8, 6).
    dpi : int, optional
        Figure DPI. Default is 100.
    save_path : str, optional
        Path to save the figure. If None, figure is not saved.
    show_grid : bool, optional
        Show vertical lines at grid point positions. Default is True.
    show_colorbar : bool, optional
        Show colorbar. Default is True.
    feature_order : np.ndarray, optional
        Indices to reorder features (genes). If None, genes are shown in original order.
        Typically obtained from sort_features_by_trajectory().

    Returns
    -------
    fig : matplotlib.figure.Figure
        The matplotlib figure object.
    axes : np.ndarray
        Array of matplotlib axes objects.

    See Also
    --------
    plot_trajectory_heatmap : General-purpose trajectory heatmap for any feature type.
    """
    # Get expression and pseudotime
    if layer == "X":
        expr = adata.X.T
    else:
        expr = adata.layers[layer].T
    if sparse.issparse(expr):
        expr = expr.toarray()

    pseudotime = adata.obs[pseudotime_col].values
    
    # Apply gene mask if provided
    if gene_mask_key is not None:
        gene_mask = adata.var[gene_mask_key].values
        expr = expr[gene_mask, :]
    
    n_genes, n_cells = expr.shape
    
    # Generate grid and compute grid cell positions
    grid, _ = _build_time_grid(
        pd.Series(pseudotime),
        numPts=numPts,
        method=grid_method,
        endpoint_mode=endpoint_mode
    )
    
    # Compute grid cell positions for visualization
    pseudotime_sorted = np.sort(pseudotime)
    grid_cell_positions = np.searchsorted(pseudotime_sorted, grid) / (n_cells - 1)
    
    # Create schematic GP-smoothed trajectory (for visualization without running GP)
    # Note: This is a schematic approximation. For real GP results, use plot_trajectory_heatmap()
    # with trajectory_gp from gp_interpolate() output.
    cell_ranks = np.arange(n_cells) / (n_cells - 1)
    expr_sorted = expr[:, np.argsort(pseudotime)]
    
    smooth_schematic = np.zeros((n_genes, len(grid)))
    for i in range(n_genes):
        raw_gene = expr_sorted[i, :]
        smooth_schematic[i, :] = np.interp(
            grid,
            cell_ranks,
            np.convolve(raw_gene, np.ones(5) / 5, mode="same")
        )
    
    # Use plot_trajectory_heatmap for actual plotting
    feature_labels = adata.var_names[gene_mask].values if gene_mask_key is not None else adata.var_names.values
    
    fig, axes = plot_trajectory_heatmap(
        trajectory_raw=expr,
        pseudotime=pseudotime,
        trajectory_gp=smooth_schematic,
        grid_values=grid,
        grid_cell_positions=grid_cell_positions,
        feature_order=feature_order,
        feature_labels=feature_labels,
        scaling=scaling,
        figsize=figsize,
        dpi=dpi,
        save_path=save_path,
        show_grid=show_grid,
        show_colorbar=show_colorbar,
        left_title=f"{grid_method}: Raw (grid lines)",
        right_title=f"{grid_method}: GP-smoothed ({len(grid)} pts)",
        suptitle=f"GP Interpolation: {grid_method} (endpoint_mode='{endpoint_mode}')"
    )
    
    return fig, axes


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
        # Grid positions in cell rank space [0, 1] for visualization
        "grid_cell_positions": pd.Series(
            np.searchsorted(np.sort(pseudotime.values), grid_values) / (len(pseudotime) - 1),
            index=grid_index,
            name="cell_rank_position"
        ),
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
