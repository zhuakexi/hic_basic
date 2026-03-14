import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import pandas as pd
import numpy as np

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from typing import Union, List, Optional
import numpy as np

def plot_pca(
    adata,
    color: str,
    PCs: Optional[Union[List[int], List[List[int]]]] = None,
    ncols: int = 2,
    background: Optional[str] = None,
    **kwargs
):
    """
    Plot PCA results from AnnData object with flexible subplot configurations.

    This function creates scatter plots of PCA components from single-cell data.
    It can generate single plots or multiple subplots showing different combinations
    of principal components.

    Parameters
    ----------
    adata : AnnData
        AnnData object containing PCA results in ``obsm['X_pca']``.
    color : str
        Column name in ``adata.obs`` to use for coloring foreground points.
    PCs : list of int or list of list of int, optional
        Specification of which principal components to plot:

        - ``None``: plots first two PCs (PC1 vs PC2).
        - ``List[int]``: plots all pairwise combinations of the specified PCs.
        - ``List[List[int]]``: plots specific pairs of PCs as subplots.
    ncols : int, optional
        Number of columns for subplot layout when multiple plots are generated
        (default 2).
    background : str, optional
        A :meth:`pandas.DataFrame.query` expression (same syntax as
        ``DataFrame.query``) that selects samples to treat as *background*.
        Background samples are rendered in light grey and drawn behind the
        foreground points.  Foreground samples (those **not** matching the
        expression) are colored normally according to *color*.  Defaults to
        ``None`` (all samples are foreground).
    **kwargs
        Additional arguments passed to ``px.scatter`` for foreground traces.

    Returns
    -------
    plotly.graph_objects.Figure
        Figure object containing the PCA plot(s).

    Examples
    --------
    >>> # Plot first two PCs
    >>> fig = plot_pca(adata, color='cell_type')

    >>> # Plot first two PCs with a grey background group
    >>> fig = plot_pca(adata, color='cell_type', background="batch == 'B'")

    >>> # Plot all pairwise combinations of PCs 1-3
    >>> fig = plot_pca(adata, color='cell_type', PCs=[0, 1, 2])

    >>> # Plot specific PC pairs in subplots
    >>> fig = plot_pca(adata, color='cell_type', PCs=[[0, 1], [0, 2], [1, 3]])
    """
    # Extract PCA results
    pca_data = pd.DataFrame(adata.obsm["X_pca"])
    pca_data.index = adata.obs_names
    pca_data.columns = [f"PC{i+1}" for i in range(pca_data.shape[1])]

    # Add annotations
    pca_data = pd.concat([pca_data, adata.obs], axis=1)

    # Split into background / foreground based on query expression
    if background is not None:
        bg_idx = pca_data.query(background).index
        bg_data = pca_data.loc[bg_idx]
        fg_data = pca_data.drop(index=bg_idx)
    else:
        bg_data = None
        fg_data = pca_data

    # Determine which PC pairs to plot
    if PCs is None:
        pc_pairs = [[0, 1]]
    elif isinstance(PCs[0], list):
        pc_pairs = PCs
    else:
        pc_pairs = []
        for i in range(len(PCs)):
            for j in range(i + 1, len(PCs)):
                pc_pairs.append([PCs[i], PCs[j]])

    n_plots = len(pc_pairs)

    def _bg_trace(df, x_col, y_col, show_legend_entry=True):
        """Build a single light-grey scatter trace for background samples."""
        return go.Scatter(
            x=df[x_col],
            y=df[y_col],
            mode="markers",
            marker=dict(color="lightgrey"),
            name="Background",
            hovertext=df.index,
            showlegend=show_legend_entry,
            legendgroup="background",
        )

    if n_plots == 1:
        pc1, pc2 = pc_pairs[0]
        x_col = f"PC{pc1+1}"
        y_col = f"PC{pc2+1}"

        if bg_data is not None and not bg_data.empty:
            # Build foreground figure first to inherit its layout/colorscale
            fg_fig = px.scatter(
                fg_data,
                x=x_col,
                y=y_col,
                color=color,
                hover_name=fg_data.index,
                **kwargs
            )
            fig = go.Figure()
            # Background drawn first (behind foreground)
            fig.add_trace(_bg_trace(bg_data, x_col, y_col))
            for trace in fg_fig.data:
                fig.add_trace(trace)
            fig.update_layout(fg_fig.layout)
        else:
            fig = px.scatter(
                fg_data,
                x=x_col,
                y=y_col,
                color=color,
                hover_name=fg_data.index,
                **kwargs
            )

        fig.update_layout(
            height=500,
            width=600,
            title=f"{x_col} vs {y_col}"
        )

    else:
        # Multiple plots – create subplots
        nrows = int(np.ceil(n_plots / ncols))

        fig = make_subplots(
            rows=nrows,
            cols=ncols,
            subplot_titles=[f"PC{pc1+1} vs PC{pc2+1}" for pc1, pc2 in pc_pairs]
        )

        bg_legend_added = False  # show "Background" legend entry only once
        for i, (pc1, pc2) in enumerate(pc_pairs):
            row = i // ncols + 1
            col = i % ncols + 1

            x_col = f"PC{pc1+1}"
            y_col = f"PC{pc2+1}"

            if bg_data is not None and not bg_data.empty:
                # Background trace first
                fig.add_trace(
                    _bg_trace(bg_data, x_col, y_col,
                               show_legend_entry=not bg_legend_added),
                    row=row, col=col
                )
                bg_legend_added = True

                # Foreground traces
                scatter_fig = px.scatter(
                    fg_data,
                    x=x_col,
                    y=y_col,
                    color=color,
                    hover_name=fg_data.index,
                    **kwargs
                )
            else:
                scatter_fig = px.scatter(
                    fg_data,
                    x=x_col,
                    y=y_col,
                    color=color,
                    hover_name=fg_data.index,
                    **kwargs
                )

            for scatter_trace in scatter_fig.data:
                fig.add_trace(scatter_trace, row=row, col=col)

            fig.update_xaxes(title_text=x_col, row=row, col=col)
            fig.update_yaxes(title_text=y_col, row=row, col=col)

        base_height = 400
        base_width = 500
        fig.update_layout(
            height=base_height * nrows,
            width=base_width * ncols,
            showlegend=True,
            title_text="PCA Components Visualization"
        )

    return fig
def plot_elbow(adata):
    fig = px.scatter(adata.uns["pca"]["variance_ratio"])
    fig.update_layout(
        showlegend = False,
        height = 500,
        width = 600,
        title = "PCA elbow plot"
    )
    return fig
def plot_umap(adata, color, key="X_umap", **kwargs):
    """
    Plot scanpy umap using plotly.
    """
    umap_res = pd.DataFrame(adata.obsm[key][:, :2])
    umap_res.index = adata.obs_names
    umap_res.columns = ["U1","U2"]
    # adding annotations
    umap_res = pd.concat([umap_res, adata.obs], axis=1)
    # plot
    if umap_res[color].dtype == "category":
        umap_res = umap_res.sort_values(color)
    fig = px.scatter(umap_res, x = "U1", y = "U2", color = color, hover_name = umap_res.index, **kwargs)
    fig.update_layout(
        height = 500,
        width = 600
    )
    return fig
def plot_tsne(adata, color):
    """
    Plot scanpy tsne using plotly.
    """
    tsne_res = pd.DataFrame(adata.obsm["X_tsne"][:, :2])
    tsne_res.index = adata.obs_names
    tsne_res.columns = ["T1","T2"]
    # adding annotations
    tsne_res = pd.concat([tsne_res, adata.obs], axis=1)
    # plot
    fig = px.scatter(tsne_res, x = "T1", y = "T2", color = color, hover_name = tsne_res.index)
    fig.update_layout(
        height = 500,
        width = 600
    )
    return fig
def plot_diffmap(adata, color, title):
    # --- create dataframe for plotly ---
    obs = adata.obs.copy()
    n_DCs = adata.obsm["X_diffmap"].shape[1]
    diffmap = pd.DataFrame(
        adata.obsm["X_diffmap"],
        index = obs.index,
        columns = [f"DC{i}" for i in range(1, n_DCs+1)]
        )
    # --- plot ---
    fig = px.scatter(
        diffmap,
        x="DC2", # DC1 is steady-state solution, which is non-informative in diffusion maps
        y="DC3",
        color=obs[color],
        title=title,
        )
    fig.update_traces(marker=dict(size=3))
    return fig