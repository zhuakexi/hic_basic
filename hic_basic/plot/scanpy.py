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
    **kwargs
):
    """
    Plot PCA results from AnnData object with flexible subplot configurations.
    
    This function creates scatter plots of PCA components from single-cell data.
    It can generate single plots or multiple subplots showing different combinations
    of principal components.
    
    Args:
        adata: AnnData object containing PCA results in obsm['X_pca']
        color: Column name in adata.obs to use for coloring points
        PCs: Specification of which principal components to plot. Can be:
            - None: plots first two PCs (PC1 vs PC2)
            - List[int]: plots all pairwise combinations of the specified PCs
            - List[List[int]]: plots specific pairs of PCs as subplots
        ncols: Number of columns for subplot layout when multiple plots are generated
        **kwargs: Additional arguments passed to px.scatter
        
    Returns:
        plotly.graph_objects.Figure: Figure object containing the PCA plot(s)
        
    Examples:
        >>> # Plot first two PCs
        >>> fig = plot_pca(adata, color='cell_type')
        
        >>> # Plot all pairwise combinations of PCs 1-3
        >>> fig = plot_pca(adata, color='cell_type', PCs=[0, 1, 2])
        
        >>> # Plot specific PC pairs in subplots
        >>> fig = plot_pca(adata, color='cell_type', PCs=[[0, 1], [0, 2], [1, 3]])
        
        >>> # Customize subplot layout
        >>> fig = plot_pca(adata, color='cell_type', PCs=[0, 1, 2, 3], ncols=3)
    """
    # Extract PCA results
    pca_data = pd.DataFrame(adata.obsm["X_pca"])
    pca_data.index = adata.obs_names
    pca_data.columns = [f"PC{i+1}" for i in range(pca_data.shape[1])]
    
    # Add annotations
    pca_data = pd.concat([pca_data, adata.obs], axis=1)
    
    # Determine which PC pairs to plot
    if PCs is None:
        # Default: plot first two PCs
        pc_pairs = [[0, 1]]
    elif isinstance(PCs[0], list):
        # List of lists: use specified pairs directly
        pc_pairs = PCs
    else:
        # List of integers: generate all pairwise combinations
        pc_pairs = []
        for i in range(len(PCs)):
            for j in range(i + 1, len(PCs)):
                pc_pairs.append([PCs[i], PCs[j]])
    
    n_plots = len(pc_pairs)
    
    # Create appropriate figure based on number of plots
    if n_plots == 1:
        # Single plot - use original behavior
        pc1, pc2 = pc_pairs[0]
        x_col = f"PC{pc1+1}"
        y_col = f"PC{pc2+1}"
        
        fig = px.scatter(
            pca_data, 
            x=x_col, 
            y=y_col, 
            color=color,
            hover_name=pca_data.index,
            **kwargs
        )
        
        fig.update_layout(
            height=500,
            width=600,
            title=f"{x_col} vs {y_col}"
        )
        
    else:
        # Multiple plots - create subplots
        nrows = int(np.ceil(n_plots / ncols))
        
        fig = make_subplots(
            rows=nrows, 
            cols=ncols,
            subplot_titles=[f"PC{pc1+1} vs PC{pc2+1}" for pc1, pc2 in pc_pairs]
        )
        
        # Create individual scatter plots and add to subplots
        for i, (pc1, pc2) in enumerate(pc_pairs):
            row = i // ncols + 1
            col = i % ncols + 1
            
            x_col = f"PC{pc1+1}"
            y_col = f"PC{pc2+1}"
            
            # Create individual scatter plot
            scatter_fig = px.scatter(
                pca_data,
                x=x_col,
                y=y_col,
                color=color,
                hover_name=pca_data.index,
                **kwargs
            )
            
            # Extract the scatter trace
            for scatter_trace in scatter_fig.data:
                fig.add_trace(scatter_trace, row=row, col=col)
            
            # Update axis labels
            fig.update_xaxes(title_text=x_col, row=row, col=col)
            fig.update_yaxes(title_text=y_col, row=row, col=col)
        
        # Calculate appropriate figure size
        base_height = 400
        base_width = 500
        height = base_height * nrows
        width = base_width * ncols
        
        fig.update_layout(
            height=height,
            width=width,
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