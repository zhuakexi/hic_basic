import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import pandas as pd
import numpy as np

def sc_pca(adata, color):
    # extract pca result
    pca_res = pd.DataFrame(adata.obsm["X_pca"][:,:2])
    pca_res.index = adata.obs_names
    pca_res.columns = ["PC1","PC2"]
    # adding annotations
    pca_res = pd.concat([pca_res,adata.obs],axis=1)
    # plot
    fig = px.scatter(pca_res, x = "PC1", y= "PC2", color = color,hover_name=pca_res.index)
    fig.update_layout(
        height = 500,
        width = 600
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
def plot_umap(adata, color):
    """
    Plot scanpy umap using plotly.
    """
    umap_res = pd.DataFrame(adata.obsm["X_umap"][:, :2])
    umap_res.index = adata.obs_names
    umap_res.columns = ["U1","U2"]
    # adding annotations
    umap_res = pd.concat([umap_res, adata.obs], axis=1)
    # plot
    fig = px.scatter(umap_res, x = "U1", y = "U2", color = color, hover_name = umap_res.index)
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