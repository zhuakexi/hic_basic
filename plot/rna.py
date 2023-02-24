from math import ceil
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import anndata as ad
import pandas as pd
import numpy as np
from .utils import filling_l2r_plotly
from .general import scatter_cols

def _module_mean(adata, gene_set, layer):
    """
    Get mean expression of a set of genes.
    """
    data = pd.DataFrame(
        adata.layers[layer],
        index = adata.obs_names,
        columns = adata.var_names
    )
    return data[gene_set].mean(axis=1)
def _traj_module_mean(adata, gene_sets, ps_col, layer="not_known"):
    """
    Get mean expression of module along trajectory. 
    Input:
        adata: anndata object, or dataframe
        gene_sets: dict, module names as keys, list of module gene names as values
        ps_col: pseudotime column name, must be in adata.obs
        layer: layer name of count matrix
    Return:
        dataframe, index as sorted sample names, columns as module names
    """
    if isinstance(adata, ad.AnnData):
        if layer == "X":
            #print("is X")
            mat = adata.to_df()
        else:
            mat = adata.layers[layer]
        #print(mat.shape)
        data = pd.DataFrame(
            mat,
            index = adata.obs_names,
            columns = adata.var_names
        )
        cell_order = adata.obs.sort_values(ps_col).index
    elif isinstance(adata, pd.DataFrame):
        data = adata
        cell_order = ps_col.sort_values().index
    means = []
    for key in gene_sets:
        valid_genes = data.columns.intersection(gene_sets[key])
        print("{}: target {} valid {} ".format(key, str(len(gene_sets[key])), str(len(valid_genes))))
        #print(key, ":", str(len(gene_sets[key])), ":",str(len(valid_genes)) )
        means.append(data[valid_genes].mean(axis=1))
    data = pd.concat(means, axis=1)
    data.columns = [key for key in gene_sets]
    data = data.loc[cell_order]
    return data
def plot_gene_module_mean(adata, gene_sets, ps_col, layer="not_known", trends=False, points=True):
    """
    Plot mean expression of module along trajectory. 
    Input:
        adata: anndata object, or dataframe
        gene_sets: dict, module names as keys, list of module gene names as values
        ps_col: pseudotime column name, must be in adata.obs
        layer: layer name of count matrix
    Return:
        plotly figure
    """
    return scatter_cols(_traj_module_mean(adata, gene_sets, ps_col, layer), trends, points)
def plot_gene_trend(adata, gene, order_col="velocity_pseudotime", additional=None):
    """
    Ploting expression trends along time.
    Must have pSmooth in uns; normalized and log1p X counts
    """
    fig = go.Figure()
    try:
        smooth_counts = np.log(adata.uns["pSmooth"].loc[gene, adata.uns["pSmooth"].columns[::-1]] + 1)
    except KeyError:
        print("No Smoothing Model Found.")
        smooth_counts = []
    try:
        gene_count = adata[adata.obs.sort_values(order_col).index,:].obs_vector(gene)
    except KeyError:
        if isinstance(additional, ad.AnnData):
            print("Using additional counts.")
            if gene not in additional.var_names:
                gene_count = []
            else:
                gene_count = additional[additional.obs.sort_values(order_col).index,:].obs_vector(gene)
        else:
            gene_count = []
    x_num = max(len(smooth_counts),len(gene_count))
    fig.add_trace(
        go.Scatter(
            y = smooth_counts,
            x = np.linspace(0, x_num, len(smooth_counts)),
            mode = "lines",
            name = gene + " GAM fitting"
        )
    )
    """fig.add_trace(
        go.Scatter(
            y = adata.uns["pCells"].loc[gene, adata.obs.sort_values("velocity_pseudotime").index],
            mode = "markers"
        )
    )"""

    fig.add_trace(
        go.Scatter(
            y = gene_count,
            x = np.linspace(0, x_num, len(gene_count)),
            mode = "markers",
            name = gene + " normalized counts"
        )
    )
    fig.update_layout(
        height = 500,
        width = 700,
        title = gene
    )
    return fig
def plot_gene_trends(data, additional, genes, order_col, ncols):
    genes = [i for i in genes if i in data.var_names]
    nrows = ceil(len(genes)/ncols)
    fig = make_subplots(rows = nrows, cols = ncols, shared_xaxes=True, shared_yaxes=True)
    for i, j, k, gene in filling_l2r_plotly(nrows, ncols, genes):
        if gene is not None:
            subfig = plot_gene_trend(data, gene, order_col, additional)
            for trace in subfig.data:
                fig.add_trace(
                    trace,
                    row = i,
                    col = j
                )
    return fig

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from pp import standard_scaler
def _plot_gene_heatmap(df, scale=True,**args):
    if scale == True:
        df = standard_scaler(df, axis=1, with_std = True)
    fig = go.Figure()
    fig.add_trace(
        go.Heatmap(
            z = df.values,
            #text = ["a" for i in range(st_df.shape[0])],
            colorscale = "RdBu_r",
            **args
        )
    )
    fig.update_layout(
        height = 700,
        width = 500
    )
    return fig