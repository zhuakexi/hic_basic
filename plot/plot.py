from itertools import repeat, zip_longest
from os import sched_rr_get_interval

import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import pandas as pd
import statsmodels.api as sm
import numpy as np

from .rna import _plot_gene_heatmap
from .hic import _plot_cdps, plot_compartment_strength
# cell cycle
## cdp scatter plot
def cdp_scatter(dat:pd.DataFrame)->go.Figure:
    fig = px.scatter(
        dat,
        x="dist", 
        y="cratio"
        )
    fig.update_layout(
        height=500,
        width=500
    )
    fig.update_xaxes(type="log")
    return fig

# cdps heatmap with cycle-phasing marker
def plot_cdps_mark(cdps,orig_annote,sample_col="sample_name",order_col="index_order",group_col="group",color_map=None,hline=False):
    """
    # plot cdps heatmap with group marker in different color
    # Input:
    #    cdps: contact decay profile, row for sample
    #    cell_annote: must contain group and ordering col
    #    sample_col: name of sample id column
    #    order_col: name of order column
    #    hline: whether to mark named distance interval
    # Output:
    #    plotly Figure obj
    """
    if color_map == None:
        # using default color map(for cell cycle)
        annote = orig_annote.assign(
            group_color = "blue"
        )
        annote.loc[annote[group_col]=="Post-M","group_color"] = "red"
        annote.loc[annote[group_col]=="Pre-M","group_color"] = "goldenrod"
        annote.loc[annote[group_col]=="early/mid-S","group_color"] = "green"
        annote.loc[annote[group_col]=="mid-S/G2","group_color"] = "lightgreen"
        annote = annote.sort_values(order_col)
    elif isinstance(color_map,dict):
        # using custom color map
        annote = orig_annote.assign(
            # grey for non assigned sample
            group_color = "grey"
        )
        for key in color_map:
            # assign group_color according to color_map
            annote.loc[
                annote[group_col] == key,
                "group_color"
            ] = color_map[key]
    else:
        raise ValueError("cdps_plot: Must provide color_map when set color_map to None.")
    order = list(annote[sample_col].values)
    fig = make_subplots(rows=2,cols=1,row_heights=[0.05,1],vertical_spacing=0.05)
    fig.add_trace(
        go.Bar(
        x = list(range(annote.shape[0])),
        y = [1 for i in range(0,annote.shape[0])],
        marker_color = annote["group_color"],
        hovertemplate='%{text}',
        text = order
        ),
        row=1,
        col=1
    )
    fig.add_trace(
        go.Heatmap(
            z = cdps.loc[order].T.values,
            hovertemplate = '%{text} '+'%{z}',
            text = [order for i in range(annote.shape[0])]
        ),
        row=2,
        col=1
    )
    fig.update_xaxes(
        row=1,
        col=1,
        visible=False
    )
    fig.update_yaxes(
        row=1,
        col=1,
        visible=False
    )
    fig.update_layout(
        height=500,
        width=1000
    )
    if hline != False:
        # mark mitotic
        fig.add_hline(
            y = 90,
            row = 2,
            col= 1
        )
        fig.add_hline(
            y = 109,
            row = 2,
            col= 1
        )
        # mark short% 
        fig.add_hline(
            y = 38,
            row = 2,
            col= 1
        )
    return fig

def add_cdps(fig, row, col, adata, sorted_obs):
    # add cdps subplot
    cdps_fig = _plot_cdps(
        adata.uns["cdps"].loc[sorted_obs.index].T
    )
    cdps_fig = cdps_fig.update_traces(
       showscale=False 
    )
    fig.add_trace(
        cdps_fig.data[0],
        row=row,
        col=col
    )
    return fig
def add_compartment_strength(fig, row, col, adata, sorted_obs):
    # add compartment strength subplot
    cs_fig = plot_compartment_strength(
        adata
        #adata.obs.sort_values("velocity_pseudotime")
    )
    for trace in cs_fig.data:
        fig.add_trace(trace,row=row,col=col)
    return fig
def add_pmUMI(fig, row, col, adata, sorted_obs):
    # add paternal maternal umi count subplot
    if "pmUMIs" not in adata.obs.columns:
        pmUMI_fig = px.scatter(
            sorted_obs,
            x=sorted_obs.index,
            y=sorted_obs["g1_UMIs"]/(sorted_obs["g1_UMIs"]+sorted_obs["g2_UMIs"]),
            #color = "seurat_clusters",
            title="genome1_ratio")
    else:
        pmUMI_fig = px.scatter(
            x = sorted_obs.index,
            y = adata.obs.loc[sorted_obs.index,"pmUMIs"],
            title = "genome1_ratio"
        )
    for trace in pmUMI_fig.data:
        fig.add_trace(
            #pmUMI_fig.data[0],
            trace,
            row = row,
            col = col
        )
    return fig
def add_contact_number(fig, row, col, adata, sorted_obs):
    # add contact number subplot
    con_num_fig = px.scatter(
        sorted_obs,
        x=sorted_obs.index,
        y=sorted_obs["pairs_c123_num"],
        #color = "seurat_clusters",
        title="c123 contact number")
    for trace in con_num_fig.data:
        fig.add_trace(
            #pmUMI_fig.data[0],
            trace,
            row = row,
            col = col
        )
    return fig
def add_intra(fig, row, col, adata, sorted_obs):
    # add intra subplot
    intra_fig = px.scatter(
        sorted_obs,
        x=sorted_obs.index,
        y=sorted_obs["intra"],
        #color = "seurat_clusters",
        title="intra contact ratio")
    for trace in intra_fig.data:
        fig.add_trace(
            #pmUMI_fig.data[0],
            trace,
            row = row,
            col = col
        )
    return fig
def add_rs(fig, row, col, adata, sorted_obs):
    rs_fig = px.scatter(
        sorted_obs,
        x=sorted_obs.index,
        y=sorted_obs["repli_score"],
        #color = "seurat_clusters",
        title="repli_score")
    for trace in rs_fig.data:
        fig.add_trace(
            #pmUMI_fig.data[0],
            trace,
            row = row,
            col = col
        )
    return fig
def add_schicluster_res(fig, row, col, adata, sorted_obs):
    # add schicluster subplot
    schicluster_fig = px.scatter(
        sorted_obs,
        x=sorted_obs.index,
        y=sorted_obs["schicluster_res"],
        #color = "seurat_clusters",
        title="schicluster_res")
    for trace in schicluster_fig.data:
        fig.add_trace(
            trace,
            row = row,
            col = col
        )
    return fig
def add_gene_heatmap(fig, row, col, adata, sorted_obs, gene_order):
    gene_fig = _plot_gene_heatmap(
        np.log(adata.uns["pCells"].loc[gene_order,sorted_obs.index]+1)
    )
    gene_fig = gene_fig.update_traces(
        showscale = False
    )
    fig.add_trace(
        gene_fig.data[0],
        row = row,
        col = col
    )
    return fig

def add_rs_config(fig, row, col, adata, sorted_obs):
    fig.update_xaxes(
        autorange=False, 
        range=[0,sorted_obs.shape[0]],
        row = row,
        col = col
    )
    return fig

def add_compartment_strength_config(fig, row, col, adata, sorted_obs):
    fig.update_xaxes(
        autorange=False, 
        range=[0,sorted_obs.shape[0]],
        row = row,
        col = col
    )
    return fig

def add_schicluster_res_config(fig, row, col, adata, sorted_obs):
    fig.update_xaxes(
        autorange=False, 
        range=[0,sorted_obs.shape[0]],
        row = row,
        col = col
    )
    return fig
def time_attr(adata, order_col="velocity_pseudotime", using=None, ascending=True, gene_order = None):
    """
    Plot the time attributes of the cells.
    Input:
        adata: AnnData object
        order_col: column in adata.obs to sort cells
        using: name of attrs to plot; list; cdps, rs, gene, pmUMI, cs
        ascending: whether to revert order
    """
    # Input:
    #  adata: AnnData object
    # plot attributes in single figure
    
    sorted_obs = adata.obs.sort_values(order_col, ascending=ascending)

    stub_trace = {
     "cdps":add_cdps,
        "rs":add_rs,
        "gene":add_gene_heatmap,
        "pmUMI":add_pmUMI,
        "cs":add_compartment_strength,
        "schicluster_res":add_schicluster_res,
    }
    stub_layout = {
        "rs" : add_rs_config,
        "cs" : add_compartment_strength_config,
        "schicluster_res" : add_schicluster_res_config,
    }
    if using == None:
        using = stub_trace.keys()
    fig = make_subplots(rows=len(using),cols=1,vertical_spacing=0.02,shared_xaxes=True)
    for i, v in enumerate(using):
        if v == "gene":
            fig = stub_trace[v](fig, i+1, 1, adata, sorted_obs, gene_order)
        else:
            fig = stub_trace[v](fig, i+1, 1, adata, sorted_obs)
    for i, v in enumerate(using):
        try:
            fig = stub_layout[v](fig, i+1, 1, adata, sorted_obs)
        except KeyError:
            continue
    # config figure
    fig.update_layout(
        height = 1500,
        width = 1000,
        title = "cell-order vs. attrs plot"
    )
    fig.update_yaxes(
        visible=False,
        row=1,
        col=1
    )
    for i in range(1,len(using)+1):
        fig.update_xaxes(
            visible=False,
            row=i,
            col=1
        )
    return fig
# plot pheatmap style gapped heatmap
def gen_gap(features:list, number:int, k:int):
    """
    Input:
        index: index of features
        m: number of gaps
        k: width of gaps
    Output:
        iterator of NA dataframes
    """
    n_features = len(features)
    for i in range(number):
        yield pd.DataFrame( 
                np.array(
                    list(repeat(np.nan, k*n_features))
                ).reshape(k, n_features),
                # fake sample_name
                index = ["gap-{fake_id}".format(fake_id = j) for j in range(i*k, i*k+k)],
                columns = features
            )
def _plot_cluster_gapped_heatmap(data, grouping, cluster_order, gap_width=1, T=True, hm_layout=_plot_cdps):
    """
    Plot gap between heatmap clusters using dataframe hack.
    TODO:
        1. handle cluster number == 1
        2. mask gap-xx xticks
        3. don't drop cluster or sample silently
    Input:
        data: m( samples ) * n( features ) dataframe
        cluster: cluster ID of each sample. will drop missing levels; pd.Series
        cluster_order: order of clusters to show on heatmap
        gap_width: each gap equals to `gap_width` samples
        T: whether to transpose data( usually in time-dependent process )
    Output:
        go.Figure
    """
    # only keep cluster mentioned in cluster_order
    data = data.loc[grouping[grouping.isin(cluster_order)].index]
    grouping = grouping.cat.set_categories(cluster_order, ordered=True)
    dfs = (v for k, v in data.groupby(grouping, sort=True))
    gaps = gen_gap(data.columns, grouping.cat.categories.shape[0] - 1, gap_width)
    intervening = (value for row in zip_longest(dfs, gaps, fillvalue=pd.DataFrame()) for value in row)
    compiled_df = pd.concat(intervening)
    if T == True:
        fig = hm_layout(compiled_df.T)
    else:
        fig = hm_layout(compiled_df)
    fig.update_layout(
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)'
    )
    return fig