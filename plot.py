import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import pandas as pd

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
# plot cdps heatmap
def plot_cdps(cdps):
    fig = go.Figure()
    fig.add_trace(
        go.Heatmap(
            z = cdps,
            y = cdps.index
            #yaxis="y"
        )
    )
    fig.update_layout(
        height = 500,
        width = 1000
    )
    fig.update_yaxes(
        type = "log",
        range = [3,8.3]
    )
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