from math import ceil

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import plotly.graph_objects as go
import cooler
import cooltools.lib.plotting

from .utils import filling_l2r_mpl

# --- Hi-C heatmap plot ---
def _plot_mat(mat,title="",vmax=500):
    norm = LogNorm(vmin=1, vmax=vmax)
    plt.figure(figsize=(11, 10))
    plt.gcf().canvas.set_window_title("Contact matrix".format())
    plt.title(title)
    return plt.imshow(
        mat,
        interpolation="none",
        #extent=[col_lo, col_hi, row_hi, row_lo],
        norm = norm,
        cmap="fall"
    )
def plot_cool(cool, title="", region="chr1",vmax=100, balance=False):
    """"
    Plot heatmap of single cooler file.
    Input:
        cool: path to cooler file.
        title: name of the plot.
        region: genome region to plot.
        vmax: max z value.
        balance: whether to load balanced cooler matrix.
    """
    clr = cooler.Cooler(cool)
    return _plot_mat(
        clr.matrix(balance=balance).fetch(region),
        title = title,
        vmax = vmax
    )

def plot_cools(cools, region, titles, ncols=3, vmax=500, height=50, width=50):
    """
    Plot heatmap of multiple cooler files.
    Input:
        cools: list of cooler file paths.
        region: genome region to plot.
        titles: title of each subplot(each sample).
        ncols: subplot cols.
        vmax: max z value.
    """
    norm = LogNorm(vmin=1, vmax=vmax)
    nrows = ceil(len(cools)/ncols)
    f, axs = plt.subplots(
        figsize=(height, width),
        nrows = nrows,
        ncols = ncols,
        sharex = True, sharey = True
    )
    for i, j, k, cool in filling_l2r_mpl(nrows, ncols, cools):
        if cool is not None:
            ax = axs[i,j]
            ax.set_title(
                titles[k])
            im = ax.matshow(
                cooler.Cooler(cool).matrix(balance=False).fetch(region),
                norm=norm,
                cmap = "fall"
            )
        else:
            ax = axs[i, j]
            ax.set_visible(False)
    plt.tight_layout()
    return plt

# --- HIC contact decay profile plot ---
# plot cdps heatmap
def _plot_cdps(cdps:pd.DataFrame):
    """
    Input:
        cdps: sample x bin contact decay profile dataframe
    """
    fig = go.Figure()
    fig.add_trace(
        go.Heatmap(
            z = cdps,
            y = cdps.index,
            colorscale="bluyl"
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
def plot_cdps(adata, index_col="velocity_pseudotime"):
    """
    Input:
        index_col: wich adata.obs col to sort sample by
    """
    cdps = adata.uns["cdps"].loc[adata.obs.sort_values(index_col).index].T
    return _plot_cdps(cdps)

# --- HIC compartment strength plot ---
def _plot_compartment_strength(cs,g1_col="g1",g2_col="g2"):
    """
    # plot compartment strength
    # cs index must be sample_name
    # cs must have g1 and g2 col
    # cs must be sorted
    """
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x = cs.index,
            y = cs[g1_col],
            mode = "markers",
            name = "g1"
        )
    )
    fig.add_trace(
        go.Scatter(
            x = cs.index,
            y = cs[g2_col],
            mode = "markers",
            name = "g2"
        )
    )
    g1line = sm.nonparametric.lowess(
        exog=list(range(cs.shape[0])),
        endog=cs[g1_col],
        frac=0.2)
    fig.add_trace(
        go.Scatter(
            x = cs.index,
            y = g1line[:,1],
            name = "g1_trend"
        )
    )
    g2line = sm.nonparametric.lowess(
        exog=list(range(cs.shape[0])),
        endog=cs[g2_col],
        frac=0.2)
    fig.add_trace(
        go.Scatter(
            x = cs.index,
            y = g2line[:,1],
            name = "g2_trend"
        )
    )
    fig.update_layout(
        height = 500,
        width = 800,
        title = "cell-order vs. compartment strength"
    )
    return fig
def plot_compartment_strength(adata, index_col="velocity_pseudotime"):
    order = adata.obs.sort_values(index_col).index
    data = pd.concat([adata.uns["g1cs_500k"].loc[order, "20"], adata.uns["g2cs_500k"].loc[order, "20"]],axis=1)
    data.columns = ["g1cs","g2cs"]
    return _plot_compartment_strength(data, "g1cs", "g2cs")