from math import ceil

import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import plotly.graph_objects as go
import cooltools.lib.plotting # import the `fall` cmap
import cooler
import numpy as np
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import bioframe

from .utils import filling_l2r_mpl, pcolormesh_45deg

# --- Hi-C heatmap plot ---
def _plot_mat_mpl(mat, title="", vmax=500, ignore_diags=True, donorm=True, cmap="fall", balancing=False):
    mat = mat.copy()
    if ignore_diags:
        np.fill_diagonal(mat, 0)
    if balancing:
        norm = LogNorm(vmax=vmax) # using 0.1 here in cooltools tutorial
    else:
        norm = LogNorm(vmin=1,vmax=vmax)
    plt.figure(figsize=(11, 10))
    plt.gcf().canvas.set_window_title("Contact matrix".format())
    plt.title(title)
    if donorm:
        return plt.imshow(
            mat,
            interpolation="none",
            #extent=[col_lo, col_hi, row_hi, row_lo],
            norm = norm,
            cmap=cmap
        )
    else:
        return plt.imshow(
            mat,
            interpolation="none",
            #extent=[col_lo, col_hi, row_hi, row_lo],
            cmap=cmap
        )
fall = [
    [0.00,'rgb(255, 255, 255)'],
    [0.10,'rgb(255, 255, 204)'],
    [0.20,'rgb(255, 237, 160)'],
    [0.30,'rgb(254, 217, 118)'],
    [0.40,'rgb(254, 178, 76)'],
    [0.50,'rgb(253, 141, 60)'],
    [0.60,'rgb(252, 78, 42)'],
    [0.70,'rgb(227, 26, 28)'],
    [0.80,'rgb(189, 0, 38)'],
    [0.90,'rgb(128, 0, 38)'],
    [1.00,'rgb(0, 0, 0)'],
]
log_fall = [
    [0,'rgb(255, 255, 255)'],
    [1/10**9,'rgb(255, 255, 204)'],
    [1/10**8,'rgb(255, 237, 160)'],
    [1/10**7,'rgb(254, 217, 118)'],
    [1/10**6,'rgb(254, 178, 76)'],
    [1/10**5,'rgb(253, 141, 60)'],
    [1/10**4,'rgb(252, 78, 42)'],
    [1/10**3,'rgb(227, 26, 28)'],
    [1/10**2,'rgb(189, 0, 38)'],
    [1/10**1,'rgb(128, 0, 38)'],
    [1/10**0,'rgb(0, 0, 0)'],
]
def _plot_mat(mat, title="", vmax=500, ignore_diags=True, donorm=True, cmap="fall", balancing=False):
    """
    TODO: make a real colorscale bar according to LogNorm
    """
    mat = mat.copy()
    if ignore_diags:
        np.fill_diagonal(mat, 0)
    if donorm:
        if balancing:
            mat = LogNorm(vmin=np.nanmin(mat) if np.nanmin(mat) > 0 else 0.0001, vmax=np.nanmax(mat), clip=True)(mat).data # balanced cooler(don't know why mat+=1 failed here)
        else:
            mat = LogNorm(vmin=1, vmax=vmax, clip=True)(mat).data
    fig = go.Figure()
    fig.add_trace(
        go.Heatmap(
            z = mat,
            colorscale=fall, # don't know why log_fall failed here
            showscale=False
        )
    )
    return fig
def plot_cool(cool, title="", region="chr1",vmax=100, balance=False, ignore_diags=True, donorm=True):
    """"
    Plot heatmap of single cooler file.
    Input:
        cool: path to cooler file.
        title: name of the plot.
        region: genome region to plot.
            "chr1" or "chr1:1000000-2000000" or ["chr1:1,000,000-2,000,000", "chr2"] or
            slice(0,-1) or [slice(0,1000000), slice(1000000,2000000)]
            Note: can only cross chrom boundaries if using slice.
        vmax: max z value.
        balance: whether to load balanced cooler matrix.
        norm: whether to use lognorm.
    """
    using_index = False
    clr = cooler.Cooler(cool)
    if isinstance(region, str):
        region = [region]
    elif isinstance(region, list):
        if isinstance(region[0], slice):
            using_index = True
        pass
    elif isinstance(region, slice):
        region = [region]
        using_index = True
    else:
        raise ValueError("region must be str or list or slice")
    if using_index:
        mat = clr.matrix(balance=balance)[region[0]] if len(region) == 1 else clr.matrix(balance=balance)[region[0], region[1]] # because can't unpack in slicer
    else:
        mat = clr.matrix(balance=balance).fetch(*region)
    return _plot_mat_mpl(
        mat,
        title = title,
        vmax = vmax,
        ignore_diags = ignore_diags,
        donorm = donorm
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
# Functions to help with plotting
from matplotlib.ticker import EngFormatter
bp_formatter = EngFormatter('b')
def format_ticks(ax, x=True, y=True, rotate=True):
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x',rotation=45)
plt.rcParams['font.size'] = 12
def plot_IS(clr,insulation_table,title,resolution=500e3,window=1.5e6,balance=True,chrom="chr2",start=10_500_000,steps=90,vmin=0.001,vmax=0.1):
    """
    Plot insulation score together with Hi-C matrix strata.
    Input:
        clr: cooler obj.
        insulation_table: insulation score table.
    """
    resolution, window = int(resolution), int(window)
    end = start + steps*window
    region = (chrom, start, end)
    norm = LogNorm(vmax=vmax, vmin=vmin)
    data = clr.matrix(balance=balance).fetch(region)
    insulation_table = pd.read_table(insulation_table)
    
    f, ax = plt.subplots(figsize=(18, 6))
    im = pcolormesh_45deg(ax, data, start=region[1], resolution=resolution, norm=norm, cmap='fall')
    ax.set_aspect(0.5)
    ax.set_ylim(0, 10*window)
    format_ticks(ax, rotate=False)
    ax.xaxis.set_visible(False)
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="1%", pad=0.1, aspect=6)
    plt.colorbar(im, cax=cax)

    insul_region = bioframe.select(insulation_table, region)

    ins_ax = divider.append_axes("bottom", size="50%", pad=0., sharex=ax)
    ins_ax.set_prop_cycle(plt.cycler("color", plt.cm.plasma(np.linspace(0,1,5))))
    ins_ax.plot(insul_region[['start', 'end']].mean(axis=1),
                insul_region['log2_insulation_score_'+str(window)],
                label=f'Window {window} bp')

    ins_ax.legend(bbox_to_anchor=(0., -1), loc='lower left', ncol=4);

    format_ticks(ins_ax, y=False, rotate=False)
    ax.set_xlim(region[1], region[2])
    plt.title(title)
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
def plot_cdps(adata, index_col="velocity_pseudotime", reverse=False):
    """
    Input:
        index_col: wich adata.obs col to sort sample by
        reverse: whether to reverse the order of the sample
    """
    if reverse:
        order = adata.obs.sort_values(index_col, ascending=False).index
    else:
        order = adata.obs.sort_values(index_col, ascending=True).index
    cdps = adata.uns["cdps"].loc[order].T
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

# --- HiC ps-curve plot ---
def _plot_ps_curve(cvd_merged, der, xlim=(1e3,1e8)):
    """
    Plot ps curve.
    Input:
        cvd_merged: 2 col dataframe
        der: array of ps-curve's derivative
    """
    f, axs = plt.subplots(
        figsize=(6.5,13),
        nrows=2,
        gridspec_kw={'height_ratios':[6,2]},
        sharex=True)
    ax = axs[0]
    ax.loglog(
        cvd_merged['s_bp'],
        cvd_merged['balanced.avg.smoothed.agg'],
        '-',
        markersize=5,
    )

    ax.set(
        ylabel='IC contact frequency',
        xlim=xlim
    )
    ax.set_aspect(1.0)
    ax.grid(lw=0.5)


    ax = axs[1]
    ax.semilogx(
        cvd_merged['s_bp'],
        der,
        alpha=0.5
    )

    ax.set(
        xlabel='separation, bp',
        ylabel='slope')

    ax.grid(lw=0.5)
    return f, ax