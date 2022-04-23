from math import ceil

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
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