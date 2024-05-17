from math import ceil

import bioframe
import cooler
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import seaborn as sns
import statsmodels.api as sm
from cooler import Cooler
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from plotly.subplots import make_subplots
from skimage import exposure
import cooltools.lib.plotting

from .utils import filling_l2r_mpl, pcolormesh_45deg, tiling_mat
from ..compartment import compartments
from ..coolstuff import cool2mat
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
    #plt.gcf().canvas.set_window_title("Contact matrix".format())
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
fruitpunch = sns.blend_palette(['white', 'red'], as_cmap=True)
num_colors = 256
colors = [fruitpunch(i) for i in range(num_colors)]
fruitpunch = ['rgb({},{},{})'.format(int(c[0]*255), int(c[1]*255), int(c[2]*255)) for c in colors]
def _plot_mat(orig_mat, title="", vmax=500, ignore_diags=True, donorm=False, cmap="fall", balancing=False, fillna=False, **args):
    """
    TODO:
        1.make a real colorscale bar according to LogNorm
        2. make real zmax
    """
    mat = orig_mat.copy()
    if isinstance(mat, pd.DataFrame):
        mat = mat.values
        if (len(orig_mat.index.levels) > 1) or (len(orig_mat.columns.levels) > 1):
            columns = None
            index = None
        else:
            columns = orig_mat.columns
            index = orig_mat.index
    else:
        columns = None
        index = None
    if ignore_diags:
        np.fill_diagonal(mat, 0)
    if donorm:
        if balancing:
            mat = LogNorm(vmin=np.nanmin(mat) if np.nanmin(mat) > 0 else 0.0001, vmax=np.nanmax(mat), clip=True)(mat).data # balanced cooler(don't know why mat+=1 failed here)
        else:
            mat = LogNorm(vmin=1, vmax=vmax, clip=True)(mat).data
    if fillna:
        # fillna with 0
        mat = np.nan_to_num(mat, nan=0)
    fig = go.Figure()
    # default_args = dict(
    #     zmax = vmax,
    # )
    # default_args.update(args)
    fig.add_trace(
        go.Heatmap(
            z = mat,
            x = columns,
            y = index,
            colorscale=fall if cmap=="fall" else cmap, # don't know why log_fall failed here
            showscale=False,
            #**default_args
            **args
        )
    )
    fig.update_layout(
        title=title
    )
    return fig
def plot_cool(coolp, title="", region="chr1",vmax=100, balance=False, ignore_diags=True, donorm=True, **args):
    """"
    Plot heatmap of single cooler file.
    Input:
        coolp: path to cooler file.
        title: name of the plot.
        region: genome region to plot.
            "chr1" or "chr1:1000000-2000000" or ["chr1:1,000,000-2,000,000", "chr2"] or
            slice(0,-1) or [slice(0,1000000), slice(1000000,2000000)]
            Note: can only cross chrom boundaries if using slice.
        vmax: max z value.
        balance: whether to load balanced cooler matrix.
        norm: whether to use lognorm.
    """
    mat = cool2mat(coolp, region, balance=balance)
    return _plot_mat(
        mat,
        title = title,
        vmax = vmax,
        ignore_diags = ignore_diags,
        donorm = donorm,
        balancing=balance,
        **args
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
def merge_track_data(clr, track_file, region):
    """
    Pick out the track data that is in the region.
    Input:
        clr: cooler obj.
        track_file: path to track file
        region: genome region to plot.
    Output:
        pd.Series; index: start of each bin, value: track value.
    """
    track_file = track_file
    track_data = pd.read_table(track_file)
    bin_s, bin_e = clr.extent(region)
    region_track_data = pd.merge(
        left = clr.bins()[bin_s:bin_e],
        right = track_data,
        how = "inner",
        on = ("chrom","start","end")
    ).copy().reset_index(drop=True)
    return region_track_data

def plot_cool_track(coolp, track_files, region, title, balance=False, **args):
    """
    Plot cooler with additional track files.
    Input:
        coolp: path to cooler file.
        track_files: dict; key: track name, value: (path to track file, column name in track file)
        region: genome region to plot.
        title: name of the plot.
        balance: whether to load balanced cooler matrix.
    """
    clr = Cooler(str(coolp))
    region_mat = clr.matrix(balance = balance).fetch(region)
    mat_index = clr.bins().fetch(region)
    region_mat = pd.DataFrame(region_mat, index=mat_index["start"], columns=mat_index["start"])
    figure = make_subplots(
        rows=len(track_files) + 1,
        cols=1,
        vertical_spacing=0.01,
        row_heights = [10] + [2] * len(track_files)
    )
    # add heatmap
    mat_fig = _plot_mat(region_mat, **args)
    figure.add_trace(
        mat_fig.data[0],
        row=1,
        col=1
    )
    # add tracks
    for i, (track_name, track_tuple) in enumerate(track_files.items(), start=2):
        track_file, colname = track_tuple
        region_track_data = merge_track_data(clr, track_file, region)
        # plotly can't handle multi-index, use start as index
        region_track_data = region_track_data.set_index("start").drop(["chrom","end"], axis=1)[colname]
        figure.add_trace(
            go.Scatter(
                x = region_track_data.index,
                y = region_track_data.values,
                name = track_name
            ),
            row = i,
            col = 1
        )
    figure.update_xaxes(
        visible=False,
        row = 1,
        col = 1
    )
    for i in range(2, len(track_files) + 2):
        figure.update_xaxes(
            visible=False,
            row = i,
            col = 1
        )
    figure.update_layout(
        paper_bgcolor = "white",
        plot_bgcolor = "white",
        height = 850,
        width = 700,
        title = title
    )
    return figure
def plot_compare_cool_track_hor(
    coolp1, coolp2, 
    IS_file1, IS_file2,
    eigs_file1, eigs_file2, 
    subplot_titles,
    region, title, balance=False, winSize=500000):
    """
    Plot cooler with additional track files
    """
    clrs = [Cooler(str(coolp1)), Cooler(str(coolp2))] 
    ISs = [pd.read_table(IS_file1), pd.read_table(IS_file2)]
    eigss = [pd.read_table(eigs_file1), pd.read_table(eigs_file2)]
    bin_s, bin_e = clrs[0].extent(region) # assume same bin
    bin_sele = clrs[0].bins()[bin_s:bin_e]
    region_ISs = [pd.merge(
            left = bin_sele,
            right = IS,
            how = "inner",
            on = ("chrom","start","end")
        ).copy().reset_index(drop=True) for IS in ISs]
    region_eigss = [pd.merge(
            left = bin_sele,
            right = eigs,
            how = "inner",
            on = ("chrom","start","end")
        ) for eigs in eigss]
    region_mats = [clr.matrix(balance = balance).fetch(region) 
                  for clr in clrs]
    figure = make_subplots(
        rows=1,
        cols=6,
        horizontal_spacing=0.01,
        column_widths = [2,2,10,2,2,10],
        subplot_titles = [None,None, subplot_titles[0],
                          None,None, subplot_titles[1]]
    )
    # add heatmap
    for i, region_mat in enumerate(region_mats):
        mat_fig = _plot_mat(region_mat)
        figure.add_trace(
            mat_fig.data[0],
            row=1,
            col=i*3 + 3
        )
    # add insulation score
    for i, region_IS in enumerate(region_ISs):
        figure.add_trace(
            go.Scatter(
                y = region_IS.index,
                x = region_IS[f"log2_insulation_score_{winSize}"],
                name = "IS",
                orientation='h',
                showlegend = True if i == 0 else False,
                mode = "lines",
                line = dict(color = "orange")
            ),
            row = 1,
            col = i*3 + 2,
        )
        # add TAD boundaries
        dat = region_IS.loc[region_IS[f"is_boundary_{winSize}"]]
        figure.add_trace(
            go.Scatter(
                y = dat.index,
                x = dat[f"log2_insulation_score_{winSize}"],
                mode = "markers",
                marker_symbol = "circle-open",
                marker = dict(color = "lightgreen"),
                orientation='h',
                name = "TAD boundaries",
                showlegend = True if i == 0 else False
            ),
            row = 1,
            col = i*3 + 2
        )
        # add boundary strengths
        figure.add_trace(
            go.Bar(
                y = dat.index,
                x = -dat[f"boundary_strength_{winSize}"],
                orientation='h',
                name = "boundary strengths",
                marker = dict(color = "purple"),
                showlegend = True if i == 0 else False
            ),
            row = 1,
            col = i*3 + 2
        )
        figure.update_yaxes(
            autorange = False,
            range = [0, len(region_IS.index)],
            row = 1,
            col = i*3 + 2
        )
    # add eigen value
    for i, region_eigs in enumerate(region_eigss):
        dat = region_eigs.copy()
        dat = dat.assign(AB = "A")
        dat.loc[dat["E1"] <0, "AB"] = "B"
        eigs_fig = px.bar(
            dat, y="start", x ="E1", color="AB",
            color_discrete_map={"A":"red","B":"blue"},
            orientation='h'
        )
        eigs_fig.update_traces(
            showlegend = True if i == 0 else False
        )
        for trace in eigs_fig.data:
            figure.add_trace(
                trace, 
                row=1,
                col=i*3 + 1
            )
    figure.update_yaxes(
        visible = False
    )
    figure.update_yaxes(
        visible=True,
        row = 1,
        col = 1
    )
    figure.update_xaxes(
        visible=False
    )
    figure.update_layout(
        paper_bgcolor = "white",
        plot_bgcolor = "white",
        height = 400,
        width = 900,
        title = title
    )
    return figure

def plot_compare_cool_pixels(coolpA, coolpB, region, outline_pixels, subplot_titles=["A","B"]):
    # extract outline pixels in the region
    clr = Cooler(str(coolpA))
    if outline_pixels is not None:
        valid_bins = clr.bins().fetch(region)[["chrom","start","end"]]
        valid_bins_reset = valid_bins.reset_index(drop=True)
        valid_pixels = pd.merge(valid_bins_reset.assign(key=1), valid_bins_reset.assign(key=1), on='key').drop('key', axis=1)
        valid_pixels.columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']
        print("Input pixels:", outline_pixels.shape[0])
        outline_pixels = pd.merge(
            left = outline_pixels,
            right = valid_pixels,
            how = "inner",
            on = ["chrom1","start1","end1","chrom2","start2","end2"]
        )
        print("Pixels in region:", outline_pixels.shape[0])
    # plot
    fig = make_subplots(
        rows=1,
        cols=2,
        subplot_titles = subplot_titles
    )
    for i, coolp in enumerate([coolpA, coolpB]):
        fig.add_trace(
            plot_cool(
                coolp,
                region = region
            ).data[0],
            row=1,
            col=i+1
        )
        if outline_pixels is not None:
            for _, row in outline_pixels.iterrows():
                start1, end1, start2, end2 = row[["start1","end1","start2","end2"]]
                if end1 - start1 != end2 - start2:
                    raise ValueError("Different region length")
                else:
                    binsize = end1 - start1
                fig.add_shape(
                    type="circle",
                    x0=start1 - 0.5 * binsize,
                    y0=start2 - 0.5 * binsize,
                    x1=end1 - 0.5 * binsize,
                    y1=end2 - 0.5 * binsize,
                    line=dict(
                        color="RoyalBlue",
                        width=1,
                    ),
                    xref=f"x{i+1}",
                    yref=f"y{i+1}"
                )
        else:
            pass
    fig.update_layout(
        margin = dict(l=0, r=10, t=30, b=0),
        height = 350,
        width = 700
    )
    return fig
def plot_tiling_compartment(coolps, eigs_files, region, title, corr=True, strip=True,
                            balance=False, cmap="RdBu", donorm=False, Enames=["E1","E1"], **args):
    """
    Plot cooler with additional track files
    """
    clrs = [Cooler(str(coolp)) for coolp in coolps]
    eigs = [pd.read_table(eigs_file) for eigs_file in eigs_files]
    bin_s, bin_e = clrs[0].extent(region)
    if corr:
        region_mat1 = compartments(
            clrs[0].matrix(balance = balance).fetch(region),
            matrixonly = True
        )
        region_mat2 = compartments(
            clrs[1].matrix(balance = balance).fetch(region),
            matrixonly = True
        )
    else:
        region_mat1 = clrs[0].matrix(balance = balance).fetch(region)
        region_mat2 = clrs[1].matrix(balance = balance).fetch(region)
    #return region_mat1, region_mat2
    region_mat = tiling_mat(
        region_mat1,
        region_mat2
    )
    pos = clrs[0].bins()[bin_s:bin_e]["start"]
    region_mat = pd.DataFrame(region_mat, index=pos, columns=pos)
    if strip:
        # strip consecutive 0s on left or right side
        for i in range(region_mat.shape[0]):
            if np.all(region_mat.iloc[i,:] == 0):
                i_s = i
            else:
                i_s = i
                break
        for i in range(region_mat.shape[0]-1, -1, -1):
            if np.all(region_mat.iloc[i,:] == 0):
                i_e = i
            else:
                i_e = i
                break
        region_mat = region_mat.iloc[i_s:i_e,i_s:i_e]
    figure = make_subplots(
        rows=2,
        cols=2,
        vertical_spacing=0.01,
        horizontal_spacing=0.01,
        row_heights = [10,2],
        column_widths = [2,10]
    )
    # add heatmap
    mat_fig = _plot_mat(
        region_mat,
        cmap=cmap,
        donorm=donorm,
        **args
        )
    figure.add_trace(
        mat_fig.data[0],
        row=1,
        col=2
    )
    # add eigen value 1
    Ename = Enames[1]
    dat = pd.merge(
        left = clrs[1].bins()[bin_s:bin_e],
        right = eigs[1],
        how = "inner",
        on = ("chrom","start","end")
    ).copy()
    dat = dat.assign(AB = "A")
    dat.loc[dat[Ename] <0, "AB"] = "B"
    eigs_fig = px.bar(
        dat, y="start", x =Ename, color="AB",
        color_discrete_map={"A":"red","B":"blue"},
        orientation='h'
    )
    eigs_fig.update_traces(
        showlegend = False
    )
    for trace in eigs_fig.data:
        figure.add_trace(trace, row=1, col=1)
    # add eigen value 2
    Ename = Enames[0]
    dat = pd.merge(
        left = clrs[0].bins()[bin_s:bin_e],
        right = eigs[0],
        how = "inner",
        on = ("chrom","start","end")
    ).copy()
    dat = dat.assign(AB = "A")
    dat.loc[dat[Ename] <0, "AB"] = "B"
    eigs_fig = px.bar(
        dat, x="start", y =Ename, color="AB",
        color_discrete_map={"A":"red","B":"blue"},
        #orientation='h'
    )
    eigs_fig.update_traces(
        showlegend = False
    )
    for trace in eigs_fig.data:
        figure.add_trace(trace, row=2,col=2)
    # layout  
    ## heatmap
    figure.update_xaxes(
        visible=False,
        row = 1,
        col = 2
    )
    figure.update_yaxes(
        visible=False,
        row = 1,
        col = 2
    )
    # eigen 2 (upper triangle)
    figure.update_xaxes(
        visible = False,
        row = 1,
        col = 1
    )
    # eigen 1 (lower triangle)
    figure.update_yaxes(
        visible = False,
        row = 2,
        col = 2
    )
    figure.update_layout(
        paper_bgcolor = "white",
        plot_bgcolor = "white",
        height = 700,
        width = 700,
        title = title
    )
    return figure
def plot_compartment(coolp, eigs_file, region, title, eigen_col="E1", strip=False, balance=False, mask_eig_na=True,
    give_mat=False, quant=0.005, minmax="min", eq_hist=False, **args):
    """
    Plot distance-normalized Hi-C correlation matrix with compartment track.
    Input:
        coolp: path to cooler file.
        eigs_file: path to eigen vector file.
        region: genome region to plot.
        title: name of the plot.
        eigen_col: which eigen vector to plot.
        strip: whether to strip the consecutive 0s in the matrix.
        balance: whether to load balanced cooler matrix.
        mask_eig_na: whether to mask the heatmap where eigen value is NA.
            Note: you need also to set fillna=False in _plot_mat to make it work.
        give_mat: whether to return the matrix.
        quant: quantile to set the bright extent.
        minmax: whether to use min or max to set the bright extent.
            work with quant.
            Note: this assumes the matrix values are symmetric around 0.
            if eq_hist is True, this will be ignored.
        eq_hist: whether to use equalized histogram.
    """
    if mask_eig_na:
        args["fillna"] = False
    clr = Cooler(str(coolp))
    eig = pd.read_table(eigs_file)
    bin_s, bin_e = clr.extent(region)
    i_s, i_e = None, None # additional index for strip

    # --- prepare heatmap --- #
    region_mat = compartments(
        clr.matrix(balance = balance).fetch(region),
        normalize= False if balance else True, # don't balance twice
        matrixonly = True
    )
    pos = clr.bins()[bin_s:bin_e]["start"]
    region_mat = pd.DataFrame(region_mat, index=pos, columns=pos)
    if strip:
        # strip consecutive 0s on left or right side
        for i in range(region_mat.shape[0]):
            if np.all(region_mat.iloc[i,:] == 0):
                i_s = i
            else:
                i_s = i
                break
        for i in range(region_mat.shape[0]-1, -1, -1):
            if np.all(region_mat.iloc[i,:] == 0):
                i_e = i
            else:
                break
        region_mat = region_mat.iloc[i_s:i_e,i_s:i_e]
    figure = make_subplots(
        rows=2,
        cols=2,
        vertical_spacing=0.01,
        horizontal_spacing=0.01,
        row_heights = [10,2],
        column_widths = [2,10]
    )
    # add eigen value 1
    dat = pd.merge(
        left = clr.bins()[bin_s:bin_e],
        right = eig,
        how = "inner",
        on = ("chrom","start","end")
    ).copy()
    dat = dat.assign(AB = "A")
    dat.loc[dat[eigen_col] <0, "AB"] = "B"
    if strip:
        dat = dat.iloc[i_s:i_e,:]
    eigs_fig = px.bar(
        dat, y="start", x =eigen_col, color="AB",
        color_discrete_map={"A":"red","B":"blue"},
        orientation='h'
    )
    eigs_fig.update_traces(
        showlegend = False
    )
    for trace in eigs_fig.data:
        figure.add_trace(trace, row=1, col=1)
    # add eigen value 2
    dat = pd.merge(
        left = clr.bins()[bin_s:bin_e],
        right = eig,
        how = "inner",
        on = ("chrom","start","end")
    ).copy()
    dat = dat.assign(AB = "A")
    dat.loc[dat[eigen_col] <0, "AB"] = "B"
    if strip:
        dat = dat.iloc[i_s:i_e,:]
    eigs_fig = px.bar(
        dat, x="start", y =eigen_col, color="AB",
        color_discrete_map={"A":"red","B":"blue"},
        #orientation='h'
    )
    eigs_fig.update_traces(
        showlegend = False
    )
    for trace in eigs_fig.data:
        figure.add_trace(trace, row=2,col=2)
    eig_dat2 = dat.copy()

    if eq_hist:
        region_mat = pd.DataFrame(
            exposure.equalize_hist(
                region_mat.values
                ),
            index = region_mat.index,
            columns = region_mat.columns
        )

    eigen_na_mask = eig_dat2.loc[dat[eigen_col].isna(),"start"].values
    if mask_eig_na:
        region_mat.loc[
            eigen_na_mask,
            :
        ] = np.nan
        region_mat.loc[
            :,
            eigen_na_mask
        ] = np.nan
    # add heatmap
    mat_fig = _plot_mat(region_mat,cmap="RdBu_r",donorm=False, **args)
    figure.add_trace(
        mat_fig.data[0],
        row=1,
        col=2
    )
    # traces
    ## heatmap
    if not eq_hist:
        dat = region_mat.values.flatten()
        dat = dat[~np.isnan(dat)]
        right = np.quantile(dat, 1-quant)
        left = np.quantile(dat, quant)
        if minmax == "min":
            bright_extent = min(abs(right), abs(left))
        elif minmax == "max":
            bright_extent = max(abs(right), abs(left))
        else:
            raise ValueError("minmax must be min or max")
        figure.update_traces(
            selector=dict(type="heatmap"),
            zmin = -bright_extent,
            zmax = bright_extent
        )
    else:
        figure.update_traces(
            selector=dict(type="heatmap"),
            zmin = 0,
            zmax = 1
        )
    ## barplot
    figure.update_traces(
        selector=dict(type="bar"),
        width=clr.binsize,
        marker=dict(
            line = dict(
                color = "purple", # looks better than black
                width = 0 # this actually tries to removes the line, 
            )
        )
    )
    # layout  
    ## heatmap
    figure.update_xaxes(
        visible=False,
        row = 1,
        col = 2
    )
    figure.update_yaxes(
        visible=False,
        row = 1,
        col = 2
    )
    # eigen 2 (upper triangle)
    figure.update_xaxes(
        visible = False,
        row = 1,
        col = 1
    )
    # eigen 1 (lower triangle)
    figure.update_yaxes(
        visible = False,
        row = 2,
        col = 2
    )
    # global layout
    figure.update_layout(
        paper_bgcolor = "white",
        plot_bgcolor = "white",
        height = 700,
        width = 700,
        title = title
    )
    if give_mat:
        return figure, region_mat
    else:
        return figure
def plot_saddle_mpl(file, title, vmin=10**(-1),vmax=10**1):
    saddle = np.load(file, allow_pickle=True)

    plt.figure(figsize=(6,6))
    norm = LogNorm(vmin=vmin, vmax=vmax)
    im = plt.imshow(
        saddle['saddledata'],
        cmap='RdBu_r',
        norm = norm
    )
    plt.xlabel("saddle category")
    plt.ylabel("saddle category")
    plt.colorbar(im, label='obs/exp', pad=0.025, shrink=0.7)
    plt.title(title)
def plot_saddle(file, title, vmin=10**(-1), vmax=10**1, **kwargs):
    saddle = np.load(file, allow_pickle=True)
    saddledata = saddle['saddledata']
    strength = saddle['saddle_strength'][10]
    # fig = go.Figure(data=go.Heatmap(
    #     z=saddledata,
    #     colorscale='RdBu_r', # 使用红蓝色彩
    #     zmin=vmin, # 设置色彩映射的最小值
    #     zmax=vmax, # 设置色彩映射的最大值
    #     colorbar=dict(title='obs/exp'), # 颜色条属性,
    #     **kwargs
    # ))
    fig = px.imshow(
        saddledata,
        color_continuous_scale='RdBu_r',
        zmin=vmin,
        zmax=vmax,
        **kwargs
    )
    fig.add_annotation(
        text=f"strength={strength:.2f}",
        x=1,
        y=0.95,
        xref="paper",
        yref="paper",
        showarrow=False,
        font=dict(
            size=16,
            color="black"
        ),
        bgcolor="white",
        opacity=0.5
    )
    # 设置图表的标题和坐标轴标签
    fig.update_layout(
        title=title,
        xaxis_title="saddle category",
        yaxis_title="saddle category",
        height = 500,
        width = 500,
        coloraxis_colorbar=dict(
            title='obs/exp',
            thicknessmode='pixels',
            thickness=15,
            lenmode='pixels',
            len=200,
            yanchor='bottom',
            y=0.5
        )
    )
    
    return fig
# --- diagonal-track plot ---
def plot_IS(IS_file):
    pass

# --- mat pileup plots ---
def pileup_IS(coolp, IS, genome="hg19", flank=800_000, **args):
    """
    Pileup TAD boundaries from a cool file.
    Input:
        coolp: cooler path(adding resolutions when mcool)
        IS: Insulation score; pd.Dataframe
        genome: reference genome. eg. hg19
        flank: flank size
    Output:
        pd.Dataframe
    """
    c = Cooler(coolp)
    # Use bioframe to fetch the genomic features from the UCSC.
    chromsizes = bioframe.fetch_chromsizes(genome)
    cens = bioframe.fetch_centromeres(genome)
    arms = bioframe.make_chromarms(chromsizes, cens)
    # this will drop the hsv and other non-standard chromosomes
    arms =  arms.set_index("chrom").loc[[i for i in c.chromnames if i.startswith("chr")]].reset_index()
    # c_expect = cooltools.expected_cis(c, view_df=arms, nproc=2, chunksize=1_000_000)
    # stack = cooltools.pileup(c, IS, view_df=arms, flank=flank,
    #                          expected_df=c_expect, expected_value_col="balanced.avg")
    stack = cooltools.pileup(c, IS, view_df=arms, flank=flank, **args)
    mtx = np.nanmean(stack, axis=2)
    mtx = pd.DataFrame(mtx).T
    return mtx
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
    end = start + steps*resolution
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

    ins_ax.set_ylim(-1.5, 1)  # Set y-axis range

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