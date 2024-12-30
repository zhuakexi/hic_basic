# pixel-wise comparison
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from cooler import Cooler
from .coolstuff import cool2pixels
def pick_band_pixels(coolp, threshold=5_000_000, **kwargs)->pd.DataFrame:
    """
    Pick (intra) pixels within a certain band distance.
    Input:
        coolp: str, path to coolp file
        threshold: int, band distance threshold
        **kwargs: additional arguments for cool2pixels
    Output:
        pixels: pixels within the threshold
            see cool2pixels for columns, with additional band_dist
    """
    pixels = cool2pixels(coolp, **kwargs)
    pixels = pixels.query('chrom1 == chrom2')
    pixels = pixels.assign(
        band_dist = (pixels["start1"] - pixels["start2"]).abs(),
    )
    pixels = pixels.query("band_dist < @threshold")
    return pixels
def compare_pixels(this_coolp, ref_coolp, threshold=5_000_000, clr_weight_name="weight")->pd.DataFrame:
    """
    Extract and align pixels from two coolp files. 
    Input:
        this_coolp: str, path to first coolp file
        ref_coolp: str, path to second coolp file
        threshold: int, band distance threshold
            only pixels within this distance will be compared
        clr_weight_name: str, column name for clr weight
            None for raw contact values
    Output:
        concated: dataframe of merged pixels
            chrom1, start1, start2, balanced_this, balanced_ref, distance
    """
    # --- get pixels --- #
    this_pixels = pick_band_pixels(this_coolp, clr_weight_name=clr_weight_name)
    ref_pixels = pick_band_pixels(ref_coolp, clr_weight_name=clr_weight_name)
    concated = pd.merge(
        this_pixels[["chrom1","start1","start2","balanced"]].rename(columns={"balanced":"balanced_this"}),
        ref_pixels[["chrom1","start1","start2","balanced"]].rename(columns={"balanced":"balanced_ref"}),
        on=["chrom1","start1","start2"],
        how="inner"
    )
    concated = concated.query('chrom1 not in ["chrX","chrY"]')
    concated = concated.assign(
        distance = (concated["start1"] - concated["start2"]).abs(),
    )
    return concated
def plot_compare_pixels(this_coolp, ref_coolp, sample_titles=["THIS","REF"], ignore_diags=0, 
    clr_weight_name="weight", lognorm=False, sample_ratio=None, debug=False):
    """
    Compare pixels between two coolp files.
    Scatter plot version.
    Input:
        this_coolp: str, path to first coolp file
        ref_coolp: str, path to second coolp file
        sample_titles: list of str, titles for each coolp
        ignore_diags: bool or int, ignore diagonal pixels
            if True, ignore 1 bin, if int, ignore n bins
        lognorm: bool, log10 transform contact values
    Output:
        fig: plotly figure
    """
    # --- get pixels --- #
    concated = compare_pixels(this_coolp, ref_coolp, clr_weight_name=clr_weight_name)
    # --- addtionally transforms --- #
    if lognorm:
        concated["balanced_this"] = np.log1p(concated["balanced_this"])
        concated["balanced_ref"] = np.log1p(concated["balanced_ref"])

    # --- plot scatter -- #
    # check quantiles
    quantiles = concated[["balanced_this","balanced_ref"]].quantile([0.95, 0.05])
    this_q95, ref_q95 = quantiles.loc[0.95]
    # get 10%x, 90%y to put annotation text
    annote_x = concated["balanced_this"].quantile(0.1)
    annote_y = concated["balanced_ref"].quantile(0.9)
    
    # ignore diags
    if ignore_diags:
        ignore_diags = ignore_diags - 1
        binsize = Cooler(this_coolp).info["bin-size"]
        if debug:
            print(f"binsize: {binsize}", f"ignored: {ignore_diags * binsize}")
        concated = concated.query("distance > @ignore_diags * @binsize")
    if sample_ratio is not None:
        concated = concated.sample(frac=sample_ratio)
    debug_kwargs = {}
    if debug:
        print(concated.head())
        concated = concated.assign(
            log_distance = np.log10((concated["start1"] - concated["start2"]).abs())
        )
        debug_kwargs["color"] = "log_distance"
    fig = px.scatter(
        concated,
        x="balanced_this",
        y="balanced_ref",
        **debug_kwargs,
        range_x = [0, this_q95], # rm outliers
        range_y = [0, ref_q95],
        labels={
            "balanced_this":sample_titles[0],
            "balanced_ref":sample_titles[1]
            },
        #hover template
        hover_name="chrom1",
        hover_data=["start1","start2"],
        #= "{text}<extra></extra>"
    )
    if not debug:
        # use black, small dots
        fig.update_traces(
            marker=dict(
                size=2,
                color="black",
                opacity=0.05,
                #line=dict(width=1, color='DarkSlateGrey')
            )
        )
    # add pearson correlation
    fig.add_annotation(
        x=annote_x,
        y=annote_y,
        text=f"R = {concated[['balanced_this', 'balanced_ref']].corr(method='pearson').iloc[0,1]:.2f}",
        showarrow=False
    )
    fig.update_layout(
        template = h.template,
        height = 500,
        width = 500
    )
    if debug:
        return fig, concated
    return fig
def histogram2d(data, x1, x2, bins=100, range=None):
    """
    Calculate 2D histogram for two columns in a DataFrame.
    Density plot version.
    Input:
        data: pd.DataFrame, data
        x1: str, column name for x axis
        x2: str, column name for y axis
        bins: int, number of bins
        range: list of list, range for x and y axis
    Output:
        histo: np.array, 2D histogram
        x_edges: np.array, x axis edges
        y_edges: np.array, y axis edges
    """
    # clean data
    data = data.dropna(subset=[x1, x2],how="any")
    # calc hist
    hist_kwargs = {
        "bins":bins,
        "density":False
    }
    if range is not None:
        hist_kwargs["range"] = range
    histo, x_edges, y_edges = np.histogram2d(
        data[x1],
        data[x2],
        **hist_kwargs
    )
    return histo, x_edges, y_edges
def plot_compare_pixels_density(this_coolp, ref_coolp, qrange=None, sample_titles=["THIS","REF"], bins=100,
    ignore_diags=0, clr_weight_name="weight", lognorm=False, vlognorm=False, retdata=False, vrange=None, label_x=0.01, label_y=0.015,
    debug=False, **kwargs):
    """
    Compare pixels between two coolp files and plot density plot.
    Input:
        this_coolp: str, path to first coolp file
        ref_coolp: str, path to second coolp file
        qrange: list of float, quantiles for x and y axis
        sample_titles: list of str, titles for each coolp
        ignore_diags: bool or int, ignore diagonal pixels
            if True, ignore 1 bin, if int, ignore n bins
        vlognorm: bool, log10 transform contact values
        clr_weight_name: col name of weight
        lognorm: bool, log10 transform of values for histogram2d
        retdata: bool, return data for further processing,
        vrange: list of list, range for x and y axis for histogram2d,
            if both qrange and vrange are set, vrange will be used
        label_x: float, relative x position for pearson correlation annotation
        label_y: float, relative y position for pearson correlation annotation
    Output:
        if retdata:
            data: pd.DataFrame, data used for plotting
        else:
            fig: plotly figure
    """
    # --- get pixels --- #
    concated = compare_pixels(this_coolp, ref_coolp, clr_weight_name=clr_weight_name)
    
    # --- ignore diags --- #
    if ignore_diags:
        ignore_diags = ignore_diags - 1
        binsize = Cooler(this_coolp).info["bin-size"]
        concated = concated.query("distance > @ignore_diags * @binsize")
    # --- check 0.95 and 0.99 quantiles for x and y --- #
    if debug:
        quantiles = concated[["balanced_this","balanced_ref"]].quantile([0.01, 0.05, 0.95, 0.99])
        print("Data quantiles: 0.01, 0.05, 0.95, 0.99")
        print(f"{sample_titles[0]}: {quantiles.loc[0.01, 'balanced_this']:.2f}, {quantiles.loc[0.05, 'balanced_this']:.2f}, {quantiles.loc[0.95, 'balanced_this']:.2f}, {quantiles.loc[0.99, 'balanced_this']:.2f}")
        print(f"{sample_titles[1]}: {quantiles.loc[0.01, 'balanced_ref']:.2f}, {quantiles.loc[0.05, 'balanced_ref']:.2f}, {quantiles.loc[0.95, 'balanced_ref']:.2f}, {quantiles.loc[0.99, 'balanced_ref']:.2f}")
    
    # --- calc pearson correlation --- #
    pearson_r = concated[['balanced_this', 'balanced_ref']].corr(method='pearson').iloc[0,1]
    
    # --- check quantiles --- #
    if qrange is not None:
        x_range = concated["balanced_this"].quantile(qrange[0])
        y_range = concated["balanced_ref"].quantile(qrange[1])
        if vrange is None:
            vrange = [x_range, y_range]
        else:
            x_range, y_range = vrange
    else:
        if vrange is not None:
            x_range, y_range = vrange
        else:
            x_range= concated["balanced_this"].min(), concated["balanced_this"].max()
            y_range= concated["balanced_ref"].min(), concated["balanced_ref"].max()
    
    # --- calc histogram2d --- #
    histo, this_edges, ref_edges = histogram2d(
        concated, "balanced_this", "balanced_ref", bins=bins, range=vrange)

    # --- plot density -- #
    # check histo transforms
    if lognorm:
        #log10
        histo = np.log10(histo + 1)
    fig = go.Figure()
    fig.add_trace(
        go.Heatmap(
            z=histo.T, # x for this/row, y for ref/col
            x=this_edges,
            y=ref_edges,
            colorscale="blues",
            #colorscale=colorscale,
            colorbar=dict(
                title="Density"
            ),
            **kwargs
        )
    )
    # add pearson correlation
    fig.add_annotation(
        # relative to x and y axis
        x=label_x,
        y=label_y,
        text=f"R = {pearson_r:.2f}",
        showarrow=False
    )
    fig.update_layout(
        height = 500,
        width = 500,
        xaxis_title=sample_titles[0],
        xaxis_range=x_range,
        yaxis_title=sample_titles[1],
        yaxis_range=y_range,
    )
    if debug:
        return fig, concated, histo
    return fig

"""
# Example
fig, concated, histo = plot_compare_pixels_density(
    coolps["Sperm"],
    coolps["Vara2019"],
    #qrange=[[0.025, 0.975], [0.025, 0.975]],
    qrange=[[0.01, 0.99], [0.01, 0.99]],
    sample_titles=["This study", "Vara2019"],
    clr_weight_name="cis_weight",
    ignore_diags=1,
    bins=200,
    lognorm=True,
    retdata=False,
    #range=[[0,0.02],[0,0.02]],
    label_x=0.005,
    label_y=0.015,
    debug=True,
    zmax=5,
    zmin=0,
)
fig.update_layout(
    template = h.template,
    title = "",
)
fig.show(renderer="png")
"""