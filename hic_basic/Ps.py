import bioframe
import cooltools
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from cooler import Cooler

"""
# example of how to use the functions
clr = cooler.Cooler('/home/cooler/hg19/hg19.cool')
cvd_merged, der = ps_curve(clr, nproc = 1)
fig = plot_ps_curve(cvd_merged)
"""
# Ps curve related functions
def get_full_chrom_arm_view(clr, genome):
    """
    Get chromosome arm regions used in clr.
    Input:
        clr: cooler object
        genome: name of clr's ref genome
    Output:
        dataframe
    """
    chromsizes = bioframe.fetch_chromsizes(genome)
    cens = bioframe.fetch_centromeres(genome)
    arms = bioframe.core.construction.add_ucsc_name_column(bioframe.make_chromarms(chromsizes,  cens))
    
    arms = arms[arms.chrom.isin(clr.chromnames)].reset_index(drop=True)
    return arms
def ps_curve(coolp, view_df=None, all_region=False, nproc=4):
    """
    Get P(s) curve of a balanced clr.
    cooltools version 0.5.1 or higher is required.
    Input:
        coolp: cooler file path
        view_df: dataframe with columns ['chrom', 'start', 'end', 'name']
        all_region: if True, return cis-exp for all regions
            else, return aggregated Ps and derivative
        nproc: number of processes to use
    Output:
        if all_region:
            cvd_smooth_agg( ps curve ), dataframe
        else:
            cvd_merged( ps curve ), der( ps curve derivative ); list of dataframe
    """
    #arms = get_full_chrom_arm_view(clr, genome)
    clr = Cooler(coolp)
    # calculate expected
    cvd_smooth_agg = cooltools.expected_cis(
        clr = clr,
        view_df = view_df, # if None, use all chromosomes

        smooth = True, # use smoothing with expanding windows
        aggregate_smoothed = True, # add a new column with the average of regions
        nproc = nproc
    )
    # adding real distance in bp
    cvd_smooth_agg['s_bp'] = cvd_smooth_agg['dist'] * clr.binsize
    # drop first 2 diagonals
    cvd_smooth_agg['balanced.avg.smoothed.agg'].loc[cvd_smooth_agg['dist'] < 2] = np.nan
    if all_region:
        return cvd_smooth_agg
    else:
        # Just take a single value for each genomic separation
        # because all regions with the same separation have the same balanced.avg.smoothed.agg
        cvd_merged = cvd_smooth_agg.drop_duplicates(subset=['dist'])[['s_bp', 'balanced.avg.smoothed.agg']]
        # Calculate derivative in log-log space
        der = np.gradient(np.log(cvd_merged['balanced.avg.smoothed.agg']),
                    np.log(cvd_merged['s_bp']))
        return cvd_merged, der
def plot_ps_curve(cvd_smooth_agg, ps_col="balanced.avg.smoothed.agg"):
    """
    Plot P(s) curve.
    Input:
        cvd_smooth_agg: output of ps_curve
        ps_col: column name of P(s) curve
    Output:
        fig: plotly figure
    """
    fig = px.scatter(
        cvd_smooth_agg,
        x = "s_bp",
        y = ps_col,
        log_x=True,
        log_y=True,
    )
    fig.update_traces(
        mode = "markers",
        marker_color = "blue",
        marker_size = 1,
        marker_opacity = 1,
    )
    fig.update_layout(
        xaxis_title = "Distances(bp)",
        yaxis_title = "Contact Probability P(s)",
        height = 500,
        width = 600,
        title = "P(s) curve",
        plot_bgcolor = "rgba(0,0,0,0)",
    )
    return fig
def plot_ps_curve_all_chromosomes(cvd_smooth_agg, ps_col="balanced.avg.smoothed", highlight_chroms=None, subset=None):
    chroms = cvd_smooth_agg["region1"].unique()
    if subset is not None:
        chroms = [chrom for chrom in chroms if chrom in subset]
    if highlight_chroms is None: # plot all chromosomes
        highlight_chroms = chroms
        fig = go.Figure()
        fig.update_xaxes(type="log")
        fig.update_yaxes(type="log")
    else: # highlight some chromosomes, others in gray
        other_chroms = [chrom for chrom in chroms if chrom not in highlight_chroms]
        other_chroms_ps = cvd_smooth_agg.loc[cvd_smooth_agg["region1"].isin(other_chroms)]
        fig = px.scatter(
            other_chroms_ps,
            x = "s_bp",
            y = ps_col,
            log_x=True,
            log_y=True,
            )
        fig.update_traces(
            mode = "markers",
            marker_color = "gray",
            marker_size = 1,
            marker_opacity = 0.5,
        )
    for chrom in highlight_chroms:
        chrom_ps = cvd_smooth_agg.loc[cvd_smooth_agg["region1"]==chrom]
        fig.add_scatter(
            x = chrom_ps["s_bp"],
            y = chrom_ps[ps_col],
            name = chrom,
            mode = "markers",
        )
    fig.update_layout(
        xaxis_title = "Distances(bp)",
        yaxis_title = "Contact Probability P(s)",
        height = 500,
        width = 600,
        title = "P(s) curves",
        plot_bgcolor = "rgba(0,0,0,0)",
    )
    return fig