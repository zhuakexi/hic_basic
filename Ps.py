from itertools import combinations

import numpy as np
import pandas as pd

import bioframe

import cooler
import cooltools
"""
# example of how to use the functions
clr = cooler.Cooler('/home/cooler/hg19/hg19.cool')
cvd_merged, der = ps_curve(clr, 'hg19', nproc = 1)
_plot_ps_curve(cvd_merged, der, (20e3, 1e8))
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
def ps_curve(clr, genome, nproc):
    """
    Get P(s) curve of a balanced clr.
    cooltools version 0.5.1 or higher is required.
    Input:
        clr: cooler object
        genome: name of clr's ref genome
    Output:
        cvd_merged( ps curve ), der( ps curve derivative ); list of dataframe
    """
    arms = get_full_chrom_arm_view(clr, genome)
    cvd_smooth_agg = cooltools.expected_cis(
        clr = clr,
        view_df = arms,
        smooth = True,
        aggregate_smoothed = True,
        nproc = nproc
    )
    cvd_smooth_agg['s_bp'] = cvd_smooth_agg['dist']* clr.binsize
    cvd_smooth_agg['balanced.avg.smoothed.agg'].loc[cvd_smooth_agg['dist'] < 2] = np.nan
    # Just take a single value for each genomic separation
    cvd_merged = cvd_smooth_agg.drop_duplicates(subset=['dist'])[['s_bp', 'balanced.avg.smoothed.agg']]
    # Calculate derivative in log-log space
    der = np.gradient(np.log(cvd_merged['balanced.avg.smoothed.agg']),
                  np.log(cvd_merged['s_bp']))
    return cvd_merged, der