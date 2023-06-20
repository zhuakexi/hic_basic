import os
import sys
from .coolstuff import gen_bins
from .hicio import parse_3dg
import numpy as np
import pandas as pd
import cooler
from scipy.sparse import triu
from sklearn.neighbors import radius_neighbors_graph
def boolean_radius_neighbor(df, min_dist=3, n_jobs=4):
    """
    generate boolean matrix from dataframe
    Input:
        df: dataframe, index is position, columns are x, y, z
        min_dist: int, minimum distance to be considered as neighbor
        n_jobs: int, number of jobs to run in parallel
    Output:
        graph: boolean matrix, True if two points are neighbors
    """
    graph = radius_neighbors_graph(
        df.values,
        min_dist,
        mode = "distance",
        metric = "minkowski",
        p = 2,
        include_self = False,
        n_jobs=n_jobs)
    return triu(graph)
def cis_proximity_graph(_3dg_path, fo, min_dist=3, genome="mm10", binsize=20000, n_jobs=4):
    """
    Genereate 0,1 cooler file of region proximity from 3dg file.
    Only cis region is considered.
    Input:
        _3dg_path: str, path to 3dg file
        fo: str, path to output cooler file
        min_dist: int, minimum distance to be considered as neighbor
        genome: str, genome version
        binsize: int, binsize
    Output:
        None
    """
    # --- load data ---
    structure = parse_3dg(_3dg_path)
    bins = gen_bins(genome, binsize)
    # --- lign bins and structure ---
    structure = structure.reset_index()
    structure.columns = ["chrom", "start", "x", "y", "z"]
    index_structure = pd.merge(
        bins, structure.reset_index(),
        left_on=["chrom","start"], right_on=["chrom","start"], how="left").dropna()[["chrom","x","y","z"]]
    # --- prepare pixels for cooler ---
    chunks = (chunk[["x","y","z"]] for _, chunk in index_structure.groupby("chrom"))
    local_pixels = (
        (
            boolean_radius_neighbor(chunk, min_dist=min_dist, n_jobs=n_jobs).nonzero(),
            pd.Series(chunk.index, index=list(range(chunk.shape[0]))) 
            )
        for chunk in chunks
        )
    pixels = (
        {
            "bin1_id" : lifter[idx[0]].values,
            "bin2_id" : lifter[idx[1]].values,
            "count" : np.ones(idx[0].shape[0])
        } for idx, lifter in local_pixels
    )
    # --- generate cooler ---
    cooler.create_cooler(
        fo,
        bins,
        pixels,
        dtypes={"count": int},
        symmetric_upper = True,
        assembly = genome,
        ordered = False, # groupby doesn't ensure order
        triucheck = True,
        dupcheck = True,
        boundscheck = True,
        ensure_sorted = True,
    )