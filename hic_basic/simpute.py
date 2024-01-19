import os
import sys
import time

import cooler
import dask
import dask.dataframe as dd
import numpy as np
import pandas as pd
from dask.delayed import delayed
from hic_basic.binnify import GenomeIdeograph
from hires_utils.hires_io import parse_3dg
from scipy.sparse import triu
from scipy.spatial.distance import euclidean
from scipy.spatial import distance_matrix
from sklearn.neighbors import radius_neighbors_graph

from .coolstuff import gen_bins

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
def calculate_distance(row1, row2):
    """
    Calculate Euclidean distance between two points.
    """
    return euclidean(row1[['x', 'y', 'z']], row2[['x', 'y', 'z']])
def cis_distance_graph(_3dg_path, fo=None, max_dist=2000000, binsize=20000, n_jobs=4):
    """
    Generate distance matrix (store in bedpe-like format) from 3dg file.
    Only cis region is considered.
    Input:
        _3dg_path: str, path to 3dg file
        fo: str, path to output parquet file, if None, return dataframe
        max_dist: int, pixels with distance larger than max_dist will be discarded
        binsize: int, binsize of input 3dg file
        n_jobs: int, number of jobs to run in parallel
    Output:
        df: n * 7 dask dataframe, (chrom1, start1, end1, chrom2, start2, end2, distance)
        None if fo is not None
    """
    # --- load data ---
    structure = parse_3dg(_3dg_path)

    # Reset index and add index columns to the DataFrame
    structure = structure.reset_index()
    structure.columns = ['chrom', 'pos', 'x', 'y', 'z']
    
    # Convert to Dask DataFrame for parallel processing
    ddf = dd.from_pandas(structure, npartitions=n_jobs)

    # Define a delayed function for processing chunks
    def process_chunk(df_chunk):
        results = []
        chrom = df_chunk['chrom'].iloc[0]
        for i in range(len(df_chunk)):
            for j in range(i, len(df_chunk)):  # Only compute upper triangle
                if abs(df_chunk['pos'].iloc[i] - df_chunk['pos'].iloc[j]) <= max_dist:
                    dist = calculate_distance(df_chunk.iloc[i], df_chunk.iloc[j])
                    results.append([chrom, df_chunk['pos'].iloc[i], df_chunk['pos'].iloc[i] + binsize, 
                                    chrom, df_chunk['pos'].iloc[j], df_chunk['pos'].iloc[j] + binsize, 
                                    dist])
        return pd.DataFrame(results, columns=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'distance'])

    # Process each chromosome using Dask
    chroms = structure['chrom'].unique()
    delayed_results = []
    for chrom in chroms:
        chrom_ddf = ddf[ddf['chrom'] == chrom].compute()  # Compute each chromosome separately
        delayed_results.append(delayed(process_chunk)(chrom_ddf))

    # Combine results and save to Parquet
    ddf_result = dd.from_delayed(delayed_results)
    ddf_result = ddf_result.repartition(npartitions=4)
    if fo is None:
        return ddf_result
    else:
        dask.config.set(scheduler='single-threaded')
        print("Saving to Parquet...")
        io_start = time.time()
        ddf_result.to_parquet(fo)
        io_end = time.time()
        print(f"IO time: {io_end - io_start}")
        return None
def cis_distance_graph_df(_3dg_path, chrom=None, genome=None, fo=None, max_dist=2000000, binsize=20000, fill=False):
    """
    Generate distance matrix (store in bedpe-like format) from 3dg file.
    Only cis region is considered.
    pandas.DataFrame version, all in memory.
    Input:
        _3dg_path: str, path to 3dg file
        chrom: str, chromosome to be processed.
            if None, process all chromosomes(highly not recommended, may cause memory error)
        genome: give genome name or length path and add pixel_id to output
            if None, pixel_id will not be added
        fo: str, path to output parquet file, if None, return dataframe
        max_dist: int, pixels with distance larger than max_dist will be discarded
        binsize: int, binsize of input 3dg file
        fill: return all possible pixels (< max_dist of course), fill missing pixels with NaN
            only valid when both chrom and genome are specified
            will return dask array
    Output:
        df: n * 7 dataframe, (chrom1, start1, end1, chrom2, start2, end2, distance)
        None if fo is not None
    """
    # --- load data ---
    structure = parse_3dg(_3dg_path) # (chr, pos): x, y, z
    if genome is not None:
        ideograph = GenomeIdeograph(genome)
        # save memory
        structure = structure.reset_index()
        structure["chr"] = pd.Categorical(
            structure["chr"],
            categories=ideograph.chromosomes.index.categories,
            ordered=True
        )
        structure = structure.set_index(["chr", "pos"])
    def process_chunk(chunk, chrom=None):
        xyz = chunk[["x","y","z"]]
        if chrom is None:
            chrom = chunk.index.get_level_values(0)[0]
            start = chunk.index.get_level_values(1)
        else:
            start = chunk.index

        dist_mat = distance_matrix(xyz.values, xyz.values)
        dist_df = pd.DataFrame(
            dist_mat,
            index=pd.Index(start,name="start1"),
            columns=pd.Index(start,name="start2")
            )

        # keep triu
        mask = np.triu(np.ones(dist_df.shape, dtype=bool), k=1)
        dist_long_df = dist_df.where(mask).stack()
        dist_long_df.name = "distance"
        dist_long_df = dist_long_df.reset_index()

        # keep distance <= max_dist
        dist_long_df = dist_long_df.loc[dist_long_df['distance'] <= max_dist]
        dist_long_df['chrom1'] = chrom
        dist_long_df['chrom2'] = chrom
        dist_long_df['end1'] = dist_long_df['start1'] + binsize
        dist_long_df['end2'] = dist_long_df['start2'] + binsize
        dist_long_df = dist_long_df[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'distance']]

        if genome is not None:
            dist_long_df = ideograph.join_pixel_id(dist_long_df, binsize, intra=True)
            dist_long_df = dist_long_df[["pixel_id", "chrom1", "start1", "end1", "chrom2", "start2", "end2", "distance"]]
        else:
            pass

        return dist_long_df
    if chrom is not None:
        structure = structure.loc[chrom]
        df = process_chunk(structure, chrom)
    else:
        df = structure.groupby(level=0).apply(process_chunk).reset_index(drop=True) # chrom1, start1, end1, chrom2, start2, end2, distance
    # --- fill to all possible pixels ---
    if fill:
        if genome is None:
            raise ValueError("fill=True requires genome to be specified")
        else:
            if chrom is None: # TODO: work with all chromosomes
                raise ValueError("fill=True requires chrom to be specified")
            else:
                start_id = ideograph.pixel_id(
                    pd.Series(
                        [chrom, 0, binsize, chrom, 0, binsize],
                        index = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]
                    ),
                    binsize=binsize
                )
                last_bin = ideograph.chromosomes.loc[chrom] - ideograph.chromosomes.loc[chrom] % binsize
                end_id = ideograph.pixel_id(
                    pd.Series(
                        [chrom, last_bin, ideograph.chromosomes.loc[chrom], chrom, last_bin, ideograph.chromosomes.loc[chrom]],
                        index = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]
                    ),
                    binsize=binsize
                )
                new_df = dd.from_array(
                    np.full(end_id - start_id + 1, np.nan),
                    columns=["distance"],
                )
                new_df.index = np.arange(start_id, end_id + 1)
                new_df["distance"] = new_df["distance"].fillna(df.set_index("pixel_id")["distance"])
                new_da = new_df["distance"].to_dask_array(lengths=True)
                df = new_da

    # --- save to Parquet ---
    if fo is None:
        return df
    else:
        print("Saving to Parquet...")
        io_start = time.time()
        df.to_parquet(fo)
        io_end = time.time()
        print(f"IO time: {io_end - io_start}")
        return None