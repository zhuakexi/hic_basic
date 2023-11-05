import numpy as np
import pandas as pd
import toolz

from hires_utils.hires_io import parse_3dg
from sklearn.neighbors import radius_neighbors_graph

def shannon_index(x):
    """
    Calculate shannon index for a row of frequency table.
    """
    x = x[x > 0]
    return -(x * np.log(x)).sum()
def intermingle(_3dg_path, min_dist=3, n_jobs=12, table=False):
    """
    Calculate intermingling metrics for each particle in 3dg file.
    Input:
        _3dg_path: str, path to 3dg file
        min_dist: int, minimum distance to be considered as neighbor
        n_jobs: int, number of jobs to run in parallel
    Output:
        intermingling_metrics: same length dataframe, index is particle index, columns are:
            chrom: str, chromosome name
            start: int, start position
            intermingling_ratio: float, ratio of intermingling
            multi_chrom_intermingling: float, shannon index of multi-chrom intermingling
            species_richness: int, number of species in the neighborhood
    """
    # --- load data ---
    structure = parse_3dg(_3dg_path)
    structure = structure.reset_index()
    structure.columns = ["chrom", "start", "x", "y", "z"]
    # --- generate radius neighbor ---
    # get sparse matrix, basically tuple of two arrays, (row, col)
    RN = radius_neighbors_graph(
        structure[["x","y","z"]].values,
        min_dist,
        mode = "distance",
        metric = "minkowski",
        p = 2,
        include_self = False,
        n_jobs=n_jobs
        ).nonzero()
    # --- count RN chromosome distribution for each particle ---
    # seq
    for_chrom_count = ((i, structure.iat[j, 0]) for i, j in zip(*RN))
    # init
    chroms = structure["chrom"].dropna().unique()
    frequency_def = dict(zip(
        chroms,
        range(len(chroms))
        ))
    init_fequency = (0,) * len(chroms)
    # binop
    def chrom_stat(acc, x):
        """
        Add value to chromosome frequency stat.
        Input:
            acc: original stat, dict, {chrom: count}, must be well initiated (has all possible chrom names as keys)
            x: input_line, (index, chromosome name), chromosome name must be in o.keys()
        Ouput:
            new_acc: updated stat
        """
        mid = list(acc)
        mid[frequency_def[x[1]]] += 1
        return tuple(mid)
    # start reduce
    per_bin_frequency = toolz.reduceby(
        lambda x: x[0], 
        chrom_stat,
        for_chrom_count,
        init_fequency
    )
    # --- generate frequency table ---
    frequency_table = pd.DataFrame(
        per_bin_frequency,
        index=chroms
    ).T
    # tidy up index
    frequency_table = pd.concat(
        [structure, frequency_table],
        axis=1,
        join="outer"
    )[frequency_table.columns]
    # --- calculate meaningful metrics ---
    # intermingling ratio
    mask = pd.DataFrame(
        [frequency_table.columns]*frequency_table.shape[0],
        columns=frequency_table.columns
        )
    mask = (mask.T == structure["chrom"]).T
    intra_near = frequency_table.where(mask, 0).sum(axis=1)
    all_near = frequency_table.sum(axis=1)
    intermingling_ratio = (all_near - intra_near) / all_near
    frequency_table_r = frequency_table.div(frequency_table.sum(axis=1, skipna=False), axis=0)
    multi_chrom_intermingling = frequency_table_r.apply(shannon_index, axis=1)
    # species richness
    species_richness = frequency_table.apply(lambda x: (x>0).sum(), axis=1)
    # prepare output
    intermingling_metrics = pd.concat(
        [
            structure[["chrom","start"]],
            intermingling_ratio.rename("intermingling_ratio"),
            multi_chrom_intermingling.rename("multi_chrom_intermingling"),
            species_richness.rename("species_richness")
        ],
        axis=1,
        join="outer"
    )
    if table:
        return intermingling_metrics, frequency_table
    else:
        return intermingling_metrics
if __name__ == "__main__":
    _3dg = "/share/home/ychi/dev/hic_basic/tests/data/intermingle.3dg"
    print(intermingle(_3dg, 0.21,table=True)[1])
    print(intermingle(_3dg, 0.21))
    print("done")