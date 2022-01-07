from python_tsp.heuristics.simulated_annealing import solve_tsp_simulated_annealing
from python_tsp.distances import euclidean_distance_matrix
from concurrent import futures
import pandas as pd
from .cycle_phasing import dis_counts
def ra_ordering(filesp, threads=24):
    """
    Random Annealing ordering of contact-decay-profiles.
    Input:
        filesp : DataFrame, must have pairs col
    Output:
        DataFrame with additional col "order_index"
    """
    with futures.ProcessPoolExecutor(threads) as pool:
        res = pool.map(dis_counts,filesp["pairs"])
    ares = list(res)
    cdps = pd.DataFrame(ares)
    cdps.columns = cdps.columns.astype("string")
    cdps.index = filesp.index
    dm = euclidean_distance_matrix(cdps.values)
    # return two ele list, [order,distance]
    order = solve_tsp_simulated_annealing(dm)
    new_filesp = filesp.assign(order_index=order[0])
    return new_filesp