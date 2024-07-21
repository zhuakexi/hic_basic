from functools import partial
from itertools import product, combinations

import numpy as np

def coord2point(coords,bases):
    """
    Input:
        coord: column array, homogeneous coordinates (dim N+1 in dim N space)
        bases: column vectors (N+1 vectors for dim N space, the last col is origin point)
    Output:
        points in space (dim N), column array
    """
    N = bases.shape[0]
    assert bases.shape[1] == N + 1, f"{bases.shape[0]} x {bases.shape[1]} is not proper shape for bases, using N * (N+1)"
    assert (coords.shape[0] == N + 1) & (coords.shape[1] == 1)
    return ( bases * coords.reshape(1,N+1) ).sum(axis=1).reshape(N,1) # using broadcasting
def coords2point(coords, bases):
    """
    Input:
        coords: column arrays, homogeneous coordinates (dim N+1 in dim N space)
        bases: column vectors (N+1 vectors for dim N space, the last col is origin point)
    Output:
        points in space (dim N), column arrays
    """
    N = coords.shape[0]
    # ndarray iterrate row by row
    points = np.concatenate(
        list(map(partial(coord2point,bases=bases), map(lambda x:x.reshape(N,1), coords.T)))
        ,axis=1)
    return points
def space_grid(bases,extent=(1,1,1),num=8):
    """
    Generate (3D) spactial grid within a box defined by some bases.
    Input:
        bases: 4 vectors. 3 from box center to its three faces. 1 the center.
        extent: size of the box.
        num: number of grids to sampling in each directon.
    Output:
        dim(num, num, num, 3) ndarray, first three dims represent relative position of the grid point, 
        last dim is the real 3D coords of that grid point.
    TODO:
        1. surpport any dim
        2. set num for each dim
    """
    if isinstance(extent, np.ndarray):
        pass
    elif hasattr(extent, "__iter__"):
        extent = np.array(extent)
    else:
        raise TypeError("extent must be np.ndarray or iterable")
    query_points = np.array(list(product(*np.linspace(-1*extent, extent,num=num).T)))
    #xyz_range = np.linspace(-1*extent, extent, num=num)
    #query_points = np.stack(np.meshgrid(*xyz_range.T), axis=-1).astype(np.float32)
    grid_coords = np.concatenate(
        [query_points, np.array([1]*query_points.shape[0]).reshape(query_points.shape[0],1)],
        axis=1).T # homogeneous coordinates of grids
    grid_points = coords2point(grid_coords, bases)
    return grid_points.T.reshape(num,num,num,3)
def grid_cords(grid):
    return grid.reshape(grid.shape[0]**3, 3)
def vector_fbclose(v1, v2, atol=1e-3):
    """
    Check if two vectors are close to each other in same direction or opposite direction
    """
    return np.isclose(v1, v2, atol=atol).all() or np.isclose(v1, -v2, atol=atol).all()
def corners2edges(corners, atol=1e-3, ret_xyz=False):
    """
    Generate edges from random order of 8 corners.
    Input:
        corners: np.ndarray or pd.DataFrame of shape (8, 3)
        atol: float, tolerance for checking if two vectors are close
        ret_xyz: whether return xyz of edges. If false, return index of corners
    Output:
        edges: list of tuple of index of corners
            if ret_xyz is True, return list of tuple of xyz of edges
    """
    if isinstance(corners, np.ndarray):
        pass
    elif isinstance(corners, pd.DataFrame):
        corners = corners.values
    else:
        raise ValueError("Input should be np.ndarray or pd.DataFrame")

    candidates = {
        pair : corners[pair[0]] - corners[pair[1]]
        for pair in combinations(range(8), 2)
    }
    candidates_likey = {}
    for pair, vector in candidates.items():
        likely = [
            vector_fbclose(vector, v)
            for v in candidates.values()
        ]
        likely = sum(likely)
        candidates_likey[pair] = likely
    edges = [pair for pair, likely in candidates_likey.items() if likely == 4]
    if ret_xyz:
        return [(corners[i], corners[j]) for i, j in edges]
    return edges