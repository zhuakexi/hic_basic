from itertools import repeat
from copy import deepcopy
import open3d as o3d
import numpy as np
import pandas as pd
from ..data import fetch_cent_chromlen
from .utils import space_grid, corners2edges
def parallel_light(plate, direction):
    """
    Adding same direction to all origin points.
    Input:
        plate: N * N * 3 array
        direction: (3,) 3D vector
    Output:
        N * N * 6 array
    """
    N = plate.shape[0]
    plate = plate.reshape(N**2, 3)
    stacked_direction = np.concatenate(list(repeat(direction.reshape(3,1),N**2)),axis=1).T
    ray = np.concatenate([plate, stacked_direction], axis=1)
    return ray.reshape(N, N, 6)
def say_cheese(plate, direction, scene):
    """
    Take a photo. The plate way!
    """
    N = int(plate.shape[0]**0.5)
    ray = parallel_light(plate, direction)
    ray = o3d.core.Tensor(ray,dtype=o3d.core.Dtype.Float32)
    return scene.cast_rays(ray)["t_hit"].numpy()
def primary_views(_3dg, ngrid=16, method="ray", keep_main=True):
    """
    Get primary views from three orthogonal faces of the nucleus' oriented bounding box.
    Input:
        _3dg: inner _3dg data structure (parsing from hickit output .3dg file)
        ngrid: number of grids in each dimension
        method: "ray" (depth from each primary direction) or "distance" (distance from grid point to mesh surface)
    Output:
        dict with keys (bases, name_of_vectors, primary_figures):
            bases: three major vectors (from center to three faces of the box) and center point vector; all are column arrays
            extent: extent of the box
            name_of_vectors: names of the three major vectors and "center"
            primary_figures: 3 * 2 figures (each figure is a 2D array), 2 direction for each major vector
    """
    result = {}

    # --- remove wandering fragments --- #
    # useful when sample have Y chromosome
    points = o3d.geometry.PointCloud()
    points.points = o3d.utility.Vector3dVector(_3dg.values)
    if keep_main:
        labels = np.array(points.cluster_dbscan(eps=5, min_points=40))
        clusters = pd.Series(labels).value_counts()
        result["clusters"] = clusters[:5]
        largest_cluster = clusters.index[0]
        result["bad_particles"] = np.where(labels!=largest_cluster)[0]
        points = points.select_by_index(result["bad_particles"], invert=True)
    else:
        result["bad_particles"] = []

    # --- generate oriented bounding box --- # 
    mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(points, 2)
    obb = mesh.get_oriented_bounding_box()
    # get bounding box
    corners = np.asarray(obb.get_box_points())
    result["obb_corners"] = corners
    result["obb_edges"] = corners2edges(corners, ret_xyz=True)
    # get box shape
    result["extent"] = obb.extent
    # get three major vectors and center point vector
    bases = np.concatenate([obb.R, obb.center.reshape(3,1)], axis=1) # generating bases for homogeneous coordinates 
    result["bases"] = bases
    grid = space_grid(
        bases,
        obb.extent*0.5, 
        ngrid)
    # naming vectors according to length
    axis_def = {0:"left-right",1:"dorsal-ventral",2:"head-tail"} # shortes to longest
    axis_names = list(map(lambda x:axis_def[x], np.argsort(obb.extent)))
    axis_names.append("center")
    result["name_of_vectors"] = axis_names
    # get 8 side view figures
    mesht = o3d.t.geometry.TriangleMesh.from_legacy(mesh)
    scene = o3d.t.geometry.RaycastingScene()
    _ = scene.add_triangles(mesht)
    if method == "distance":
        grid = o3d.core.Tensor(grid, dtype=o3d.core.Dtype.Float32)
        signed_distance = scene.compute_signed_distance(grid)
        signed_distance = signed_distance.numpy()
        primary_figures = [
            [signed_distance[0,:,:], signed_distance[-1,:,:]],
            [signed_distance[:,0,:], signed_distance[:,-1,:]],
            [signed_distance[:,:,0], signed_distance[:,:,-1]],
        ]
        result["primary_figures"] = primary_figures
    elif method == "ray":
        ray_directions = (obb.R * obb.extent.reshape(1,3)*0.5)
        res = map(
            say_cheese,
            [grid[0,:,:,:], grid[-1,:,:,:], grid[:,0,:,:], grid[:,-1,:,:], grid[:,:,0,:], grid[:,:,-1,:]],
            [ray_directions[:, 0], (-1)*ray_directions[:, 0], ray_directions[:, 1], (-1)*ray_directions[:, 1], ray_directions[:, 2], (-1)*ray_directions[:, 2]],
            repeat(scene)
        )
        primary_figures = []
        for i in range(3):
            tmp = []
            for j in range(2):
                tmp.append(next(res))
            primary_figures.append(tmp)
        result["primary_figures"] = primary_figures
    return result
    # obb.R * obb.extent.reshape(1,3) * 0.5 # vectors with proper lengths
def _3dg2mesh(_3dg, method="convex_hull", watertight=False, **args):
    """
    Convert 3dg data structure to mesh.
    Input:
        _3dg: _3dg data structure (parsing from hickit output .3dg file)
        method: method to generate mesh, one of ["convex_hull", "alpha_shape", "poisson"]
        watertight: whether to force watertight
    Output:
        mesh
    TODO: make poisson and alpha_shape watertight
    """
    # --- generate point cloud ---
    xyz = _3dg.values # (N, 3)
    points = o3d.geometry.PointCloud()
    points.points = o3d.utility.Vector3dVector(xyz)
    # --- build mesh from point cloud ---
    default_args = {
        "alpha_shape" : {"alpha" : 2},
        "convex_hull" : {},
        "poisson" : {"depth":9, "n_threads":8},
    }
    if method not in default_args:
        raise ValueError(f"method {method} not supported")
    args = {**default_args[method], **args}
    if method == "poisson":
        if watertight:
            raise ValueError("poisson method does not support forcing watertight")
        points.estimate_normals()
        points.orient_normals_consistent_tangent_plane(100)
        mesh, _ = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(points, **args)
        mesh = mesh.simplify_quadric_decimation(10000)
    elif method == "convex_hull":
        mesh, _ = points.compute_convex_hull()
    elif method == "alpha_shape":
        mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(points, **args)
        if watertight:
            while not mesh.is_watertight():
                if args["alpha"] > 100:
                    raise ValueError("Fail to generate watertight mesh")
                args["alpha"] += 1
                print("alpha shape not watertight, increase alpha to %d" % args["alpha"])
                mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(points, **args)
    else:
        raise ValueError(f"method {method} not supported")
    if method != "convex_hull":
        mesh.remove_degenerate_triangles()
        mesh.remove_duplicated_triangles()
        mesh.remove_duplicated_vertices()
        mesh.remove_non_manifold_edges()
    return mesh
def _keep_largest_mesh(mesh):
    """
    Keep the largest connected mesh.
    Input:
        mesh
    Output:
        mesh
    """
    triangle_clusters, cluster_n_triangles, cluster_area \
        = mesh.cluster_connected_triangles()
    triangle_clusters = np.asarray(triangle_clusters)
    cluster_n_triangles = np.asarray(cluster_n_triangles)
    cluster_area = np.asarray(cluster_area)
    # only keep cluster with most triangles
    edge_to_remove = \
        cluster_n_triangles[triangle_clusters] < max(cluster_n_triangles)
    mesh1 = deepcopy(mesh)
    mesh1.remove_triangles_by_mask(edge_to_remove)
    mesh1 = mesh1.remove_unreferenced_vertices()
    return mesh1
def calc_depth(_3dg):
    """
    Calculate the depth (from nuclear membrane) of a chromatin bin.
    Input:
        inner _3dg data structure (parsing from hickit output .3dg file)
    Output:
        same dataframe with new column 'depth'
    """
    xyz = _3dg.values # (N, 3)
    mesh = _3dg2mesh(_3dg, method="alpha_shape", watertight=False, alpha=2)
    mesh = _keep_largest_mesh(mesh) # try to remove nucleolus holes
    # transform to implicit representation
    mesh = o3d.t.geometry.TriangleMesh.from_legacy(mesh)
    scene = o3d.t.geometry.RaycastingScene()
    _ = scene.add_triangles(mesh)
    # query depth
    query_points = o3d.core.Tensor(xyz, dtype=o3d.core.Dtype.Float32)
    res = scene.compute_signed_distance(query_points)
    #eturn res
    return _3dg.assign(depth=res.numpy())
def radial_position(_3dg:pd.DataFrame):
    """
    Calculate the radial position of each chromatin bin.
    Input:
        inner _3dg data structure (parsing from hickit output .3dg file)
    Output:
        same dataframe with new column 'radial_position'
    """
    xyz = _3dg.values # (N, 3)
    ref_pos = xyz.mean(axis=0)
    mean_radius = np.sqrt(((xyz - ref_pos)**2).sum(axis=1)).mean()
    radial_positions = np.sqrt(((xyz - ref_pos)**2).sum(axis=1)) / mean_radius
    return _3dg.assign(
        radial_position=radial_positions
        )
def _3dg_volume(_3dg, method="convex_hull", **args):
    """
    Compute volume of 3dg data structure.
    Input:
        _3dg: _3dg data structure (parsing from hickit output .3dg file)
        method: method to generate mesh, one of ["convex_hull", "alpha_shape", "poisson"]
    Output:
        volume
    """
    mesh = _3dg2mesh(_3dg, method, watertight=True, **args)
    mesh = _keep_largest_mesh(mesh)
    return mesh.get_volume()
def _radius_of_gyration(xyz):
    """
    Compute radius of gyration of a set of points.
    Input:
        xyz: (N, 3) array
    Output:
        radius of gyration
    """
    xyz = xyz - xyz.mean(axis=0)
    return np.sqrt((xyz**2).sum(axis=1).mean())
def _radius_of_gyration_matrix(xyz):
    """
    Compute pairwise radius of gyration matrix of a set of points.
    Input:
        xyz: (N, 3) array
    Output:
        (N, N) array
    """
    N = xyz.shape[0]
    rg = np.zeros((N, N))
    for i in range(N):
        for j in range(i, N):
            rg[i, j] = rg[j, i] = _radius_of_gyration(xyz[i:j+1])
    return rg
def _3dg_radius_of_gyration(_3dg):
    """
    Compute radius of gyration of 3dg data structure.
    Input:
        _3dg: _3dg data structure (parsing from hickit output .3dg file)
    Output:
        radius of gyration
    """
    xyz = _3dg.values # (N, 3)
    xyz = xyz - xyz.mean(axis=0)
    return np.sqrt((xyz**2).sum(axis=1).mean())
def _3dg_chrom_radius_of_gyration(_3dg):
    """
    Compute radius of gyration of 3dg data structure for each chromosome.
    Input:
        _3dg: _3dg data structure (parsing from hickit output .3dg file)
    Output:
        n chroms * 1 array
    """
    chroms = _3dg.index.get_level_values("chr").unique()
    return pd.Series(
        list(map(
            lambda x:_3dg_radius_of_gyration(_3dg.loc[x]),
            chroms)),
        index=chroms)
def append_centelo_to_index(df, genome="mm10", dis=2e6, p=False, q=True):
    """
    Append centromere and telomere information to the index.
    Input:
        df: inner _3dg data structure (see parse_3dg in hires_utils)
        genome: genome name, use this to determine centromere and telomere positions
        dis: distance to centromere/telomere
        p: whether to include p arm
            For acrocentric chromosomes, set p to False
        q: whether to include q arm
    Output:
        dat: data structure with centromere and telomere information
    """
    dat = df.assign(centelo="Other")
    centelo = fetch_cent_chromlen(genome)
    for chrom in dat.index.get_level_values(0).unique():
        cent_start, cent_end = centelo.loc[chrom, ["start","end"]]
        # arm1 paracentric
        if p:
            dat.loc[(chrom, cent_start-dis):(chrom, cent_start), "centelo"] = "Cent1"
        else:
            pass
        if q:
            # arm2 paracentric
            dat.loc[(chrom, cent_end):(chrom, cent_end+dis), "centelo"] = "Cent2"
        else:
            pass

        start, end = 0, centelo.loc[chrom, "chrom_length"]
        # arm1 near telomere
        if p:
            dat.loc[(chrom, start):(chrom, start+dis), "centelo"] = "Telo1"
        else:
            pass
        # arm2 near telomere
        if q:
            dat.loc[(chrom, end-dis):(chrom, end), "centelo"] = "Telo2"
        else:
            pass
    dat = dat.set_index("centelo", append=True)
    return dat
def C2T_diff(chunk):
    """
    Compute the difference between centromere and telomere.
    Use this to treat various centromere and telomere circumstances.
    """
    chrom = chunk.index.get_level_values(0).unique()
    if len(chrom) > 1:
        raise ValueError("Chunk contains multiple chromosomes")
    chrom = chrom[0]
    centelo = chunk.index.get_level_values(1)
    diffs = {}
    if ("Cent1" in centelo) and ("Telo1" in centelo):
        diffs["C2T_diff1"] = chunk.loc[chrom, "Cent1"] - chunk.loc[(chrom, "Telo1")]
    if ("Cent2" in centelo) and ("Telo2" in centelo):
        diffs["C2T_diff2"] = chunk.loc[chrom, "Cent2"] - chunk.loc[(chrom, "Telo2")]
    if len(diffs) == 0:
        diffs["C2T_diff"] = pd.Series(pd.NA, index=chunk.columns)
    res = pd.concat(diffs, axis=1).T
    res.index = pd.MultiIndex.from_product(
        [[chrom], res.index]
    )
    return res
def cent2telo_vector(_3dg, genome="mm10", dis=2e6, p=False, q=True, show_chroms=False):
    """
    Compute the C2T vector.
    Input:
        _3dg: inner _3dg data structure (see parse_3dg in hires_utils)
        genome: genome name, use this to determine centromere and telomere positions
        dis: distance to centromere/telomere
        p: whether to include p arm
            For acrocentric chromosomes, set p to False
        q: whether to include q arm
        show_chroms: whether to return chromosomal C2T vectors
    Output:
        if show_chroms:
            dataframe of chromosomal C2T vectors
        else:
            normed_c2t_L: float of normalized C2T vector length
                normally in the range of 0.00x
    Usage:
        _3dg = parse_3dg(
            _3dg_path
        )
        c2t_vec = cent2telo_vector(_3dg, genome="mm10", dis=2e6, p=False, q=True)
    """
    c2t_vector = append_centelo_to_index(
        _3dg,
        genome,
        dis=dis,
        p=p, q=q
        ).groupby(# chrom, centelo
            level=[0,2], observed=True
        ).mean()
    res = []
    for group, chunk in c2t_vector.groupby(
        level = 0, observed=True
    ):
        chunk_res = C2T_diff(chunk)
        if chunk_res.isna().all().all():
            continue
        res.append(chunk_res)
    res = pd.concat(
        res,
        axis=0
    )
    if show_chroms:
        return res
    c2t_vector = res.sum()
    c2t_L = np.sqrt((c2t_vector**2).sum())
    normed_c2t_L = c2t_L / _3dg.shape[0]
    return normed_c2t_L