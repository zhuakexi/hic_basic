from itertools import repeat
import open3d as o3d
import numpy as np
from .utils import space_grid
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
def primary_views(_3dg, ngrid=16, method="distance"):
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
    TODO: adding ray casting method
    """
    result = {}
    mesh = _3dg2mesh(_3dg)
    obb = mesh.get_oriented_bounding_box()
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
def _3dg2mesh(_3dg):
    """
    Convert 3dg data structure to mesh.
    Input:
        inner _3dg data structure (parsing from hickit output .3dg file)
    Output:
        mesh
    """
    xyz = _3dg.values # (N, 3)
    # generate point cloud
    points = o3d.geometry.PointCloud()
    points.points = o3d.utility.Vector3dVector(xyz)
    # build mesh from point cloud
    mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(points, 2)
    return mesh
def calc_depth(_3dg):
    """
    Calculate the depth (from nuclear membrane) of a chromatin bin.
    Input:
        inner _3dg data structure (parsing from hickit output .3dg file)
    Output:
        same dataframe with new column 'depth'
    """
    xyz = _3dg.values # (N, 3)
    mesh = _3dg2mesh(_3dg)
    # transform to implicit representation
    mesh = o3d.t.geometry.TriangleMesh.from_legacy(mesh)
    scene = o3d.t.geometry.RaycastingScene()
    _ = scene.add_triangles(mesh)
    # query depth
    query_points = o3d.core.Tensor(xyz, dtype=o3d.core.Dtype.Float32)
    res = scene.compute_signed_distance(query_points)
    #eturn res
    return _3dg.assign(depth=res.numpy())