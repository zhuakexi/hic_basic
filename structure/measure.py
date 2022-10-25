import open3d as o3d
def get_bounding_box(_3dg):
    """
    Get the bounding box of the nucleus.
    Input:
        inner _3dg data structure (parsing from hickit output .3dg file)
    Output:
        three major vectors (from center to three faces of the box) and center point vector; all are column arrays
    """
    mesh = _3dg2mesh(_3dg)
    pass
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