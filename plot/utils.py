import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
def filling_l2r_plotly(rows, cols, features):
    """
    Helper to iterate within row-cols. 
    Plotly_flavor, starts with 1.
    """
    for i in range(rows):
        for j in range(cols):
            k = i * cols + j
            try:
                feature = features[k]
            except IndexError:
                feature = None
            yield i+1, j+1, k, feature
def filling_l2r_mpl(rows, cols, features):
    """
    Helper to iterate within row-cols. 
    Mpl_flavor, starts with 0.
    """
    for i in range(rows):
        for j in range(cols):
            k = i * cols + j
            try:
                feature = features[k]
            except IndexError:
                feature = None
            yield i, j, k, feature
def list2colorlist(celltypes):
    """
    Tranform a data list to a color list, ready for all color argument.
    TODO:
        custom mapper;
        custom color sequence
    """
    # prepare mapper
    color_discrete_sequence = px.colors.qualitative.Plotly
    cat = tuple(set(celltypes))
    mapper = pd.Series(
        index = cat,
        data = [color_discrete_sequence[i % len(color_discrete_sequence)] for i, key in enumerate(cat)],
        name = "color"
    )
    color = [mapper[i] for i in celltypes]
    return color
def add_cat_marker(fig, data, catcol="cell_type", ypos=1):
    """
    Adding colored scatter trace to figure to mark x categories.
    Input:
        fig: plotly figure
        data: dataframe storing category info, must have same x index with fig data
        catcol: category column name, must be in adata.obs
        ypos: y position of marker points
    Output:
        fig with new traces, each cat one trace
    """
    data = data.copy()
    data = data.assign(color = list2colorlist(data[catcol]))
    for name, dat in data.groupby(catcol):
        fig.add_trace(
            go.Scatter(
                x = dat.index,
                y = [ypos for i in dat.index],
                mode = "markers",
                marker = dict(
                    color = dat["color"],
                    size = 3
                ),
                legendgroup = catcol,
                name = name
            )
    )
def tiling_mat(A, Ref, adjust = True):
    """
    Tiling 2 matrix. Rising light to the lighter one.
    Input:
        A: matrix to show on upper right.
        Ref: matrix to show on lower left.
        adjust: adjust the color scale to the stronger one. Note: maybe overflow if brightness of the two matrix are too different.
    Output:
        A tiled matrix.
    """
    #m = np.tril(Ref/lighter) + np.triu(A)
    As,Rs = A.sum(),Ref.sum()
    if adjust:
        if As < Rs:
            m = np.tril(Ref) + np.triu(A*(Rs/As)) # value overflow without parentheses in (Rs/As)
        if As > Rs:
            m = np.tril(Ref*(As/Rs)) + np.triu(A)
    else:
        m = np.tril(Ref) + np.triu(A)
    return m
def pcolormesh_45deg(ax, matrix_c, start=0, resolution=1, *args, **kwargs):
    """
    Helper to plot a matrix with 45 degree angle.
    """
    start_pos_vector = [start+resolution*i for i in range(len(matrix_c)+1)]
    import itertools
    n = matrix_c.shape[0]
    t = np.array([[1, 0.5], [-1, 0.5]])
    matrix_a = np.dot(np.array([(i[1], i[0])
                                for i in itertools.product(start_pos_vector[::-1],
                                                           start_pos_vector)]), t)
    x = matrix_a[:, 1].reshape(n + 1, n + 1)
    y = matrix_a[:, 0].reshape(n + 1, n + 1)
    im = ax.pcolormesh(x, y, np.flipud(matrix_c), *args, **kwargs)
    im.set_rasterized(True)
    return im
def mat_coarsen(mat, coarseness):
    """
    # coarsen matrix to lower resolution
    # Input:
    #    mat: matrix
    # Output:
    #    matrix
    """
    shape = np.array(mat.shape, dtype=float)
    new_shape = coarseness * np.ceil(shape/ coarseness).astype(int)
    # zero-padded array
    zp_mat = np.zeros(new_shape)
    zp_mat[:mat.shape[0],:mat.shape[1]] = mat
    temp = zp_mat.reshape(
        (   zp_mat.shape[0] // coarseness,
            coarseness,
            zp_mat.shape[1] // coarseness,
            coarseness
        )
    )
    zp_mat_c = np.sum(temp,axis=(1,3))
    return zp_mat_c