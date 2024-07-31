import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
def filling_l2r_plotly(rows, cols, features):
    """
    Helper to iterate within row-cols.
    Plotly_flavor, starts with 1.
    Return:
        row, col, index, feature
    """
    features_iter = iter(features)
    for i in range(rows):
        for j in range(cols):
            k = i * cols + j
            try:
                feature = next(features_iter)
            except StopIteration:
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
def tiling_mat(A_orig, Ref_orig, adjust = True, ignore_diags = True):
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
    A, Ref = A_orig.copy(), Ref_orig.copy()
    if ignore_diags:
        np.fill_diagonal(A,0)
        np.fill_diagonal(Ref,0)
        As, Rs = np.nansum(A),np.nansum(Ref)
    else:
        As, Rs = np.nansum(A),np.nansum(Ref)
    if adjust:
        if As < Rs:
            m = np.tril(Ref) + np.triu(A*(Rs/As)) # value overflow without parentheses in (Rs/As)
        elif As > Rs:
            m = np.tril(Ref*(As/Rs)) + np.triu(A)
        else:
            print("Warning: two matrix have same brightness, this seldom happens.")
            m = np.tril(Ref) + np.triu(A)
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
def spread_text(text_list, track_width=5, fold=2):
    """
    Spread text list to a several tracks.
    Input:
        text_list: list of text to spread
        track_width: width of each track
        fold: number of tracks
    Output:
        list of position shift for each text
    TODO:
        add random shift
    """
    shift_dict = dict(zip(
        list(range(fold)),
        [i*track_width for i in range(fold)]
    ))
    return [shift_dict[i%fold] for i in range(len(text_list))]
def hex_to_rgb(hex_color):
    """
    将十六进制颜色字符串转换为 RGB 格式。
    
    参数:
        hex_color (str): 十六进制颜色字符串，例如 "#4C78A8"
    
    返回:
        tuple: 三元组 (r, g, b)，其中 r, g, b 分别为红、绿、蓝通道的值 (0-255)
    """
    # 移除字符串前的 '#' 然后分解成三个长度为 2 的字符串
    hex_color = hex_color.lstrip('#')
    # 每两个字符代表一个颜色通道的十六进制值
    return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))
import io
from PIL import Image
from scipy.ndimage import rotate
def plotly_fig2array(fig):
    #convert a Plotly fig to  a RGB-array
    #fig_bytes = fig.to_image(format="png", height = 1600, width = 1600, scale=4)
    fig_bytes = fig.to_image(format="png", height = 2000, width = 2000, scale=4)
    buf = io.BytesIO(fig_bytes)
    img = Image.open(buf)
    return np.asarray(img)