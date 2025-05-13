import colorsys
import io
import os
import re

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from scipy.ndimage import rotate
from PIL import Image


### --- plot utils --- ###


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
def sub_genome_mat(mat, keep_regions):
    """
    Keep only the regions in the keep_regions list.
    Input:
        mat: index and columns are [chrom, start]
        keep_regions: list of regions to keep
            regions are tuples of (start, end)
    Output:
        mat: matrix with only the regions in keep_regions
    """
    # drop rows
    row_dropped_dfs = []
    for region in keep_regions:
        row_dropped_dfs.append(
            mat.loc[region[0]:region[1],:]
        )
    row_dropped = pd.concat(row_dropped_dfs,axis=0)
    # drop columns
    col_dropped_dfs = []
    for region in keep_regions:
        col_dropped_dfs.append(
            row_dropped.loc[:,region[0]:region[1]]
        )
    col_dropped = pd.concat(col_dropped_dfs,axis=1)
    return col_dropped
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


### --- colors --- ###


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
def hex_split(hex_color):
    """
    解析十六进制颜色代码，并返回 R、G、B 和 A 分量。

    :param hex_color: 十六进制颜色代码，可以是六位 (#RRGGBB) 或八位 (#RRGGBBAA)
    :return: 一个元组 (R, G, B, A)，其中 A 是可选的透明度分量
    """
    # 去掉颜色代码前面的 '#' 符号
    hex_color = hex_color.lstrip('#')

    # 判断颜色代码是否有 Alpha 通道
    if len(hex_color) == 8:
        # 8 位颜色代码
        r = int(hex_color[0:2], 16)
        g = int(hex_color[2:4], 16)
        b = int(hex_color[4:6], 16)
        a = int(hex_color[6:8], 16)
    elif len(hex_color) == 6:
        # 6 位颜色代码，默认 Alpha 通道为不透明
        r = int(hex_color[0:2], 16)
        g = int(hex_color[2:4], 16)
        b = int(hex_color[4:6], 16)
        a = 255  # 默认 Alpha 通道为不透明
    else:
        raise ValueError("Invalid hex color code. It must be 6 or 8 characters long.")

    return r, g, b, a

def compute_hue_diff(h1_deg, h2_deg):
    """
    计算两个色相之间的最短差值（考虑环形特性）
    """
    diff = h2_deg - h1_deg
    if diff > 180:
        diff -= 360
    elif diff < -180:
        diff += 360
    return diff

def expand_colors(colors, n_interpolate=1):
    """
    在相邻颜色之间插入 n_interpolate 个新颜色，使颜色过渡自然。
    
    参数:
        colors (List[tuple]): 输入颜色序列，每个颜色是 (r, g, b) 的 RGB 元组，范围 0~1。
        n_interpolate (int): 每对颜色之间插入的新颜色数量。
    
    返回:
        List[tuple]: 扩展后的颜色序列。
    """
    if not colors:
        return []
    
    expanded = [colors[0]]  # 保留第一个颜色

    for i in range(len(colors) - 1):
        c1 = colors[i]
        c2 = colors[i + 1]

        # 转换为 HSV 颜色空间
        h1, s1, v1 = colorsys.rgb_to_hsv(*c1)
        h2, s2, v2 = colorsys.rgb_to_hsv(*c2)

        h1_deg = h1 * 360
        h2_deg = h2 * 360
        dh = compute_hue_diff(h1_deg, h2_deg)

        for j in range(n_interpolate):
            t = (j + 1) / (n_interpolate + 1)  # 插值系数，从 0 到 1

            # 插值色相
            h_deg = h1_deg + t * dh
            h = h_deg / 360  # 回到 [0,1] 范围

            # 插值饱和度和明度
            s = s1 * (1 - t) + s2 * t
            v = v1 * (1 - t) + v2 * t

            # 降低中间点的饱和度，使过渡更柔和
            s *= (1 - t * (1 - t))

            # 转换回 RGB
            r, g, b = colorsys.hsv_to_rgb(h, s, v)
            # transform to int
            r, g, b = int(r), int(g), int(b)
            expanded.append((r, g, b))

        # 添加当前对的第二个颜色
        expanded.append(c2)

    return expanded

def plotly_fig2array(fig):
    #convert a Plotly fig to  a RGB-array
    #fig_bytes = fig.to_image(format="png", height = 1600, width = 1600, scale=4)
    fig_bytes = fig.to_image(format="png", height = 2000, width = 2000, scale=4)
    buf = io.BytesIO(fig_bytes)
    img = Image.open(buf)
    return np.asarray(img)

def plot_color_sequence(colors, output_file=None):
    """
    Plot a color sequence as a 1-row heatmap.
    
    Args:
        colors (List[str] or List[tuple]): List of colors in hex or RGB format.
        output_file (str, optional): Path to save the plot as an image. If None, the plot is shown.
    """
    # Convert hex colors to RGB if necessary
    if isinstance(colors[0], str):
        colors = [hex_to_rgb(color) for color in colors]
    
    # Create a 1-row heatmap
    fig = go.Figure(
        data=go.Heatmap(
            z=[[i for i in range(len(colors))]],  # Dummy data
            colorscale=[(i / (len(colors) - 1), f"rgb({r}, {g}, {b})")
                        for i, (r, g, b) in enumerate(colors)],
            showscale=False
        )
    )
    fig.update_layout(
        height=100,
        width=800,
        margin=dict(l=10, r=10, t=10, b=10),
        xaxis=dict(showticklabels=False),
        yaxis=dict(showticklabels=False)
    )
    
    if output_file:
        fig.write_image(output_file)
        return output_file
    else:
        return fig


### --- manuscript notebook compile --- ###
def get_fig_outprefix()->str:
    """
    Get the output prefix for the current figure.
    Change directory name and this will get the right prefix.
    Fig1a -> output/Fig.1a_
    Fig1Sb -> output/Extended_Data_Fig.1b_
    """
    tstring = os.getcwd()
    name = re.search(
        r'[^/]+$',
        tstring
    )
    if name is not None:
        name = name.group(0)
        comps = re.search(
            r"Fig(\d+)(S?)([a-zA-Z]*)",
            name
        )
        if comps is not None:
            fig = comps.group(1)
            sup = comps.group(2)
            let = comps.group(3)
        else:
            raise ValueError("Cannot parse figure name")
    else:
        raise ValueError("No figures found in path")
    if sup == "S":
        figpr = f"output/Extended_Data_Fig.{fig}{let}_"
    else:
        figpr = f"output/Fig.{fig}{let}_"
    return figpr