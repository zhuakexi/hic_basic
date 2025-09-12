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


### --- axes --- ###
def clip_axis_range(
    fig,
    x_lower=True,
    x_upper=True,
    y_lower=False,
    y_upper=False,
    main_fig=None,
    row=None,
    col=None
    ):
    """
    Clip plotly axis range to exactly match data.

    x- and y-axis, lower and upper bound are both free to choose.
    Iterates through data traces contained in figure, extracts min and max, and sets the respective axes range bounds to these.

    Parameters
    ----------
    fig: plotly.graph_objs._figure.Figure
        The figure to update.
    x_lower: bool
        Clip x axis lower range. Defautls to True.
    x_upper: bool
        Clip x axis upper range. Defautls to True.
    y_lower: bool
        Clip y axis lower range. Defautls to False.
    y_upper: bool
        Clip y axis upper range. Defautls to False.
    main_fig:
        If not None, change main_fig axes according to fig traces.
        Use this when you are make subplots.
    row:
        If main_fig is not None, specify the row of the subplot where you want to change axes.
    col:
        If main_fig is not None, specify the column of the subplot where you want to change axes.

    Returns
    -------
    plotly.graph_objs._figure.Figure
        Plotly graph with updated axis ranges.

    FROM: https://community.plotly.com/t/set-axis-range-to-match-data/93168/3
    """
    x_vals_lower, x_vals_upper, y_vals_lower, y_vals_upper = [], [], [], []
    
    for _trace in fig["data"]:
        if _trace["x"] is not None:  # Trace has not empty "x" data
            x_vals_lower.append(min(_trace["x"]))
            x_vals_upper.append(max(_trace["x"]))

        if _trace["y"] is not None:  # Trace has not empty "y" data
            y_vals_lower.append(min(_trace["y"]))
            y_vals_upper.append(max(_trace["y"]))

    x_vals_lower = min(x_vals_lower) if len(x_vals_lower) > 0 else None
    x_vals_upper = max(x_vals_upper) if len(x_vals_upper) > 0 else None
    y_vals_lower = min(y_vals_lower) if len(y_vals_lower) > 0 else None
    y_vals_upper = max(y_vals_upper) if len(y_vals_upper) > 0 else None

    if main_fig is None:
        fig = fig.update_xaxes(
            autorangeoptions={
                "minallowed":[None, x_vals_lower][x_lower],
                "maxallowed":[None, x_vals_upper][x_upper]
            }
        )
        
        fig = fig.update_yaxes(
            autorangeoptions={
                "minallowed":[None, y_vals_lower][y_lower],
                "maxallowed":[None, y_vals_upper][y_upper]
            }       
        )
        return fig
    else:
        main_fig = main_fig.update_xaxes(
            autorangeoptions={
                "minallowed":[None, x_vals_lower][x_lower],
                "maxallowed":[None, x_vals_upper][x_upper]
            },
            row=row,
            col=col
        )
        
        main_fig = main_fig.update_yaxes(
            autorangeoptions={
                "minallowed":[None, y_vals_lower][y_lower],
                "maxallowed":[None, y_vals_upper][y_upper]
            },
            row=row,
            col=col
        )
        return main_fig



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
def rgb_to_hex(rgb_color):
    """
    将 RGB 格式转换为十六进制颜色字符串。
    
    参数:
        r (int): 红色通道的值 (0-255)
        g (int): 绿色通道的值 (0-255)
        b (int): 蓝色通道的值 (0-255)
    
    返回:
        str: 十六进制颜色字符串，例如 "#4C78A8"
    """
    r, g, b = rgb_color
    # 确保每个通道的值在 0-255 之间，并转换为整数
    r = max(0, min(255, int(r)))
    g = max(0, min(255, int(g)))
    b = max(0, min(255, int(b)))
    
    # 使用大写字母格式化为两位十六进制字符串
    return "#{:02X}{:02X}{:02X}".format(r, g, b)
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

def interpolate_color_rgb_linear(orig, target, n, cap=1):
    """
    Interpolate between two colors in RGB space.
    Input:
        orig: original color (r, g, b)
        target: target color (r, g, b)
        n: number of colors to generate
        cap: 0<cap<1, ensure the last interpolated color is not too close to the target color
            near 0, almost same as orig
            near 1, interpolated color cover the whole range from orig to target
    Output:
        list of n interpolated colors
    """
    assert isinstance(orig, tuple) and len(orig) == 3, "orig must be a tuple of (r, g, b)"
    assert isinstance(target, tuple) and len(target) == 3, "target must be a tuple of (r, g, b)"
    assert isinstance(n, int) and n >= 0, "n must be a non negative integer"
    assert 0 < cap < 1, "cap must be between 0 and 1"
    if n == 0:
        return []
    # determine the step size for each color channel according to cap
    step = [(channel_value - orig[channel]) * cap / n for channel, channel_value in enumerate(target)]
    # generate the interpolated colors
    colors = []
    for i in range(n):
        # calculate the new color value
        new_color = tuple(
            max(0, min(255, round(orig[channel] + step[channel] * (i + 1))))
            for channel in range(3)
        )
        colors.append(new_color)
    return colors
def expand_color_sequence(colors, n, cap=0.75):
    """
    Expand color sequence by create n interpolated colors for each original color.
    The created colors are "whiter" version of the original colors.
    Group the original colors (orig_colors) in pairs, and within each pair,
    arrange the interpolated colors in opposite directions.
    For example: blue → light blue → light red → red.
    Input:
        colors: list of original colors
        n: number of interpolated colors for each original color
    Output:
        list of expanded colors
    """
    white = (255, 255, 255)
    interpolated_colors = [
        interpolate_color_rgb_linear(colors[i], white, n, cap=cap)
        for i in range(len(colors))
    ]
    expanded_color_sequence = []
    for i in range(len(colors)):
        if i % 2 == 0:
            # even index, add original color and interpolated colors
            expanded_color_sequence.append(colors[i])
            expanded_color_sequence.extend(interpolated_colors[i])
        else:
            # odd index, add reversed interpolated colors and original color
            expanded_color_sequence.extend(interpolated_colors[i][::-1])
            expanded_color_sequence.append(colors[i])
    # remove duplicates
    expanded_color_sequence = list(dict.fromkeys(expanded_color_sequence))
    return expanded_color_sequence
def expand_colors(colors, n_interpolate=1):
    """
    在相邻颜色之间插入 n_interpolate 个新颜色，使颜色过渡自然。
    
    参数:
        colors (List[tuple]): 输入颜色序列，每个颜色是 (r, g, b) 的 RGB 元组，范围 0~1。
        n_interpolate (int): 每对颜色之间插入的新颜色数量。
    
    返回:
        List[tuple]: 扩展后的颜色序列，每个颜色为 (r, g, b) 的整数元组，范围 0~255。
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
            h_deg = (h1_deg + t * dh) % 360  # 保证色相在 0~360 度之间
            h = h_deg / 360  # 归一化到 [0, 1]

            # 插值饱和度和明度
            s = s1 * (1 - t) + s2 * t
            v = v1 * (1 - t) + v2 * t

            # 降低中间点的饱和度，使过渡更柔和
            s *= (1 - t * (1 - t))

            # 限制 s 和 v 在 [0, 1] 范围内
            s = max(0.0, min(1.0, s))
            v = max(0.0, min(1.0, v))

            # 转换回 RGB
            r, g, b = colorsys.hsv_to_rgb(h, s, v)

            # 限制 RGB 分量在 [0, 1] 之间
            r = max(0.0, min(1.0, r))
            g = max(0.0, min(1.0, g))
            b = max(0.0, min(1.0, b))

            # 转换为 [0, 255] 的整数
            r_int = int(r * 255)
            g_int = int(g * 255)
            b_int = int(b * 255)

            # 确保最终值在 [0, 255] 之间
            r_int = max(0, min(255, r_int))
            g_int = max(0, min(255, g_int))
            b_int = max(0, min(255, b_int))

            expanded.append((r_int, g_int, b_int))

        # 添加当前对的第二个颜色
        expanded.append(tuple(int(c * 255) for c in c2))

    return expanded
def plotly_fig2array(fig):
    #convert a Plotly fig to  a RGB-array
    #fig_bytes = fig.to_image(format="png", height = 1600, width = 1600, scale=4)
    fig_bytes = fig.to_image(format="png", height = 2000, width = 2000, scale=4)
    buf = io.BytesIO(fig_bytes)
    img = Image.open(buf)
    return np.asarray(img)


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