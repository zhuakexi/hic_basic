"""
Slingshot 绘图模块 - 基于 plotly 的交互式可视化

提供以下函数：
- plot_slingshot: 主绘图函数，绘制细胞散点图（按 pseudotime 着色）+ 主曲线
- plot_clusters: 绘制细胞聚类（按细胞类型着色）
- plot_curves: 仅绘制主曲线
"""
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from typing import Optional, List, Union

from hic_basic.plot.plotly_template import apply_template, base_size, font_base_size


def _save_figure(fig, save_path, show):
    """
    保存或显示图形（仅 PNG）
    
    Parameters
    ----------
    fig : plotly.graph_objects.Figure
        Plotly 图形对象
    save_path : str
        保存路径（.png）
    show : bool
        是否显示图形
    """
    if save_path:
        if save_path.endswith('.png'):
            fig.write_image(save_path, scale=2)
        else:
            # 默认保存为 PNG
            fig.write_image(save_path + '.png', scale=2)
        print(f"Figure saved to {save_path}")
    
    if show:
        fig.show()


def plot_slingshot(
    slingshot_obj,
    title: str = "Slingshot Trajectory",
    width: int = 700,
    height: int = 600,
    point_size: int = 5,
    point_opacity: float = 0.6,
    curve_width: int = 4,
    curve_opacity: float = 0.8,
    colorscale: str = "Viridis",
    show_legend: bool = True,
    save_path: Optional[str] = None,
    show: bool = True
):
    """
    绘制 Slingshot 轨迹推断结果
    
    主绘图函数：绘制细胞散点图（按 pseudotime 着色）并叠加所有主曲线
    
    Parameters
    ----------
    slingshot_obj : Slingshot
        Slingshot 对象（需已调用 fit()）
    title : str, optional
        图形标题
    width : int, optional
        图形宽度（像素）
    height : int, optional
        图形高度（像素）
    point_size : int, optional
        散点大小
    point_opacity : float, optional
        散点透明度
    curve_width : int, optional
        曲线线宽
    curve_opacity : float, optional
        曲线透明度
    colorscale : str, optional
        颜色映射名称（plotly 支持的颜色映射）
    show_legend : bool, optional
        是否显示图例
    save_path : str, optional
        保存路径（如 "output.png"）
    show : bool, optional
        是否显示图形
    
    Returns
    -------
    fig : plotly.graph_objects.Figure
        Plotly 图形对象
    
    Examples
    --------
    >>> from hic_basic.TI.slingshot import Slingshot
    >>> from hic_basic.TI.slingshot.plot import plot_slingshot
    >>> 
    >>> # 初始化并拟合
    >>> slingshot = Slingshot(data=adata, celltype_key='celltype', start_node='SSC')
    >>> slingshot.fit(num_epochs=10)
    >>> 
    >>> # 绘图
    >>> fig = plot_slingshot(slingshot, title="Cell Trajectory")
    >>> fig.show()
    """
    # 检查是否已拟合
    if slingshot_obj.curves is None or len(slingshot_obj.curves) == 0:
        raise ValueError("Slingshot object has not been fitted yet. Call slingshot.fit() first.")
    
    # 计算统一 pseudotime
    pseudotime = _compute_unified_pseudotime(slingshot_obj)
    
    # 创建图形
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=("Cells (colored by pseudotime)", "Cells + Principal Curves"),
        specs=[[{"type": "scatter"}, {"type": "scatter"}]]
    )
    
    # 左侧：细胞按 pseudotime 着色
    scatter_left = go.Scatter(
        x=slingshot_obj.data[:, 0],
        y=slingshot_obj.data[:, 1],
        mode='markers',
        marker=dict(
            size=point_size,
            color=pseudotime,
            colorscale=colorscale,
            opacity=point_opacity,
            colorbar=dict(title="Pseudotime", x=0.45)
        ),
        name="Cells",
        hovertemplate="UMAP1: %{x:.2f}<br>UMAP2: %{y:.2f}<br>Pseudotime: %{marker.color:.3f}<extra></extra>"
    )
    fig.add_trace(scatter_left, row=1, col=1)
    
    # 右侧：细胞（灰色）+ 主曲线
    # 灰色散点
    scatter_right = go.Scatter(
        x=slingshot_obj.data[:, 0],
        y=slingshot_obj.data[:, 1],
        mode='markers',
        marker=dict(
            size=point_size,
            color='lightgray',
            opacity=point_opacity * 0.5
        ),
        name="Cells",
        showlegend=False
    )
    fig.add_trace(scatter_right, row=1, col=2)
    
    # 添加主曲线
    curve_names = _get_curve_names(slingshot_obj)
    for l_idx, curve in enumerate(slingshot_obj.curves):
        s_interp, p_interp, order = curve.unpack_params()
        curve_name = curve_names[l_idx] if l_idx < len(curve_names) else f"Lineage {l_idx}"
        
        curve_trace = go.Scatter(
            x=p_interp[order, 0],
            y=p_interp[order, 1],
            mode='lines',
            line=dict(width=curve_width),
            opacity=curve_opacity,
            name=curve_name,
            hovertemplate=f"{curve_name}<br>UMAP1: %{{x:.2f}}<br>UMAP2: %{{y:.2f}}<extra></extra>"
        )
        fig.add_trace(curve_trace, row=1, col=2)
    
    # 标记起始节点
    start_center = slingshot_obj.cluster_centres[slingshot_obj.start_node]
    start_marker = go.Scatter(
        x=[start_center[0]],
        y=[start_center[1]],
        mode='markers',
        marker=dict(
            size=15,
            color='red',
            symbol='star',
            line=dict(color='black', width=2)
        ),
        name="Start Node",
        showlegend=show_legend
    )
    fig.add_trace(start_marker, row=1, col=2)
    
    # 更新布局
    fig.update_layout(
        title=title,
        width=width,
        height=height,
        showlegend=show_legend,
        legend=dict(x=1.02, y=1, xanchor="left", yanchor="top")
    )
    
    # 应用 hic_basic 模板
    fig = apply_template(fig, title=None, width=width, height=height)
    
    fig.update_xaxes(title_text="UMAP 1", row=1, col=1)
    fig.update_yaxes(title_text="UMAP 2", row=1, col=1)
    fig.update_xaxes(title_text="UMAP 1", row=1, col=2)
    fig.update_yaxes(title_text="UMAP 2", row=1, col=2)
    
    # 保存或显示
    _save_figure(fig, save_path, show)
    
    return fig


def plot_clusters(
    slingshot_obj,
    title: str = "Slingshot Clusters",
    width: int = 600,
    height: int = 600,
    point_size: int = 6,
    point_opacity: float = 0.8,
    save_path: Optional[str] = None,
    show: bool = True
):
    """
    绘制细胞聚类（按细胞类型着色）
    
    Parameters
    ----------
    slingshot_obj : Slingshot
        Slingshot 对象
    title : str, optional
        图形标题
    width : int, optional
        图形宽度
    height : int, optional
        图形高度
    point_size : int, optional
        散点大小
    point_opacity : float, optional
        散点透明度
    save_path : str, optional
        保存路径
    show : bool, optional
        是否显示图形
    
    Returns
    -------
    fig : plotly.graph_objects.Figure
        Plotly 图形对象
    """
    # 获取细胞类型标签
    if slingshot_obj.string_labels:
        labels = [slingshot_obj.int_to_celltype.get(i, str(i)) 
                  for i in slingshot_obj.cluster_label_indices]
        label_list = list(slingshot_obj.int_to_celltype.values())
    else:
        labels = slingshot_obj.cluster_label_indices
        label_list = list(range(slingshot_obj.num_clusters))
    
    # 创建图形
    fig = go.Figure()
    
    # 为每个细胞类型添加散点
    for label in label_list:
        mask = np.array(labels) == label
        scatter = go.Scatter(
            x=slingshot_obj.data[mask, 0],
            y=slingshot_obj.data[mask, 1],
            mode='markers',
            marker=dict(size=point_size, opacity=point_opacity),
            name=str(label),
            hovertemplate=f"{label}<br>UMAP1: %{{x:.2f}}<br>UMAP2: %{{y:.2f}}<extra></extra>"
        )
        fig.add_trace(scatter)
    
    # 标记起始节点
    start_center = slingshot_obj.cluster_centres[slingshot_obj.start_node]
    start_marker = go.Scatter(
        x=[start_center[0]],
        y=[start_center[1]],
        mode='markers',
        marker=dict(
            size=15,
            color='red',
            symbol='star',
            line=dict(color='black', width=2)
        ),
        name="Start Node"
    )
    fig.add_trace(start_marker)
    
    # 更新布局
    fig.update_layout(
        title=title,
        width=width,
        height=height,
        xaxis_title="UMAP 1",
        yaxis_title="UMAP 2",
        legend=dict(x=1.02, y=1, xanchor="left", yanchor="top")
    )
    
    # 应用 hic_basic 模板
    fig = apply_template(fig, title=None, width=width, height=height)
    
    _save_figure(fig, save_path, show)
    
    return fig


def plot_curves(
    slingshot_obj,
    title: str = "Slingshot Principal Curves",
    width: int = 700,
    height: int = 600,
    curve_width: int = 4,
    curve_opacity: float = 0.9,
    show_cells: bool = True,
    save_path: Optional[str] = None,
    show: bool = True
):
    """
    绘制主曲线（可选显示细胞散点）
    
    Parameters
    ----------
    slingshot_obj : Slingshot
        Slingshot 对象
    title : str, optional
        图形标题
    width : int, optional
        图形宽度
    height : int, optional
        图形高度
    curve_width : int, optional
        曲线线宽
    curve_opacity : float, optional
        曲线透明度
    show_cells : bool, optional
        是否显示细胞散点
    save_path : str, optional
        保存路径
    show : bool, optional
        是否显示图形
    
    Returns
    -------
    fig : plotly.graph_objects.Figure
        Plotly 图形对象
    """
    if slingshot_obj.curves is None or len(slingshot_obj.curves) == 0:
        raise ValueError("Slingshot object has not been fitted yet. Call slingshot.fit() first.")
    
    fig = go.Figure()
    
    # 添加细胞散点（灰色背景）
    if show_cells:
        scatter = go.Scatter(
            x=slingshot_obj.data[:, 0],
            y=slingshot_obj.data[:, 1],
            mode='markers',
            marker=dict(
                size=4,
                color='lightgray',
                opacity=0.4
            ),
            name="Cells",
            showlegend=False
        )
        fig.add_trace(scatter)
    
    # 添加主曲线
    curve_names = _get_curve_names(slingshot_obj)
    for l_idx, curve in enumerate(slingshot_obj.curves):
        s_interp, p_interp, order = curve.unpack_params()
        curve_name = curve_names[l_idx] if l_idx < len(curve_names) else f"Lineage {l_idx}"
        
        curve_trace = go.Scatter(
            x=p_interp[order, 0],
            y=p_interp[order, 1],
            mode='lines',
            line=dict(width=curve_width),
            opacity=curve_opacity,
            name=curve_name
        )
        fig.add_trace(curve_trace)
    
    # 标记起始节点
    start_center = slingshot_obj.cluster_centres[slingshot_obj.start_node]
    start_marker = go.Scatter(
        x=[start_center[0]],
        y=[start_center[1]],
        mode='markers',
        marker=dict(
            size=15,
            color='red',
            symbol='star',
            line=dict(color='black', width=2)
        ),
        name="Start Node"
    )
    fig.add_trace(start_marker)
    
    # 更新布局
    fig.update_layout(
        title=title,
        width=width,
        height=height,
        xaxis_title="UMAP 1",
        yaxis_title="UMAP 2",
        legend=dict(x=1.02, y=1, xanchor="left", yanchor="top")
    )
    
    # 应用 hic_basic 模板
    fig = apply_template(fig, title=None, width=width, height=height)
    
    _save_figure(fig, save_path, show)
    
    return fig


def _compute_unified_pseudotime(slingshot_obj) -> np.ndarray:
    """
    计算统一的 pseudotime（所有细胞）
    
    Parameters
    ----------
    slingshot_obj : Slingshot
        Slingshot 对象
    
    Returns
    -------
    pseudotime : np.ndarray
        统一 pseudotime 数组
    """
    pseudotime = np.zeros(slingshot_obj.data.shape[0])
    
    for l_idx, lineage in enumerate(slingshot_obj.lineages):
        curve = slingshot_obj.curves[l_idx]
        cell_mask = np.logical_or.reduce(
            np.array([slingshot_obj.cluster_label_indices == k for k in lineage])
        )
        if hasattr(curve, 'pseudotimes_interp'):
            pseudotime[cell_mask] = curve.pseudotimes_interp[cell_mask]
    
    return pseudotime


def _get_curve_names(slingshot_obj) -> List[str]:
    """
    获取曲线名称列表
    
    Parameters
    ----------
    slingshot_obj : Slingshot
        Slingshot 对象
    
    Returns
    -------
    names : List[str]
        曲线名称列表
    """
    names = []
    for lineage in slingshot_obj.lineages:
        if slingshot_obj.string_labels:
            # 使用细胞类型名称
            cluster_names = [slingshot_obj.int_to_celltype.get(c, str(c)) for c in lineage.clusters]
            names.append(" → ".join(cluster_names[:4]) + ("..." if len(cluster_names) > 4 else ""))
        else:
            # 使用整数索引
            names.append(f"Lineage {lineage.clusters}")
    return names
