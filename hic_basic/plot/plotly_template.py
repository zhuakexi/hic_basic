"""
hic_basic Plotly 模板 - 统一可视化风格

基于 spg 项目模板的 hic_basic 版本
"""
import plotly.graph_objects as go

# 柔和、干净的配色方案
soft_clean_palette = [
    "#5B8DB8",
    "#E3A06F",
    "#8FB7A1",
    "#C08CB4",
    "#7CA7B5",
    "#C7B07B",
    "#A7B7D8",
    "#E2B6A7",
    "#B6CBB5",
    "#D5B7CC",
]

# 基础尺寸
base_size = 500
font_base_size = 12
line_base_size = 2
ticklen = 6

# 颜色设置
title_color = 'rgba(0,0,0,0)'  # 透明标题
font_color = "black"
line_color = "black"

# 创建模板
template = go.layout.Template()

# 字体
template.layout.font = dict(
    size=font_base_size,
    color=font_color
)

# 配色方案
template.layout.colorway = soft_clean_palette

# 图形尺寸
template.layout.height = base_size + 100
template.layout.width = base_size + 100

# 标题
template.layout.title = dict(
    text="TITLE",
    font=dict(
        size=font_base_size,
        color=title_color
    ),
    xanchor="center",
    x=0.5
)

# 背景
template.layout.plot_bgcolor = 'rgba(0,0,0,0)'
template.layout.paper_bgcolor = 'rgba(0,0,0,0)'

# X 轴
template.layout.xaxis = dict(
    showline=True,
    linecolor=line_color,
    linewidth=line_base_size,
    tickwidth=line_base_size,
    ticklen=ticklen,
    showgrid=False,
    zeroline=False,
    ticks="outside"
)

# Y 轴
template.layout.yaxis = dict(
    showline=True,
    linecolor=line_color,
    linewidth=line_base_size,
    tickwidth=line_base_size,
    ticklen=ticklen,
    showgrid=False,
    zeroline=False,
    ticks="outside"
)

# 颜色条
template.layout.coloraxis = dict(
    colorbar=dict(
        len=0.5,
        yanchor="bottom",
        y=0.5,
        xanchor="left",
        x=1.0,
        showticklabels=True,
        thickness=20,
    )
)


def apply_template(fig, title=None, width=None, height=None):
    """
    应用 hic_basic 模板到图形
    
    Parameters
    ----------
    fig : plotly.graph_objects.Figure
        Plotly 图形对象
    title : str, optional
        标题（可选覆盖）
    width : int, optional
        宽度（可选覆盖）
    height : int, optional
        高度（可选覆盖）
    
    Returns
    -------
    fig : plotly.graph_objects.Figure
        应用模板后的图形
    """
    # 获取模板的 JSON
    template_json = template.to_plotly_json()
    
    # 更新布局
    fig.update_layout(
        template=template_json,
        title=dict(
            text=title if title else "TITLE",
            font=dict(size=font_base_size, color=title_color),
            xanchor="center",
            x=0.5
        )
    )
    
    if width:
        fig.update_layout(width=width)
    
    if height:
        fig.update_layout(height=height)
    
    return fig
