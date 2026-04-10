import pandas as pd
import plotly.express as px
def c_u_filter(meta, ccol, ucol, c_threshold, u_threshold):
    """
    Filter out low-quality cells.
    Input:
        ccol: column name of contact number
        ucol: column name of umis
        c_threshold: [low, high] for contacts, use None for unbounded (e.g., [None, 2000] means only upper limit)
        u_threshold: [low, high] for umis, use None for unbounded (e.g., [None, 2000] means only upper limit)
    Output:
        new df with qc-passed rows
    """
    io_filter = (~meta[ccol].isna())&(~meta[ucol].isna())&(meta[ccol]>1)&(meta[ucol]>1)
    
    # Build contact filter
    if c_threshold[0] is not None:
        c_filter_low = meta[ccol] > c_threshold[0]
    else:
        c_filter_low = True
    
    if c_threshold[1] is not None:
        c_filter_high = meta[ccol] < c_threshold[1]
    else:
        c_filter_high = True
    
    # Build umi filter
    if u_threshold[0] is not None:
        u_filter_low = meta[ucol] > u_threshold[0]
    else:
        u_filter_low = True
    
    if u_threshold[1] is not None:
        u_filter_high = meta[ucol] < u_threshold[1]
    else:
        u_filter_high = True
    
    metric_filter = c_filter_low & c_filter_high & u_filter_low & u_filter_high
    qc = meta.loc[io_filter & metric_filter]
    print("Filter out %d samples with key-NA." % (~io_filter).sum())
    print("Filter out %d samples with bad contacts or umis." % (~metric_filter).sum())
    print("Keep %d samples from %d samples" % (qc.shape[0], meta.shape[0]))
    return qc
def plot_qc(meta, ccol="contacts", ucol="umis", c_threshold=[150_000, 1_000_000], u_threshold=[10_000, 600_000],**kwargs):
    fig = px.scatter(
        meta,
        x=ccol,
        y=ucol,
        **kwargs
        )
    
    # 计算中位数
    ccol_median = meta[ccol].median()
    ucol_median = meta[ucol].median()
    
    # 添加中位数线
    fig.add_hline(y=ucol_median, line_dash="dash", line_color="red", 
                  annotation_text=f"Median: {ucol_median:.1f}", 
                  annotation_position="right")
    fig.add_vline(x=ccol_median, line_dash="dash", line_color="red",
                  annotation_text=f"Median: {ccol_median:.1f}", 
                  annotation_position="top")
    
    # 添加阈值线
    fig.add_hline(y=u_threshold[0])
    fig.add_hline(y=u_threshold[1])
    fig.add_vline(x=c_threshold[0])
    fig.add_vline(x=c_threshold[1])
    
    fig.update_layout(
        height = 700,
        width = 700
    )
    return fig
def plot_qc(meta, ccol="contacts", ucol="umis", c_threshold=[150_000, 1_000_000], u_threshold=[10_000, 600_000],**kwargs):
    fig = px.scatter(
        meta,
        x=ccol,
        y=ucol,
        **kwargs
        )
    
    # 计算中位数
    ccol_median = meta[ccol].median()
    ucol_median = meta[ucol].median()
    
    # 添加中位数线 - 使用淡蓝色
    fig.add_hline(y=ucol_median, line_dash="dash", line_color="lightblue", 
                  annotation_text=f"Median: {ucol_median:.1f}", 
                  annotation_position="top right")  # 横线上方右侧
    fig.add_vline(x=ccol_median, line_dash="dash", line_color="lightblue",
                  annotation_text=f"Median: {ccol_median:.1f}", 
                  annotation_position="bottom right")  # 竖线右侧底部
    
    # 添加阈值线
    fig.add_hline(y=u_threshold[0])
    fig.add_hline(y=u_threshold[1])
    fig.add_vline(x=c_threshold[0])
    fig.add_vline(x=c_threshold[1])
    
    fig.update_layout(
        height = 700,
        width = 700
    )
    return fig