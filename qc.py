import pandas as pd
import plotly.express as px
def basic_filter(meta):
    """
    Filter out low-quality cells.
        1. must be diploid
        2. 1M > contact number > 150k; 0.6M > umis > 10k
    Male: ypercent > 0.0005
    Input:
        meta: dataframe
    Output:
        new df with qc-passed rows
    """
    io_filter = (~meta["pairs_c123_num"].isna())&(~meta["umis"].isna())&(meta["pairs_c123_num"]>1)&(meta["umis"]>1)
    ploidy_filter = (meta["ypercent"] != -1) & (meta["biasedX_score"] != -1) & (meta["hap_score"] != -1) & (meta["hap_score"] < 0.5)
    metric_filter = (meta["contacts"] > 150_000) & (meta["contacts"] < 1_000_000) & (meta["umis"] < 600_000) & (meta["umis"] > 10_000)
    qc = meta.loc[io_filter & ploidy_filter & metric_filter]
    print("Filter out %d samples with key-NA." % (~io_filter).sum())
    print("Filter out %d haploid samples." % (~ploidy_filter).sum())
    print("Filter out %d samples with bad contacts or umis." % (~metric_filter).sum())
    print("Keep %d samples from %d samples" % (qc.shape[0], meta.shape[0]))
    print("  with %d male: %d female." % ((qc["ypercent"] > 0.0005).sum(), (qc["ypercent"] <= 0.0005).sum()))
    return qc
def plot_qc(meta):
    fig = px.scatter(meta,x="contacts",y="umis")
    fig.update_layout(
        height = 700,
        width = 700
    )
    fig.add_hline(10_000)
    fig.add_hline(600_000)
    fig.add_vline(150_000)
    fig.add_vline(1_000_000)
    return fig