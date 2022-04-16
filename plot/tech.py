from plotly.subplots import make_subplots
import plotly.graph_objects as go

from .utils  import filling_l2r_plotly
def tech_box(annote, grouping="condition", size=200):
    """
    DNA and RNA tech metrics.
    Input:
        annote: sample * feature dataframe
        grouping: col name to group by
        size: subplot size
    """
    features = ["contacts","umis","rna_ratio","con_per_reads","umis_per_reads"]
    fig = make_subplots(rows=2,cols=3,subplot_titles=features)
    for row, col, _, feature in filling_l2r_plotly(2,3,features):
        if isinstance(feature, str):
            fig.add_trace(
                go.Box(
                    x = annote["condition"],
                    y = annote[feature],
                    name = feature
                ),
                row = row,
                col = col
            )
    fig.update_layout(
        height = 2*size,
        width = 3*size
    )
    fig.update_layout(
        showlegend=False
    )
    return fig