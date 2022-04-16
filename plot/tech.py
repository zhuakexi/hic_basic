from math import ceil
from plotly.subplots import make_subplots
import plotly.graph_objects as go

from .utils  import filling_l2r_plotly
def tech_box(annote, grouping="condition", features=["defaults"], ncols=3, size=200):
    """
    DNA and RNA tech metrics.
    Input:
        annote: sample * feature dataframe
        grouping: col name to group by
        features: features to be plot. ["defaults", ...]
        ncols: number of cols of subplots
        size: subplot size
    """
    real_features = []
    for i in features:
        if i == "defaults":
            real_features.extend(["contacts","umis","rna_ratio","con_per_reads","umis_per_reads"])
        else:
            real_features.append(i)
    
    nrows = ceil(len(real_features)/ncols)
    fig = make_subplots(rows=nrows, cols=ncols, subplot_titles=features)
    for row, col, _, feature in filling_l2r_plotly(nrows, ncols, features):
        if feature is not None:
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
        height = nrows*size,
        width = ncols*size
    )
    fig.update_layout(
        showlegend=False
    )
    return fig