import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from .utils import filling_l2r_plotly
def plot_eigenfaces(pca,rows=2,cols=2,fsize=16):
    fig = make_subplots(rows=rows,cols=cols)
    for row, col, i, feature in filling_l2r_plotly(rows, cols, pca.components_):
        fig.add_trace(
            go.Heatmap(
                z = feature.reshape(fsize,-1)
            ),
            row = row,
            col = col
        )
    fig.update_traces(
        showscale = False
    )
    fig.update_xaxes(visible=False)
    fig.update_yaxes(visible=False)
    fig.update_layout(
        height = 500,
        width = 800
    )
    return fig
def plot_elbow(pca):
    fig = px.scatter(pca.explained_variance_ratio_)
    fig.update_layout(
        height = 500,
        width = 500,
        title = "PCA elbow plot"
    )
    return fig