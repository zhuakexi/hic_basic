import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import pandas as pd

# cell cycle
## cdp scatter plot
def cdp_scatter(dat:pd.DataFrame)->go.Figure:
    fig = px.scatter(
        dat,
        x="dist", 
        y="cratio"
        )
    fig.update_layout(
        height=500,
        width=500
    )
    fig.update_xaxes(type="log")
    return fig