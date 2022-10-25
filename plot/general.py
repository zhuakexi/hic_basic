# --- part 1 Set/Venn ---
from upsetplot import from_contents, UpSet
import numpy as np
import pandas as pd

import plotly.express as px
import plotly.graph_objects as go
import statsmodels.api as sm

def plot_upsetplot(sets, show_counts):
    """
    Upset plot like in R.
    Input:
        sets: dict of set
    Output:
        upset figure
    """
    return UpSet(from_contents(sets), subset_size="count", show_counts=show_counts)
def upset_plot_getter(data, keys):
    """
    Get aggregate data from upsetplot internal data format (boolen multiindex).
    Input:
        data: dict of set
        keys: set names to intersect
    Output:
        (count, set)
    """
    data = from_contents(data)
    data = data.sort_index()
    index = np.array([True for i in data.index])
    if len(keys) == 0:
        return None
    elif len(keys) == 1:
        index = [False for i in data.index.names]
        index[data.index.names.index(keys[0])] = True
        tset = data.loc[tuple(index)].values
        count = len(tset)
    else:
        for key in keys:
            index = index & data.index.get_level_values(key).values
        tset = data.loc[index].values
        count = len(tset)
        #count = np.sum(index)
    return count, tset

# --- part 2 scatter ---
def scatter_cols(data, points=None, trends=[]):
    """
    Multi-traces scatter plot.
    Input:
        data: x as index, traces as cols
        trend: plot lowess trends
        point: plot scatter points
    Return:
        go.Figure
    """
    if points is None:
        # plot point scatter for all cols by default.
        points = data.columns
    fig = go.Figure()
    for col in points:
        fig.add_trace(
            go.Scatter(
                x = data.index,
                y = data[col],
                name = col
            )
        )
    for col in trends:
        fig.add_trace(
            go.Scatter(
                x = data.index,
                y = sm.nonparametric.lowess(
                    exog = list(range(data.shape[0])),
                    endog = data[col],
                    frac = 0.2
                )[:, 1],
                name = col + "_trend"
            )
        )
    fig.update_layout(
        height = 500,
        width = 800
    )
    return fig

# --- part N 3D ---
def plot_points(array,**args):
    """
    Wrapping points (represent in numpy column arrays) to dataframe(treate x, y, z as features so in shape N * 3) and plot.
    Input:
        array: column array
    Output:
        dataframe, x, y, z as columns
    """
    data = array.T
    data = pd.DataFrame(data, columns="x y z".split())
    fig = px.scatter_3d(data, x="x",y="y",z="z",**args)
    fig.update_layout(
        height = 500,
        width = 500
    )
    return fig