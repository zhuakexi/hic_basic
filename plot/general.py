# --- part 1 Set/Venn ---
from upsetplot import from_contents, UpSet
import numpy as np

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
def scatter_cols(data, trends=False, points=True):
    """
    Multi-traces scatter plot.
    Input:
        data: x as index, traces as cols
        trend: plot lowess trends
        point: plot scatter points
    Return:
        go.Figure
    """
    fig = go.Figure()
    if points:
        for col in data:
            fig.add_trace(
                go.Scatter(
                    x = data.index,
                    y = data[col],
                    name = col
                )
            )
    if trends:
        for col in data:
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