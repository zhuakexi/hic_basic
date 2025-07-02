from itertools import cycle

from .general import template
import plotly.express as px
import plotly.graph_objects as go

def plot_bar(df, val_col="value", group_by="group", color_map=None, x_sep=" "):
    """
    Plot bar plot with point and error bar.
    Input:
        df: DataFrame with columns [group_by, val_col]
        val_col: column name for the value to plot
        group_by: column name to group by
        color_map: dict mapping group names to colors
        x_sep: separator of groupby names since there might be multiple. e.g. "patient + time"
    Output:
        fig: plotly figure object
    """
    groups = list(df.groupby(group_by).groups.keys())
    # --- prepare bar data --- #
    bar_data = df.groupby(group_by)[val_col].mean()
    bar_error = df.groupby(group_by)[val_col].std()
    if color_map is None:
        color_map = dict(zip(groups,px.colors.qualitative.G10))
    elif isinstance(color_map, list):
        # if color_map is a list, convert it to a dict
        color_map = dict(zip(groups, cycle(color_map)))
    elif isinstance(color_map, dict):
        # if color_map is a dict, make sure it has all groups
        missing_groups = set(groups) - set(color_map.keys())
        if missing_groups:
            raise ValueError(f"color_map is missing groups: {missing_groups}")
    else:
        raise ValueError("color_map must be a list or a dict")
    # --- plot --- #
    fig = go.Figure()
    for group in groups:
        if isinstance(group, tuple):
            group_name = [str(s) for s in group]
            group_name = x_sep.join(group_name)
        fig.add_trace(
            go.Bar(
                x=[group_name],
                y=[bar_data[group]],
                name = group_name,
                error_y=dict(
                    type='data',
                    array=[bar_error[group]],
                    visible = True
                    ),
                marker_color = color_map[group],
                marker_line_color = "black",
                marker_line_width = 1
        )
        )
        strip_data = df.set_index(group_by).sort_index().loc[
            group,
            val_col
            ]
        strip_fig = px.strip(
            x = [group_name] * len(strip_data),
            y = strip_data,
        )
        strip_fig.update_traces(
            marker = dict(
                color = "black",
                size = 6,
                symbol = "circle-open"
            )
        )
        for trace in strip_fig.data:
            fig.add_trace(trace)
    fig.update_layout(
        template = template,
        plot_bgcolor = "white",
        height = 500,
        width = 500
    )
    return fig