import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.stats import norm
from .utils import hex_to_rgb

### --- line plot with error band --- ###

def smooth_with_ci(df, x, y, frac=0.3, ci=0.95):
    """
    Apply LOESS smoothing to time-series data and estimate confidence intervals.

    Parameters
    ----------
    x : array-like
        Independent variable (e.g., time points).
    y : array-like
        Dependent variable (e.g., gene expression values).
    frac : float, optional
        Fraction of data used for each local regression (0 < frac <= 1).
    ci : float, optional
        Confidence interval level (0 < ci < 1, e.g., 0.95 for 95% CI).

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns:
        - 'x': Original x values
        - 'y': Original y values
        - 'y_smooth': Smoothed y values
        - 'y_upper': Upper bound of CI
        - 'y_lower': Lower bound of CI
    """
    assert x in df.columns and y in df.columns, \
        f"Columns '{x}' and '{y}' must exist in the DataFrame."
    # Validate CI
    if ci > 1:
        ci = ci / 100.0  # e.g., 95 -> 0.95
    if not (0 < ci < 1):
        raise ValueError("ci must be between 0 and 1.")

    df = df.sort_values(by=x)
    x, y = df[x].values, df[y].values
    x = np.array(x, dtype=float)
    y = np.array(y, dtype=float)

    # Apply LOESS smoothing
    smoothed = lowess(y, x, frac=frac, return_sorted=False)

    # Residuals and std deviation
    residuals = y - smoothed
    sigma = np.nanstd(residuals, ddof=1)  # avoid division bias

    # z-score for CI
    z = norm.isf((1 - ci) / 2)

    upper = smoothed + z * sigma
    lower = smoothed - z * sigma

    return pd.DataFrame({
        'x': x,
        'y': y,
        'y_smooth': smoothed,
        'y_upper': upper,
        'y_lower': lower
    })
def rgb_grayer(color, white_ratio=0.7, opacity=0.3):
    """
    Convert a color to a lighter, more transparent version by blending with white
    and reducing opacity.

    Parameters:
        color: A list of three integers representing RGB values (0-255).

    Returns:
        A string in "rgba(r, g, b, a)" format with adjusted color and opacity.
    """
    # Blend the original color with white to make it lighter
    # WHITE_RATIO = 0.7 means 70% white, 30% original color
    # WHITE_RATIO = 0.7
    # Reduce opacity for a more transparent effect (0.3 = 30% opacity)
    # OPACITY = 0.3

    # Calculate new color by linear interpolation toward white
    new_color = [
        int(c * (1 - white_ratio) + 255 * white_ratio)  # 255 is white
        for c in color
    ]

    return "rgba({},{},{},{})".format(*new_color, opacity)
def add_error_band(fig, df, orig_color=px.colors.qualitative.Plotly[0], debug=False):
    """
    Input:
        df: must have x, upper, lower cols
        orig_color: color of original line, will use a grayer version
    """
    df = df.loc[df["x"].notnull() & df["y_upper"].notnull() & df["y_lower"].notnull()]
    fig.add_trace(
        go.Scatter(
            x = df["x"],
            y = df["y_upper"],
            mode = "lines",
            line = dict(
                width = 0 if not debug else 5,
            ),
            showlegend = False
        )
    )
    fig.add_trace(
        go.Scatter(
            x = df["x"],
            y = df["y_lower"],
            mode = "lines",
            line = dict(
                width = 0 if not debug else 5,
            ),
            fillcolor=rgb_grayer(hex_to_rgb(orig_color)),
            fill = "tonexty",
            showlegend = False
        )
    )
    return fig
def lines(df, x, y, color_discrete_map=None, with_error=False, frac=0.3, ci=0.95, debug=False, gray_bg=False, **kwargs):
    # fig = px.line(
    #     df,
    #     x = x,
    #     y = y
    # )
    fig = go.Figure()
    if isinstance(y, str):
        y_list = [y]
    elif isinstance(y, list):
        y_list = y
    else:
        raise ValueError("y must be str or list of str")
    for y in y_list:
        color = color_discrete_map[y] if color_discrete_map and y in color_discrete_map else None
        error_df = smooth_with_ci(
            df,
            x = x,
            y = y,
            frac = frac,
            ci = ci
        )
        fig.add_trace(
            go.Scatter(
                x = error_df["x"],
                y = error_df["y_smooth"],
                mode = "lines+markers",
                line = dict(
                    color = color
                ),
                marker = dict(
                    color = color
                ),
                name = y,
            )
        )
        if with_error:
            fig = add_error_band(
                fig,
                error_df,
                orig_color=color if not gray_bg else "#808080",
                debug=debug
            )
    fig.update_layout(
        title = "",
        # xaxis_title = "Time (H)",
        # yaxis_title = "HSV centromere / other odds ratio",
        # yaxis_range = [0,1.2],
        height = 400,
        width = 800,
        margin = dict(l=40, r=20, t=20, b=40),
        **kwargs
    )
    return fig
def line(df, x, y, color=None, color_discrete_map=None, with_error=False, 
                  frac=0.3, ci=0.95, debug=False, gray_bg=False, 
                  hover_data=None, hover_format=None, **kwargs):
    """
    Plot a single trend line with optional error bands and point coloring.
    
    Creates a line plot for a single variable with optional confidence intervals.
    Points can be colored according to a specified column in the dataframe.
    
    Args:
        df (pd.DataFrame): Input data.
        x (str): Column name for x-axis values.
        y (str): Column name for y-axis values.
        color (str, optional): Column name for point coloring. Points will be 
                              colored according to this column's values. 
                              Defaults to None (uniform color).
        color_discrete_map (dict, optional): Color mapping for points when 
                                            `color` is specified. Keys should be 
                                            values from the `color` column.
                                            Defaults to None.
        with_error (bool, optional): Whether to add confidence interval bands.
                                    Defaults to False.
        frac (float, optional): Fraction of data to use for smoothing (LOESS).
                                Defaults to 0.3.
        ci (float, optional): Confidence interval level (0-1). Defaults to 0.95.
        debug (bool, optional): Enable debug mode for error bands. 
                                Defaults to False.
        gray_bg (bool, optional): Use gray color for error bands. Defaults to False.
        **kwargs: Additional keyword arguments passed to fig.update_layout().
    
    Returns:
        plotly.graph_objs.Figure: Plotly figure object.
    
    Examples:
        >>> line(df, x='time', y='value', color='group')
        >>> line(df, x='x', y='y', with_error=True, color_discrete_map={'A': 'red', 'B': 'blue'})
    """
    fig = go.Figure()
    
    # Smooth the data
    error_df = smooth_with_ci(df, x=x, y=y, frac=frac, ci=ci)
    
    # Determine line color
    line_color = None
    if color_discrete_map and y in color_discrete_map:
        line_color = color_discrete_map[y]
    elif color_discrete_map and len(color_discrete_map) == 1:
        line_color = list(color_discrete_map.values())[0]
    
    # Prepare hover template
    def build_hover_template(group_name=None):
        """Build hover template based on available data."""
        template_parts = [f"{x}: %{{x}}", f"{y}: %{{y}}"]
        
        if color and group_name:
            template_parts.append(f"{color}: {group_name}")
        
        if hover_data:
            # Add additional hover data columns
            for col_idx, col in enumerate(hover_data):
                if col in df.columns and col != color:
                    format_spec = hover_format.get(col, '') if hover_format else ''
                    if format_spec:
                        # Fix: escape curly braces properly for Plotly template
                        template_parts.append(f"{col}: %{{customdata[{col_idx}]:{format_spec}}}")
                    else:
                        template_parts.append(f"{col}: %{{customdata[{col_idx}]}}")
        
        template_parts.append("<extra></extra>")
        return "<br>".join(template_parts)
    
    # Add the main line trace
    fig.add_trace(
        go.Scatter(
            x=error_df["x"],
            y=error_df["y_smooth"],
            mode="lines",
            line=dict(color=line_color, width=2),
            name=y,
            showlegend=True,
            hovertemplate=build_hover_template()
        )
    )
    
    # Add scatter points with optional coloring
    if color and color in df.columns:
        color_groups = df.groupby(color)
        
        for group_name, group_data in color_groups:
            point_color = None
            if color_discrete_map and group_name in color_discrete_map:
                point_color = color_discrete_map[group_name]
            elif line_color:
                point_color = line_color
            
            # Prepare custom data for hover
            custom_data = None
            if hover_data:
                hover_cols = [col for col in hover_data if col in group_data.columns and col != color]
                custom_data = group_data[hover_cols].values
            
            fig.add_trace(
                go.Scatter(
                    x=group_data[x],
                    y=group_data[y],
                    mode="markers",
                    marker=dict(
                        color=point_color,
                        size=8,
                        line=dict(width=1, color='white')
                    ),
                    name=f"{y} ({group_name})",
                    showlegend=True,
                    legendgroup=group_name,
                    hovertemplate=build_hover_template(group_name),
                    customdata=custom_data
                )
            )
    else:
        # Add uniform scatter points
        custom_data = None
        if hover_data:
            hover_cols = [col for col in hover_data if col in df.columns]
            customdata = df[hover_cols].values
        
        fig.add_trace(
            go.Scatter(
                x=df[x],
                y=df[y],
                mode="markers",
                marker=dict(
                    color=line_color,
                    size=8,
                    line=dict(width=1, color='white')
                ),
                name=y,
                showlegend=False,
                hovertemplate=build_hover_template(),
                customdata=custom_data
            )
        )
    
    # Add error bands if requested
    if with_error:
        fig = add_error_band(
            fig,
            error_df,
            orig_color=line_color if not gray_bg else "#808080",
            debug=debug
        )
    
    # Update layout
    fig.update_layout(
        title="",
        height=400,
        width=800,
        margin=dict(l=40, r=20, t=20, b=40),
        **kwargs
    )
    
    return fig