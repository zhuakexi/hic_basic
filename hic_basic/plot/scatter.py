import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.stats import norm
from .utils import hex_to_rgb

### --- line plot with error band --- ###

def smooth_with_ci(df, x, y, frac=0.3, ci=0.95, skip_smooth=False, error_type="ci", sem_col=None):
    """
    Apply LOESS smoothing to one-row-per-x curve data and estimate ribbons.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame containing one y value for each x value. Repeated
        observations must be aggregated upstream before calling this function.
    x : str
        Column name for independent variable (e.g., time points).
    y : str
        Column name for dependent variable (e.g., aggregated mean values).
    frac : float, optional
        Fraction of data used for each local regression (0 < frac <= 1).
        Default is 0.3.
    ci : float, optional
        Confidence interval level (0 < ci < 1, e.g., 0.95 for 95% CI).
        Default is 0.95.
    skip_smooth : bool, optional
        If True, bypass LOESS smoothing and use raw y values directly.
        Default is False.
    error_type : str, optional
        Type of error estimation: "ci" for pooled residual CI, "sem" for
        per-point SEM. Default is "ci". When "sem", requires sem_col to be
        specified.
    sem_col : str, optional
        Column name containing pre-computed SEM values for each x value.
        Required when error_type="sem". SEM is not computed internally.

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns:
        - 'x': Original x values (sorted)
        - 'y': Original y values
        - 'y_smooth': Smoothed y values (or raw y if skip_smooth=True)
        - 'y_upper': Upper bound of error band
        - 'y_lower': Lower bound of error band
        - 'y_sem': SEM values (if error_type="sem")

    Notes
    -----
    This function assumes df already describes a single curve. It does not
    group, bin, or aggregate repeated observations by x.
    """
    assert x in df.columns and y in df.columns, \
        f"Columns '{x}' and '{y}' must exist in the DataFrame."
    # Validate CI
    if ci > 1:
        ci = ci / 100.0  # e.g., 95 -> 0.95
    if not (0 < ci < 1):
        raise ValueError("ci must be between 0 and 1.")
    
    # Validate error_type
    if error_type not in ("ci", "sem"):
        raise ValueError("error_type must be 'ci' or 'sem'.")
    
    # Validate sem_col when error_type="sem"
    if error_type == "sem":
        if sem_col is None:
            raise ValueError("sem_col must be specified when error_type='sem'.")
        if not isinstance(sem_col, str):
            raise ValueError(
                "smooth_with_ci() accepts a single y column, so sem_col must be a single column name string."
            )
        if sem_col not in df.columns:
            raise ValueError(f"Column '{sem_col}' must exist in the DataFrame when error_type='sem'.")

    df = df.sort_values(by=x)
    x_arr = df[x].values
    y_arr = df[y].values
    x_arr = np.array(x_arr, dtype=float)
    y_arr = np.array(y_arr, dtype=float)

    # Apply LOESS smoothing or use raw values
    if skip_smooth:
        smoothed = y_arr  # no smoothing, use raw y
    else:
        smoothed = lowess(y_arr, x_arr, frac=frac, return_sorted=False)

    # z-score for CI/SEM multiplier
    z = norm.isf((1 - ci) / 2)

    # Compute error bounds based on error_type
    if error_type == "ci":
        # Pooled residual CI (constant width)
        residuals = y_arr - smoothed
        sigma = np.nanstd(residuals, ddof=1)  # avoid division bias
        upper = smoothed + z * sigma
        lower = smoothed - z * sigma
        sem_values = None
    else:  # error_type == "sem"
        # Per-bin SEM (varying width from pre-computed column)
        sem_values = df[sem_col].values
        upper = smoothed + z * sem_values
        lower = smoothed - z * sem_values

    result = pd.DataFrame({
        'x': x_arr,
        'y': y_arr,
        'y_smooth': smoothed,
        'y_upper': upper,
        'y_lower': lower
    })
    
    if sem_values is not None:
        result['y_sem'] = sem_values
    
    return result
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
def add_error_band(fig, df, orig_color=px.colors.qualitative.Plotly[0], debug=False, error_type="ci"):
    """
    Add error visualization to a figure.
    
    Parameters
    ----------
    fig : plotly.graph_objects.Figure
        Figure to add error to.
    df : pandas.DataFrame
        DataFrame with x, y_upper, y_lower columns.
    orig_color : str
        Color of original line, will use a grayer version for the ribbon.
    debug : bool
        If True, show ribbon outline.
    error_type : str
        "ci" for pooled CI ribbon, "sem" for varying-width ribbon from
        pre-computed SEM values.
    
    Returns
    -------
    plotly.graph_objects.Figure
        Figure with error added.
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
            showlegend = False,
            name = f"{error_type.upper()} upper"
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
            showlegend = False,
            name = f"{error_type.upper()} ribbon"
        )
    )
    return fig
def _normalize_lines_sem_cols(y_list, sem_col, y_is_list, error_type):
    """
    Normalize sem_col for lines().

    Returns one sem column name per y series when error_type="sem".
    For other error types, returns a placeholder list without validation.
    """
    if error_type != "sem":
        return [sem_col] * len(y_list)

    if sem_col is None:
        raise ValueError("sem_col must be specified when error_type='sem'.")

    if isinstance(sem_col, str):
        return [sem_col] * len(y_list)

    if isinstance(sem_col, (list, tuple)):
        if not y_is_list:
            raise ValueError(
                "When y is a single column name, sem_col must be a single column name string."
            )
        if len(sem_col) != len(y_list):
            raise ValueError("sem_col must have the same length as y when both are lists.")
        if not all(isinstance(col, str) for col in sem_col):
            raise ValueError("Each sem_col entry must be a column name string.")
        return list(sem_col)

    raise ValueError("sem_col must be a column name string or a list of column name strings.")
def lines(df, x, y, color_discrete_map=None, with_error=False, frac=0.3, ci=0.95, 
          debug=False, gray_bg=False, skip_smooth=False, error_type="ci", sem_col=None, **kwargs):
    """
    Plot multiple line trends with optional error visualization.

    Args:
        df (pd.DataFrame): Input curve data with one row per x value.
        x (str): Column name for x-axis values.
        y (str or list): Column name(s) for y-axis values. Each y series must
            already be aligned one-to-one with x.
        color_discrete_map (dict, optional): Color mapping for each y column.
        with_error (bool): Whether to add error visualization.
        frac (float): Fraction of data for LOESS smoothing.
        ci (float): Confidence interval level (0-1).
        debug (bool): Enable debug mode for error bands.
        gray_bg (bool): Use gray color for error bands.
        skip_smooth (bool): If True, bypass LOESS and plot raw y values.
        error_type (str): "ci" for pooled CI ribbon, "sem" for ribbon from
            pre-computed SEM values at each x.
        sem_col (str or list[str]): Column name containing pre-computed SEM
            values for each x value, or a list aligned one-to-one with y when
            y is a list. Required when error_type="sem".
        **kwargs: Additional arguments passed to fig.update_layout().

    Returns:
        plotly.graph_objects.Figure

    Notes:
        lines() does not aggregate repeated observations. If your raw data has
        multiple rows per x value, compute summary statistics upstream and pass
        the aggregated curve plus an optional sem_col.
        When y is a list and error_type="sem", sem_col may be a single shared
        SEM column name or a list aligned one-to-one with y.
    """
    fig = go.Figure()
    if isinstance(y, str):
        y_list = [y]
        y_is_list = False
    elif isinstance(y, list):
        y_list = y
        y_is_list = True
    else:
        raise ValueError("y must be str or list of str")
    sem_cols = _normalize_lines_sem_cols(y_list, sem_col, y_is_list, error_type)

    for y, y_sem_col in zip(y_list, sem_cols):
        color = color_discrete_map[y] if color_discrete_map and y in color_discrete_map else None
        # Use default Plotly color if no color specified
        error_color = color if color else px.colors.qualitative.Plotly[0]
        
        error_df = smooth_with_ci(
            df,
            x = x,
            y = y,
            frac = frac,
            ci = ci,
            skip_smooth = skip_smooth,
            error_type = error_type,
            sem_col = y_sem_col
        )
        
        trace_kwargs = {
            "x": error_df["x"],
            "y": error_df["y_smooth"],
            "mode": "lines+markers",
            "line": dict(color=color),
            "marker": dict(color=color),
            "name": y,
        }

        fig.add_trace(go.Scatter(**trace_kwargs))

        if with_error:
            fig = add_error_band(
                fig,
                error_df,
                orig_color=error_color if not gray_bg else "#808080",
                debug=debug,
                error_type=error_type
            )
    
    # Handle title separately to avoid conflict with kwargs
    layout_args = {
        "height": 400,
        "width": 800,
        "margin": dict(l=40, r=20, t=20, b=40),
    }
    if "title" not in kwargs:
        layout_args["title"] = ""
    layout_args.update(kwargs)
    fig.update_layout(**layout_args)
    return fig
def line(df, x, y, color=None, color_discrete_map=None, with_error=False,
                  frac=0.3, ci=0.95, debug=False, gray_bg=False,
                  hover_data=None, hover_format=None,
                  skip_smooth=False, error_type="ci", sem_col=None, **kwargs):
    """
    Plot a single trend line with optional error bands and point coloring.

    Creates a line plot for a single variable with optional confidence intervals.
    Points can be colored according to a specified column in the dataframe.

    Args:
        df (pd.DataFrame): Input curve data with one row per x value.
        x (str): Column name for x-axis values.
        y (str): Column name for y-axis values. Must already be one-to-one with x.
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
        hover_data (list, optional): Additional columns to show in hover tooltip.
        hover_format (dict, optional): Format specifiers for hover columns.
        skip_smooth (bool, optional): If True, bypass LOESS and plot raw y values.
                                      Defaults to False.
        error_type (str, optional): "ci" for pooled CI ribbon, "sem" for
                                    ribbon from pre-computed SEM values at
                                    each x. Defaults to "ci".
        sem_col (str, optional): Column name with pre-computed SEM values for
                                 each x value. Required when error_type="sem".
        **kwargs: Additional keyword arguments passed to fig.update_layout().

    Returns:
        plotly.graph_objs.Figure: Plotly figure object.

    Notes:
        line() does not aggregate repeated observations. If your raw data has
        multiple rows per x value, aggregate it upstream so that each x maps to
        exactly one y, then pass SEM explicitly through sem_col if needed.

    Examples:
        >>> line(curve_df, x='time', y='value', color='group')
        >>> line(curve_df, x='x', y='y', with_error=True, color_discrete_map={'A': 'red', 'B': 'blue'})
        >>> line(agg_df, x='bin', y='mean', with_error=True, error_type="sem", sem_col="sem")
    """
    fig = go.Figure()

    if error_type == "sem" and isinstance(sem_col, (list, tuple)):
        raise ValueError(
            "line() accepts a single y column, so sem_col must be a single column name string."
        )

    # Smooth the data
    error_df = smooth_with_ci(df, x=x, y=y, frac=frac, ci=ci, 
                               skip_smooth=skip_smooth, error_type=error_type, sem_col=sem_col)
    
    # Determine line color (use default if None)
    line_color = None
    if color_discrete_map and y in color_discrete_map:
        line_color = color_discrete_map[y]
    elif color_discrete_map and len(color_discrete_map) == 1:
        line_color = list(color_discrete_map.values())[0]
    # Use default Plotly color if no color specified
    error_color = line_color if line_color else px.colors.qualitative.Plotly[0]
    
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
    
    main_trace_kwargs = {
        "x": error_df["x"],
        "y": error_df["y_smooth"],
        "mode": "lines",
        "line": dict(color=line_color, width=2),
        "name": y,
        "showlegend": True,
        "hovertemplate": build_hover_template()
    }

    fig.add_trace(go.Scatter(**main_trace_kwargs))
    
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
    
    if with_error:
        fig = add_error_band(
            fig,
            error_df,
            orig_color=error_color if not gray_bg else "#808080",
            debug=debug,
            error_type=error_type
        )
    
    # Update layout (handle title separately to avoid conflict with kwargs)
    layout_args = {
        "height": 400,
        "width": 800,
        "margin": dict(l=40, r=20, t=20, b=40),
    }
    if "title" not in kwargs:
        layout_args["title"] = ""
    layout_args.update(kwargs)
    fig.update_layout(**layout_args)
    
    return fig
