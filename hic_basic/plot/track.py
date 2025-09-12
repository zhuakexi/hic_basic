
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

### --- helpers --- ###

def sort_index_chrom(df, chromosomes):
    if isinstance(chromosomes, pd.DataFrame):
        pass
    else:
        try:
            chromosomes = chromosomes.data
        except AttributeError:
            raise ValueError("chromosomes must be a dataframe or Feature obj")
    index_frame = df.index.to_frame()
    index_frame["chrom"] = index_frame["chrom"].astype(
        chromosomes.index.dtype
    )
    df.index = pd.MultiIndex.from_frame(index_frame)
    df = df.sort_index()
    return df

### --- genome scatter --- ###

def plot_genome_scatter(genome_data, chromosomes=None, value_cols=None, shift_unit=0.05, shift_threshold=100,
                        chorm_label_shift=0.05, y_extra_expansion=0.1, range_y=None, **kwargs):
    """
    Scatter genome wide data. Sepaarate chromosomes by vertical lines and add chromosome names.
    TODO: add real positions, now only add 1,2,3,4... as x axis.
    Input:
        genome_data: pd.DataFrame, index is multiindex with level 0 as chromosome names 
            and level 1 as start position. First column is the data to plot.
            Will sort the index first. To plot chromosome
            according to your desired order, set the index level 0 as ordered category.
        shift_unit: float, default 0.05, the unit of shift when two chromosomes are too close.
        shift_threshold: int, default 100, the threshold to determine if two chromosomes are too close.
        chorm_label_shift: float, default 0.05, the shift of chromosome name from the top of the plot.
        y_extra_expansion: float, default 0.1, the extra expansion (relative to min, max) of y axis.
        range_y: tuple or list, default None, the y axis range.
    Output:
        fig: plotly.graph_objs._figure.Figure
    """
    assert chromosomes is not None, "chromosomes must be provided"
    x_extra_expansion = 0

    genome_data = sort_index_chrom(genome_data, chromosomes)

    if value_cols is None:
        value_cols = genome_data.columns
    dat = genome_data.assign( # add simple id as x axis
        id = list(range(len(genome_data))),
    )
    fig = px.scatter(
        dat,
        y = value_cols,
        x = "id",
        **kwargs
    )
    fig.update_layout(
        height = 500,
        width = 2000,
        plot_bgcolor = "white"
    )
    fig = add_chrom_grid(
        genome_data,
        fig,
        fig_type="scatter",
        x_extra_expansion=x_extra_expansion,
        y_extra_expansion=y_extra_expansion,
        value_cols=value_cols,
        range_y=range_y
    )
    return fig

### --- genome heatmap --- ###
def plot_genome_heatmap(genome_data, chromosomes=None, **kwargs):
    """
    Heatmap genome wide data. Sepaarate chromosomes by vertical lines and add chromosome names.
    Input:
        genome_data: pd.DataFrame, index is multiindex with level 0 as chromosome names
            columns are samples. Will sort the index first.
        chromosomes: pd.DataFrame or Feature, the chromosome information used to sort the genome_data.
    Output:
        fig: plotly.graph_objs._figure.Figure
    """
    assert chromosomes is not None, "chromosomes must be provided"
    df = genome_data.copy()
    df = sort_index_chrom(df, chromosomes)

    fig = go.Figure()
    fig.add_trace(
        go.Heatmap(
            z=df.T.values,
            **kwargs
        )
    )
    fig.update_layout(
        height = 500,
        width = 600
    )
    fig = add_chrom_grid(
        df,
        fig,
        fig_type="heatmap",
        x_extra_expansion=0,
        y_extra_expansion=0, # for heatmap, y axis is sample index, so the best is y_extra_expansion=0
    )
    return fig

def add_chrom_grid(df, fig, fig_type="scatter", x_extra_expansion=0, y_extra_expansion=0.1, value_cols=None, range_y=None):
    """
    Add chromosome grid lines to an existing genome track figure.
    Input:
        df: pd.DataFrame, the genome data used to plot the figure. Index is multiindex with level 0 as chromosome names 
            and level 1 as start position. First column is the data to plot.
            1. For scatter plot, use first value column to determine y axis title.
            2. For heatmap, df is `pos * samples`, use "Sample" as y axis title, each column is a sample.
        fig: plotly.graph_objs._figure.Figure, the figure to add grid lines to.
        fig_type: str, default "scatter", type of the figure, currently only support "scatter", "heatmap".
        x_extra_expansion: float, default 0, the extra expansion (relative to min, max) of x axis.
        y_extra_expansion: float, default 0.1, the extra expansion (relative to min, max) of y axis.
            For heatmap, y axis is sample index, so the best is y_extra_expansion=0.
        range_y: tuple or list, default None, the y axis range.
    Output:
        fig: plotly.graph_objs._figure.Figure, the figure with added grid lines.
    """
    assert fig_type in ["scatter", "heatmap"], "fig_type must be 'scatter' or 'heatmap'"
    assert isinstance(
        df.index.get_level_values("chrom"),
        pd.CategoricalIndex), "df index level 0 must be categorical dtype"

    dat = pd.DataFrame(
        index = df.index
    )
    dat = dat.assign( # add simple id as x axis
        id = list(range(len(dat))),
        chrom = dat.index.get_level_values(0),
        chromint = dat.index.get_level_values(0).codes
    )
    x_start, x_end = dat["id"].min(), dat["id"].max()
    x_expand = x_end - x_start
    x_start, x_end = x_start - x_expand * x_extra_expansion, x_end + x_expand * x_extra_expansion
    if range_y is None:
        if fig_type == "scatter":
            if value_cols is None:
                value_cols = df.columns
            y_start, y_end = df[value_cols].min().min(), df[value_cols].max().max()
        elif fig_type == "heatmap":
            y_start, y_end = 0 - 0.5, len(df.columns) - 0.5
    else:
        y_start, y_end = range_y
    y_expand = y_end - y_start
    y_start, y_end = y_start - y_expand * y_extra_expansion, y_end + y_expand * y_extra_expansion

    # Detect chromosome boundaries by checking changes in 'chromint' (chromosome integer codes)
    # A new chromosome starts when the difference between consecutive 'chromint' values is greater than 0
    chrom_left_boundaries = dat[dat["chromint"].diff().fillna(1) > 0]

    # Calculate midpoints between the start of each chromosome and the start of the next one
    # This helps position vertical lines between chromosomes
    mids = pd.Series(chrom_left_boundaries["id"].tolist() + [x_end]).rolling(2).mean().dropna()

    # Assign midpoints to the chromosome boundary DataFrame for reference
    chrom_left_boundaries = chrom_left_boundaries.assign(
        mid = mids.values  # Store midpoints for annotation or line positioning
    )

    # Add vertical separation lines between chromosomes
    # print(y_start, y_end)
    for i, row in chrom_left_boundaries.iterrows():
        start, mid = row["id"], row["mid"]
        if row["id"] == x_start:
            # Skip adding a line at the very first chromosome's start (already marked by the y-axis)
            continue
        fig.add_shape(
            type='line',
            x0=start, y0=y_start,  # Start at the bottom of the plot
            x1=start, y1=y_end,    # End at the top of the plot
            line=dict(
                color='black',     # Line color
                width=1            # Line thickness
            )
        )
    
    # Update axis ranges to accommodate the added grid lines
    xaxis = dict(
        title = "Chromosome",
        showline=True,
        linecolor="black",
        mirror=True,
        range = [x_start, x_end],
        ticks = "",
        tickvals = chrom_left_boundaries["mid"].tolist(),
        ticktext = chrom_left_boundaries["chrom"].tolist()
    )
    yaxis = dict(
        title = value_cols[0] if fig_type == "scatter" else "Sample",
        showline=True,
        linecolor="black",
        mirror=True,
        range = [y_start, y_end]
    )
    fig.update_layout(
        xaxis = xaxis,
        yaxis = yaxis
    )
    return fig