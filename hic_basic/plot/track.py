
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
                        chorm_label_shift=0.05, y_extra_expansion=0.1, **kwargs):
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
    Output:
        fig: plotly.graph_objs._figure.Figure
    """
    assert chromosomes is not None, "chromosomes must be provided"
    x_extra_expansion = 0

    genome_data = sort_index_chrom(genome_data, chromosomes)

    dat = genome_data.sort_index()
    if value_cols is None:
        value_cols = dat.columns
    dat = dat.assign( # add simple id as x axis
        id = list(range(len(dat))),
        chrom = dat.index.get_level_values(0),
        chromint = dat.index.get_level_values(0).codes
    )
    x_start, x_end = dat["id"].min(), dat["id"].max()
    x_expand = x_end - x_start
    x_start, x_end = x_start - x_expand * x_extra_expansion, x_end + x_expand * x_extra_expansion
    y_start, y_end = dat[value_cols].min().min(), dat[value_cols].max().max()
    y_expand = y_end - y_start
    y_start, y_end = y_start - y_expand * y_extra_expansion, y_end + y_expand * y_extra_expansion
    fig = px.scatter(
        dat,
        y = value_cols,
        x = "id",
        **kwargs
    )
    chrom_left_boundaries = dat[dat["chromint"].diff().fillna(1) > 0]
    mids = pd.Series(chrom_left_boundaries["id"].tolist() + [x_end]).rolling(2).mean().dropna()
    chrom_left_boundaries = chrom_left_boundaries.assign(
        mid = mids.values
    )
    # add chromosome separation lines and names
    for i, row in chrom_left_boundaries.iterrows():
        start, mid = row["id"], row["mid"]
        if row["id"] == x_start:
            # yaxis itself is a line, so no need to add line at x=0
            continue
        fig.add_shape(
            type='line',
            x0=start, y0=y_start,
            x1=start, y1=y_end,
            line=dict(
                color='black',
                width = 1
                )
            )
        # fig.add_annotation(
        #     x=mid,
        #     y=line_y1 + chorm_label_shift - yshift,
        #     text=chroms[i],
        #     showarrow=False,
        #     yshift=10
        #     )
    fig.update_layout(
        height = 500,
        width = 2000,
        plot_bgcolor = "white",
        xaxis = dict(
            title = "Chromosome",
            showline=True,
            linecolor="black",
            mirror=True,
            range = [x_start, x_end],
            ticks = "",
            tickvals = chrom_left_boundaries["mid"].tolist(),
            ticktext = chrom_left_boundaries["chrom"].tolist()
        ),
        yaxis = dict(
            title = value_cols[0],
            showline=True,
            linecolor="black",
            mirror=True,
            range = [y_start, y_end]
        )
    )
    #fig.show(renderer="png")
    return fig