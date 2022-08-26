import pandas as pd
import plotly.express as px
def filling_l2r_plotly(rows, cols, features):
    """
    Helper to iterate within row-cols. 
    Plotly_flavor, starts with 1.
    """
    for i in range(rows):
        for j in range(cols):
            k = i * cols + j
            try:
                feature = features[k]
            except IndexError:
                feature = None
            yield i+1, j+1, k, feature
def filling_l2r_mpl(rows, cols, features):
    """
    Helper to iterate within row-cols. 
    Mpl_flavor, starts with 0.
    """
    for i in range(rows):
        for j in range(cols):
            k = i * cols + j
            try:
                feature = features[k]
            except IndexError:
                feature = None
            yield i, j, k, feature
def list2colorlist(celltypes):
    """
    Tranform a data list to a color list, ready for all color argument.
    TODO:
        custom mapper;
        custom color sequence
    """
    # prepare mapper
    color_discrete_sequence = px.colors.qualitative.Plotly
    cat = tuple(set(celltypes))
    mapper = pd.Series(
        index = cat,
        data = [color_discrete_sequence[i % len(color_discrete_sequence)] for i, key in enumerate(cat)],
        name = "color"
    )
    color = [mapper[i] for i in celltypes]
    return color