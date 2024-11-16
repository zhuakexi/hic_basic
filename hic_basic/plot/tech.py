from math import ceil
from plotly.subplots import make_subplots
import plotly.graph_objects as go

from .general import ColorDiscreteMap, color_discrete_sequence
from .utils  import filling_l2r_plotly
from ..wet.afbb import formal_feature_names

def tech_box(annote, grouping="condition", features=["defaults"], ncols=3, size=200,
    color_discrete_sequence=color_discrete_sequence, use_label=True, labels:dict=None, **kwargs)->go.Figure:
    """
    DNA and RNA tech metrics.
    Input:
        annote: sample * feature dataframe
        grouping: col name to group by
            use all if grouping is None
        features: features to be plot. ["defaults", ...]
        ncols: number of cols of subplots
        size: subplot size
        color_discrete_sequence: cycling pallete for each group
        use_label: use a more formal name for each feature
            formal names are defined in wet.afbb.formal_feature_names
        labels: customized formal feature names
    TODO:
        set group order
    """
    if grouping is not None:
        cdmap = ColorDiscreteMap(color_discrete_sequence)
        cdmap.extend(annote[grouping].unique())
        
    real_features = []
    for i in features:
        if i == "defaults":
            real_features.extend(["contacts","umis","rna_ratio","con_per_reads","umis_per_reads"])
        else:
            real_features.append(i)
    
    nrows = ceil(len(real_features)/ncols)
    if use_label:
        # translate subplot titles
        labels = formal_feature_names if labels is None else {**formal_feature_names, **labels}
        subplot_titles = [
            labels[feature] if feature in labels else feature
            for feature in real_features
            ]
    else:
        subplot_titles = real_features
    fig = make_subplots(rows=nrows, cols=ncols, subplot_titles=subplot_titles)
    for row, col, _, feature in filling_l2r_plotly(nrows, ncols, real_features):
        if feature is not None:# else skip
            for group, chunk in annote.groupby(grouping) if grouping is not None else [(None, annote)]:
                if group is None:
                    # use single trace style
                    per_trace_kwargs = {"name":feature,"marker_color":"black","fillcolor":"gray"}
                else:
                    # use separate fillcolor for each group
                    per_trace_kwargs = {"name":group,"marker_color":"black","fillcolor":cdmap[group]}
                fig.add_trace(
                    go.Box(
                        y = chunk[feature],
                        showlegend = False,
                        **per_trace_kwargs
                    ),
                    row = row,
                    col = col
                )
            if grouping is not None:
                fig.update_xaxes(
                    title_text = "Group",
                    row = row,
                    col = col
                )
            else: # if value of x axis is useless, hide it
                fig.update_xaxes(
                    ticks="",
                    row = row,
                    col = col
                )
                
    fig.update_traces(
        showlegend = True,row=1,col=1
    )
    fig.update_layout(
        height = nrows*size,
        width = ncols*size,
    )
    return fig