def plot_gene_trend(adata, gene, order_col="velocity_pseudotime", additional=None):
    """
    Ploting expression trends along time.
    Must have pSmooth in uns; normalized and log1p X counts
    """
    fig = go.Figure()
    try:
        smooth_counts = np.log(adata.uns["pSmooth"].loc[gene, adata.uns["pSmooth"].columns[::-1]] + 1)
    except KeyError:
        print("No Smoothing Model Found.")
        smooth_counts = []
    try:
        gene_count = adata[adata.obs.sort_values(order_col).index,:].obs_vector(gene)
    except KeyError:
        if isinstance(additional, ad.AnnData):
            print("Using additional counts.")
            if gene not in additional.var_names:
                gene_count = []
            else:
                gene_count = additional[additional.obs.sort_values(order_col).index,:].obs_vector(gene)
        else:
            gene_count = []
    x_num = max(len(smooth_counts),len(gene_count))
    fig.add_trace(
        go.Scatter(
            y = smooth_counts,
            x = np.linspace(0, x_num, len(smooth_counts)),
            mode = "lines",
            name = gene + " GAM fitting"
        )
    )
    """fig.add_trace(
        go.Scatter(
            y = adata.uns["pCells"].loc[gene, adata.obs.sort_values("velocity_pseudotime").index],
            mode = "markers"
        )
    )"""

    fig.add_trace(
        go.Scatter(
            y = gene_count,
            x = np.linspace(0, x_num, len(gene_count)),
            mode = "markers",
            name = gene + " normalized counts"
        )
    )
    fig.update_layout(
        height = 500,
        width = 700,
        title = gene
    )
    return fig
def plot_gene_trends(data, additional, genes, order_col, ncols):
    genes = [i for i in genes if i in data.var_names]
    nrows = ceil(len(genes)/ncols)
    fig = make_subplots(rows = nrows, cols = ncols, shared_xaxes=True, shared_yaxes=True)
    for i, j, k, gene in filling_l2r_plotly(nrows, ncols, genes):
        if gene is not None:
            subfig = plot_gene_trend(data, gene, order_col, additional)
            for trace in subfig.data:
                fig.add_trace(
                    trace,
                    row = i,
                    col = j
                )
    return fig