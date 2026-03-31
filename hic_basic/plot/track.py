
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from typing import Dict, List, Optional, Union, Tuple
from pathlib import Path
from .gene_track_parser import (
    parse_annotation_file,
    filter_genes_in_region,
    Gene,
    Transcript,
    Interval
)

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


### --- Gene track --- ###


def _split_exon_by_cds_utr(
    exon: Interval,
    cds_intervals: List[Interval],
    utr5_intervals: List[Interval],
    utr3_intervals: List[Interval]
) -> List[Tuple[str, Interval]]:
    """
    Split an exon into CDS/UTR/exon parts based on overlapping features.

    Returns list of tuples: (feature_type, interval)
    """
    result = []
    exon_start, exon_end = exon.start, exon.end

    # Collect all sub-intervals with their types
    intervals_with_type = []

    for cds in cds_intervals:
        overlap_start = max(exon_start, cds.start)
        overlap_end = min(exon_end, cds.end)
        if overlap_start <= overlap_end:
            intervals_with_type.append(('cds', overlap_start, overlap_end))

    for utr in utr5_intervals + utr3_intervals:
        overlap_start = max(exon_start, utr.start)
        overlap_end = min(exon_end, utr.end)
        if overlap_start <= overlap_end:
            intervals_with_type.append(('utr', overlap_start, overlap_end))

    # If no CDS/UTR overlaps, it's just an exon
    if not intervals_with_type:
        result.append(('exon', exon))
        return result

    # Sort by start position
    intervals_with_type.sort(key=lambda x: x[1])

    # Merge overlapping intervals and track types
    merged = []
    for feat_type, start, end in intervals_with_type:
        if merged and merged[-1][0] == feat_type and merged[-1][2] >= start - 1:
            # Merge with previous
            merged[-1] = (feat_type, merged[-1][1], max(merged[-1][2], end))
        else:
            merged.append((feat_type, start, end))

    # Create intervals
    for feat_type, start, end in merged:
        interval = Interval(start=start, end=end, attrs=exon.attrs)
        result.append((feat_type, interval))

    return result


def _get_transcript_row(
    transcript: Transcript,
    y_position: float,
    utr_height: float = 0.30,
    exon_height: float = 0.40,
    cds_height: float = 0.60,
    backbone_y: float = None
) -> List[dict]:
    """
    Create shapes for a single transcript row.

    Returns list of plotly shapes (rectangles for exons/CDS/UTR, lines for introns).
    """
    shapes = []

    # Sort exons by start position
    sorted_exons = sorted(transcript.exons, key=lambda e: e.start)

    if not sorted_exons:
        return shapes

    # Determine backbone y position (middle of track)
    if backbone_y is None:
        backbone_y = y_position + 0.5

    # Draw intron backbone (thin line from transcript start to end)
    if len(sorted_exons) > 1:
        # Draw backbone segments between exons
        for i in range(len(sorted_exons) - 1):
            intron_start = sorted_exons[i].end
            intron_end = sorted_exons[i + 1].start
            shapes.append(dict(
                type='line',
                x0=intron_start,
                x1=intron_end,
                y0=backbone_y,
                y1=backbone_y,
                line=dict(color='black', width=0.5),
                xref='x',
                yref='y'
            ))
            # Add strand arrow on intron
            _add_strand_arrow(
                shapes, intron_start, intron_end, backbone_y,
                transcript.strand, arrow_type='intron'
            )

    # Process each exon
    for exon in sorted_exons:
        # Split exon into CDS/UTR parts
        parts = _split_exon_by_cds_utr(
            exon, transcript.cds, transcript.utr5, transcript.utr3
        )

        for feat_type, interval in parts:
            # Determine height and color based on feature type
            if feat_type == 'cds':
                height = cds_height
                color = 'darkblue'
                y0 = y_position + (1 - height) / 2
            elif feat_type == 'utr':
                height = utr_height
                color = 'lightblue'
                y0 = y_position + (1 - height) / 2
            else:  # exon
                height = exon_height
                color = 'skyblue'
                y0 = y_position + (1 - height) / 2

            # Draw exon rectangle
            shapes.append(dict(
                type='rect',
                x0=interval.start,
                x1=interval.end,
                y0=y0,
                y1=y0 + height,
                fillcolor=color,
                line=dict(color='black', width=0.5),
                xref='x',
                yref='y'
            ))

            # Add strand arrow on exon
            if feat_type == 'exon':
                _add_strand_arrow(
                    shapes, interval.start, interval.end, y0 + height / 2,
                    transcript.strand, arrow_type='exon'
                )

    return shapes


def _add_strand_arrow(
    shapes: List[dict],
    start: int,
    end: int,
    y: float,
    strand: str,
    arrow_type: str = 'exon',
    arrow_y_offset: float = 0.15  # Y offset in track coordinates (not genomic)
):
    """Add strand direction arrow to shapes list.
    
    TODO: Implement proper strand direction visualization.
    Current issues:
    - Arrow direction is not correctly rendered
    - Solid triangle arrows are not visually appealing
    - Consider implementing fishbone-style arrows instead
    
    Parameters
    ----------
    shapes : List[dict]
        List to append shapes to.
    start : int
        Start position (genomic coordinate).
    end : int
        End position (genomic coordinate).
    y : float
        Y position in track coordinates (typically 0-1 per track).
    strand : str
        '+' or '-' for strand direction.
    arrow_type : str
        'exon' or 'intron' - determines arrow style.
    arrow_y_offset : float
        Vertical offset for arrow head in track coordinates. Default 0.15.
    """
    # TODO: Temporarily disabled arrow rendering until proper implementation
    # See issue: arrow direction and styling needs improvement
    return
    
    # --- Original arrow code (commented out) ---
    # if strand not in ('+', '-'):
    #     return
    #
    # # Calculate arrow position (middle of segment)
    # mid = (start + end) / 2
    # segment_length = end - start
    #
    # # Skip if segment is too short
    # if segment_length < 100:
    #     return
    #
    # # Arrow x-length in genomic coordinates
    # arrow_x_length = min(0.15 * segment_length, 500)
    #
    # if strand == '+':
    #     # Arrow pointing right
    #     if arrow_type == 'intron':
    #         # Triangle arrowhead using path
    #         shapes.append(dict(
    #             type='path',
    #             path=f'M {mid - arrow_x_length/2} {y} L {mid + arrow_x_length/2} {y} L {mid} {y - arrow_y_offset} Z',
    #             fillcolor='black',
    #             line=dict(width=0),
    #             xref='x',
    #             yref='y'
    #         ))
    #     else:
    #         # Simple arrow on exon - horizontal line
    #         shapes.append(dict(
    #             type='line',
    #             x0=mid - arrow_x_length/2,
    #             x1=mid + arrow_x_length/2,
    #             y0=y,
    #             y1=y,
    #             line=dict(color='black', width=2),
    #             xref='x',
    #             yref='y'
    #         ))
    #         # Arrow head - diagonal lines with limited y offset
    #         shapes.append(dict(
    #             type='line',
    #             x0=mid + arrow_x_length/2,
    #             x1=mid + arrow_x_length/4,
    #             y0=y,
    #             y1=y - arrow_y_offset,
    #             line=dict(color='black', width=2),
    #             xref='x',
    #             yref='y'
    #         ))
    #         shapes.append(dict(
    #             type='line',
    #             x0=mid + arrow_x_length/2,
    #             x1=mid + arrow_x_length/4,
    #             y0=y,
    #             y1=y + arrow_y_offset,
    #             line=dict(color='black', width=2),
    #             xref='x',
    #             yref='y'
    #         ))
    # else:  # strand == '-'
    #     # Arrow pointing left
    #     if arrow_type == 'intron':
    #         shapes.append(dict(
    #             type='path',
    #             path=f'M {mid + arrow_x_length/2} {y} L {mid - arrow_x_length/2} {y} L {mid} {y - arrow_y_offset} Z',
    #             fillcolor='black',
    #             line=dict(width=0),
    #             xref='x',
    #             yref='y'
    #         ))
    #     else:
    #         shapes.append(dict(
    #             type='line',
    #             x0=mid + arrow_x_length/2,
    #             x1=mid - arrow_x_length/2,
    #             y0=y,
    #             y1=y,
    #             line=dict(color='black', width=2),
    #             xref='x',
    #             yref='y'
    #         ))
    #         shapes.append(dict(
    #             type='line',
    #             x0=mid - arrow_x_length/2,
    #             x1=mid - arrow_x_length/4,
    #             y0=y,
    #             y1=y - arrow_y_offset,
    #             line=dict(color='black', width=2),
    #             xref='x',
    #             yref='y'
    #         ))
    #         shapes.append(dict(
    #             type='line',
    #             x0=mid - arrow_x_length/2,
    #             x1=mid - arrow_x_length/4,
    #             y0=y,
    #             y1=y + arrow_y_offset,
    #             line=dict(color='black', width=2),
    #             xref='x',
    #             yref='y'
    #         ))


def _get_gene_label(gene: Gene) -> str:
    """Get gene label according to priority rules."""
    if gene.gene_name:
        return gene.gene_name
    # Try Name attribute
    if gene.attrs.get('Name'):
        return gene.attrs['Name']
    # Fall back to gene_id
    return gene.gene_id


def _get_transcript_label(transcript: Transcript) -> str:
    """Get transcript label according to priority rules."""
    if transcript.transcript_name:
        return transcript.transcript_name
    # Try Name attribute
    if transcript.attrs.get('Name'):
        return transcript.attrs['Name']
    # Fall back to transcript_id
    return transcript.transcript_id


def plot_gene_track(
    annotation_file: Union[str, Path],
    region: str,
    mode: str = "all_transcripts",
    title: Optional[str] = None,
    height: int = 400,
    width: int = 1000,
    utr_height: float = 0.30,
    exon_height: float = 0.40,
    cds_height: float = 0.60,
    row_spacing: float = 0.2,
    show_labels: bool = True,
    label_font_size: int = 10,
) -> go.Figure:
    """
    Plot gene annotation track in genome browser style.

    This function creates a Plotly figure showing gene structures including
    exons, introns, CDS, and UTRs with strand direction indicators.

    Parameters
    ----------
    annotation_file : str or Path
        Path to GTF or GFF3 annotation file.
    region : str
        Genomic region to display, e.g., "chr1:1000000-2000000".
    mode : str, optional
        Display mode. Options:
        - "all_transcripts": Show all transcripts for each gene
        - "canonical": Show only canonical/longest transcript per gene
        - "gene_only": Show collapsed gene models (one row per gene)
        Default is "all_transcripts".
    title : str, optional
        Plot title. If None, uses region as title.
    height : int, optional
        Figure height in pixels. Default is 400.
    width : int, optional
        Figure width in pixels. Default is 1000.
    utr_height : float, optional
        Relative height of UTR boxes (0-1). Default is 0.30.
    exon_height : float, optional
        Relative height of exon boxes (0-1). Default is 0.40.
    cds_height : float, optional
        Relative height of CDS boxes (0-1). Default is 0.60.
    row_spacing : float, optional
        Spacing between transcript rows. Default is 0.2.
    show_labels : bool, optional
        Whether to show gene/transcript labels. Default is True.
    label_font_size : int, optional
        Font size for labels. Default is 10.

    Returns
    -------
    plotly.graph_objects.Figure
        Plotly figure with gene track visualization.

    Examples
    --------
    >>> fig = plot_gene_track(
    ...     annotation_file="annotations.gtf",
    ...     region="chr1:1000000-2000000",
    ...     mode="all_transcripts",
    ...     title="Gene Track Example"
    ... )
    >>> fig.show()
    
    Notes
    -----
    TODO: Strand direction arrows are currently disabled pending improved
    implementation. Planned improvements include:
    - Fishbone-style arrows instead of solid triangles
    - Proper direction rendering for both + and - strands
    - Better visual integration with exon/intron structures
    """
    # Parse region
    if ':' in region:
        seqid, coords = region.split(':', 1)
        if '-' in coords:
            start, end = map(int, coords.replace(',', '').split('-'))
        else:
            start, end = 0, int(coords.replace(',', ''))
    else:
        seqid = region
        start, end = 0, 1000000  # Default range

    # Parse annotation file
    genes = parse_annotation_file(annotation_file)

    # Filter genes in region
    genes_in_region = filter_genes_in_region(genes, seqid, start, end)

    if not genes_in_region:
        # Create empty figure with message
        fig = go.Figure()
        fig.add_annotation(
            text=f"No genes found in {region}",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig

    # Collect transcripts to display
    transcripts_to_show = []
    genes_to_label = []

    for gene_id, gene in sorted(genes_in_region.items(), key=lambda x: x[1].start):
        if mode == "gene_only":
            # Select representative transcript
            if gene.transcripts:
                # Prefer longest CDS transcript
                cds_transcripts = [t for t in gene.transcripts if t.cds]
                if cds_transcripts:
                    rep = max(cds_transcripts, key=lambda t: sum(c.end - c.start for c in t.cds))
                else:
                    rep = max(gene.transcripts, key=lambda t: sum(e.end - e.start for e in t.exons))
                transcripts_to_show.append((gene, rep))
                genes_to_label.append(gene)
        else:
            # Show all transcripts
            for transcript in gene.transcripts:
                transcripts_to_show.append((gene, transcript))
            if mode == "all_transcripts":
                genes_to_label.append(gene)

    if not transcripts_to_show:
        fig = go.Figure()
        fig.add_annotation(
            text=f"No transcripts found in {region}",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig

    # Create figure
    fig = go.Figure()

    # Calculate y-axis range
    n_rows = len(transcripts_to_show)
    y_max = n_rows * (1 + row_spacing)

    # Add shapes for each transcript
    for i, (gene, transcript) in enumerate(transcripts_to_show):
        y_position = i * (1 + row_spacing)

        # Add transcript shapes
        shapes = _get_transcript_row(
            transcript,
            y_position=y_position,
            utr_height=utr_height,
            exon_height=exon_height,
            cds_height=cds_height
        )
        for shape in shapes:
            fig.add_shape(shape)

        # Add label
        if show_labels:
            if mode == "all_transcripts":
                label = _get_transcript_label(transcript)
            else:
                label = _get_gene_label(gene)

            fig.add_annotation(
                x=transcript.start,
                y=y_position + 1 + row_spacing / 2,
                text=label,
                xref='x',
                yref='y',
                showarrow=False,
                xanchor='left',
                font=dict(size=label_font_size, color='black'),
                align='left'
            )

    # Update layout
    fig.update_layout(
        title=title if title else f"Gene Track: {region}",
        height=height,
        width=width,
        plot_bgcolor='white',
        xaxis=dict(
            title='Genomic Position (bp)',
            showline=True,
            linecolor='black',
            mirror=True,
            range=[start, end]
        ),
        yaxis=dict(
            title='Transcripts',
            showline=True,
            linecolor='black',
            mirror=True,
            range=[-0.5, y_max],
            showticklabels=False
        ),
        margin=dict(l=50, r=50, t=50, b=50)
    )

    return fig


def plot_gene_track_with_heatmap(
    annotation_file: Union[str, Path],
    region: str,
    heatmap_data: Optional[pd.DataFrame] = None,
    heatmap_title: Optional[str] = None,
    gene_track_height: int = 300,
    heatmap_height: int = 400,
    **kwargs
) -> go.Figure:
    """
    Plot gene track below a heatmap (e.g., Hi-C contact matrix).

    Parameters
    ----------
    annotation_file : str or Path
        Path to GTF or GFF3 annotation file.
    region : str
        Genomic region to display.
    heatmap_data : pd.DataFrame, optional
        Heatmap data (e.g., contact matrix). If None, only gene track is shown.
    heatmap_title : str, optional
        Title for heatmap subplot.
    gene_track_height : int, optional
        Height of gene track subplot. Default is 300.
    heatmap_height : int, optional
        Height of heatmap subplot. Default is 400.
    **kwargs
        Additional arguments passed to plot_gene_track.

    Returns
    -------
    plotly.graph_objects.Figure
        Combined figure with heatmap and gene track.
    """
    from plotly.subplots import make_subplots

    # Create subplots
    if heatmap_data is not None:
        fig = make_subplots(
            rows=2,
            cols=1,
            row_heights=[0.6, 0.4],
            vertical_spacing=0.02,
            shared_xaxes=True
        )

        # Add heatmap
        fig.add_trace(
            go.Heatmap(
                z=heatmap_data.values,
                x=heatmap_data.columns,
                y=heatmap_data.index,
                colorscale='RdBu_r'
            ),
            row=1, col=1
        )
        fig.update_xaxes(title_text='', row=1, col=1)
        fig.update_yaxes(title_text='Position', row=1, col=1)
    else:
        fig = make_subplots(
            rows=1,
            cols=1
        )

    # Add gene track
    gene_fig = plot_gene_track(
        annotation_file=annotation_file,
        region=region,
        **kwargs
    )

    # Add gene track shapes and annotations
    row = 2 if heatmap_data is not None else 1
    for shape in gene_fig.layout.shapes:
        fig.add_shape(shape, row=row, col=1)
    for annotation in gene_fig.layout.annotations:
        fig.add_annotation(annotation, row=row, col=1)

    fig.update_xaxes(title_text='Genomic Position (bp)', row=row, col=1)
    fig.update_yaxes(showticklabels=False, row=row, col=1)

    fig.update_layout(
        title=heatmap_title if heatmap_title else f"Gene Track: {region}",
        height=gene_track_height + heatmap_height,
        plot_bgcolor='white'
    )

    return fig