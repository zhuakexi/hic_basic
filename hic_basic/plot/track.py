
"""
Track module for genomic data visualization.

This module provides classes and functions for visualizing genomic data
as tracks that can be composed into multi-panel figures.

Example usage:
    >>> from hic_basic.plot.track import EigenTrack, ISTrack, BigWigTrack, GeneTrack, TrackComposer
    >>> composer = TrackComposer.for_trapezoid_with_tracks(n_tracks=3)
    >>> composer.add_track(trapezoid, row=1, col=1)
    >>> composer.add_track(bw_track, row=2, col=1)
    >>> composer.add_track(gene_track, row=3, col=1)
    >>> fig = composer.build(region="chr1:1000000-5000000")
"""

import ast
import gzip
import io
import json
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Union, Tuple, Any

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
from cooler import Cooler
from .gene_track_parser import (
    parse_annotation_file,
    filter_genes_in_region,
    Gene,
    Transcript,
    Interval
)


# =============================================================================
# Track Configuration
# =============================================================================


@dataclass
class TrackConfig:
    """Configuration for a track's position and styling.

    Attributes:
        row: Row position in the subplot grid (1-indexed).
        col: Column position in the subplot grid (1-indexed).
        height_ratio: Relative height ratio for the track's row.
        title: Optional title for the track.
        y_range: Optional (min, max) tuple for y-axis range.
        style_kwargs: Additional styling keyword arguments.
    """
    row: int
    col: int
    height_ratio: float = 1.0
    title: Optional[str] = None
    y_range: Optional[Tuple[float, float]] = None
    style_kwargs: Dict[str, Any] = field(default_factory=dict)


# =============================================================================
# Track Base Class
# =============================================================================


class Track(ABC):
    """Abstract base class for genomic tracks.

    Tracks are modular visualization components that can be composed into
    multi-panel figures. Each track is responsible for loading its own data
    and rendering itself onto a plotly figure.

    Subclasses must implement:
        - load_data(): Load and prepare data for visualization
        - render(): Draw the track onto a figure at specified position

    Example:
        >>> track = EigenTrack("sample.vecs")
        >>> track.load_data(clr, "chr1:1000000-5000000")
        >>> track.render(fig, row=2, col=1)
    """

    def __init__(self, name: str, config: Optional[TrackConfig] = None):
        """Initialize a Track.

        Args:
            name: Name identifier for the track.
            config: Optional TrackConfig for position and styling.
        """
        self.name = name
        self.config = config
        self._data: Optional[pd.DataFrame] = None
        self._clr: Optional[Cooler] = None

    @abstractmethod
    def load_data(self, **kwargs):
        """Load and prepare data for the track.

        Args:
            **kwargs: Data loading parameters (e.g., cooler, region).

        Returns:
            self: Returns self for method chaining.
        """
        pass

    @abstractmethod
    def render(self, fig, row: int, col: int, **kwargs):
        """Render the track onto a figure.

        Args:
            fig: Plotly figure to add traces to.
            row: Row position in subplot grid.
            col: Column position in subplot grid.
            **kwargs: Additional rendering parameters.

        Returns:
            fig: Modified figure with track added.
        """
        pass

    def get_default_config(self) -> TrackConfig:
        """Get default configuration for this track type.

        Returns:
            TrackConfig: Default configuration.
        """
        return TrackConfig(row=1, col=1, height_ratio=1.0)

    @property
    def data(self) -> Optional[pd.DataFrame]:
        """Get the track's data."""
        return self._data


# =============================================================================
# EigenTrack - Eigenvector/Compartment Track
# =============================================================================


class EigenTrack(Track):
    """Eigenvector track for A/B compartment visualization.

    Displays eigenvector values as bar plots, colored by compartment
    assignment (A positive, B negative).

    Attributes:
        data_path: Path to eigenvector file (bedGraph or TSV format).
        eigen_col: Column name containing eigenvector values.
        cmap: Color mapping for A/B compartments.
        orientation: Bar orientation ('v' for vertical, 'h' for horizontal).

    Example:
        >>> eigen = EigenTrack("sample.vecs", eigen_col="E1")
        >>> eigen.load_data(clr, "chr1:1000000-5000000")
        >>> eigen.render(fig, row=2, col=1)
    """

    def __init__(
        self,
        data_path: Union[str, Path],
        eigen_col: str = "E1",
        cmap: Optional[Dict[str, str]] = None,
        orientation: str = "v",
        config: Optional[TrackConfig] = None
    ):
        """Initialize EigenTrack.

        Args:
            data_path: Path to eigenvector file.
            eigen_col: Column name for eigenvector values. Default 'E1'.
            cmap: Color map for compartments. Default {'A': 'red', 'B': 'blue'}.
            orientation: Bar orientation. 'v' for vertical, 'h' for horizontal.
            config: Optional TrackConfig.
        """
        super().__init__(name="eigen", config=config)
        self.data_path = Path(data_path)
        self.eigen_col = eigen_col
        self.cmap = cmap or {"A": "rgb(178,24,43)", "B": "rgb(33,102,172)"}
        self.orientation = orientation

    def load_data(self, clr: Cooler, region: str) -> "EigenTrack":
        """Load eigenvector data aligned with cooler bins.

        Args:
            clr: Cooler object for bin information.
            region: Genomic region (e.g., "chr1:1000000-5000000").

        Returns:
            self: Returns self for method chaining.
        """
        self._clr = clr
        bin_s, bin_e = clr.extent(region)
        bins = clr.bins()[bin_s:bin_e]

        eig = pd.read_table(self.data_path)
        self._data = pd.merge(
            left=bins,
            right=eig,
            how="inner",
            on=("chrom", "start", "end")
        ).assign(
            AB=lambda df: np.where(df[self.eigen_col] < 0, "B", "A")
        )
        return self

    def render(self, fig, row: int, col: int, **kwargs) -> go.Figure:
        """Render eigenvector track as bar plot.

        Args:
            fig: Plotly figure to add traces to.
            row: Row position in subplot grid.
            col: Column position in subplot grid.
            **kwargs: Additional parameters (width, showlegend, etc.).

        Returns:
            fig: Modified figure with track added.
        """
        if self._data is None:
            raise ValueError("Data not loaded. Call load_data() first.")

        bar_width = self._data["start"].diff().mode()[0]
        default_kwargs = dict(
            showlegend=False,
            width=bar_width,
            marker=dict(line=dict(width=0, color="purple"))
        )
        default_kwargs.update(kwargs)

        if self.orientation == "h":
            trace = go.Bar(
                y=self._data["start"],
                x=self._data[self.eigen_col],
                marker_color=self._data["AB"].map(self.cmap),
                orientation="h",
                **default_kwargs
            )
        else:
            trace = go.Bar(
                x=self._data["start"],
                y=self._data[self.eigen_col],
                marker_color=self._data["AB"].map(self.cmap),
                orientation="v",
                **default_kwargs
            )

        fig.add_trace(trace, row=row, col=col)
        self._update_axes(fig, row, col, **kwargs)
        return fig

    def _update_axes(self, fig, row: int, col: int, **kwargs):
        """Update axis styling for the track.

        Args:
            fig: Plotly figure.
            row: Row position.
            col: Column position.
            **kwargs: Additional axis parameters.
        """
        y_range = kwargs.get("y_range", [-0.5, 0.5])
        y_kwargs = dict(
            showline=True,
            linecolor="black",
            tickvals=[-0.5, 0.5],
            ticks="",
            range=y_range
        )

        if self.orientation == "h":
            fig.update_xaxes(**y_kwargs, row=row, col=col)
            fig.update_yaxes(showline=False, row=row, col=col)
        else:
            fig.update_yaxes(**y_kwargs, row=row, col=col)
            fig.update_xaxes(showline=False, row=row, col=col)


# =============================================================================
# ArcTrack - BEDPE links as arc connections
# =============================================================================


BEDPE_COORD_COLUMNS = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]


def _parse_track_region(region: str) -> Tuple[str, float, float]:
    """Parse a simple genomic region string into chrom, start, end."""
    if ":" not in region:
        return region, 0.0, float("inf")
    chrom, coords = region.split(":", 1)
    if "-" not in coords:
        return chrom, 0.0, float(coords.replace(",", ""))
    start, end = coords.replace(",", "").split("-", 1)
    return chrom, float(start), float(end)


def _read_bedpe_like(data: Union[str, Path, pd.DataFrame]) -> pd.DataFrame:
    """Read a BEDPE-like table with optional header and score/info columns."""
    if isinstance(data, pd.DataFrame):
        df = data.copy()
    else:
        data_path = Path(data)
        opener = gzip.open if data_path.suffix == ".gz" else open
        with opener(data_path, "rt") as handle:
            content = "".join(line for line in handle if not line.lstrip().startswith("#"))
        df = pd.read_table(io.StringIO(content), header=None)
        if df.empty:
            return pd.DataFrame(columns=BEDPE_COORD_COLUMNS)

        first_row = [str(value).strip() for value in df.iloc[0, :6].tolist()]
        has_header = first_row == BEDPE_COORD_COLUMNS
        if has_header:
            header = [str(value).strip() for value in df.iloc[0].tolist()]
            df = df.iloc[1:].reset_index(drop=True)
            df.columns = header
        else:
            names = BEDPE_COORD_COLUMNS.copy()
            if df.shape[1] >= 7:
                names.append("score")
            if df.shape[1] >= 8:
                names.append("info")
            names.extend([f"extra_{i}" for i in range(df.shape[1] - len(names))])
            df.columns = names

    if set(BEDPE_COORD_COLUMNS).issubset(df.columns):
        return df

    if all(col in df.columns for col in range(6)):
        rename = {idx: name for idx, name in enumerate(BEDPE_COORD_COLUMNS)}
        if 6 in df.columns:
            rename[6] = "score"
        if 7 in df.columns:
            rename[7] = "info"
        return df.rename(columns=rename)

    missing = set(BEDPE_COORD_COLUMNS) - set(df.columns)
    raise ValueError(f"BEDPE data is missing required columns: {sorted(missing)}")


def _parse_line_info(value: Any) -> Dict[str, Any]:
    """Parse a Plotly line style dict from JSON or Python literal text."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return {}
    if isinstance(value, dict):
        return dict(value)
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return {}
    for parser in (json.loads, ast.literal_eval):
        try:
            parsed = parser(text)
            if isinstance(parsed, dict):
                return parsed
        except (ValueError, SyntaxError, TypeError, json.JSONDecodeError):
            continue
    raise ValueError(f"info must be a dict/JSON/Python-literal line style, got: {value!r}")


class ArcTrack(Track):
    """Arc track for BEDPE connections on a genomic axis.

    The track reads BEDPE coordinates and renders each intra-region pair as a
    Plotly line arc. The ``score`` column controls default line width after
    min/max scaling. The ``info`` column may contain a JSON or Python literal
    dictionary of Plotly ``scatter.line`` attributes, for example
    ``{"color": "red", "width": 4, "dash": "dot"}``; these attributes
    override score-derived defaults.

    Parameters
    ----------
    data : str, Path, or pandas.DataFrame
        BEDPE-like table. Required columns are chrom1, start1, end1, chrom2,
        start2, and end2. Optional columns score and info are recognized.
    score_col : str, optional
        Column used to scale arc line width. Default is "score".
    info_col : str, optional
        Column with Plotly line style dictionaries. Default is "info".
    line : dict, optional
        Default Plotly scatter line attributes.
    score_width_range : tuple of float, optional
        Width range used when scaling score. Default is (1.0, 6.0).
    n_points : int, optional
        Number of points sampled per arc. Default is 60.
    arc_height_scale : float, optional
        Multiplier for semicircle arc height. Default is 1.0.
    display_name : str, optional
        Track name shown in hover/legend context. Default is "arc".
    """

    def __init__(
        self,
        data: Union[str, Path, pd.DataFrame],
        score_col: str = "score",
        info_col: str = "info",
        line: Optional[Dict[str, Any]] = None,
        score_width_range: Tuple[float, float] = (1.0, 6.0),
        n_points: int = 60,
        arc_height_scale: float = 1.0,
        display_name: str = "arc",
        config: Optional[TrackConfig] = None,
    ):
        super().__init__(name=display_name, config=config)
        self.data_source = data
        self.score_col = score_col
        self.info_col = info_col
        self.line = line or {"color": "#4c78a8", "width": 1.5}
        self.score_width_range = score_width_range
        self.n_points = n_points
        self.arc_height_scale = arc_height_scale
        self.display_name = display_name
        self._region: Optional[Tuple[str, float, float]] = None
        self._max_arc_height: float = 1.0

    def load_data(self, region: str, **kwargs) -> "ArcTrack":
        """Load and filter BEDPE arcs for a genomic region.

        Args:
            region: Genomic region such as "chr1:0-5000000".
            **kwargs: Ignored, accepted for TrackComposer compatibility.

        Returns:
            self: Returns self for method chaining.
        """
        chrom, region_start, region_end = _parse_track_region(region)
        df = _read_bedpe_like(self.data_source)
        if df.empty:
            self._data = df.copy()
            self._region = (chrom, region_start, region_end)
            self._max_arc_height = 1.0
            return self

        for coord_col in ["start1", "end1", "start2", "end2"]:
            df[coord_col] = pd.to_numeric(df[coord_col], errors="coerce")
        df = df.dropna(subset=["start1", "end1", "start2", "end2"]).copy()

        mask = (
            (df["chrom1"].astype(str) == str(chrom))
            & (df["chrom2"].astype(str) == str(chrom))
            & (df["end1"] >= region_start)
            & (df["start1"] <= region_end)
            & (df["end2"] >= region_start)
            & (df["start2"] <= region_end)
        )
        df = df.loc[mask].reset_index(drop=True)
        if df.empty:
            self._data = df
            self._region = (chrom, region_start, region_end)
            self._max_arc_height = 1.0
            return self

        df["x1"] = (df["start1"] + df["end1"]) / 2
        df["x2"] = (df["start2"] + df["end2"]) / 2
        df["span"] = (df["x2"] - df["x1"]).abs()
        df = df[df["span"] > 0].reset_index(drop=True)
        if df.empty:
            self._data = df
            self._region = (chrom, region_start, region_end)
            self._max_arc_height = 1.0
            return self

        if self.score_col in df.columns:
            score = pd.to_numeric(df[self.score_col], errors="coerce")
            min_width, max_width = self.score_width_range
            finite_score = score[np.isfinite(score)]
            if len(finite_score) > 0 and finite_score.max() > finite_score.min():
                width = min_width + (score - finite_score.min()) / (finite_score.max() - finite_score.min()) * (max_width - min_width)
                df["_score_width"] = width.fillna(min_width)
            elif len(finite_score) > 0:
                df["_score_width"] = (min_width + max_width) / 2
            else:
                df["_score_width"] = self.line.get("width", min_width)
        else:
            df["_score_width"] = self.line.get("width", self.score_width_range[0])

        if self.info_col in df.columns:
            df["_line_info"] = df[self.info_col].apply(_parse_line_info)
        else:
            df["_line_info"] = [{} for _ in range(len(df))]

        self._max_arc_height = max(float(df["span"].max()) * self.arc_height_scale / 2, 1.0)
        self._data = df
        self._region = (chrom, region_start, region_end)
        return self

    def render(self, fig, row: int, col: int, **kwargs) -> go.Figure:
        """Render BEDPE links as arcs."""
        if self._data is None:
            raise ValueError("Data not loaded. Call load_data() first.")

        if self._data.empty:
            fig.update_yaxes(showticklabels=False, range=[0, 1], row=row, col=col)
            return fig

        theta = np.linspace(0, np.pi, self.n_points)
        custom_cols = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]
        if self.score_col in self._data.columns:
            custom_cols.append(self.score_col)

        for idx, arc in self._data.iterrows():
            x_left = min(float(arc["x1"]), float(arc["x2"]))
            x_right = max(float(arc["x1"]), float(arc["x2"]))
            half_span = (x_right - x_left) / 2
            center = (x_left + x_right) / 2
            height = half_span * self.arc_height_scale
            x = center - half_span * np.cos(theta)
            y = height * np.sin(theta)

            line = dict(self.line)
            line["width"] = float(arc["_score_width"])
            line.update(arc["_line_info"])

            customdata = np.tile(arc[custom_cols].to_numpy(dtype=object), (len(x), 1))
            score_hover = ""
            if self.score_col in self._data.columns:
                score_hover = f"<br>{self.score_col}=%{{customdata[6]}}"

            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode="lines",
                    line=line,
                    name=self.display_name,
                    showlegend=False,
                    customdata=customdata,
                    hovertemplate=(
                        "arc "
                        "%{customdata[0]}:%{customdata[1]}-%{customdata[2]} -> "
                        "%{customdata[3]}:%{customdata[4]}-%{customdata[5]}"
                        + score_hover
                        + "<extra></extra>"
                    ),
                ),
                row=row,
                col=col,
            )

        if self._region is not None:
            _, region_start, region_end = self._region
            if np.isfinite(region_end):
                fig.update_xaxes(range=[region_start, region_end], row=row, col=col)

        fig.update_xaxes(showline=True, linecolor="black", row=row, col=col)
        fig.update_yaxes(
            showticklabels=False,
            showgrid=False,
            zeroline=False,
            range=[0, self._max_arc_height * 1.08],
            row=row,
            col=col,
        )
        return fig


def plot_arc_track(
    data: Union[str, Path, pd.DataFrame],
    region: str,
    title: Optional[str] = None,
    width: int = 1000,
    height: int = 260,
    output_path: Optional[Union[str, Path]] = None,
    **arc_kwargs,
) -> go.Figure:
    """Plot BEDPE links as a standalone genomic arc track.

    Parameters
    ----------
    data : str, Path, or pandas.DataFrame
        BEDPE-like table accepted by :class:`ArcTrack`. The first six columns
        are chrom1, start1, end1, chrom2, start2, and end2. Optional score
        controls default arc width, and optional info supplies Plotly line
        attributes as a dict/JSON/Python-literal string.
    region : str
        Genomic region to show, for example "chr1:0-5000000".
    title : str, optional
        Figure title. If None, a title is generated from region.
    width, height : int, optional
        Figure size in pixels. Defaults are 1000 by 260.
    output_path : str or Path, optional
        PNG output path. If provided, the figure is written to this path.
    **arc_kwargs
        Additional keyword arguments passed to :class:`ArcTrack`.

    Returns
    -------
    plotly.graph_objects.Figure
        The generated arc track figure.
    """
    from plotly.subplots import make_subplots

    track = ArcTrack(data, **arc_kwargs)
    track.load_data(region)
    fig = make_subplots(rows=1, cols=1)
    track.render(fig, row=1, col=1)
    fig.update_layout(
        title=title or f"Arc track: {region}",
        width=width,
        height=height,
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(l=50, r=20, t=45, b=35),
    )
    if output_path is not None:
        fig.write_image(str(output_path), width=width, height=height)
    return fig


# =============================================================================
# BigWigTrack - Continuous Signal Track from BigWig Files
# =============================================================================


class BigWigTrack(Track):
    """BigWig signal track for continuous genomic data visualization.

    Displays signal from bigWig files as a bar plot.

    Attributes:
        data_path: Path to bigWig file.
        bin_size: Bin size for aggregating signal.
        agg_method: Aggregation method ('mean', 'max', 'min', 'sum').
        color: Bar color.
        orientation: Bar orientation ('v' for vertical, 'h' for horizontal).
        display_name: Name to display on the track label.

    Example:
        >>> bw_track = BigWigTrack("H3K27ac.bw", bin_size=1000, display_name="H3K27ac")
        >>> bw_track.load_data("chr1:1000000-5000000")
        >>> bw_track.render(fig, row=2, col=1)
    """

    def __init__(
        self,
        data_path: Union[str, Path],
        bin_size: int = 1000,
        agg_method: str = "mean",
        color: str = "#1f77b4",
        orientation: str = "v",
        display_name: Optional[str] = None,
        config: Optional[TrackConfig] = None
    ):
        """Initialize BigWigTrack.

        Args:
            data_path: Path to bigWig file.
            bin_size: Bin size in base pairs. Default 1000.
            agg_method: Aggregation method. Default 'mean'.
            color: Bar color. Default '#1f77b4'.
            orientation: Bar orientation. 'v' for vertical, 'h' for horizontal.
            display_name: Name to display on track label. Default is filename.
            config: Optional TrackConfig.
        """
        # Extract filename for display name if not provided
        if display_name is None:
            display_name = Path(data_path).name
        super().__init__(name=display_name, config=config)
        self.data_path = Path(data_path)
        self.bin_size = bin_size
        self.agg_method = agg_method
        self.color = color
        self.orientation = orientation
        self.display_name = display_name

    def load_data(self, region: str, chrom_sizes: Optional[Dict[str, int]] = None) -> "BigWigTrack":
        """Load bigWig signal data for a genomic region.

        Args:
            region: Genomic region (e.g., "chr1:1000000-5000000").
            chrom_sizes: Optional dict of chromosome sizes. If not provided,
                will attempt to read from bigWig file.

        Returns:
            self: Returns self for method chaining.

        Raises:
            ImportError: If pyBigWig is not installed.
        """
        try:
            import pyBigWig
        except ImportError:
            raise ImportError(
                "pyBigWig is required for BigWigTrack. "
                "Install with: pip install pyBigWig"
            )

        # Parse region
        if ':' in region:
            chrom, coords = region.split(':', 1)
            if '-' in coords:
                start, end = map(int, coords.replace(',', '').split('-'))
            else:
                start, end = 0, int(coords.replace(',', ''))
        else:
            raise ValueError("Region must be in format 'chrom:start-end'")

        # Open bigWig file and read signal
        with pyBigWig.open(str(self.data_path)) as bw:
            # Get chromosome size if not provided
            if chrom_sizes is None:
                chrom_sizes = bw.chroms()

            if chrom not in chrom_sizes:
                raise ValueError(f"Chromosome {chrom} not found in bigWig file")

            chrom_size = chrom_sizes[chrom]
            end = min(end, chrom_size)

            # Read signal in bins
            bin_starts = list(range(start, end, self.bin_size))
            bin_ends = [min(e, end) for e in range(start + self.bin_size, end + self.bin_size, self.bin_size)]

            signals = []
            positions = []

            for bin_start, bin_end in zip(bin_starts, bin_ends):
                try:
                    values = bw.values(chrom, bin_start, bin_end)
                    # Filter out NaN values
                    values = [v for v in values if not np.isnan(v)]
                    if values:
                        if self.agg_method == "mean":
                            sig = np.mean(values)
                        elif self.agg_method == "max":
                            sig = np.max(values)
                        elif self.agg_method == "min":
                            sig = np.min(values)
                        elif self.agg_method == "sum":
                            sig = np.sum(values)
                        else:
                            sig = np.mean(values)
                        signals.append(sig)
                        positions.append(bin_start)
                except Exception:
                    # Skip bins with no coverage
                    signals.append(0)
                    positions.append(bin_start)

            self._data = pd.DataFrame({
                "start": positions,
                "signal": signals
            })

        return self

    def render(self, fig, row: int, col: int, **kwargs) -> go.Figure:
        """Render bigWig signal as bar plot.

        Args:
            fig: Plotly figure to add traces to.
            row: Row position in subplot grid.
            col: Column position in subplot grid.
            **kwargs: Additional parameters (width, showlegend, etc.).

        Returns:
            fig: Modified figure with track added.
        """
        if self._data is None:
            raise ValueError("Data not loaded. Call load_data() first.")

        bar_width = self.bin_size
        default_kwargs = dict(
            showlegend=False,
            width=bar_width,
            marker=dict(color=self.color, line=dict(width=0))
        )
        default_kwargs.update(kwargs)

        if self.orientation == "h":
            trace = go.Bar(
                y=self._data["start"],
                x=self._data["signal"],
                orientation="h",
                **default_kwargs
            )
        else:
            trace = go.Bar(
                x=self._data["start"],
                y=self._data["signal"],
                orientation="v",
                **default_kwargs
            )

        fig.add_trace(trace, row=row, col=col)
        self._update_axes(fig, row, col, **kwargs)
        return fig

    def _update_axes(self, fig, row: int, col: int, **kwargs):
        """Update axis styling for bigWig track.

        Args:
            fig: Plotly figure.
            row: Row position.
            col: Column position.
            **kwargs: Additional axis parameters.
        """
        y_range = kwargs.get("y_range", None)

        if self.orientation == "h":
            x_kwargs = dict(
                showline=True,
                linecolor="black",
                title="Signal"
            )
            if y_range:
                x_kwargs["range"] = y_range
            fig.update_xaxes(**x_kwargs, row=row, col=col)
            fig.update_yaxes(showline=False, row=row, col=col)
        else:
            # Use empty title and add annotation for horizontal track label
            y_kwargs = dict(
                showline=True,
                linecolor="black",
                title="",  # Empty title
                showticklabels=True,
                tickfont=dict(size=9)
            )
            if y_range:
                y_kwargs["range"] = y_range
            fig.update_yaxes(**y_kwargs, row=row, col=col)
            fig.update_xaxes(showline=False, row=row, col=col)
            
            # Add horizontal track name annotation at left side, vertically centered
            # Use domain coordinates for y axis (0-1 within the subplot)
            yref_domain = f"y{row} domain" if row > 1 else "y domain"
            
            fig.add_annotation(
                text=self.display_name,
                xref="paper",  # Use paper coordinates for x (0-1 across entire figure)
                yref=yref_domain,
                x=0.02,  # 2% from left edge of figure
                y=0.5,  # Vertically centered in track (domain coordinates)
                showarrow=False,
                textangle=0,  # Horizontal text
                xanchor="left",
                yanchor="middle",
                font=dict(size=9, color="black"),
            )


# =============================================================================
# ISTrack - Insulation Score Track
# =============================================================================


class ISTrack(Track):
    """Insulation Score track for TAD boundary visualization.

    Displays insulation score as a line plot.

    Attributes:
        data_path: Path to insulation score file.
        is_col: Column name containing insulation score values.
        line_color: Color of the line plot.
        orientation: Line orientation ('v' for vertical, 'h' for horizontal).

    Example:
        >>> is_track = ISTrack("sample.IS.tsv", is_col="IS_score")
        >>> is_track.load_data(clr, "chr1:1000000-5000000")
        >>> is_track.render(fig, row=3, col=1)
    """

    def __init__(
        self,
        data_path: Union[str, Path],
        is_col: str = "IS",
        line_color: str = "black",
        orientation: str = "v",
        config: Optional[TrackConfig] = None
    ):
        """Initialize ISTrack.

        Args:
            data_path: Path to insulation score file.
            is_col: Column name for insulation score. Default 'IS'.
            line_color: Line color. Default 'black'.
            orientation: Line orientation. 'v' for vertical, 'h' for horizontal.
            config: Optional TrackConfig.
        """
        super().__init__(name="IS", config=config)
        self.data_path = Path(data_path)
        self.is_col = is_col
        self.line_color = line_color
        self.orientation = orientation

    def load_data(self, clr: Cooler, region: str) -> "ISTrack":
        """Load insulation score data aligned with cooler bins.

        Args:
            clr: Cooler object for bin information.
            region: Genomic region.

        Returns:
            self: Returns self for method chaining.
        """
        self._clr = clr
        bin_s, bin_e = clr.extent(region)
        bins = clr.bins()[bin_s:bin_e]

        is_data = pd.read_table(self.data_path)
        self._data = pd.merge(
            left=bins,
            right=is_data,
            how="inner",
            on=("chrom", "start", "end")
        )
        return self

    def render(self, fig, row: int, col: int, **kwargs) -> go.Figure:
        """Render insulation score track as line plot.

        Args:
            fig: Plotly figure to add traces to.
            row: Row position in subplot grid.
            col: Column position in subplot grid.
            **kwargs: Additional parameters (line width, etc.).

        Returns:
            fig: Modified figure with track added.
        """
        if self._data is None:
            raise ValueError("Data not loaded. Call load_data() first.")

        default_kwargs = dict(
            mode="lines",
            line=dict(color=self.line_color, width=1)
        )
        default_kwargs.update(kwargs)

        if self.orientation == "h":
            trace = go.Scatter(
                y=self._data["start"],
                x=self._data[self.is_col],
                **default_kwargs
            )
        else:
            trace = go.Scatter(
                x=self._data["start"],
                y=self._data[self.is_col],
                **default_kwargs
            )

        fig.add_trace(trace, row=row, col=col)
        self._update_axes(fig, row, col, **kwargs)
        return fig

    def _update_axes(self, fig, row: int, col: int, **kwargs):
        """Update axis styling for the track.

        Args:
            fig: Plotly figure.
            row: Row position.
            col: Column position.
            **kwargs: Additional axis parameters.
        """
        y_kwargs = dict(
            showline=True,
            linecolor="black",
            ticks=""
        )

        if self.orientation == "h":
            y_kwargs["side"] = "top"
            fig.update_xaxes(**y_kwargs, row=row, col=col)
        else:
            fig.update_yaxes(**y_kwargs, row=row, col=col)


# =============================================================================
# HeatmapTrack - Hi-C Contact Matrix Heatmap
# =============================================================================


class HeatmapTrack(Track):
    """Heatmap track for Hi-C contact matrix visualization.

    Displays contact frequency matrix as a heatmap.

    Attributes:
        cool_path: Path to cooler file.
        balance: Whether to use balanced matrix.
        marks: List of marks to display on the heatmap.
        heatmap_kwargs: Additional heatmap parameters.

    Example:
        >>> heatmap = HeatmapTrack("sample.cool", balance=True)
        >>> heatmap.load_data("chr1:1000000-5000000")
        >>> heatmap.render(fig, row=1, col=1)
    """

    def __init__(
        self,
        cool_path: Union[str, Path],
        balance: bool = True,
        marks: Optional[List[Dict]] = None,
        config: Optional[TrackConfig] = None,
        **heatmap_kwargs
    ):
        """Initialize HeatmapTrack.

        Args:
            cool_path: Path to cooler file.
            balance: Whether to use balanced matrix. Default True.
            marks: List of mark dictionaries. Each dict should contain:
                - chrom1, start1, end1: First genomic position
                - chrom2, start2, end2: Second genomic position
                - color: Optional, color of the mark (default: red)
                - shape: Optional, shape type ("circle", "rect") (default: circle)
                - width: Optional, line width (default: 2)
                - label: Optional, text label for the mark
            config: Optional TrackConfig.
            **heatmap_kwargs: Additional parameters for go.Heatmap.
        """
        super().__init__(name="heatmap", config=config)
        self.cool_path = Path(cool_path)
        self.balance = balance
        self.marks = marks or []
        self.heatmap_kwargs = heatmap_kwargs
        self._matrix: Optional[np.ndarray] = None
        self._bins: Optional[pd.DataFrame] = None

    def load_data(self, region: str, clr: Optional[Cooler] = None) -> "HeatmapTrack":
        """Load contact matrix from cooler file.

        Args:
            region: Genomic region.
            clr: Optional Cooler object. If not provided, will create from cool_path.

        Returns:
            self: Returns self for method chaining.
        """
        if clr is not None:
            self._clr = clr
        else:
            self._clr = Cooler(str(self.cool_path))
        self._matrix = self._clr.matrix(balance=self.balance).fetch(region)
        self._bins = self._clr.bins().fetch(region)
        return self

    def render(self, fig, row: int, col: int, **kwargs) -> go.Figure:
        """Render contact matrix as heatmap.

        Args:
            fig: Plotly figure to add traces to.
            row: Row position in subplot grid.
            col: Column position in subplot grid.
            **kwargs: Additional parameters for go.Heatmap.

        Returns:
            fig: Modified figure with track added.
        """
        if self._matrix is None:
            raise ValueError("Data not loaded. Call load_data() first.")

        merged_kwargs = {**self.heatmap_kwargs, **kwargs}
        bin_centers = (self._bins["start"] + self._bins["end"]) / 2
        trace = go.Heatmap(
            z=self._matrix,
            x=bin_centers,
            y=bin_centers,
            **merged_kwargs
        )
        fig.add_trace(trace, row=row, col=col)
        fig.update_xaxes(
            range=[float(self._bins["start"].min()), float(self._bins["end"].max())],
            row=row,
            col=col,
        )
        fig.update_yaxes(
            range=[float(self._bins["start"].min()), float(self._bins["end"].max())],
            row=row,
            col=col,
        )

        # Render marks if present
        if self.marks:
            mark_track = MarkTrack(
                marks=self.marks,
                default_shape="circle",
                default_color="red",
                default_width=2
            )
            mark_track._clr = getattr(self, '_clr', None)
            mark_track._bins = self._bins
            mark_track.render(fig, row=row, col=col)

        return fig

    @classmethod
    def from_matrix(
        cls,
        matrix: np.ndarray,
        bins: pd.DataFrame,
        config: Optional[TrackConfig] = None,
        marks: Optional[List[Dict]] = None,
        **heatmap_kwargs
    ) -> "HeatmapTrack":
        """Create HeatmapTrack from pre-computed matrix.

        Args:
            matrix: Contact matrix as numpy array.
            bins: DataFrame with bin positions (chrom, start, end).
            config: Optional TrackConfig.
            marks: List of mark dictionaries.
            **heatmap_kwargs: Additional parameters for go.Heatmap.

        Returns:
            HeatmapTrack: New instance with pre-computed data.
        """
        track = cls.__new__(cls)
        Track.__init__(track, name="heatmap", config=config)
        track._matrix = matrix
        track._bins = bins
        track.marks = marks or []
        track.heatmap_kwargs = heatmap_kwargs
        return track


# =============================================================================
# MarkTrack - Mark Specific Pixels on Heatmap
# =============================================================================


class MarkTrack(Track):
    """Mark specific pixels or regions on a Hi-C heatmap.

    Supports marking position pairs with various shapes (circle, rectangle, line).
    Marks are rendered as outlines on top of heatmap tracks.

    Attributes:
        marks: List of dictionaries specifying position pairs and styling.
        default_shape: Default shape for marks ("circle", "rect", "line").
        default_color: Default color for marks.
        default_width: Default line width for marks.

    Example:
        >>> mark = MarkTrack(
        ...     marks=[
        ...         {"chrom1": "chr1", "start1": 100000,
        ...          "chrom2": "chr1", "start2": 200000,
        ...          "color": "red", "shape": "circle"}
        ...     ]
        ... )
        >>> mark.load_data(clr=clr, region="chr1:0-500000")
        >>> mark.render(fig, row=1, col=1)
    """

    def __init__(
        self,
        marks: Optional[List[Dict]] = None,
        default_shape: str = "circle",
        default_color: str = "red",
        default_width: int = 2,
        config: Optional[TrackConfig] = None,
    ):
        """Initialize MarkTrack.

        Args:
            marks: List of mark dictionaries. Each dict should contain:
                - chrom1, start1, end1: First genomic position
                - chrom2, start2, end2: Second genomic position
                - color: Optional, color of the mark (default: red)
                - shape: Optional, shape type ("circle", "rect") (default: circle)
                - width: Optional, line width (default: 2)
                - label: Optional, text label for the mark
            default_shape: Default shape for marks.
            default_color: Default color for marks.
            default_width: Default line width for marks.
            config: Optional TrackConfig.
        """
        super().__init__(name="mark", config=config)
        self.marks = marks or []
        self.default_shape = default_shape
        self.default_color = default_color
        self.default_width = default_width
        self._clr: Optional[Cooler] = None
        self._bins: Optional[pd.DataFrame] = None

    def load_data(
        self,
        clr: Optional[Cooler] = None,
        region: Optional[str] = None,
        bins: Optional[pd.DataFrame] = None
    ) -> "MarkTrack":
        """Load cooler data for coordinate conversion.

        Args:
            clr: Cooler object.
            region: Genomic region (used to filter marks).
            bins: Optional bins DataFrame. If not provided, will fetch from clr.

        Returns:
            self: Returns self for method chaining.
        """
        if clr is not None:
            self._clr = clr
        if bins is not None:
            self._bins = bins
        elif self._clr is not None:
            self._bins = self._clr.bins().fetch(region) if region else self._clr.bins()[:]
        return self

    def render(self, fig, row: int, col: int, **kwargs) -> go.Figure:
        """Render marks as shapes on the figure.

        Args:
            fig: Plotly figure to add shapes to.
            row: Row position in subplot grid.
            col: Column position in subplot grid.
            **kwargs: Additional parameters (ignored).

        Returns:
            fig: Modified figure with marks added.
        """
        if self._clr is None and self._bins is None:
            raise ValueError("Data not loaded. Call load_data() first.")

        for mark in self.marks:
            # Get genomic coordinates
            chrom1 = mark.get("chrom1")
            start1 = mark.get("start1")
            end1 = mark.get("end1")
            chrom2 = mark.get("chrom2")
            start2 = mark.get("start2")
            end2 = mark.get("end2")

            # Draw the provided intervals directly. This keeps multi-bin loop
            # anchors aligned with guide lines and avoids half-bin shifts from
            # treating bin starts as heatmap cell centers.
            if self._bins is not None:
                region_chroms = set(self._bins["chrom"].astype(str))
                if str(chrom1) not in region_chroms or str(chrom2) not in region_chroms:
                    continue

                x0 = float(start1)
                x1 = float(end1)
                y0 = float(start2)
                y1 = float(end2)
                if x1 <= x0 or y1 <= y0:
                    continue

                # Shape positioning
                shape = mark.get("shape", self.default_shape)
                color = mark.get("color", self.default_color)
                width = mark.get("width", self.default_width)

                if shape == "circle":
                    fig.add_shape(
                        type="circle",
                        x0=x0, y0=y0, x1=x1, y1=y1,
                        line=dict(color=color, width=width),
                        row=row,
                        col=col,
                    )
                elif shape == "rect":
                    fig.add_shape(
                        type="rect",
                        x0=x0, y0=y0,
                        x1=x1, y1=y1,
                        line=dict(color=color, width=width),
                        row=row,
                        col=col,
                    )
                elif shape == "line":
                    fig.add_shape(
                        type="line",
                        x0=x0, y0=y0,
                        x1=x1, y1=y1,
                        line=dict(color=color, width=width),
                        row=row,
                        col=col,
                    )

                # Add label if provided
                if "label" in mark:
                    # Position label outside the shape
                    # Offset to top-right of the mark
                    center_x = (x0 + x1) / 2
                    center_y = (y0 + y1) / 2
                    
                    # Calculate offset based on interval size
                    offset_x = max(x1 - x0, y1 - y0) * 0.75
                    offset_y = max(x1 - x0, y1 - y0) * 0.75
                    
                    # Position label to the right and above the mark
                    label_x = center_x + offset_x
                    label_y = center_y + offset_y
                    
                    fig.add_annotation(
                        x=label_x,
                        y=label_y,
                        text=mark["label"],
                        showarrow=False,
                        font=dict(color=color, size=11),
                        row=row,
                        col=col,
                    )

        return fig

    @classmethod
    def from_pixel_ids(
        cls,
        pixel_ids: List[int],
        clr: Cooler,
        binsize: int,
        **kwargs
    ) -> "MarkTrack":
        """Create MarkTrack from pixel IDs.

        Args:
            pixel_ids: List of pixel IDs to mark.
            clr: Cooler object for coordinate conversion.
            binsize: Bin size in base pairs.
            **kwargs: Additional arguments passed to MarkTrack.

        Returns:
            MarkTrack: New instance with marks at specified pixel IDs.
        """
        from ..genome import GenomeIdeograph

        genome = GenomeIdeograph.from_cooler(clr)
        marks = []

        for pixel_id in pixel_ids:
            # Convert pixel_id to genomic coordinates
            temp_df = pd.DataFrame({"pixel_id": [pixel_id]})
            result = genome.join_positions(temp_df, binsize)

            mark = {
                "chrom1": result["chrom1"].iloc[0],
                "start1": result["start1"].iloc[0],
                "end1": result["end1"].iloc[0],
                "chrom2": result["chrom2"].iloc[0],
                "start2": result["start2"].iloc[0],
                "end2": result["end2"].iloc[0],
            }
            marks.append(mark)

        return cls(marks=marks, **kwargs)


# =============================================================================
# TrapezoidHeatmapTrack - Rotated Hi-C Heatmap for Compartment Visualization
# =============================================================================


class TrapezoidHeatmapTrack(Track):
    """Trapezoid heatmap track for compartment visualization.

    Displays Hi-C contact matrix rotated 45 degrees to show compartment
    patterns in a trapezoidal shape.

    Attributes:
        cool_path: Path to cooler file.
        vec_file: Path to eigenvector file.
        balance: Whether to use balanced matrix.
        strata: Number of rows to keep in the trapezoid.
        cmap: Colormap name.
        quant: Quantile for color scale clipping.
        fillna: Whether to fill NaN values.

    Example:
        >>> trapezoid = TrapezoidHeatmapTrack(
        ...     "sample.cool", "eigvec.bedgraph", strata=30
        ... )
        >>> trapezoid.load_data("chr1:1000000-5000000")
        >>> trapezoid.render(fig, row=1, col=1)
    """

    def __init__(
        self,
        cool_path: Union[str, Path],
        vec_file: Union[str, Path],
        balance: bool = True,
        strata: int = 20,
        cmap: str = "RdBu_r",
        quant: float = 0.01,
        minmax: str = "min",
        eigen_col: str = "E1",
        fillna: bool = False,
        strip: bool = False,
        mask_eig_na: bool = True,
        config: Optional[TrackConfig] = None,
        **heatmap_kwargs
    ):
        """Initialize TrapezoidHeatmapTrack.

        Args:
            cool_path: Path to cooler file.
            vec_file: Path to eigenvector file.
            balance: Whether to use balanced matrix. Default True.
            strata: Number of rows in trapezoid. Default 20.
            cmap: Colormap name. Default 'RdBu_r'.
            quant: Quantile for color scale clipping. Default 0.01.
            minmax: How to compute color limits. Default 'min'.
            eigen_col: Column name for eigenvector values. Default 'E1'.
            fillna: Whether to fill NaN values. Default False.
            strip: Whether to strip consecutive zeros. Default False.
            mask_eig_na: Whether to mask NaN eigenvalues. Default True.
            config: Optional TrackConfig.
            **heatmap_kwargs: Additional parameters for go.Heatmap.
        """
        super().__init__(name="trapezoid", config=config)
        self.cool_path = Path(cool_path)
        self.vec_file = Path(vec_file)
        self.balance = balance
        self.strata = strata
        self.cmap = cmap
        self.quant = quant
        self.minmax = minmax
        self.eigen_col = eigen_col
        self.fillna = fillna
        self.strip = strip
        self.mask_eig_na = mask_eig_na
        self.heatmap_kwargs = heatmap_kwargs
        self._matrix: Optional[np.ndarray] = None
        self._bins: Optional[pd.DataFrame] = None
        self._extent: Optional[Union[float, Tuple[float, float]]] = None

    def load_data(self, region: str, clr: Optional[Cooler] = None) -> "TrapezoidHeatmapTrack":
        """Load and process trapezoid heatmap data.

        Args:
            region: Genomic region.
            clr: Optional Cooler object. If not provided, will create from cool_path.

        Returns:
            self: Returns self for method chaining.
        """
        from scipy.ndimage import rotate
        from hic_basic.compartment import compartments

        # Always use cool_path to ensure fresh file handle
        self._clr = Cooler(str(self.cool_path))

        bin_s, bin_e = self._clr.extent(region)
        bins = self._clr.bins()[bin_s:bin_e]

        # Load contact matrix and compute compartment correlation matrix
        mat = self._clr.matrix(balance=self.balance).fetch(region)
        cmat = compartments(mat, matrixonly=True)

        # Process matrix for trapezoid
        min_v = np.nanmin(cmat)
        cmat_filled = np.nan_to_num(cmat, nan=min_v - 100)
        r_cmat = rotate(cmat_filled, 45, cval=np.nan, order=1)
        r_cmat[r_cmat < min_v - 0.5] = np.nan

        # Extract trapezoid region
        df = pd.DataFrame(r_cmat)
        N = df.shape[0]
        if N % 2 == 0:
            trapezoid = df.iloc[N // 2: N // 2 + self.strata, :].values
        else:
            trapezoid = df.iloc[N // 2 + 1: N // 2 + self.strata, :].values

        # Compute extent for color scale
        extent = np.nanpercentile(np.abs(trapezoid), (1 - self.quant) * 100)

        self._matrix = trapezoid
        self._bins = bins
        self._extent = extent

        return self

    def render(self, fig, row: int, col: int, **kwargs) -> go.Figure:
        """Render trapezoid heatmap.

        Args:
            fig: Plotly figure to add traces to.
            row: Row position in subplot grid.
            col: Column position in subplot grid.
            **kwargs: Additional parameters for go.Heatmap.

        Returns:
            fig: Modified figure with track added.
        """
        if self._matrix is None:
            raise ValueError("Data not loaded. Call load_data() first.")

        merged_kwargs = {**self.heatmap_kwargs, **kwargs}

        # Set default zmin/zmax from extent
        if isinstance(self._extent, tuple):
            if "zmin" not in merged_kwargs:
                merged_kwargs["zmin"] = self._extent[0]
            if "zmax" not in merged_kwargs:
                merged_kwargs["zmax"] = self._extent[1]
        else:
            if "zmin" not in merged_kwargs:
                merged_kwargs["zmin"] = -self._extent
            if "zmax" not in merged_kwargs:
                merged_kwargs["zmax"] = self._extent

        trace = go.Heatmap(
            z=self._matrix,
            x=self._bins["start"],
            y=self._bins["start"],
            colorscale=self.cmap,
            **merged_kwargs
        )
        fig.add_trace(trace, row=row, col=col)
        return fig

    @classmethod
    def from_matrix(
        cls,
        matrix: np.ndarray,
        bins: pd.DataFrame,
        extent: Union[float, Tuple[float, float]],
        config: Optional[TrackConfig] = None,
        cmap: str = "RdBu_r",
        **heatmap_kwargs
    ) -> "TrapezoidHeatmapTrack":
        """Create TrapezoidHeatmapTrack from pre-computed matrix.

        Args:
            matrix: Trapezoid matrix as numpy array.
            bins: DataFrame with bin positions.
            extent: Color scale extent (single value or (min, max) tuple).
            config: Optional TrackConfig.
            cmap: Colormap name. Default 'RdBu_r'.
            **heatmap_kwargs: Additional parameters for go.Heatmap.

        Returns:
            TrapezoidHeatmapTrack: New instance with pre-computed data.
        """
        track = cls.__new__(cls)
        Track.__init__(track, name="trapezoid", config=config)
        track._matrix = matrix
        track._bins = bins
        track._extent = extent
        track.cmap = cmap
        track.heatmap_kwargs = heatmap_kwargs
        return track


# =============================================================================
# Genome-wide visualization functions
# =============================================================================


DEFAULT_PHASED_READ_COUNT_COLORS = {
    "ref": "#d62728",
    "alt": "#1f77b4",
    "other": "#7f7f7f",
}


def _write_plotly_png(fig, output_path: Union[str, Path], width: int = 1200, height: int = 700) -> None:
    """Write a Plotly figure to PNG, with a simple matplotlib fallback."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        fig.write_image(str(output_path), width=width, height=height)
    except Exception as exc:  # pragma: no cover - only used when kaleido/chrome is unavailable
        import matplotlib.pyplot as plt

        fallback_title = getattr(fig.layout.title, "text", None) or output_path.name
        fallback_text = f"{fallback_title}\n\nPlotly image export failed:\n{exc}"
        fallback_fig, ax = plt.subplots(figsize=(width / 150, height / 150))
        ax.axis("off")
        ax.text(0.5, 0.5, fallback_text, ha="center", va="center", wrap=True)
        fallback_fig.savefig(output_path, dpi=150, bbox_inches="tight")
        plt.close(fallback_fig)


def _normalize_region(region):
    if region is None:
        return None
    if isinstance(region, str):
        chrom, coords = region.split(":", 1)
        start, end = coords.replace(",", "").split("-", 1)
        return chrom, float(start), float(end)
    if isinstance(region, tuple) and len(region) == 3:
        chrom, start, end = region
        return chrom, float(start), float(end)
    raise ValueError("region must be None, 'chrom:start-end', or (chrom, start, end)")


def _chromosome_sort_value(chromosome):
    label = str(chromosome)
    core = label[3:] if label.lower().startswith("chr") else label
    core_lower = core.lower()
    special = {"x": 23, "y": 24, "m": 25, "mt": 25}
    if core_lower in special:
        return (0, special[core_lower], "")

    digits = []
    for char in core:
        if char.isdigit():
            digits.append(char)
        else:
            break
    if digits:
        number = int("".join(digits))
        suffix = core[len(digits):].lower()
        return (0, number, suffix)
    return (1, 0, core_lower)


def _sort_chromosome_position_table(df: pd.DataFrame) -> pd.DataFrame:
    chrom_order = sorted(
        pd.unique(df["chromosome"].astype(str)),
        key=_chromosome_sort_value,
    )
    chrom_key = pd.Categorical(df["chromosome"].astype(str), categories=chrom_order, ordered=True)
    return (
        df.assign(_chromosome_order=chrom_key)
        .sort_values(["_chromosome_order", "position"])
        .drop(columns="_chromosome_order")
        .reset_index(drop=True)
    )


def _normalize_phased_read_count_table(phased_read_counts, region=None):
    df = phased_read_counts.copy().reset_index(drop=True)
    if "chrom" in df.columns and "chromosome" not in df.columns:
        df["chromosome"] = df["chrom"]
    if "pos" in df.columns and "position" not in df.columns:
        df["position"] = df["pos"]

    missing = {"chromosome", "position", "ref", "alt", "other"} - set(df.columns)
    if missing:
        raise ValueError(f"phased_read_counts is missing required columns: {sorted(missing)}")

    if "start_pos" not in df.columns:
        df["start_pos"] = df["position"]
    if "end_pos" not in df.columns:
        df["end_pos"] = df["position"]

    for col in ["position", "start_pos", "end_pos", "ref", "alt", "other"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=["position", "start_pos", "end_pos"]).reset_index(drop=True)

    parsed_region = _normalize_region(region)
    if parsed_region is not None:
        chrom, start, end = parsed_region
        mask = (
            (df["chromosome"].astype(str) == str(chrom))
            & (df["end_pos"] >= start)
            & (df["start_pos"] <= end)
        )
        df = df.loc[mask].reset_index(drop=True)

    if len(df) == 0:
        raise ValueError("No phased_read_counts rows overlap the requested region")

    df = _sort_chromosome_position_table(df)
    single_chrom = df["chromosome"].astype(str).nunique() == 1
    if single_chrom:
        df["x_plot"] = df["position"]
        x_title = "genomic position (bp)"
    else:
        df["x_plot"] = np.arange(len(df))
        x_title = "genome-wide bin index"
    return df, single_chrom, x_title


def _add_genome_wide_chrom_boundaries(fig, df, single_chrom, annotation_y=1.08):
    if single_chrom:
        return fig

    chrom_starts = df.groupby("chromosome", sort=False)["x_plot"].min()
    chrom_ends = df.groupby("chromosome", sort=False)["x_plot"].max()
    for start in chrom_starts.iloc[1:]:
        fig.add_vline(x=float(start), line_width=1, line_color="black", opacity=0.5)
    for chrom in chrom_starts.index:
        mid = (chrom_starts[chrom] + chrom_ends[chrom]) / 2
        fig.add_annotation(
            x=mid,
            y=annotation_y,
            xref="x",
            yref="paper",
            text=str(chrom),
            showarrow=False,
            font=dict(size=18),
        )
    return fig


def plot_genome_phased_read_count_track(
    phased_read_counts,
    region=None,
    output_path: Optional[Union[str, Path]] = None,
    title: str = "phased_read_count",
    count_cols: Tuple[str, ...] = ("ref", "alt", "other"),
    colors: Optional[Dict[str, str]] = None,
    y_cap_quantile: float = 0.995,
    width: int = 1400,
    height: int = 420,
):
    """
    Plot a genome-wide phased read-count track.

    This function is intended for internal phased-read count tables rather than
    standard signal files such as bigWig. The input table must contain
    chromosome/position coordinates and phased count columns.

    Parameters
    ----------
    phased_read_counts : pandas.DataFrame
        Table with chromosome, position, ref, alt, and other columns. The
        aliases chrom and pos are accepted. Optional start_pos and end_pos
        columns are used for interval hover text and region filtering.
    region : str or tuple, optional
        Region to plot, either "chr:start-end" or (chrom, start, end). If None,
        all rows are shown with a genome-wide bin index when multiple
        chromosomes are present.
    output_path : str or Path, optional
        PNG path. If provided, the figure is written to this path.
    title : str, optional
        Figure title.
    count_cols : tuple of str, optional
        Count columns to draw. Default is ("ref", "alt", "other").
    colors : dict, optional
        Mapping from count column to trace color.
    y_cap_quantile : float, optional
        Quantile used to cap the displayed y-axis. Default 0.995.
    width, height : int, optional
        PNG export size in pixels when output_path is provided.

    Returns
    -------
    plotly.graph_objects.Figure
        The generated read-count track figure.
    """
    df, single_chrom, x_title = _normalize_phased_read_count_table(phased_read_counts, region=region)
    missing_count_cols = set(count_cols) - set(df.columns)
    if missing_count_cols:
        raise ValueError(f"phased_read_counts is missing count columns: {sorted(missing_count_cols)}")

    values = df[list(count_cols)].to_numpy(dtype=float).ravel()
    finite_values = values[np.isfinite(values)]
    if len(finite_values) == 0:
        y_cap = 1.0
    else:
        y_cap = max(float(np.quantile(finite_values, y_cap_quantile)), 1.0)
    n_clipped = int(np.sum(finite_values > y_cap))

    colors = {**DEFAULT_PHASED_READ_COUNT_COLORS, **(colors or {})}
    customdata = df[["chromosome", "start_pos", "end_pos"]].to_numpy()
    hover_base = (
        "chrom=%{customdata[0]}<br>"
        "start=%{customdata[1]}<br>"
        "end=%{customdata[2]}<br>"
    )

    fig = go.Figure()
    for col in count_cols:
        fig.add_trace(
            go.Scattergl(
                x=df["x_plot"],
                y=df[col],
                mode="lines+markers",
                line=dict(color=colors.get(col, "#333333"), width=1.6),
                marker=dict(size=5),
                name=col,
                customdata=customdata,
                hovertemplate=hover_base + f"{col}=%{{y}}<extra></extra>",
            )
        )

    _add_genome_wide_chrom_boundaries(fig, df, single_chrom)
    fig.update_xaxes(title=dict(text=x_title, font=dict(size=20)), tickfont=dict(size=18))
    fig.update_yaxes(
        title=dict(text="read count", font=dict(size=20)),
        tickfont=dict(size=18),
        range=[0, y_cap * 1.05],
    )
    fig.update_xaxes(showline=True, linewidth=1.4, linecolor="black", mirror=True, showgrid=False, zeroline=False, ticks="outside")
    fig.update_yaxes(showline=True, linewidth=1.4, linecolor="black", mirror=True, showgrid=False, zeroline=False, ticks="outside")
    cap_label = f"p{y_cap_quantile * 100:g}"
    fig.add_annotation(
        text=f"count y-axis capped at {cap_label}={y_cap:.1f}; clipped points={n_clipped}",
        x=0.5,
        y=0.94,
        xref="paper",
        yref="paper",
        showarrow=False,
        font=dict(size=16),
        bgcolor="rgba(255,255,255,0.85)",
    )
    fig.update_layout(
        title=title,
        title_font=dict(size=28),
        font=dict(size=18),
        width=width,
        height=height,
        plot_bgcolor="white",
        paper_bgcolor="white",
        legend=dict(orientation="h", yanchor="bottom", y=-0.32, font=dict(size=16)),
    )

    if output_path is not None:
        _write_plotly_png(fig, output_path, width=width, height=height)
    return fig


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


# =============================================================================
# GeneTrack - Gene Annotation Track
# =============================================================================


class GeneTrack(Track):
    """Gene annotation track for displaying gene structures from GTF/GFF3 files.

    Displays gene models including exons, introns, CDS, and UTRs.

    Attributes:
        annotation_file: Path to GTF or GFF3 annotation file.
        mode: Display mode ('all_transcripts', 'canonical', 'gene_only').
        utr_height: Relative height of UTR boxes (0-1).
        exon_height: Relative height of exon boxes (0-1).
        cds_height: Relative height of CDS boxes (0-1).
        row_spacing: Spacing between transcript rows.
        show_labels: Whether to show gene/transcript labels.
        label_font_size: Font size for labels.

    Example:
        >>> gene_track = GeneTrack("annotations.gtf", mode="all_transcripts")
        >>> gene_track.load_data("chr1:1000000-5000000")
        >>> gene_track.render(fig, row=3, col=1)
    """

    def __init__(
        self,
        annotation_file: Union[str, Path],
        mode: str = "all_transcripts",
        utr_height: float = 0.30,
        exon_height: float = 0.40,
        cds_height: float = 0.60,
        row_spacing: float = 0.2,
        show_labels: bool = True,
        label_font_size: int = 10,
        config: Optional[TrackConfig] = None
    ):
        """Initialize GeneTrack.

        Args:
            annotation_file: Path to GTF or GFF3 annotation file.
            mode: Display mode. Default 'all_transcripts'.
            utr_height: Relative height of UTR boxes. Default 0.30.
            exon_height: Relative height of exon boxes. Default 0.40.
            cds_height: Relative height of CDS boxes. Default 0.60.
            row_spacing: Spacing between transcript rows. Default 0.2.
            show_labels: Whether to show labels. Default True.
            label_font_size: Font size for labels. Default 10.
            config: Optional TrackConfig.
        """
        super().__init__(name="gene", config=config)
        self.annotation_file = Path(annotation_file)
        self.mode = mode
        self.utr_height = utr_height
        self.exon_height = exon_height
        self.cds_height = cds_height
        self.row_spacing = row_spacing
        self.show_labels = show_labels
        self.label_font_size = label_font_size
        self._genes: Optional[Dict[str, Gene]] = None
        self._shapes: List[dict] = []
        self._annotations: List[dict] = []

    def load_data(self, region: str) -> "GeneTrack":
        """Load gene annotation data for a genomic region.

        Args:
            region: Genomic region (e.g., "chr1:1000000-5000000").

        Returns:
            self: Returns self for method chaining.
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
            start, end = 0, 1000000

        # Parse annotation file
        self._genes = parse_annotation_file(self.annotation_file)

        # Filter genes in region
        self._genes = filter_genes_in_region(self._genes, seqid, start, end)

        # Generate shapes based on mode
        self._shapes = []
        self._annotations = []

        if not self._genes:
            return self

        # Collect transcripts to display
        transcripts_to_show = []
        genes_to_label = []

        for gene_id, gene in sorted(self._genes.items(), key=lambda x: x[1].start):
            if self.mode == "gene_only":
                # Select representative transcript (longest CDS)
                if gene.transcripts:
                    rep_tx = max(gene.transcripts, key=lambda t: sum(e.end - e.start for e in t.cds) if t.cds else 0)
                    transcripts_to_show.append((gene, rep_tx))
                    genes_to_label.append(gene)
            elif self.mode == "canonical":
                # Select longest transcript
                if gene.transcripts:
                    rep_tx = max(gene.transcripts, key=lambda t: sum(e.end - e.start for e in t.exons))
                    transcripts_to_show.append((gene, rep_tx))
                    genes_to_label.append(gene)
            else:  # all_transcripts
                for tx in gene.transcripts:
                    transcripts_to_show.append((gene, tx))
                genes_to_label.append(gene)

        # Generate shapes for each transcript
        # All transcripts share the same y range [0, 1]
        # Transcripts are drawn overlapping (all at y=0)
        for i, (gene, tx) in enumerate(transcripts_to_show):
            # All transcripts at y_position=0
            # Backbone at y=0.5 (center of track)
            shapes = _get_transcript_row(
                tx,
                y_position=0.0,
                utr_height=self.utr_height,
                exon_height=self.exon_height,
                cds_height=self.cds_height,
                backbone_y=0.5  # Fixed: backbone in center of track [0, 1]
            )
            self._shapes.extend(shapes)

        # Add gene labels with staggered rows to avoid overlap
        # Labels are placed above the transcript structures
        # Use multiple rows (y=1.2, 1.5, 1.8, etc.) to stagger overlapping labels
        # Start well above the transcript structures (which are at y=0 to y=1)
        label_rows = []  # List of (y_position, end_position) for staggered rows
        label_row_height = 0.25  # Height between label rows
        base_label_y = 1.25  # Starting y position for first label row (well above transcripts at y=0-1)
        
        for i, (gene, tx) in enumerate(transcripts_to_show):
            if self.show_labels and gene in genes_to_label:
                label = _get_gene_label(gene)
                # Estimate label width in base pairs
                # Use wider estimate to avoid visual overlap
                label_width_bp = len(label) * 1500  # More conservative estimate
                
                label_start = tx.start
                label_end = tx.start + label_width_bp
                
                # Find a row that doesn't overlap with existing labels
                # Start from the lowest row and go up
                row_found = False
                for row_idx, (row_y, row_end) in enumerate(label_rows):
                    # Check if this label overlaps with the last label in this row
                    if label_start > row_end + 15000:  # 15kb padding between labels in same row
                        # Can use this row
                        label_y = base_label_y + row_idx * label_row_height
                        self._annotations.append(dict(
                            x=tx.start,
                            y=label_y,
                            text=label,
                            showarrow=False,
                            font=dict(size=self.label_font_size),
                            xref='x',
                            yref='y'
                        ))
                        label_rows[row_idx] = (row_y, label_end)
                        row_found = True
                        break
                
                if not row_found:
                    # Create a new row above all existing rows
                    new_row_y = base_label_y + len(label_rows) * label_row_height
                    self._annotations.append(dict(
                        x=tx.start,
                        y=new_row_y,
                        text=label,
                        showarrow=False,
                        font=dict(size=self.label_font_size),
                        xref='x',
                        yref='y'
                    ))
                    label_rows.append((new_row_y, label_end))

        # Store y-axis max for render method
        # Account for label rows
        if label_rows:
            max_label_y = max(y for y, _ in label_rows)
            self._y_max = max(1.0, max_label_y + 0.3)  # Add padding above labels
        else:
            self._y_max = 1.0

        return self

    def render(self, fig, row: int, col: int, **kwargs) -> go.Figure:
        """Render gene track as shapes and annotations.

        Args:
            fig: Plotly figure to add traces to.
            row: Row position in subplot grid.
            col: Column position in subplot grid.
            **kwargs: Additional parameters.

        Returns:
            fig: Modified figure with track added.
        """
        if self._genes is None:
            raise ValueError("Data not loaded. Call load_data() first.")

        # Check if no genes found
        if not self._genes:
            fig.add_annotation(
                text="No genes in region",
                xref='x', yref='y',
                x=0.5, y=0.5,
                showarrow=False,
                row=row, col=col
            )
            return fig

        # Add shapes
        for shape in self._shapes:
            fig.add_shape(shape, row=row, col=col)

        # Add annotations
        for annotation in self._annotations:
            fig.add_annotation(annotation, row=row, col=col)

        # Update axes
        fig.update_xaxes(showline=True, linecolor='black', row=row, col=col)
        # Set y-axis range based on number of transcripts
        y_max = getattr(self, '_y_max', 1.0)
        # Add some padding above transcripts for labels
        fig.update_yaxes(showticklabels=False, range=[0, y_max + 0.5], row=row, col=col)

        return fig


# =============================================================================
# Helper functions for gene track visualization
# =============================================================================


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
    
    Parameters
    ----------
    y_position : float
        Bottom y position of the track row (0-1 scale).
    utr_height, exon_height, cds_height : float
        Relative heights of each feature type (0-1 scale).
    backbone_y : float, optional
        Y position for intron backbone line. If None, calculated as center
        of the track row based on max feature height.
    """
    shapes = []

    # Sort exons by start position
    sorted_exons = sorted(transcript.exons, key=lambda e: e.start)

    if not sorted_exons:
        return shapes

    # Determine backbone y position (middle of track row)
    # Track row spans from y_position to y_position + 1
    # Backbone should be in the center: y_position + 0.5
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
                # Center CDS in the track row
                y0 = y_position + (1 - height) / 2
            elif feat_type == 'utr':
                height = utr_height
                color = 'lightblue'
                # Center UTR in the track row
                y0 = y_position + (1 - height) / 2
            else:  # exon
                height = exon_height
                color = 'skyblue'
                # Center exon in the track row
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
