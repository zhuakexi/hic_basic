
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

    Example:
        >>> bw_track = BigWigTrack("H3K27ac.bw", bin_size=1000)
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
        config: Optional[TrackConfig] = None
    ):
        """Initialize BigWigTrack.

        Args:
            data_path: Path to bigWig file.
            bin_size: Bin size in base pairs. Default 1000.
            agg_method: Aggregation method. Default 'mean'.
            color: Bar color. Default '#1f77b4'.
            orientation: Bar orientation. 'v' for vertical, 'h' for horizontal.
            config: Optional TrackConfig.
        """
        super().__init__(name="bigwig", config=config)
        self.data_path = Path(data_path)
        self.bin_size = bin_size
        self.agg_method = agg_method
        self.color = color
        self.orientation = orientation

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
            y_kwargs = dict(
                showline=True,
                linecolor="black",
                title="Signal"
            )
            if y_range:
                y_kwargs["range"] = y_range
            fig.update_yaxes(**y_kwargs, row=row, col=col)
            fig.update_xaxes(showline=False, row=row, col=col)


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
        trace = go.Heatmap(
            z=self._matrix,
            x=self._bins["start"],
            y=self._bins["start"],
            **merged_kwargs
        )
        fig.add_trace(trace, row=row, col=col)

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

            # Convert to bin coordinates
            if self._bins is not None:
                # Find bin indices
                bin1_mask = (
                    (self._bins["chrom"] == chrom1) &
                    (self._bins["start"] >= start1) &
                    (self._bins["end"] <= end1)
                )
                bin2_mask = (
                    (self._bins["chrom"] == chrom2) &
                    (self._bins["start"] >= start2) &
                    (self._bins["end"] <= end2)
                )

                bin1_idx = self._bins[bin1_mask].index.min()
                bin2_idx = self._bins[bin2_mask].index.min()

                # Get bin centers for shape positioning
                bin1_start = self._bins.loc[bin1_idx, "start"]
                bin1_end = self._bins.loc[bin1_idx, "end"]
                bin2_start = self._bins.loc[bin2_idx, "start"]
                bin2_end = self._bins.loc[bin2_idx, "end"]

                binsize1 = bin1_end - bin1_start
                binsize2 = bin2_end - bin2_start

                # Shape positioning
                shape = mark.get("shape", self.default_shape)
                color = mark.get("color", self.default_color)
                width = mark.get("width", self.default_width)

                if shape == "circle":
                    # Circle centered on the pixel
                    x0 = bin1_start - 0.5 * binsize1
                    y0 = bin2_start - 0.5 * binsize2
                    x1 = bin1_end - 0.5 * binsize1
                    y1 = bin2_end - 0.5 * binsize2

                    fig.add_shape(
                        type="circle",
                        x0=x0, y0=y0, x1=x1, y1=y1,
                        line=dict(color=color, width=width),
                        xref=f"x{row}",
                        yref=f"y{col}"
                    )
                elif shape == "rect":
                    # Rectangle outlining the pixel
                    fig.add_shape(
                        type="rect",
                        x0=bin1_start, y0=bin2_start,
                        x1=bin1_end, y1=bin2_end,
                        line=dict(color=color, width=width),
                        xref=f"x{row}",
                        yref=f"y{col}"
                    )
                elif shape == "line":
                    # Line connecting two points
                    fig.add_shape(
                        type="line",
                        x0=bin1_start, y0=bin2_start,
                        x1=bin1_end, y1=bin2_end,
                        line=dict(color=color, width=width),
                        xref=f"x{row}",
                        yref=f"y{col}"
                    )

                # Add label if provided
                if "label" in mark:
                    # Position label outside the shape
                    # Offset to top-right of the mark
                    center_x = (bin1_start + bin1_end) / 2
                    center_y = (bin2_start + bin2_end) / 2
                    
                    # Calculate offset based on binsize
                    offset_x = max(binsize1, binsize2) * 0.75
                    offset_y = max(binsize1, binsize2) * 0.75
                    
                    # Position label to the right and above the mark
                    label_x = center_x + offset_x
                    label_y = center_y + offset_y
                    
                    fig.add_annotation(
                        x=label_x,
                        y=label_y,
                        text=mark["label"],
                        showarrow=False,
                        font=dict(color=color, size=11),
                        xref=f"x{row}",
                        yref=f"y{col}"
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
        
        # Store y-axis max for render method (always 1.0 for single row)
        self._y_max = 1.0

        # Add gene labels
        for i, (gene, tx) in enumerate(transcripts_to_show):
            if self.show_labels and gene in genes_to_label:
                label = _get_gene_label(gene)
                # Position label above the transcript row
                y = float(i)  # Same y as transcript row
                self._annotations.append(dict(
                    x=tx.start,
                    y=y + self.exon_height + 0.05,
                    text=label,
                    showarrow=False,
                    font=dict(size=self.label_font_size),
                    xref='x',
                    yref='y'
                ))

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