"""
Track Composer for assembling multi-track genomic figures.

This module provides the TrackComposer class for composing multiple tracks
into a single cohesive figure. It handles layout configuration and track
orchestration.

Example usage:
    >>> from hic_basic.plot.track import EigenTrack, HeatmapTrack, BigWigTrack, GeneTrack
    >>> from hic_basic.plot.track_composer import TrackComposer
    >>>
    >>> composer = TrackComposer.for_trapezoid_with_tracks(n_tracks=3)
    >>> composer.add_track(trapezoid, row=1, col=1)
    >>> composer.add_track(bw_track, row=2, col=1)
    >>> composer.add_track(gene_track, row=3, col=1)
    >>> fig = composer.build(region="chr1:1000000-5000000", title="Hi-C Analysis")
"""

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Dict, List, Optional, Tuple, Any

from plotly.subplots import make_subplots
from plotly.graph_objects import Figure

if TYPE_CHECKING:
    from cooler import Cooler
    from .track import Track, EigenTrack


@dataclass
class LayoutConfig:
    """Layout configuration for multi-track figures.

    Attributes:
        rows: Number of rows in the subplot grid.
        cols: Number of columns in the subplot grid.
        row_heights: Relative heights for each row. Default is equal heights.
        col_widths: Relative widths for each column. Default is equal widths.
        vertical_spacing: Spacing between rows as fraction of figure height.
        horizontal_spacing: Spacing between columns as fraction of figure width.
        shared_xaxes: Whether to share x-axes across rows.
        shared_yaxes: Whether to share y-axes across columns.
        **kwargs: Additional layout parameters passed to make_subplots.

    Example:
        >>> config = LayoutConfig(
        ...     rows=3, cols=2,
        ...     row_heights=[10, 2, 2],
        ...     col_widths=[2, 10],
        ...     vertical_spacing=0.02
        ... )
    """
    rows: int
    cols: int
    row_heights: Optional[List[float]] = None
    col_widths: Optional[List[float]] = None
    vertical_spacing: float = 0.02
    horizontal_spacing: float = 0.02
    shared_xaxes: bool = False
    shared_yaxes: bool = False
    kwargs: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        """Set default values for row_heights and col_widths."""
        if self.row_heights is None:
            self.row_heights = [1.0] * self.rows
        if self.col_widths is None:
            self.col_widths = [1.0] * self.cols


class TrackComposer:
    """Composer for assembling multiple tracks into a single figure.

    TrackComposer manages the layout and rendering of multiple tracks
    in a coordinated figure. It handles subplot creation, track positioning,
    and final figure assembly.

    Attributes:
        layout_config: LayoutConfig defining the subplot grid.
        tracks: Dictionary mapping (row, col) positions to Track objects.
        fig: The assembled plotly figure (after build()).

    Example:
        >>> composer = TrackComposer(LayoutConfig(rows=2, cols=2))
        >>> composer.add_track(heatmap_track, row=1, col=2)
        >>> composer.add_track(eigen_track, row=2, col=2)
        >>> fig = composer.build(region="chr1:1000000-5000000")
    """

    def __init__(self, layout_config: LayoutConfig):
        """Initialize TrackComposer.

        Args:
            layout_config: LayoutConfig defining the subplot grid.
        """
        self.layout_config = layout_config
        self.tracks: Dict[Tuple[int, int], "Track"] = {}
        self._fig: Optional[Figure] = None

    def add_track(
        self,
        track: "Track",
        row: int,
        col: int,
        height_ratio: Optional[float] = None
    ) -> "TrackComposer":
        """Add a track to the composer at specified position.

        Args:
            track: Track instance to add.
            row: Row position (1-indexed).
            col: Column position (1-indexed).
            height_ratio: Optional height ratio override for the track's row.

        Returns:
            self: Returns self for method chaining.

        Example:
            >>> eigen = EigenTrack("sample.vecs")
            >>> composer.add_track(eigen, row=2, col=1)
        """
        from .track import TrackConfig

        if height_ratio is not None:
            self.layout_config.row_heights[row - 1] = height_ratio

        track.config = TrackConfig(row=row, col=col, height_ratio=1.0)
        self.tracks[(row, col)] = track
        return self

    def build(
        self,
        region: str,
        title: str = "",
        clr: Optional["Cooler"] = None,
        **layout_kwargs
    ) -> Figure:
        """Build the complete figure with all tracks.

        Args:
            region: Genomic region for all tracks (e.g., "chr1:1000000-5000000").
            title: Figure title.
            clr: Optional Cooler object shared across tracks.
            **layout_kwargs: Additional layout parameters.

        Returns:
            Figure: Assembled plotly figure.

        Raises:
            ValueError: If no tracks have been added.

        Example:
            >>> fig = composer.build(
            ...     region="chr1:1000000-5000000",
            ...     title="Hi-C Analysis"
            ... )
        """
        from cooler import Cooler

        if not self.tracks:
            raise ValueError("No tracks added. Use add_track() to add tracks first.")

        # Create subplot grid
        self._fig = make_subplots(
            rows=self.layout_config.rows,
            cols=self.layout_config.cols,
            row_heights=self.layout_config.row_heights,
            column_widths=self.layout_config.col_widths,
            vertical_spacing=self.layout_config.vertical_spacing,
            horizontal_spacing=self.layout_config.horizontal_spacing,
            shared_xaxes=self.layout_config.shared_xaxes,
            shared_yaxes=self.layout_config.shared_yaxes,
            **self.layout_config.kwargs
        )

        # Load data and render each track
        for (row, col), track in self.tracks.items():
            # Use delayed import to avoid circular dependency
            from .track import ArcTrack, EigenTrack, HeatmapTrack, BigWigTrack, GeneTrack, TrapezoidHeatmapTrack

            # Skip loading if data is already loaded (pre-loaded tracks)
            if track._data is not None or (hasattr(track, '_matrix') and track._matrix is not None):
                # Data already loaded, just render
                track.render(self._fig, row=row, col=col)
                continue

            # Load data with appropriate parameters
            if isinstance(track, EigenTrack):
                # EigenTrack requires clr parameter
                if clr is not None:
                    track.load_data(clr=clr, region=region)
                else:
                    raise ValueError(
                        f"EigenTrack at ({row},{col}) requires a 'clr' parameter. "
                        "Please provide 'clr' when calling build()."
                    )
            elif isinstance(track, HeatmapTrack):
                # HeatmapTrack can use clr or cool_path
                if clr is not None:
                    track.load_data(region=region, clr=clr)
                else:
                    track.load_data(region=region)
            elif isinstance(track, TrapezoidHeatmapTrack):
                # TrapezoidHeatmapTrack can use clr or cool_path
                if clr is not None:
                    track.load_data(region=region, clr=clr)
                else:
                    track.load_data(region=region)
            elif isinstance(track, BigWigTrack):
                # BigWigTrack loads from region string
                track.load_data(region=region)
            elif isinstance(track, GeneTrack):
                # GeneTrack loads from region string
                track.load_data(region=region)
            elif isinstance(track, ArcTrack):
                # ArcTrack loads BEDPE-like links from region string
                track.load_data(region=region)
            elif clr is not None:
                # Use shared cooler if provided
                track.load_data(clr=clr, region=region)
            elif hasattr(track, 'cool_path'):
                # Track has its own cooler path
                track.load_data(region=region)
            else:
                # Generic load_data call
                track.load_data(region=region)

            # Render track to figure
            track.render(self._fig, row=row, col=col)

        # Update layout
        self._fig.update_layout(
            title=title,
            plot_bgcolor="white",
            **layout_kwargs
        )

        return self._fig

    @classmethod
    def for_tiling_compartment(
        cls,
        row_heights: Optional[List[float]] = None,
        col_widths: Optional[List[float]] = None,
        **kwargs
    ) -> "TrackComposer":
        """Create composer for tiling compartment layout.

        Creates a 2x2 grid optimized for displaying two triangular
        contact matrices with eigenvector tracks on left and bottom.

        Args:
            row_heights: Optional override for row heights. Default [10, 1].
            col_widths: Optional override for column widths. Default [1, 10].
            **kwargs: Additional LayoutConfig parameters.

        Returns:
            TrackComposer: Configured composer instance.

        Example:
            >>> composer = TrackComposer.for_tiling_compartment()
            >>> composer.add_track(heatmap, row=1, col=2)
            >>> composer.add_track(eigen_left, row=1, col=1)
            >>> composer.add_track(eigen_bottom, row=2, col=2)
        """
        return cls(LayoutConfig(
            rows=2,
            cols=2,
            row_heights=row_heights or [10, 1],
            col_widths=col_widths or [1, 10],
            vertical_spacing=0.01,
            horizontal_spacing=0.01,
            **kwargs
        ))

    @classmethod
    def for_trapezoid_2tracks(
        cls,
        row_heights: Optional[List[float]] = None,
        **kwargs
    ) -> "TrackComposer":
        """Create composer for trapezoid heatmap with two tracks.

        Creates a 3-row vertical layout for trapezoid compartment
        heatmap with two eigenvector tracks below.

        Args:
            row_heights: Optional override for row heights. Default [6, 2, 2].
            **kwargs: Additional LayoutConfig parameters.

        Returns:
            TrackComposer: Configured composer instance.

        Example:
            >>> composer = TrackComposer.for_trapezoid_2tracks()
            >>> composer.add_track(trapezoid, row=1, col=1)
            >>> composer.add_track(eigen1, row=2, col=1)
            >>> composer.add_track(eigen2, row=3, col=1)
        """
        return cls(LayoutConfig(
            rows=3,
            cols=1,
            row_heights=row_heights or [6, 2, 2],
            col_widths=[1],
            vertical_spacing=0.05,
            shared_xaxes=False,
            **kwargs
        ))

    @classmethod
    def for_trapezoid_with_tracks(
        cls,
        n_tracks: int,
        trapezoid_height: float = 6,
        track_heights: Optional[List[float]] = None,
        vertical_spacing: float = 0.05,
        **kwargs
    ) -> "TrackComposer":
        """Create composer for trapezoid heatmap with arbitrary tracks.

        Creates a vertical layout with:
        - Row 1: Trapezoid compartment heatmap
        - Rows 2 to n_tracks+1: Additional tracks (bigWig, gene, etc.)

        Args:
            n_tracks: Number of additional tracks below the trapezoid.
            trapezoid_height: Height ratio for trapezoid row. Default 6.
            track_heights: List of height ratios for each track.
                If None, all tracks have equal height of 2.
            vertical_spacing: Spacing between rows. Default 0.05.
            **kwargs: Additional LayoutConfig parameters.

        Returns:
            TrackComposer: Configured composer instance.

        Example:
            >>> composer = TrackComposer.for_trapezoid_with_tracks(
            ...     n_tracks=3,
            ...     trapezoid_height=6,
            ...     track_heights=[2, 2, 3]  # gene track taller
            ... )
            >>> composer.add_track(trapezoid, row=1, col=1)
            >>> composer.add_track(bw_track1, row=2, col=1)
            >>> composer.add_track(bw_track2, row=3, col=1)
            >>> composer.add_track(gene_track, row=4, col=1)
        """
        if track_heights is None:
            track_heights = [2.0] * n_tracks
        elif len(track_heights) != n_tracks:
            raise ValueError(
                f"track_heights length ({len(track_heights)}) must match "
                f"n_tracks ({n_tracks})"
            )

        row_heights = [trapezoid_height] + track_heights

        return cls(LayoutConfig(
            rows=1 + n_tracks,
            cols=1,
            row_heights=row_heights,
            col_widths=[1],
            vertical_spacing=vertical_spacing,
            shared_xaxes=False,
            **kwargs
        ))

    @classmethod
    def for_heatmap_with_tracks(
        cls,
        n_tracks: int,
        heatmap_height: float = 10,
        track_height: float = 2,
        **kwargs
    ) -> "TrackComposer":
        """Create composer for heatmap with multiple tracks below.

        Creates a vertical layout with heatmap on top and specified
        number of tracks below.

        Args:
            n_tracks: Number of tracks below the heatmap.
            heatmap_height: Relative height of heatmap row. Default 10.
            track_height: Relative height of each track row. Default 2.
            **kwargs: Additional LayoutConfig parameters.

        Returns:
            TrackComposer: Configured composer instance.

        Example:
            >>> composer = TrackComposer.for_heatmap_with_tracks(n_tracks=3)
            >>> composer.add_track(heatmap, row=1, col=1)
            >>> composer.add_track(track1, row=2, col=1)
            >>> composer.add_track(track2, row=3, col=1)
            >>> composer.add_track(track3, row=4, col=1)
        """
        row_heights = [heatmap_height] + [track_height] * n_tracks
        return cls(LayoutConfig(
            rows=n_tracks + 1,
            cols=1,
            row_heights=row_heights,
            vertical_spacing=0.02,
            shared_xaxes=True,
            **kwargs
        ))

    @classmethod
    def for_cool_with_tracks(
        cls,
        n_tracks: int,
        heatmap_height: float = 10,
        track_heights: Optional[List[float]] = None,
        **kwargs
    ) -> "TrackComposer":
        """Create composer for Hi-C contact heatmap (with optional tiling) + tracks.

        Creates a vertical layout with:
        - Row 1: Hi-C contact matrix heatmap (can be tiling of two samples)
        - Rows 2 to n_tracks+1: Additional tracks (BigWig, Gene, Eigen, IS, etc.)

        This is similar to for_heatmap_with_tracks() but with different default
        spacing optimized for Hi-C contact matrix visualization.

        Args:
            n_tracks: Number of tracks below the heatmap.
            heatmap_height: Relative height of heatmap row. Default 10.
            track_heights: List of height ratios for each track.
                If None, all tracks have equal height of 2.
            **kwargs: Additional LayoutConfig parameters.

        Returns:
            TrackComposer: Configured composer instance.

        Example:
            >>> composer = TrackComposer.for_cool_with_tracks(
            ...     n_tracks=3,
            ...     heatmap_height=10,
            ...     track_heights=[2, 2, 3]  # gene track taller
            ... )
            >>> composer.add_track(heatmap, row=1, col=1)
            >>> composer.add_track(bw_track1, row=2, col=1)
            >>> composer.add_track(bw_track2, row=3, col=1)
            >>> composer.add_track(gene_track, row=4, col=1)
        """
        if track_heights is None:
            track_heights = [2.0] * n_tracks
        elif len(track_heights) != n_tracks:
            raise ValueError(
                f"track_heights length ({len(track_heights)}) must match "
                f"n_tracks ({n_tracks})"
            )

        row_heights = [heatmap_height] + track_heights

        return cls(LayoutConfig(
            rows=1 + n_tracks,
            cols=1,
            row_heights=row_heights,
            vertical_spacing=0.02,
            shared_xaxes=True,
            **kwargs
        ))

    @property
    def fig(self) -> Optional[Figure]:
        """Get the assembled figure.

        Returns:
            The assembled plotly figure, or None if build() has not been called.
        """
        return self._fig

    def get_track(self, row: int, col: int) -> Optional["Track"]:
        """Get track at specified position.

        Args:
            row: Row position.
            col: Column position.

        Returns:
            Track at position, or None if no track exists.
        """
        return self.tracks.get((row, col))

    def list_tracks(self) -> List[Tuple[int, int, str]]:
        """List all tracks with their positions.

        Returns:
            List of (row, col, track_name) tuples.
        """
        return [(row, col, track.name) for (row, col), track in self.tracks.items()]
