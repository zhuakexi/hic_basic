# Track architecture exports
from .hic import (
    plot_tiling_compartment,
    plot_compartment_trapezoid_2tracks,
    plot_compartment_trapezoid,
    plot_trapezoid_with_tracks,
    plot_cool_with_tracks,
    # Legacy functions (deprecated, kept for historical reference)
    _plot_tiling_compartment_legacy,
    _plot_compartment_trapezoid_2tracks_legacy,
)

# Track classes for multi-track composition
from .track import (
    Track,
    TrackConfig,
    ArcTrack,
    plot_arc_track,
    StateTrack,
    resolve_state_palette,
    EigenTrack,
    ISTrack,
    HeatmapTrack,
    BigWigTrack,
    GeneTrack,
    TrapezoidHeatmapTrack,
    MarkTrack,
    plot_genome_phased_read_count_track,
)

# Track composer
from .track_composer import TrackComposer, LayoutConfig
