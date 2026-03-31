# Track architecture exports
from .hic import (
    plot_tiling_compartment,
    plot_compartment_trapezoid_2tracks,
    plot_compartment_trapezoid,
    plot_trapezoid_with_tracks,
    # Legacy functions (deprecated, kept for historical reference)
    _plot_tiling_compartment_legacy,
    _plot_compartment_trapezoid_2tracks_legacy,
)

# Track classes for multi-track composition
from .track import (
    Track,
    TrackConfig,
    EigenTrack,
    ISTrack,
    HeatmapTrack,
    BigWigTrack,
    GeneTrack,
    TrapezoidHeatmapTrack,
)

# Track composer
from .track_composer import TrackComposer, LayoutConfig