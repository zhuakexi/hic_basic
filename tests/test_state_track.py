from pathlib import Path

import pytest
from plotly.subplots import make_subplots

from hic_basic.plot import StateTrack
from hic_basic.plot.hic import plot_cool_with_tracks
from hic_basic.plot.track import resolve_state_palette
from hic_basic.plot.track_composer import LayoutConfig, TrackComposer


DATA_DIR = Path(__file__).parent / "data"
EARLY_BED = DATA_DIR / "synthetic_chromhmm_early.bed"
LATE_BED = DATA_DIR / "synthetic_chromhmm_late.bed"
REGION = "chr1:50000-550000"


def test_state_track_loads_clips_and_merges_adjacent_segments():
    track = StateTrack(EARLY_BED, palette="chromhmm_15").load_data(REGION)

    assert list(track.data.columns) == [
        "chrom", "start", "end", "state", "color", "resolved_color"
    ]
    assert track.data.iloc[0]["start"] == 50000
    assert track.data.iloc[-1]["end"] == 550000
    enhancer = track.data[track.data["state"] == "4_Strong_Enhancer"]
    assert enhancer[["start", "end"]].values.tolist() == [[100000, 240000]]
    insulator = track.data[track.data["state"] == "8_Insulator"].iloc[0]
    assert insulator["resolved_color"] == "#5f6570"


def test_state_track_renders_gaps_and_deduplicates_shared_palette_legend():
    composer = TrackComposer(LayoutConfig(
        rows=2,
        cols=1,
        row_heights=[1, 1],
        vertical_spacing=0.08,
        shared_xaxes=True,
    ))
    composer.add_track(
        StateTrack(EARLY_BED, palette="chromhmm_15", display_name="early"),
        row=1,
        col=1,
    )
    composer.add_track(
        StateTrack(LATE_BED, palette="chromhmm_15", display_name="late"),
        row=2,
        col=1,
    )
    fig = composer.build(region=REGION)

    visible_legend_names = [trace.name for trace in fig.data if trace.showlegend]
    assert len(visible_legend_names) == len(set(visible_legend_names))
    assert "1_Active_Promoter" in visible_legend_names
    assert tuple(fig.layout.xaxis2.range) == (50000, 550000)


def test_state_track_empty_region_renders_background():
    track = StateTrack(EARLY_BED).load_data("chr3:0-100000")
    fig = make_subplots(rows=1, cols=1)
    track.render(fig, row=1, col=1)

    assert track.data.empty
    assert len(fig.layout.shapes) == 1
    assert fig.layout.shapes[0].fillcolor == "#eeeeee"


def test_state_track_integrates_below_cool_heatmap():
    fig = plot_cool_with_tracks(
        coolp1=DATA_DIR / "synthetic_cool_with_tracks.cool",
        tracks=[StateTrack(EARLY_BED, palette="chromhmm_15")],
        region="chr1:0-500000",
        balance=False,
        track_heights=[1],
    )

    assert fig.data[0].type == "heatmap"
    assert any(trace.name == "1_Active_Promoter" for trace in fig.data[1:])
    assert tuple(fig.layout.xaxis.range) == tuple(fig.layout.xaxis2.range)


def test_state_track_rejects_overlapping_intervals(tmp_path):
    bed = tmp_path / "overlap.bed"
    bed.write_text("chr1\t0\t20\tA\nchr1\t10\t30\tB\n")

    with pytest.raises(ValueError, match="must not overlap"):
        StateTrack(bed).load_data("chr1:0-100")


def test_state_track_rejects_rows_missing_state(tmp_path):
    bed = tmp_path / "missing_state.bed"
    bed.write_text("chr1\t0\t20\tA\nchr1\t20\t30\n")

    with pytest.raises(ValueError, match="missing required fields"):
        StateTrack(bed).load_data("chr1:0-100")


def test_custom_palette_falls_back_to_auto_colors():
    colors = resolve_state_palette(["A", "B"], {"A": "#112233"})
    assert colors["A"] == "#112233"
    assert colors["B"].startswith("#")


def test_chromhmm_18_palette_uses_roadmap_colors():
    colors = resolve_state_palette(
        ["6_TxWk", "14_TssBiv", "18_Quies"],
        "chromhmm_18",
    )
    assert colors == {
        "6_TxWk": "#006400",
        "14_TssBiv": "#cd5c5c",
        "18_Quies": "#ffffff",
    }
