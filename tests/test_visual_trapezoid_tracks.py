#!/usr/bin/env python
"""
Visual test for trapezoid Hi-C with multi-track visualization.

This script generates PNG images for visual inspection of the 
plot_trapezoid_with_tracks function.

Usage:
    python tests/test_visual_trapezoid_tracks.py

Output:
    tests/output/visual/trapezoid_tracks_*.png
"""
import os
import sys
from pathlib import Path

# Add paths
sys.path.insert(0, "/share/home/ychi/dev/hic_basic")
sys.path.insert(0, "/share/home/ychi/dev/hires_utils")

from hic_basic.plot import (
    plot_trapezoid_with_tracks,
    BigWigTrack,
    GeneTrack,
    EigenTrack,
    ISTrack
)
from hic_basic.plot.track_composer import TrackComposer

# Test data paths
DATA_DIR = Path(__file__).parent / "data"
OUTPUT_DIR = Path(__file__).parent / "output" / "visual"

COOL_FILE = DATA_DIR / "visual_test.cool"
EIGS_FILE = DATA_DIR / "visual_test_eigs.bedgraph"
BW_FILE1 = DATA_DIR / "visual_test_signal1.bw"
BW_FILE2 = DATA_DIR / "visual_test_signal2.bw"
GTF_FILE = DATA_DIR / "visual_test.gtf"

REGION = "chr1:0-5000000"  # Full chr1


def test_basic_trapezoid_with_tracks():
    """Test basic trapezoid + tracks layout."""
    print("Test 1: Basic trapezoid + 2 tracks (BigWig + Gene)")
    
    tracks = [
        BigWigTrack(str(BW_FILE1), bin_size=50000, color="#1f77b4"),
        GeneTrack(str(GTF_FILE), mode="all_transcripts", show_labels=True)
    ]
    
    fig = plot_trapezoid_with_tracks(
        coolp=str(COOL_FILE),
        vec_file=str(EIGS_FILE),
        region=REGION,
        tracks=tracks,
        track_heights=[2, 3],
        trapezoid_height=None,  # Auto-calculate
        strata=25,
        title="Trapezoid Hi-C + BigWig + Gene Track",
        balance=False,
        total_height=700,
        width=800,
        auto_scale=True  # Enable auto-scaling
    )
    
    out_path = OUTPUT_DIR / "trapezoid_tracks_basic.png"
    fig.write_image(str(out_path))
    print(f"  Saved: {out_path}")


def test_trapezoid_with_multi_bigwig():
    """Test trapezoid with multiple BigWig tracks."""
    print("Test 2: Trapezoid + 3 BigWig tracks")
    
    tracks = [
        BigWigTrack(str(BW_FILE1), bin_size=50000, color="#1f77b4"),
        BigWigTrack(str(BW_FILE2), bin_size=50000, color="#ff7f0e"),
        GeneTrack(str(GTF_FILE), mode="gene_only", show_labels=True)
    ]
    
    fig = plot_trapezoid_with_tracks(
        coolp=str(COOL_FILE),
        vec_file=str(EIGS_FILE),
        region=REGION,
        tracks=tracks,
        track_heights=[2, 2, 2],
        trapezoid_height=None,  # Auto-calculate
        strata=30,
        title="Trapezoid Hi-C + Multi BigWig + Gene",
        balance=False,
        total_height=800,
        width=900,
        auto_scale=True
    )
    
    out_path = OUTPUT_DIR / "trapezoid_tracks_multi_bw.png"
    fig.write_image(str(out_path))
    print(f"  Saved: {out_path}")


def test_trapezoid_large_strata():
    """Test trapezoid with larger strata value."""
    print("Test 3: Trapezoid with large strata (50)")
    
    tracks = [
        BigWigTrack(str(BW_FILE1), bin_size=50000, color="#1f77b4"),
    ]
    
    fig = plot_trapezoid_with_tracks(
        coolp=str(COOL_FILE),
        vec_file=str(EIGS_FILE),
        region=REGION,
        tracks=tracks,
        track_heights=[2],
        trapezoid_height=10,
        strata=50,
        title="Trapezoid Hi-C (strata=50)",
        balance=False,
        total_height=600,
        width=800
    )
    
    out_path = OUTPUT_DIR / "trapezoid_large_strata.png"
    fig.write_image(str(out_path))
    print(f"  Saved: {out_path}")


def test_trapezoid_gene_only():
    """Test trapezoid with gene track only."""
    print("Test 4: Trapezoid + Gene track (all transcripts)")
    
    tracks = [
        GeneTrack(str(GTF_FILE), mode="all_transcripts", show_labels=True)
    ]
    
    fig = plot_trapezoid_with_tracks(
        coolp=str(COOL_FILE),
        vec_file=str(EIGS_FILE),
        region=REGION,
        tracks=tracks,
        track_heights=[4],
        trapezoid_height=6,
        strata=25,
        title="Trapezoid Hi-C + Gene Annotation",
        balance=False,
        total_height=700,
        width=900
    )
    
    out_path = OUTPUT_DIR / "trapezoid_gene_only.png"
    fig.write_image(str(out_path))
    print(f"  Saved: {out_path}")


def test_trapezoid_custom_colors():
    """Test trapezoid with custom colormap."""
    print("Test 5: Trapezoid with custom colormap")
    
    tracks = [
        BigWigTrack(str(BW_FILE1), bin_size=50000, color="#2ca02c"),
        BigWigTrack(str(BW_FILE2), bin_size=50000, color="#d62728"),
    ]
    
    fig = plot_trapezoid_with_tracks(
        coolp=str(COOL_FILE),
        vec_file=str(EIGS_FILE),
        region=REGION,
        tracks=tracks,
        track_heights=[2, 2],
        trapezoid_height=6,
        strata=25,
        title="Trapezoid Hi-C (Custom Colors)",
        balance=False,
        cmap="viridis",  # Different colormap
        total_height=700,
        width=800
    )
    
    out_path = OUTPUT_DIR / "trapezoid_custom_cmap.png"
    fig.write_image(str(out_path))
    print(f"  Saved: {out_path}")


def test_trapezoid_tall():
    """Test tall figure layout."""
    print("Test 6: Tall figure layout")
    
    tracks = [
        BigWigTrack(str(BW_FILE1), bin_size=50000, color="#1f77b4"),
        BigWigTrack(str(BW_FILE2), bin_size=50000, color="#ff7f0e"),
        GeneTrack(str(GTF_FILE), mode="all_transcripts", show_labels=True),
    ]
    
    fig = plot_trapezoid_with_tracks(
        coolp=str(COOL_FILE),
        vec_file=str(EIGS_FILE),
        region=REGION,
        tracks=tracks,
        track_heights=[2, 2, 4],
        trapezoid_height=8,
        strata=30,
        title="Trapezoid Hi-C (Tall Layout)",
        balance=False,
        total_height=1000,
        width=800
    )
    
    out_path = OUTPUT_DIR / "trapezoid_tall_layout.png"
    fig.write_image(str(out_path))
    print(f"  Saved: {out_path}")


def test_trapezoid_wide():
    """Test wide figure layout."""
    print("Test 7: Wide figure layout")
    
    tracks = [
        BigWigTrack(str(BW_FILE1), bin_size=50000, color="#1f77b4"),
        GeneTrack(str(GTF_FILE), mode="gene_only", show_labels=True),
    ]
    
    fig = plot_trapezoid_with_tracks(
        coolp=str(COOL_FILE),
        vec_file=str(EIGS_FILE),
        region=REGION,
        tracks=tracks,
        track_heights=[2, 2],
        trapezoid_height=6,
        strata=25,
        title="Trapezoid Hi-C (Wide Layout)",
        balance=False,
        total_height=500,
        width=1200
    )
    
    out_path = OUTPUT_DIR / "trapezoid_wide_layout.png"
    fig.write_image(str(out_path))
    print(f"  Saved: {out_path}")


if __name__ == "__main__":
    print("=" * 60)
    print("Visual Test: Trapezoid Hi-C with Multi-Track Visualization")
    print("=" * 60)
    
    # Ensure output directory exists
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Check test data exists
    if not COOL_FILE.exists():
        print(f"Error: Test data not found: {COOL_FILE}")
        print("Run: python tests/data/create_visual_test_data.py")
        sys.exit(1)
    
    # Run all tests
    test_basic_trapezoid_with_tracks()
    test_trapezoid_with_multi_bigwig()
    test_trapezoid_large_strata()
    test_trapezoid_gene_only()
    test_trapezoid_custom_colors()
    test_trapezoid_tall()
    test_trapezoid_wide()
    
    print("\n" + "=" * 60)
    print(f"All visual tests completed!")
    print(f"Output directory: {OUTPUT_DIR}")
    print("=" * 60)
