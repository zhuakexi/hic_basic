import random

import numpy as np
import pysam
import pytest

from hic_basic.atac import (
    calculate_coverage,
    calculate_fragment_lengths,
    get_tss_regions,
    normalize_coverage,
    plot_fragment_length_distribution,
)


@pytest.fixture
def atac_inputs(tmp_path):
    bam_path = tmp_path / "reads.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 100}],
    }

    with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
        for name, start, length in (
            ("boundary_read", 0, 1),
            ("read1", 48, 5),
            ("read2", 50, 3),
        ):
            read = pysam.AlignedSegment()
            read.query_name = name
            read.query_sequence = "A" * length
            read.flag = 0
            read.reference_id = 0
            read.reference_start = start
            read.mapping_quality = 60
            read.cigar = ((0, length),)
            read.query_qualities = pysam.qualitystring_to_array("I" * length)
            bam.write(read)
    pysam.index(str(bam_path))

    chrom_sizes = tmp_path / "chrom.sizes"
    chrom_sizes.write_text("chr1\t100\n")
    tss_bed = tmp_path / "tss.bed"
    tss_bed.write_text("chr1\t50\t51\nchr1\t70\t71\n")
    return bam_path, tss_bed, chrom_sizes


def test_calculate_coverage_uses_indexed_bam_and_fills_zeros(atac_inputs):
    bam_path, tss_bed, chrom_sizes = atac_inputs
    regions = get_tss_regions(tss_bed, chrom_sizes, window_size=2)

    coverage, distances = calculate_coverage(bam_path, regions, read_length=5)

    np.testing.assert_array_equal(
        coverage,
        np.array([[1, 1, 2, 2, 2], [0, 0, 0, 0, 0]]),
    )
    assert distances == [-2, -1, 0, 1, 2]


def test_calculate_coverage_sampling_is_reproducible(atac_inputs):
    bam_path, tss_bed, chrom_sizes = atac_inputs
    regions = get_tss_regions(tss_bed, chrom_sizes, window_size=2)
    random.seed(7)

    coverage, distances = calculate_coverage(
        bam_path,
        regions,
        read_length=5,
        sampling=1,
    )

    assert coverage.shape == (1, 5)
    np.testing.assert_array_equal(coverage[0], [0, 0, 0, 0, 0])
    assert distances == [-2, -1, 0, 1, 2]


def test_calculate_coverage_pads_chromosome_boundary(atac_inputs, tmp_path):
    bam_path, _, chrom_sizes = atac_inputs
    tss_bed = tmp_path / "boundary_tss.bed"
    tss_bed.write_text("chr1\t0\t1\nchr1\t50\t51\n")
    regions = get_tss_regions(tss_bed, chrom_sizes, window_size=2)

    coverage, distances = calculate_coverage(bam_path, regions, read_length=5)

    np.testing.assert_array_equal(coverage[0], [0, 0, 1, 0, 0])
    assert coverage.shape == (2, 5)
    assert distances == [-2, -1, 0, 1, 2]


def test_normalize_coverage_greenleaf():
    coverage = np.array([[1, 1, 2, 2, 2], [0, 0, 0, 0, 0]])

    normalized = normalize_coverage(coverage, read_length=5)

    np.testing.assert_allclose(
        normalized,
        np.array([[4 / 3, 4 / 3, 8 / 3, 8 / 3, 8 / 3], [0, 0, 0, 0, 0]]),
    )


def test_normalize_coverage_per_million_reads(atac_inputs):
    bam_path, _, _ = atac_inputs
    coverage = np.array([[1, 2]])

    normalized = normalize_coverage(
        coverage,
        read_length=5,
        greenleaf_norm=False,
        bam_path=bam_path,
    )

    np.testing.assert_allclose(normalized, [[1_000_000 / 3, 2_000_000 / 3]])


def test_calculate_fragment_lengths_and_plot(atac_inputs):
    bam_path, _, _ = atac_inputs

    fragment_lengths = calculate_fragment_lengths(bam_path)
    figure = plot_fragment_length_distribution(np.array([50, 100, 150]))

    np.testing.assert_array_equal(fragment_lengths, [0, 0, 0])
    assert len(figure.data) == 1
