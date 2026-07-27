"""Optional real-data ATAC coverage test.

Run explicitly with:
HIC_BASIC_ATAC_BAM=/path/sample.bam \
HIC_BASIC_ATAC_TSS_BED=/path/tss.bed \
HIC_BASIC_ATAC_CHROM_SIZES=/path/chrom.sizes \
./run_in_env.sh pytest -m integration -q tests/test_atac_integration.py
"""

import os
from pathlib import Path

import numpy as np
import pytest

from hic_basic.atac import calculate_coverage, get_tss_regions, normalize_coverage


REQUIRED_ENV = (
    "HIC_BASIC_ATAC_BAM",
    "HIC_BASIC_ATAC_TSS_BED",
    "HIC_BASIC_ATAC_CHROM_SIZES",
)


@pytest.mark.integration
@pytest.mark.slow
@pytest.mark.skipif(
    not all(os.environ.get(name) for name in REQUIRED_ENV),
    reason="set HIC_BASIC_ATAC_BAM, HIC_BASIC_ATAC_TSS_BED, and HIC_BASIC_ATAC_CHROM_SIZES",
)
def test_real_bam_tss_coverage():
    bam_path, tss_bed, chrom_sizes = (Path(os.environ[name]) for name in REQUIRED_ENV)
    for path in (bam_path, tss_bed, chrom_sizes):
        if not path.exists():
            pytest.skip(f"integration input does not exist: {path}")

    regions = get_tss_regions(tss_bed, chrom_sizes, window_size=2000)
    coverage, distances = calculate_coverage(
        bam_path,
        regions,
        read_length=150,
        sampling=100,
    )
    normalized = normalize_coverage(coverage, read_length=150)

    assert coverage.shape == (min(100, len(regions)), 4001)
    assert len(distances) == coverage.shape[1]
    assert normalized.shape == coverage.shape
    assert np.isfinite(normalized).all()
