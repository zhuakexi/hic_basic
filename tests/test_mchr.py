"""Tests for hic_basic.mchr._3dg_to_xr reference alignment.

Covers the haploid/diploid reference handling:
- phased (diploid) structures must use a diploid reference (e.g. mm10_dip)
- a haploid reference with a phased structure warns early and then fails
  with an actionable mismatch error instead of silently producing NaN
- unphased structures keep working with haploid references
"""
import warnings
from pathlib import Path

import pytest

from hic_basic.mchr import _3dg_to_xr
from hic_basic.hicio import get_ref_dir

REPRESENTATIVE_3DG = Path(__file__).parent / "data" / "representative_structure.3dg"


def _write_3dg(fp, rows):
    """Write a minimal .3dg file (chrom, start, x, y, z per line)."""
    with open(fp, "w") as f:
        for chrom, start, x, y, z in rows:
            f.write(f"{chrom}\t{start}\t{x}\t{y}\t{z}\n")
    return fp


def test_phased_structure_haploid_ref_warns_then_raises():
    """A phased structure against a haploid reference must warn at
    validation time and then raise a mismatch error naming the diploid
    reference instead of silently dropping all coordinates."""
    try:
        get_ref_dir()
    except (FileNotFoundError, TypeError):
        pytest.skip("No reference data available")

    with pytest.warns(UserWarning, match="diploid reference"):
        with pytest.raises(ValueError, match="diploid reference"):
            _3dg_to_xr(REPRESENTATIVE_3DG, "s", genome="mm10", binsize=20000)


def test_phased_structure_diploid_ref_keeps_homologs():
    """A phased structure aligned to a diploid reference keeps both
    homologous copies (chr1(mat)/chr1(pat)) with their own coordinates."""
    try:
        get_ref_dir()
    except (FileNotFoundError, TypeError):
        pytest.skip("No reference data available")

    out = _3dg_to_xr(REPRESENTATIVE_3DG, "s", genome="mm10_dip", binsize=20000)
    coords = out["3dg"].sel(sample_name="s")
    chroms = set(coords["chrom"].values)
    assert "chr1(mat)" in chroms
    assert "chr1(pat)" in chroms
    # both homologs carry non-null coordinates (5 bins x 3 features each)
    assert int(coords.sel(chrom=["chr1(mat)", "chr1(pat)"]).notnull().sum()) == 30


def test_unphased_structure_haploid_ref_ok(tmp_path):
    """An unphased structure against a haploid reference works without any
    ploidy warning."""
    try:
        get_ref_dir()
    except (FileNotFoundError, TypeError):
        pytest.skip("No reference data available")

    fp = _write_3dg(
        tmp_path / "unphased.3dg",
        [
            ("chr1", 0, 0.0, 0.0, 0.0),
            ("chr1", 20000, 1.0, 0.4, 0.2),
            ("chr1", 40000, 1.8, 1.1, 0.5),
            ("chr2", 0, 0.5, 0.1, 0.3),
            ("chr2", 20000, 2.0, 1.5, 0.7),
        ],
    )
    with warnings.catch_warnings(record=True) as recorded:
        warnings.simplefilter("always")
        _3dg_to_xr(fp, "s", genome="mm10", binsize=20000)
    assert not any(issubclass(w.category, UserWarning) for w in recorded)


def test_mismatched_assembly_raises(tmp_path):
    """A structure whose chromosome names do not match the reference at all
    raises a mismatch error pointing at the reference (non-ploidy branch)."""
    try:
        get_ref_dir()
    except (FileNotFoundError, TypeError):
        pytest.skip("No reference data available")

    fp = _write_3dg(
        tmp_path / "mismatch.3dg",
        [
            ("chrUn_KI270442v1", 0, 0.0, 0.0, 0.0),
            ("chrUn_KI270442v1", 20000, 1.0, 0.4, 0.2),
        ],
    )
    with pytest.raises(ValueError, match="reference genome"):
        _3dg_to_xr(fp, "s", genome="mm10", binsize=20000)


def test_no_genome_skips_alignment():
    """Without a reference genome the structure is converted as-is."""
    out = _3dg_to_xr(REPRESENTATIVE_3DG, "s")
    assert "3dg" in out.data_vars
    # without a genome, parse_3dg's original index names are kept
    assert list(out["3dg"].dims) == ["sample_name", "chr", "pos", "features"]
