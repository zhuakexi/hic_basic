"""Tests for hic_basic.data module.

This test suite covers the data module's main functions:
- dupref_annote: Annotate diploid genome bins with haploid reference
- fetch_cent_chromlen: Get centromere regions and chromosome lengths
- chromosomes: Get chromosome lengths for various genomes
- fetch_centromeres: Get centromere regions
"""
import os
from pathlib import Path

import pandas as pd
import pytest

from hic_basic.genome import GenomeIdeograph
from hic_basic.data import (
    dupref_annote, 
    fetch_cent_chromlen, 
    chromosomes,
    fetch_centromeres,
    mm10,
    GRCh38,
)
from hic_basic.hicio import get_ref_dir


def test_dupref_annote():
    """Test the dupref_annote function using a reference CpG file.

    The test will be skipped if we cannot locate the expected reference file.
    """
    try:
        ref_dir = Path(get_ref_dir())
    except (FileNotFoundError, TypeError):
        pytest.skip("No reference data available for dupref_annote test")

    # Look for reference CpG files in common locations
    candidates = [
        ref_dir / "color" / "hg19.cpg.20k.txt", 
        ref_dir / "color" / "hg19.cpg.1m.txt", 
        ref_dir / "hg19.cpg.20k.txt",
        ref_dir / "hg19.cpg.1m.txt"
    ]
    ref_fp = None
    for c in candidates:
        if c.exists():
            ref_fp = c
            break
    if ref_fp is None:
        pytest.skip("No hg19 CpG reference file found in reference dir")

    genomes = GenomeIdeograph("hg19_dip")
    bins = genomes.bins(20e3, bed=True)
    bins = bins.set_index(["chrom", "start"])
    bins = bins.drop(columns=["end"]) if "end" in bins.columns else bins

    ref = pd.read_table(ref_fp, names=["chrom", "pos", "CpG"], index_col=["chrom", "pos"])

    annoted = dupref_annote(bins, ref)
    # expect that the annotated frame gained the reference columns
    assert annoted.shape[1] >= bins.shape[1] + ref.shape[1]
    # expect not all values are missing
    assert annoted.isna().sum().max() < 0.5 * annoted.shape[0]


def test_chromosomes():
    """Test chromosomes function for various genome assemblies."""
    # Test basic mm10
    chroms_mm10 = chromosomes("mm10")
    assert isinstance(chroms_mm10, pd.DataFrame)
    assert "length" in chroms_mm10.columns
    assert "chr1" in chroms_mm10.index
    assert chroms_mm10.loc["chr1", "length"] > 0
    
    # Test GRCh38
    chroms_grch38 = chromosomes("GRCh38")
    assert isinstance(chroms_grch38, pd.DataFrame)
    assert "chr1" in chroms_grch38.index
    
    # Test diploid genome
    chroms_mm10_dip = chromosomes("mm10_dip")
    assert isinstance(chroms_mm10_dip, pd.DataFrame)
    # Should have haplotype suffixes
    assert any("(mat)" in str(idx) or "(pat)" in str(idx) for idx in chroms_mm10_dip.index)
    
    # Test ordered output
    chroms_ordered = chromosomes("mm10", order=True)
    assert isinstance(chroms_ordered, pd.DataFrame)
    # Check that chr2 comes before chr10 (natural ordering)
    chr_list = chroms_ordered.index.tolist()
    if "chr2" in chr_list and "chr10" in chr_list:
        assert chr_list.index("chr2") < chr_list.index("chr10")


def test_fetch_centromeres():
    """Test fetch_centromeres for supported genomes."""
    # Test mm10
    cent_mm10 = fetch_centromeres("mm10")
    assert isinstance(cent_mm10, pd.DataFrame)
    assert all(col in cent_mm10.columns for col in ["chrom", "start", "end"])
    assert len(cent_mm10) > 0
    assert cent_mm10["start"].dtype in [int, "int64", "Int64"]
    assert cent_mm10["end"].dtype in [int, "int64", "Int64"]
    
    # Test hg19_dip
    cent_hg19_dip = fetch_centromeres("hg19_dip")
    assert isinstance(cent_hg19_dip, pd.DataFrame)
    assert all(col in cent_hg19_dip.columns for col in ["chrom", "start", "end"])
    assert len(cent_hg19_dip) > 0


def test_fetch_cent_chromlen():
    """Test fetch_cent_chromlen across a few genome names.

    This validates that we get centromere positions and chromosome lengths.
    """
    # Test mm10
    cent_chromlen = fetch_cent_chromlen("mm10")
    assert isinstance(cent_chromlen, pd.DataFrame)
    assert "chrY" in cent_chromlen.index
    assert all(col in cent_chromlen.columns for col in ["start", "end", "chrom_length"])
    # Verify that centromere start < end
    assert (cent_chromlen["start"] < cent_chromlen["end"]).all()
    # Verify that centromere end < chromosome length
    assert (cent_chromlen["end"] < cent_chromlen["chrom_length"]).all()

    # Test GRCh38
    cent_chromlen = fetch_cent_chromlen("GRCh38")
    assert isinstance(cent_chromlen, pd.DataFrame)
    assert "chrY" in cent_chromlen.index
    assert all(col in cent_chromlen.columns for col in ["start", "end", "chrom_length"])

    # Test diploid variant should include haplotypes like 'chr1(mat)'
    cent_chromlen = fetch_cent_chromlen("mm10_dip")
    assert isinstance(cent_chromlen, pd.DataFrame)
    assert any("chr1" in idx for idx in cent_chromlen.index)
    # Should have both maternal and paternal haplotypes
    assert any("(mat)" in str(idx) for idx in cent_chromlen.index)
    assert any("(pat)" in str(idx) for idx in cent_chromlen.index)
    assert all(col in cent_chromlen.columns for col in ["start", "end", "chrom_length"])


def test_mm10_ref_genome():
    """Test the mm10 RefGenome object.
    
    Validates that mm10 is properly configured with chromosomes feature.
    """
    # Test that mm10 is a RefGenome object
    from hic_basic.data import RefGenome
    assert isinstance(mm10, RefGenome)
    assert mm10.name == "mm10"
    
    # Test that chromosomes feature is accessible
    chroms = mm10.chromosomes
    assert chroms is not None
    
    # Test that chromosome data is loaded
    chrom_data = chroms.data
    assert isinstance(chrom_data, pd.DataFrame)
    assert "length" in chrom_data.columns
    assert len(chrom_data) > 0
    
    # Verify common mouse chromosomes are present
    assert "chr1" in chrom_data.index
    assert "chrX" in chrom_data.index
    assert "chrY" in chrom_data.index
    
    # Verify chromosome lengths are positive integers
    assert (chrom_data["length"] > 0).all()
    assert chrom_data["length"].dtype in [int, "int64", "Int64"]
    
    # Test chromosome alias also works (creates new instance but same data)
    chrom_alias = mm10.chromosome
    assert chrom_alias.data.equals(mm10.chromosomes.data)


def test_GRCh38_ref_genome():
    """Test the GRCh38 RefGenome object.
    
    Validates that GRCh38 is properly configured with chromosomes feature.
    """
    # Test that GRCh38 is a RefGenome object
    from hic_basic.data import RefGenome
    assert isinstance(GRCh38, RefGenome)
    assert GRCh38.name == "GRCh38"
    
    # Test that chromosomes feature is accessible
    chroms = GRCh38.chromosomes
    assert chroms is not None
    
    # Test that chromosome data is loaded
    chrom_data = chroms.data
    assert isinstance(chrom_data, pd.DataFrame)
    assert "length" in chrom_data.columns
    assert len(chrom_data) > 0
    
    # Verify common human chromosomes are present
    assert "chr1" in chrom_data.index
    assert "chrX" in chrom_data.index
    assert "chrY" in chrom_data.index
    
    # Verify chromosome lengths are positive integers
    assert (chrom_data["length"] > 0).all()
    assert chrom_data["length"].dtype in [int, "int64", "Int64"]
    
    # Test chromosome alias also works (creates new instance but same data)
    chrom_alias = GRCh38.chromosome
    assert chrom_alias.data.equals(GRCh38.chromosomes.data)
    
    # Verify hg38 alias exists
    from hic_basic.data import hg38
    assert hg38 is GRCh38


def test_ref_genome_cache_dir():
    """Test that RefGenome objects have proper cache directories."""
    # Verify cache directories exist
    assert mm10.cache_dir.exists()
    assert GRCh38.cache_dir.exists()
    
    # Verify cache directory structure
    assert "mm10" in str(mm10.cache_dir)
    assert "GRCh38" in str(GRCh38.cache_dir)
