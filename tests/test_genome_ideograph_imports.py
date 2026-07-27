import importlib

import pytest

from hic_basic.genome import GenomeIdeograph


def test_modules_import_genome_ideograph_from_new_location():
    di = importlib.import_module("hic_basic.DI")
    simpute = importlib.import_module("hic_basic.impute.simpute")

    assert di.GenomeIdeograph is GenomeIdeograph
    assert simpute.GenomeIdeograph is GenomeIdeograph
    assert "chr1" in di.GenomeIdeograph("mm10").chromosomes.index
    assert "chr1" in simpute.GenomeIdeograph("mm10").chromosomes.index


def test_legacy_binnify_import_warns():
    with pytest.warns(DeprecationWarning, match="hic_basic.genome"):
        from hic_basic.binnify import GenomeIdeograph as LegacyGenomeIdeograph

    assert LegacyGenomeIdeograph is GenomeIdeograph
