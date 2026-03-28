"""
Crossover detection from single-cell sequencing data.

This module implements methods for detecting meiotic crossover events
from single-cell whole-genome sequencing data.

Submodules:
    - detector: Lu et al. (2012) Science method for haploid cells (sperm)
    - detector_mii: HMM for sister-chromatid composition in diploid cells (mII)
    - pipeline: End-to-end pipelines for both methods

References:
    1. Lu et al. (2012). Probing Meiotic Recombination and Aneuploidy
       of Single Sperm Cells by Whole-Genome Sequencing. Science.
       DOI: 10.1126/science.1229112

    2. See issues/mII_crossover/alg.md for sister-chromatid HMM method.

Example:
    >>> # Haploid cell (sperm) analysis
    >>> from hic_basic.crossover import detect_crossovers_end2end
    >>> results = detect_crossovers_end2end(
    ...     sam_file="sample.bam",
    ...     ref_snp_file="snps.tsv",
    ...     output_dir="output/",
    ...     n_jobs=4
    ... )

    >>> # Diploid cell (mII) analysis
    >>> from hic_basic.crossover import detect_sister_chromatid_states_end2end
    >>> results = detect_sister_chromatid_states_end2end(
    ...     bam_file="mII_cell.bam",
    ...     snp_file="snps.vcf",
    ...     output_dir="output/"
    ... )
"""

from .detector import (
    HMMCrossoverDetector,
    CrossoverDetector,
    HaplotypeBlock,
    Crossover,
    CrossoverResults,
)
from .detector_mii import (
    SisterChromatidHMM,
    SisterChromatidBlock,
    SisterChromatidResults,
    estimate_emissions_from_counts,
    estimate_transitions_rough,
    rough_classification,
)
from .pipeline import (
    detect_crossovers_end2end,
    process_sample_parallel,
)
from .pipeline_mii import (
    detect_sister_chromatid_states_end2end,
    save_results_mii,
    bin_by_snp,
    convert_unit_coordinates_to_snp,
)

__all__ = [
    # Lu et al. (2012) - Haploid cells
    "HMMCrossoverDetector",
    "CrossoverDetector",
    "HaplotypeBlock",
    "Crossover",
    "CrossoverResults",
    "detect_crossovers_end2end",
    "process_sample_parallel",
    # Sister chromatid - Diploid cells (mII)
    "SisterChromatidHMM",
    "SisterChromatidBlock",
    "SisterChromatidResults",
    "detect_sister_chromatid_states_end2end",
    "save_results_mii",
    "bin_by_snp",
    "convert_unit_coordinates_to_snp",
    "estimate_emissions_from_counts",
    "estimate_transitions_rough",
    "rough_classification",
]
