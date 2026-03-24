"""
Crossover detection from single-cell sequencing data.

This module implements the Lu et al. (2012) Science method for detecting
meiotic crossover events from single sperm whole-genome sequencing data.

Reference:
    Lu S, Zong C, Fan W, et al. Probing Meiotic Recombination and Aneuploidy
    of Single Sperm Cells by Whole-Genome Sequencing. Science. 2012;338(6114):1627-1630.
    DOI: 10.1126/science.1229112

Key Features:
    - Two-step crossover detection (rough estimation + HMM refinement)
    - Parallel processing with intermediate caching
    - Fault-tolerant pipeline with flag files
    - Integration with hic_basic.sam and hic_basic.phasing modules

Example:
    >>> from hic_basic.crossover import detect_crossovers_end2end
    >>> results = detect_crossovers_end2end(
    ...     sam_file="sample.bam",
    ...     ref_snp_file="snps.tsv",
    ...     output_dir="output/",
    ...     n_jobs=4
    ... )
"""

from .detector import (
    HMMCrossoverDetector,
    CrossoverDetector,
    HaplotypeBlock,
    Crossover,
    CrossoverResults,
)
from .pipeline import (
    detect_crossovers_end2end,
    process_sample_parallel,
)

__all__ = [
    # Core classes
    "HMMCrossoverDetector",
    "CrossoverDetector",
    "HaplotypeBlock",
    "Crossover",
    "CrossoverResults",
    # Pipeline functions
    "detect_crossovers_end2end",
    "process_sample_parallel",
]
