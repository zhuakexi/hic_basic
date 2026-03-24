"""
End-to-end crossover detection pipeline with parallel processing.

This module provides the main entry point for running crossover detection
from BAM/SAM files with fault-tolerant caching and parallel execution.

Key Features:
    - Step-by-step intermediate caching
    - Flag files for completion tracking
    - Resume capability after interruption
    - Parallel processing with shared memory
    - Integration with hic_basic.sam and hic_basic.phasing

Pipeline Steps:
    1. BAM/SAM → Marked Alleles (extract alleles at SNP positions)
    2. Marked Alleles → Allele Counts (aggregate per SNP)
    3. Allele Counts → VCF (genotype calling)
    4. VCF + Phasing → Haplotype Blocks (HMM inference)
    5. Haplotype Blocks → Crossovers (transition detection)

Example:
    >>> from hic_basic.crossover import detect_crossovers_end2end
    >>> results = detect_crossovers_end2end(
    ...     sam_file="sample.bam",
    ...     ref_snp_file="snps.tsv",
    ...     phasing_file="phasing.tsv",
    ...     output_dir="output/",
    ...     n_jobs=4,
    ...     cache_dir="/tmp/cache"
    ... )
"""

import logging
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from .cache import CacheManager, SharedMemoryManager
from .detector import CrossoverDetector, CrossoverResults
from .io import (
    load_snp_data,
    load_phasing_data,
    load_phased_snp_from_hic_basic,
    save_results,
    merge_snp_with_phasing,
    extract_genotype_arrays,
    validate_input_data,
)

logger = logging.getLogger(__name__)


@dataclass
class PipelineStep:
    """
    Definition of a pipeline step.

    Attributes:
        name: Step name
        input_files: List of input file patterns
        output_file: Output file pattern
        description: Step description
    """
    name: str
    input_files: List[str]
    output_file: str
    description: str


PIPELINE_STEPS = [
    PipelineStep(
        name="mark_alleles",
        input_files=["sam_file", "ref_snp_file"],
        output_file="{sample_id}.marked_alleles.parquet",
        description="Extract alleles from BAM/SAM at SNP positions"
    ),
    PipelineStep(
        name="count_alleles",
        input_files=["{sample_id}.marked_alleles.parquet"],
        output_file="{sample_id}.allele_counts.parquet",
        description="Aggregate allele counts per SNP"
    ),
    PipelineStep(
        name="call_genotypes",
        input_files=["{sample_id}.allele_counts.parquet", "ref_snp_file"],
        output_file="{sample_id}.vcf.parquet",
        description="Call genotypes from allele counts"
    ),
    PipelineStep(
        name="infer_haplotypes",
        input_files=["{sample_id}.vcf.parquet", "phasing_file"],
        output_file="{sample_id}.haplotype_blocks.parquet",
        description="Infer haplotype blocks using HMM"
    ),
    PipelineStep(
        name="detect_crossovers",
        input_files=["{sample_id}.haplotype_blocks.parquet"],
        output_file="{sample_id}.crossovers.parquet",
        description="Detect crossover positions from haplotype transitions"
    ),
]


def _process_sample_worker(
    sample_id: str,
    sam_file: str,
    ref_snp_file: str,
    phasing_file: str,
    output_dir: str,
    cache_dir: str,
    window_size: int,
    min_fragment_size: int,
    error_rate: float,
    shared_snp_dict: Optional[str] = None
) -> Tuple[str, Dict[str, Any]]:
    """
    Worker function for processing a single sample.

    This function is designed to be called in a separate process.
    It uses cached intermediate results and flag files for fault tolerance.

    Args:
        sample_id: Sample identifier
        sam_file: Path to BAM/SAM file
        ref_snp_file: Path to reference SNP file
        phasing_file: Path to phasing information file
        output_dir: Output directory
        cache_dir: Cache directory
        window_size: HMM window size
        min_fragment_size: Minimum fragment size
        error_rate: Error rate for HMM
        shared_snp_dict: Optional key to shared SNP dictionary

    Returns:
        Tuple of (sample_id, result dictionary)
    """
    from hic_basic.phasing import sam_mark_alleles, sam_count_alleles

    try:
        cache = CacheManager(cache_dir)
        sample_params = {
            "sample_id": sample_id,
            "sam_file": sam_file,
            "ref_snp_file": ref_snp_file,
            "window_size": window_size,
            "min_fragment_size": min_fragment_size,
            "error_rate": error_rate,
        }
        sample_key = cache.compute_key(sample_params)

        logger.info(f"Processing sample {sample_id} (key: {sample_key[:8]}...)")

        result = {
            "sample_id": sample_id,
            "success": False,
            "steps_completed": [],
            "output_files": {},
            "error": None,
        }

        marked_alleles_file = cache._get_data_path(sample_key, "mark_alleles").with_suffix('.parquet')
        if not cache.exists(sample_key, "mark_alleles"):
            logger.info(f"Step 1/5: Marking alleles for {sample_id}")
            marked_allele_tmp = marked_alleles_file.with_suffix('.tsv.gz')
            sam_mark_alleles(
                sam_file=sam_file,
                phased_snp_file=ref_snp_file,
                outfile=str(marked_allele_tmp),
                show_progress=False,
                verbose=0
            )
            allele_df = pd.read_csv(
                str(marked_allele_tmp),
                names=['chrom', 'pos', 'allele', 'read_name', 'relative_pos'],
                usecols=[0, 1, 2, 3, 4]
            )
            allele_df.to_parquet(str(marked_alleles_file), index=False)
            marked_allele_tmp.unlink()
            cache.set_flag(sample_key, "mark_alleles", {"n_reads": len(allele_df)})
            logger.info(f"  Saved {len(allele_df)} marked alleles")
        else:
            logger.info(f"Step 1/5: Using cached marked alleles")
        result["steps_completed"].append("mark_alleles")

        allele_counts_file = cache._get_data_path(sample_key, "count_alleles").with_suffix('.parquet')
        if not cache.exists(sample_key, "count_alleles"):
            logger.info(f"Step 2/5: Counting alleles for {sample_id}")
            marked_alleles_df = pd.read_parquet(str(marked_alleles_file))
            allele_count_df = sam_count_alleles(
                allele_df=marked_alleles_df,
                ref_snp_file=ref_snp_file,
                gt_strategy="no_conflict"
            )
            allele_count_df.to_parquet(str(allele_counts_file), index=False)
            cache.set_flag(sample_key, "count_alleles", {"n_snps": len(allele_count_df)})
            logger.info(f"  Saved {len(allele_count_df)} allele counts")
        else:
            logger.info(f"Step 2/5: Using cached allele counts")
        result["steps_completed"].append("count_alleles")

        vcf_file = cache._get_data_path(sample_key, "call_genotypes").with_suffix('.parquet')
        if not cache.exists(sample_key, "call_genotypes"):
            logger.info(f"Step 3/5: Calling genotypes for {sample_id}")
            allele_count_df = pd.read_parquet(str(allele_counts_file))
            ref_snp_df = load_phased_snp_from_hic_basic(ref_snp_file)
            vcf_df = pd.merge(
                ref_snp_df,
                allele_count_df,
                how="left",
                on=["chromosome", "position"]
            )
            vcf_df['sperm_genotype'] = vcf_df['gt'].fillna('N')
            vcf_df.to_parquet(str(vcf_file), index=False)
            cache.set_flag(sample_key, "call_genotypes", {"n_snps": len(vcf_df)})
            logger.info(f"  Saved {len(vcf_df)} genotypes")
        else:
            logger.info(f"Step 3/5: Using cached genotypes")
        result["steps_completed"].append("call_genotypes")

        phasing_df = load_phasing_data(phasing_file)
        vcf_df = pd.read_parquet(str(vcf_file))
        merged_df = merge_snp_with_phasing(vcf_df, phasing_df)

        logger.info(f"Step 4/5: Inferring haplotypes for {sample_id}")
        detector = CrossoverDetector(
            window_size=window_size,
            min_fragment_size=min_fragment_size,
            error_rate=error_rate
        )

        results_dict = {}
        for chrom in merged_df['chromosome'].unique():
            chrom_data = merged_df[merged_df['chromosome'] == chrom].copy()
            if len(chrom_data) == 0:
                continue

            positions, sperm_genotypes, haplotype_1, haplotype_2 = extract_genotype_arrays(chrom_data)

            states, crossover_positions = detector.hmm_detector.detect(
                positions=positions,
                sperm_genotypes=sperm_genotypes,
                haplotype_1=haplotype_1,
                haplotype_2=haplotype_2
            )

            blocks = detector._states_to_blocks(chrom, positions, states)
            results_dict[chrom] = {
                "blocks": blocks,
                "crossovers": crossover_positions,
                "states": states.tolist(),
                "positions": positions.tolist(),
            }

        haplotype_blocks_list = []
        for chrom, data in results_dict.items():
            for block in data["blocks"]:
                haplotype_blocks_list.append({
                    "chromosome": chrom,
                    "start_pos": block.start_pos,
                    "end_pos": block.end_pos,
                    "haplotype": block.haplotype,
                    "confidence": block.confidence,
                    "num_snps": block.num_snps,
                })

        haplotype_blocks_df = pd.DataFrame(haplotype_blocks_list)
        haplotype_blocks_file = cache._get_data_path(sample_key, "infer_haplotypes").with_suffix('.parquet')
        haplotype_blocks_df.to_parquet(str(haplotype_blocks_file), index=False)
        cache.set_flag(sample_key, "infer_haplotypes", {"n_blocks": len(haplotype_blocks_df)})
        logger.info(f"  Saved {len(haplotype_blocks_df)} haplotype blocks")
        result["steps_completed"].append("infer_haplotypes")

        logger.info(f"Step 5/5: Detecting crossovers for {sample_id}")
        crossovers_list = []
        for chrom, data in results_dict.items():
            blocks = data["blocks"]
            crossover_positions = data["crossovers"]

            for i, pos in enumerate(crossover_positions):
                if i < len(blocks) - 1:
                    from_hap = blocks[i].haplotype
                    to_hap = blocks[i + 1].haplotype
                    confidence = (blocks[i].confidence + blocks[i + 1].confidence) / 2
                    resolution = (blocks[i + 1].start_pos - blocks[i].end_pos) / 1000
                else:
                    from_hap = 'U'
                    to_hap = 'U'
                    confidence = 0.5
                    resolution = 0

                crossovers_list.append({
                    "chromosome": chrom,
                    "position": pos,
                    "start_pos": blocks[i].end_pos if i < len(blocks) else positions[0],
                    "end_pos": blocks[i + 1].start_pos if i < len(blocks) - 1 else positions[-1],
                    "from_haplotype": from_hap,
                    "to_haplotype": to_hap,
                    "confidence": confidence,
                    "resolution_kb": resolution,
                })

        crossovers_df = pd.DataFrame(crossovers_list)
        crossovers_file = cache._get_data_path(sample_key, "detect_crossovers").with_suffix('.parquet')
        crossovers_df.to_parquet(str(crossovers_file), index=False)
        cache.set_flag(sample_key, "detect_crossovers", {"n_crossovers": len(crossovers_df)})
        logger.info(f"  Saved {len(crossovers_df)} crossovers")
        result["steps_completed"].append("detect_crossovers")

        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        blocks_out = output_path / f"{sample_id}_haplotype_blocks.tsv"
        haplotype_blocks_df.to_csv(blocks_out, sep='\t', index=False)
        result["output_files"]["haplotype_blocks"] = str(blocks_out)

        if len(crossovers_df) > 0:
            crossovers_out = output_path / f"{sample_id}_crossovers.tsv"
            crossovers_df.to_csv(crossovers_out, sep='\t', index=False)
            result["output_files"]["crossovers"] = str(crossovers_out)

        stats = {
            "total_crossovers": len(crossovers_df),
            "total_blocks": len(haplotype_blocks_df),
            "mean_confidence": haplotype_blocks_df['confidence'].mean() if len(haplotype_blocks_df) > 0 else 0,
            "mean_resolution_kb": crossovers_df['resolution_kb'].mean() if len(crossovers_df) > 0 else 0,
        }
        stats_out = output_path / f"{sample_id}_statistics.tsv"
        pd.DataFrame([stats]).to_csv(stats_out, sep='\t', index=False)
        result["output_files"]["statistics"] = str(stats_out)
        result["statistics"] = stats

        result["success"] = True
        logger.info(f"Completed processing sample {sample_id}")

        return sample_id, result

    except Exception as e:
        logger.error(f"Error processing sample {sample_id}: {e}")
        return sample_id, {
            "sample_id": sample_id,
            "success": False,
            "error": str(e),
            "steps_completed": [],
            "output_files": {},
        }


def process_sample_parallel(
    sample_ids: List[str],
    sam_files: Dict[str, str],
    ref_snp_file: str,
    phasing_file: str,
    output_dir: str,
    n_jobs: int = 4,
    cache_dir: Optional[str] = None,
    window_size: int = 10,
    min_fragment_size: int = 300000,
    error_rate: float = 0.05,
) -> Dict[str, Dict[str, Any]]:
    """
    Process multiple samples in parallel.

    Args:
        sample_ids: List of sample IDs to process
        sam_files: Dictionary mapping sample_id to BAM/SAM file path
        ref_snp_file: Path to reference SNP file
        phasing_file: Path to phasing information file
        output_dir: Output directory
        n_jobs: Number of parallel processes
        cache_dir: Cache directory (default: output_dir/.cache)
        window_size: HMM window size
        min_fragment_size: Minimum fragment size
        error_rate: Error rate for HMM

    Returns:
        Dictionary mapping sample_id to result dictionary
    """
    if cache_dir is None:
        cache_dir = str(Path(output_dir) / ".cache")

    Path(cache_dir).mkdir(parents=True, exist_ok=True)
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    logger.info(f"Processing {len(sample_ids)} samples with {n_jobs} processes")

    all_results = {}

    with ProcessPoolExecutor(max_workers=n_jobs) as executor:
        futures = {}

        for sample_id in sample_ids:
            if sample_id not in sam_files:
                logger.warning(f"Sample {sample_id} not found in sam_files, skipping")
                continue

            sam_file = sam_files[sample_id]

            future = executor.submit(
                _process_sample_worker,
                sample_id=sample_id,
                sam_file=sam_file,
                ref_snp_file=ref_snp_file,
                phasing_file=phasing_file,
                output_dir=output_dir,
                cache_dir=cache_dir,
                window_size=window_size,
                min_fragment_size=min_fragment_size,
                error_rate=error_rate,
            )
            futures[future] = sample_id

        for future in as_completed(futures):
            sample_id, result = future.result()
            all_results[sample_id] = result

            if result["success"]:
                logger.info(f"✓ {sample_id}: {result['statistics'].get('total_crossovers', 0)} crossovers detected")
            else:
                logger.error(f"✗ {sample_id}: {result.get('error', 'Unknown error')}")

    return all_results


def detect_crossovers_end2end(
    sam_file: str,
    ref_snp_file: str,
    phasing_file: str,
    output_dir: str,
    sample_id: Optional[str] = None,
    n_jobs: int = 1,
    cache_dir: Optional[str] = None,
    window_size: int = 10,
    min_fragment_size: int = 300000,
    error_rate: float = 0.05,
    chromosomes: Optional[List[str]] = None,
) -> Union[CrossoverResults, Dict[str, CrossoverResults]]:
    """
    End-to-end crossover detection from BAM/SAM file.

    This is the main entry point for the crossover detection pipeline.
    It processes BAM/SAM files through all steps:
    1. Extract alleles at SNP positions
    2. Aggregate allele counts
    3. Call genotypes
    4. Infer haplotype blocks using HMM
    5. Detect crossover positions

    The pipeline uses intermediate caching with flag files for fault tolerance.
    If interrupted, it can resume from the last completed step.

    Args:
        sam_file: Path to input BAM/SAM file
        ref_snp_file: Path to reference SNP file (TSV or parquet)
                      Expected columns: chromosome, position, haplotype_1, haplotype_2
        phasing_file: Path to phasing information file
                      Expected columns: chromosome, position, haplotype_1, haplotype_2
        output_dir: Output directory for results
        sample_id: Sample ID (default: derived from sam_file name)
        n_jobs: Number of parallel processes (default: 1)
        cache_dir: Cache directory (default: output_dir/.cache)
        window_size: SNP window size for rough estimation
        min_fragment_size: Minimum fragment size in bp
        error_rate: SNP calling error rate for HMM
        chromosomes: Optional list of chromosomes to analyze (default: all)

    Returns:
        CrossoverResults object (single sample) or dictionary of results (multiple samples)

    Raises:
        FileNotFoundError: If input files don't exist
        ValueError: If input data is invalid

    Example:
        Single sample:
        >>> results = detect_crossovers_end2end(
        ...     sam_file="sperm01.bam",
        ...     ref_snp_file="snps.tsv",
        ...     phasing_file="phasing.tsv",
        ...     output_dir="output/",
        ...     n_jobs=4
        ... )
        >>> print(f"Detected {len(results.crossovers)} crossovers")

        Multiple samples:
        >>> sam_files = {
        ...     "S01": "sperm01.bam",
        ...     "S02": "sperm02.bam",
        ... }
        >>> all_results = detect_crossovers_end2end(
        ...     sam_file=None,
        ...     sam_files=sam_files,
        ...     ref_snp_file="snps.tsv",
        ...     phasing_file="phasing.tsv",
        ...     output_dir="output/",
        ...     n_jobs=4
        ... )
    """
    if sample_id is None:
        sample_id = Path(sam_file).stem

    if not Path(sam_file).exists():
        raise FileNotFoundError(f"SAM/BAM file not found: {sam_file}")
    if not Path(ref_snp_file).exists():
        raise FileNotFoundError(f"Reference SNP file not found: {ref_snp_file}")
    if not Path(phasing_file).exists():
        raise FileNotFoundError(f"Phasing file not found: {phasing_file}")

    if cache_dir is None:
        cache_dir = str(Path(output_dir) / ".cache")

    logger.info(f"Starting crossover detection pipeline")
    logger.info(f"  Sample: {sample_id}")
    logger.info(f"  SAM file: {sam_file}")
    logger.info(f"  Reference SNPs: {ref_snp_file}")
    logger.info(f"  Phasing: {phasing_file}")
    logger.info(f"  Output: {output_dir}")
    logger.info(f"  Parallel jobs: {n_jobs}")

    if n_jobs > 1:
        logger.info("Using parallel processing mode")
        sam_files = {sample_id: sam_file}
        results = process_sample_parallel(
            sample_ids=[sample_id],
            sam_files=sam_files,
            ref_snp_file=ref_snp_file,
            phasing_file=phasing_file,
            output_dir=output_dir,
            n_jobs=n_jobs,
            cache_dir=cache_dir,
            window_size=window_size,
            min_fragment_size=min_fragment_size,
            error_rate=error_rate,
        )

        sample_result = results.get(sample_id, {})
        if not sample_result.get("success", False):
            raise RuntimeError(f"Pipeline failed: {sample_result.get('error', 'Unknown error')}")

        detector = CrossoverDetector(
            window_size=window_size,
            min_fragment_size=min_fragment_size,
            error_rate=error_rate
        )

        crossover_results = CrossoverResults(sperm_id=sample_id)

        blocks_file = Path(output_dir) / f"{sample_id}_haplotype_blocks.tsv"
        if blocks_file.exists():
            blocks_df = pd.read_csv(blocks_file, sep='\t')
            for _, row in blocks_df.iterrows():
                crossover_results.haplotype_blocks.append(HaplotypeBlock(
                    chromosome=row['chromosome'],
                    start_pos=row['start_pos'],
                    end_pos=row['end_pos'],
                    haplotype=row['haplotype'],
                    confidence=row['confidence'],
                    num_snps=row['num_snps'],
                ))

        crossovers_file = Path(output_dir) / f"{sample_id}_crossovers.tsv"
        if crossovers_file.exists():
            crossovers_df = pd.read_csv(crossovers_file, sep='\t')
            for _, row in crossovers_df.iterrows():
                crossover_results.crossovers.append(Crossover(
                    sperm_id=sample_id,
                    chromosome=row['chromosome'],
                    position=row['position'],
                    start_pos=row['start_pos'],
                    end_pos=row['end_pos'],
                    from_haplotype=row['from_haplotype'],
                    to_haplotype=row['to_haplotype'],
                    confidence=row['confidence'],
                    resolution_kb=row['resolution_kb'],
                ))

        stats_file = Path(output_dir) / f"{sample_id}_statistics.tsv"
        if stats_file.exists():
            stats_df = pd.read_csv(stats_file, sep='\t')
            crossover_results.statistics = stats_df.iloc[0].to_dict()

        return crossover_results

    else:
        logger.info("Using sequential processing mode")

        from hic_basic.phasing import sam_mark_alleles, sam_count_alleles

        cache = CacheManager(cache_dir)
        params = {
            "sample_id": sample_id,
            "sam_file": sam_file,
            "ref_snp_file": ref_snp_file,
            "window_size": window_size,
            "min_fragment_size": min_fragment_size,
            "error_rate": error_rate,
        }
        key = cache.compute_key(params)

        Path(output_dir).mkdir(parents=True, exist_ok=True)

        marked_alleles_file = cache._get_data_path(key, "mark_alleles").with_suffix('.parquet')
        if not cache.exists(key, "mark_alleles"):
            logger.info("Step 1/5: Marking alleles...")
            marked_allele_tmp = Path(output_dir) / f"{sample_id}.marked_alleles.tmp.tsv.gz"
            sam_mark_alleles(
                sam_file=sam_file,
                phased_snp_file=ref_snp_file,
                outfile=str(marked_allele_tmp),
                show_progress=False,
                verbose=0
            )
            allele_df = pd.read_csv(
                str(marked_allele_tmp),
                names=['chrom', 'pos', 'allele', 'read_name', 'relative_pos'],
                usecols=[0, 1, 2, 3, 4]
            )
            allele_df.to_parquet(str(marked_alleles_file), index=False)
            marked_allele_tmp.unlink()
            cache.set_flag(key, "mark_alleles", {"n_reads": len(allele_df)})
        else:
            logger.info("Step 1/5: Using cached marked alleles")

        allele_counts_file = cache._get_data_path(key, "count_alleles").with_suffix('.parquet')
        if not cache.exists(key, "count_alleles"):
            logger.info("Step 2/5: Counting alleles...")
            marked_alleles_df = pd.read_parquet(str(marked_alleles_file))
            allele_count_df = sam_count_alleles(
                allele_df=marked_alleles_df,
                ref_snp_file=ref_snp_file,
                gt_strategy="no_conflict"
            )
            allele_count_df.to_parquet(str(allele_counts_file), index=False)
            cache.set_flag(key, "count_alleles", {"n_snps": len(allele_count_df)})
        else:
            logger.info("Step 2/5: Using cached allele counts")

        vcf_file = cache._get_data_path(key, "call_genotypes").with_suffix('.parquet')
        if not cache.exists(key, "call_genotypes"):
            logger.info("Step 3/5: Calling genotypes...")
            allele_count_df = pd.read_parquet(str(allele_counts_file))
            ref_snp_df = load_phased_snp_from_hic_basic(ref_snp_file)
            vcf_df = pd.merge(
                ref_snp_df,
                allele_count_df,
                how="left",
                on=["chromosome", "position"]
            )
            vcf_df['sperm_genotype'] = vcf_df['gt'].fillna('N')
            vcf_df.to_parquet(str(vcf_file), index=False)
            cache.set_flag(key, "call_genotypes", {"n_snps": len(vcf_df)})
        else:
            logger.info("Step 3/5: Using cached genotypes")

        logger.info("Step 4/5: Running HMM inference...")
        detector = CrossoverDetector(
            window_size=window_size,
            min_fragment_size=min_fragment_size,
            error_rate=error_rate
        )

        vcf_df = pd.read_parquet(str(vcf_file))
        phasing_df = load_phasing_data(phasing_file)

        if chromosomes is not None:
            vcf_df = vcf_df[vcf_df['chromosome'].isin(chromosomes)]
            phasing_df = phasing_df[phasing_df['chromosome'].isin(chromosomes)]

        merged_df = merge_snp_with_phasing(vcf_df, phasing_df)
        detector.snp_data = merged_df.rename(columns={'sperm_genotype': 'sperm_genotype'})
        detector.phasing_data = phasing_df

        results = detector.detect_crossovers(sample_id)

        logger.info("Step 5/5: Saving results...")
        output_files = save_results(results, output_dir)

        cache.set_flag(key, "complete", {
            "n_crossovers": len(results.crossovers),
            "n_blocks": len(results.haplotype_blocks),
            "output_files": output_files,
        })

        logger.info(f"Pipeline completed: {len(results.crossovers)} crossovers detected")

        return results


def get_pipeline_status(
    sam_file: str,
    ref_snp_file: str,
    cache_dir: str,
    sample_id: Optional[str] = None,
    window_size: int = 10,
    min_fragment_size: int = 300000,
    error_rate: float = 0.05,
) -> Dict[str, bool]:
    """
    Get the completion status of pipeline steps.

    Args:
        sam_file: Path to BAM/SAM file
        ref_snp_file: Path to reference SNP file
        cache_dir: Cache directory
        sample_id: Sample ID
        window_size: HMM window size
        min_fragment_size: Minimum fragment size
        error_rate: Error rate for HMM

    Returns:
        Dictionary mapping step name to completion status
    """
    if sample_id is None:
        sample_id = Path(sam_file).stem

    cache = CacheManager(cache_dir)
    params = {
        "sample_id": sample_id,
        "sam_file": sam_file,
        "ref_snp_file": ref_snp_file,
        "window_size": window_size,
        "min_fragment_size": min_fragment_size,
        "error_rate": error_rate,
    }
    key = cache.compute_key(params)

    status = {}
    for step in PIPELINE_STEPS:
        status[step.name] = cache.is_flagged(key, step.name)

    return status


def resume_pipeline(
    sam_file: str,
    ref_snp_file: str,
    phasing_file: str,
    output_dir: str,
    sample_id: Optional[str] = None,
    cache_dir: Optional[str] = None,
    **kwargs
) -> CrossoverResults:
    """
    Resume an interrupted pipeline from the last completed step.

    Args:
        sam_file: Path to BAM/SAM file
        ref_snp_file: Path to reference SNP file
        phasing_file: Path to phasing file
        output_dir: Output directory
        sample_id: Sample ID
        cache_dir: Cache directory
        **kwargs: Additional arguments passed to detect_crossovers_end2end

    Returns:
        CrossoverResults object
    """
    if cache_dir is None:
        cache_dir = str(Path(output_dir) / ".cache")

    status = get_pipeline_status(
        sam_file=sam_file,
        ref_snp_file=ref_snp_file,
        cache_dir=cache_dir,
        sample_id=sample_id,
    )

    logger.info("Pipeline status:")
    for step, completed in status.items():
        symbol = "✓" if completed else "○"
        logger.info(f"  {symbol} {step}")

    if all(status.values()):
        logger.info("All steps completed, loading results...")

        if sample_id is None:
            sample_id = Path(sam_file).stem

        results = CrossoverResults(sperm_id=sample_id)

        blocks_file = Path(output_dir) / f"{sample_id}_haplotype_blocks.tsv"
        if blocks_file.exists():
            blocks_df = pd.read_csv(blocks_file, sep='\t')
            for _, row in blocks_df.iterrows():
                results.haplotype_blocks.append(HaplotypeBlock(
                    chromosome=row['chromosome'],
                    start_pos=row['start_pos'],
                    end_pos=row['end_pos'],
                    haplotype=row['haplotype'],
                    confidence=row['confidence'],
                    num_snps=row['num_snps'],
                ))

        crossovers_file = Path(output_dir) / f"{sample_id}_crossovers.tsv"
        if crossovers_file.exists():
            crossovers_df = pd.read_csv(crossovers_file, sep='\t')
            for _, row in crossovers_df.iterrows():
                results.crossovers.append(Crossover(
                    sperm_id=sample_id,
                    chromosome=row['chromosome'],
                    position=row['position'],
                    start_pos=row['start_pos'],
                    end_pos=row['end_pos'],
                    from_haplotype=row['from_haplotype'],
                    to_haplotype=row['to_haplotype'],
                    confidence=row['confidence'],
                    resolution_kb=row['resolution_kb'],
                ))

        stats_file = Path(output_dir) / f"{sample_id}_statistics.tsv"
        if stats_file.exists():
            stats_df = pd.read_csv(stats_file, sep='\t')
            results.statistics = stats_df.iloc[0].to_dict()

        return results

    logger.info("Resuming pipeline from last completed step...")
    return detect_crossovers_end2end(
        sam_file=sam_file,
        ref_snp_file=ref_snp_file,
        phasing_file=phasing_file,
        output_dir=output_dir,
        sample_id=sample_id,
        cache_dir=cache_dir,
        **kwargs
    )
