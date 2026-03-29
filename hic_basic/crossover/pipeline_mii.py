"""
End-to-end pipeline for sister chromatid composition inference.

This module provides the main entry point for running sister chromatid
state inference from BAM/BCF files with intermediate caching.

Pipeline Steps:
    1. BAM/BCF → Read Counts (extract ref/alt/other counts at SNP positions)
    2. Read Counts → Binning (optional: merge adjacent SNPs into units)
    3. Binned Data → Rough Classification (initial state labels)
    4. Rough Classification → Parameter Estimation (emission & transition)
    5. Parameters + Binned Data → HMM Inference (Viterbi decoding)
    6. HMM States → Blocks (contiguous segment grouping)
    7. Unit-level coordinates → Original SNP coordinates (conversion)
"""

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import numpy as np
import pandas as pd

from .cache import CacheManager
from .detector_mii import (
    SisterChromatidBlock,
    SisterChromatidHMM,
    SisterChromatidResults,
    estimate_emissions_from_counts,
    estimate_transitions_rough,
    rough_classification,
)

logger = logging.getLogger(__name__)


def bin_by_snp(read_counts_df: pd.DataFrame, n_snps: int = 10) -> pd.DataFrame:
    """
    Merge adjacent SNPs into units based on SNP count.
    
    Args:
        read_counts_df: DataFrame with columns [chromosome or chrom, pos or position, ref, alt, other]
        n_snps: Number of SNPs per unit (default: 10)
    
    Returns:
        DataFrame with aggregated units, columns: [chromosome, start_pos, end_pos, ref, alt, other, n_snps, position]
    """
    if len(read_counts_df) == 0:
        return read_counts_df
    
    # Normalize column names
    df = read_counts_df.copy()
    if 'chrom' in df.columns and 'chromosome' not in df.columns:
        df['chromosome'] = df['chrom']
    if 'pos' in df.columns and 'position' not in df.columns:
        df['position'] = df['pos']
    
    # Create unit IDs
    unit_ids = np.arange(len(df)) // n_snps
    df['unit_id'] = unit_ids
    
    # Aggregate
    aggregated = df.groupby('unit_id').agg({
        'chromosome': 'first',
        'position': ['first', 'last'],
        'ref': 'sum',
        'alt': 'sum',
        'other': 'sum',
    }).reset_index()
    
    # Flatten column names
    aggregated.columns = ['unit_id', 'chromosome', 'start_pos', 'end_pos', 'ref', 'alt', 'other']
    aggregated['n_snps'] = n_snps
    
    # Use start_pos as position for HMM
    aggregated['position'] = aggregated['start_pos']
    
    return aggregated


def convert_unit_coordinates_to_snp(unit_results: SisterChromatidResults, 
                                     original_read_counts: pd.DataFrame,
                                     n_snps: int = 10) -> SisterChromatidResults:
    """
    Convert unit-level state assignments back to original SNP coordinates.
    
    Args:
        unit_results: Results from HMM at unit level
        original_read_counts: Original SNP-level read counts
        n_snps: Number of SNPs per unit
    
    Returns:
        SisterChromatidResults with SNP-level states
    """
    # Normalize column names
    if 'pos' in original_read_counts.columns and 'position' not in original_read_counts.columns:
        original_read_counts = original_read_counts.copy()
        original_read_counts['position'] = original_read_counts['pos']
    
    # Expand states from units to SNPs
    snp_states = []
    snp_positions = []
    
    for i, state in enumerate(unit_results.states):
        # Get SNPs in this unit
        start_idx = i * n_snps
        end_idx = min((i + 1) * n_snps, len(original_read_counts))
        
        for j in range(start_idx, end_idx):
            if j < len(original_read_counts):
                snp_states.append(state)
                snp_positions.append(original_read_counts.iloc[j]['position'])
    
    # Recalculate blocks at SNP level
    snp_blocks = []
    n_snps_total = len(snp_states)
    if n_snps_total > 0:
        i = 0
        while i < n_snps_total:
            current_state = snp_states[i]
            start_pos = snp_positions[i]
            start_idx = i
            
            while i < n_snps_total and snp_states[i] == current_state:
                i += 1
            
            end_pos = snp_positions[i - 1]
            num_snps = i - start_idx
            
            block = SisterChromatidBlock(
                chromosome=unit_results.chromosome,
                start_pos=start_pos,
                end_pos=end_pos,
                state=current_state,
                confidence=1.0,
                num_snps=num_snps
            )
            snp_blocks.append(block)
    
    # Create new results object
    snp_results = SisterChromatidResults(
        sample_id=unit_results.sample_id,
        chromosome=unit_results.chromosome,
        states=snp_states,
        positions=snp_positions,
        blocks=snp_blocks,
        statistics=unit_results.statistics
    )
    
    return snp_results


def extract_read_counts(
    bam_file: str,
    snp_file: str,
    output_file: Optional[str] = None,
    min_mapq: int = 20,
    min_baseq: int = 20
) -> pd.DataFrame:
    """
    Extract read counts at SNP positions from BAM file.

    Args:
        bam_file: Path to BAM/BCF file
        snp_file: Path to SNP file (VCF or TSV with positions)
        output_file: Optional output file path
        min_mapq: Minimum mapping quality
        min_baseq: Minimum base quality

    Returns:
        DataFrame with columns: chromosome, position, ref, alt, other
    """
    import pysam

    logger.info(f"Extracting read counts from {bam_file}")

    # Load SNP positions
    if snp_file.endswith('.vcf') or snp_file.endswith('.vcf.gz'):
        snps = []
        with pysam.VariantFile(snp_file) as vcf:
            for record in vcf:
                snps.append({
                    'chromosome': record.chrom,
                    'position': record.pos,
                    'ref': record.ref,
                    'alt': record.alts[0] if record.alts else '.',
                })
        snp_df = pd.DataFrame(snps)
    else:
        snp_df = pd.read_csv(snp_file, sep='\t')

    # Count reads at each SNP
    counts = []
    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        for _, row in snp_df.iterrows():
            chrom = row['chromosome']
            pos = row['position']

            ref_count = 0
            alt_count = 0
            other_count = 0

            try:
                for pileupcolumn in bam.pileup(
                    chrom,
                    start=pos - 1,
                    end=pos,
                    stepper='samtools',
                    min_mapping_quality=min_mapq
                ):
                    if pileupcolumn.pos != pos - 1:
                        continue

                    for pileupread in pileupcolumn.pileups:
                        if pileupread.is_del or pileupread.is_refskip:
                            continue

                        base = pileupread.alignment.query_sequence[pileupread.query_position]
                        qual = pileupread.alignment.query_qualities[pileupread.query_position]

                        if qual < min_baseq:
                            continue

                        ref_allele = row.get('ref', '')
                        alt_allele = row.get('alt', '')

                        if len(ref_allele) == 1 and len(alt_allele) == 1:
                            if base.upper() == ref_allele.upper():
                                ref_count += 1
                            elif base.upper() == alt_allele.upper():
                                alt_count += 1
                            else:
                                other_count += 1
            except ValueError:
                # Chromosome not found in BAM
                pass

            counts.append({
                'chromosome': chrom,
                'position': pos,
                'ref': ref_count,
                'alt': alt_count,
                'other': other_count,
            })

    counts_df = pd.DataFrame(counts)

    if output_file:
        counts_df.to_csv(output_file, sep='\t', index=False)
        logger.info(f"Saved read counts to {output_file}")

    return counts_df


def _process_sample_mii(
    sample_id: str,
    bam_file: str,
    snp_file: str,
    output_dir: str,
    cache_dir: str,
    chromosome: Optional[str] = None,
    transition_params: str = "default",  # "default" | "estimated"
    alpha: float = 0.01,
    beta: float = 0.02,
    gamma: float = 0.001,
    rho: float = 0.1,
    depth_threshold: int = 2,
    other_threshold: float = 0.4,
    balance_threshold: float = 0.25,
    binning_method: Optional[str] = None,
    binning_snp_number: int = 250,
) -> Dict[str, Any]:
    """
    Process a single sample for sister chromatid inference.

    Args:
        sample_id: Sample identifier
        bam_file: Path to BAM file
        snp_file: Path to SNP file
        output_dir: Output directory
        cache_dir: Cache directory
        chromosome: Optional chromosome to analyze
        transition_params: Transition parameter source ("default" | "estimated")
        alpha, beta, gamma, rho: Transition parameters (used when transition_params="default")
        depth_threshold: Minimum depth for confident call (default: 2)
        other_threshold: Maximum other fraction for confident call (default: 0.4)
        balance_threshold: Threshold for classifying as balanced (default: 0.25)
                          PP: ref_frac > 0.75, MM: alt_frac > 0.75
        binning_method: Binning method (None for single-SNP, "by_snp" for merging by SNP count)
        binning_snp_number: Number of SNPs per unit (default: 250)

    Returns:
        Dictionary with results and statistics
    """
    cache = CacheManager(cache_dir)

    # Compute cache key
    params = {
        "sample_id": sample_id,
        "bam_file": bam_file,
        "snp_file": snp_file,
        "chromosome": chromosome,
        "alpha": alpha,
        "beta": beta,
        "gamma": gamma,
        "rho": rho,
    }
    key = cache.compute_key(params)

    result = {
        "sample_id": sample_id,
        "success": False,
        "steps_completed": [],
        "output_files": {},
        "error": None,
        "results": None,
    }

    try:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # Step 1: Extract read counts
        counts_file = output_path / f"{sample_id}_read_counts.tsv"
        if chromosome:
            counts_file = output_path / f"{sample_id}_{chromosome}_read_counts.tsv"

        if not cache.exists(key, "extract_counts"):
            logger.info(f"Step 1/4: Extracting read counts for {sample_id}")
            counts_df = extract_read_counts(bam_file, snp_file, str(counts_file))
            cache.save(counts_df, key, "extract_counts")
            cache.set_flag(key, "extract_counts", {"n_snps": len(counts_df)})
        else:
            logger.info("Step 1/4: Using cached read counts")
            counts_df = cache.load(key, "extract_counts")

        if chromosome:
            counts_df = counts_df[counts_df['chromosome'] == chromosome]

        if len(counts_df) == 0:
            raise ValueError("No SNPs found after filtering")

        result["steps_completed"].append("extract_counts")

        # Step 1.5: Binning (optional)
        original_counts_df = counts_df.copy()  # Keep original for coordinate conversion
        if binning_method == "by_snp":
            logger.info(f"Step 1.5/5: Binning SNPs into units (n_snps={binning_snp_number})")
            counts_df = bin_by_snp(counts_df, n_snps=binning_snp_number)
            logger.info(f"  Binned {len(original_counts_df)} SNPs into {len(counts_df)} units")
            result["steps_completed"].append("binning")
        elif binning_method is None:
            logger.info("Step 1.5/5: Using single-SNP mode (no binning)")
        else:
            logger.warning(f"Unknown binning_method: {binning_method}, using single-SNP mode")

        # Step 2: Rough classification
        if not cache.exists(key, "rough_classification"):
            logger.info("Step 2/4: Performing rough classification")
            rough_labels = rough_classification(
                counts_df,
                depth_threshold=depth_threshold,
                other_threshold=other_threshold,
                balance_threshold=balance_threshold
            )
            cache.save(rough_labels, key, "rough_classification")
            cache.set_flag(key, "rough_classification", {"n_snps": len(rough_labels)})
        else:
            logger.info("Step 2/4: Using cached rough classification")
            rough_labels = cache.load(key, "rough_classification")

        result["steps_completed"].append("rough_classification")

        # Step 3: Estimate parameters
        if not cache.exists(key, "estimate_params"):
            logger.info("Step 3/4: Estimating HMM parameters")
            emission_params = estimate_emissions_from_counts(counts_df, rough_labels)
            
            # Use default or estimated transition parameters
            if transition_params == "estimated":
                logger.info("  Using estimated transition parameters from rough classification")
                alpha_est, beta_est, gamma_est, rho_est = estimate_transitions_rough(rough_labels)
            else:
                logger.info(f"  Using default transition parameters: alpha={alpha}, beta={beta}, gamma={gamma}, rho={rho}")
                alpha_est, beta_est, gamma_est, rho_est = alpha, beta, gamma, rho

            params_dict = {
                "emission_params": emission_params,
                "alpha": alpha_est,
                "beta": beta_est,
                "gamma": gamma_est,
                "rho": rho_est,
                "transition_params_source": transition_params,
            }
            cache.save(params_dict, key, "estimate_params")
            cache.set_flag(key, "estimate_params", params_dict)
        else:
            logger.info("Step 3/4: Using cached parameters")
            params_dict = cache.load(key, "estimate_params")
            emission_params = params_dict["emission_params"]
            alpha_est = params_dict.get("alpha", alpha)
            beta_est = params_dict.get("beta", beta)
            gamma_est = params_dict.get("gamma", gamma)
            rho_est = params_dict.get("rho", rho)

        result["steps_completed"].append("estimate_params")

        # Step 4: HMM inference
        if not cache.exists(key, "hmm_inference"):
            logger.info("Step 4/4: Running HMM inference")
            hmm = SisterChromatidHMM(
                alpha=alpha_est,
                beta=beta_est,
                gamma=gamma_est,
                rho=rho_est,
                emission_params=emission_params
            )

            positions = counts_df['position'].values
            results = hmm.decode_with_blocks(
                read_counts=counts_df,
                chromosome=chromosome or "all",
                positions=positions,
                sample_id=sample_id
            )

            # Convert unit-level coordinates back to SNP coordinates if binning was used
            if binning_method == "by_snp":
                logger.info(f"  Converting {len(results.states)} units back to {len(original_counts_df)} SNPs")
                results = convert_unit_coordinates_to_snp(
                    results,
                    original_counts_df,
                    n_snps=binning_snp_number
                )

            cache.save(results, key, "hmm_inference")
            cache.set_flag(key, "hmm_inference", {
                "n_snps": len(results.states),
                "n_blocks": len(results.blocks),
            })
        else:
            logger.info("Step 4/4: Using cached HMM inference")
            results = cache.load(key, "hmm_inference")

        result["steps_completed"].append("hmm_inference")
        result["results"] = results

        # Save output files
        blocks_file = output_path / f"{sample_id}_blocks.tsv"
        if chromosome:
            blocks_file = output_path / f"{sample_id}_{chromosome}_blocks.tsv"

        blocks_df = pd.DataFrame([
            {
                'chromosome': b.chromosome,
                'start_pos': b.start_pos,
                'end_pos': b.end_pos,
                'state': b.state,
                'confidence': b.confidence,
                'num_snps': b.num_snps,
            }
            for b in results.blocks
        ])
        blocks_df.to_csv(blocks_file, sep='\t', index=False)
        result["output_files"]["blocks"] = str(blocks_file)

        states_file = output_path / f"{sample_id}_states.tsv"
        if chromosome:
            states_file = output_path / f"{sample_id}_{chromosome}_states.tsv"

        states_df = pd.DataFrame({
            'chromosome': chromosome or "all",
            'position': results.positions,
            'state': results.states,
        })
        states_df.to_csv(states_file, sep='\t', index=False)
        result["output_files"]["states"] = str(states_file)

        stats_file = output_path / f"{sample_id}_statistics.json"
        if chromosome:
            stats_file = output_path / f"{sample_id}_{chromosome}_statistics.json"

        import json
        with open(stats_file, 'w') as f:
            json.dump(results.statistics, f, indent=2)
        result["output_files"]["statistics"] = str(stats_file)

        result["success"] = True
        logger.info(f"Completed processing sample {sample_id}")

    except Exception as e:
        logger.error(f"Error processing sample {sample_id}: {e}")
        result["error"] = str(e)

    return result


def detect_sister_chromatid_states_end2end(
    bam_file: str,
    snp_file: str,
    output_dir: str,
    sample_id: Optional[str] = None,
    chromosomes: Optional[List[str]] = None,
    cache_dir: Optional[str] = None,
    transition_params: str = "default",  # "default" | "estimated"
    alpha: float = 0.01,
    beta: float = 0.02,
    gamma: float = 0.001,
    rho: float = 0.1,
    depth_threshold: int = 2,
    other_threshold: float = 0.4,
    balance_threshold: float = 0.25,
    binning_method: Optional[str] = None,
    binning_snp_number: int = 250,
) -> Union[SisterChromatidResults, Dict[str, SisterChromatidResults]]:
    """
    End-to-end sister chromatid state inference from BAM file.

    This is the main entry point for the mII crossover detection pipeline.
    It processes BAM files through all steps:
    1. Extract read counts at SNP positions
    2. Binning (optional): merge adjacent SNPs into units
    3. Rough classification for parameter initialization
    4. Estimate emission and transition parameters
    5. HMM inference with Viterbi decoding
    6. Group states into contiguous blocks
    7. Convert unit-level coordinates back to SNP coordinates (if binning was used)

    Args:
        bam_file: Path to input BAM file
        snp_file: Path to SNP file (VCF or TSV)
        output_dir: Output directory for results
        sample_id: Sample ID (default: derived from bam_file name)
        chromosomes: Optional list of chromosomes to analyze (default: all)
        cache_dir: Cache directory (default: output_dir/.cache)
        transition_params: Transition parameter source ("default" | "estimated")
        alpha: PP<->PM transition probability (default: 0.01, used when transition_params="default")
        beta: PM->PP/MM transition probability (default: 0.02, used when transition_params="default")
        gamma: probability of entering U state (default: 0.001, used when transition_params="default")
        rho: probability of leaving U state (default: 0.1, used when transition_params="default")
        depth_threshold: Minimum depth for confident call (default: 2)
        other_threshold: Maximum other fraction for confident call (default: 0.4)
        balance_threshold: Threshold for classifying as balanced (default: 0.25)
                          PP: ref_frac > 0.75, MM: alt_frac > 0.75
        binning_method: Binning method (None for single-SNP mode, "by_snp" for merging by SNP count)
        binning_snp_number: Number of SNPs per unit (default: 250)

    Returns:
        SisterChromatidResults object (single chromosome) or dictionary (multiple chromosomes)

    Example:
        >>> # Default parameters (fixed transition matrix)
        >>> results = detect_sister_chromatid_states_end2end(
        ...     bam_file="mII_cell.bam",
        ...     snp_file="snps.vcf",
        ...     output_dir="output/"
        ... )

        >>> # Estimated transition parameters from data
        >>> results = detect_sister_chromatid_states_end2end(
        ...     bam_file="mII_cell.bam",
        ...     snp_file="snps.vcf",
        ...     output_dir="output/",
        ...     transition_params="estimated"
        ... )

        >>> # Custom binning (250 SNPs per unit)
        >>> results = detect_sister_chromatid_states_end2end(
        ...     bam_file="mII_cell.bam",
        ...     snp_file="snps.vcf",
        ...     output_dir="output/",
        ...     binning_method="by_snp",
        ...     binning_snp_number=250
        ... )
    """
    if sample_id is None:
        sample_id = Path(bam_file).stem

    if not Path(bam_file).exists():
        raise FileNotFoundError(f"BAM file not found: {bam_file}")
    if not Path(snp_file).exists():
        raise FileNotFoundError(f"SNP file not found: {snp_file}")

    if cache_dir is None:
        cache_dir = str(Path(output_dir) / ".cache")

    logger.info(f"Starting sister chromatid inference pipeline")
    logger.info(f"  Sample: {sample_id}")
    logger.info(f"  BAM file: {bam_file}")
    logger.info(f"  SNPs: {snp_file}")
    logger.info(f"  Output: {output_dir}")

    if chromosomes is None:
        # Process all chromosomes together
        result = _process_sample_mii(
            sample_id=sample_id,
            bam_file=bam_file,
            snp_file=snp_file,
            output_dir=output_dir,
            cache_dir=cache_dir,
            chromosome=None,
            transition_params=transition_params,
            alpha=alpha,
            beta=beta,
            gamma=gamma,
            rho=rho,
            depth_threshold=depth_threshold,
            other_threshold=other_threshold,
            balance_threshold=balance_threshold,
            binning_method=binning_method,
            binning_snp_number=binning_snp_number,
        )

        if not result["success"]:
            raise RuntimeError(f"Pipeline failed: {result.get('error', 'Unknown error')}")

        return result["results"]

    else:
        # Process each chromosome separately
        all_results = {}
        for chrom in chromosomes:
            logger.info(f"Processing chromosome {chrom}")
            result = _process_sample_mii(
                sample_id=sample_id,
                bam_file=bam_file,
                snp_file=snp_file,
                output_dir=output_dir,
                cache_dir=cache_dir,
                chromosome=chrom,
                transition_params=transition_params,
                alpha=alpha,
                beta=beta,
                gamma=gamma,
                rho=rho,
                depth_threshold=depth_threshold,
                other_threshold=other_threshold,
                balance_threshold=balance_threshold,
                binning_method=binning_method,
                binning_snp_number=binning_snp_number,
            )

            if result["success"]:
                all_results[chrom] = result["results"]
                logger.info(f"  ✓ {chrom}: {len(result['results'].blocks)} blocks")
            else:
                logger.error(f"  ✗ {chrom}: {result.get('error', 'Unknown error')}")

        return all_results


def save_results_mii(
    results: Union[SisterChromatidResults, Dict[str, SisterChromatidResults]],
    output_dir: str,
    sample_id: str
) -> Dict[str, str]:
    """
    Save sister chromatid inference results to files.

    Args:
        results: Results object or dictionary
        output_dir: Output directory
        sample_id: Sample identifier

    Returns:
        Dictionary mapping output type to file path
    """
    from pathlib import Path

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    output_files = {}

    if isinstance(results, dict):
        for chrom, result in results.items():
            chrom_output = save_results_mii(result, output_dir, f"{sample_id}_{chrom}")
            output_files.update(chrom_output)
    else:
        # Save blocks
        blocks_file = output_path / f"{sample_id}_blocks.tsv"
        blocks_df = pd.DataFrame([
            {
                'chromosome': b.chromosome,
                'start_pos': b.start_pos,
                'end_pos': b.end_pos,
                'state': b.state,
                'confidence': b.confidence,
                'num_snps': b.num_snps,
            }
            for b in results.blocks
        ])
        blocks_df.to_csv(blocks_file, sep='\t', index=False)
        output_files['blocks'] = str(blocks_file)

        # Save states
        states_file = output_path / f"{sample_id}_states.tsv"
        states_df = pd.DataFrame({
            'chromosome': results.chromosome,
            'position': results.positions,
            'state': results.states,
        })
        states_df.to_csv(states_file, sep='\t', index=False)
        output_files['states'] = str(states_file)

        # Save statistics
        import json
        stats_file = output_path / f"{sample_id}_statistics.json"
        with open(stats_file, 'w') as f:
            json.dump(results.statistics, f, indent=2)
        output_files['statistics'] = str(stats_file)

    return output_files
