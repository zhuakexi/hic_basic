"""
Input/output operations for crossover detection.

This module handles:
- Loading SNP data from various formats (TSV, VCF, parquet)
- Loading phasing information
- Saving crossover detection results
- Integration with hic_basic.phasing module formats
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
import pysam

from .detector import Crossover, CrossoverResults, HaplotypeBlock

logger = logging.getLogger(__name__)


def load_snp_data(
    snp_file: str,
    sample_id: Optional[str] = None,
    chromosome: Optional[str] = None
) -> pd.DataFrame:
    """
    Load SNP genotype data from various formats.

    Supported formats:
        - TSV/CSV: Tab or comma-separated with columns
        - VCF: Variant Call Format
        - Parquet: Columnar storage format

    Args:
        snp_file: Path to SNP data file
        sample_id: Sample ID to extract (for VCF with multiple samples)
        chromosome: Optional chromosome to filter

    Returns:
        DataFrame with columns: chromosome, position, snp_id, ref_allele,
                               alt_allele, sperm_genotype, diploid_genotype

    Raises:
        ValueError: If file format is not supported
        FileNotFoundError: If file does not exist
    """
    snp_path = Path(snp_file)

    if not snp_path.exists():
        raise FileNotFoundError(f"SNP file not found: {snp_file}")

    suffix = snp_path.suffix.lower()

    if suffix in ['.tsv', '.txt', '.csv']:
        df = _load_snp_tsv(snp_file, chromosome)
    elif suffix == '.vcf':
        df = _load_snp_vcf(snp_file, sample_id, chromosome)
    elif suffix == '.parquet':
        df = _load_snp_parquet(snp_file, chromosome)
    elif suffix == '.gz':
        if str(snp_path).endswith('.tsv.gz'):
            df = _load_snp_tsv(snp_file, chromosome)
        elif str(snp_path).endswith('.vcf.gz'):
            df = _load_snp_vcf(snp_file, sample_id, chromosome)
        else:
            raise ValueError(f"Unsupported compressed format: {snp_file}")
    else:
        raise ValueError(f"Unsupported file format: {suffix}")

    logger.info(f"Loaded {len(df)} SNPs from {snp_file}")
    return df


def _load_snp_tsv(
    snp_file: str,
    chromosome: Optional[str] = None
) -> pd.DataFrame:
    """Load SNP data from TSV/CSV file."""
    df = pd.read_csv(snp_file, sep='\t' if '\t' in open(snp_file).readline() else ',')

    required_cols = ['chromosome', 'position', 'sperm_genotype']
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    if chromosome is not None:
        df = df[df['chromosome'] == chromosome]

    df = df.sort_values(['chromosome', 'position']).reset_index(drop=True)

    if 'snp_id' not in df.columns:
        df['snp_id'] = df['chromosome'] + ':' + df['position'].astype(str)

    if 'ref_allele' not in df.columns:
        df['ref_allele'] = '.'
    if 'alt_allele' not in df.columns:
        df['alt_allele'] = '.'
    if 'diploid_genotype' not in df.columns:
        df['diploid_genotype'] = './.'

    return df


def _load_snp_vcf(
    vcf_file: str,
    sample_id: Optional[str] = None,
    chromosome: Optional[str] = None
) -> pd.DataFrame:
    """Load SNP data from VCF file."""
    records = []

    with pysam.VariantFile(vcf_file) as vcf:
        samples = list(vcf.header.samples)

        if sample_id is None:
            if len(samples) == 1:
                sample_id = samples[0]
            else:
                raise ValueError("Multiple samples in VCF, please specify sample_id")

        if sample_id not in samples:
            raise ValueError(f"Sample {sample_id} not found in VCF. Available: {samples}")

        for record in vcf.fetch(chromosome) if chromosome else vcf:
            gt = record.samples[sample_id].get('GT')
            if gt is None:
                sperm_gt = 'N'
            else:
                sperm_gt = '/'.join(str(a) if a is not None else '.' for a in gt)

            records.append({
                'chromosome': record.chrom,
                'position': record.pos,
                'snp_id': record.id or f"{record.chrom}:{record.pos}",
                'ref_allele': record.ref,
                'alt_allele': ','.join(record.alts) if record.alts else '.',
                'sperm_genotype': sperm_gt,
                'diploid_genotype': './.',
            })

    df = pd.DataFrame(records)

    if chromosome is not None:
        df = df[df['chromosome'] == chromosome]

    df = df.sort_values(['chromosome', 'position']).reset_index(drop=True)
    return df


def _load_snp_parquet(
    parquet_file: str,
    chromosome: Optional[str] = None
) -> pd.DataFrame:
    """Load SNP data from parquet file."""
    df = pd.read_parquet(parquet_file)

    required_cols = ['chromosome', 'position', 'sperm_genotype']
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    if chromosome is not None:
        df = df[df['chromosome'] == chromosome]

    df = df.sort_values(['chromosome', 'position']).reset_index(drop=True)

    if 'snp_id' not in df.columns:
        df['snp_id'] = df['chromosome'] + ':' + df['position'].astype(str)

    return df


def load_phasing_data(
    phasing_file: str,
    chromosome: Optional[str] = None
) -> pd.DataFrame:
    """
    Load phasing information from file.

    Supported formats:
        - TSV/CSV: Tab or comma-separated
        - Parquet: Columnar storage
        - From hic_basic.phasing module output

    Args:
        phasing_file: Path to phasing data file
        chromosome: Optional chromosome to filter

    Returns:
        DataFrame with columns: chromosome, position, snp_id,
                               haplotype_1, haplotype_2

    Raises:
        ValueError: If file format is not supported
        FileNotFoundError: If file does not exist
    """
    phasing_path = Path(phasing_file)

    if not phasing_path.exists():
        raise FileNotFoundError(f"Phasing file not found: {phasing_file}")

    suffix = phasing_path.suffix.lower()

    if suffix in ['.tsv', '.txt', '.csv', '.gz']:
        df = pd.read_csv(
            phasing_file,
            sep='\t' if str(phasing_path).endswith('.tsv') or '\t' in open(phasing_file).readline() else ','
        )
    elif suffix == '.parquet':
        df = pd.read_parquet(phasing_file)
    else:
        raise ValueError(f"Unsupported phasing file format: {suffix}")

    required_cols = ['chromosome', 'position', 'haplotype_1', 'haplotype_2']
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        if {'chromosome', 'position', 'ref', 'alt'}.issubset(df.columns):
            df['haplotype_1'] = df['ref']
            df['haplotype_2'] = df['alt']
        else:
            raise ValueError(f"Missing required columns: {missing}")

    if chromosome is not None:
        df = df[df['chromosome'] == chromosome]

    df = df.sort_values(['chromosome', 'position']).reset_index(drop=True)

    if 'snp_id' not in df.columns:
        df['snp_id'] = df['chromosome'] + ':' + df['position'].astype(str)

    logger.info(f"Loaded phasing for {len(df)} SNPs from {phasing_file}")
    return df


def load_phased_snp_from_hic_basic(
    phased_snp_file: str,
    chromosome: Optional[str] = None
) -> pd.DataFrame:
    """
    Load phased SNP data from hic_basic.phasing module format.

    This function reads the output format from hic_basic.phasing.parse_phased_snp
    and converts it to a DataFrame suitable for crossover detection.

    Args:
        phased_snp_file: Path to phased SNP file (TSV or parquet)
        chromosome: Optional chromosome to filter

    Returns:
        DataFrame with columns: chromosome, position, snp_id,
                               haplotype_1, haplotype_2
    """
    if phased_snp_file.endswith('.parquet'):
        df = pd.read_parquet(phased_snp_file)
    else:
        df = pd.read_csv(phased_snp_file, sep='\t')

    if len(df.columns) >= 4:
        df = df.iloc[:, :4]
        df.columns = ['chromosome', 'position', 'haplotype_1', 'haplotype_2']

    if chromosome is not None:
        df = df[df['chromosome'] == chromosome]

    df = df.sort_values(['chromosome', 'position']).reset_index(drop=True)
    df['snp_id'] = df['chromosome'] + ':' + df['position'].astype(str)

    return df


def save_results(
    results: CrossoverResults,
    output_dir: str,
    prefix: Optional[str] = None
) -> Dict[str, str]:
    """
    Save crossover detection results to files.

    Args:
        results: CrossoverResults object
        output_dir: Output directory path
        prefix: Optional prefix for output files

    Returns:
        Dictionary mapping output type to file path
    """
    from pathlib import Path

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    sample_id = results.sperm_id
    file_prefix = prefix if prefix else sample_id

    output_files = {}

    blocks_file = output_path / f"{file_prefix}_haplotype_blocks.tsv"
    blocks_df = pd.DataFrame([
        {
            'chromosome': b.chromosome,
            'start_pos': b.start_pos,
            'end_pos': b.end_pos,
            'haplotype': b.haplotype,
            'confidence': b.confidence,
            'num_snps': b.num_snps,
        }
        for b in results.haplotype_blocks
    ])
    if len(blocks_df) > 0:
        blocks_df.to_csv(blocks_file, sep='\t', index=False)
        output_files['haplotype_blocks'] = str(blocks_file)
        logger.info(f"Saved haplotype blocks to {blocks_file}")

    crossovers_file = output_path / f"{file_prefix}_crossovers.tsv"
    crossovers_df = pd.DataFrame([
        {
            'sperm_id': c.sperm_id,
            'chromosome': c.chromosome,
            'position': c.position,
            'start_pos': c.start_pos,
            'end_pos': c.end_pos,
            'from_haplotype': c.from_haplotype,
            'to_haplotype': c.to_haplotype,
            'confidence': c.confidence,
            'resolution_kb': c.resolution_kb,
        }
        for c in results.crossovers
    ])
    if len(crossovers_df) > 0:
        crossovers_df.to_csv(crossovers_file, sep='\t', index=False)
        output_files['crossovers'] = str(crossovers_file)
        logger.info(f"Saved crossovers to {crossovers_file}")

    stats_file = output_path / f"{file_prefix}_statistics.tsv"
    stats_df = pd.DataFrame([results.statistics])
    stats_df.to_csv(stats_file, sep='\t', index=False)
    output_files['statistics'] = str(stats_file)
    logger.info(f"Saved statistics to {stats_file}")

    return output_files


def merge_snp_with_phasing(
    snp_df: pd.DataFrame,
    phasing_df: pd.DataFrame,
    on_cols: Optional[List[str]] = None
) -> pd.DataFrame:
    """
    Merge SNP data with phasing information.

    Args:
        snp_df: SNP genotype DataFrame
        phasing_df: Phasing information DataFrame
        on_cols: Columns to merge on (default: ['chromosome', 'position', 'snp_id'])

    Returns:
        Merged DataFrame with both genotype and phasing information
    """
    if on_cols is None:
        on_cols = ['chromosome', 'position', 'snp_id']

    available_cols = [c for c in on_cols if c in snp_df.columns and c in phasing_df.columns]

    if len(available_cols) < 2:
        raise ValueError(f"Insufficient common columns for merge: {available_cols}")

    merged = snp_df.merge(
        phasing_df,
        on=available_cols,
        how='left',
        suffixes=('', '_phasing')
    )

    logger.info(f"Merged {len(merged)} SNPs with phasing information")
    return merged


def extract_genotype_arrays(
    df: pd.DataFrame
) -> Tuple[np.ndarray, List[str], List[str], List[str]]:
    """
    Extract arrays needed for HMM from merged DataFrame.

    Args:
        df: Merged DataFrame with sperm_genotype, haplotype_1, haplotype_2 columns

    Returns:
        Tuple of (positions, sperm_genotypes, haplotype_1, haplotype_2)
    """
    positions = df['position'].values
    sperm_genotypes = df['sperm_genotype'].fillna('N').tolist()
    haplotype_1 = df['haplotype_1'].fillna('N').tolist()
    haplotype_2 = df['haplotype_2'].fillna('N').tolist()

    return positions, sperm_genotypes, haplotype_1, haplotype_2


def validate_input_data(
    snp_df: pd.DataFrame,
    phasing_df: Optional[pd.DataFrame] = None
) -> Tuple[bool, List[str]]:
    """
    Validate input data for crossover detection.

    Args:
        snp_df: SNP genotype DataFrame
        phasing_df: Optional phasing DataFrame

    Returns:
        Tuple of (is_valid, list of issues)
    """
    issues = []

    required_snp_cols = ['chromosome', 'position', 'sperm_genotype']
    missing_snp = [c for c in required_snp_cols if c not in snp_df.columns]
    if missing_snp:
        issues.append(f"Missing SNP columns: {missing_snp}")

    if len(snp_df) == 0:
        issues.append("SNP DataFrame is empty")

    if 'position' in snp_df.columns and snp_df['position'].isnull().any():
        issues.append("SNP positions contain null values")

    if phasing_df is not None:
        required_phasing_cols = ['chromosome', 'position', 'haplotype_1', 'haplotype_2']
        missing_phasing = [c for c in required_phasing_cols if c not in phasing_df.columns]
        if missing_phasing:
            issues.append(f"Missing phasing columns: {missing_phasing}")

        if 'chromosome' in snp_df.columns and 'chromosome' in phasing_df.columns:
            common_chroms = set(snp_df['chromosome'].unique()) & set(phasing_df['chromosome'].unique())
            if len(common_chroms) == 0:
                issues.append("No common chromosomes between SNP and phasing data")

    if len(issues) == 0:
        logger.info("Input data validation passed")

    return len(issues) == 0, issues
