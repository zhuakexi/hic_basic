import random
from typing import Optional, Union, Tuple, List

import pysam
import numpy as np
import plotly.graph_objects as go
from pybedtools import BedTool
from tqdm import tqdm
from .plot.general import _histogram

def calculate_fragment_lengths(
    bam_path: str,
    sampling: Optional[Union[float, int]] = None
) -> np.ndarray:
    """
    Calculate fragment lengths from a BAM file with optional sampling.

    Parameters:
        bam_path (str): Path to the BAM file.
        sampling (Optional[Union[float, int]]): 
            If float, sample this proportion of reads (0 < sampling <= 1).
            If int, sample this number of reads.
            If None, no sampling (default).

    Returns:
        np.ndarray: Array of fragment lengths.
    """
    fragment_lengths = []

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        total_reads = bam.count()
        print(f"Total reads in BAM file: {total_reads}")
        bam.reset()

        # Determine sampling strategy
        sample_size = total_reads
        sample_ratio = 1.0
        use_sampling = False

        if sampling is not None:
            if isinstance(sampling, float):
                if not (0 < sampling <= 1):
                    raise ValueError("sampling ratio must be between 0 and 1")
                sample_ratio = sampling
                sample_size = int(total_reads * sample_ratio)
                use_sampling = True
            elif isinstance(sampling, int):
                if sampling < 0:
                    raise ValueError("sampling size must be a non-negative integer")
                sample_size = sampling
                use_sampling = True
            else:
                raise TypeError("sampling must be a float or an int")

        # Apply sampling and process reads
        if use_sampling and isinstance(sampling, int):
            # Reservoir Sampling for fixed number of reads
            sample_pool = []
            for i, read in enumerate(tqdm(bam, desc="Sampling reads", total=total_reads)):
                if i < sample_size:
                    sample_pool.append(read)
                else:
                    r = random.randint(0, i)
                    if r < sample_size:
                        sample_pool[r] = read

            # Process sampled reads
            for read in tqdm(sample_pool, desc="Processing sampled reads", total=len(sample_pool)):
                if not read.is_unmapped and not read.mate_is_unmapped:
                    fragment_length = abs(read.template_length)
                    fragment_lengths.append(fragment_length)

        else:
            # Probability sampling or no sampling
            for read in tqdm(bam, desc="Processing reads", total=total_reads):
                if not use_sampling or random.random() < sample_ratio:
                    if not read.is_unmapped and not read.mate_is_unmapped:
                        fragment_length = abs(read.template_length)
                        fragment_lengths.append(fragment_length)

    return np.array(fragment_lengths)


def plot_fragment_length_distribution(
    fragment_lengths: np.ndarray,
    **kwargs
) -> go.Figure:
    """
    Plot fragment length distribution using Plotly.

    Input:
        fragment_lengths (np.ndarray): Array of fragment lengths.
        **kwargs: Additional keyword arguments for histogram plotting.

    Output:
        go.Figure: Plotly figure object.
    """
    fig = _histogram(
        fragment_lengths,
        **kwargs
    )
    fig.update_layout(
        title="ATAC-seq Fragment Length Distribution",
        xaxis_title="Fragment Length (bp)",
        yaxis_title="Count",
        bargap=0.1,
        showlegend=False
    )
    return fig


# def calculate_tss_enrichment(
#     bam_path: str,
#     gtf_path: str,
#     tss_window: int = 1000,
#     gene_body_window: int = 3000,
#     bed_cache: Optional[str] = None,
#     chrom_size_file = "hg38.chrom.sizes"
# ) -> tuple[float, float, float]:
#     """
#     Calculate TSS region and gene body region read depths.

#     Input:
#         bam_path (str): Path to the BAM file.
#         gtf_path (str): Path to the GTF annotation file.
#         tss_window (int): Window size around TSS (default: 1000 bp).
#         gene_body_window (int): Window size for gene body region (default: 3000 bp).

#     Output:
#         tuple[float, float, float]: (TSS_depth, gene_body_depth, tss_score)
#     """
#     # Extract TSS positions from GTF
#     tss_bed = BedTool(gtf_path).filter(lambda x: x[2] == "transcript").each(
#         lambda x: [x[0], int(x[3]) -1, int(x[3]), x[-1]]
#     )
#     tss_bed = tss_bed.saveas(bed_cache) if bed_cache else tss_bed

#     # Extend TSS and gene body regions
#     tss_flank = tss_bed.flank(l=tss_window, r=0, s=True, g=chrom_size_file)
#     gene_body = tss_bed.slop(b=gene_body_window, l=0, r=0, s=True, g=chrom_size_file)

#     # Calculate read depths
#     tss_coverage = pysam.depth("-b", tss_flank.fn, bam_path).splitlines()
#     gene_body_coverage = pysam.depth("-b", gene_body.fn, bam_path).splitlines()

#     tss_depth = np.mean([int(line.split()[2]) for line in tss_coverage])
#     gene_body_depth = np.mean([int(line.split()[2]) for line in gene_body_coverage])
#     tss_score = tss_depth / gene_body_depth if gene_body_depth > 0 else 0.0

#     return tss_depth, gene_body_depth, tss_score

def get_tss_regions(tss_bed_path=None, chrom_size_path=None, window_size: int = 2000):
    """
    Expands TSS regions by a specified window size (default: 2000 bp upstream and downstream).

    Input
       tss_bed_path (str): Path to the TSS BED file (1-based coordinates).
       chrom_size_path (str): Path to the chromosome sizes file (2-column format).
       window_size (int): Number of bases to extend on each side of the TSS.

    Output
       BedTool: Expanded TSS regions as a pybedtools BedTool object.
    """
    valid_chroms = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
    tss_bed = BedTool(tss_bed_path)
    tss_bed = tss_bed.filter(lambda x: x.chrom in valid_chroms)
    tss_ext = tss_bed.slop(b=window_size, g=chrom_size_path)
    return tss_ext

import numpy as np
import random
import pysam
from typing import Tuple, List

def calculate_coverage(
    bam_path: str, 
    tss_regions, 
    read_length: int, 
    sampling=None
):
    """
    Computes coverage across TSS regions using pysam.depth, with optional sampling.

    Input
       bam_path (str): Path to the BAM file (must be indexed).
       tss_regions (BedTool): Expanded TSS regions (output of `get_tss_regions`).
       read_length (int): Length of sequencing reads (used for alignment centering).
       sampling (float or int, optional): Sampling strategy.
           - If float: Fraction of regions to sample (0 < value <= 1).
           - If int: Exact number of regions to sample.
           - If None: Process all regions.

    Output
       Tuple[np.ndarray, List[float]]:
           coverage_array (np.ndarray): 2D array of coverage values (shape: [n_TSS, bins]).
           distances (List[float]): Distance from TSS for each bin (in base pairs).
    """
    # Validate sampling parameter
    if sampling is not None:
        total_regions = len(tss_regions)
        if isinstance(sampling, float):
            if not (0 < sampling <= 1):
                raise ValueError("sampling ratio must be between 0 and 1")
            n_samples = int(total_regions * sampling)
        elif isinstance(sampling, int):
            if sampling < 0:
                raise ValueError("sampling size must be a non-negative integer")
            n_samples = sampling
        else:
            raise TypeError("sampling must be a float or an int")

        # Ensure sample count is valid
        n_samples = min(n_samples, total_regions)
        n_samples = max(n_samples, 0)

        # Randomly sample regions
        sampled_regions = random.sample(list(tss_regions), n_samples)
    else:
        sampled_regions = tss_regions

    # Process sampled regions
    coverage_list = []
    for interval in tqdm(sampled_regions, desc="Calculating coverage", total=len(sampled_regions)):
        chrom, start, end = interval.chrom, int(interval.start), int(interval.end)
        depth = pysam.depth("-b", "-", bam_path, input=f"{chrom}\t{start}\t{end}")
        depths = [int(line.split()[2]) for line in depth.strip().splitlines()]
        coverage_list.append(depths)

    # Convert to numpy array
    coverage_array = np.array(coverage_list)
    distances = np.linspace(-2000, 2000, coverage_array.shape[1]).tolist()
    return coverage_array, distances
def normalize_coverage(coverage_array: np.ndarray, read_length: int, greenleaf_norm: bool = True) -> np.ndarray:
    """
    Normalizes coverage using Greenleaf-style or per million mapped reads normalization.

    Input:
       coverage_array (np.ndarray): 2D coverage array (shape: [n_TSS, bins]).
       read_length (int): Read length (used for centering adjustment).
       greenleaf_norm (bool): Whether to apply Greenleaf-style normalization (default: True).

    Output:
       np.ndarray: Normalized coverage array.
    """
    if greenleaf_norm:
        edge_bins = int(100 / (2 * 2000 / coverage_array.shape[1]))
        means = np.mean(coverage_array, axis=0)
        avg_noise = (sum(means[:edge_bins]) + sum(means[-edge_bins:])) / (2 * edge_bins)
        normalized = coverage_array / avg_noise
    else:
        total_mapped = pysam.AlignmentFile(bam_path, "rb").count()
        normalized = coverage_array / (total_mapped / 1e6)
    return normalized

def process_tss_enrichment(
    bam_path: str,
    tss_bed_path: str,
    chrom_size_path: str,
    read_length: int,
    window_size: int = 2000,
    greenleaf_norm: bool = True
) -> Tuple[np.ndarray, float]:
    """
    Full TSS enrichment analysis pipeline: region expansion, coverage calculation, and normalization.

    Input:
       bam_path (str): Path to the BAM file.
       tss_bed_path (str): Path to the TSS BED file.
       chrom_size_path (str): Path to the chromosome sizes file.
       read_length (int): Sequencing read length.
       window_size (int): TSS region extension size (default: 2000).
       greenleaf_norm (bool): Use Greenleaf-style normalization (default: True).

    Output:
       Tuple[np.ndarray, float]:
           normalized_array (np.ndarray): Normalized coverage array (shape: [n_TSS, bins]).
           max_signal (float): Maximum signal value across all TSS regions.
    """
    tss_ext = get_tss_regions(tss_bed_path, chrom_size_path, window_size)
    coverage_array, _ = calculate_coverage(bam_path, tss_ext, read_length)
    normalized_array = normalize_coverage(coverage_array, read_length, greenleaf_norm)
    max_signal = np.max(np.mean(normalized_array, axis=0))
    return normalized_array, max_signal

