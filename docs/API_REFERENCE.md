# hic_basic API Reference

*Auto-generated API documentation for public functions*

Total modules with public functions: 76

---

## Table of Contents

- [DI](#di) (18 functions)
- [Ps](#ps) (5 functions)
- [TAD.IS](#tad-is) (2 functions)
- [TAD.tad](#tad-tad) (8 functions)
- [ana](#ana) (12 functions)
- [atac](#atac) (6 functions)
- [basic](#basic) (5 functions)
- [binnify](#binnify) (1 functions)
- [calculate](#calculate) (4 functions)
- [cli.download](#cli-download) (2 functions)
- [cli.render](#cli-render) (3 functions)
- [cli.version](#cli-version) (2 functions)
- [clustering](#clustering) (1 functions)
- [compare_pixels](#compare_pixels) (5 functions)
- [compartment](#compartment) (8 functions)
- [coolstuff](#coolstuff) (32 functions)
- [cycle_phasing](#cycle_phasing) (2 functions)
- [data](#data) (18 functions)
- [embedding](#embedding) (5 functions)
- [feature](#feature) (10 functions)
- [genome](#genome) (11 functions)
- [heuristic_ordering](#heuristic_ordering) (1 functions)
- [hicbr](#hicbr) (6 functions)
- [hicio](#hicio) (31 functions)
- [impute.schicluster](#impute-schicluster) (3 functions)
- [impute.simpute](#impute-simpute) (9 functions)
- [inter_contact](#inter_contact) (3 functions)
- [interpolate](#interpolate) (3 functions)
- [lap_embedding](#lap_embedding) (3 functions)
- [mchr](#mchr) (7 functions)
- [metrics](#metrics) (4 functions)
- [nagano_cycle_phasing](#nagano_cycle_phasing) (6 functions)
- [paracalc](#paracalc) (2 functions)
- [phasing](#phasing) (22 functions)
- [pileup](#pileup) (6 functions)
- [pipeline.rule](#pipeline-rule) (8 functions)
- [plot.bar](#plot-bar) (1 functions)
- [plot.general](#plot-general) (8 functions)
- [plot.hic](#plot-hic) (18 functions)
- [plot.pca](#plot-pca) (2 functions)
- [plot.plot](#plot-plot) (15 functions)
- [plot.render](#plot-render) (27 functions)
- [plot.rna](#plot-rna) (3 functions)
- [plot.scanpy](#plot-scanpy) (5 functions)
- [plot.scatter](#plot-scatter) (6 functions)
- [plot.tech](#plot-tech) (1 functions)
- [plot.track](#plot-track) (4 functions)
- [plot.utils](#plot-utils) (20 functions)
- [pp](#pp) (1 functions)
- [pseudotime.TI](#pseudotime-ti) (1 functions)
- [pymol.color_centelo](#pymol-color_centelo) (2 functions)
- [pymol.glow](#pymol-glow) (1 functions)
- [pymol.sele_cent](#pymol-sele_cent) (2 functions)
- [pymol.sele_telo](#pymol-sele_telo) (2 functions)
- [pymol.utils](#pymol-utils) (2 functions)
- [scAB_embedding](#scab_embedding) (9 functions)
- [scripts.downsra](#scripts-downsra) (1 functions)
- [scripts.rescue_jupyter](#scripts-rescue_jupyter) (2 functions)
- [sequence](#sequence) (11 functions)
- [shuffle](#shuffle) (4 functions)
- [spectral_ordering](#spectral_ordering) (1 functions)
- [structure.align](#structure-align) (1 functions)
- [structure.measure](#structure-measure) (15 functions)
- [structure.pileup](#structure-pileup) (16 functions)
- [structure.utils](#structure-utils) (6 functions)
- [territory](#territory) (3 functions)
- [utils](#utils) (8 functions)
- [wet.adtools](#wet-adtools) (7 functions)
- [wet.afbb](#wet-afbb) (18 functions)
- [wet.exp_record](#wet-exp_record) (10 functions)
- [wet.gam](#wet-gam) (1 functions)
- [wet.meta_trick](#wet-meta_trick) (3 functions)
- [wet.paracalc](#wet-paracalc) (6 functions)
- [wet.qc](#wet-qc) (3 functions)
- [wet.rawcheck](#wet-rawcheck) (7 functions)
- [wet.utils](#wet-utils) (1 functions)

---


## DI

**File:** `DI.py`

**Public functions:** 18


### `N_partitions(max_dist=2000000.0, binsize=20000.0, chrom_mean_length=100000000.0, chunksize=100000)`


Pick a number of partitions suitable for the data.
Input:
    max_dist: max distance in impuation
    binsize: resolution of the data
    chrom_mean_lengths: mean length of chromosomes
    chunksize: number of rows suitable in memory, 100e3 ~ 6MB
Output:
    number of partitions


**Source:** Line 104 in [DI.py](DI.py#L104)


---


### `gen_design_matrix(imputed_paths, sample_names=None, fo=None, binsize=20000, max_dist=None)`


Horizontal concat all imputed data, with pixel_id as index
Input:
    imputed_paths: list of imputed data paths
    sample_names: list of sample names, if None, derive from imputed_paths
    fo: output parquet file, write to fo if not None(will trigger computation)
        if None, return dask dataframe
    binsize: resolution of the data
    max_dist: max distance in impuation
Output:
    dd.DataFrame: n * m dask dataframe
        (chrom, start1, start2, sample1, sample2, ...)


**Source:** Line 425 in [DI.py](DI.py#L425)


---


### `get_sample_name(imputed_path)`


Get sample name from file with many possible suffix.
Input:
    imputed_path: path to imputed data
Output:
    sample name


**Source:** Line 416 in [DI.py](DI.py#L416)


---


### `imputed_reader(imputed_path, sample_id, genome, binsize, max_dist)`


No documentation available.


**Source:** Line 274 in [DI.py](DI.py#L274)


---


### `imputed_reader_strID(imputed_path, sample_id, binsize=20000, max_dist=None)`


Read in imputed data, zscore normalize, and add pixel_id.
Note: this function treat imputation as intra-chromosomal.
It assumes all row has the same chrom1, chrom2. It then do groupby only on band,
and code pixel_id as chrom1:start1-start2.
Input:
    imputed_path: path to imputed data.
        It should be a parquet file with columns chrom1, start1, start2, value,
        or chrom, start1, start2, value. Value col name is not assumed.
        More columns are allowed and ignored.
    sample_id: sample id, used as output column name
    binsize: resolution of the data
    max_dist: max distance in impuation
Output:
    dd.Series: zscore normalized distance


**Source:** Line 372 in [DI.py](DI.py#L372)


---


### `mannwhitneyu_row(row, m1, m2)`


No documentation available.


**Source:** Line 88 in [DI.py](DI.py#L88)


---


### `mannwhitneyu_row_labeled(row, groupA, groupB)`


No documentation available.


**Source:** Line 467 in [DI.py](DI.py#L467)


---


### `mannwhitneyu_test_pandas(task_name, combined_df, groupA, groupB, max_3d_dist, withdf=False)`


No documentation available.


**Source:** Line 472 in [DI.py](DI.py#L472)


---


### `multiple_testing_correction(pvalues, correction_type='FDR')`


Consistent with R - print
correct_pvalues_for_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05,
                                      0.069, 0.07, 0.071, 0.09, 0.1])
from https://github.com/CoBiG2/cobig_misc_scripts/blob/master/FDR.py


**Source:** Line 50 in [DI.py](DI.py#L50)


---


### `normalize_band(df, chrom=None, max_dist=2000000, binsize=20000)`


Input:
    df: bedpe like dask dataframe, with columns chrom1, start1, end1, chrom2, start2, end2, value
        don't assume value name
    chrom: if not None, assume df has only intra pixels
    max_dist: max distance in impuation, used to estimiate partitions
    binsize: resolution of the data, used to estimiate partitions
Output:
    df: same as input, with value normalized by band
    TODO: if max_dist or binsize is None, df is partitioned by chrom1, chrom2


**Source:** Line 13 in [DI.py](DI.py#L13)


---


### `pileup_impute(imputes, genome, binsize, chrom, start, end, zscore=False, symmetric=True)`


No documentation available.


**Source:** Line 580 in [DI.py](DI.py#L580)


---


### `plot_DI(genome, chrom, start, end, binsize, DI)`


No documentation available.


**Source:** Line 627 in [DI.py](DI.py#L627)


---


### `plot_compare_impute(imputesA, imputesB, genome, binsize, chrom, start, end, zscore=True, subplot_titles=('A', 'B', 'A-B'))`


No documentation available.


**Source:** Line 642 in [DI.py](DI.py#L642)


---


### `process_file(imputed_path, chrom, genome, max_dist, binsize, chunksize=10000)`


No documentation available.


**Source:** Line 122 in [DI.py](DI.py#L122)


---


### `project_DI(genome, chrom, start, end, binsize, DI)`


No documentation available.


**Source:** Line 566 in [DI.py](DI.py#L566)


---


### `read_pvalues(pval_path, chrom)`


No documentation available.


**Source:** Line 526 in [DI.py](DI.py#L526)


---


### `simpleDiff(groupA, groupB, chrom='chr1', genome='mm10', fo=None, max_3d_dist=5, max_dist=2000000, binsize=20000, chunksize=10000, cache=None)`


Input:
    groupA: imputed .parquet file path for groupA
    groupB: imputed .parquet file paths for groupB
    fo: output parquet file
    filt_fdr: filter pixels with fdr > filt_fdr, if None, no filtering
    n_jobs: number of jobs
Output:
    df: n * (6+3) dask dataframe
        (chrom1, start1, end1, chrom2, start2, end2, meanA, meanB, diff, pvalue, qvalue)


**Source:** Line 210 in [DI.py](DI.py#L210)


---


### `simpleDiff_postprocess(pvalue_files, chroms, genome='mm10', binsize=20000, filt_fdr=0.05, topDI=0.2, fo=None)`


Input:
    pvalue_file: pvalue file path
        meanA, meanB, diff, pvalue
    filt_fdr: filter pixels with fdr > filt_fdr
        if None, no filtering
    topDI: mark abs diff > topDI as topDI
        if topDI is None, don't mark
    fo: output parquet file
Output:
    df: n * (6+3) dask dataframe
        (chrom1, start1, end1, chrom2, start2, end2, meanA, meanB, diff, pvalue, qvalue)


**Source:** Line 530 in [DI.py](DI.py#L530)


---


## Ps

**File:** `Ps.py`

**Public functions:** 5


### `get_full_chrom_arm_view(clr, genome)`


Get chromosome arm regions used in clr.
Input:
    clr: cooler object
    genome: name of clr's ref genome
Output:
    dataframe


**Source:** Line 16 in [Ps.py](Ps.py#L16)


---


### `plot_ps_curve(cvd_smooth_agg, ps_col='balanced.avg.smoothed.agg')`


Plot P(s) curve.
Input:
    cvd_smooth_agg: output of ps_curve
    ps_col: column name of P(s) curve
Output:
    fig: plotly figure


**Source:** Line 82 in [Ps.py](Ps.py#L82)


---


### `plot_ps_curve_all_chromosomes(cvd_smooth_agg, ps_col='balanced.avg.smoothed', highlight_chroms=None, subset=None)`


No documentation available.


**Source:** Line 137 in [Ps.py](Ps.py#L137)


---


### `plot_ps_curves(cvd_smooth_aggs, batches, subset=None, color_map=None)`


No documentation available.


**Source:** Line 113 in [Ps.py](Ps.py#L113)


---


### `ps_curve(coolp, view_df=None, all_region=False, nproc=4, clr_weight_name='weight', drop2diags=False)`


Get P(s) curve of a balanced clr.
cooltools version 0.5.1 or higher is required.
Input:
    coolp: cooler file path
    view_df: dataframe with columns ['chrom', 'start', 'end', 'name']
    all_region: if True, return cis-exp for all regions
        else, return aggregated Ps and derivative
    nproc: number of processes to use
    clr_weight_name: weight column name
    **kwargs: other parameters for cooltools.expected_cis
Output:
    if all_region:
        cvd_smooth_agg( ps curve ), dataframe
    else:
        cvd_merged( ps curve ), der( ps curve derivative ); list of dataframe


**Source:** Line 31 in [Ps.py](Ps.py#L31)


---


## TAD.IS

**File:** `TAD/IS.py`

**Public functions:** 2


### `IS2blocks(borders)`


No documentation available.


**Source:** Line 169 in [TAD/IS.py](TAD/IS.py#L169)


---


### `find_peak_prominence(arr, max_dist=None)`


Find the local maxima of an array and their prominence.
The prominence of a peak is defined as the maximal difference between the
height of the peak and the lowest point in the range until a higher peak.

Parameters
----------
arr : array_like
max_dist : int
    If specified, the distance to the adjacent higher peaks is limited
    by `max_dist`.

Returns
-------
loc_max_poss : numpy.array
    The positions of local maxima of a given array.

proms : numpy.array
    The prominence of the detected maxima.


**Source:** Line 3 in [TAD/IS.py](TAD/IS.py#L3)


---


## TAD.tad

**File:** `TAD/tad.py`

**Public functions:** 8


### `IS(R, length)`


No documentation available.


**Source:** Line 55 in [TAD/tad.py](TAD/tad.py#L55)


---


### `bestco(F, reso, length, delta, min_val, max_val)`


No documentation available.


**Source:** Line 70 in [TAD/tad.py](TAD/tad.py#L70)


---


### `corate(A, n, time)`


No documentation available.


**Source:** Line 22 in [TAD/tad.py](TAD/tad.py#L22)


---


### `part_zero(F: np.ndarray, core: int, window: float, reso: int, delta, length, min_val, max_val)`


No documentation available.


**Source:** Line 151 in [TAD/tad.py](TAD/tad.py#L151)


---


### `silhou(R, pos)`


No documentation available.


**Source:** Line 44 in [TAD/tad.py](TAD/tad.py#L44)


---


### `tad(F, core: int=1, reso: int=40, min_val: int=600, max_val: int=1000, window_size_raw: int=8000, delta_raw: int=100)`


Input:
    F: input contact matrix : np.ndarray
    reso: resolution of the matrix (kbp)
    min, max: tad mean size (kbp)
    split: window size (kbp) 
    core: treads used


**Source:** Line 215 in [TAD/tad.py](TAD/tad.py#L215)


---


### `task(i)`


No documentation available.


**Source:** Line 167 in [TAD/tad.py](TAD/tad.py#L167)


---


### `zero(R, t, length, delta)`


No documentation available.


**Source:** Line 101 in [TAD/tad.py](TAD/tad.py#L101)


---


## ana

**File:** `ana.py`

**Public functions:** 12


### `commit(self, description=None)`


Create a new commit snapshot of the current data state.

Parameters
----------
description : str, optional
    Human-readable description of the commit. If None, uses commit ID.
    Must be unique and cannot be purely numeric.
    
Returns
-------
int
    Commit ID of the newly created commit.
    
Raises
------
ValueError
    If description already exists or is purely numeric.
    
Notes
-----
Each commit creates physical copies of data files and schema files.
Uses circular file naming to prevent infinite filename number growth.


**Source:** Line 640 in [ana.py](ana.py#L640)


---


### `convert_to_nullable_bool(val)`


No documentation available.


**Source:** Line 513 in [ana.py](ana.py#L513)


---


### `create_empty_json_file(path)`


Create an empty JSON file with empty dictionary content.


**Source:** Line 199 in [ana.py](ana.py#L199)


---


### `data(self)`


Retrieve the current structured data as a pandas DataFrame with correct dtypes.

Returns
-------
pandas.DataFrame
    Current data state with restored dtypes. Returns empty DataFrame if no data exists.
    
Notes
-----
Uses schema file to restore all pandas dtypes including categorical, string,
nullable integers, etc. Falls back to standard CSV reading if schema is missing.


**Source:** Line 605 in [ana.py](ana.py#L605)


---


### `data_path(self)`


Path to the CSV file storing structured data.


**Source:** Line 226 in [ana.py](ana.py#L226)


---


### `get_schema(self)`


Get the current data schema information.

Returns
-------
dict
    Current schema information including column dtypes.


**Source:** Line 933 in [ana.py](ana.py#L933)


---


### `list_commits(self)`


Display all available commits in chronological order.


**Source:** Line 922 in [ana.py](ana.py#L922)


---


### `obj(self)`


Retrieve the current object storage data as a dictionary.

Returns
-------
dict
    Current object storage state. Returns empty dict if no objects exist.


**Source:** Line 626 in [ana.py](ana.py#L626)


---


### `obj_path(self)`


Path to the JSON file storing unstructured object data.


**Source:** Line 236 in [ana.py](ana.py#L236)


---


### `revert(self, commit_repr)`


Revert the system state to a previous commit version.

Parameters
----------
commit_repr : str or int
    Commit identifier.
    
Returns
-------
int
    Commit ID of the new commit created after revert operation.


**Source:** Line 863 in [ana.py](ana.py#L863)


---


### `schema_path(self)`


Path to the JSON schema file storing dtype information.


**Source:** Line 231 in [ana.py](ana.py#L231)


---


### `update(self, new_data, key=None, description=None, force=False)`


Update data storage and automatically commit changes.

Parameters
----------
new_data : pandas.DataFrame, pandas.Series, dict, or list
    New data to add to the database.
key : str, optional
    Specifies the column name (for dict/Series) or object key (for list).
description : str, optional
    Custom description for the automatic commit.
force : bool, optional, default False
    If True, overwrite existing data for the given key.


**Source:** Line 741 in [ana.py](ana.py#L741)


---


## atac

**File:** `atac.py`

**Public functions:** 6


### `calculate_coverage(bam_path: str, tss_regions, read_length: int, sampling=None)`


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


**Source:** Line 175 in [atac.py](atac.py#L175)


---


### `calculate_fragment_lengths(bam_path: str, sampling: Optional[Union[float, int]]=None)`

**Returns:** `np.ndarray`


Calculate fragment lengths from a BAM file with optional sampling.

Parameters:
    bam_path (str): Path to the BAM file.
    sampling (Optional[Union[float, int]]): 
        If float, sample this proportion of reads (0 < sampling <= 1).
        If int, sample this number of reads.
        If None, no sampling (default).

Returns:
    np.ndarray: Array of fragment lengths.


**Source:** Line 11 in [atac.py](atac.py#L11)


---


### `get_tss_regions(tss_bed_path=None, chrom_size_path=None, window_size: int=2000)`


Expands TSS regions by a specified window size (default: 2000 bp upstream and downstream).

Input
   tss_bed_path (str): Path to the TSS BED file (1-based coordinates).
   chrom_size_path (str): Path to the chromosome sizes file (2-column format).
   window_size (int): Number of bases to extend on each side of the TSS.

Output
   BedTool: Expanded TSS regions as a pybedtools BedTool object.


**Source:** Line 152 in [atac.py](atac.py#L152)


---


### `normalize_coverage(coverage_array: np.ndarray, read_length: int, greenleaf_norm: bool=True)`

**Returns:** `np.ndarray`


Normalizes coverage using Greenleaf-style or per million mapped reads normalization.

Input:
   coverage_array (np.ndarray): 2D coverage array (shape: [n_TSS, bins]).
   read_length (int): Read length (used for centering adjustment).
   greenleaf_norm (bool): Whether to apply Greenleaf-style normalization (default: True).

Output:
   np.ndarray: Normalized coverage array.


**Source:** Line 233 in [atac.py](atac.py#L233)


---


### `plot_fragment_length_distribution(fragment_lengths: np.ndarray)`

**Returns:** `go.Figure`


Plot fragment length distribution using Plotly.

Input:
    fragment_lengths (np.ndarray): Array of fragment lengths.
    **kwargs: Additional keyword arguments for histogram plotting.

Output:
    go.Figure: Plotly figure object.


**Source:** Line 84 in [atac.py](atac.py#L84)


---


### `process_tss_enrichment(bam_path: str, tss_bed_path: str, chrom_size_path: str, read_length: int, window_size: int=2000, greenleaf_norm: bool=True)`

**Returns:** `Tuple[np.ndarray, float]`


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


**Source:** Line 255 in [atac.py](atac.py#L255)


---


## basic

**File:** `basic.py`

**Public functions:** 5


### `add_contact_describe(adata, inplace=True)`


Get cell's contact decay statistics, defined in Nagano2017.
Input:
    adata: AnnData object, must have uns["cdps"]
Output:
    adata will extra obs cols


**Source:** Line 7 in [basic.py](basic.py#L7)


---


### `add_pmUMIs(adata, inplace=True)`


Add cell's g1 UMI ratio.
Input:
    adata: AnnData object, must have ["g1_UMIs","g2_UMIs"]
Output:
    adata will extra obs cols


**Source:** Line 25 in [basic.py](basic.py#L25)


---


### `chrom_cov(pairs: pd.DataFrame)`

**Returns:** `pd.DataFrame`


Calculate the contact coverage of each chromosome.
Input:
    pairs: DataFrame, must have ["chr1","chr2"] in columns
Output:
    DataFrame with index as chromosome and columns as:
        intra: intra-chromosomal contacts
        tot: total contacts
Note:
    total contacts equals to:
        (res["tot"].sum() + res["intra"].sum()) / 2


**Source:** Line 39 in [basic.py](basic.py#L39)


---


### `pairs_coverage(pairs_fp, genome: str=None, binsize: int=1000000, flavor: str='hickit', sub: str=None)`

**Returns:** `pd.DataFrame`


Calculate coverage of pairs in bins.
Input:
    pairs_fp: str, path to pairs file
    genome: str, genome name
    binsize: int, size of bins
    flavor: str, flavor of bins
    sub: only calculate coverage of a subset of genome
Output:
    coverage: coverage of pairs in each bin


**Source:** Line 82 in [basic.py](basic.py#L82)


---


### `pairs_coverages(pairs_fps, sample_names, n_jobs=8, genome: str=None, binsize: int=1000000, flavor: str='hickit', sub=None)`

**Returns:** `pd.DataFrame`


Multiple pairs_coverage
Input:
    pairs_fps: list of str, paths to pairs files
    sample_names: list of str, sample names
    genome, binsize, flavor: see pairs_coverage
Output:
    n_bins * n_samples DataFrame


**Source:** Line 131 in [basic.py](basic.py#L131)


---


## binnify

**File:** `binnify.py`

**Public functions:** 1


### `bed_center(bed_df)`


Get center of each bin
Input:
    bed_df: DataFrame, bed-like with columns 'chrom', 'start', 'end'
Output:
    DataFrame with additional column 'center'


**Source:** Line 11 in [binnify.py](binnify.py#L11)


---


## calculate

**File:** `calculate.py`

**Public functions:** 4


### `decorator(func)`


No documentation available.


**Source:** Line 43 in [calculate.py](calculate.py#L43)


---


### `mt(module_name: str)`


A decorator factory to create parallelized function decorators.

This factory returns a decorator that parallelizes the function it decorates,
providing features like file skipping, path templating, progress bars, and result handling.

Args:
    module_name (str): The name of the module where the function to be decorated is located.

Returns:
    A decorator function that can be applied to a function.


**Source:** Line 28 in [calculate.py](calculate.py#L28)


---


### `task_wrapper(module_name, func_name)`


Wrapper function to import and execute a function from a module.

Parameters:
    module_name (str): Name of the module to import from
    func_name (str): Name of the function to execute
    *args: Positional arguments to pass to the function
    **kwargs: Keyword arguments to pass to the function
    
Returns:
    Result of the executed function


**Source:** Line 10 in [calculate.py](calculate.py#L10)


---


### `wrapper(input_data: pd.DataFrame, force=False, nproc=None, concat=False, new_cols=None)`


Parallelized function wrapper with input/output pattern support.

The function uses a pattern-based system for input and output file handling:
- Input sources: input_cols, input_pattern, input_pattern1, input_pattern2, etc.
- Output sources: output_cols, output_pattern, output_pattern1, output_pattern2, etc.

Parameters are processed in order: input_cols first, then input_patterns.
Patterns are formatted using {sample_name} placeholder which is replaced with
the DataFrame index value for each row.

Parameters:
    input_data (pd.DataFrame): Input DataFrame containing data to process
    force (bool): Whether to overwrite existing output files (default: False)
    nproc (int): Number of parallel processes (default: all CPUs)
    concat (bool): If True, concatenate Series results into a DataFrame
    new_cols (list): Column names for pattern-based output files (non-concat mode only)
    *args: Additional positional arguments passed to the decorated function
    **kwargs: Additional keyword arguments including input/output patterns
    
Pattern Parameters (passed as kwargs):
    input_cols: Column name(s) containing input paths (single or list)
    output_cols: Column name(s) containing output paths (single or list)
    input_pattern: Template for input file paths
    input_pattern1, input_pattern2, ...: Additional input patterns
    output_pattern: Template for output file paths  
    output_pattern1, output_pattern2, ...: Additional output patterns
    
Returns:
    If concat=True and results are pandas Series: DataFrame with concatenated results
    Otherwise: DataFrame with pattern-based output paths and execution status


**Source:** Line 45 in [calculate.py](calculate.py#L45)


---


## cli.download

**File:** `cli/download.py`

**Public functions:** 2


### `add_arguments(subparser)`


No documentation available.


**Source:** Line 38 in [cli/download.py](cli/download.py#L38)


---


### `run(args)`


No documentation available.


**Source:** Line 71 in [cli/download.py](cli/download.py#L71)


---


## cli.render

**File:** `cli/render.py`

**Public functions:** 3


### `add_arguments(subparser)`


No documentation available.


**Source:** Line 102 in [cli/render.py](cli/render.py#L102)


---


### `render_task(args)`


Render a single structure.


**Source:** Line 21 in [cli/render.py](cli/render.py#L21)


---


### `run(args)`


No documentation available.


**Source:** Line 161 in [cli/render.py](cli/render.py#L161)


---


## cli.version

**File:** `cli/version.py`

**Public functions:** 2


### `add_arguments(subparser)`


No documentation available.


**Source:** Line 3 in [cli/version.py](cli/version.py#L3)


---


### `run()`


No documentation available.


**Source:** Line 5 in [cli/version.py](cli/version.py#L5)


---


## clustering

**File:** `clustering.py`

**Public functions:** 1


### `do_kmeans(adata: ad.AnnData, layer: str)`

**Returns:** `ad.AnnData`


Input:
    adata
    layer: which layer to operate kmeans on
Return:
    adata with new obs col.


**Source:** Line 20 in [clustering.py](clustering.py#L20)


---


## compare_pixels

**File:** `compare_pixels.py`

**Public functions:** 5


### `compare_pixels(this_coolp, ref_coolp, threshold=5000000, clr_weight_name='weight')`

**Returns:** `pd.DataFrame`


Extract and align pixels from two coolp files. 
Input:
    this_coolp: str, path to first coolp file
    ref_coolp: str, path to second coolp file
    threshold: int, band distance threshold
        only pixels within this distance will be compared
    clr_weight_name: str, column name for clr weight
        None for raw contact values
Output:
    concated: dataframe of merged pixels
        chrom1, start1, start2, balanced_this, balanced_ref, distance


**Source:** Line 25 in [compare_pixels.py](compare_pixels.py#L25)


---


### `histogram2d(data, x1, x2, bins=100, range=None)`


Calculate 2D histogram for two columns in a DataFrame.
Density plot version.
Input:
    data: pd.DataFrame, data
    x1: str, column name for x axis
    x2: str, column name for y axis
    bins: int, number of bins
    range: list of list, range for x and y axis
Output:
    histo: np.array, 2D histogram
    x_edges: np.array, x axis edges
    y_edges: np.array, y axis edges


**Source:** Line 140 in [compare_pixels.py](compare_pixels.py#L140)


---


### `pick_band_pixels(coolp, threshold=5000000)`

**Returns:** `pd.DataFrame`


Pick (intra) pixels within a certain band distance.
Input:
    coolp: str, path to coolp file
    threshold: int, band distance threshold
    **kwargs: additional arguments for cool2pixels
Output:
    pixels: pixels within the threshold
        see cool2pixels for columns, with additional band_dist


**Source:** Line 7 in [compare_pixels.py](compare_pixels.py#L7)


---


### `plot_compare_pixels(this_coolp, ref_coolp, sample_titles=['THIS', 'REF'], ignore_diags=0, clr_weight_name='weight', lognorm=False, sample_ratio=None, debug=False)`


Compare pixels between two coolp files.
Scatter plot version.
Input:
    this_coolp: str, path to first coolp file
    ref_coolp: str, path to second coolp file
    sample_titles: list of str, titles for each coolp
    ignore_diags: bool or int, ignore diagonal pixels
        if True, ignore 1 bin, if int, ignore n bins
    lognorm: bool, log10 transform contact values
Output:
    fig: plotly figure


**Source:** Line 53 in [compare_pixels.py](compare_pixels.py#L53)


---


### `plot_compare_pixels_density(this_coolp, ref_coolp, qrange=None, sample_titles=['THIS', 'REF'], bins=100, ignore_diags=0, clr_weight_name='weight', lognorm=False, vlognorm=False, retdata=False, vrange=None, label_x=0.01, label_y=0.015, debug=False)`


Compare pixels between two coolp files and plot density plot.
Input:
    this_coolp: str, path to first coolp file
    ref_coolp: str, path to second coolp file
    qrange: list of float, quantiles for x and y axis
    sample_titles: list of str, titles for each coolp
    ignore_diags: bool or int, ignore diagonal pixels
        if True, ignore 1 bin, if int, ignore n bins
    vlognorm: bool, log10 transform contact values
    clr_weight_name: col name of weight
    lognorm: bool, log10 transform of values for histogram2d
    retdata: bool, return data for further processing,
    vrange: list of list, range for x and y axis for histogram2d,
        if both qrange and vrange are set, vrange will be used
    label_x: float, relative x position for pearson correlation annotation
    label_y: float, relative y position for pearson correlation annotation
Output:
    if retdata:
        data: pd.DataFrame, data used for plotting
    else:
        fig: plotly figure


**Source:** Line 170 in [compare_pixels.py](compare_pixels.py#L170)


---


## compartment

**File:** `compartment.py`

**Public functions:** 8


### `AB_block_ends(orig_data: pd.DataFrame, min_winSize=3, binsize=1000000)`


Generate AB bed-like reference from cooltools vecs dataframe.
Input:
    orig_data: chrom, start, end, E1
Output:
    chrom, start, end, AB, mean_E1; one row per block


**Source:** Line 217 in [compartment.py](compartment.py#L217)


---


### `call_compartment(df: pd.DataFrame, norm: bool=True)`

**Returns:** `tuple`


No documentation available.


**Source:** Line 208 in [compartment.py](compartment.py#L208)


---


### `compartments(M, normalize=True, matrixonly=False, nPCs=2)`


A/B compartment analysis

Perform a PCA-based A/B compartment analysis on a normalized, single
chromosome contact map. The results are two vectors whose values (negative
or positive) should presumably correlate with the presence of 'active'
vs. 'inert' chromatin.

Parameters
----------
M : array_like
    The input, normalized contact map. Must be a single chromosome.
normalize : bool
    Whether to normalize the matrix beforehand.
matrixonly : bool
    Whether to return only the matrix used for the PCA (useful in compartment plotting).
nPCs : int
    Number of principal components to compute.

Returns
-------
pandas.DataFrame
    A DataFrame with the first two principal components as columns and
    the genomic coordinates as index.


**Source:** Line 146 in [compartment.py](compartment.py#L146)


---


### `distance_law_from_mat(matrix, indices=None, log_bins=True, base=1.1)`


Compute distance law as a function of the genomic coordinate aka P(s).
Bin length increases exponentially with distance if log_bins is True. Works
on dense and sparse matrices. Less precise than the one from the pairs.
Parameters
----------
matrix : numpy.array or scipy.sparse.coo_matrix
    Hi-C contact map of the chromosome on which the distance law is
    calculated.
indices : None or numpy array
    List of indices on which to compute the distance law. For example
    compartments or expressed genes.
log_bins : bool
    Whether the distance law should be computed on exponentially larger
    bins.
Returns
-------
numpy array of floats :
    The start index of each bin.
numpy array of floats :
    The distance law computed per bin on the diagonal


**Source:** Line 3 in [compartment.py](compartment.py#L3)


---


### `extrude_full_zero(mat: pd.DataFrame)`

**Returns:** `pd.DataFrame`


No documentation available.


**Source:** Line 57 in [compartment.py](compartment.py#L57)


---


### `normalize_dense(M, norm='SCN', order=1, iterations=40)`


Apply one of the many normalization types to input dense
matrix. Will also apply any callable norms such as a user-made
or a lambda function.
NOTE: Legacy function for dense maps

Parameters
----------
M : 2D numpy array of floats
norm : str
    Normalization procedure to use. Can be one of "SCN",
    "mirnylib", "frag" or "global". Can also be a user-
    defined function.
order : int
    Defines the type of vector norm to use. See numpy.linalg.norm
    for details.
iterations : int
    Iterations parameter when using an iterative normalization
    procedure.

Returns
-------
2D numpy array of floats :
    Normalized dense matrix.


**Source:** Line 72 in [compartment.py](compartment.py#L72)


---


### `sc_compartment_strength(cool_f, outfile=None, eigen_track=None)`


Calculate compartment strength from a single-cell eigenvector track.
Input:
    cool_f: cooler file path
    eigen_track: eigenvector track file path, e.g. CpG/GC count track.
    outfile: output file path, if None, return a list of compartment strength values.


**Source:** Line 255 in [compartment.py](compartment.py#L255)


---


### `upper_to_symmetry(mat: pd.DataFrame)`

**Returns:** `pd.DataFrame`


No documentation available.


**Source:** Line 64 in [compartment.py](compartment.py#L64)


---


## coolstuff

**File:** `coolstuff.py`

**Public functions:** 32


### `cli_IS(coolp, output, windowsizes, balanced=True, append_raw_scores=True, threads=8, conda_env=None, cwd=None, force=False, verbose=0)`


No documentation available.


**Source:** Line 1120 in [coolstuff.py](coolstuff.py#L1120)


---


### `cli_balance(coolp, threads=8, force=False, name='weight', cis_only=False, trans_only=False, min_nnz=None, min_count=None, mad_max=None, ignore_diags=None, conda_env=None, cwd=None, verbose=0)`


Balance a cooler matrix.
Input:
    coolp: cooler file path
    threads: number of threads
    force: whether to overwrite existing balanced matrix
    name: name of the balanced matrix
    cis_only: balance only cis contacts
    trans_only: balance only trans contacts
    min_nnz: minimum number of non-zero contacts
    min_count: minimum number of contacts
    mad_max: mad_max filter, lower values are more stringent
        NOTE: this dominates all other filtering options
    conda_env: conda environment to run in
    cwd: working directory
    verbose: verbosity level, default 0. When > 0, print message when skipping execution.
Output:
    True if successful, False otherwise


**Source:** Line 971 in [coolstuff.py](coolstuff.py#L971)


---


### `cli_compartment(coolp, phasing_track, outprefix, view=None, conda_env=None, cwd=None, force=False, verbose=0)`


No documentation available.


**Source:** Line 1053 in [coolstuff.py](coolstuff.py#L1053)


---


### `cli_downsample(coolp, output, count=100000000.0, cis_count=None, fraction=None, threads=8, force=False, conda_env=None, cwd=None, verbose=0)`


No documentation available.


**Source:** Line 951 in [coolstuff.py](coolstuff.py#L951)


---


### `cli_expected(coolp, output, balanced=False, view=None, ignore_diags=1, conda_env=None, cwd=None, threads=8, force=False, verbose=0)`


No documentation available.


**Source:** Line 1180 in [coolstuff.py](coolstuff.py#L1180)


---


### `cli_mergecool(incools, outcool, force=False, conda_env=None, skip_blank=True, cwd=None, batch_size=1000, verbose=0)`


Merge cool files with same indices to get a consensus heatmap.
Handles large numbers of input files by batching them.

Input:
    incools: list of .cool file paths
    outcool: output .cool file path
    conda_env: (optional) conda environment to run in
    cwd: (optional) working directory
    batch_size: number of files to merge in each batch
    verbose: verbosity level, default 0. When > 0, print message when skipping execution.


**Source:** Line 894 in [coolstuff.py](coolstuff.py#L894)


---


### `cli_pairs2cool(filei, fileo, sizef, binsize, force=False, conda_frontend='', verbose=0)`


Generate .cool file from 4DN .pairs file
Input:
    filei: input pairs file
    fileo: output .cool file
    sizef: chrom size file, 2col(chrom:str, size:int) tsv
    binsize: binsize as you wish
    conda_frontend: use conda to run cmd in conda env
        eg. "conda run -n cooler"
    verbose: verbosity level, default 0. When > 0, print message when skipping execution.


**Source:** Line 822 in [coolstuff.py](coolstuff.py#L822)


---


### `cli_pileup(coolp, feature, output, format='BED', view=None, expected=None, flank=100000, store_snips=False, conda_env=None, cwd=None, threads=8)`


No documentation available.


**Source:** Line 1142 in [coolstuff.py](coolstuff.py#L1142)


---


### `cli_saddle(coolp, eigv, expected, outprefix, view, bin_method='quantile', qrange=[0.02, 0.98], vrange=[-1, 1], nbins=50, conda_env=None, cwd=None, force=False, debug=False, verbose=0)`


Run cooltools saddle.
Input:
    coolp: cooler file path
    eigv: eigenvector file path
    expected: expected file path
    outprefix: output prefix
    view: genome view file path
    bin_method: binning method, either 'quantile' or 'value'.
        If 'quantile', qrange should be provided.
        If 'value', vrange should be provided.
    qrange: quantile range
    vrange: value range
    nbins: number of bins
    conda_env: conda environment to run in
    cwd: working directory
    force: whether to overwrite existing files
    verbose: verbosity level, default 0. When > 0, print message when skipping execution.


**Source:** Line 1072 in [coolstuff.py](coolstuff.py#L1072)


---


### `cli_zoomify(coolp, output, resolutions=[20000, 40000, 100000, 500000, 1000000], force=False, threads=8, conda_env=None, cwd=None, verbose=0)`


No documentation available.


**Source:** Line 1039 in [coolstuff.py](coolstuff.py#L1039)


---


### `cool2mat(cool, region: Union[str, List[str], slice, List[slice]], balance: bool=False, min_nz=None, bad_bins=None)`

**Returns:** `pd.DataFrame`


Fetch matrix from cooler file with proper index and columns.
Input:
    cool: cooler file path
    region: genome region to fetch.
        "chr1" or "chr1:1000000-2000000":
            ucsc styl, get intra-chrom matrix
        ["chr1:1,000,000-2,000,000", "chr2"]:
            list of region string to get inter-chrom matrix
        slice(0,-1):
            slicer syntax, treat the whole genome as a big matrix,get intra-chrom matrix
            Note: can only cross chrom boundaries when using slicer syntax
            For example slice(0,-1) for whole genome.
        [slice(0,1000000), slice(1000000,2000000)]:
            inter-chrom matrix using slicer syntax, can get cross-chrom matrix
        [[(r1_chrom1,start1),(r1_chrom2,end1)], [(r2_chrom1,start2),(r2_chrom2,end2)]]:
            list of 2 standard region to get cross-chrom matrix in human readable format
    balance: balance or not, bool
    min_nz: minimum number of non-zero pixels to keep, bin with less pixels will be set to nan
        if None, no filtering
    bad_bins: list of bad bins to be set to nan, a dataframe with columns ["chrom","start"]
Output:
    pd.DataFrame
    region1 -> numpy 1st axis; numpy/pandas row; pandas index;
    region2 -> numpy 2nd axis; numpy/pandas column; pandas columns;
TODO:
    min_nz and bad_bins support for separate index(region1) and columns(region2) filtering for inter-chrom matrix


**Source:** Line 132 in [coolstuff.py](coolstuff.py#L132)


---


### `cool2pixels(coolp: str, clr_weight_name: bool=None)`

**Returns:** `pd.DataFrame`


Fetch all pixels from a cooler file.
Input:
    coolp: path to cooler file
    clr_weight_name: name of weight column
        if None, use raw count
Output:
    pixels: DataFrame of pixels, [bin1_id, bin2_id, count, balanced]


**Source:** Line 113 in [coolstuff.py](coolstuff.py#L113)


---


### `cool_uri_to_path(cool_uri)`


Extract file path from a cool URI.

Examples:
    - "xx.cool" -> "xx.cool"
    - "xx.mcool::/resolutions/20000" -> "xx.mcool"


**Source:** Line 884 in [coolstuff.py](coolstuff.py#L884)


---


### `cools2scool(cools, scool_path)`


Generate scool from cooler files.
TODO:
    don't readin all cools together.
Input:
    cools: dict of name-cooler file, using first file's bins as default bins
    scool_path: outputfile name


**Source:** Line 616 in [coolstuff.py](coolstuff.py#L616)


---


### `gam2cool(gam_filepat, bin_file, outfile, force=False, intra_only=True)`


Transform gam result to cool format.
Use intra_only when binsize is small.
Input:
    gam_filepat: path pattern to gam file, use * to represent any string,
        typically chrom names
    bin_file: path to bin file created by bedtools in the pipeline
    outfile: path to output cool file
    force: whether to overwrite existing file
    intra_only: only keep intra-chrom
    **kwargs: other args for cooler.create_cooler
Output:
    outfile


**Source:** Line 547 in [coolstuff.py](coolstuff.py#L547)


---


### `gen_bins(genome, binsize)`


Generate "bins" for cooler api.
Input:
    genome: RefGenome obj or Chromosomes obj or chromsizes dataframe with columns ["length"]
    binsize: bin size; int


**Source:** Line 243 in [coolstuff.py](coolstuff.py#L243)


---


### `gen_config(sample_table, cfg, assembly='mm10', bwa='/share/home/ychi/data/genome/GRCm38/bwa_index/mm10.*', chrom_size='/share/home/ychi/data/genome/GRCm38/mm10.len.forCool.tsv', template='/shareb/ychi/ana/distiller-nf/project.yml')`


Generate distiller.nf project.yaml file from bubble_flow sample_table.csv.
Input:
    sample_table: any csv, index as sample_name, using sample-name and task_dirp to infer DNA reads files
    cfg: output project.yml file path
    assembly: genome assembly name, default mm10
    template: template project.yml file path


**Source:** Line 722 in [coolstuff.py](coolstuff.py#L722)


---


### `hic_pileup(scool, grouping, cache_pattern='{}.pileup.cool', mergebuf=1000000.0)`


Generate pileup cool file for each cluster.
To simply pileup all cells, use `grouping = pd.Series(1, index=cells)`
Input:
    scool: cooler's scool file path
    grouping: dict, key is cluster name, value is list of cell names
    cache_pattern: path to store resulting .cool file; str with {}
Output:
    write cool according to cache_pattern


**Source:** Line 636 in [coolstuff.py](coolstuff.py#L636)


---


### `iter_pixels(clr, chunksize=1000000.0, join=True)`


No documentation available.


**Source:** Line 108 in [coolstuff.py](coolstuff.py#L108)


---


### `mt_pairs2cool(sample_table, outdir, pairs_col='pairs_c123', suffix='.cool', nproc=16)`


TODO: Will still print first stdout and stderr message now. Fix this.


**Source:** Line 859 in [coolstuff.py](coolstuff.py#L859)


---


### `pairs2cool(pairs_path, coolpath, sizef, binsize)`


Generate .cool file from 4DN .pairs file


**Source:** Line 303 in [coolstuff.py](coolstuff.py#L303)


---


### `pairs2cool_worker(sample, sample_table, outdir, pairs_col, suffix, kwargs)`


This function is used to run pairs2cool in parallel.
TODO: hide output


**Source:** Line 843 in [coolstuff.py](coolstuff.py#L843)


---


### `pairs2scool(pairs_paths, coolpath, sizef, binsize)`


Generate .scool file from multiple 4DN .pairs file
Input:
    pairs_path: pairs file path; dict
        {cell_name: pairs_path}
    coolpath: output .scool file path; str
    sizef: chrom size file path, 2col(chrom:str, size:int) tsv without header
    binsize: bin size; int


**Source:** Line 592 in [coolstuff.py](coolstuff.py#L592)


---


### `parse_gam(gam_file, pixel=False)`


Input:
    gam_file: path to gam file
        temporarily, only support .npz file
    pixel: a cooler's joint pixel format
Output:
    mat: pandas.DataFrame
        if not pixel, a DataFrame with 2-level-multiindex (chrom, start) index and columns


**Source:** Line 475 in [coolstuff.py](coolstuff.py#L475)


---


### `peek_gam(gam_file)`


Peek into gam file to get some basic info.
Input:
    gam_file: path to gam file
Output:
    dict of basic info
        chrom1: most common chrom1
        chrom2: most common chrom2


**Source:** Line 457 in [coolstuff.py](coolstuff.py#L457)


---


### `pixels2cool(pixels: pd.DataFrame, bins: pd.DataFrame, outfile, columns: list=['count'], dtypes: dict={'count': int})`


Dump pixels to cool file.
Input:
    pixels: pixels dataframe or iterator of pixels dataframes,
        [chrom1, start1, (end1), chrom2, start2, (end2), vcol1, vcol2, ...]
    bins: [chrom, start, end]
    outfile: output cool file path
    columns: value columns to be stored
    dtypes: data type for each value column
    **kwargs: other args for cooler.create_cooler
Output:
    outfile
TODO: deal when bin1_id or bin2_id in pixels


**Source:** Line 313 in [coolstuff.py](coolstuff.py#L313)


---


### `pos2id(clr, chrom, pos)`


Transform a genomic position to a bin ID in the cooler matrix.

Parameters
----------
clr : cooler.Cooler or str
    A cooler object or path to cooler file.
chrom : str
    Chromosome name.
pos : int
    Genomic position.

Returns
-------
int
    The bin ID corresponding to the right boundary of the genomic position.


**Source:** Line 45 in [coolstuff.py](coolstuff.py#L45)


---


### `region2slice(clr, region)`


Transform a genomic region to a slice of a cooler matrix.
Region can cross chromosomes.

Parameters
----------
clr : cooler.Cooler or str
    A cooler object.
region : see hic_basic.genome.Region

Returns
-------
slice
    A slice object that can be used to index the cooler matrix.

Examples
--------
>>> import cooler
>>> clr = cooler.Cooler("example.mcool::/resolutions/10000")
>>> region = "chr1:0-1000000" # [('chr1', 0), ('chr1', 100000000)], or [('chr1', 200000000), ('chr2', 50749379)]
>>> islice = region2slice(clr, region)


**Source:** Line 67 in [coolstuff.py](coolstuff.py#L67)


---


### `reset_cool_bins(coolpin, coolpout, new_bins=None, genome='GRCh38', chunksize=1000000.0)`


Reset bins of cool file. Can also use this to filter out some bins.
Input:
    coolpin: input cool file path
    coolpout: output cool file path
    bins: new bins, if None, generate from genome
    genome: genome name, default "GRCh38"
    chunksize: chunksize for reading pixels


**Source:** Line 655 in [coolstuff.py](coolstuff.py#L655)


---


### `schicluster2cool(filesp, fo, genome, binsize, ordered=True, ensure_sorted=True)`


Dump schicluster imputed matrix to .cool file.
Input:
    filesp: All hdf5 file paths for 1 cell, chrom as index; pd.Series
    fo: output cool uri
    genome: eg. "mm10"
    binsize: binsize of original schicluster matrix
    **args: other args for cooler.create_cooler
Output:
    fo
    write to file and return file path
Example:
```
    from pathlib import Path
    import pandas as pd
    from hic_basic.data import mm10
    from hic_basic.coolstuff import schicluster2cool
    
    
    sample = "2021110112"
    ddir = "/shareb/ychi/repo/embryo_integrate/schicluster/imputed_matrix/100000/"

    chroms = [chrom for chrom in mm10.chromosomes.data.index if chrom not in ["chrX","chrY"]]
    filesp = pd.Series(
        data = [
            Path(ddir)/ "{chrom}".format(chrom=chrom) / "{sample}_{chrom}_pad1_std1_rp0.5_sqrtvc.hdf5".format(sample=sample, chrom=chrom) 
            for chrom in chroms
        ],
        index = chroms,
        name = sample)
    
    schicluster2cool(filesp, "test.cool","mm10",100e3)
```


**Source:** Line 404 in [coolstuff.py](coolstuff.py#L404)


---


### `stat_cool(names, ddir, genome, show=True)`


Statistic technical details.
Input:
    names: list of cooler paths
    ddir: distiller-nf directory
    genome: genome name eg. "mm10"
Output:
    pd.DataFrame


**Source:** Line 763 in [coolstuff.py](coolstuff.py#L763)


---


### `str2slice(clr, region: str)`


Convert region string to slice object.
Input:
    clr: cooler file object
    region: region string, eg. "chr1:1,000,000-2,000,000"
Output:
    slice object


**Source:** Line 98 in [coolstuff.py](coolstuff.py#L98)


---


## cycle_phasing

**File:** `cycle_phasing.py`

**Public functions:** 2


### `dis_counts(pairs_fp: str, filter_str=None)`


Get cell's intra contact's distribution in Nagano2017's window
Using customized .pairs parser. Work for 11 column table only
Input:
    pairs_fp: path to .pairs file
    filter_str: expr string to filter pairs file
Output:
    counts: pd.Series of contact distance distribution


**Source:** Line 31 in [cycle_phasing.py](cycle_phasing.py#L31)


---


### `window_count(distances: pd.Series, win_num)`

**Returns:** `pd.Series`


Count distribution of distance array with Nagano2017's window.
If distances is None, return nan Series.
Input:
    distances: distance array
    win_num: number of windows
Output:
    window_count: distribution of distances in Nagano


**Source:** Line 7 in [cycle_phasing.py](cycle_phasing.py#L7)


---


## data

**File:** `data.py`

**Public functions:** 18


### `Jichang2022_embryo_gene_module()`


No documentation available.


**Source:** Line 743 in [data.py](data.py#L743)


---


### `XijinGe2017_embryo_gene_module()`


Gene clusters from `Exploratory bioinformatics investigation reveals importance of “junk” DNA in early embryo development`.
Note:
    has excel gene-name to Date problem
Output:
    dict of list.
    A: Oocyte-specific
    B: 2-4 transient
    C: oocyte-16C
    D: 2-16C transient strong
    E: 2-16 suppressed
    F: 2-16C mild
    G: 2-4-8 graduate activation
    H: Blastocyst strong activation (starts at 8c)
    I: Blastocyst specific A
    J: Blastocyst specific B


**Source:** Line 752 in [data.py](data.py#L752)


---


### `callback(data)`


No documentation available.


**Source:** Line 265 in [data.py](data.py#L265)


---


### `chromosomes(genome, order=False)`


Get chromosome lengths.
Input:
    genome: genome name or chromosome length table; string or pd.DataFrame;
        Now support hg19, hg19_dip, GRCh38, GRCh38_dip, mm10, mm10_dip
    order: if True, order chromosomes in a natural way (chr1, chr2 ...),
        if list order chromosomes by the list given;
Return:
    pandas.DataFrame


**Source:** Line 481 in [data.py](data.py#L481)


---


### `download_geo(geo_id, outdir, method='ftp', timeout=30, retries=3, block_size=8192)`


Download GEO data with enhanced error handling, retries, larger block size, passive mode for FTP, and download resumption.

Input:
    geo_id: GEO ID; string
    outdir: output directory; string
    method: download method; string; "ftp" or "https"
    timeout: timeout for downloads in seconds; int
    retries: number of retries for failed downloads; int
    block_size: initial block size for downloading large files; int


**Source:** Line 210 in [data.py](data.py#L210)


---


### `download_sra(accessions, outdir, method='ena', timeout=30, retries=3, block_size=8192)`


Download sequencing data by SRA run accession(s).

Input:
    accessions:
        - single accession string (e.g., "SRR1234567")
        - comma-separated string (e.g., "SRR1,SRR2")
        - list/tuple of accessions
        Note: reading a run-list text file is handled by the CLI; pass a list here.
    outdir: output directory
    method:
        - "ena": use ENA filereport to resolve FASTQ URLs and download via HTTPS/FTP with resume and retries
        - "sra-tools": use prefetch + fasterq-dump if SRA Toolkit is available in PATH
    timeout: request timeout in seconds
    retries: retry count for transient errors
    block_size: chunk size for streaming downloads


**Source:** Line 331 in [data.py](data.py#L331)


---


### `dupref_annote(bins, ref)`


Annotate diploid genome with haploid ref
Input:
    bins: pd.DataFrame, 2-level index
    ref: pd.DataFrame, 2-level index
Output:
    pd.DataFrame, 3dg with ref


**Source:** Line 452 in [data.py](data.py#L452)


---


### `fetch_TSS(gnames, TSS, name_col='gene_name')`


Fetch TSS regions for given gene names.
Input:
    gnames: genename list
    TSS: genome name or TSS reference table; string or pd.DataFrame
    name_col: column name of gene name in TSS table; string
Output:
    chromname, TSS id, TSS sites; pd.Series


**Source:** Line 701 in [data.py](data.py#L701)


---


### `fetch_cent_chromlen(genome)`


Get centromere regions and chromosome lengths.
Only work with mm10 for now.
TODO: generalize to other genomes
Input:
    fp (str): path to the gap file
Returns:
    dict(chrom) of list[start, end, chrom_length]


**Source:** Line 543 in [data.py](data.py#L543)


---


### `fetch_centromeres(genome)`


Get centromere regions.
Returns:
    pandas.DataFrame


**Source:** Line 521 in [data.py](data.py#L521)


---


### `get_chrom_data(genome, order=False)`


Get chromosome data from various genome representations.

This helper function provides a unified interface to extract chromosome
length data from different genome input formats.

Args:
    genome: Genome specification in one of these forms:
        - RefGenome object (e.g., mm10, GRCh38)
        - Chromosomes object from feature module
        - pandas DataFrame with chromosome lengths
        - str: genome name (e.g., "mm10") or path to chromosome sizes file
    order (bool): If True, ensure chromosomes are ordered naturally
        (chr1, chr2, ..., chr10, ..., chrX, chrY).

Returns:
    pd.DataFrame: DataFrame with chromosome names as index and at least
        a 'length' column containing chromosome lengths.

Raises:
    ValueError: If genome cannot be resolved or has invalid format.
    FileNotFoundError: If genome string refers to non-existent file.

Examples:
    >>> from hic_basic.data import mm10, get_chrom_data
    >>> # From RefGenome object
    >>> chrom_df = get_chrom_data(mm10)
    >>> # From genome name string
    >>> chrom_df = get_chrom_data("mm10")
    >>> # From file path
    >>> chrom_df = get_chrom_data("/path/to/chrom.sizes")


**Source:** Line 125 in [data.py](data.py#L125)


---


### `id2name(genelist, genome)`


Convert gene IDs to gene names.
Input:
    genelist: list of gene IDs
    genome: genome name or gene ID to gene name mapping dataframe; string or pd.DataFrame
        gene ID as index when using a dataframe
Output:
    gene names; list


**Source:** Line 643 in [data.py](data.py#L643)


---


### `mouse_GO_cell_cycle()`


Mouse cell cycle genes with term: `G2/M transition of mitotic cell cycle` and `G1/S transition of mitotic cell cycle`.
Output:
    dict of list


**Source:** Line 789 in [data.py](data.py#L789)


---


### `mouse_cell_cycle()`


Cell cycle gene sets in mouse. Used in seurat or scanpy g1/s, g2/m score.
Output:
    dict of list.


**Source:** Line 779 in [data.py](data.py#L779)


---


### `name2id(genelist, genome, pickid=False)`


Convert gene names to gene IDs.
Input:
    genelist: list of gene names
    genome: genome name or gene ID to gene name mapping dataframe (name 2 id df will have non-unique index); string or pd.DataFrame
        gene id as index when using a dataframe
    pickid: if True, pick 1 gene ID for genes with multiple IDs, else return NA for genes with multiple IDs
Output:
    gene IDs; list


**Source:** Line 664 in [data.py](data.py#L664)


---


### `register_chromosomes(self, size_file)`


Register chromosome sizes for this genome.

Args:
    size_file (str | pathlib.Path): Path to a two-column
        chromosome sizes file.


**Source:** Line 80 in [data.py](data.py#L80)


---


### `register_feature(self, feature_name, factory)`


Register a feature factory.

Args:
    feature_name (str): Attribute name for the feature.
    factory (Callable[[], Any]): Function that builds the feature.


**Source:** Line 60 in [data.py](data.py#L60)


---


### `split_links(v)`


No documentation available.


**Source:** Line 391 in [data.py](data.py#L391)


---


## embedding

**File:** `embedding.py`

**Public functions:** 5


### `NN_spectral_embedding(symdist, k_nn=7)`


Non-linear dimensionality reduction of cdps curve
Input:
    symdist: symmetric distance matrix
    k_nn: number of neighbors when build NN graph
Output:
    2nd 3rd smallest eigenvectors of the graph Laplacian
eg:
    from hic_basic.metrics import distance_pca_euclid
    from hic_basic.embedding import NN_spectral_embedding
    from hic_basic.pseudotime.TI import angle_pseudotime

    dm = distance_pca_euclid(df.values, 6)
    dm = pd.DataFrame(
        dm,
        index = sub.index,
        columns = sub.index
    )
    spec = NN_spectral_embedding(dm)
    pt = angle_pseudotime(0,spec) # 0 is the root cell


**Source:** Line 120 in [embedding.py](embedding.py#L120)


---


### `band_seg_svd(cell_list, res=100000.0, segL=10000000.0, dist=10000000.0, dim=50)`


Doing svd on cell matrix band(leg distance < dist) segments.
Input:
    cell_list: file list of cell matrix.
    res: resolution of matrix(binsize in bp).
    segL: length of segment(in bp).
    dist: distance of band, only pixels within the near-diagonal band are considered.
    dim: number of components to keep.
Return:
    list of svd results, same length as nseg.


**Source:** Line 69 in [embedding.py](embedding.py#L69)


---


### `band_svd(cell_list, res, dist=10000000, dim=50)`


Doing svd on cell matrix band(leg distance < dist).
TODO: using api compatible to cools.


**Source:** Line 35 in [embedding.py](embedding.py#L35)


---


### `pca_rep(mat, n_components, with_std=False, retmodel=False)`


Return pca representation of data.
Input:
    mat: m * n matrix
    n_components: number of components to keep
    with_std: whether doing std scaling. No consensus but 
        usually don't scale it in omic-biology.
    retmodel: whether return all related model objects.
Output:
    m * n_components matrix


**Source:** Line 14 in [embedding.py](embedding.py#L14)


---


### `plot_elbow(pca: PCA, title='TITLE')`


No documentation available.


**Source:** Line 167 in [embedding.py](embedding.py#L167)


---


## feature

**File:** `feature.py`

**Public functions:** 10


### `compile(self)`


Process standardized genomic data formats (e.g., GTF, VCF) and convert them into simplified intermediate formats like BED or CSV.

This abstract method serves two primary purposes:
1. Data Transformation: Extracts core features from complex genomic annotations and restructures them into standardized intermediate formats for downstream analysis.
2. Caching Mechanism: Generates reference files (stored in `self.fn`) that can be directly consumed by external tools like pybedtools, BEDTools, or other bioinformatics pipelines.

The implementation must be customized by subclasses to handle specific genome assemblies and data formats.
When processing new genomes, this method pre-processes the raw data to create optimized intermediate representations that improve computational efficiency for subsequent operations.

Set the `self._is_compiled` attribute to True after successful compilation.


**Source:** Line 84 in [feature.py](feature.py#L84)


---


### `compile(self, size_file=None, force=False)`


No documentation available.


**Source:** Line 124 in [feature.py](feature.py#L124)


---


### `compile(self, gtf_path, force=False)`


Compile TSS regions from GTF file.


**Source:** Line 200 in [feature.py](feature.py#L200)


---


### `compile(self, fasta, force=False)`


No documentation available.


**Source:** Line 247 in [feature.py](feature.py#L247)


---


### `data(self, fn=None)`


Lazy-loaded data property.


**Source:** Line 76 in [feature.py](feature.py#L76)


---


### `get_tss_region_from_gtf(gtf_path, cache_file=None)`


Get TSS region from GTF file.
Input:
    gtf_path: path to the GTF file; string
    cache_file: path to the cache file; string or None
Output:
    TSS table; pd.DataFrame with columns:
        gene_id, gene_name, transcript_id, seqname, txStart, source
    Same gene will have multiple rows if it has multiple transcripts.


**Source:** Line 170 in [feature.py](feature.py#L170)


---


### `load(self, fn=None)`


Load the feature data to inner data structure from the compiled file.

Sets the `self.data` attribute to the loaded data.


**Source:** Line 99 in [feature.py](feature.py#L99)


---


### `load(self, fn=None)`


Load chromosome sizes from compiled file.


**Source:** Line 132 in [feature.py](feature.py#L132)


---


### `load(self)`


Load TSS regions from compiled file.


**Source:** Line 227 in [feature.py](feature.py#L227)


---


### `load(self)`


Load CpG regions from compiled file.


**Source:** Line 276 in [feature.py](feature.py#L276)


---


## genome

**File:** `genome.py`

**Public functions:** 11


### `append_bins(self, df: pd.DataFrame, binsize: int, chr1_col: str='chr1', pos1_col: str='pos1', chr2_col: str=None, pos2_col: str=None, flavor: str='hickit')`


Append bins to dataframe.

Input:
    df: input dataframe
    binsize: binsize, see GenomeIdeograph.bins
    chr1_col: column name for first chromosome
    pos1_col: column name for first position
    chr2_col: column name for second chromosome (optional)
    pos2_col: column name for second position (optional)
    flavor: flavor of the binning, see GenomeIdeograph.bins
    
Output:
    df_e: dataframe with new bin columns
        For pairs mode: adds "chrom1","start1","end1","chrom2","start2","end2"
        For single mode: adds "chrom1","start1","end1"


**Source:** Line 410 in [genome.py](genome.py#L410)


---


### `bins(self, binsize: int, bed=False, order=False, flavor='hickit', within_chromosome=True)`


Get binned reference(IntervalIndex version)
Input:
    binsize: int
    bed: bool, if True, return bed format
    order: bool, if True, sort by chromosome and start position
    flavor: str, binning flavor
    within_chromosome: bool, if True, bins are confined within chromosomes;
        if False, bins can span multiple chromosomes
Output:
    If within_chromosome=True: intervals of each bin(dict of IntervalIndex)
    If within_chromosome=False: list of regions, each as [(start_chrom, start_pos), (end_chrom, end_pos)]


**Source:** Line 354 in [genome.py](genome.py#L354)


---


### `breaks(self, binsize: int, flavor='hickit', within_chromosome=True)`


Get binned reference(int version)
Input:
    binsize: int
    flavor: str, "hickit", "bedtools", "cooler_compat", "default", or "nearest_boundary"
    within_chromosome: bool, if True, each bin is confined within a chromosome;
        if False, bins can span multiple chromosomes with boundary handling specified by flavor
Return:
    If within_chromosome=True: breaks of bins(dict of list)
    If within_chromosome=False: list of regions, each as [(start_chrom, start_pos), (end_chrom, end_pos)]
Note:
    For hickit-flavored binning, the size of the last bin 
    >= 0.5 * binsize and < 1.5 * binsize. For bedtools-flavored binning, 
    the size of the last bin > 0 and <= binsize.
    For cooler_compat-flavored binning, the size of the last bin
    > binsize is trimmed to binsize.
    When within_chromosome is False, flavor can be:
        - "default": no boundary attachment
        - "nearest_boundary": snap to nearest chromosome boundary


**Source:** Line 236 in [genome.py](genome.py#L236)


---


### `coarsen_grouper(self, binsize1: int, binsize2: int, flavor: str='hickit')`

**Returns:** `pd.Series`


Get grouper to map from binsize1 to a bigger binsize2.
Input:
    binsize1: int
    binsize2: int
    flavor: str, "hickit" or "bedtools"
Output:
    pd.Series, 2-level index (chrom, start) from binsize1 as index,
        tuple (chrom, start) from binsize2 as value


**Source:** Line 478 in [genome.py](genome.py#L478)


---


### `index_merge_hom(df)`


Cancel ploidiness in index. You can then merge homologous chromosomes use functions like groupby.
This will change index inplace.
Input:
    df: DataFrame with MultiIndex ["chrom", "start"], chom is ploidified, like "chr1(mat)"
Output:
    a dataframe with modified level 0 index value, like "chr1(mat)" -> "chr1"


**Source:** Line 643 in [genome.py](genome.py#L643)


---


### `join_pixel_id(self, df, binsize, intra=True, filep=None)`


Add a 'pixel_id' column to the input DataFrame in a vectorized manner.
Input:
    df: DataFrame, bedpe-like with columns 'chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'
    binsize: int, size of the bins used in pixelation
    intra: bool, if True, process intra-chromosomal pixels
Output:
    DataFrame with an additional 'pixel_id' column


**Source:** Line 533 in [genome.py](genome.py#L533)


---


### `join_positions(self, df, binsize, intra=True)`


Add columns 'chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2' to the input DataFrame in a vectorized manner.
Input:
    df: DataFrame, bedpe-like with column 'pixel_id'
    binsize: int, size of the bins used in pixelation
Output:
    DataFrame with additional columns 'chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'


**Source:** Line 562 in [genome.py](genome.py#L562)


---


### `new_start(row, chrom_intervals=None)`


No documentation available.


**Source:** Line 489 in [genome.py](genome.py#L489)


---


### `parse_ucsc_region(region_str)`


Parses a UCSC-style region string and returns a tuple of (chromosome, start, end).

Parameters:
    region_str (str): A UCSC-style region string like 'chr1:2,345-6,789' or 'chrX:2345-6789'.

Returns:
    tuple: A tuple containing the chromosome (str), start position (int), and end position (int).

Raises:
    ValueError: If the input string is not in a valid UCSC region format.


**Source:** Line 12 in [genome.py](genome.py#L12)


---


### `pixel_id(self, row, binsize: int, intra=True)`


Get pixel_id from row
Input:
    row: pd.Series, assume bedpe format
    binsize: int
Output:
    int


**Source:** Line 511 in [genome.py](genome.py#L511)


---


### `sort_chrom(df, genome)`


Sort dataframe by natural chromosome order.
TODO:
    1. ask whether to rm chromint after sorting
    2. ask whether to reset index
Input:
    df: must have a column named "chrom"
Output:
    df: sorted by natural chromosome order


**Source:** Line 662 in [genome.py](genome.py#L662)


---


## heuristic_ordering

**File:** `heuristic_ordering.py`

**Public functions:** 1


### `ra_ordering(filesp, threads=24)`


Random Annealing ordering of contact-decay-profiles.
Input:
    filesp : DataFrame, must have pairs col
Output:
    DataFrame with additional col "order_index"


**Source:** Line 6 in [heuristic_ordering.py](heuristic_ordering.py#L6)


---


## hicbr

**File:** `hicbr.py`

**Public functions:** 6


### `dump(self, h5ad_file: str, output_dir: str, layers: Optional[List[str]]=None)`

**Returns:** `Dict[str, Any]`


Convert h5ad file to generic file formats.

Parameters
----------
h5ad_file : str
    Path to input h5ad file
output_dir : str
    Directory to save converted files
layers : list, optional
    List of layers to convert. If None, converts all supported layers.
    
Returns
-------
dict
    Metadata about the conversion process


**Source:** Line 87 in [hicbr.py](hicbr.py#L87)


---


### `files_to_anndata(input_dir: str, layers: Optional[List[str]]=None)`

**Returns:** `Any`


Load converted files back into an AnnData object.

Parameters
----------
input_dir : str
    Directory containing converted files
layers : list, optional
    List of layers to load. If None, loads all available layers.
    
Returns
-------
anndata.AnnData
    Reconstructed AnnData object


**Source:** Line 777 in [hicbr.py](hicbr.py#L777)


---


### `h5ad_to_files(h5ad_file: str, output_dir: str, layers: Optional[List[str]]=None)`

**Returns:** `Dict[str, Any]`


Convert h5ad file to generic file formats.

Parameters
----------
h5ad_file : str
    Path to input h5ad file
output_dir : str
    Directory to save converted files
layers : list, optional
    List of layers to convert. If None, converts all supported layers.
    
Returns
-------
dict
    Metadata about the conversion process
TODO:
    Fix: pd.Index type in obs with col name xxx now will be dumped as xxx_mask, xxx_values in csv


**Source:** Line 752 in [hicbr.py](hicbr.py#L752)


---


### `load(self, input_dir: str, layers: Optional[List[str]]=None)`

**Returns:** `Any`


Load converted files back into an AnnData object.

Parameters
----------
input_dir : str
    Directory containing converted files
layers : list, optional
    List of layers to load. If None, loads all available layers.
    
Returns
-------
anndata.AnnData
    Reconstructed AnnData object


**Source:** Line 155 in [hicbr.py](hicbr.py#L155)


---


### `register_dump_handler(self, layer_name: str, handler_func)`


Register a custom dump handler for a specific layer.

Parameters
----------
layer_name : str
    Name of the layer to register handler for
handler_func : callable
    Function that handles dumping this layer


**Source:** Line 61 in [hicbr.py](hicbr.py#L61)


---


### `register_load_handler(self, layer_name: str, handler_func)`


Register a custom load handler for a specific layer.

Parameters
----------
layer_name : str
    Name of the layer to register handler for
handler_func : callable
    Function that handles loading this layer


**Source:** Line 74 in [hicbr.py](hicbr.py#L74)


---


## hicio

**File:** `hicio.py`

**Public functions:** 31


### `align_to_df(primary_views)`


Extract general dataframe from primary_view function output.
Input:
    primary_views: complex data structure list of dict
Output:
    dataframe


**Source:** Line 941 in [hicio.py](hicio.py#L941)


---


### `combine_10x(annoted_file: dict)`

**Returns:** `pd.DataFrame`


Extract data from files and combine them into a expression matrix.
Input:
    annoted_file: {"matrix":matrix_file, "genes":gene, "barcodes":barcodes}
Output:
    a cell * gene matrix


**Source:** Line 636 in [hicio.py](hicio.py#L636)


---


### `decode_sam_flag(flag: int)`

**Returns:** `Dict[str, bool]`


Decode SAM bitwise FLAG into human-readable components.

Args:
    flag: Integer FLAG value from SAM field 2
    
Returns:
    Dictionary with flag component names as keys and boolean values
    
Reference: SAM v1.6 specification


**Source:** Line 526 in [hicio.py](hicio.py#L526)


---


### `divide_name(filename)`


No documentation available.


**Source:** Line 37 in [hicio.py](hicio.py#L37)


---


### `dump_json(obj, filep)`


Dump to json shortcut.


**Source:** Line 976 in [hicio.py](hicio.py#L976)


---


### `dump_pickle(obj, filep)`


Dump to pickle shortcut.


**Source:** Line 994 in [hicio.py](hicio.py#L994)


---


### `flush(self)`


No documentation available.


**Source:** Line 28 in [hicio.py](hicio.py#L28)


---


### `get_chrom_contact_counts(dump_dir)`


Get per-chromosome-per-phase contacts counts of samples.
Input:
    dump_dir: rd/dump generated by hires_utils
Output:
    pd.DataFrame


**Source:** Line 930 in [hicio.py](hicio.py#L930)


---


### `get_ref_dir()`


Return a path to src's ref


**Source:** Line 1000 in [hicio.py](hicio.py#L1000)


---


### `load_cool(cool, root='/')`


Load simple cooler file as sparsematrix

Parameters
----------
cool : str
    Path to the input .cool file.

Returns
-------
mat : scipy coo_matrix
    Hi-C contact map in COO format.
frags : pandas DataFrame
    Table of bins matching the matrix.
chroms : pandas DataFrame
    Table of chromosome informations.


**Source:** Line 719 in [hicio.py](hicio.py#L719)


---


### `load_json(filep)`


Load from json file shortcut.


**Source:** Line 982 in [hicio.py](hicio.py#L982)


---


### `load_pickle(filep)`


Load from pickle file shortcut.


**Source:** Line 988 in [hicio.py](hicio.py#L988)


---


### `match_10x(file_list: list, retind=False)`

**Returns:** `dict`


Check if input file names contains a 10x result.
Input:
    file_list: list of file path
    retind: return index of matched file instead of file path
Output:
    {"matrix":matrix_file, "genes":gene_file, "barcodes":barcodes_file}


**Source:** Line 603 in [hicio.py](hicio.py#L603)


---


### `matra(file)`


Read umi-tools output to anndata object.


**Source:** Line 919 in [hicio.py](hicio.py#L919)


---


### `parse_bed(bed_path: str)`

**Returns:** `pd.DataFrame`


Load and parse BED format files, automatically assigning BED12 column names.

Inputs:
    bed_path (str): Path to the BED file.
    
Returns:
    pd.DataFrame: Parsed BED data with columns named according to BED12 specification:
        - chrom (str): Chromosome name.
        - start (int): Start coordinate (0-based).
        - end (int): End coordinate (exclusive).
        - name (str, optional): Feature name.
        - score (float, optional): Score value (0-1000).
        - strand (str, optional): Strand ('+' or '-').
        - thickStart (int, optional): Thick start coordinate.
        - thickEnd (int, optional): Thick end coordinate.
        - itemRgb (str, optional): RGB color (e.g., '255,0,0').
        - blockCount (int, optional): Number of blocks.
        - blockSizes (str, optional): Comma-separated block sizes.
        - blockStarts (str, optional): Comma-separated block start positions.
        
Notes:
    - Automatically detects the number of columns in the BED file
    - Assigns BED12 column names in order: chrom, start, end, name, score, 
      strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts
    - If the file has fewer than 12 columns, only the first N columns are named
    - If the file has more than 12 columns, only the first 12 are read
    - Performs automatic type conversion for numeric columns
    - Coordinates are 0-based and end-exclusive by BED convention


**Source:** Line 47 in [hicio.py](hicio.py#L47)


---


### `parse_gff(file, ID=False, Name=False)`


Parsing gff file. Read all in memory. Extract ID(gene), and Name if set true.


**Source:** Line 248 in [hicio.py](hicio.py#L248)


---


### `parse_gtf(file, ID=True, name=True, mgi_id=True)`


Parsing gtf file. Read all in memory. Extract gene_id to df if set true.
Input:
    file: path to .gtf file.
    ID: whether to parse gene_id from attribute string.
    name: whether to parse gene_name from attribute string.
    **args: arguments passed to pd.read_table
Return:
    pd.DataFrame


**Source:** Line 220 in [hicio.py](hicio.py#L220)


---


### `parse_hicluster_res(embed, sample_table)`


Read schicluster embedding result into dataframe.
Input:
    embed: schicluster concat-cell output hdf5 file
    sample_table: input sample file of schicluster pipeline or list of sample names
Examples:
    from hic_basic.io import parse_hicluster_res
    from hic_basic.scAB_embedding import do_umap
    hicluster_res = parse_hicluster_res(embed, sample_table)
    umap_res = do_umap(hicluster_res.values)


**Source:** Line 754 in [hicio.py](hicio.py#L754)


---


### `parse_pairs(filename: str)`

**Returns:** `pd.DataFrame`


read from 4DN's standard .pairs format
compatible with all hickit originated pairs-like format


**Source:** Line 125 in [hicio.py](hicio.py#L125)


---


### `parse_pairs_like(filename: str)`

**Returns:** `pd.DataFrame`


read from 4DN's standard .pairs format
compatible with all hickit originated pairs-like format


**Source:** Line 175 in [hicio.py](hicio.py#L175)


---


### `parse_sam_line(sam_line: Union[str, List[str]])`

**Returns:** `pd.Series`


Parse a SAM format line into a pandas Series with standardized fields and optional tags.

This function parses a SAM (Sequence Alignment/Map) format record, extracting both
the 11 mandatory fields and any optional TAG:TYPE:VALUE fields according to the
SAM specification v1.6+.

Args:
    sam_line: A SAM format line as a string or pre-split list of strings.
              If string, it will be split by tabs.
              
Returns:
    pd.Series: A pandas Series where:
               - Index includes: 'QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR',
                 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL' (mandatory fields)
               - Additional indices for any optional tags (e.g., 'NM', 'MD', 'AS')
               - Values are appropriately typed (int, float, str) based on SAM type codes
               
Raises:
    ValueError: If the input doesn't contain at least 11 fields (mandatory SAM fields)
    
Notes:
    - POS, PNEXT: 1-based coordinates (0 indicates unmapped)
    - FLAG: Bitwise flag, returned as integer for further decoding
    - Optional tags follow pattern: TAG:TYPE:VALUE
    - TYPE codes: i (int), f (float), Z (string), A (char), H (hex), B (byte array)
    - For 'B' type (byte arrays), values are returned as list of appropriate type
    
Examples:
    >>> sam_line = "read1   16      chr1    100     255     30M     =       200     50      ACGT    IIII    NM:i:1  MD:Z:30"
    >>> series = parse_sam_line(sam_line)
    >>> series['QNAME']  # Returns 'read1'
    >>> series['FLAG']   # Returns 16
    >>> series['NM']     # Returns 1
    >>> series['MD']     # Returns '30'


**Source:** Line 340 in [hicio.py](hicio.py#L340)


---


### `parse_vcf(vcf_path, region, return_alleles=False)`


Extracts genotype data for all samples within a specified genomic region from a VCF file.
Optionally returns genotypes in allele form (e.g., 'GA/GA' instead of '1/1').

Parameters:
    vcf_path (str): Path to the VCF file.
    chromosome (str): Target chromosome (e.g., '1', 'chr1').
    start (int): Start position (1-based).
    end (int): End position (1-based).
    return_alleles (bool): If True, returns genotype in allele form. Default is False.

Returns:
    pd.DataFrame: A DataFrame where each row represents a variant and each column represents a sample.
                  Genotypes are returned either as indices (e.g., '1/1') or as allele strings (e.g., 'GA/GA').


**Source:** Line 272 in [hicio.py](hicio.py#L272)


---


### `read_10x(fp, gene_column=1, cell_column=0)`

**Returns:** `pd.DataFrame`


Read 10x-flavor matrix market file.
Input:
    fp: file path
    gene_column: use this column from gene/feature matrix
    cell_column: use this column from cell matrix
Output:
    a cell * gene matrix


**Source:** Line 645 in [hicio.py](hicio.py#L645)


---


### `read_expr(path, sep=None)`

**Returns:** `pd.DataFrame`


Read expression matrix from file.
Don't guarantee OV(obs*var) or VO(var*obs) output. It just read in the file.
Input:
    path: path to matrix file
    sep: separator, if None will infer from file extension
Output:
    expression matrix


**Source:** Line 557 in [hicio.py](hicio.py#L557)


---


### `read_h5_ds(group)`


Read .h5 bottom level group into dataframe.
Input:
    group: h5py.Group; must have same length datasets.
Output:
    pd.DataFrame


**Source:** Line 704 in [hicio.py](hicio.py#L704)


---


### `read_meta(fp)`


Read general metadata, take care of sample_name.


**Source:** Line 892 in [hicio.py](hicio.py#L892)


---


### `read_umi_tools(path, sep=None)`

**Returns:** `pd.DataFrame`


Read umi_tools long-form output matrix.
A wrapper for read_expr.
By default umi_tools file store matrix in VO format, this function will transpose it to OV.
See read_expr for more details.
Input:
    path: path to matrix file
    sep: separator, if None will infer from file extension
Output:
    obs * var matrix


**Source:** Line 906 in [hicio.py](hicio.py#L906)


---


### `schicluster2mat(filei)`


Read in schicluster impute .hdf5 or np.savez_compressed .npz file to sparse matrix.

This function supports two file formats:
1. HDF5 files with standard schicluster sparse matrix structure
2. NPZ files saved using np.savez_compressed with 'data' key

Args:
    filei: Path to input file. Supported formats:
           - .hdf5: schicluster imputed HDF5 file
           - .npz: Compressed numpy archive with sparse matrix components

Returns:
    csr_matrix: Scipy sparse matrix in CSR format

Raises:
    ValueError: If file format is not supported or required keys are missing
    OSError: If file cannot be opened or read


**Source:** Line 781 in [hicio.py](hicio.py#L781)


---


### `schiclusterDir2mat(hiclusterdir)`


Read in schicluster imputed .hdf5 files (from a dir) to sparse matrix.
Input:
    hicluster: schicluster imputed .hdf5 file dir. e.g. f"/shareb/ychi/repo/sperm22_schicluster/imputed_matrix/100000/{chrom}/"
Output:
    sparse matrix


**Source:** Line 849 in [hicio.py](hicio.py#L849)


---


### `write(self, _)`


No documentation available.


**Source:** Line 25 in [hicio.py](hicio.py#L25)


---


### `write_triplet(sparseM, filep, max_coo=False, zipping=True)`


Write sparse matrix to triplet txt, assume symmetry.
For scipy doesn't give a way to write triplet format directly.
Input:
    sparseM: scipy sparse matrix
    filesp: output file path
    max_coo: if True, write largest possible coordinate at first row
Output:
    indi, indj, data; separated by tab


**Source:** Line 861 in [hicio.py](hicio.py#L861)


---


## impute.schicluster

**File:** `impute/schicluster.py`

**Public functions:** 3


### `schicluster_imputation_for_mat(mat, alpha=0.05, kernel_size=3, sigma=2, if_convolve=True)`


No documentation available.


**Source:** Line 33 in [impute/schicluster.py](impute/schicluster.py#L33)


---


### `schicluster_impute(coolp, region, fo, max_dist: Optional[int]=20000000, alpha: float=0.05, kernel_size: int=3, sigma: int=2, if_convolve: bool=True)`


Impute Hi-C matrix using schicluster method from cooler file.
Input:
    coolp: cooler file
    region: region to impute, normally a chromosome,
        otherwise maybe incompatible with downstream analysis.
        See cool2mat's region argument
    fo: output file, str
    max_dist: max linear distance, int
    alpha: alpha, float
    kernel_size: kernel size, int
    sigma: sigma, float
    if_convolve: if convolve, bool
Output:
    write imputed matrix to fo


**Source:** Line 49 in [impute/schicluster.py](impute/schicluster.py#L49)


---


### `solve_rwr_inverse(stoch_matrix, alpha=0.05)`


No documentation available.


**Source:** Line 10 in [impute/schicluster.py](impute/schicluster.py#L10)


---


## impute.simpute

**File:** `impute/simpute.py`

**Public functions:** 9


### `DMimpute(_3dg_path, fo, genome='GRCh38', binsize=20000, agg_dip=False, flavor='cooler_compat', max_dist=None)`


Generate distance matrix (store in cooler) from 3dg file.
Only cis region is considered.
Input:
    _3dg_path: str, path to 3dg file
    genome: give genome name or length path and add pixel_id to output
        if None, pixel_id will not be added
    fo: str, path to output parquet file
    max_dist: int, pixels with distance larger than max_dist will be discarded
        if None, all pixels will be returned
    binsize: int, binsize of input 3dg file
    agg_dip: if True, aggregate diploid structure
Output:
    fo filepath


**Source:** Line 402 in [impute/simpute.py](impute/simpute.py#L402)


---


### `boolean_radius_neighbor(df, min_dist=3, n_jobs=4, pixels=True)`


generate boolean matrix from dataframe
Input:
    df: dataframe, index is position, columns are x, y, z
    min_dist: int, minimum distance to be considered as neighbor
    n_jobs: int, number of jobs to run in parallel
    pixels: bool, if True, return pixels, otherwise return graph
Output:
    if pixels:
        dataframe: bedpe-like dataframe, columns are chrom1, start1, chrom2, start2, count
    else:
        graph: CSR format boolean matrix, True if two points are neighbors


**Source:** Line 21 in [impute/simpute.py](impute/simpute.py#L21)


---


### `chrom_dist_mat_long(chunk)`


No documentation available.


**Source:** Line 240 in [impute/simpute.py](impute/simpute.py#L240)


---


### `cis_distance_graph(_3dg_path, fo=None, genome=None, max_dist=2000000, binsize=20000, n_jobs=None)`


Generate distance matrix (store in bedpe-like format) from 3dg file.
Only cis region is considered.
Input:
    _3dg_path: str, path to 3dg file
    fo: str, path to output parquet file, if None, return dataframe
    genome: give genome name or length path and add pixel_id to output
        if None, pixel_id will not be added
    max_dist: int, pixels with distance larger than max_dist will be discarded
    binsize: int, binsize of input 3dg file
    n_jobs: int, number of jobs to run in parallel
Output:
    df: n * 7 dask dataframe, (chrom1, start1, end1, chrom2, start2, end2, distance)
    None if fo is not None


**Source:** Line 216 in [impute/simpute.py](impute/simpute.py#L216)


---


### `cis_distance_graph_df(_3dg_path, chrom=None, genome=None, fo=None, max_dist=2000000, binsize=20000, fill=False)`


Generate distance matrix (store in bedpe-like format) from 3dg file.
TODO: fix max_dist, now is using 3d distance mistakenly!
Only cis region is considered.
pandas.DataFrame version, all in memory.
Input:
    _3dg_path: str, path to 3dg file
    chrom: str, chromosome to be processed.
        if None, process all chromosomes(highly not recommended, may cause memory error)
    genome: give genome name or length path and add pixel_id to output
        if None, pixel_id will not be added
    fo: str, path to output parquet file, if None, return dataframe
    max_dist: int, pixels with distance larger than max_dist will be discarded
        if None, all pixels will be returned
    binsize: int, binsize of input 3dg file
    fill: return all possible pixels (< max_dist of course), fill missing pixels with NaN
        only valid when both chrom and genome are specified
        will return dask array
Output:
    df: n * 7 dataframe, (chrom1, start1, end1, chrom2, start2, end2, distance)
    None if fo is not None


**Source:** Line 281 in [impute/simpute.py](impute/simpute.py#L281)


---


### `cis_proximity_graph(_3dg_path, fo, min_dist=3, genome='mm10', binsize=20000, agg_dip=False, flavor='hickit', n_jobs=4)`


Genereate 0,1 cooler file of region proximity from 3dg file.
Only cis region is considered.
Input:
    _3dg_path: str, path to 3dg file
    fo: str, path to output cooler file
    min_dist: int, minimum distance to be considered as neighbor
    genome: str, genome version
    binsize: int, binsize
    agg_dip: if True, aggregate diploid structure
    flavor: str, flavor of binnify, see GenomeIdeograph.bins
Output:
    None


**Source:** Line 122 in [impute/simpute.py](impute/simpute.py#L122)


---


### `parse_3dg_dask(_3dg_path)`


Parse 3dg file to dask dataframe.
Input:
    _3dg_path: str, path to 3dg file
Output:
    structure: dask dataframe, index is (chr, pos), columns are x, y, z


**Source:** Line 194 in [impute/simpute.py](impute/simpute.py#L194)


---


### `process_chunk(chunk, chrom=None)`


No documentation available.


**Source:** Line 317 in [impute/simpute.py](impute/simpute.py#L317)


---


### `process_chunk(chunk, chrom=None)`


No documentation available.


**Source:** Line 421 in [impute/simpute.py](impute/simpute.py#L421)


---


## inter_contact

**File:** `inter_contact.py`

**Public functions:** 3


### `cell_sig(file_path: str)`

**Returns:** `pd.Series`


No documentation available.


**Source:** Line 32 in [inter_contact.py](inter_contact.py#L32)


---


### `expected_inter(pairs: pd.DataFrame)`

**Returns:** `pd.DataFrame`


No documentation available.


**Source:** Line 6 in [inter_contact.py](inter_contact.py#L6)


---


### `pick_triu(matrix: pd.DataFrame)`

**Returns:** `pd.DataFrame`


No documentation available.


**Source:** Line 19 in [inter_contact.py](inter_contact.py#L19)


---


## interpolate

**File:** `interpolate.py`

**Public functions:** 3


### `correct_pseudotime(adata, pseudotime_col, n_pcs=15, inplace=True)`


Correct pseudotime of adata to make it linearly correlated with a distance metrix.
Inspired by Alpert et al. (2018)
Input:
    adata: AnnData object, must have adata.obs[pseudotime_col]
    pseudotime_col: pseudotime column name
Output:
    new adata object


**Source:** Line 90 in [interpolate.py](interpolate.py#L90)


---


### `gaussian_interpolate(adata, pseudotime_col, obs_using, winSz=0.1, numPts=200, inplace=True)`


Interpolation of all concerned traits in adata.
Input:
    adata: AnnData object, must have adata.obs[pseudotime_col]
    pseudotime_col: pseudotime column name
    obs_using: list of obs column names to be interpolated
    winSz: window size of the interpolation
    numPts: number of desired interpolated points
Output:
    new adata object numPts * adata.shape[1]


**Source:** Line 63 in [interpolate.py](interpolate.py#L63)


---


### `gaussian_weights(x, y)`


No documentation available.


**Source:** Line 20 in [interpolate.py](interpolate.py#L20)


---


## lap_embedding

**File:** `lap_embedding.py`

**Public functions:** 3


### `circle_arctan(cordinate: pd.Series)`

**Returns:** `float`


No documentation available.


**Source:** Line 14 in [lap_embedding.py](lap_embedding.py#L14)


---


### `get_distance_matrix(cell_cdps: pd.DataFrame)`

**Returns:** `pd.DataFrame`


No documentation available.


**Source:** Line 7 in [lap_embedding.py](lap_embedding.py#L7)


---


### `spectral_ordering(distance_matrix: pd.DataFrame)`

**Returns:** `pd.DataFrame`


No documentation available.


**Source:** Line 31 in [lap_embedding.py](lap_embedding.py#L31)


---


## mchr

**File:** `mchr.py`

**Public functions:** 7


### `DI(self, groupA_samples, groupB_samples, max_linear_dist=2000000, max_3d_dist=5, fdrcut=0.05, chrom_blacklist: list=['chrX', 'chrY'], n_jobs=None, mem: int=4)`

**Returns:** `pd.DataFrame`


Simplediff method to compute differential interactions.
Input:
    groupA_samples: sample names of group A
    groupB_samples: sample names of group B
    max_linear_dist: max linear distance to consider
    max_3d_dist: max 3D distance to consider, if any of groupA and groupB have mean distance
        larger than this value in pixel, this pixel will be filtered out
    fdrcut: FDR cutoff, if None, no multiple testing correction and will give all valid pixels
        if set, only pixels with FDR < fdrcut will be returned
    chrom_blacklist: chromosomes to exclude
    n_jobs: number of threads for parallel computing
    mem: memory limit for each worker(GB)
Output:
    DI: differential interactions, index and columns are 3-level multiindex
    (chrom, start1, start2)


**Source:** Line 358 in [mchr.py](mchr.py#L358)


---


### `DM(self, samples: list, region_pair: list, min_samples=None, n_jobs=None)`

**Returns:** `pd.DataFrame`


Compute distance matrix.
Input:
    samples: sample names to compute
    region_pair: genome region to compute DM, can include inter-chromosome regions,
        can cross chromosomes.
        [[(chrom1,pos1),(chrom1,pos2)],[(chrom2,pos1),(chrom2,pos2)]
    proximity: if not None, calculate proximity instead of distance
        count pixels within proximity
    n_jobs: number of threads for parallel computing
        if None, use a simple in-memory method
        else, use dask
Output:
    DM: distance matrix, index and columns are 2-level multiindex
    (chrom, start)


**Source:** Line 310 in [mchr.py](mchr.py#L310)


---


### `IS(self, samples: list, region_pair: list, w: int)`

**Returns:** `pd.DataFrame`


Compute insulation score for each sample.
NOTE: this method is not recommended for large regions.
    Will be replaced by a more efficient method.
Input:
    samples: sample names to compute
    region_pair: genome region to compute IS, can only use intra-chromosome regions here
    w: window size
Output:
    IS: insulation score, index is 2-level multiindex
    (chrom, start), columns are sample names


**Source:** Line 277 in [mchr.py](mchr.py#L277)


---


### `PM(self, samples: list, region_pair: list, proximity: int, min_samples=None, n_jobs=None)`

**Returns:** `pd.DataFrame`


Compute proximity matrix.
Input:
    samples: sample names to compute
    region_pair: genome region to compute PM, can include inter-chromosome regions,
        can cross chromosomes.
        [[(chrom1,pos1),(chrom1,pos2)],[(chrom2,pos1),(chrom2,pos2)]
    proximity: count pixels within proximity
    n_jobs: number of threads for parallel computing
        if None, use a simple in-memory method
        else, use dask
Output:
    PM: proximity matrix, index and columns are 2-level multiindex
    (chrom, start)


**Source:** Line 334 in [mchr.py](mchr.py#L334)


---


### `from_netcdf(cls, _3dg_ds_fp: str, ddir=None, _3dg_files=None, samples=None, genome=None, binsize=None, flavor=None)`


Initialize a Mchr object from a netcdf file.
Input:
    _3dg_ds_fp: netcdf file path
Output:
    Mchr object


**Source:** Line 226 in [mchr.py](mchr.py#L226)


---


### `multiDM(self, samples: list, region_pair: list)`

**Returns:** `xr.DataArray`


Compute distance matrix for each single sample.
NOTE: this method is not recommended for large regions.
A better way is to aggregate DM to 1D feature in a split-apply-combine manner.
Input:
    samples: sample names to compute
    region_pair: genome region to compute DM, can include inter-chromosome regions,
        can cross chromosomes.
        [[(chrom1,pos1),(chrom1,pos2)],[(chrom2,pos1),(chrom2,pos2)]
    n_jobs: number of threads for parallel computing
        if None, use a simple in-memory method
        else, use dask
Output:
    xarray DataArray, sample_name as a dimension


**Source:** Line 254 in [mchr.py](mchr.py#L254)


---


### `multiple_testing_correction(pvalues, correction_type='FDR')`


Consistent with R - print
correct_pvalues_for_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05,
                                      0.069, 0.07, 0.071, 0.09, 0.1])
from https://github.com/CoBiG2/cobig_misc_scripts/blob/master/FDR.py


**Source:** Line 21 in [mchr.py](mchr.py#L21)


---


## metrics

**File:** `metrics.py`

**Public functions:** 4


### `distance_pca_euclid(mat: np.ndarray, n_components: int, with_std=False)`

**Returns:** `np.ndarray`


Euclid distances in PCA-space.
Input:
    mat: m( sample ) * n( feature ) matrix
    n_components: number of PCA components used to calculate distance.
Output:
    m * m distance matrix


**Source:** Line 24 in [metrics.py](metrics.py#L24)


---


### `kernel_c_rbf(mat)`


No documentation available.


**Source:** Line 51 in [metrics.py](metrics.py#L51)


---


### `kernel_pca_euclid(mat: np.ndarray, n_components: int)`

**Returns:** `np.ndarray`


Affinity based on PCA-euclid-distance
Input:
    mat: m( sample ) * n( feature ) matrix
    n_components: number of PCA components used to calculate distance.
Output:
    m * m similarity matrix


**Source:** Line 39 in [metrics.py](metrics.py#L39)


---


### `pairwise_DKL(df)`


compute pairwise Kullback-Leibler Divergence.
resulting D'(i,j) = D'(j,i) = D_KL(i,j) + D_KL(j,i).
Input:
    df: sample x feature matrix
Output:
    matrix: symdist


**Source:** Line 55 in [metrics.py](metrics.py#L55)


---


## nagano_cycle_phasing

**File:** `nagano_cycle_phasing.py`

**Public functions:** 6


### `Nagano_ordering(metrics: pd.DataFrame)`


No documentation available.


**Source:** Line 91 in [nagano_cycle_phasing.py](nagano_cycle_phasing.py#L91)


---


### `contact_describe(cell_name: str, c1=1, p1=2, c2=3, p2=4)`

**Returns:** `pd.Series`


No documentation available.


**Source:** Line 11 in [nagano_cycle_phasing.py](nagano_cycle_phasing.py#L11)


---


### `order_sample(ordering_values)`


No documentation available.


**Source:** Line 140 in [nagano_cycle_phasing.py](nagano_cycle_phasing.py#L140)


---


### `order_sample_within_group(chunk, shift)`


No documentation available.


**Source:** Line 133 in [nagano_cycle_phasing.py](nagano_cycle_phasing.py#L133)


---


### `schic_spectral_embedding(df)`


non-linear dimensionality reduction of cdps curve
Input:
    df: sample x feature matrix
Output:
    2nd 3rd smallest eigenvectors of the graph Laplacian


**Source:** Line 164 in [nagano_cycle_phasing.py](nagano_cycle_phasing.py#L164)


---


### `smooth_cluster(cdps, annote, cluster_col='named_cluster', batchsize=5, pca_n=6)`


# fine-tuned group assignment using k-means clustering voting
Input:
    cdps: contact decay profiles
    annote: df that records original group assignment, must using sample names as index
    cluster_col: colname in annote that records originmal group assignment
    batchsize: number of samples in each voting batch
    pca_n: number of components used in distance calculation
Output:
    annote with 2 new col: km_label km_label_group


**Source:** Line 60 in [nagano_cycle_phasing.py](nagano_cycle_phasing.py#L60)


---


## paracalc

**File:** `paracalc.py`

**Public functions:** 2


### `mt(func)`


No documentation available.


**Source:** Line 22 in [paracalc.py](paracalc.py#L22)


---


### `wrapper(iterable)`


No documentation available.


**Source:** Line 24 in [paracalc.py](paracalc.py#L24)


---


## phasing

**File:** `phasing.py`

**Public functions:** 22


### `check_allele(entry, rs, v, append_features=None, min_baseq=20, verbose=0)`


Check if a SAM entry contains a specific allele at a given position.

Args:
    entry (list): A list representing a SAM format entry where:
                  entry[2] is the reference sequence name,
                  entry[5] is the CIGAR string,
                  entry[9] is the sequence,
                  entry[10] is the base quality string.
    rs (int): Reference start position (0-based).
    v (list): A list containing SNP information [position, snp_id, ref_allele, alt_allele].
    append_features (list, optional): List of additional features to append to the result.
    min_baseq (int, optional): Minimum base quality threshold. Defaults to 20.
    verbose (int, optional): Verbosity level. Defaults to 0.

Returns:
    list or None: A list containing [chrom, pos, allele, read_name, relative_pos] and additional features if the allele matches,
                  otherwise None.


**Source:** Line 406 in [phasing.py](phasing.py#L406)


---


### `convert_alignment_to_entry(alignment)`


Convert a pysam.AlignedSegment object to a list in the same format as a SAM entry.

Args:
    alignment (pysam.AlignedSegment): An alignment object from pysam.

Returns:
    list: A list representing the alignment in the same format as a SAM entry.


**Source:** Line 474 in [phasing.py](phasing.py#L474)


---


### `count_chrom_phased(filename: str)`

**Returns:** `pd.Series`


Count number of phased (with separable SNP) segments per chromosome from a .seg.gz file.
Input:
    filename: Path to the .seg.gz file.
Output:
    A pandas Series with multi-index (chromosome, phasing) and counts.


**Source:** Line 874 in [phasing.py](phasing.py#L874)


---


### `count_lines_in_file(file_path)`


Count the number of lines in a file efficiently.

Args:
    file_path (str): Path to the file.

Returns:
    int: Number of lines in the file.


**Source:** Line 505 in [phasing.py](phasing.py#L505)


---


### `count_phased_seg_perbin(filename, genome='mm10', binsize=1000000.0)`


Count the number of phased segments per bin of a given genome.
Input:
    filename: str, path to the .seg file
    genome: str, genome name (default: "mm10")
    binsize: int, size of the bins (default: 1e6)
Output:
    pandas DataFrame with columns: "chrom", "start", "phase", "count"


**Source:** Line 892 in [phasing.py](phasing.py#L892)


---


### `do_sam2vcf_allele_count(sam_file, output, ref_snp_file=None, marked_allele_file=None, force=False)`


Given a SAM file and a reference SNP file, generate a parquet file with allele counts at each SNP.

{pipeline}[Gamete SNP Phasing Pipeline][1.count alleles]:
    `do_sam2vcf_allele_count` for each gamete SAM file to get allele counts at each SNP position.

Parameters:
    sam_file (str): Path to the SAM file, can be sam.gz, bam, or sam.
    ref_snp_file (str): Path to the reference SNP file. Can be tab-delimited (.tsv, .txt, etc.)
                       or parquet format (.parquet). Expected columns: chrom, pos, ref_allele, alt_allele.
    marked_allele_file (str, optional): Path to the intermediate marked allele file.
    output (str): Path to the output parquet file.
    force (bool, optional): Whether to overwrite existing output files. Defaults to False.
Output (str): Path to the output parquet file.
Note:
    The output is aligned to ref_snp_file and includes all columns from the reference SNP table
    plus allele count fields and genotype (e.g., input_ref/input_alt, ref/alt/other, gt).


**Source:** Line 793 in [phasing.py](phasing.py#L793)


---


### `entry_align_pos(entry)`


Calculate reference and query start/end positions based on CIGAR string.

Args:
    entry (list): A list representing a SAM format entry where:
                  entry[1] is the FLAG field,
                  entry[3] is the reference start position,
                  entry[5] is the CIGAR string.

Returns:
    tuple: A tuple containing (rs, re, qs, qe) where:
           rs: Reference start position (0-based)
           re: Reference end position
           qs: Query start position
           qe: Query end position


**Source:** Line 334 in [phasing.py](phasing.py#L334)


---


### `find_intv(a: List, x: int)`

**Returns:** `int`


Find the interval containing a given value using binary search.

Args:
    a: List of intervals
    x: Value to search for
    
Returns:
    Index of the interval containing x, or -1 if not found


**Source:** Line 105 in [phasing.py](phasing.py#L105)


---


### `find_ovlp(a: List, st: int, en: int)`

**Returns:** `List`


Find intervals overlapping with a given range.

Args:
    a: List of intervals [start, end, ...]
    st: Start of the query range
    en: End of the query range
    
Returns:
    List of overlapping intervals


**Source:** Line 146 in [phasing.py](phasing.py#L146)


---


### `get_chrom_hap_score(dump_dir)`


Get chrom_hap_score of samples.
Input:
    dump_dir: rd/dump generated by hires_utils
Output:
    pd.Series


**Source:** Line 857 in [phasing.py](phasing.py#L857)


---


### `get_file_type(file_path)`


Determine the file type by reading the file content.

Args:
    file_path (str): Path to the file to check.

Returns:
    str: The file type ('SAM', 'BAM', 'SAMgz', or 'unknown').


**Source:** Line 233 in [phasing.py](phasing.py#L233)


---


### `hap_score_gb(df)`


No documentation available.


**Source:** Line 854 in [phasing.py](phasing.py#L854)


---


### `index_end(a: List, sorted: bool=True)`

**Returns:** `None`


Add index information for efficient overlap finding.

Args:
    a: List of intervals [start, end, ...]
    sorted: Whether the intervals are already sorted


**Source:** Line 71 in [phasing.py](phasing.py#L71)


---


### `is_sam_content(text, max_lines=20)`


检查文本内容是否符合 SAM 格式特征

Args:
    text (str): 要检查的文本内容
    max_lines (int): 最多检查的行数

Returns:
    bool: 如果是 SAM 格式返回 True，否则返回 False


**Source:** Line 300 in [phasing.py](phasing.py#L300)


---


### `merge(a: List, sorted: bool=True)`

**Returns:** `None`


Merge overlapping intervals.

Args:
    a: List of intervals to merge, each interval is [start, end]
    sorted: Whether the intervals are already sorted


**Source:** Line 42 in [phasing.py](phasing.py#L42)


---


### `mt_count_chrom_phased(filename: str)`

**Returns:** `pd.Series`


No documentation available.


**Source:** Line 889 in [phasing.py](phasing.py#L889)


---


### `parse_cigar(cigar: str)`

**Returns:** `List[Tuple[int, str]]`


Parse CIGAR string into operations and lengths.

Args:
    cigar: CIGAR string (e.g., "100M5D10I")
    
Returns:
    List of (length, operation) tuples (e.g., [(100, 'M'), (5, 'D'), (10, 'I')])


**Source:** Line 219 in [phasing.py](phasing.py#L219)


---


### `parse_phased_snp(phased_snp_file: str)`

**Returns:** `dict`


Parse a phased SNP file into a dictionary for interval-based lookups.

Args:
    phased_snp_file (str): Path to the phased SNP file. Can be tab-delimited (.tsv, .txt, etc.)
                          or parquet format (.parquet). Expected columns: chrom, pos, ref_allele, alt_allele.

Returns:
    dict: A dictionary where keys are chromosome names and values are lists of 
          [start, end, ref_allele, alt_allele, index] intervals for efficient overlap finding.


**Source:** Line 172 in [phasing.py](phasing.py#L172)


---


### `read_ref_snp_file(ref_snp_file)`


Read a reference SNP file, supporting both TSV and parquet formats.

Args:
    ref_snp_file (str): Path to the reference SNP file. Can be tab-delimited (.tsv, .txt, etc.)
                       or parquet format (.parquet).

Returns:
    pandas.DataFrame: A DataFrame with columns 'chrom', 'pos', 'input_ref', 'input_alt'.


**Source:** Line 670 in [phasing.py](phasing.py#L670)


---


### `sam_count_alleles(allele_df, ref_snp_file, gt_strategy='max_allele')`


Count reference, alternative, and other alleles for each chromosome and site.

This function takes a DataFrame containing genetic data and calculates the count
of reference alleles, alternative alleles, and other alleles for each unique
chromosome and position combination.
TODO: Better GT strategies. For example, is 0/1 better than 1/4?

Args:
    allele_df (pandas.DataFrame): A DataFrame with columns 'chrom', 'pos', 'allele',
                           'input_ref', and 'input_alt'. Each row represents
                           an allele observation at a specific genomic location.
    ref_snp_file (str): Path to the reference SNP file. Can be tab-delimited or parquet format.
    gt_strategy (str): Strategy for generating the GT column. 
                      "max_allele": Choose the allele with the highest count (ref or alt) and convert to ATGC
                      "no_conflict": Choose the non-zero allele if one is 0 and the other > 0, otherwise NA

Returns:
    pandas.DataFrame: A DataFrame with columns 'chrom', 'pos', 'ref', 'alt', 'other', 'gt'.
                      Each row represents a unique chromosome and position combination
                      with the counts of reference alleles ('ref'), alternative alleles ('alt'),
                      other alleles ('other'), and the genotype ('gt').


**Source:** Line 698 in [phasing.py](phasing.py#L698)


---


### `sam_mark_alleles(sam_file, phased_snp_file, outfile=None, append_features=None, min_mapq=20, min_baseq=20, show_progress=False, verbose=0)`


Mark alleles in SAM/BAM file based on phased SNP information.

Args:
    sam_file (str): Path to the input SAM or BAM file.
    phased_snp_file (str): Path to the phased SNP file. Can be tab-delimited (.tsv, .txt, etc.)
                          or parquet format (.parquet). Expected columns: chrom, pos, ref_allele, alt_allele.
    outfile (str, optional): Path to the output file. If None, returns a DataFrame. Defaults to None.
    append_features (list, optional): List of additional features to append to the output. Defaults to None.
    min_mapq (int, optional): Minimum mapping quality threshold. Defaults to 20.
    min_baseq (int, optional): Minimum base quality threshold. Defaults to 20.
    show_progress (bool, optional): Whether to show a progress bar. Defaults to False.
    verbose (int, optional): Verbosity level. Defaults to 0.

Returns:
    str or pandas.DataFrame: Path to the output file if outfile is specified, otherwise a DataFrame containing marked alleles.


**Source:** Line 518 in [phasing.py](phasing.py#L518)


---


### `sort(a: List)`

**Returns:** `None`


Sort intervals based on their coordinates.

Args:
    a: List of intervals to sort. Each interval can be a number or [start, end] pair.


**Source:** Line 25 in [phasing.py](phasing.py#L25)


---


## pileup

**File:** `pileup.py`

**Public functions:** 6


### `TAD_pileup_strength(mat, strip=3)`

**Returns:** `float`


Calculate TAD strength from pileup matrix.
To remove influence from diagonal, we calculate region around domain loop.
The size of region is determined by strip.
Input:
    mat: 2D numpy array
    strip: min distance to diagonal, count by pixels
Output:
    strength: float


**Source:** Line 192 in [pileup.py](pileup.py#L192)


---


### `add_diag_law(raw_mat, binsize=20000, power=0.25)`


Add artificial diagonal law to the matrix.


**Source:** Line 35 in [pileup.py](pileup.py#L35)


---


### `asymmetric_pileup(coolp, refs, expand, expected=None, binsize=None, power=0.25, give_snips=False, balance=False, show_progress=True)`


Fetch any regions from cooler, iterate chrom by chrom.
TODO: support inter-chrom regions
Input:
    coolp: path to cooler file
    refs: list of (chrom1, start1, end1, chrom2, start2, end2)
    expand: treat input refs as point, expand by expand bp
    expected: path to expected file, will give OBS/EXP if provided
    binsize: if provided, will use this to give index and columns for output matrix
    power: power to add to the diagonal, use this to visualize near-diagonal features such as TADs
    give_snips: return the snips used to calculate the pileup
    balance: use balanced matrix
Output:
    mat: pileup matrix


**Source:** Line 112 in [pileup.py](pileup.py#L112)


---


### `block_pileup(coolp, refs, expected=None, power=0.25, give_snips=False, balance=False, debug=False)`


Fetch regions from cooler, iterate chrom by chrom.
Input:
    coolp: path to cooler file
    refs: list of (chrom, start, end)
    expected: expected dataframe, will give OBS/EXP if provided
    power: power to add to the diagonal, use this to visualize TADs
    give_snips: return the snips used to calculate the pileup
    balance: use balanced matrix


**Source:** Line 45 in [pileup.py](pileup.py#L45)


---


### `cool2mat_OE(coolp, chrom, expected, balance=False)`


Fetch cooler matrix and calculate OE.
Input:
    coolp: path to cooler
    chrom: chromosome


**Source:** Line 9 in [pileup.py](pileup.py#L9)


---


### `kernel_strength(kernel, matrix, retmat=False)`


Calculate the signal strength of matrix according to the kernel.

The function computes the signal strength by averaging the values of the kernel
where the mask is True, and dividing by the mean of the values of the False mask positions.

Input:
    kernel (np.ndarray): A 2D numpy array representing the kernel. 0 for signal positions, 1 for background positions.
    Size of kernel must be smaller than the matrix size. Only the region defined by the mask is considered.
    matrix (np.ndarray): A 2D numpy array representing the matrix from which the signal strength is calculated.
retmat (bool): If True, returns the matrix values corresponding to the kernel mask region.

Returns:
float: The signal strength.


**Source:** Line 212 in [pileup.py](pileup.py#L212)


---


## pipeline.rule

**File:** `pipeline/rule.py`

**Public functions:** 8


### `bubble_flow_touched()`


This function returns a dictionary that maps the file paths pattern for
all the files that are touched by the bubble_flow pipeline.


**Source:** Line 116 in [pipeline/rule.py](pipeline/rule.py#L116)


---


### `check_exist(path)`


Check if a file exists.
Tolerate the pd.NA, None, math.nan, np.nan


**Source:** Line 96 in [pipeline/rule.py](pipeline/rule.py#L96)


---


### `check_input(file_pat)`


No documentation available.


**Source:** Line 42 in [pipeline/rule.py](pipeline/rule.py#L42)


---


### `check_key(file_pat, key_to_check)`


Check one key against other keys in the file pattern.
For exmaple, if key_to_check is "sample", pick sample
that has all files that are expanted from file_pat using
other keys in kwargs.
Input:
    file_pat: file pattern
    key_to_check: key to check, one of the keys in kwargs, str
    kwargs: other keys, dict
Output:
    list of values for key_to_check


**Source:** Line 59 in [pipeline/rule.py](pipeline/rule.py#L59)


---


### `json_reader(file)`


No documentation available.


**Source:** Line 270 in [pipeline/rule.py](pipeline/rule.py#L270)


---


### `safe_input(file_pat, with_key=False)`


Expand file_pat to see if the input files exist, if not, skip.
Input:
    file_pat: file pattern
    with_key: return outkeys or not, bool
    kwargs: key word arguments, to expand file_pat with wildcards
Output:
    list of input files or (list of input files, list of outkeys)
    filtered by the existence of the input files


**Source:** Line 9 in [pipeline/rule.py](pipeline/rule.py#L9)


---


### `symlink_files(sample, task_dirp, target_dir, ref, omit=None, filepats=None)`


Create symbolic links from a old task directory to a new task directory.
Use this to partially run the pipeline.
Input:
    sample: the sample name
    task_dirp: the old task directory
    target_dir: the new task directory
    ref: the reference genome, used to create the count matrix file pattern
    omit: the pipeline steps to omit, choose from
        ["DNA", "RNA", "hic_mapped", "seg", "pairs_0", "pairs_c1",
         "pairs_c12", "pairs_c123", "dip", "_3dg", "_3dg_c", "count_matrix"]
        if None, use all steps
    filepats: key-value pairs of file patterns created by the pipeline,
        if None, use the file patterns created by the bubble_flow
        key: the pipeline step name
        value: the file pattern, must contain {sample_name} and {task_dirp} as placeholders
Return:
    None


**Source:** Line 201 in [pipeline/rule.py](pipeline/rule.py#L201)


---


### `unlink_files(sample, target_dir, ref, subset=None)`


No documentation available.


**Source:** Line 256 in [pipeline/rule.py](pipeline/rule.py#L256)


---


## plot.bar

**File:** `plot/bar.py`

**Public functions:** 1


### `plot_bar(df, val_col='value', group_by='group', color_map=None, x_sep=' ')`


Plot bar plot with point and error bar.
Input:
    df: DataFrame with columns [group_by, val_col]
    val_col: column name for the value to plot
    group_by: column name to group by
    color_map: dict mapping group names to colors
    x_sep: separator of groupby names since there might be multiple. e.g. "patient + time"
Output:
    fig: plotly figure object


**Source:** Line 8 in [plot/bar.py](plot/bar.py#L8)


---


## plot.general

**File:** `plot/general.py`

**Public functions:** 8


### `extend(self, trace_names)`


No documentation available.


**Source:** Line 24 in [plot/general.py](plot/general.py#L24)


---


### `get(self, trace_name=None)`


No documentation available.


**Source:** Line 19 in [plot/general.py](plot/general.py#L19)


---


### `plot_colorbar(vmin, vmax, colorscale='Viridis', title='TITLE')`


Plot colorbar.
It is achieved by invisible scatter plot with colorbar.
Input:
    vmin: min value
    vmax: max value
    colorscale: color scale
    title: colorbar title
    colorbar_kwargs: other arguments for coloraxis.colorbar
Output:
    fig


**Source:** Line 74 in [plot/general.py](plot/general.py#L74)


---


### `plot_ols(data, x, y, title)`


Plot linear regression of dataframe


**Source:** Line 281 in [plot/general.py](plot/general.py#L281)


---


### `plot_points(array)`


Wrapping points (represent in numpy column arrays) to dataframe(treate x, y, z as features so in shape N * 3) and plot.
Input:
    array: column array
Output:
    dataframe, x, y, z as columns


**Source:** Line 316 in [plot/general.py](plot/general.py#L316)


---


### `plot_upsetplot(sets, show_counts)`


Upset plot like in R.
Input:
    sets: dict of set
Output:
    upset figure


**Source:** Line 128 in [plot/general.py](plot/general.py#L128)


---


### `scatter_cols(data, points=None, trends=[])`


Multi-traces scatter plot.
Input:
    data: x as index, traces as cols
    trend: plot lowess trends
    point: plot scatter points
Return:
    go.Figure


**Source:** Line 241 in [plot/general.py](plot/general.py#L241)


---


### `upset_plot_getter(data, keys)`


Get aggregate data from upsetplot internal data format (boolen multiindex).
Input:
    data: dict of set
    keys: set names to intersect
Output:
    (count, set)


**Source:** Line 137 in [plot/general.py](plot/general.py#L137)


---


## plot.hic

**File:** `plot/hic.py`

**Public functions:** 18


### `add_eig_track(fig, orig_data, eig_col=None, y_kwargs={}, row=None, col=None)`


Add eigen vector track to a plotly figure or subplot.
Input:
    fig: figure to add trace.
    orig_data: dataframe or series; [chrom, start, eig_col]
    eig_col: column name of eigen vector.
    y_kwargs: kwargs for y axis.
    row: row index of subplot to add trace, set to None if not subplot.
    col: col index of subplot to add trace, set to None if not subplot.
Output:
    fig: figure with eigen track added.


**Source:** Line 638 in [plot/hic.py](plot/hic.py#L638)


---


### `format_ticks(ax, x=True, y=True, rotate=True)`


No documentation available.


**Source:** Line 1272 in [plot/hic.py](plot/hic.py#L1272)


---


### `merge_track_data(clr, track_file, region)`


Pick out the track data that is in the region.
Input:
    clr: cooler obj.
    track_file: path to track file
    region: genome region to plot.
Output:
    pd.Series; index: start of each bin, value: track value.


**Source:** Line 359 in [plot/hic.py](plot/hic.py#L359)


---


### `pileup_IS(coolp, IS, genome='hg19', flank=800000)`


Pileup TAD boundaries from a cool file.
Input:
    coolp: cooler path(adding resolutions when mcool)
    IS: Insulation score; pd.Dataframe
    genome: reference genome. eg. hg19
    flank: flank size
Output:
    pd.Dataframe


**Source:** Line 1244 in [plot/hic.py](plot/hic.py#L1244)


---


### `plot_CM(cm, grid=True)`


Plot contact matrix.
Input:
    cm: pd.DataFrame; contact matrix.
    grid: whether to add chrom boundaries.


**Source:** Line 196 in [plot/hic.py](plot/hic.py#L196)


---


### `plot_IS(IS_file)`


No documentation available.


**Source:** Line 1176 in [plot/hic.py](plot/hic.py#L1176)


---


### `plot_IS(clr, insulation_table, title, resolution=500000.0, window=1500000.0, balance=True, chrom='chr2', start=10500000, steps=90, vmin=0.001, vmax=0.1)`


Plot insulation score together with Hi-C matrix strata.
Input:
    clr: cooler obj.
    insulation_table: insulation score table.


**Source:** Line 1281 in [plot/hic.py](plot/hic.py#L1281)


---


### `plot_cdps(adata, index_col='velocity_pseudotime', reverse=False)`


Input:
    index_col: wich adata.obs col to sort sample by
    reverse: whether to reverse the order of the sample


**Source:** Line 1351 in [plot/hic.py](plot/hic.py#L1351)


---


### `plot_compare_cool_pixels(coolpA, coolpB, region, outline_pixels, subplot_titles=['A', 'B'])`


No documentation available.


**Source:** Line 574 in [plot/hic.py](plot/hic.py#L574)


---


### `plot_compare_cool_track_hor(coolp1, coolp2, IS_file1, IS_file2, eigs_file1, eigs_file2, subplot_titles, region, title, balance=False, winSize=500000)`


Plot cooler with additional track files


**Source:** Line 441 in [plot/hic.py](plot/hic.py#L441)


---


### `plot_compartment(coolp, eigs_file, region, title, eigen_col='E1', strip=False, balance=False, mask_eig_na=True, give_mat=False, quant=0.005, minmax='min', cmap='RdBu_r', eq_hist=False)`


Plot distance-normalized Hi-C correlation matrix with compartment track.
Input:
    coolp: path to cooler file.
    eigs_file: path to eigen vector file.
    region: genome region to plot.
    title: name of the plot.
    eigen_col: which eigen vector to plot.
    strip: whether to strip the consecutive 0s in the matrix.
    balance: whether to load balanced cooler matrix.
    mask_eig_na: whether to mask the heatmap where eigen value is NA.
        Note: you need also to set fillna=False in _plot_mat to make it work.
    give_mat: whether to return the matrix.
    quant: quantile to set the bright extent.
    minmax: whether to use min or max to set the bright extent.
        work with quant.
        Note: this assumes the matrix values are symmetric around 0.
        if eq_hist is True, this will be ignored.
    eq_hist: whether to use equalized histogram.


**Source:** Line 874 in [plot/hic.py](plot/hic.py#L874)


---


### `plot_compartment_strength(adata, index_col='velocity_pseudotime')`


No documentation available.


**Source:** Line 1417 in [plot/hic.py](plot/hic.py#L1417)


---


### `plot_cool(coolp, title='', region='chr1', vmax=100, balance=False, ignore_diags=False, grid=True, donorm=False)`


"
Plot heatmap of single cooler file.
Input:
    coolp: path to cooler file.
    title: name of the plot.
    region: genome region to plot.
        "chr1" or "chr1:1000000-2000000" or ["chr1:1,000,000-2,000,000", "chr2"] or
        slice(0,-1) or [slice(0,1000000), slice(1000000,2000000)]
        Note: can only cross chrom boundaries if using slice.
    vmax: max z value.
    balance: whether to load balanced cooler matrix.
    norm: whether to use lognorm.


**Source:** Line 235 in [plot/hic.py](plot/hic.py#L235)


---


### `plot_cool_track(coolp, track_files, region, title, balance=False)`


Plot cooler with additional track files.
Input:
    coolp: path to cooler file.
    track_files: dict; key: track name, value: (path to track file, column name in track file)
    region: genome region to plot.
    title: name of the plot.
    balance: whether to load balanced cooler matrix.


**Source:** Line 380 in [plot/hic.py](plot/hic.py#L380)


---


### `plot_cools(cools, regions, ncols=3, subplot_titles=None, sub_height=100, sub_width=100, margin=None)`


Plot multiple regions in multiple cooler files.
Iterate over regions first, then coolers. Will skip if either is None.
Input:
    cools: list of cooler file paths.
    regions: list of genome regions to plot.
    ncols: subplot cols.
    sub_height: height of each subplot.
    sub_width: width of each subplot.
    margin: dict; margin of the whole figure, keys: l, r, t, b.
        If None, use default value.
    kwargs: kwargs for _plot_mat.


**Source:** Line 293 in [plot/hic.py](plot/hic.py#L293)


---


### `plot_saddle(file, title, vmin=10 ** (-1), vmax=10 ** 1, strength_square_size=10)`


Plot saddle plot.
Input:
    file: path to saddle file, usually ends with .saddledump.npz.
    title: title of the plot.
    vmin: min z value.
    vmax: max z value.
    strength_square_size: which strength to use, different square size to calculate AA, BB... will give different strength.
    kwargs: kwargs for px.imshow.
Output:
    plot


**Source:** Line 1102 in [plot/hic.py](plot/hic.py#L1102)


---


### `plot_saddle_mpl(file, title, vmin=10 ** (-1), vmax=10 ** 1)`


No documentation available.


**Source:** Line 1088 in [plot/hic.py](plot/hic.py#L1088)


---


### `plot_tiling_compartment(coolps, eigs_files, region, title, corr=True, strip=True, balance=False, cmap='RdBu', donorm=False, Enames=['E1', 'E1'], eig_y_kwargs={})`


Plot 2 coolers in a upper-lower manner with compartment track data.
First cooler will be on the lower triangle, second cooler will be on the upper triangle.
Will ignore diagonal because their values can't be determined.
Input:
    coolps: list of cooler file paths.
    eigs_files: list of eigen vector files.
    region: genome region to plot.
    title: name of the plot.
    corr: whether to plot correlation matrix.
    strip: whether to strip the consecutive 0s in the matrix.
    balance: whether to load balanced cooler matrix.
    cmap: color map.
    donorm: whether to use lognorm.
    Enames: col names of eigen
Output:
    plotly figure.
TODO: rewrite with add_eig_track


**Source:** Line 696 in [plot/hic.py](plot/hic.py#L696)


---


## plot.pca

**File:** `plot/pca.py`

**Public functions:** 2


### `plot_eigenfaces(pca, rows=2, cols=2, fsize=16)`


No documentation available.


**Source:** Line 5 in [plot/pca.py](plot/pca.py#L5)


---


### `plot_elbow(pca)`


No documentation available.


**Source:** Line 25 in [plot/pca.py](plot/pca.py#L25)


---


## plot.plot

**File:** `plot/plot.py`

**Public functions:** 15


### `add_cdps(fig, row, col, adata, sorted_obs)`


No documentation available.


**Source:** Line 124 in [plot/plot.py](plot/plot.py#L124)


---


### `add_compartment_strength(fig, row, col, adata, sorted_obs)`


No documentation available.


**Source:** Line 138 in [plot/plot.py](plot/plot.py#L138)


---


### `add_compartment_strength_config(fig, row, col, adata, sorted_obs)`


No documentation available.


**Source:** Line 255 in [plot/plot.py](plot/plot.py#L255)


---


### `add_contact_number(fig, row, col, adata, sorted_obs)`


No documentation available.


**Source:** Line 170 in [plot/plot.py](plot/plot.py#L170)


---


### `add_gene_heatmap(fig, row, col, adata, sorted_obs, gene_order)`


No documentation available.


**Source:** Line 232 in [plot/plot.py](plot/plot.py#L232)


---


### `add_intra(fig, row, col, adata, sorted_obs)`


No documentation available.


**Source:** Line 186 in [plot/plot.py](plot/plot.py#L186)


---


### `add_pmUMI(fig, row, col, adata, sorted_obs)`


No documentation available.


**Source:** Line 147 in [plot/plot.py](plot/plot.py#L147)


---


### `add_rs(fig, row, col, adata, sorted_obs)`


No documentation available.


**Source:** Line 202 in [plot/plot.py](plot/plot.py#L202)


---


### `add_rs_config(fig, row, col, adata, sorted_obs)`


No documentation available.


**Source:** Line 246 in [plot/plot.py](plot/plot.py#L246)


---


### `add_schicluster_res(fig, row, col, adata, sorted_obs)`


No documentation available.


**Source:** Line 217 in [plot/plot.py](plot/plot.py#L217)


---


### `add_schicluster_res_config(fig, row, col, adata, sorted_obs)`


No documentation available.


**Source:** Line 264 in [plot/plot.py](plot/plot.py#L264)


---


### `cdp_scatter(dat: pd.DataFrame)`

**Returns:** `go.Figure`


No documentation available.


**Source:** Line 15 in [plot/plot.py](plot/plot.py#L15)


---


### `gen_gap(features: list, number: int, k: int)`


Input:
    index: index of features
    m: number of gaps
    k: width of gaps
Output:
    iterator of NA dataframes


**Source:** Line 332 in [plot/plot.py](plot/plot.py#L332)


---


### `plot_cdps_mark(cdps, orig_annote, sample_col='sample_name', order_col='index_order', group_col='group', color_map=None, hline=False)`


# plot cdps heatmap with group marker in different color
# Input:
#    cdps: contact decay profile, row for sample
#    cell_annote: must contain group and ordering col
#    sample_col: name of sample id column
#    order_col: name of order column
#    hline: whether to mark named distance interval
# Output:
#    plotly Figure obj


**Source:** Line 29 in [plot/plot.py](plot/plot.py#L29)


---


### `time_attr(adata, order_col='velocity_pseudotime', using=None, ascending=True, gene_order=None)`


Plot the time attributes of the cells.
Input:
    adata: AnnData object
    order_col: column in adata.obs to sort cells
    using: name of attrs to plot; list; cdps, rs, gene, pmUMI, cs
    ascending: whether to revert order


**Source:** Line 272 in [plot/plot.py](plot/plot.py#L272)


---


## plot.render

**File:** `plot/render.py`

**Public functions:** 27


### `apply_rotation(df, R)`


No documentation available.


**Source:** Line 207 in [plot/render.py](plot/render.py#L207)


---


### `apply_rotation_commands(df, commands_str)`


No documentation available.


**Source:** Line 166 in [plot/render.py](plot/render.py#L166)


---


### `centelo_relpos(positions_2col, genome, dupref=False)`


No documentation available.


**Source:** Line 477 in [plot/render.py](plot/render.py#L477)


---


### `clip_b_pymol(_3dg, b_factor, png, clip=0, slab=2, cmap='magenta green, all, 0.005, 0.02', turn_cmd='', tmpdir=None, conda='pymol')`


Render clip view, color by b factor.
Input:
    _3dg: _3dg file path
    b_factor: a 3 column tsv (chrom pos b_factor) without header
    png: output png file path
    clip: clip position 0 for middle, negative for more back, positive for more front
    cmap: pymol color map
    turn_cmd: PyMOL turn command to adjust camera orientation
    tmpdir: directory to save intermediate files
    conda: conda environment name, None for no conda
    **args: transparent to threedg_to_cif


**Source:** Line 389 in [plot/render.py](plot/render.py#L389)


---


### `clip_centelo_pymol(_3dg_file, png, genome='mm10', clip=0, slab=2, tmpdir=None, cif_name=None, dupref=False, conda='pymol', turn_cmd='')`


Render clip view, color centromere-telomere.
Input:
    _3dg_file: _3dg file path
    png: output png file path
    genome: genome name, used to fetch centromere position
    tmpdir: directory to save intermediate files
    cif_name: cif file path to dump, if None, a random name will be generated and removed after rendering
    dupref: whether to annote diploid genome with haploid reference
    conda: conda environment name, None for no conda
    turn_cmd: PyMOL turn command to adjust camera orientation
    **args: transparent to threedg_to_cif


**Source:** Line 540 in [plot/render.py](plot/render.py#L540)


---


### `clip_single_centelo_pymol(_3dg_file, png, target_chroms=['chrX', 'chrY'], genome='mm10', clip=0, slab=2, tmpdir=None, cif_name=None, dupref=False, conda='pymol', turn_cmd='')`


Render clip view, color centromere-telomere of target chromosomes. Other chromosomes are set to 0.5/white.
Input:
    _3dg_file: _3dg file path
    png: output png file path
    target_chroms: chromosomes to color; list
    genome: genome name, used to fetch centromere position
    clip: clip position 0 for middle, negative for more back, positive for more front
    tmpdir: directory to save intermediate files
    cif_name: cif file path to dump, if None, a random name will be generated and removed after rendering
    dupref: whether to annote diploid genome with haploid reference
    turn_cmd: PyMOL turn command to adjust camera orientation
    **args: transparent to threedg_to_cif


**Source:** Line 614 in [plot/render.py](plot/render.py#L614)


---


### `clip_single_territory_pymol(_3dg_file, png, target_chroms=['chrX', 'chrY'], clip=0, slab=2, tmpdir=None, conda=None, turn_cmd='')`


Render clip view, only color target chromosomes.
Input:
    _3dg_file: _3dg file path
    png: output png file path
    target_chroms: chromosomes to color; list
    clip: clip position 0 for middle, negative for more back, positive for more front
    tmpdir: directory to save intermediate files
    conda: conda environment name, None for no conda
    turn_cmd: PyMOL turn command to adjust camera orientation
    **args: transparent to threedg_to_cif


**Source:** Line 437 in [plot/render.py](plot/render.py#L437)


---


### `clip_territory_pymol(_3dg, png, clip=0, slab=2, tmpdir=None, conda='pymol', turn_cmd='')`


Render clip view, color by each chromosome.
Input:
    _3dg: _3dg file path
    png: output png file path
    clip: clip position 0 for middle, negative for more back, positive for more front
    tmpdir: directory to save intermediate files
    conda: conda environment name, None for no conda
    turn_cmd: PyMOL turn command to adjust camera orientation
    **args: transparent to threedg_to_cif


**Source:** Line 346 in [plot/render.py](plot/render.py#L346)


---


### `dump_chimera_links(links, outfile, rgb: str='255,255,255', radius=0.1, atol=0.001)`


Dump links to Chimera marker file (.cmm file).
Input:
    links: list of tuple of np.ndarray of shape (3,)
        [(corner1, corner2), ...]
    outfile: str, output file name
    rgb: str, color in rgb format
    radius: float, radius of marker
    atol: float, tolerance for checking if two points are close
Output:
    None


**Source:** Line 38 in [plot/render.py](plot/render.py#L38)


---


### `gen_cif(self, _3dg)`


Rewrite the intermediate cif file with the given kwargs.


**Source:** Line 258 in [plot/render.py](plot/render.py#L258)


---


### `gen_script(self)`


Rewrite the intermediate pymol script with the given kwargs.


**Source:** Line 264 in [plot/render.py](plot/render.py#L264)


---


### `get_rotation_commands(target_vector, retvec=False)`


Calculate PyMOL turn commands to rotate the camera to face a given direction.

Input:
    target_vector: A list or numpy array representing the unit vector towards which the camera should be oriented.
    retvec: If True, return the rotation degrees along x, y, z axis as well.
Output:
    A list of string containing a sequence of 'turn' commands for PyMOL to adjust the camera's orientation accordingly.


**Source:** Line 99 in [plot/render.py](plot/render.py#L99)


---


### `get_rotation_matrix(axis, theta)`


No documentation available.


**Source:** Line 189 in [plot/render.py](plot/render.py#L189)


---


### `highlight_surface_b_pymol(_3dg, b_factor, chain, png, cmap='magenta green, chain {}, 0.005, 0.02', turn_cmd='', tmpdir=None, conda='pymol')`


Render surface, color by b factor, highlight only one chain, other chains are transparent.
Input:
    _3dg: _3dg file path
    b_factor: a 3 column tsv (chrom pos b_factor) without header
    chain: chain to highlight
    png: output png file path
    cmap: pymol color map
    turn_cmd: PyMOL turn command to adjust camera orientation
    tmpdir: directory to save intermediate files
    conda: conda environment name, None for no conda
    **args: transparent to threedg_to_cif


**Source:** Line 411 in [plot/render.py](plot/render.py#L411)


---


### `parse_commands(commands_str)`


No documentation available.


**Source:** Line 167 in [plot/render.py](plot/render.py#L167)


---


### `pymol_primary_views(_3dg_file, view='lr', targets=None)`


Rotate a 3dg structure relative to primary axis.
Input:
    _3dg_file: 3dg file of a single-cell structure
    view: view to look at in pymol, ["lr", "ht", "dv"], means looking towards             left-right, head-tail or dorsal-ventral axis.
    targets: flip according to this.
Return:
    _3dg: a rotated 3dg structure suitable for pymol rendering


**Source:** Line 660 in [plot/render.py](plot/render.py#L660)


---


### `render(self)`


No documentation available.


**Source:** Line 281 in [plot/render.py](plot/render.py#L281)


---


### `render_clip_b_primary_view(_3dg_file, b_factor, outpng, view='lr', targets=None)`


Rendering clip primary views of a single-cell structure with pymol, colored by territory.
Designed for a none-glomerulus nucleus.
Input:
    _3dg_file: 3dg file of a single-cell structure
    b_factor: reference b factor
    outpng: output png file
    view: view to render, ["lr", "ht", "dv"], means looking towards             left-right, head-tail or dorsal-ventral axis.
Output:
    outpng path; and wrote outpng file


**Source:** Line 797 in [plot/render.py](plot/render.py#L797)


---


### `render_clip_centelo_primary_view(_3dg_file, outpng, genome, clip, slab, tmpdir, view='lr', targets=None)`


Color primary views of a single-cell structure with relative dist to centromeres and telomeres.
Designed for a none-glomerulus nucleus.
Input:
    _3dg_file: 3dg file of a single-cell structure
    outpng: output png file
    genome: genome name of the structure
    tmpdir: temporary directory
    view: view to render, ["lr", "ht", "dv"], means looking towards             left-right, head-tail or dorsal-ventral axis.
    targets: flip according to this.
Output:
    outpng path; and wrote outpng file


**Source:** Line 818 in [plot/render.py](plot/render.py#L818)


---


### `render_clip_primary_view(_3dg_file, outpng, view='lr', targets=None)`


Rendering clip primary views of a single-cell structure with pymol, colored by territory.
Designed for a none-glomerulus nucleus.
Input:
    _3dg_file: 3dg file of a single-cell structure
    outpng: output png file
    view: view to render, ["lr", "ht", "dv"], means looking towards             left-right, head-tail or dorsal-ventral axis.
Output:
    outpng path; and wrote outpng file


**Source:** Line 778 in [plot/render.py](plot/render.py#L778)


---


### `render_surface_centelo_primary_view(_3dg_file, outpng, genome, tmpdir, view='lr', targets=None)`


Color primary views of a single-cell structure with relative dist to centromeres and telomeres.
Designed for a none-glomerulus nucleus.
Input:
    _3dg_file: 3dg file of a single-cell structure
    outpng: output png file
    genome: genome name of the structure
    tmpdir: temporary directory
    view: view to render, ["lr", "ht", "dv"], means looking towards             left-right, head-tail or dorsal-ventral axis.
    targets: flip according to this.
Output:
    outpng path; and wrote outpng file


**Source:** Line 735 in [plot/render.py](plot/render.py#L735)


---


### `render_surface_primary_view(_3dg_file, outpng, tmpdir, view='lr', targets=None)`


Rendering primary views of a single-cell structure with pymol.
Designed for a none-glomerulus nucleus.
Input:
    _3dg_file: 3dg file of a single-cell structure
    outpng: output png file
    tmpdir: temporary directory
    view: view to render, ["lr", "ht", "dv"], means looking towards             left-right, head-tail or dorsal-ventral axis.
Output:
    outpng path; and wrote outpng file


**Source:** Line 714 in [plot/render.py](plot/render.py#L714)


---


### `search_point(p, points, atol=0.001)`


Iterate over points and check if p is close to any point in points.
Input:
    p: flattened np.ndarray
    points: dict of flattened np.ndarray
    atol: float, tolerance for checking if two points are close
Output:
    If close, return True and index of point.
    If not, return False and None.
NOTE: this is not very efficient. Use KDTree for large number of points.


**Source:** Line 22 in [plot/render.py](plot/render.py#L22)


---


### `single_centelo_relpos(positions_2col, genome, target_chroms, dupref=False)`


No documentation available.


**Source:** Line 577 in [plot/render.py](plot/render.py#L577)


---


### `surface_b_pymol(_3dg, b_factor, png, cmap='magenta green, all, 0.005, 0.02', tmpdir=None, conda='pymol', turn_cmd='')`


Render surface, color by b factor.
Input:
    _3dg: _3dg file path
    b_factor: a 3 column tsv (chrom pos b_factor) without header
    png: output png file path
    cmap: pymol color map
    tmpdir: directory to save intermediate files
    conda: conda environment name, None for no conda
    turn_cmd: PyMOL turn command to adjust camera orientation
    **args: transparent to threedg_to_cif


**Source:** Line 369 in [plot/render.py](plot/render.py#L369)


---


### `surface_centelo_pymol(_3dg_file, png, genome='mm10', tmpdir=None, cif_name=None, dupref=False, conda='pymol', turn_cmd='')`


Render clip view, color centromere-telomere.
Input:
    _3dg_file: _3dg file path
    png: output png file path
    genome: genome name, used to fetch centromere position
    tmpdir: directory to save intermediate files
    cif_name: cif file path to dump, if None, a random name will be generated and removed after rendering
    dupref: whether to annote diploid genome with haploid reference
    turn_cmd: PyMOL turn command to adjust camera orientation
    **args: transparent to threedg_to_cif


**Source:** Line 506 in [plot/render.py](plot/render.py#L506)


---


### `surface_territory_pymol(_3dg, png, tmpdir=None, conda='pymol', dcmap='dip_c', turn_cmd='', cif_name=None)`


Render surface, color by each chromosome.
Input:
    _3dg: _3dg file path
    png: output png file path
    dcmap: discrete color map.
        dip_c: each chromosome has a different rainbow color
    tmpdir: directory to save intermediate files
    conda: conda environment name, None for no conda
    turn_cmd: PyMOL turn command to adjust camera orientation


**Source:** Line 288 in [plot/render.py](plot/render.py#L288)


---


## plot.rna

**File:** `plot/rna.py`

**Public functions:** 3


### `plot_gene_module_mean(adata, gene_sets, ps_col, layer='not_known', trends=False, points=True)`


Plot mean expression of module along trajectory. 
Input:
    adata: anndata object, or dataframe
    gene_sets: dict, module names as keys, list of module gene names as values
    ps_col: pseudotime column name, must be in adata.obs
    layer: layer name of count matrix
Return:
    plotly figure


**Source:** Line 57 in [plot/rna.py](plot/rna.py#L57)


---


### `plot_gene_trend(adata, gene, order_col='velocity_pseudotime', additional=None)`


Ploting expression trends along time.
Must have pSmooth in uns; normalized and log1p X counts


**Source:** Line 69 in [plot/rna.py](plot/rna.py#L69)


---


### `plot_gene_trends(data, additional, genes, order_col, ncols)`


No documentation available.


**Source:** Line 121 in [plot/rna.py](plot/rna.py#L121)


---


## plot.scanpy

**File:** `plot/scanpy.py`

**Public functions:** 5


### `plot_diffmap(adata, color, title)`


No documentation available.


**Source:** Line 193 in [plot/scanpy.py](plot/scanpy.py#L193)


---


### `plot_elbow(adata)`


No documentation available.


**Source:** Line 150 in [plot/scanpy.py](plot/scanpy.py#L150)


---


### `plot_pca(adata, color: str, PCs: Optional[Union[List[int], List[List[int]]]]=None, ncols: int=2)`


Plot PCA results from AnnData object with flexible subplot configurations.

This function creates scatter plots of PCA components from single-cell data.
It can generate single plots or multiple subplots showing different combinations
of principal components.

Args:
    adata: AnnData object containing PCA results in obsm['X_pca']
    color: Column name in adata.obs to use for coloring points
    PCs: Specification of which principal components to plot. Can be:
        - None: plots first two PCs (PC1 vs PC2)
        - List[int]: plots all pairwise combinations of the specified PCs
        - List[List[int]]: plots specific pairs of PCs as subplots
    ncols: Number of columns for subplot layout when multiple plots are generated
    **kwargs: Additional arguments passed to px.scatter
    
Returns:
    plotly.graph_objects.Figure: Figure object containing the PCA plot(s)
    
Examples:
    >>> # Plot first two PCs
    >>> fig = plot_pca(adata, color='cell_type')
    
    >>> # Plot all pairwise combinations of PCs 1-3
    >>> fig = plot_pca(adata, color='cell_type', PCs=[0, 1, 2])
    
    >>> # Plot specific PC pairs in subplots
    >>> fig = plot_pca(adata, color='cell_type', PCs=[[0, 1], [0, 2], [1, 3]])
    
    >>> # Customize subplot layout
    >>> fig = plot_pca(adata, color='cell_type', PCs=[0, 1, 2, 3], ncols=3)


**Source:** Line 14 in [plot/scanpy.py](plot/scanpy.py#L14)


---


### `plot_tsne(adata, color)`


Plot scanpy tsne using plotly.


**Source:** Line 177 in [plot/scanpy.py](plot/scanpy.py#L177)


---


### `plot_umap(adata, color, key='X_umap')`


Plot scanpy umap using plotly.


**Source:** Line 159 in [plot/scanpy.py](plot/scanpy.py#L159)


---


## plot.scatter

**File:** `plot/scatter.py`

**Public functions:** 6


### `add_error_band(fig, df, orig_color=px.colors.qualitative.Plotly[0], debug=False)`


Input:
    df: must have x, upper, lower cols
    orig_color: color of original line, will use a grayer version


**Source:** Line 93 in [plot/scatter.py](plot/scatter.py#L93)


---


### `build_hover_template(group_name=None)`


Build hover template based on available data.


**Source:** Line 229 in [plot/scatter.py](plot/scatter.py#L229)


---


### `line(df, x, y, color=None, color_discrete_map=None, with_error=False, frac=0.3, ci=0.95, debug=False, gray_bg=False, hover_data=None, hover_format=None)`


Plot a single trend line with optional error bands and point coloring.

Creates a line plot for a single variable with optional confidence intervals.
Points can be colored according to a specified column in the dataframe.

Args:
    df (pd.DataFrame): Input data.
    x (str): Column name for x-axis values.
    y (str): Column name for y-axis values.
    color (str, optional): Column name for point coloring. Points will be 
                          colored according to this column's values. 
                          Defaults to None (uniform color).
    color_discrete_map (dict, optional): Color mapping for points when 
                                        `color` is specified. Keys should be 
                                        values from the `color` column.
                                        Defaults to None.
    with_error (bool, optional): Whether to add confidence interval bands.
                                Defaults to False.
    frac (float, optional): Fraction of data to use for smoothing (LOESS).
                            Defaults to 0.3.
    ci (float, optional): Confidence interval level (0-1). Defaults to 0.95.
    debug (bool, optional): Enable debug mode for error bands. 
                            Defaults to False.
    gray_bg (bool, optional): Use gray color for error bands. Defaults to False.
    **kwargs: Additional keyword arguments passed to fig.update_layout().

Returns:
    plotly.graph_objs.Figure: Plotly figure object.

Examples:
    >>> line(df, x='time', y='value', color='group')
    >>> line(df, x='x', y='y', with_error=True, color_discrete_map={'A': 'red', 'B': 'blue'})


**Source:** Line 179 in [plot/scatter.py](plot/scatter.py#L179)


---


### `lines(df, x, y, color_discrete_map=None, with_error=False, frac=0.3, ci=0.95, debug=False, gray_bg=False)`


No documentation available.


**Source:** Line 125 in [plot/scatter.py](plot/scatter.py#L125)


---


### `rgb_grayer(color, white_ratio=0.7, opacity=0.3)`


Convert a color to a lighter, more transparent version by blending with white
and reducing opacity.

Parameters:
    color: A list of three integers representing RGB values (0-255).

Returns:
    A string in "rgba(r, g, b, a)" format with adjusted color and opacity.


**Source:** Line 69 in [plot/scatter.py](plot/scatter.py#L69)


---


### `smooth_with_ci(df, x, y, frac=0.3, ci=0.95)`


Apply LOESS smoothing to time-series data and estimate confidence intervals.

Parameters
----------
x : array-like
    Independent variable (e.g., time points).
y : array-like
    Dependent variable (e.g., gene expression values).
frac : float, optional
    Fraction of data used for each local regression (0 < frac <= 1).
ci : float, optional
    Confidence interval level (0 < ci < 1, e.g., 0.95 for 95% CI).

Returns
-------
pandas.DataFrame
    DataFrame with columns:
    - 'x': Original x values
    - 'y': Original y values
    - 'y_smooth': Smoothed y values
    - 'y_upper': Upper bound of CI
    - 'y_lower': Lower bound of CI


**Source:** Line 11 in [plot/scatter.py](plot/scatter.py#L11)


---


## plot.tech

**File:** `plot/tech.py`

**Public functions:** 1


### `tech_box(annote, grouping='condition', features=['defaults'], ncols=3, size=200, color_discrete_sequence=color_discrete_sequence, use_label=True, labels: dict=None)`

**Returns:** `go.Figure`


DNA and RNA tech metrics.
Input:
    annote: sample * feature dataframe
    grouping: col name to group by
        use all if grouping is None
    features: features to be plot. ["defaults", ...]
    ncols: number of cols of subplots
    size: subplot size
    color_discrete_sequence: cycling pallete for each group
    use_label: use a more formal name for each feature
        formal names are defined in wet.afbb.formal_feature_names
    labels: customized formal feature names
TODO:
    set group order


**Source:** Line 9 in [plot/tech.py](plot/tech.py#L9)


---


## plot.track

**File:** `plot/track.py`

**Public functions:** 4


### `add_chrom_grid(df, fig, fig_type='scatter', x_extra_expansion=0, y_extra_expansion=0.1, value_cols=None, range_y=None)`


Add chromosome grid lines to an existing genome track figure.
Input:
    df: pd.DataFrame, the genome data used to plot the figure. Index is multiindex with level 0 as chromosome names 
        and level 1 as start position. First column is the data to plot.
        1. For scatter plot, use first value column to determine y axis title.
        2. For heatmap, df is `pos * samples`, use "Sample" as y axis title, each column is a sample.
    fig: plotly.graph_objs._figure.Figure, the figure to add grid lines to.
    fig_type: str, default "scatter", type of the figure, currently only support "scatter", "heatmap".
    x_extra_expansion: float, default 0, the extra expansion (relative to min, max) of x axis.
    y_extra_expansion: float, default 0.1, the extra expansion (relative to min, max) of y axis.
        For heatmap, y axis is sample index, so the best is y_extra_expansion=0.
    range_y: tuple or list, default None, the y axis range.
Output:
    fig: plotly.graph_objs._figure.Figure, the figure with added grid lines.


**Source:** Line 111 in [plot/track.py](plot/track.py#L111)


---


### `plot_genome_heatmap(genome_data, chromosomes=None)`


Heatmap genome wide data. Sepaarate chromosomes by vertical lines and add chromosome names.
Input:
    genome_data: pd.DataFrame, index is multiindex with level 0 as chromosome names
        columns are samples. Will sort the index first.
    chromosomes: pd.DataFrame or Feature, the chromosome information used to sort the genome_data.
Output:
    fig: plotly.graph_objs._figure.Figure


**Source:** Line 77 in [plot/track.py](plot/track.py#L77)


---


### `plot_genome_scatter(genome_data, chromosomes=None, value_cols=None, shift_unit=0.05, shift_threshold=100, chorm_label_shift=0.05, y_extra_expansion=0.1, range_y=None)`


Scatter genome wide data. Sepaarate chromosomes by vertical lines and add chromosome names.
TODO: add real positions, now only add 1,2,3,4... as x axis.
Input:
    genome_data: pd.DataFrame, index is multiindex with level 0 as chromosome names 
        and level 1 as start position. First column is the data to plot.
        Will sort the index first. To plot chromosome
        according to your desired order, set the index level 0 as ordered category.
    shift_unit: float, default 0.05, the unit of shift when two chromosomes are too close.
    shift_threshold: int, default 100, the threshold to determine if two chromosomes are too close.
    chorm_label_shift: float, default 0.05, the shift of chromosome name from the top of the plot.
    y_extra_expansion: float, default 0.1, the extra expansion (relative to min, max) of y axis.
    range_y: tuple or list, default None, the y axis range.
Output:
    fig: plotly.graph_objs._figure.Figure


**Source:** Line 26 in [plot/track.py](plot/track.py#L26)


---


### `sort_index_chrom(df, chromosomes)`


No documentation available.


**Source:** Line 8 in [plot/track.py](plot/track.py#L8)


---


## plot.utils

**File:** `plot/utils.py`

**Public functions:** 20


### `add_cat_marker(fig, data, catcol='cell_type', ypos=1)`


Adding colored scatter trace to figure to mark x categories.
Input:
    fig: plotly figure
    data: dataframe storing category info, must have same x index with fig data
    catcol: category column name, must be in adata.obs
    ypos: y position of marker points
Output:
    fig with new traces, each cat one trace


**Source:** Line 160 in [plot/utils.py](plot/utils.py#L160)


---


### `clip_axis_range(fig, x_lower=True, x_upper=True, y_lower=False, y_upper=False, main_fig=None, row=None, col=None)`


Clip plotly axis range to exactly match data.

x- and y-axis, lower and upper bound are both free to choose.
Iterates through data traces contained in figure, extracts min and max, and sets the respective axes range bounds to these.

Parameters
----------
fig: plotly.graph_objs._figure.Figure
    The figure to update.
x_lower: bool
    Clip x axis lower range. Defautls to True.
x_upper: bool
    Clip x axis upper range. Defautls to True.
y_lower: bool
    Clip y axis lower range. Defautls to False.
y_upper: bool
    Clip y axis upper range. Defautls to False.
main_fig:
    If not None, change main_fig axes according to fig traces.
    Use this when you are make subplots.
row:
    If main_fig is not None, specify the row of the subplot where you want to change axes.
col:
    If main_fig is not None, specify the column of the subplot where you want to change axes.

Returns
-------
plotly.graph_objs._figure.Figure
    Plotly graph with updated axis ranges.

FROM: https://community.plotly.com/t/set-axis-range-to-match-data/93168/3


**Source:** Line 15 in [plot/utils.py](plot/utils.py#L15)


---


### `compute_hue_diff(h1_deg, h2_deg)`


计算两个色相之间的最短差值（考虑环形特性）


**Source:** Line 397 in [plot/utils.py](plot/utils.py#L397)


---


### `expand_color_sequence(colors, n, cap=0.75)`


Expand color sequence by create n interpolated colors for each original color.
The created colors are "whiter" version of the original colors.
Group the original colors (orig_colors) in pairs, and within each pair,
arrange the interpolated colors in opposite directions.
For example: blue → light blue → light red → red.
Input:
    colors: list of original colors
    n: number of interpolated colors for each original color
Output:
    list of expanded colors


**Source:** Line 439 in [plot/utils.py](plot/utils.py#L439)


---


### `expand_colors(colors, n_interpolate=1)`


在相邻颜色之间插入 n_interpolate 个新颜色，使颜色过渡自然。

参数:
    colors (List[tuple]): 输入颜色序列，每个颜色是 (r, g, b) 的 RGB 元组，范围 0~1。
    n_interpolate (int): 每对颜色之间插入的新颜色数量。

返回:
    List[tuple]: 扩展后的颜色序列，每个颜色为 (r, g, b) 的整数元组，范围 0~255。


**Source:** Line 470 in [plot/utils.py](plot/utils.py#L470)


---


### `filling_l2r_mpl(rows, cols, features)`


Helper to iterate within row-cols. 
Mpl_flavor, starts with 0.


**Source:** Line 130 in [plot/utils.py](plot/utils.py#L130)


---


### `filling_l2r_plotly(rows, cols, features)`


Helper to iterate within row-cols.
Plotly_flavor, starts with 1.
Return:
    row, col, index, feature


**Source:** Line 114 in [plot/utils.py](plot/utils.py#L114)


---


### `get_fig_outprefix()`

**Returns:** `str`


Get the output prefix for the current figure.
Change directory name and this will get the right prefix.
Fig1a -> output/Fig.1a_
Fig1Sb -> output/Extended_Data_Fig.1b_


**Source:** Line 550 in [plot/utils.py](plot/utils.py#L550)


---


### `hex_split(hex_color)`


解析十六进制颜色代码，并返回 R、G、B 和 A 分量。

:param hex_color: 十六进制颜色代码，可以是六位 (#RRGGBB) 或八位 (#RRGGBBAA)
:return: 一个元组 (R, G, B, A)，其中 A 是可选的透明度分量


**Source:** Line 369 in [plot/utils.py](plot/utils.py#L369)


---


### `hex_to_rgb(hex_color)`


将十六进制颜色字符串转换为 RGB 格式。

参数:
    hex_color (str): 十六进制颜色字符串，例如 "#4C78A8"

返回:
    tuple: 三元组 (r, g, b)，其中 r, g, b 分别为红、绿、蓝通道的值 (0-255)


**Source:** Line 335 in [plot/utils.py](plot/utils.py#L335)


---


### `interpolate_color_rgb_linear(orig, target, n, cap=1)`


Interpolate between two colors in RGB space.
Input:
    orig: original color (r, g, b)
    target: target color (r, g, b)
    n: number of colors to generate
    cap: 0<cap<1, ensure the last interpolated color is not too close to the target color
        near 0, almost same as orig
        near 1, interpolated color cover the whole range from orig to target
Output:
    list of n interpolated colors


**Source:** Line 408 in [plot/utils.py](plot/utils.py#L408)


---


### `list2colorlist(celltypes)`


Tranform a data list to a color list, ready for all color argument.
TODO:
    custom mapper;
    custom color sequence


**Source:** Line 143 in [plot/utils.py](plot/utils.py#L143)


---


### `mat_coarsen(mat, coarseness)`


# coarsen matrix to lower resolution
# Input:
#    mat: matrix
# Output:
#    matrix


**Source:** Line 257 in [plot/utils.py](plot/utils.py#L257)


---


### `pcolormesh_45deg(ax, matrix_c, start=0, resolution=1)`


Helper to plot a matrix with 45 degree angle.


**Source:** Line 241 in [plot/utils.py](plot/utils.py#L241)


---


### `plot_color_sequence(colors, output_file=None)`


Plot a color sequence as a 1-row heatmap.

Args:
    colors (List[str] or List[tuple]): List of colors in hex or RGB format.
    output_file (str, optional): Path to save the plot as an image. If None, the plot is shown.


**Source:** Line 299 in [plot/utils.py](plot/utils.py#L299)


---


### `plotly_fig2array(fig)`


No documentation available.


**Source:** Line 540 in [plot/utils.py](plot/utils.py#L540)


---


### `rgb_to_hex(rgb_color)`


将 RGB 格式转换为十六进制颜色字符串。

参数:
    r (int): 红色通道的值 (0-255)
    g (int): 绿色通道的值 (0-255)
    b (int): 蓝色通道的值 (0-255)

返回:
    str: 十六进制颜色字符串，例如 "#4C78A8"


**Source:** Line 349 in [plot/utils.py](plot/utils.py#L349)


---


### `spread_text(text_list, track_width=5, fold=2)`


Spread text list to a several tracks.
Input:
    text_list: list of text to spread
    track_width: width of each track
    fold: number of tracks
Output:
    list of position shift for each text
TODO:
    add random shift


**Source:** Line 279 in [plot/utils.py](plot/utils.py#L279)


---


### `sub_genome_mat(mat, keep_regions)`


Keep only the regions in the keep_regions list.
Input:
    mat: index and columns are [chrom, start]
    keep_regions: list of regions to keep
        regions are tuples of (start, end)
Output:
    mat: matrix with only the regions in keep_regions


**Source:** Line 216 in [plot/utils.py](plot/utils.py#L216)


---


### `tiling_mat(A_orig, Ref_orig, adjust=True, ignore_diags=True)`


Tiling 2 matrix. Rising light to the lighter one.
Input:
    A: matrix to show on upper right.
    Ref: matrix to show on lower left.
    adjust: adjust the color scale to the stronger one. Note: maybe overflow if brightness of the two matrix are too different.
Output:
    A tiled matrix.


**Source:** Line 187 in [plot/utils.py](plot/utils.py#L187)


---


## pp

**File:** `pp.py`

**Public functions:** 1


### `standard_scaler(df, axis=0, with_std=False)`


Return scaled dataframe.
Input:
    df: pandas dataframe
    axis: 0 or 1. 0 to scale each col(col as feature), 1 to scale each row(row as feature).
    with_std: whether doing std scaling. No consensus but 
        usually don't scale it in omic-biology.
Output:
    scaled dataframe.


**Source:** Line 77 in [pp.py](pp.py#L77)


---


## pseudotime.TI

**File:** `pseudotime/TI.py`

**Public functions:** 1


### `angle_pseudotime(root, embedding, direction=1, colname='pseudotime')`

**Returns:** `pd.DataFrame`


Transform 2D embedding to pseudotime by naive circular ordering.
Input:
    root: set what sample/angle to 0.
        str: name of sample, must be in data.index; float: angle value
    df: 2d embedding of samples; dataframe
    direction: 1 for clockwise, -1 for anticlockwise
Output:
    0 - 1 pseudotime, dataframe with new column


**Source:** Line 33 in [pseudotime/TI.py](pseudotime/TI.py#L33)


---


## pymol.color_centelo

**File:** `pymol/color_centelo.py`

**Public functions:** 2


### `color_centelo(obj, genome)`


No documentation available.


**Source:** Line 71 in [pymol/color_centelo.py](pymol/color_centelo.py#L71)


---


### `fetch_cent_chromlen(genome)`


Get centromere regions and chromosome lengths.
Only work with mm10 for now.
TODO: generalize to other genomes
Input:
    fp (str): path to the gap file
Returns:
    dict(chrom) of list[start, end, chrom_length]


**Source:** Line 17 in [pymol/color_centelo.py](pymol/color_centelo.py#L17)


---


## pymol.glow

**File:** `pymol/glow.py`

**Public functions:** 1


### `glow(colors, selections)`


No documentation available.


**Source:** Line 15 in [pymol/glow.py](pymol/glow.py#L15)


---


## pymol.sele_cent

**File:** `pymol/sele_cent.py`

**Public functions:** 2


### `fetch_centromeres(genome)`


Get centromere regions.
Only work with mm10 for now.
TODO: generalize to other genomes
Input:
    fp (str): path to the gap file
Returns:
    dict(chrom) of list[start, end]


**Source:** Line 27 in [pymol/sele_cent.py](pymol/sele_cent.py#L27)


---


### `sele_cent(obj, genome, flank=1)`


No documentation available.


**Source:** Line 62 in [pymol/sele_cent.py](pymol/sele_cent.py#L62)


---


## pymol.sele_telo

**File:** `pymol/sele_telo.py`

**Public functions:** 2


### `fetch_telomeres(genome)`


Get centromere regions.
Returns:
    dict(chrom) of list[left centromere, right centromere] of list[start, end]


**Source:** Line 27 in [pymol/sele_telo.py](pymol/sele_telo.py#L27)


---


### `sele_telo(obj, genome, flank=1)`


No documentation available.


**Source:** Line 61 in [pymol/sele_telo.py](pymol/sele_telo.py#L61)


---


## pymol.utils

**File:** `pymol/utils.py`

**Public functions:** 2


### `fetch_centromeres(genome)`


Get centromere regions.
Only work with mm10 for now.
TODO: generalize to other genomes
Input:
    fp (str): path to the gap file
Returns:
    dict(chrom) of list[start, end]


**Source:** Line 18 in [pymol/utils.py](pymol/utils.py#L18)


---


### `get_ref_dir()`


Return a path to src's ref


**Source:** Line 11 in [pymol/utils.py](pymol/utils.py#L11)


---


## scAB_embedding

**File:** `scAB_embedding.py`

**Public functions:** 9


### `CpG_stat(acc, x)`


Add value to CpG sum and neighbor count.
Input:
    acc: original stat, 2-tuple, (CpG_sum, n_neighbors)
    x: input_line, (index, neighbor CpG value)
Ouput:
    new_acc: updated stat
Note: if nan, will be ignored


**Source:** Line 186 in [scAB_embedding.py](scAB_embedding.py#L186)


---


### `calc_color2(filesp, file_col, color_file, binsize, rank_norm=True, merge_haplotypes=True, dropXY=True, dupref=False, col_thresh=0.9, row_thresh=0.9, color_fill=True, threads=24)`


Input:
    filesp: dataframe, must have file_col col
    file_col: column in filesp that stores pairs file path
    color_file: ref bed file that stores linear CpG density
    binsize: binsize of color_file
    merge_haplotypes: whether treat maternal paternal chromosome differently,
        when False, chromosomes must be like "chr1(mat)"
    dropXY: only output autosome result
    col_thresh: [0,1] larger is stricter, 0 to keep all
        keeping cols that at leat *ratio* of samples have nonNA value
    row_thresh: [0,1] larger is stricter, 0 to keep all
        after col cleaning, keeping samples that have nonNA value for at leat *ratio* of all attrs
    threads: number of cores using
Note:
    when binsize is small, valid col and row thresh drop dramatically, it's better to use color3 at that time
TODO: use 2-layer column


**Source:** Line 267 in [scAB_embedding.py](scAB_embedding.py#L267)


---


### `color2(pairsf, color_file, bin_size, merge_haplotypes=True, dupref=False, dropXY=True, custom_color_file=False)`


Calculate 3D CpG density of each chromosome bin
Input:
    pairsf: 4DN's pairs file
    color_file: reference linear CpG density
    bin_size: binsize of color_file(and the result output)
    merge_haplotypes: whether treat maternal paternal chromosome differently,
        when False, chromosomes must be like "chr1(mat)"
    dupref: whether to duplicate reference CpG density for both homologous chromosomes
        In practice, this is implemented by omitting the suffix in homologous chromosome name when querying color_file
    dropXY: only output autosome result


**Source:** Line 22 in [scAB_embedding.py](scAB_embedding.py#L22)


---


### `do_umap(data, ndims=30)`


No documentation available.


**Source:** Line 308 in [scAB_embedding.py](scAB_embedding.py#L308)


---


### `fill_color(data, color_file, grt_full=True)`


# Fill in NA using reference CpG value
# This make sense because spacing CpG is highly correlated with Sequencing CpG
# Input:
#    grt_full: guarantee that output dataframe is NA free,
#        this will drop all col's that can't be filled by ref color_file


**Source:** Line 244 in [scAB_embedding.py](scAB_embedding.py#L244)


---


### `get_bin_locus(pos, size)`


No documentation available.


**Source:** Line 20 in [scAB_embedding.py](scAB_embedding.py#L20)


---


### `s_color2(_3dg: pd.DataFrame, color_file, min_dist=3, rank_norm=False, n_jobs=12, dupref=True)`

**Returns:** `pd.DataFrame`


Calculate intermingling metrics for each particle in 3dg file.
Input:
    _3dg: dataframe, result of hic_basic.hicio parse_3dg
    color_file: reference linear CpG density
    min_dist: int, minimum distance to be considered as neighbor
    rank_norm: bool, whether to rank normalize output
    n_jobs: int, number of jobs to run in parallel
Output:
    dataframe with same index, and a scAB column


**Source:** Line 112 in [scAB_embedding.py](scAB_embedding.py#L112)


---


### `scAB_embedding(data, ndims=30, dorank=False)`


No documentation available.


**Source:** Line 337 in [scAB_embedding.py](scAB_embedding.py#L337)


---


### `stack_dict(ares, sample_name: None, col_thresh=0.9, row_thresh=0.9)`


# stack list of dict, drop bad rows and columns
# first clean bad cols, then clean bad rows
# Input:
#    ares: list of dict
#    sample_name: name list of samples
#    col_thresh: [0,1] larger is stricter, 0 to keep all
#        keeping cols that at leat *ratio* of samples have nonNA value
#    row_thresh: [0,1] larger is stricter, 0 to keep all
#        after col cleaning, keeping samples that have nonNA value for at leat *ratio* of all attrs


**Source:** Line 226 in [scAB_embedding.py](scAB_embedding.py#L226)


---


## scripts.downsra

**File:** `scripts/downsra.py`

**Public functions:** 1


### `newname(url)`


No documentation available.


**Source:** Line 23 in [scripts/downsra.py](scripts/downsra.py#L23)


---


## scripts.rescue_jupyter

**File:** `scripts/rescue_jupyter.py`

**Public functions:** 2


### `find_code(text)`


No documentation available.


**Source:** Line 8 in [scripts/rescue_jupyter.py](scripts/rescue_jupyter.py#L8)


---


### `print_code(source)`


No documentation available.


**Source:** Line 20 in [scripts/rescue_jupyter.py](scripts/rescue_jupyter.py#L20)


---


## sequence

**File:** `sequence.py`

**Public functions:** 11


### `add_suffix_to_fastq(input_file: str, output_file: str, skip_incomplete=True)`

**Returns:** `None`


Process a FASTQ file by appending "/1" or "/2" (MGI/DNBSEQ flavor) to read IDs if not already present (illumina_flavor).

For Illumina reads, the prefix of the read ID is determined according to the comment part (e.g. `/1` for `1:N:0:1`).

Parameters:
input_file (str): Path to the input FASTQ file (can be gzip-compressed if ending with .gz).
output_file (str): Path to the output FASTQ file (will be gzip-compressed if input is compressed).

Returns:
None


**Source:** Line 342 in [sequence.py](sequence.py#L342)


---


### `compare_fastq_pairs(fq1_path, fq2_path, showstats=False)`


Compare paired-end FASTQ files and return statistics.

Args:
    fq1_path (str): Path to FASTQ file 1 (R1)
    fq2_path (str): Path to FASTQ file 2 (R2)
    
Returns:
    dict: Statistics including:
        'total_fq1' (int): Total reads in FQ1
        'total_fq2' (int): Total reads in FQ2
        'valid_pairs' (int): Valid paired reads
        'invalid_pairs' (int): Reads with mismatched ends
        'missing_in_fq1' (int): Reads in FQ2 missing from FQ1
        'missing_in_fq2' (int): Reads in FQ1 missing from FQ2


**Source:** Line 89 in [sequence.py](sequence.py#L89)


---


### `count_CpG(bed_df: pd.DataFrame, fasta: str)`

**Returns:** `pd.DataFrame`


Calculate CpG ratio for each bin in bed_df.
Input:
    bed_df: pd.DataFrame, bed file with columns chrom, start, end
    fasta: str, path to the fasta file
Output:
    pd.DataFrame, bed file with additional column CpG_ratio
Note:
    1. omit all Ns in the sequence
    2. theoretically, max CpG ratio is 0.5 here.
    3. only bins with CpG ratio >= 0 are returned


**Source:** Line 147 in [sequence.py](sequence.py#L147)


---


### `count_cg(seq)`


No documentation available.


**Source:** Line 166 in [sequence.py](sequence.py#L166)


---


### `count_fastq_worker(file_path, form='s')`


No documentation available.


**Source:** Line 13 in [sequence.py](sequence.py#L13)


---


### `extract_first_n_reads(input_file: str, output_file: str, n: int)`

**Returns:** `None`


Extract the first N reads from a FASTQ file.

Reads from the input FASTQ file and writes the first N reads to the output file.
The compression format of the output matches the input (compressed if input ends with .gz).

Parameters:
input_file (str): Path to the input FASTQ file (can be gzip-compressed if ending with .gz).
output_file (str): Path to the output FASTQ file.
n (int): Number of reads to extract.

Returns:
None


**Source:** Line 425 in [sequence.py](sequence.py#L425)


---


### `parse_fastq_id(line)`


Parse a FASTQ ID line into (prefix, end) tuple or None if invalid.

Args:
    line (str): The FASTQ ID line (e.g., '@M06168:... 1:N:0:1' or '@E250109998L1C001R00300005452/1')
    
Returns:
    tuple or None: (prefix, end) if valid, else None.
        prefix: The read identifier prefix (e.g., '@M06168:...' or '@E250109998L1C001R00300005452').
        end: The read end identifier (e.g., '1' or '2').

Notes
-----
Supports two formats:
1. Illumina format:
    `@<instrument>:<run_id>:<flowcell_id>:<lane>:<tile>:<x_pos>:<y_pos> <read>:<is_filtered>:<control_number>:<index_seq>`
    Example: `@M06168:86:000000000-LMWTD:1:1101:14802:1380 1:N:0:1`

2. DNBSEQ format:
    `@<unique_identifier>/<read_end>`
    Example: `@E250109998L1C001R00300005452/1`


**Source:** Line 20 in [sequence.py](sequence.py#L20)


---


### `plot_search_primers(result, output_file=None)`


Visualize the primer search result using bar plots with 4 subplots.

Args:
    result (pd.DataFrame): Primer search result with columns for each form.
    output_file (str, optional): Path to save the plot as an HTML file. If None, the plot is shown.


**Source:** Line 292 in [sequence.py](sequence.py#L292)


---


### `print_stats(stats)`


Print comparison statistics in a human-readable format.


**Source:** Line 138 in [sequence.py](sequence.py#L138)


---


### `read_fastq_ids(file_path)`


Read FASTQ file (supports gzipped) and return a dictionary of {prefix: end}.

This function correctly identifies ID lines by tracking the 4-line FASTQ structure.

Args:
    file_path (str): Path to the FASTQ file (can be gzipped)
    
Returns:
    dict: {prefix: end} mapping for all valid ID lines


**Source:** Line 61 in [sequence.py](sequence.py#L61)


---


### `search_primers(fq_fp, primers, sample_n=None, sample_r=None, seed=None, keep_hits=0)`


Search primers in fastq file and return distribution of primer start positions.
Will search forward, reverse, complement and reverse complement.

Input:
    fq_fp: str, fastq file path
    primers: list of str, primers to search
        If list, search all primers as one.
    sample_n: int, optional, number of reads to sample
    sample_r: float, optional, fraction of reads to sample (0 < sample_r <= 1)
    seed: int, optional, random seed for reproducibility
    keep_hits: keep number of hits for further analysis, default 0

Output:
    if keep_hits == 0:
        result (pd.DataFrame): Primer search result with columns for each form.
    if keep_hits > 0:
        result,
        hits (dict of list): sampled hits for each primer form.


**Source:** Line 181 in [sequence.py](sequence.py#L181)


---


## shuffle

**File:** `shuffle.py`

**Public functions:** 4


### `AB_mix(obsA: list, obsB: list, size_of_each_group=100, fraction=0.5, N=1000, seed=None)`


Shuffle two groups. Will first random sample each group N times, then exchange fraction of elements for each pair.
Both group-pair and mix-pair are returned, all pairs are non-overlapping.
This is used to calculate distribution of group-to-group and mix-to-mix differences.
Input:
    obsA: list of elements
    obsB: list of elements
        Note: obsA and obsB better have same length(TODO: check this in code)
    fraction: fraction of elements to switch
    size_of_each_group: size of each group selected
    N: number of groups
    seed: random seed
Output:
    groupA, groupB, mix1, mix2


**Source:** Line 48 in [shuffle.py](shuffle.py#L48)


---


### `exchange_group(groupA, groupB, fraction=0.5)`


Swith fraction of elements of groupA with groupB.
Note:
    1. groupA and groupB must have same length
    2. result groups have no intersection
Input:
    groupA: list of elements
    groupB: list of elements
    fraction: fraction of elements to switch
Output:
    mix_groupA, mix_groupB


**Source:** Line 27 in [shuffle.py](shuffle.py#L27)


---


### `get_shape(nested_list)`


Get shape of a nested list with any depth


**Source:** Line 2 in [shuffle.py](shuffle.py#L2)


---


### `random_group(group, size_of_each_group, N)`


Random sample each group N times.
Input:
    group: list of elements
    size_of_each_group: size of each group selected
    N: number of groups
Output:
    list of random groups


**Source:** Line 10 in [shuffle.py](shuffle.py#L10)


---


## spectral_ordering

**File:** `spectral_ordering.py`

**Public functions:** 1


### `spectral_ordering(filesp, threads=24)`


Input:
    filesp : DataFrame, must have pairs col
Output:
    DataFrame with additional col "order_index"


**Source:** Line 28 in [spectral_ordering.py](spectral_ordering.py#L28)


---


## structure.align

**File:** `structure/align.py`

**Public functions:** 1


### `align_structures(df1: pd.DataFrame, df2: pd.DataFrame, agg_chrom=True, ret3d: bool=False)`

**Returns:** `np.ndarray`


Compute the homogeneous transformation matrix or transformed coordinates to align structure 2 to structure 1.

Args:
    df1 (pd.DataFrame): DataFrame with MultiIndex (chr, pos) and columns ['x', 'y', 'z']
    df2 (pd.DataFrame): DataFrame with MultiIndex (chr, pos) and columns ['x', 'y', 'z']
    agg_chrom (bool, optional): Whether to aggregate points by chromosome (chr) first. Defaults to True.
    ret3d (bool, optional): If True, return transformed coordinates instead of the transformation matrix. Defaults to False.

Returns:
    np.ndarray: 
        - If ret3d=False: 4x4 homogeneous transformation matrix M
        - If ret3d=True: Transformed coordinates of structure 2 (n_points × 3)

Raises:
    ValueError: If input DataFrames have mismatched indices or columns

Example:
    # Return transformation matrix
    M = align_structures(df1, df2)
    
    # Return transformed coordinates (aggregated by chromosome)
    Q_aligned = align_structures(df1, df2, ret3d=True, agg_chrom=True)
    
    # Verify equivalence
    assert np.allclose(Q_aligned, (M @ np.hstack([df2.values, np.ones((len(df2),1))]).T).T[:,:3])


**Source:** Line 4 in [structure/align.py](structure/align.py#L4)


---


## structure.measure

**File:** `structure/measure.py`

**Public functions:** 15


### `C2T_diff(chunk)`


Compute the difference between centromere and telomere.
Use this to treat various centromere and telomere circumstances.


**Source:** Line 449 in [structure/measure.py](structure/measure.py#L449)


---


### `angle_between_vectors(vector_a, vector_b)`


Calculate the angle between two 3D vectors.

Input:
    vector_a (np.array): First 3D vector.
    vector_b (np.array): Second 3D vector.
Output:
    angle_degrees (float): The angle between the two vectors in degrees.


**Source:** Line 552 in [structure/measure.py](structure/measure.py#L552)


---


### `append_centelo_to_index(df, genome='mm10', dis=2000000.0, p=False, q=True)`


Append centromere and telomere information to the index.
Input:
    df: inner _3dg data structure (see parse_3dg in hires_utils)
    genome: genome name, use this to determine centromere and telomere positions
    dis: distance to centromere/telomere
    p: whether to include p arm
        For acrocentric chromosomes, set p to False
    q: whether to include q arm
Output:
    dat: data structure with centromere and telomere information


**Source:** Line 408 in [structure/measure.py](structure/measure.py#L408)


---


### `append_pm_to_index(df)`


Append paternal and maternal information to the index.
Input:
    df: inner _3dg data structure (see parse_3dg in hires_utils)
Output:
    dat: data structure with paternal and maternal information
TODO: A more general way to determine paternal/maternal rather than name matching


**Source:** Line 521 in [structure/measure.py](structure/measure.py#L521)


---


### `calc_depth(_3dg)`


Calculate the depth (from nuclear membrane) of a chromatin bin.
Input:
    inner _3dg data structure (parsing from hickit output .3dg file)
Output:
    same dataframe with new column 'depth'


**Source:** Line 199 in [structure/measure.py](structure/measure.py#L199)


---


### `cent2telo_vector(_3dg, genome='mm10', dis=2000000.0, p=False, q=True, show_chroms=False)`


Compute the C2T vector.
Input:
    _3dg: inner _3dg data structure (see parse_3dg in hires_utils)
    genome: genome name, use this to determine centromere and telomere positions
    dis: distance threshold to define regions near centromeres and telomeres
    p: whether to include p arm
        For acrocentric chromosomes, set p to False
    q: whether to include q arm
    show_chroms: whether to return chromosomal C2T vectors
Output:
    if show_chroms:
        res: dataframe of chromosomal C2T vectors
        c2t_vector: sum of all chromosomal C2T vectors
    else:
        normed_c2t_L: float of normalized C2T vector length
            normally in the range of 0.00x
Usage:
    _3dg = parse_3dg(
        _3dg_path
    )
    c2t_vec = cent2telo_vector(_3dg, genome="mm10", dis=2e6, p=False, q=True)


**Source:** Line 471 in [structure/measure.py](structure/measure.py#L471)


---


### `fetch_strong_violates(_3dg_f, pairs_f, genome, binsize, d_thres=4)`


Fetch strong constraint violations from pairs and _3dg data.
Input:
    _3dg_f: path to _3dg file
    pairs_f: path to pairs file
    genome: genome assembly identifier (e.g., 'mm10')
    binsize: Chromosome bin size in base pairs.
    d_thres: Distance threshold for defining strong violations.
Output:
    strong_violates: DataFrame containing strong constraint violations.


**Source:** Line 335 in [structure/measure.py](structure/measure.py#L335)


---


### `mt_strong_violation_ratio(sample_table, outfile, force=False, pairs_col='pairs_c12', _3dg_col='20k_g_struct1', genome=None, binsize=20000.0, d_thres=4, n_jobs=16)`


Process multiple samples in parallel to calculate strong constraint violations.
TODO: make explicit cache

Input:
    sample_table (pd.DataFrame): DataFrame containing sample metadata with columns 
                               '20k_g_struct1' and 'pairs_c12' for file paths.
    genome (str): Genome assembly identifier (e.g., 'mm10').
    binsize (int): Chromosome bin size in base pairs.
    d_thres (int): Distance threshold for defining strong violations.
    n_jobs (int): Number of parallel processes to use.
    outfile (str): Path to output JSON file storing results.

Output:
    None. Results are saved directly to the specified outfile.


**Source:** Line 360 in [structure/measure.py](structure/measure.py#L360)


---


### `parallel_light(plate, direction)`


Adding same direction to all origin points.
Input:
    plate: N * N * 3 array
    direction: (3,) 3D vector
Output:
    N * N * 6 array


**Source:** Line 20 in [structure/measure.py](structure/measure.py#L20)


---


### `perpendicular_unit_vector(A, B)`


Calculate a unit vector that is perpendicular to both A and B,
ensuring the resulting coordinate system (A, B, C) follows the right-hand rule.

Input:
    A (np.array): First 3D vector.
    B (np.array): Second 3D vector.

Output:
    C_normalized (np.array): Unit vector perpendicular to both A and B.


**Source:** Line 585 in [structure/measure.py](structure/measure.py#L585)


---


### `pm_vector(_3dg)`


Compute the paternal/maternal vector.
Assume chromosome names are in the format of "chr1(pat)/chr1(mat)".
Input:
    _3dg: inner _3dg data structure (see parse_3dg in hires_utils)
Output:
    pm_vector: paternal centroid point to maternal centroid point vector (mat - pat)


**Source:** Line 535 in [structure/measure.py](structure/measure.py#L535)


---


### `primary_views(_3dg, ngrid=16, method='ray', keep_main=True)`


Get primary views from three orthogonal faces of the nucleus' oriented bounding box.
Input:
    _3dg: inner _3dg data structure (parsing from hickit output .3dg file)
    ngrid: number of grids in each dimension
    method: "ray" (depth from each primary direction) or "distance" (distance from grid point to mesh surface)
Output:
    dict with keys (bases, name_of_vectors, primary_figures):
        bases: three major vectors (from center to three faces of the box) and center point vector; all are column arrays
        extent: extent of the box
        name_of_vectors: names of the three major vectors and "center"
        primary_figures: 3 * 2 figures (each figure is a 2D array), 2 direction for each major vector


**Source:** Line 42 in [structure/measure.py](structure/measure.py#L42)


---


### `radial_position(_3dg: pd.DataFrame)`


Calculate the radial position of each chromatin bin.
Input:
    inner _3dg data structure (parsing from hickit output .3dg file)
Output:
    same dataframe with new column 'radial_position'


**Source:** Line 219 in [structure/measure.py](structure/measure.py#L219)


---


### `restraint_violations(pairs, _3dg, genome, binsize, flavor='hickit')`


Calculate distances between pairs, this is used to evaluate the extent of how reconstructed structure
violates the input contact restraints.
Input:
    pairs: pairs data structure (parsing from hickit output .pairs file)
    _3dg: inner _3dg data structure (see parse_3dg in hires_utils)
    genome: genome name
    binsize: binsize used in the hic data
    flavor: flavor of the bins, see GenomeIdeograph
Output:
    pairs with new column 'distance'


**Source:** Line 300 in [structure/measure.py](structure/measure.py#L300)


---


### `say_cheese(plate, direction, scene)`


Take a photo. The plate way!


**Source:** Line 34 in [structure/measure.py](structure/measure.py#L34)


---


## structure.pileup

**File:** `structure/pileup.py`

**Public functions:** 16


### `consecutive_slice_pileup_bf(primary_views, targets, sample_table, features, step=2, upper=6, lower=-6, axis='ht')`


No documentation available.


**Source:** Line 739 in [structure/pileup.py](structure/pileup.py#L739)


---


### `consecutive_slice_pileup_bf_parallel(primary_views, targets, sample_table, features, step=2, upper=6, lower=-6, axis='ht', prefix=None, threads=8)`


No documentation available.


**Source:** Line 798 in [structure/pileup.py](structure/pileup.py#L798)


---


### `echo_file(file_path)`


No documentation available.


**Source:** Line 408 in [structure/pileup.py](structure/pileup.py#L408)


---


### `mix_layer(x_voxel, axis='ht', near=1)`


Mix the layer of the x_voxel with the near layer


**Source:** Line 897 in [structure/pileup.py](structure/pileup.py#L897)


---


### `pc_bin_feature(xyzf: pd.DataFrame, step: int, plot: bool=True)`

**Returns:** `pd.DataFrame`


Binnify feature in primary-coordinate(3D) space.
Input:
    xyzf: first 3 col must be ["ht","dv","lr"], index must be ["chrom","start"]
    step: bin size
Output:
    3d-binned feature
    if plot, return df with index ["ht","dv","lr"]
    else, return df with index ["chrom","start"]


**Source:** Line 192 in [structure/pileup.py](structure/pileup.py#L192)


---


### `pc_binagg_feature(xyzf, step, grouping=['ht', 'dv', 'lr'], agg=dict(density=('particle', 'sum')), min_particles={}, plot=True, global_min_particles=10)`


Binnify feature and aggregate whithin in primary-coordinate(3D) space.
Input:
    xyzf: first 3 col must be ["ht","dv","lr"]
    step: bin size
    min_particles: min particles in a bin to be valid, for each feature
    global_min_particles: default min particles, use this if feature not in min_particles
        set to None to disable all min_particles
Output:
    3d-binned feature


**Source:** Line 241 in [structure/pileup.py](structure/pileup.py#L241)


---


### `pileup_bf(primary_views, targets, sample_table, features, pb_features=None, agg={'density': ('particle', 'sum')}, sub=None, step=2, grouping=['ht', 'dv'], min_particles={}, global_min_particle=10)`


Pileup binned features.
Input:
    primary_views: to get three bases of each sample; structured data.
    targets: additional direction for each base. df.
    sample_table: to get path of structure(3dg) file. df.
    features: bulk features, df:=bin*feature
    pb_features: dict of df, df:=bin*sample
    agg: aggregation method; dict; find key in pb_features or columns of features.
    sub: epression for primary views, used for slicing.
    step: bin size
    grouping: bin grouping
    min_particles: min particles in a bin to be valid, for each feature
    global_min_particles: default min particles, use this if feature not in min_particles
        set to None to disable all min_particles
Output:
    bfs


**Source:** Line 269 in [structure/pileup.py](structure/pileup.py#L269)


---


### `pileup_bf_mt(sample_names: list, primary_views: dict, targets: pd.DataFrame, _3dg_files: Union[pd.DataFrame, List[str]], features: pd.DataFrame, outfile: Union[str, Path]=None, pb_features: Union[dict, pd.DataFrame]=None, pb_readers: dict=None, agg={'density': ('particle', 'sum')}, sub=None, step=2, grouping=['ht', 'dv', 'lr'], min_particles={}, global_min_particle=10, cache_dir: Union[str, Path]=None, n_threads=8)`


Pileup binned features.
Multi-threaded version.
Note: output not compatible with single-threade version.
Input:
    sample_names: list of sample names
    primary_views: to get three bases of each sample; structured data.
    targets: additional direction for each base.
        treat first 3 cols as ht, dv, lr directions
    _3dg_files: path of structure(3dg) file.
    features: bulk features, df:=bin*feature
    outfile: output file path
    pb_features:
        features that are unique for each sample
        dict of df, df:=bin*sample
        df: samples * features, each element is path to single feature file
    pb_readers:
        reader of feature files := dict(pb_feature name = reader function)
        reader function give a 2-level index series := (chr,pos):value
        will use a default reader for keys in pb_features that 
    agg: aggregation method; find key in pb_features or columns of features.
    sub: epression for primary views, used for slicing.
    step: bin size
    grouping: bin grouping
    min_particles: min particles in a bin to be valid, for each feature
    global_min_particles: default min particles, use this if feature not in min_particles
        set to None to disable all min_particles
    cache_dir: cache directory, required. Intermediate file names are not unique, isolated just by sample name.
    n_threads: number of threads
Output:
    bfs


**Source:** Line 410 in [structure/pileup.py](structure/pileup.py#L410)


---


### `pileup_bf_mt_task(sample_name: str, primary_view: dict, directions, _3dg_file, features, cache_file, pb_features=None, pb_readers=None, agg=None, sub=None, step=2, grouping=['ht', 'dv', 'lr'], min_particles={}, global_min_particle=10)`


No documentation available.


**Source:** Line 348 in [structure/pileup.py](structure/pileup.py#L348)


---


### `pileup_bf_parallel_task(primary_views, targets, sample_table, features, step, grouping, sub, axis, prefix, key)`


No documentation available.


**Source:** Line 768 in [structure/pileup.py](structure/pileup.py#L768)


---


### `project_back(sample_name: str, primary_view: dict, directions: pd.Series, _3dg_file: pd.DataFrame, spatial_features: pd.DataFrame, cache_file, step=2)`


Project back spatial feature to primary coordinates.
Input:
    sample_name: sample name
    primary_view: primary view
    directions: directions
    _3dg_file: 3dg file path
    spatial_features: spatial feature
        (ht, dv, lr) [feature1, feature2, ...]
    cache_file: cache file
    step: bin size


**Source:** Line 549 in [structure/pileup.py](structure/pileup.py#L549)


---


### `project_back_mt(sample_names: list, primary_views: dict, targets: pd.DataFrame, _3dg_files: Union[pd.DataFrame, List[str]], spatial_features: pd.DataFrame, outfile: Union[str, Path]=None, step=2, cache_dir: Union[str, Path]=None, n_threads=8)`


Project spatial features back to single-cell genome.
Input:
    sample_names: list of sample names
    primary_views: to get three bases of each sample; structured data.
    targets: additional direction for each base.
        treat first 3 cols as ht, dv, lr directions
    _3dg_files: path of structure(3dg) file.
    spatial_features: mean spatial features, df:=(ht,dv,lr)*features
    outfile: output file path
    step: bin size
    cache_dir: cache directory, required. Intermediate file names are not unique, isolated just by sample name.
    n_threads: number of threads
Output:
    bfs


**Source:** Line 630 in [structure/pileup.py](structure/pileup.py#L630)


---


### `project_back_mt_task(sample_name: str, primary_view: dict, directions: pd.Series, _3dg_file: pd.DataFrame, spatial_features: pd.DataFrame, cache_file, step=2)`


Project back spatial feature to primary coordinates.
A wrapper for project_back.
Input:
    sample_name: sample name
    primary_view: primary view
    directions: directions
    _3dg_file: 3dg file path
    spatial_features: spatial feature
        (ht, dv, lr) [feature1, feature2, ...]
    cache_file: cache file
    step: bin size


**Source:** Line 594 in [structure/pileup.py](structure/pileup.py#L594)


---


### `sig_primary_coords(primary_views, sig, _3dg)`


Transform xyz coordinates to primary coordinates.
Input:
    primary_views: primary_views[sample]
    sig: 3 iterable, ht,dv, lr
    _3dg: xyz of sample
Output:
    new _3dg; columns as head_tail coordinates, 
        dorsal_ventral coordinates,
        left_right coordinates; range 0-1


**Source:** Line 55 in [structure/pileup.py](structure/pileup.py#L55)


---


### `spatial_binnify(left: tuple=(-70, -50, -40), right: tuple=(80, 60, 50), step: int=2, plot: bool=True)`

**Returns:** `pd.DataFrame`


Generate spatial grid in primary-coordinate(3D) space.
Input:
    xyzf: first 3 col must be ["ht","dv","lr"], index must be ["chrom","start"]
    step: bin size
    plot: for plot compatibility, use left bound of each bin
Output:
    3d-binned feature
    if plot, return df with index ["ht","dv","lr"]
    else, return df with index ["chrom","start"]


**Source:** Line 15 in [structure/pileup.py](structure/pileup.py#L15)


---


### `voxelize(_3dg, primary_view, target, step=2, plot=True, observed=True, lr90=False)`


No documentation available.


**Source:** Line 835 in [structure/pileup.py](structure/pileup.py#L835)


---


## structure.utils

**File:** `structure/utils.py`

**Public functions:** 6


### `coord2point(coords, bases)`


Input:
    coord: column array, homogeneous coordinates (dim N+1 in dim N space)
    bases: column vectors (N+1 vectors for dim N space, the last col is origin point)
Output:
    points in space (dim N), column array


**Source:** Line 6 in [structure/utils.py](structure/utils.py#L6)


---


### `coords2point(coords, bases)`


Input:
    coords: column arrays, homogeneous coordinates (dim N+1 in dim N space)
    bases: column vectors (N+1 vectors for dim N space, the last col is origin point)
Output:
    points in space (dim N), column arrays


**Source:** Line 18 in [structure/utils.py](structure/utils.py#L18)


---


### `corners2edges(corners, atol=0.001, ret_xyz=False)`


Generate edges from random order of 8 corners.
Input:
    corners: np.ndarray or pd.DataFrame of shape (8, 3)
    atol: float, tolerance for checking if two vectors are close
    ret_xyz: whether return xyz of edges. If false, return index of corners
Output:
    edges: list of tuple of index of corners
        if ret_xyz is True, return list of tuple of xyz of edges


**Source:** Line 67 in [structure/utils.py](structure/utils.py#L67)


---


### `grid_cords(grid)`


No documentation available.


**Source:** Line 60 in [structure/utils.py](structure/utils.py#L60)


---


### `space_grid(bases, extent=(1, 1, 1), num=8)`


Generate (3D) spactial grid within a box defined by some bases.
Input:
    bases: 4 vectors. 3 from box center to its three faces. 1 the center.
    extent: size of the box.
    num: number of grids to sampling in each directon.
Output:
    dim(num, num, num, 3) ndarray, first three dims represent relative position of the grid point, 
    last dim is the real 3D coords of that grid point.
TODO:
    1. surpport any dim
    2. set num for each dim


**Source:** Line 32 in [structure/utils.py](structure/utils.py#L32)


---


### `vector_fbclose(v1, v2, atol=0.001)`


Check if two vectors are close to each other in same direction or opposite direction


**Source:** Line 62 in [structure/utils.py](structure/utils.py#L62)


---


## territory

**File:** `territory.py`

**Public functions:** 3


### `chrom_stat(acc, x)`


Add value to chromosome frequency stat.
Input:
    acc: original stat, dict, {chrom: count}, must be well initiated (has all possible chrom names as keys)
    x: input_line, (index, chromosome name), chromosome name must be in o.keys()
Ouput:
    new_acc: updated stat


**Source:** Line 54 in [territory.py](territory.py#L54)


---


### `intermingle(_3dg, min_dist=3, n_jobs=12, table=False)`


Calculate intermingling metrics for each particle in 3dg file.
Input:
    _3dg: dataframe, result of hic_basic.hicio parse_3dg
    min_dist: int, minimum distance to be considered as neighbor
    n_jobs: int, number of jobs to run in parallel
Output:
    intermingling_metrics: same length dataframe, index is particle index, columns are:
        chrom: str, chromosome name
        start: int, start position
        intermingling_ratio: float, ratio of intermingling
        multi_chrom_intermingling: float, shannon index of multi-chrom intermingling
        species_richness: int, number of species in the neighborhood


**Source:** Line 13 in [territory.py](territory.py#L13)


---


### `shannon_index(x)`


Calculate shannon index for a row of frequency table.


**Source:** Line 7 in [territory.py](territory.py#L7)


---


## utils

**File:** `utils.py`

**Public functions:** 8


### `argnatsort(array)`


No documentation available.


**Source:** Line 34 in [utils.py](utils.py#L34)


---


### `binnify(chromsizes, binsize)`


Divide a genome into evenly sized bins.

Parameters
----------
chromsizes : Series
    pandas Series indexed by chromosome name with chromosome lengths in bp.
binsize : int
    size of bins in bp

Returns
-------
bins : :py:class:`pandas.DataFrame`
    Dataframe with columns: ``chrom``, ``start``, ``end``.


**Source:** Line 93 in [utils.py](utils.py#L93)


---


### `gen_fileis(sample_table, dir_path, str_pat)`


Check and generate fileis(input file tables storing file path) for downstream calculating.
Input:
    sample_table: meta, sample_names as index, try to ensure all samples have their input file.
    dir_path: where to find those input files.
    str_pat: string pattern, gen exact file name. using {} to represent sample_name.
Output:
    pd.DataFrame


**Source:** Line 128 in [utils.py](utils.py#L128)


---


### `mouse_id2name(id_lists, ref_file, multi=False)`


Convert moust gene_IDs to gene_names.
Input:
    id_list: gene_ID list or list of gene_ID lists, do this because ref parsing is slow
        e.g. [["ENSG0000012345", "ENSG0000012346"], ["ENSG0000012347"]]
    ref_file: reference file
    multi: whether id_lists has multiple lists
Output:
    list of gene_names


**Source:** Line 172 in [utils.py](utils.py#L172)


---


### `natsort_key(s, _NS_REGEX=re.compile('(\\d+)', re.U))`


No documentation available.


**Source:** Line 32 in [utils.py](utils.py#L32)


---


### `read_chromsizes(filepath_or, name_patterns=('^chr[0-9]+$', '^chr[XY]$', '^chrM$'), all_names=False)`


Parse a ``<db>.chrom.sizes`` or ``<db>.chromInfo.txt`` file from the UCSC
database, where ``db`` is a genome assembly name.

note: sep = "       ", can't have header.

Parameters
----------
filepath_or : str or file-like
    Path or url to text file, or buffer.
name_patterns : sequence, optional
    Sequence of regular expressions to capture desired sequence names.
    Each corresponding set of records will be sorted in natural order.
all_names : bool, optional
    Whether to return all contigs listed in the file. Default is
    ``False``.

Returns
-------
:py:class:`pandas.Series`
    Series of integer bp lengths indexed by sequence name.

References
----------
* `UCSC assembly terminology <http://genome.ucsc.edu/FAQ/FAQdownloads.html#download9>`_
* `GRC assembly terminology <https://www.ncbi.nlm.nih.gov/grc/help/definitions>`_


**Source:** Line 40 in [utils.py](utils.py#L40)


---


### `two_sets(ref, check, warning=False)`


Check two sample list.
Input:
    ref: reference list, iterable.
    check: list to be checked, iterable.
Output:
    [missing, extra]; list of sets


**Source:** Line 8 in [utils.py](utils.py#L8)


---


### `z_score(array, mean=0, std=1)`


No documentation available.


**Source:** Line 194 in [utils.py](utils.py#L194)


---


## wet.adtools

**File:** `wet/adtools.py`

**Public functions:** 7


### `add_obsm(adata, data, obsm_key)`


Adding obsm to anndata.
Input:
    data: low-dim representation of samples. samples x feature
    obsm_key: name of the representation


**Source:** Line 327 in [wet/adtools.py](wet/adtools.py#L327)


---


### `clean_velocyto_names(adata, using_dup: list)`


No documentation available.


**Source:** Line 72 in [wet/adtools.py](wet/adtools.py#L72)


---


### `create_adata_layers(ws, funcs, debug=False)`


Create adata object form the layers
Input:
    ws: dict of layer name and dataframe
    funcs: functions to access layers
Output:
    adata object
TODO: figure out what cause the "Index of obs must match index of X." error


**Source:** Line 266 in [wet/adtools.py](wet/adtools.py#L266)


---


### `expand_df(target_df, ref_df, count_matrix=False, fill_na_value=0)`


expand target_df to shape(with same index, columns) of ref_df, fill with 0
Input:
    target_df: df to be transformed; sparse df
    ref_df: axis reference, must contain all rows and cols of target_df; sparse df
    count_matrix: if True, use sparse matrix, else use dense matrix


**Source:** Line 98 in [wet/adtools.py](wet/adtools.py#L98)


---


### `fix_velocyto_names(adata, annote, verbose=False)`


Fix obs.index of velocyto output.
e.g. 220110_embryo_8EZD5:2021123110 -> 2021123110.


**Source:** Line 63 in [wet/adtools.py](wet/adtools.py#L63)


---


### `gen_adata(qc, cache_dir, rewrite=[], debug=False)`


Generate anndata object from the QC data.
Input:
    qc: QC dataframe
    cache_dir: directory to store the generated adata object
    rewrite: list of dataframe name to rewrite

    --- layers ---
    expr: expression matrix
        Generate single expression matrix that contains all (has valid task directory) qc passed samples from qc file.
        qc file must contain "task_dirp" col and use sample ID as index.
        fix sample ID with "_" for umi_tools limitation. 
    g1: phase 0 expression matrix
    g2: phase 1 expression matrix
    velo: velocity matrix
    --- uns ---
    cdps: contact decay profile
    g1cs: phase 0 compartment strength
    g2cs: phase 1 compartment strength
    --- obs ---
    g1_UMIs: phase 0 UMIs
    g2_UMIs: phase 1 UMIs
    pm: phase 0 1 interaction(contacts)
    rs: repli_score
        Input:
            outdir: directory to store the repli_score intermediate files
    collect_hour: collection hour parsed from group name
    cell_type: cell type parsed from group name
TODO: for samples without results in some uns, add proper nan values


**Source:** Line 345 in [wet/adtools.py](wet/adtools.py#L345)


---


### `gen_cs(fileis)`


Generate compartment strength table from pre-computed files.


**Source:** Line 87 in [wet/adtools.py](wet/adtools.py#L87)


---


## wet.afbb

**File:** `wet/afbb.py`

**Public functions:** 18


### `add_extra(annote)`


Calculate contact_per_reads, umis_per_reads, rna_ratio.
Handling zero-division and NA values(and set to -1).
Output:
    new df with extra cols.


**Source:** Line 347 in [wet/afbb.py](wet/afbb.py#L347)


---


### `add_info(annote, rd)`


Add sample information from JSON files to an existing annotation DataFrame.

Args:
    annote: Existing annotation DataFrame
    rd: Directory containing JSON files with additional sample information
    
Returns:
    Combined DataFrame with original annotations and new sample information


**Source:** Line 219 in [wet/afbb.py](wet/afbb.py#L219)


---


### `add_mapping(annote, threads=16)`


Adding mapping and primary mapping rate of DNA library.
Input:
    old: sam version compatible to old pipeline. task_dirp/sam/xxx.aln.sam.gz


**Source:** Line 382 in [wet/afbb.py](wet/afbb.py#L382)


---


### `add_pairs(annote, task_dirp)`


try add pairs file path, by default using:
    pairs_0/*sample*.pairs.gz
    pairs_c1/*sample*.c1.pairs.gz
    pairs_c12/*sample*.c12.pairs.gz
    pairs_c123/*sample*.c123.pairs.gz
    dip/*sample*.dip.pairs.gz
assign None if file doesn't exist
Input:
    annote: must use sample name as index
Output:
    anntoe with extra cols


**Source:** Line 126 in [wet/afbb.py](wet/afbb.py#L126)


---


### `add_pairs_num(annote, threads, pairs_types=['pairs_c1', 'pairs_c12', 'pairs_c123', 'dip'])`


count contact number of pairs file in parallel
Input:
    annote: must have all pairs_type cols
Output:
    annote with extra cols, pairs_type + "_num" by default


**Source:** Line 56 in [wet/afbb.py](wet/afbb.py#L56)


---


### `add_umis(annote, umip)`


Adding umi stat to annote.
Accept a single umi count matrix path or an iterable (list/tuple) of umi count matrix paths.
Missing files will be skipped with a warning. If no valid files are found, assign 0 for both
"umis" and "genes" for all samples.

Input:
    umip: umi count matrix path or iterable of paths


**Source:** Line 233 in [wet/afbb.py](wet/afbb.py#L233)


---


### `check_RNA(task_dirp)`


No documentation available.


**Source:** Line 433 in [wet/afbb.py](wet/afbb.py#L433)


---


### `count_pairs(filename: str)`

**Returns:** `int`


No documentation available.


**Source:** Line 43 in [wet/afbb.py](wet/afbb.py#L43)


---


### `get_info(log_d)`


Extract and combine information from all JSON files in a directory.

Processes JSON files in chronological order, merging their contents.
Each JSON file should contain sample information as key-value pairs.

Args:
    log_d: Directory containing JSON files with sample information
    
Returns:
    DataFrame where columns are sample IDs and rows are attributes


**Source:** Line 196 in [wet/afbb.py](wet/afbb.py#L196)


---


### `listdir_bytime(folder)`


List all files in a directory sorted by modification time (oldest first).

Args:
    folder: Path to directory to scan
    
Returns:
    List of full file paths sorted by modification time


**Source:** Line 179 in [wet/afbb.py](wet/afbb.py#L179)


---


### `mapping_rate(bam, threads=4)`


No documentation available.


**Source:** Line 76 in [wet/afbb.py](wet/afbb.py#L76)


---


### `multi_re(key_array, re_array, base_array)`


Search according to regular expressions. Used in log file info extraction.
Print finding result and return finding as dict.
Input:
    key_array: feature name of each re target.
    re_s_array: re strings.
    base_array: divide by elements to get real ratio.
Output:
    searching result as dict.


**Source:** Line 91 in [wet/afbb.py](wet/afbb.py#L91)


---


### `multi_re_(logfile)`


No documentation available.


**Source:** Line 103 in [wet/afbb.py](wet/afbb.py#L103)


---


### `pick_useful(annote, features=['basic'])`


Pick useful cols from meta dataframe.
Input:
    annote: meta dataframe
    features: list of features to pick, default words are "DNA" and "DNA_RNA"
Output:
    new df with useful cols.


**Source:** Line 476 in [wet/afbb.py](wet/afbb.py#L476)


---


### `real_file(i)`


No documentation available.


**Source:** Line 84 in [wet/afbb.py](wet/afbb.py#L84)


---


### `task_stat(task_dirp, ref=None, threads=32)`

**Returns:** `pd.DataFrame`


Get all infomation about this bubble_flow task in the form of a dataframe.
The meta dataframe is constructed by adding additional columns to "contacts_info.csv" generated by bubble_flow.
TODO:
    make decent dtypes for every column step by step.


**Source:** Line 440 in [wet/afbb.py](wet/afbb.py#L440)


---


### `update_values(base: dict, record: str)`

**Returns:** `dict`


Update a base dictionary with values from a JSON record file.

For each key in the record, if the key exists in base, update the nested dictionary.
If the key doesn't exist, add the entire record to base.

Args:
    base: Base dictionary to be updated
    record: Path to JSON file containing key-value pairs to merge

Returns:
    Updated base dictionary


**Source:** Line 156 in [wet/afbb.py](wet/afbb.py#L156)


---


### `zcount(filename: str, target: str)`

**Returns:** `int`


No documentation available.


**Source:** Line 34 in [wet/afbb.py](wet/afbb.py#L34)


---


## wet.exp_record

**File:** `wet/exp_record.py`

**Public functions:** 10


### `add_cell_type(annote)`


Add cell_type col to input df according to group name.
Input must have "group" col.


**Source:** Line 215 in [wet/exp_record.py](wet/exp_record.py#L215)


---


### `add_group_hour(annote, col='collect_hour')`


Add additional order index col to input df indicating hour post fertilizing time when the sample was collected.
Input:
    annote: sample-feature table, must have a "group" col.
    col: name of new column.


**Source:** Line 203 in [wet/exp_record.py](wet/exp_record.py#L203)


---


### `add_group_order(annote)`


Add additional order index col indicating sample collection time order to input df.
Input must have a "group" col.


**Source:** Line 189 in [wet/exp_record.py](wet/exp_record.py#L189)


---


### `d2(time_str)`


No documentation available.


**Source:** Line 126 in [wet/exp_record.py](wet/exp_record.py#L126)


---


### `drange(a, b, seats)`


No documentation available.


**Source:** Line 12 in [wet/exp_record.py](wet/exp_record.py#L12)


---


### `find_day(group)`


Find day mark.


**Source:** Line 147 in [wet/exp_record.py](wet/exp_record.py#L147)


---


### `find_time_str(group)`


Find time substring in sample's group_string.


**Source:** Line 133 in [wet/exp_record.py](wet/exp_record.py#L133)


---


### `get_cell_type(grp)`


No documentation available.


**Source:** Line 220 in [wet/exp_record.py](wet/exp_record.py#L220)


---


### `parse_group_string(group_string, last_mouse)`


Parsing the cell meta data record from experiment.
Input:
    group_string: string to be parsed, should be in yaml format.
    last_mouse: largest mouse number in past experiments


**Source:** Line 21 in [wet/exp_record.py](wet/exp_record.py#L21)


---


### `plot_real_pseudo_time(adata, ps_col='velocity_pseudotime', hour=True)`


Plot experiment time with pseudotime.
(adata plot; scatter plot;)
Input:
    adata: AnnData object with "group" col in obs.
    ps_col: obs col storing pseudotime to plot real time with. 
    hour: sort by real experiment hour; False to use relative ordering.


**Source:** Line 232 in [wet/exp_record.py](wet/exp_record.py#L232)


---


## wet.gam

**File:** `wet/gam.py`

**Public functions:** 1


### `add_gam_mapping(sample_table_fp, ana_home)`

**Returns:** `pd.DataFrame`


Adding mapping and primary mapping rate to meta table.
Input:
    sample_table_fp: path to the meta table.
    ana_home: analysis home directory, check ana_home/logs/mapping.log
        for the mapping rate.
Output:
    sample_table: the original meta table with mapping rate added.


**Source:** Line 34 in [wet/gam.py](wet/gam.py#L34)


---


## wet.meta_trick

**File:** `wet/meta_trick.py`

**Public functions:** 3


### `last_mouse(df)`


Find last batch number.


**Source:** Line 36 in [wet/meta_trick.py](wet/meta_trick.py#L36)


---


### `merge_meta()`


Merging metadata object.
Input:
    DataFrame object seperate with , .


**Source:** Line 4 in [wet/meta_trick.py](wet/meta_trick.py#L4)


---


### `rmsd_filter(meta, reso='20k', rmsd_thres=1, samples=True, rmsd=False)`


No documentation available.


**Source:** Line 19 in [wet/meta_trick.py](wet/meta_trick.py#L19)


---


## wet.paracalc

**File:** `wet/paracalc.py`

**Public functions:** 6


### `count_MP(filename: str)`

**Returns:** `int`


No documentation available.


**Source:** Line 62 in [wet/paracalc.py](wet/paracalc.py#L62)


---


### `gen_PM_interactions(filesp, threads=32)`


MP: Maternal-Paternal interaction, e.g. chr1(mat)-chr1(pat), chr1(mat)-chr2(pat)
MMPP: Maternal-Maternal/Paternal-Paternal interaction, e.g. chr1(mat)-chr2(mat)
All_inter_contacts: MP + MMPP is all inter-chromosomal contacts.


**Source:** Line 80 in [wet/paracalc.py](wet/paracalc.py#L80)


---


### `gen_cdps(filesp, threads=32, range_dtype='string', pairs_col='pairs_c123', filter_str=None)`


No documentation available.


**Source:** Line 45 in [wet/paracalc.py](wet/paracalc.py#L45)


---


### `gen_cool(qc, cached, threads=16, sizef='/share/home/ychi/data/genome/GRCm38/mm10.chrom.sizes', binsize=40000, cool_addr='cool_paths.json')`


Generate cooler files from qc file.
Using mid-files of compartment strength pipeline is better.


**Source:** Line 15 in [wet/paracalc.py](wet/paracalc.py#L15)


---


### `gen_repli_score(qc, cached, threads=32, ref='/share/home/ychi/data/genome/GRCm38/mm10_repli_chip.wig')`


No documentation available.


**Source:** Line 110 in [wet/paracalc.py](wet/paracalc.py#L110)


---


### `safe_repli_score(pairsf, ref, outf)`


No documentation available.


**Source:** Line 102 in [wet/paracalc.py](wet/paracalc.py#L102)


---


## wet.qc

**File:** `wet/qc.py`

**Public functions:** 3


### `c_u_filter(meta, ccol, ucol, c_threshold, u_threshold)`


Filter out low-quality cells.
Input:
    ccol: column name of contact number
    ucol: column name of umis
    c_threshold: [low, high] for contacts
    u_threshold: [low, high] for contacts
Output:
    new df with qc-passed rows


**Source:** Line 3 in [wet/qc.py](wet/qc.py#L3)


---


### `plot_qc(meta, ccol='contacts', ucol='umis', c_threshold=[150000, 1000000], u_threshold=[10000, 600000])`


No documentation available.


**Source:** Line 21 in [wet/qc.py](wet/qc.py#L21)


---


### `plot_qc(meta, ccol='contacts', ucol='umis', c_threshold=[150000, 1000000], u_threshold=[10000, 600000])`


No documentation available.


**Source:** Line 52 in [wet/qc.py](wet/qc.py#L52)


---


## wet.rawcheck

**File:** `wet/rawcheck.py`

**Public functions:** 7


### `calc_md5(input_files, file_labels=None, result_file=None, force=False, threads=8)`


Calculate md5 checksums of input files.
Input:
    input_files: list of files to calculate md5 checksums
    file_labels: rename files in md5 checksum results, for example: /dir/file -> file
    result_file: file to write md5 checksum results
Output:
    True if all md5 checksums are OK, False otherwise


**Source:** Line 45 in [wet/rawcheck.py](wet/rawcheck.py#L45)


---


### `calculate_md5_for_file(file_path)`


Calculate md5 checksum for a file.
This give exactly the same result as md5sum command.
Input:
    file_path: path to the file
Output:
    file_path: path to the file
    md5: md5 checksum of the file


**Source:** Line 7 in [wet/rawcheck.py](wet/rawcheck.py#L7)


---


### `gen_sample_table(task_dirp, R1_file='_R1.fq.gz', R2_file='_R2.fq.gz')`


Create sample table for a raw directory.
Input:
    task_dirp: directory containing samples
    R1_file: R1 file suffix
    R2_file: R2 file suffix
Output:
    sample_table will be written to task_dirp/../sample_table.tsv


**Source:** Line 214 in [wet/rawcheck.py](wet/rawcheck.py#L214)


---


### `ls_md5_files(download_dir)`


Find all md5.txt files in download_dir.
Input:
    download_dir: directory containing downloaded files
Output:
    md5_lines: list of tuples, (parent_dir, md5, filename)


**Source:** Line 22 in [wet/rawcheck.py](wet/rawcheck.py#L22)


---


### `ls_samples(download_dir, show=False)`


Find all samples in download_dir. Samples are defined as the directories containing fastq files.
Input:
    download_dir: directory containing downloaded files
    show: print samples if True
Output:
    samples: list of sample names


**Source:** Line 134 in [wet/rawcheck.py](wet/rawcheck.py#L134)


---


### `split_download_dir(download_dir, task_dirps, sub='Rawdata')`


Split download_dir into multiple directories, each containing multiple samples and 1 md5 file for each sample.
Input:
    download_dir: directory containing downloaded files
    task_dirps: dictionary, {task_dirp: [sample_names]}
    sub: only consider samples that have this substring in their path
Output:
    None


**Source:** Line 161 in [wet/rawcheck.py](wet/rawcheck.py#L161)


---


### `verify_md5(download_dir, result_file, threads=8)`


Verify md5 checksums of downloaded files.
Input:
    download_dir: directory containing downloaded files
    result_file: file to write md5 checksum results
Output:
    True if all md5 checksums are OK, False otherwise


**Source:** Line 93 in [wet/rawcheck.py](wet/rawcheck.py#L93)


---


## wet.utils

**File:** `wet/utils.py`

**Public functions:** 1


### `check_input(filesp, cols)`


Check if the input files are valid.
Input:
    filesp: dataframe storing file paths.
    cols: cols to check; list or string.


**Source:** Line 4 in [wet/utils.py](wet/utils.py#L4)


---
