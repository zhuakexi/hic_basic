# hic_basic API Reference

*Auto-generated API documentation for public functions*

Total modules with public functions: 94

---

## Table of Contents

- [DI](#di) (18 functions)
- [Ps](#ps) (5 functions)
- [TAD.IS](#tad-is) (2 functions)
- [TAD.tad](#tad-tad) (7 functions)
- [TI.circle](#ti-circle) (1 functions)
- [TI.slingshot.plot](#ti-slingshot-plot) (3 functions)
- [TI.slingshot.utils](#ti-slingshot-utils) (5 functions)
- [api_query](#api_query) (13 functions)
- [atac](#atac) (6 functions)
- [basic](#basic) (5 functions)
- [binnify](#binnify) (1 functions)
- [calculate](#calculate) (2 functions)
- [cli.ana_describe](#cli-ana_describe) (2 functions)
- [cli.api](#cli-api) (2 functions)
- [cli.download](#cli-download) (2 functions)
- [cli.render](#cli-render) (3 functions)
- [cli.render3d](#cli-render3d) (2 functions)
- [cli.version](#cli-version) (2 functions)
- [cli.view3d](#cli-view3d) (2 functions)
- [clustering](#clustering) (1 functions)
- [compare_pixels](#compare_pixels) (5 functions)
- [compartment](#compartment) (8 functions)
- [coolstuff](#coolstuff) (32 functions)
- [crossover.detector_mii](#crossover-detector_mii) (5 functions)
- [crossover.io](#crossover-io) (7 functions)
- [crossover.pipeline](#crossover-pipeline) (4 functions)
- [crossover.pipeline_mii](#crossover-pipeline_mii) (8 functions)
- [cycle_phasing](#cycle_phasing) (2 functions)
- [data.file_formats](#data-file_formats) (2 functions)
- [embedding](#embedding) (5 functions)
- [genome](#genome) (4 functions)
- [heuristic_ordering](#heuristic_ordering) (1 functions)
- [hicbr](#hicbr) (2 functions)
- [hicio](#hicio) (29 functions)
- [impute.schicluster](#impute-schicluster) (3 functions)
- [impute.schicluster_gpu](#impute-schicluster_gpu) (10 functions)
- [impute.simpute](#impute-simpute) (6 functions)
- [inter_contact](#inter_contact) (3 functions)
- [interpolate](#interpolate) (9 functions)
- [lap_embedding](#lap_embedding) (3 functions)
- [mchr](#mchr) (1 functions)
- [metrics](#metrics) (4 functions)
- [nagano_cycle_phasing](#nagano_cycle_phasing) (6 functions)
- [paracalc](#paracalc) (1 functions)
- [phasing](#phasing) (11 functions)
- [pileup](#pileup) (6 functions)
- [pipeline.rule](#pipeline-rule) (8 functions)
- [plot.bar](#plot-bar) (1 functions)
- [plot.gene_track_parser](#plot-gene_track_parser) (2 functions)
- [plot.general](#plot-general) (6 functions)
- [plot.hic](#plot-hic) (22 functions)
- [plot.pca](#plot-pca) (2 functions)
- [plot.plot](#plot-plot) (15 functions)
- [plot.plotly_template](#plot-plotly_template) (1 functions)
- [plot.render](#plot-render) (21 functions)
- [plot.rna](#plot-rna) (3 functions)
- [plot.scanpy](#plot-scanpy) (5 functions)
- [plot.scatter](#plot-scatter) (5 functions)
- [plot.tech](#plot-tech) (1 functions)
- [plot.track](#plot-track) (9 functions)
- [plot.utils](#plot-utils) (20 functions)
- [pp](#pp) (1 functions)
- [pymol.color_centelo](#pymol-color_centelo) (2 functions)
- [pymol.glow](#pymol-glow) (1 functions)
- [pymol.sele_cent](#pymol-sele_cent) (2 functions)
- [pymol.sele_telo](#pymol-sele_telo) (2 functions)
- [pymol.utils](#pymol-utils) (2 functions)
- [sam](#sam) (7 functions)
- [scAB_embedding](#scab_embedding) (8 functions)
- [scripts.downsra](#scripts-downsra) (1 functions)
- [scripts.rescue_jupyter](#scripts-rescue_jupyter) (2 functions)
- [sequence](#sequence) (10 functions)
- [shuffle](#shuffle) (4 functions)
- [spectral_ordering](#spectral_ordering) (1 functions)
- [structure.align](#structure-align) (1 functions)
- [structure.measure](#structure-measure) (15 functions)
- [structure.pileup](#structure-pileup) (16 functions)
- [structure.utils](#structure-utils) (6 functions)
- [territory](#territory) (2 functions)
- [utils](#utils) (8 functions)
- [visualize.annotations](#visualize-annotations) (2 functions)
- [visualize.config](#visualize-config) (1 functions)
- [visualize.html](#visualize-html) (2 functions)
- [visualize.render](#visualize-render) (1 functions)
- [visualize.structure](#visualize-structure) (2 functions)
- [wet.adtools](#wet-adtools) (7 functions)
- [wet.afbb](#wet-afbb) (18 functions)
- [wet.exp_record](#wet-exp_record) (6 functions)
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


### `ps_curve(coolp, view_df=None, all_region=False, nproc=4, clr_weight_name='weight', drop2diags=False, **kwargs)`


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

**Public functions:** 7


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


### `zero(R, t, length, delta)`


No documentation available.


**Source:** Line 101 in [TAD/tad.py](TAD/tad.py#L101)


---


## TI.circle

**File:** `TI/circle.py`

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


**Source:** Line 33 in [TI/circle.py](TI/circle.py#L33)


---


## TI.slingshot.plot

**File:** `TI/slingshot/plot.py`

**Public functions:** 3


### `plot_clusters(slingshot_obj, title: str='Slingshot Clusters', width: int=600, height: int=600, point_size: int=6, point_opacity: float=0.8, save_path: Optional[str]=None, show: bool=True)`


绘制细胞聚类（按细胞类型着色）

Parameters
----------
slingshot_obj : Slingshot
    Slingshot 对象
title : str, optional
    图形标题
width : int, optional
    图形宽度
height : int, optional
    图形高度
point_size : int, optional
    散点大小
point_opacity : float, optional
    散点透明度
save_path : str, optional
    保存路径
show : bool, optional
    是否显示图形

Returns
-------
fig : plotly.graph_objects.Figure
    Plotly 图形对象


**Source:** Line 210 in [TI/slingshot/plot.py](TI/slingshot/plot.py#L210)


---


### `plot_curves(slingshot_obj, title: str='Slingshot Principal Curves', width: int=700, height: int=600, curve_width: int=4, curve_opacity: float=0.9, show_cells: bool=True, save_path: Optional[str]=None, show: bool=True)`


绘制主曲线（可选显示细胞散点）

Parameters
----------
slingshot_obj : Slingshot
    Slingshot 对象
title : str, optional
    图形标题
width : int, optional
    图形宽度
height : int, optional
    图形高度
curve_width : int, optional
    曲线线宽
curve_opacity : float, optional
    曲线透明度
show_cells : bool, optional
    是否显示细胞散点
save_path : str, optional
    保存路径
show : bool, optional
    是否显示图形

Returns
-------
fig : plotly.graph_objects.Figure
    Plotly 图形对象


**Source:** Line 306 in [TI/slingshot/plot.py](TI/slingshot/plot.py#L306)


---


### `plot_slingshot(slingshot_obj, title: str='Slingshot Trajectory', width: int=700, height: int=600, point_size: int=5, point_opacity: float=0.6, curve_width: int=4, curve_opacity: float=0.8, colorscale: str='Viridis', show_legend: bool=True, save_path: Optional[str]=None, show: bool=True)`


绘制 Slingshot 轨迹推断结果

主绘图函数：绘制细胞散点图（按 pseudotime 着色）并叠加所有主曲线

Parameters
----------
slingshot_obj : Slingshot
    Slingshot 对象（需已调用 fit()）
title : str, optional
    图形标题
width : int, optional
    图形宽度（像素）
height : int, optional
    图形高度（像素）
point_size : int, optional
    散点大小
point_opacity : float, optional
    散点透明度
curve_width : int, optional
    曲线线宽
curve_opacity : float, optional
    曲线透明度
colorscale : str, optional
    颜色映射名称（plotly 支持的颜色映射）
show_legend : bool, optional
    是否显示图例
save_path : str, optional
    保存路径（如 "output.png"）
show : bool, optional
    是否显示图形

Returns
-------
fig : plotly.graph_objects.Figure
    Plotly 图形对象

Examples
--------
>>> from hic_basic.TI.slingshot import Slingshot
>>> from hic_basic.TI.slingshot.plot import plot_slingshot
>>>
>>> # 初始化并拟合
>>> slingshot = Slingshot(data=adata, celltype_key='celltype', start_node='SSC')
>>> slingshot.fit(num_epochs=10)
>>>
>>> # 绘图
>>> fig = plot_slingshot(slingshot, title="Cell Trajectory")
>>> fig.show()


**Source:** Line 42 in [TI/slingshot/plot.py](TI/slingshot/plot.py#L42)


---


## TI.slingshot.utils

**File:** `TI/slingshot/utils.py`

**Public functions:** 5


### `infer_cluster_label_indices(cluster_labels)`


推断聚类标签索引

参数:
    cluster_labels: 字符串标签列表或整数聚类分配列表

返回:
    聚类分配的 numpy 数组


**Source:** Line 61 in [TI/slingshot/utils.py](TI/slingshot/utils.py#L61)


---


### `isint(x)`


检查是否为整数


**Source:** Line 26 in [TI/slingshot/utils.py](TI/slingshot/utils.py#L26)


---


### `isstr(x)`


检查是否为字符串


**Source:** Line 41 in [TI/slingshot/utils.py](TI/slingshot/utils.py#L41)


---


### `mahalanobis(X1, X2, S1, S2)`


计算两个向量之间的马氏距离

参数:
    X1: 第一个向量 (n,)
    X2: 第二个向量 (n,)
    S1: 第一个协方差矩阵 (n, n)
    S2: 第二个协方差矩阵 (n, n)

返回:
    马氏距离 (标量)


**Source:** Line 8 in [TI/slingshot/utils.py](TI/slingshot/utils.py#L8)


---


### `scale_to_range(x, a=0, b=1)`


将数组缩放到指定范围

参数:
    x: 输入数组
    a: 目标最小值
    b: 目标最大值

返回:
    缩放后的数组


**Source:** Line 46 in [TI/slingshot/utils.py](TI/slingshot/utils.py#L46)


---


## api_query

**File:** `api_query.py`

**Public functions:** 13


### `entry_to_json(entry: Dict[str, Any], include_resolved_source: bool=False)`

**Returns:** `Dict[str, Any]`


Return a stable JSON-serializable entry payload.


**Source:** Line 192 in [api_query.py](api_query.py#L192)


---


### `first_doc_line(docstring: str)`

**Returns:** `str`


Extract the first non-empty line from a docstring.


**Source:** Line 48 in [api_query.py](api_query.py#L48)


---


### `iter_entries(api_index: Index)`

**Returns:** `Iterable[Dict[str, Any]]`


Yield normalized function entries from an API index.


**Source:** Line 24 in [api_query.py](api_query.py#L24)


---


### `list_module_functions(api_index: Index, module_query: str)`

**Returns:** `List[Dict[str, Any]]`


List functions for exact or partial module matches.


**Source:** Line 65 in [api_query.py](api_query.py#L65)


---


### `list_modules(api_index: Index)`

**Returns:** `List[Dict[str, Any]]`


Return module names and public function counts.


**Source:** Line 57 in [api_query.py](api_query.py#L57)


---


### `load_api_index(json_path: Optional[Path]=None)`

**Returns:** `Index`


Load an API index from a path or from package resources.


**Source:** Line 16 in [api_query.py](api_query.py#L16)


---


### `normalize_entry(module: str, func: Dict[str, Any])`

**Returns:** `Dict[str, Any]`


Return a function entry with stable derived fields.


**Source:** Line 31 in [api_query.py](api_query.py#L31)


---


### `resolve_source(entry: Dict[str, Any])`

**Returns:** `Dict[str, Any]`


Resolve source location in the installed environment when possible.


**Source:** Line 167 in [api_query.py](api_query.py#L167)


---


### `resolve_symbol(api_index: Index, symbol: str)`

**Returns:** `Tuple[str, List[Dict[str, Any]]]`


Resolve a short function name or fully qualified hic_basic symbol.


**Source:** Line 145 in [api_query.py](api_query.py#L145)


---


### `score_entry(entry: Dict[str, Any], tokens: List[str])`

**Returns:** `int`


Score an entry for query tokens.


**Source:** Line 119 in [api_query.py](api_query.py#L119)


---


### `search_by_keyword(api_index: Index, keyword: str, case_sensitive: bool=False)`

**Returns:** `List[Dict[str, Any]]`


Search for functions by keyword in docstrings.


**Source:** Line 89 in [api_query.py](api_query.py#L89)


---


### `search_by_name(api_index: Index, query: str, case_sensitive: bool=False)`

**Returns:** `List[Dict[str, Any]]`


Search for functions by name.


**Source:** Line 78 in [api_query.py](api_query.py#L78)


---


### `search_entries(api_index: Index, query: str, limit: int=10)`

**Returns:** `List[Dict[str, Any]]`


Rank API entries using deterministic keyword matching.


**Source:** Line 101 in [api_query.py](api_query.py#L101)


---


## atac

**File:** `atac.py`

**Public functions:** 6


### `calculate_coverage(bam_path: str, tss_regions, read_length: int, sampling=None)`


Compute per-base coverage across TSS regions with optional sampling.

Input
   bam_path (str): Path to a coordinate-sorted, indexed BAM file.
   tss_regions (BedTool): Expanded TSS regions (output of `get_tss_regions`).
   read_length (int): Sequencing read length. Retained for API compatibility.
   sampling (float or int, optional): Sampling strategy.
       - If float: Fraction of regions to sample (0 < value <= 1).
       - If int: Exact number of regions to sample.
       - If None: Process all regions.

Output
   Tuple[np.ndarray, List[int]]:
       coverage_array (np.ndarray): Integer array with shape (n_TSS, region_width).
           Uncovered positions are represented by zeros.
       distances (List[int]): Per-base offsets from the center of each region.

Notes
   Regions clipped at chromosome boundaries are padded with zeros to the
   widest selected region. The BAM file is opened once for the complete
   operation, avoiding one `pysam.depth` process per region.


**Source:** Line 170 in [atac.py](atac.py#L170)


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


### `normalize_coverage(coverage_array: np.ndarray, read_length: int, greenleaf_norm: bool=True, bam_path: Optional[str]=None)`

**Returns:** `np.ndarray`


Normalizes coverage using Greenleaf-style or per million mapped reads normalization.

Input:
   coverage_array (np.ndarray): 2D coverage array (shape: [n_TSS, bins]).
   read_length (int): Sequencing read length. Retained for API compatibility.
   greenleaf_norm (bool): Whether to apply Greenleaf-style normalization (default: True).
   bam_path (str, optional): Indexed BAM used to calculate mapped-read scaling
       when `greenleaf_norm` is False.

Output:
   np.ndarray: Normalized coverage array.


**Source:** Line 277 in [atac.py](atac.py#L277)


---


### `plot_fragment_length_distribution(fragment_lengths: np.ndarray, **kwargs)`

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


**Source:** Line 318 in [atac.py](atac.py#L318)


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


**Source:** Line 29 in [binnify.py](binnify.py#L29)


---


## calculate

**File:** `calculate.py`

**Public functions:** 2


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


### `task_wrapper(module_name, func_name, *args, **kwargs)`


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


## cli.ana_describe

**File:** `cli/ana_describe.py`

**Public functions:** 2


### `add_arguments(subparser)`


No documentation available.


**Source:** Line 12 in [cli/ana_describe.py](cli/ana_describe.py#L12)


---


### `run(args)`


No documentation available.


**Source:** Line 32 in [cli/ana_describe.py](cli/ana_describe.py#L32)


---


## cli.api

**File:** `cli/api.py`

**Public functions:** 2


### `add_arguments(subparser)`


No documentation available.


**Source:** Line 23 in [cli/api.py](cli/api.py#L23)


---


### `run(args)`


No documentation available.


**Source:** Line 47 in [cli/api.py](cli/api.py#L47)


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


## cli.render3d

**File:** `cli/render3d.py`

**Public functions:** 2


### `add_arguments(subparser)`


No documentation available.


**Source:** Line 29 in [cli/render3d.py](cli/render3d.py#L29)


---


### `run(args)`


No documentation available.


**Source:** Line 91 in [cli/render3d.py](cli/render3d.py#L91)


---


## cli.version

**File:** `cli/version.py`

**Public functions:** 2


### `add_arguments(subparser)`


No documentation available.


**Source:** Line 3 in [cli/version.py](cli/version.py#L3)


---


### `run(*args, **kwargs)`


No documentation available.


**Source:** Line 5 in [cli/version.py](cli/version.py#L5)


---


## cli.view3d

**File:** `cli/view3d.py`

**Public functions:** 2


### `add_arguments(subparser)`


No documentation available.


**Source:** Line 16 in [cli/view3d.py](cli/view3d.py#L16)


---


### `run(args)`


No documentation available.


**Source:** Line 75 in [cli/view3d.py](cli/view3d.py#L75)


---


## clustering

**File:** `clustering.py`

**Public functions:** 1


### `do_kmeans(adata: ad.AnnData, layer: str, **k)`

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


### `pick_band_pixels(coolp, threshold=5000000, **kwargs)`

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


### `plot_compare_pixels_density(this_coolp, ref_coolp, qrange=None, sample_titles=['THIS', 'REF'], bins=100, ignore_diags=0, clr_weight_name='weight', lognorm=False, vlognorm=False, retdata=False, vrange=None, label_x=0.01, label_y=0.015, debug=False, **kwargs)`


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


### `gam2cool(gam_filepat, bin_file, outfile, force=False, intra_only=True, **kwargs)`


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


### `gen_config(sample_table, cfg, assembly='mm10', bwa='/share/home/ychi/data/genome/GRCm38/bwa_index/mm10.*', chrom_size='/share/home/ychi/data/genome/GRCm38/mm10.len.forCool.tsv', template='/shareb/ychi/ana/distiller-nf/project.yml', **args)`


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


### `mt_pairs2cool(sample_table, outdir, pairs_col='pairs_c123', suffix='.cool', nproc=16, **kwargs)`


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


### `pixels2cool(pixels: pd.DataFrame, bins: pd.DataFrame, outfile, columns: list=['count'], dtypes: dict={'count': int}, **kwargs)`


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


### `schicluster2cool(filesp, fo, genome, binsize, ordered=True, ensure_sorted=True, **args)`


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


## crossover.detector_mii

**File:** `crossover/detector_mii.py`

**Public functions:** 5


### `estimate_emissions_from_counts(read_counts: pd.DataFrame, rough_labels: pd.Series)`

**Returns:** `Dict[str, Tuple[float, float, float]]`


Estimate emission probabilities from read counts and rough labels.

Args:
    read_counts: DataFrame with columns 'ref', 'alt', 'other'
    rough_labels: Series with rough state labels ('PP', 'PM', 'MM', 'U')

Returns:
    Dictionary mapping state name to (p_ref, p_alt, p_other)


**Source:** Line 463 in [crossover/detector_mii.py](crossover/detector_mii.py#L463)


---


### `estimate_transitions_rough(rough_labels: pd.Series, pseudocount: float=1.0)`

**Returns:** `Tuple[float, float, float, float]`


Estimate transition probabilities from rough labels.

Self transitions are included in the row totals so the returned parameters
represent per-adjacent-observation probabilities, not conditional
probabilities after a non-self jump has already occurred. The returned rho
is the per-target U->biological-state probability used by the HMM; the
total U exit probability is 3*rho.

Args:
    rough_labels: Series with rough state labels
    pseudocount: Pseudocount for smoothing (deprecated, kept for API compatibility)

Returns:
    Tuple of (alpha, beta, gamma, rho)


**Source:** Line 512 in [crossover/detector_mii.py](crossover/detector_mii.py#L512)


---


### `get_default_transition_params()`

**Returns:** `Dict[str, float]`


Get default transition matrix parameters.

These parameters represent a conservative model with:
- Low transition rates between biological states (PP, PM, MM)
- Rare entry into U state
- Moderate exit from U state

Returns:
    Dictionary with keys: alpha, beta, gamma, rho


**Source:** Line 33 in [crossover/detector_mii.py](crossover/detector_mii.py#L33)


---


### `rough_classification(read_counts: pd.DataFrame, depth_threshold: int=2, other_threshold: float=0.4, balance_threshold: float=0.25)`

**Returns:** `pd.Series`


Perform rough classification of SNPs based on read counts.

Args:
    read_counts: DataFrame with columns 'ref', 'alt', 'other'
    depth_threshold: Minimum total depth for confident call (default: 2)
    other_threshold: Maximum other fraction for confident call (default: 0.4)
    balance_threshold: Threshold for classifying as balanced (PM) (default: 0.25)
                      PP: ref_frac > 0.75, MM: alt_frac > 0.75, PM: |ref-alt| < 0.25.
                      Ratio-intermediate observations are labeled AMBIG, not U.

Returns:
    Series with rough state labels. U is reserved for low-depth or high-other
    observations; AMBIG marks ratio-intermediate observations.


**Source:** Line 607 in [crossover/detector_mii.py](crossover/detector_mii.py#L607)


---


### `simulate_mii_crossover_counts(chromosome: str='chr1', chromosome_length: int=100000000, snp_spacing: int=1000, crossover_rate_per_mb: float=0.02, mean_depth: float=5.0, sequencing_error_rate: float=0.005, depth_cv: float=0.55, local_depth_window_snps: int=100, random_state: Optional[Union[int, np.random.Generator]]=None)`

**Returns:** `MiiCrossoverSimulation`


Simulate SNP read counts for mII sister-chromatid crossover inference.

The simulated mII signal is the aggregate of two independent chromatid
haplotype paths. Each chromatid starts as either P or M, switches haplotype
at crossover events, and contributes equally to the observed read pool.
Read counts are generated with WGS-like local depth variation, Poisson
sampling, and random sequencing errors.

Args:
    chromosome: Chromosome name assigned to simulated SNPs.
    chromosome_length: Chromosome length in bp.
    snp_spacing: Distance between adjacent SNPs in bp.
    crossover_rate_per_mb: Per-chromatid crossover probability per Mb.
    mean_depth: Mean total read depth per SNP.
    sequencing_error_rate: Per-read probability of observing the wrong base.
    depth_cv: Coefficient of variation for local depth multipliers.
    local_depth_window_snps: Number of SNPs sharing a local depth bias.
    random_state: Optional seed or numpy Generator for reproducibility.

Returns:
    MiiCrossoverSimulation with read-count DataFrame, crossover positions,
    and parameter metadata.


**Source:** Line 662 in [crossover/detector_mii.py](crossover/detector_mii.py#L662)


---


## crossover.io

**File:** `crossover/io.py`

**Public functions:** 7


### `extract_genotype_arrays(df: pd.DataFrame)`

**Returns:** `Tuple[np.ndarray, List[str], List[str], List[str]]`


Extract arrays needed for HMM from merged DataFrame.

Args:
    df: Merged DataFrame with sperm_genotype, haplotype_1, haplotype_2 columns

Returns:
    Tuple of (positions, sperm_genotypes, haplotype_1, haplotype_2)


**Source:** Line 381 in [crossover/io.py](crossover/io.py#L381)


---


### `load_phased_snp_from_hic_basic(phased_snp_file: str, chromosome: Optional[str]=None)`

**Returns:** `pd.DataFrame`


Load phased SNP data from hic_basic.phasing module format.

This function reads the output format from hic_basic.phasing.parse_phased_snp
and converts it to a DataFrame suitable for crossover detection.

Args:
    phased_snp_file: Path to phased SNP file (TSV or parquet)
    chromosome: Optional chromosome to filter

Returns:
    DataFrame with columns: chromosome, position, snp_id,
                           haplotype_1, haplotype_2


**Source:** Line 238 in [crossover/io.py](crossover/io.py#L238)


---


### `load_phasing_data(phasing_file: str, chromosome: Optional[str]=None)`

**Returns:** `pd.DataFrame`


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


**Source:** Line 176 in [crossover/io.py](crossover/io.py#L176)


---


### `load_snp_data(snp_file: str, sample_id: Optional[str]=None, chromosome: Optional[str]=None)`

**Returns:** `pd.DataFrame`


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


**Source:** Line 24 in [crossover/io.py](crossover/io.py#L24)


---


### `merge_snp_with_phasing(snp_df: pd.DataFrame, phasing_df: pd.DataFrame, on_cols: Optional[List[str]]=None)`

**Returns:** `pd.DataFrame`


Merge SNP data with phasing information.

Args:
    snp_df: SNP genotype DataFrame
    phasing_df: Phasing information DataFrame
    on_cols: Columns to merge on (default: ['chromosome', 'position', 'snp_id'])

Returns:
    Merged DataFrame with both genotype and phasing information


**Source:** Line 346 in [crossover/io.py](crossover/io.py#L346)


---


### `save_results(results: CrossoverResults, output_dir: str, prefix: Optional[str]=None)`

**Returns:** `Dict[str, str]`


Save crossover detection results to files.

Args:
    results: CrossoverResults object
    output_dir: Output directory path
    prefix: Optional prefix for output files

Returns:
    Dictionary mapping output type to file path


**Source:** Line 274 in [crossover/io.py](crossover/io.py#L274)


---


### `validate_input_data(snp_df: pd.DataFrame, phasing_df: Optional[pd.DataFrame]=None)`

**Returns:** `Tuple[bool, List[str]]`


Validate input data for crossover detection.

Args:
    snp_df: SNP genotype DataFrame
    phasing_df: Optional phasing DataFrame

Returns:
    Tuple of (is_valid, list of issues)


**Source:** Line 401 in [crossover/io.py](crossover/io.py#L401)


---


## crossover.pipeline

**File:** `crossover/pipeline.py`

**Public functions:** 4


### `detect_crossovers_end2end(sam_file: str, ref_snp_file: str, phasing_file: str, output_dir: str, sample_id: Optional[str]=None, n_jobs: int=1, cache_dir: Optional[str]=None, window_size: int=10, min_fragment_size: int=300000, error_rate: float=0.05, chromosomes: Optional[List[str]]=None)`

**Returns:** `Union[CrossoverResults, Dict[str, CrossoverResults]]`


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


**Source:** Line 429 in [crossover/pipeline.py](crossover/pipeline.py#L429)


---


### `get_pipeline_status(sam_file: str, ref_snp_file: str, cache_dir: str, sample_id: Optional[str]=None, window_size: int=10, min_fragment_size: int=300000, error_rate: float=0.05)`

**Returns:** `Dict[str, bool]`


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


**Source:** Line 693 in [crossover/pipeline.py](crossover/pipeline.py#L693)


---


### `process_sample_parallel(sample_ids: List[str], sam_files: Dict[str, str], ref_snp_file: str, phasing_file: str, output_dir: str, n_jobs: int=4, cache_dir: Optional[str]=None, window_size: int=10, min_fragment_size: int=300000, error_rate: float=0.05)`

**Returns:** `Dict[str, Dict[str, Any]]`


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


**Source:** Line 353 in [crossover/pipeline.py](crossover/pipeline.py#L353)


---


### `resume_pipeline(sam_file: str, ref_snp_file: str, phasing_file: str, output_dir: str, sample_id: Optional[str]=None, cache_dir: Optional[str]=None, **kwargs)`

**Returns:** `CrossoverResults`


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


**Source:** Line 738 in [crossover/pipeline.py](crossover/pipeline.py#L738)


---


## crossover.pipeline_mii

**File:** `crossover/pipeline_mii.py`

**Public functions:** 8


### `bin_by_snp(read_counts_df: pd.DataFrame, n_snps: int=10)`

**Returns:** `pd.DataFrame`


Merge adjacent SNPs into units based on SNP count.

Args:
    read_counts_df: DataFrame with columns [chromosome or chrom, pos or position, ref, alt, other]
    n_snps: Number of SNPs per unit (default: 10)

Returns:
    DataFrame with aggregated units, columns: [chromosome, start_pos, end_pos, ref, alt, other, n_snps, position]


**Source:** Line 55 in [crossover/pipeline_mii.py](crossover/pipeline_mii.py#L55)


---


### `convert_unit_coordinates_to_snp(unit_results: SisterChromatidResults, original_read_counts: pd.DataFrame, n_snps: int=10)`

**Returns:** `SisterChromatidResults`


Convert unit-level state assignments back to original SNP coordinates.

Args:
    unit_results: Results from HMM at unit level
    original_read_counts: Original SNP-level read counts
    n_snps: Number of SNPs per unit

Returns:
    SisterChromatidResults with SNP-level states


**Source:** Line 1865 in [crossover/pipeline_mii.py](crossover/pipeline_mii.py#L1865)


---


### `detect_sister_chromatid_states_end2end(bam_file: str, snp_file: str, output_dir: str, sample_id: Optional[str]=None, chromosomes: Optional[List[str]]=None, cache_dir: Optional[str]=None, transition_params: Union[str, Dict[str, float]]='default', emission_params: Optional[Dict[str, Tuple[float, float, float]]]=None, alpha: float=0.01, beta: float=0.02, gamma: float=0.001, rho: float=0.1, depth_threshold: int=2, other_threshold: float=0.4, balance_threshold: float=0.25, binning_method: Optional[str]=None, binning_snp_number: int=250, transition_param_unit: str='step', estimate_use: Union[str, List[str]]='all', transition_min_segment_bins: int=3, transition_u_gamma_prior: Optional[float]=None, transition_u_gamma_factor: float=0.1, transition_u_rho_prior: float=0.32, u_emission_other_prior: float=0.05, estimate_only: bool=False, params_vis_prefix: Optional[Union[str, Path]]=None, verbose: int=0)`

**Returns:** `Union[SisterChromatidResults, Dict[str, SisterChromatidResults], Dict[str, Any]]`


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
    transition_params: Transition parameter source ("default" | "estimated") or alpha/beta/gamma/rho dict
    emission_params: Optional specified emission probabilities by PP/PM/MM/U label
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
    transition_param_unit: Unit for specified transition parameters:
        "step", "mb", or "snp". Estimated parameters are always per-step.
    estimate_use: Chromosomes used for parameter estimation when chromosomes
        is provided. Use "all" for all extracted chromosomes, or a list of
        chromosome names.
    transition_min_segment_bins: Raw rough-label segments shorter than this
        many bins are marked as gaps and skipped for transition estimation.
    transition_u_gamma_prior: Optional fixed gamma prior used for U entry
        when binning_method="by_snp" and transition_params="estimated". If
        None, uses transition_u_gamma_factor * min(alpha, beta).
    transition_u_gamma_factor: Factor for deriving the by-SNP U-entry prior
        from estimated alpha/beta when transition_u_gamma_prior is None.
    transition_u_rho_prior: By-SNP prior for U -> each biological state.
    u_emission_other_prior: By-SNP fixed U-emission other probability
        used when emission_params is not specified.
    estimate_only: If True, stop after parameter estimation and return parameter diagnostics
    params_vis_prefix: Prefix for parameter diagnostic PNG files when verbose >= 1.
        Writes rough_label_classification, rough_decision_distributions,
        rough_segment_lengths, rough_label_filtering, emission_matrix, and
        transition_matrix PNGs. During full decoding it also writes
        final_assigned_label_count PNGs.
    verbose: If >=1, print parameter estimation diagnostics

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


**Source:** Line 2334 in [crossover/pipeline_mii.py](crossover/pipeline_mii.py#L2334)


---


### `estimate_mii_hmm_parameters(read_counts_df: pd.DataFrame, emission_params: Optional[Dict[str, Tuple[float, float, float]]]=None, transition_params: Union[str, Dict[str, float]]='default', alpha: float=0.01, beta: float=0.02, gamma: float=0.001, rho: float=0.1, depth_threshold: int=2, other_threshold: float=0.4, balance_threshold: float=0.25, binning_method: Optional[str]=None, binning_snp_number: int=250, transition_param_unit: str='step', transition_min_segment_bins: int=3, transition_u_gamma_prior: Optional[float]=None, transition_u_gamma_factor: float=0.1, transition_u_rho_prior: float=0.32, u_emission_other_prior: float=0.05, params_vis_prefix: Optional[Union[str, Path]]=None, verbose: int=0)`

**Returns:** `Dict[str, Any]`


Estimate or resolve mII HMM parameters from read-count data.

Args:
    read_counts_df: DataFrame with chromosome, position, ref, alt, other.
    emission_params: Optional specified emission probabilities by label.
    transition_params: "default", "estimated", or a dict with alpha/beta/gamma/rho.
    alpha, beta, gamma, rho: Fixed transition parameters used for "default".
    depth_threshold: Minimum depth for rough classification.
    other_threshold: Maximum other fraction for rough classification.
    balance_threshold: Ref/alt balance threshold for rough classification.
    binning_method: None or "by_snp".
    binning_snp_number: Number of SNPs per bin when binning_method="by_snp".
    transition_param_unit: Unit for caller-specified transition parameters:
        "step" for per HMM observation, "mb" for per-Mb rates, or "snp" for
        per-SNP probabilities. Estimated parameters are always per-step.
    transition_min_segment_bins: Raw rough-label segments shorter than this
        many bins are marked as gaps and skipped for transition estimation.
    transition_u_gamma_prior: Optional fixed gamma prior used for U entry
        when binning_method="by_snp" and transition_params="estimated". If
        None, uses transition_u_gamma_factor * min(alpha, beta).
    transition_u_gamma_factor: Factor for deriving the by-SNP U-entry prior
        from estimated alpha/beta when transition_u_gamma_prior is None.
    transition_u_rho_prior: By-SNP prior for U -> each biological state.
    u_emission_other_prior: By-SNP fixed U-emission other probability
        used when emission_params is None. Ref and alt share the remaining
        probability equally.
    params_vis_prefix: Prefix for diagnostic PNG files. When verbose >= 1,
        writes rough_label_classification, rough_decision_distributions,
        rough_segment_lengths, rough_label_filtering, emission_matrix, and
        transition_matrix PNGs. During full decoding it also writes
        final_assigned_label_count PNGs.
    verbose: If >=1, print diagnostics and, when params_vis_prefix is set, write PNGs.

Returns:
    Dictionary containing binned counts, rough labels, emission parameters,
    transition parameters, transition matrix, and diagnostic tables.


**Source:** Line 1592 in [crossover/pipeline_mii.py](crossover/pipeline_mii.py#L1592)


---


### `extract_read_counts(bam_file: str, snp_file: str, output_file: Optional[str]=None, min_mapq: int=20, min_baseq: int=20)`

**Returns:** `pd.DataFrame`


Extract read counts at SNP positions from BAM file.

Args:
    bam_file: Path to BAM/BCF file
    snp_file: Path to SNP file (VCF or TSV with positions)
    output_file: Optional output file path
    min_mapq: Minimum mapping quality
    min_baseq: Minimum base quality

Returns:
    DataFrame with columns: chromosome, position, ref, alt, other


**Source:** Line 1937 in [crossover/pipeline_mii.py](crossover/pipeline_mii.py#L1937)


---


### `plot_mii_assignment_track(pipeline_output: Any, region: Optional[Union[str, List[Tuple[str, int]], Region, Tuple[Any, ...]]]=None, output_path: Optional[Union[str, Path]]=None, title: str='mII_final_assignment', width: int=1400, height: int=360)`


Plot the mII final assignment label track.

Parameters
----------
pipeline_output
    Pipeline output as a DataFrame, tabular file, output directory, pipeline
    result dictionary, or SisterChromatidResults. The resolved table must
    contain chromosome, state, and position or start_pos/end_pos columns.
region
    Region to plot using hic_basic standard representations, such as
    "chr1:1000000-2000000", [(chrom, start), (chrom, end)], or Region.
    If None, all rows are shown.
output_path
    Optional PNG path. If provided, the Plotly figure is written to this path.
title
    Figure title.
width, height
    PNG export size in pixels when output_path is provided.

Returns
-------
plotly.graph_objects.Figure
    The generated label track figure.


**Source:** Line 1000 in [crossover/pipeline_mii.py](crossover/pipeline_mii.py#L1000)


---


### `plot_mii_final_assignment(pipeline_output: Any, region: Optional[Union[str, List[Tuple[str, int]], Region, Tuple[Any, ...]]]=None, output_path: Optional[Union[str, Path]]=None, title: str='Final assigned labels: genome_wide_label + genome_wide_count', width: int=1400, height: int=760)`


Plot final mII HMM labels together with ref/alt/other counts.

Parameters
----------
pipeline_output
    Pipeline output as a DataFrame, tabular file, output directory, pipeline
    result dictionary. The resolved table must contain chromosome, position
    (or start_pos/end_pos), state, ref, alt, and other columns.
region
    Region to plot using hic_basic standard representations, such as
    "chr1:1000000-2000000", [(chrom, start), (chrom, end)], or Region.
    If None, all rows are shown.
output_path
    Optional PNG path. If provided, the Plotly figure is written to this path.
title
    Figure title.
width, height
    PNG export size in pixels when output_path is provided.

Returns
-------
plotly.graph_objects.Figure
    The generated Plotly figure.


**Source:** Line 1093 in [crossover/pipeline_mii.py](crossover/pipeline_mii.py#L1093)


---


### `save_results_mii(results: Union[SisterChromatidResults, Dict[str, SisterChromatidResults]], output_dir: str, sample_id: str)`

**Returns:** `Dict[str, str]`


Save sister chromatid inference results to files.

Args:
    results: Results object or dictionary
    output_dir: Output directory
    sample_id: Sample identifier

Returns:
    Dictionary mapping output type to file path


**Source:** Line 2641 in [crossover/pipeline_mii.py](crossover/pipeline_mii.py#L2641)


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


## data.file_formats

**File:** `data/file_formats.py`

**Public functions:** 2


### `get_format(key: str)`

**Returns:** `dict | None`


Look up a file format definition by key.

Parameters
----------
key : str
    Format key (e.g., ``"cool"``, ``"bam"``).

Returns
-------
dict or None
    Format definition dict, or None if not found.


**Source:** Line 142 in [data/file_formats.py](data/file_formats.py#L142)


---


### `list_formats()`

**Returns:** `list[str]`


Return sorted list of all registered format keys.


**Source:** Line 158 in [data/file_formats.py](data/file_formats.py#L158)


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
    from hic_basic.TI.circle import angle_pseudotime

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


## genome

**File:** `genome.py`

**Public functions:** 4


### `chrom_is_ploidified(chrom)`


Detect ploidified chromosome names, like "chr1(mat)", "chr1(pat)", "chr1a".
Mirrors the suffixes handled by chrom_rm_suffix, anchored at the name end
so scaffold/contig names with internal underscores are not flagged.
Input:
    chrom: pd.Index or pd.Series of chromosome names
Output:
    bool: True if any name carries a ploidy suffix


**Source:** Line 648 in [genome.py](genome.py#L648)


---


### `index_merge_hom(df)`


Cancel ploidiness in index. You can then merge homologous chromosomes use functions like groupby.
This will change index inplace.
Input:
    df: DataFrame with MultiIndex ["chrom", "start"], chom is ploidified, like "chr1(mat)"
Output:
    a dataframe with modified level 0 index value, like "chr1(mat)" -> "chr1"


**Source:** Line 659 in [genome.py](genome.py#L659)


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


### `sort_chrom(df, genome)`


Sort dataframe by natural chromosome order.
TODO:
    1. ask whether to rm chromint after sorting
    2. ask whether to reset index
Input:
    df: must have a column named "chrom"
Output:
    df: sorted by natural chromosome order


**Source:** Line 678 in [genome.py](genome.py#L678)


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

**Public functions:** 2


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


## hicio

**File:** `hicio.py`

**Public functions:** 29


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


### `parse_gtf(file, ID=True, name=True, mgi_id=True, **args)`


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


## impute.schicluster_gpu

**File:** `impute/schicluster_gpu.py`

**Public functions:** 10


### `gpu_processing_worker(input_queue, output_queue, processing_params, verbosity=0)`


GPU worker thread: impute matrices sequentially and hand them to the
saver queue.

Runs as a single thread in the parent process because CUDA cannot be
used in forked subprocesses; GPU work is inherently sequential (one
device, one stream of batched solves), so one worker suffices.

``processing_params`` keys: ``method`` ('torch_ver' or 'numpy_ver'),
``device`` ('gpu' or 'cpu'), ``if_convolve``, ``if_rwr``,
``ignore_diag``, ``norm``, ``oldway``, and the three dtypes
(``convolve_dtype`` / ``kernel_dtype`` / ``rwr_dtype``).


**Source:** Line 393 in [impute/schicluster_gpu.py](impute/schicluster_gpu.py#L393)


---


### `hic_imputation_multiprocess(coolps, outfiles, chroms=None, resolution=10000, band_width=4000000, ignore_diag=True, norm=True, verbosity=0, nproc_loading=1, nproc_saving=5, method='torch_ver', device='gpu', oldway=False, if_convolve=True, if_rwr=True, convolve_dtype=None, kernel_dtype=None, rwr_dtype=None)`


Multi-process Hi-C imputation using a producer-consumer pipeline.

Layout: ``nproc_loading`` loader processes -> 1 GPU worker thread ->
``nproc_saving`` saver processes.  Queues are ``JoinableQueue`` with
capacity ``3 x n_workers`` to buffer while a slower stage catches up;
``None`` poison pills shut each stage down in order after all work is
drained.

Args:
    coolps (list): paths to cool files (one per task).
    outfiles (list): output ``.npz`` paths, one per task.
    chroms (list or str): chromosome per task, or a single chromosome
        applied to all tasks.
    resolution (int): resolution of the Hi-C data (bp).
    band_width (int): maximum cis interaction distance (bp).
    ignore_diag (bool): whether to zero the main diagonal.
    norm (bool): whether to full-sum-normalize.
    verbosity (int): 0 = silent, >=1 = pipeline progress logs.
    nproc_loading (int): number of loader processes.
    nproc_saving (int): number of saver processes.
    method (str): ``'torch_ver'`` (GPU path) or ``'numpy_ver'`` (CPU,
        in-thread).
    device (str): ``'gpu'`` or ``'cpu'`` (only for ``'torch_ver'``).
    oldway (bool): legacy diagonal +1 mode for the CPU path.
    if_convolve (bool): whether to apply Gaussian convolution.
    if_rwr (bool): whether to run the RWR step.
    convolve_dtype: convolution output dtype; torch dtype for
        ``'torch_ver'`` (float32 default), numpy dtype for
        ``'numpy_ver'`` (float64 default).
    kernel_dtype: Gaussian kernel dtype (float).
    rwr_dtype: RWR solve dtype.


**Source:** Line 538 in [impute/schicluster_gpu.py](impute/schicluster_gpu.py#L538)


---


### `load_single_hic_matrix(coolp, chrom='chr1', resolution=10000, band_width=4000000, verbosity=0)`


Load a single Hi-C matrix from a cool file into the pipeline format.

Args:
    coolp (str): path to the cool file.
    chrom (str): chromosome (or region) to extract.
    resolution (int): resolution of the Hi-C data (bp).
    band_width (int): maximum cis interaction distance (bp); ``None``
        disables the upper-band constraint.
    verbosity (int): 0 = silent, >=1 = log loading time.

Returns:
    tuple: ``(mat_X, params)`` where ``mat_X`` has shape
    ``(1, 1, H, W)`` and ``params`` carries ``chrom``, ``resolution``,
    ``band_width``, ``upper_band`` and ``coolp``.


**Source:** Line 319 in [impute/schicluster_gpu.py](impute/schicluster_gpu.py#L319)


---


### `matrix_loading_worker(input_queue, output_queue, verbosity=0)`


Loader process: fetch ``(coolp, chrom, resolution, band_width, outfile)``
jobs and push ``(mat_X, params, outfile)`` tuples to the processing queue.

A ``None`` item on the input queue is a poison pill that stops the
worker.


**Source:** Line 362 in [impute/schicluster_gpu.py](impute/schicluster_gpu.py#L362)


---


### `save_single_result(imputed_mat_X, params, outfile, verbosity=0)`


Save one imputed matrix as a compressed CSR ``.npz``.

The stored arrays are ``data`` / ``indices`` / ``indptr`` / ``shape``;
load with ``csr_matrix((z["data"], z["indices"], z["indptr"]),
shape=z["shape"])``.  The stored matrix is the imputed, normalized,
band-limited upper triangle (plus main diagonal).

Args:
    imputed_mat_X (np.ndarray): imputed matrix, shape ``(H, W)``.
    params (dict): pipeline params (informational; not stored).
    outfile (str): output ``.npz`` path.
    verbosity (int): 0 = silent, >=1 = per-stage timing logs.


**Source:** Line 505 in [impute/schicluster_gpu.py](impute/schicluster_gpu.py#L505)


---


### `saving_worker(input_queue, verbosity=0)`


Saver process: convert imputed matrices to CSR and write compressed
``.npz`` files.

A ``None`` item on the input queue is a poison pill that stops the
worker.


**Source:** Line 477 in [impute/schicluster_gpu.py](impute/schicluster_gpu.py#L477)


---


### `schicluster_imputation_cpu(mat, alpha=0.05, kernel_size=3, sigma=2, if_convolve=True, ignore_diag=True, oldway=False, norm=True, upper_band=None, convolve_dtype=np.float64, kernel_dtype=np.float64, rwr_dtype=np.float64, if_rwr=True)`


Schicluster imputation of one contact matrix, NumPy reference.

Same parameter surface as the GPU version, so the two paths can be
cross-validated.  Note: with ``ignore_diag=True`` the input ``mat`` is
modified in place (main diagonal zeroed) — pass a copy if the input must
be preserved.

Args:
    mat (np.ndarray): input contact matrix, shape ``(H, W)``.
    alpha (float): RWR restart parameter.
    kernel_size (int): size of the Gaussian kernel.
    sigma (float): standard deviation of the Gaussian kernel.
    if_convolve (bool): whether to apply Gaussian convolution.
    ignore_diag (bool): whether to zero the main diagonal (before and
        after imputation).
    oldway (bool): legacy diagonal +1 mode (differs from the production
        mode when convolution has placed non-zero values on the
        sub/super-diagonal; do not compare against the reference data).
    norm (bool): whether to full-sum-normalize the result.
    upper_band (int): if not None, keep only the main diagonal and the
        upper ``upper_band`` super-diagonals, zero the rest.
    convolve_dtype (np.dtype): dtype of the convolution output (can be
        int for research purposes, e.g. ``np.int32`` in production).
    kernel_dtype (np.dtype): dtype of the Gaussian kernel (should be
        float for proper convolution).
    rwr_dtype (np.dtype): dtype for the RWR solve (float64 default).
    if_rwr (bool): whether to run the RWR step.

Returns:
    np.ndarray: imputed matrix, shape ``(H, W)``.


**Source:** Line 115 in [impute/schicluster_gpu.py](impute/schicluster_gpu.py#L115)


---


### `schicluster_imputation_gpu(mat, alpha=0.05, kernel_size=3, sigma=2, if_convolve=True, ignore_diag=True, norm=True, upper_band=None, verbosity=0, device=None, convolve_dtype=torch.float32, kernel_dtype=torch.float32, rwr_dtype=torch.float32, if_rwr=True)`


Schicluster imputation of a batch of contact matrices, GPU version.

Matches the CPU reference (:func:`schicluster_imputation_cpu`) step for
step, but operates on a batch with ``torch.nn.functional.conv2d`` and a
batched ``torch.linalg.solve``.

Args:
    mat: batch of input matrices, shape ``(batch_size, 1, H, W)``
        (ndarray or tensor).
    alpha (float): RWR restart parameter.
    kernel_size (int): size of the Gaussian kernel.
    sigma (float): standard deviation of the Gaussian kernel.
    if_convolve (bool): whether to apply Gaussian convolution.
    ignore_diag (bool): whether to zero the main diagonal (before and
        after imputation).
    norm (bool): whether to full-sum-normalize the result.
    upper_band (int): if not None, keep only the main diagonal and the
        upper ``upper_band`` super-diagonals, zero the rest.
    verbosity (int): 0 = silent, >=1 = per-stage timing logs.
    device (torch.device): device to run on; defaults to CUDA if
        available, else CPU.
    convolve_dtype (torch.dtype): dtype of the convolution output
        (``torch.int32`` in production).
    kernel_dtype (torch.dtype): dtype of the Gaussian kernel (float).
    rwr_dtype (torch.dtype): dtype for the RWR solve (float32 in
        production).
    if_rwr (bool): whether to run the RWR step.

Returns:
    torch.Tensor: imputed matrices on CPU, shape ``(B, 1, H, W)``.


**Source:** Line 199 in [impute/schicluster_gpu.py](impute/schicluster_gpu.py#L199)


---


### `solve_rwr_inverse(stoch_matrix, alpha=0.05)`


Solve the Random Walk with Restart (RWR) equations by matrix inverse.

Solves ``(I - (1-alpha) * P^T) * R = alpha * I`` for the per-source
restart distributions and symmetrizes the result.

Args:
    stoch_matrix (np.ndarray): column-stochastic matrix ``P``, shape
        ``(N, N)``.
    alpha (float): restart probability.

Returns:
    np.ndarray: symmetrized RWR solution, shape ``(N, N)``.


**Source:** Line 48 in [impute/schicluster_gpu.py](impute/schicluster_gpu.py#L48)


---


### `solve_rwr_inverse_gpu(stoch_matrix, alpha=0.05, rwr_dtype=torch.float32, device=None)`


GPU version of :func:`solve_rwr_inverse`, batched over ``(B, C, N, N)``.

A single batched ``torch.linalg.solve`` replaces the per-matrix
``scipy.linalg.solve`` of the CPU reference; this is the main speedup of
the GPU path.

Args:
    stoch_matrix: column-stochastic matrix or tensor, shape
        ``(batch_size, channels, N, N)``.
    alpha (float): restart probability.
    rwr_dtype (torch.dtype): arithmetic dtype for the solve.
    device (torch.device): device to run on; defaults to CUDA if
        available, else CPU.

Returns:
    torch.Tensor: symmetrized RWR solution, shape ``(B, C, N, N)``.


**Source:** Line 73 in [impute/schicluster_gpu.py](impute/schicluster_gpu.py#L73)


---


## impute.simpute

**File:** `impute/simpute.py`

**Public functions:** 6


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

**Public functions:** 9


### `compute_center_of_mass(trajectory, grid_values)`


Compute center of mass (weighted average time position) for each feature.

The center of mass is computed as:
    com = sum(grid_values * |trajectory|) / sum(|trajectory|)

Using absolute values to handle both positive and negative trajectories.

Parameters
----------
trajectory : np.ndarray
    Trajectory matrix of shape (n_features, n_grid_points).
grid_values : np.ndarray
    Grid point values of shape (n_grid_points,).

Returns
-------
com : np.ndarray
    Center of mass for each feature, shape (n_features,).


**Source:** Line 238 in [interpolate.py](interpolate.py#L238)


---


### `compute_continuous_peak_time(trajectory, grid_values)`


Compute continuous peak time for each feature using parabolic interpolation.

For each feature, fits a parabola through the peak grid point and its two
neighbors, then computes the vertex position as the continuous peak time.

Parameters
----------
trajectory : np.ndarray
    Trajectory matrix of shape (n_features, n_grid_points).
    Each row is the GP-smoothed trajectory for one feature.
grid_values : np.ndarray
    Grid point values of shape (n_grid_points,). The actual time/pseudotime
    values at each grid point.

Returns
-------
peak_time : np.ndarray
    Continuous peak time for each feature, shape (n_features,).
    If a feature has no clear peak (monotonic or flat), returns the grid
    position of the maximum value.
peak_value : np.ndarray
    The interpolated peak value at the continuous peak time.


**Source:** Line 144 in [interpolate.py](interpolate.py#L144)


---


### `correct_pseudotime(adata, pseudotime_col, n_pcs=15, inplace=True)`


Correct pseudotime of adata to make it linearly correlated with a distance metrix.
Inspired by Alpert et al. (2018)
Input:
    adata: AnnData object, must have adata.obs[pseudotime_col]
    pseudotime_col: pseudotime column name
Output:
    new adata object


**Source:** Line 1495 in [interpolate.py](interpolate.py#L1495)


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


**Source:** Line 1230 in [interpolate.py](interpolate.py#L1230)


---


### `gp_interpolate(adata, pseudotime_col, layer='X', gene_mask_key=None, numPts=14, grid_method='linspace', grid_custom=None, endpoint_mode='half-interval', alpha=1e-06, length_scale=0.1, length_scale_bounds=(0.01, 1.0), n_restarts_optimizer=5, normalize_y=True, result_key='gp_interpolate', inplace=True)`


Interpolate single-trajectory expression dynamics with Gaussian Processes.

Parameters
----------
adata : anndata.AnnData
    AnnData object containing expression values and pseudotime annotations.
pseudotime_col : str
    Column in ``adata.obs`` containing numeric pseudotime values.
layer : str, optional
    Expression source. Use ``"X"`` for ``adata.X`` or the name of an
    entry in ``adata.layers``.
gene_mask_key : str, optional
    Boolean column in ``adata.var`` used to select genes for fitting.
numPts : int, optional
    Number of interpolation points on normalized pseudotime. Default is 14.
    Ignored when ``grid_method='optimal'`` (auto-determined) or
    ``grid_method='custom'`` (uses ``grid_custom``).
grid_method : str, optional
    Strategy for constructing the interpolation time grid. Options:
    - "linspace": evenly spaced points (default)
    - "quantile": quantile-based spacing (equal cells per bin)
    - "density": density-adaptive spacing (more points in dense regions)
    - "optimal": use optimal binning to determine numPts automatically
    - "custom": use user-provided grid (requires ``grid_custom``)
grid_custom : array-like, optional
    Custom grid values in range [0, 1], strictly increasing.
    Only used when ``grid_method='custom'``.
endpoint_mode : str, optional
    How to handle endpoints. Default is "half-interval". Options:
    - "half-interval": endpoints at [half_interval, 1-half_interval]
      to avoid sparse cells at t=0 and t=1 (recommended)
    - "full": endpoints at [0, 1] (legacy behavior)
alpha : float, optional
    Fixed observation noise variance passed to the GP regressor.
length_scale : float, optional
    Initial RBF length scale.
length_scale_bounds : tuple, optional
    Bounds for RBF length-scale optimization.
n_restarts_optimizer : int, optional
    Number of optimizer restarts for hyperparameter fitting.
normalize_y : bool, optional
    Whether to normalize the target values before fitting.
result_key : str, optional
    Key used to store outputs in ``adata.uns``.
inplace : bool, optional
    If ``True``, write results to ``adata`` directly. Otherwise fit on a
    copy and return the copy.

Returns
-------
anndata.AnnData
    AnnData with GP interpolation results stored in ``adata.uns[result_key]``.


**Source:** Line 1347 in [interpolate.py](interpolate.py#L1347)


---


### `plot_gp_heatmap(adata, pseudotime_col, layer='X', gene_mask_key=None, numPts=14, grid_method='linspace', endpoint_mode='half-interval', scaling='z_score', figsize=(8, 6), dpi=100, save_path=None, show_grid=True, show_colorbar=True, feature_order=None)`


Plot a single heatmap for GP interpolation results.

This is a convenience wrapper around plot_trajectory_heatmap() for AnnData-based
RNA/gene expression data. For more control or for non-RNA data (e.g., loops),
use plot_trajectory_heatmap() directly.

Parameters
----------
adata : anndata.AnnData
    AnnData object with expression data and pseudotime.
pseudotime_col : str
    Column name in adata.obs containing pseudotime values.
layer : str, optional
    Expression layer to use. Default is "X".
gene_mask_key : str, optional
    Boolean column in adata.var for gene selection.
numPts : int, optional
    Number of grid points. Default is 14.
grid_method : str, optional
    Grid construction method. Default is "linspace".
    Options: "linspace", "quantile", "density", "optimal", "custom".
endpoint_mode : str, optional
    Endpoint handling mode. Default is "half-interval".
    Options: "half-interval", "full".
scaling : str, optional
    Color scaling. Default is "z_score".
    Options: "z_score", "zero_one".
figsize : tuple, optional
    Figure size (width, height) in inches. Default is (8, 6).
dpi : int, optional
    Figure DPI. Default is 100.
save_path : str, optional
    Path to save the figure. If None, figure is not saved.
show_grid : bool, optional
    Show vertical lines at grid point positions. Default is True.
show_colorbar : bool, optional
    Show colorbar. Default is True.
feature_order : np.ndarray, optional
    Indices to reorder features (genes). If None, genes are shown in original order.
    Typically obtained from sort_features_by_trajectory().

Returns
-------
fig : matplotlib.figure.Figure
    The matplotlib figure object.
axes : np.ndarray
    Array of matplotlib axes objects.

See Also
--------
plot_trajectory_heatmap : General-purpose trajectory heatmap for any feature type.


**Source:** Line 890 in [interpolate.py](interpolate.py#L890)


---


### `plot_trajectory_heatmap(adata, pseudotime_col, layer='X', gene_mask_key=None, numPts=14, grid_methods=None, endpoint_mode='half-interval', scaling='z_score', figsize=(10, 12), dpi=100)`


Internal function: Visualize and compare multiple grid methods side-by-side.

This is a debugging utility, not intended for production use.
Use plot_gp_heatmap() for single-method visualization.

Parameters
----------
adata : anndata.AnnData
    AnnData object with expression data and pseudotime.
pseudotime_col : str
    Column name in adata.obs containing pseudotime values.
layer : str, optional
    Expression layer to use. Default is "X".
gene_mask_key : str, optional
    Boolean column in adata.var for gene selection.
numPts : int, optional
    Number of grid points. Default is 14.
grid_methods : list, optional
    List of grid methods to compare. Default is ["linspace", "quantile", "density"].
endpoint_mode : str, optional
    Endpoint handling mode. Default is "half-interval".
scaling : str, optional
    Color scaling: "z_score" or "zero_one". Default is "z_score".
figsize : tuple, optional
    Figure size (width, height) in inches. Default is (10, 12).
dpi : int, optional
    Figure DPI. Default is 100.

Returns
-------
fig : matplotlib.figure.Figure
    The matplotlib figure object.
axes : np.ndarray
    Array of matplotlib axes objects.


**Source:** Line 341 in [interpolate.py](interpolate.py#L341)


---


### `plot_trajectory_heatmap(trajectory_raw, pseudotime, trajectory_gp=None, grid_values=None, grid_cell_positions=None, feature_order=None, feature_labels=None, scaling='z_score', cmap='RdBu_r', vmin=None, vmax=None, figsize=(8, 6), dpi=100, save_path=None, show_grid=True, grid_line_color='black', grid_line_style='--', grid_line_alpha=0.7, show_colorbar=True, left_title='Raw (cell space)', right_title='GP-smoothed (grid space)', suptitle=None)`


Plot trajectory heatmap for features (genes, loops, etc.) with GP-smoothed trajectories.

This is a context-agnostic visualization function that works with any feature type
(RNA genes, chromatin loops, etc.) as long as the data is provided as numpy arrays.

Parameters
----------
trajectory_raw : np.ndarray or pd.DataFrame
    Raw trajectory matrix of shape (n_features, n_cells). Each row is a feature's
    trajectory across cells. Can also be a DataFrame with feature names as index.
pseudotime : np.ndarray or pd.Series
    Pseudotime values for cells, shape (n_cells,). Used to sort cells and compute
    grid line positions.
trajectory_gp : np.ndarray, optional
    GP-smoothed trajectory matrix of shape (n_features, n_grid_points).
    If None, only the raw trajectory panel is shown.
grid_values : np.ndarray, optional
    Grid point values in pseudotime space, shape (n_grid_points,).
    Used for right panel x-axis labels and to compute grid_cell_positions if needed.
grid_cell_positions : np.ndarray, optional
    Grid point positions in cell rank space [0, 1], shape (n_grid_points,).
    If None but grid_values is provided, computed from pseudotime.
feature_order : np.ndarray, optional
    Indices to reorder features (rows). If None, features are shown in original order.
    Typically obtained from sort_features_by_trajectory().
feature_labels : list or np.ndarray, optional
    Labels for features (rows). If None, uses integer indices.
scaling : str, optional
    Color scaling method. Options:
    - "z_score": Z-score normalization per feature
    - "zero_one": Min-max normalization to [0, 1]
cmap : str or matplotlib.colors.Colormap, optional
    Colormap for the heatmap. Default is "RdBu_r" (diverging).
vmin, vmax : float, optional
    Color scale limits. If None, computed from data based on scaling.
figsize : tuple, optional
    Figure size (width, height) in inches. Default is (8, 6).
dpi : int, optional
    Figure DPI. Default is 100.
save_path : str, optional
    Path to save the figure. If None, figure is not saved.
show_grid : bool, optional
    Show vertical lines at grid point positions. Default is True.
grid_line_color : str, optional
    Color of grid lines. Default is "black".
grid_line_style : str, optional
    Line style of grid lines. Default is "--".
grid_line_alpha : float, optional
    Alpha transparency of grid lines. Default is 0.7.
show_colorbar : bool, optional
    Show colorbar. Default is True.
left_title : str, optional
    Title for left panel (raw data). Default is "Raw (cell space)".
right_title : str, optional
    Title for right panel (GP-smoothed). Default is "GP-smoothed (grid space)".
suptitle : str, optional
    Figure super title. If None, auto-generated from parameters.

Returns
-------
fig : matplotlib.figure.Figure
    The matplotlib figure object.
axes : np.ndarray
    Array of matplotlib axes objects (1 or 2 elements depending on trajectory_gp).

Examples
--------
>>> # RNA/gene trajectory visualization
>>> gp_result = gp_interpolate(adata, pseudotime_col="pt")
>>> gp_mean = gp_result.uns["gp_interpolate"]["mean"].values
>>> grid_vals = gp_result.uns["gp_interpolate"]["time_grid"].values
>>> sorted_idx, _ = sort_features_by_trajectory(gp_mean, grid_vals)
>>> fig, axes = plot_trajectory_heatmap(
...     trajectory_raw=adata.X.T,
...     pseudotime=adata.obs["pt"],
...     trajectory_gp=gp_mean,
...     grid_values=grid_vals,
...     feature_order=sorted_idx
... )

>>> # Loop trajectory visualization (same interface)
>>> gp_result = gp_interpolate(loop_adata, pseudotime_col="pt")
>>> gp_mean = gp_result.uns["gp_interpolate"]["mean"].values
>>> sorted_idx, _ = sort_features_by_trajectory(gp_mean, grid_vals)
>>> fig, axes = plot_trajectory_heatmap(
...     trajectory_raw=loop_adata.X.T,
...     pseudotime=loop_adata.obs["pt"],
...     trajectory_gp=gp_mean,
...     grid_values=grid_vals,
...     feature_order=sorted_idx
... )


**Source:** Line 613 in [interpolate.py](interpolate.py#L613)


---


### `sort_features_by_trajectory(trajectory, grid_values, sort_keys=('peak_time', 'com', '-peak_value'))`


Sort features based on trajectory dynamics characteristics.

Parameters
----------
trajectory : np.ndarray
    Trajectory matrix of shape (n_features, n_grid_points).
grid_values : np.ndarray
    Grid point values of shape (n_grid_points,).
sort_keys : tuple of str
    List of sort keys in priority order. Supported keys:
    - "peak_time": continuous peak time (parabolic interpolation)
    - "com": center of mass
    - "peak_value": peak value at continuous peak time
    - "-peak_value": negative peak value (for descending sort)
    - "max_grid_idx": discrete grid index of maximum (fallback)

Returns
-------
sorted_indices : np.ndarray
    Indices that would sort the features.
sort_values : pd.DataFrame
    DataFrame with columns for each computed metric.


**Source:** Line 274 in [interpolate.py](interpolate.py#L274)


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

**Public functions:** 1


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

**Public functions:** 1


### `mt(func)`


No documentation available.


**Source:** Line 22 in [paracalc.py](paracalc.py#L22)


---


## phasing

**File:** `phasing.py`

**Public functions:** 11


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


**Source:** Line 68 in [phasing.py](phasing.py#L68)


---


### `count_chrom_phased(filename: str, *args)`

**Returns:** `pd.Series`


Count number of phased (with separable SNP) segments per chromosome from a .seg.gz file.
Input:
    filename: Path to the .seg.gz file.
Output:
    A pandas Series with multi-index (chromosome, phasing) and counts.


**Source:** Line 492 in [phasing.py](phasing.py#L492)


---


### `count_phased_seg_perbin(filename, genome='mm10', binsize=1000000.0)`


Count the number of phased segments per bin of a given genome.
Input:
    filename: str, path to the .seg file
    genome: str, genome name (default: "mm10")
    binsize: int, size of the bins (default: 1e6)
Output:
    pandas DataFrame with columns: "chrom", "start", "phase", "count"


**Source:** Line 510 in [phasing.py](phasing.py#L510)


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


**Source:** Line 411 in [phasing.py](phasing.py#L411)


---


### `get_chrom_hap_score(dump_dir)`


Get chrom_hap_score of samples.
Input:
    dump_dir: rd/dump generated by hires_utils
Output:
    pd.Series


**Source:** Line 475 in [phasing.py](phasing.py#L475)


---


### `hap_score_gb(df)`


No documentation available.


**Source:** Line 472 in [phasing.py](phasing.py#L472)


---


### `mt_count_chrom_phased(filename: str)`

**Returns:** `pd.Series`


No documentation available.


**Source:** Line 507 in [phasing.py](phasing.py#L507)


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


**Source:** Line 21 in [phasing.py](phasing.py#L21)


---


### `read_ref_snp_file(ref_snp_file)`


Read a reference SNP file, supporting both TSV and parquet formats.

Args:
    ref_snp_file (str): Path to the reference SNP file. Can be tab-delimited (.tsv, .txt, etc.)
                       or parquet format (.parquet).

Returns:
    pandas.DataFrame: A DataFrame with columns 'chrom', 'pos', 'input_ref', 'input_alt'.


**Source:** Line 288 in [phasing.py](phasing.py#L288)


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


**Source:** Line 316 in [phasing.py](phasing.py#L316)


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


**Source:** Line 136 in [phasing.py](phasing.py#L136)


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


### `check_input(file_pat, **kwargs)`


No documentation available.


**Source:** Line 42 in [pipeline/rule.py](pipeline/rule.py#L42)


---


### `check_key(file_pat, key_to_check, **kwargs)`


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


### `safe_input(file_pat, with_key=False, **kwargs)`


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


## plot.gene_track_parser

**File:** `plot/gene_track_parser.py`

**Public functions:** 2


### `filter_genes_in_region(genes: Dict[str, Gene], seqid: str, start: int, end: int, padding: int=0)`

**Returns:** `Dict[str, Gene]`


Filter genes that overlap with a genomic region.

Parameters
----------
genes : Dict[str, Gene]
    Dictionary of genes.
seqid : str
    Chromosome/sequence ID.
start : int
    Region start position.
end : int
    Region end position.
padding : int, optional
    Padding to add to region boundaries. Default is 0.

Returns
-------
Dict[str, Gene]
    Filtered dictionary of genes.


**Source:** Line 392 in [plot/gene_track_parser.py](plot/gene_track_parser.py#L392)


---


### `parse_annotation_file(filepath: Union[str, Path])`

**Returns:** `Dict[str, Gene]`


Parse GTF or GFF3 annotation file.

Parameters
----------
filepath : str or Path
    Path to annotation file.

Returns
-------
Dict[str, Gene]
    Dictionary of gene_id -> Gene objects.


**Source:** Line 353 in [plot/gene_track_parser.py](plot/gene_track_parser.py#L353)


---


## plot.general

**File:** `plot/general.py`

**Public functions:** 6


### `plot_colorbar(vmin, vmax, colorscale='Viridis', title='TITLE', **colorbar_kwargs)`


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


### `plot_ols(data, x, y, title, **kwargs)`


Plot linear regression of dataframe


**Source:** Line 281 in [plot/general.py](plot/general.py#L281)


---


### `plot_points(array, **args)`


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

**Public functions:** 22


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


**Source:** Line 653 in [plot/hic.py](plot/hic.py#L653)


---


### `format_ticks(ax, x=True, y=True, rotate=True)`


No documentation available.


**Source:** Line 1111 in [plot/hic.py](plot/hic.py#L1111)


---


### `merge_track_data(clr, track_file, region)`


Pick out the track data that is in the region.
Input:
    clr: cooler obj.
    track_file: path to track file
    region: genome region to plot.
Output:
    pd.Series; index: start of each bin, value: track value.


**Source:** Line 374 in [plot/hic.py](plot/hic.py#L374)


---


### `pileup_IS(coolp, IS, genome='hg19', flank=800000, **args)`


Pileup TAD boundaries from a cool file.
Input:
    coolp: cooler path(adding resolutions when mcool)
    IS: Insulation score; pd.Dataframe
    genome: reference genome. eg. hg19
    flank: flank size
Output:
    pd.Dataframe


**Source:** Line 1083 in [plot/hic.py](plot/hic.py#L1083)


---


### `plot_CM(cm, grid=True, **kwargs)`


Plot contact matrix.
Input:
    cm: pd.DataFrame; contact matrix.
    grid: whether to add chrom boundaries.


**Source:** Line 199 in [plot/hic.py](plot/hic.py#L199)


---


### `plot_IS(IS_file)`


No documentation available.


**Source:** Line 1015 in [plot/hic.py](plot/hic.py#L1015)


---


### `plot_IS(clr, insulation_table, title, resolution=500000.0, window=1500000.0, balance=True, chrom='chr2', start=10500000, steps=90, vmin=0.001, vmax=0.1)`


Plot insulation score together with Hi-C matrix strata.
Input:
    clr: cooler obj.
    insulation_table: insulation score table.


**Source:** Line 1120 in [plot/hic.py](plot/hic.py#L1120)


---


### `plot_cdps(adata, index_col='velocity_pseudotime', reverse=False)`


Input:
    index_col: wich adata.obs col to sort sample by
    reverse: whether to reverse the order of the sample


**Source:** Line 1190 in [plot/hic.py](plot/hic.py#L1190)


---


### `plot_compare_cool_pixels(coolpA, coolpB, region, outline_pixels, subplot_titles=['A', 'B'])`


No documentation available.


**Source:** Line 589 in [plot/hic.py](plot/hic.py#L589)


---


### `plot_compare_cool_track_hor(coolp1, coolp2, IS_file1, IS_file2, eigs_file1, eigs_file2, subplot_titles, region, title, balance=False, winSize=500000)`


Plot cooler with additional track files


**Source:** Line 456 in [plot/hic.py](plot/hic.py#L456)


---


### `plot_compartment(coolp, eigs_file, region, title, eigen_col='E1', strip=False, balance=False, mask_eig_na=True, give_mat=False, quant=0.005, minmax='min', cmap='RdBu_r', eq_hist=False, **args)`


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


**Source:** Line 713 in [plot/hic.py](plot/hic.py#L713)


---


### `plot_compartment_strength(adata, index_col='velocity_pseudotime')`


No documentation available.


**Source:** Line 1256 in [plot/hic.py](plot/hic.py#L1256)


---


### `plot_compartment_trapezoid(coolp, vec_file, region, balance=True, strata=20, title='TITLE', quant=0.01, minmax='min', cmap='RdBu_r', eigen_col='E1', fillna=False, strip=False, mask_eig_na=True, add_track=False, **kwargs)`


Plot compartment heatmap as a trapezoid shape by rotating the square contact matrix 45 degrees.

This function extracts the compartment signal from a cooler file and eigenvector track,
rotates the contact matrix by 45 degrees, and displays a trapezoidal region showing the
compartment pattern.

Parameters
----------
coolp : str or Path
    Path to cooler file.
vec_file : str or Path
    Path to eigenvector file (e.g., from cooltools compartment calling).
region : str
    Genome region to plot, e.g., "chr1" or "chr1:1000000-2000000".
balance : bool, optional
    Whether to use balanced cooler matrix. Default is True.
strata : int, optional
    Number of rows to keep in the trapezoid (height of the trapezoid). Default is 20.
title : str, optional
    Plot title. Default is "TITLE".
quant : float, optional
    Quantile for color scale clipping. Default is 0.01.
minmax : str or None, optional
    How to compute color limits. Default is "min".
cmap : str, optional
    Colormap name. Default is "RdBu_r".
eigen_col : str, optional
    Column name for eigenvector values. Default is "E1".
fillna : bool, optional
    Whether to fill NaN values in the matrix. Default is False.
strip : bool, optional
    Whether to strip consecutive zeros. Default is False.
mask_eig_na : bool, optional
    Whether to mask NaN eigenvalues. Default is True.
add_track : bool, optional
    Whether to add eigenvector track below the heatmap. Default is False.
**kwargs : dict
    Additional keyword arguments passed to go.Heatmap.

Returns
-------
plotly.graph_objects.Figure
    Plotly figure with trapezoid compartment heatmap.

Examples
--------
>>> fig = plot_compartment_trapezoid(
...     coolp="sample.cool",
...     vec_file="eigvec.bedgraph",
...     region="chr1",
...     strata=50,
...     title="chr1 compartment"
... )
>>> fig.show()


**Source:** Line 1309 in [plot/hic.py](plot/hic.py#L1309)


---


### `plot_compartment_trapezoid_2tracks(coolp0: Union[str, Path], vec_file1: Union[str, Path], vec_file2: Union[str, Path], region: str, balance: bool=True, strata: int=20, title: str='TITLE', quant: float=0.01, minmax: str='min', cmap: str='RdBu_r', eigen_col: str='E1', fillna: bool=False, strip: bool=False, mask_eig_na: bool=True, vec_file0: Optional[Union[str, Path]]=None, coolp1: Optional[Union[str, Path]]=None, coolp2: Optional[Union[str, Path]]=None)`

**Returns:** `go.Figure`


Plot compartment trapezoid heatmap with two eigenvector tracks.

This function creates a 3-row subplot layout using the Track architecture:
- Row 1: Trapezoid compartment heatmap
- Row 2: Eigenvector track 1 (e.g., sample 1)
- Row 3: Eigenvector track 2 (e.g., sample 2 or reference)

Args:
    coolp0: Path to primary cooler file.
    vec_file1: Path to first eigenvector file.
    vec_file2: Path to second eigenvector file.
    region: Genome region to plot (e.g., "chr1").
    balance: Whether to use balanced cooler matrix. Default True.
    strata: Number of rows to keep in the trapezoid. Default 20.
    title: Plot title. Default "TITLE".
    quant: Quantile for color scale clipping. Default 0.01.
    minmax: How to compute color limits. Default "min".
    cmap: Colormap name. Default "RdBu_r".
    eigen_col: Column name for eigenvector values. Default "E1".
    fillna: Whether to fill NaN values in the matrix. Default False.
    strip: Whether to strip consecutive zeros. Default False.
    mask_eig_na: Whether to mask NaN eigenvalues. Default True.
    vec_file0: Path to eigenvector file for primary cooler. Default None (uses vec_file1).
    coolp1: Path to second cooler file for eigenvector 1. Default None (uses coolp0).
    coolp2: Path to third cooler file for eigenvector 2. Default None (uses coolp0).

Returns:
    go.Figure: Plotly figure with trapezoid heatmap and two eigenvector tracks.

Example:
    >>> fig = plot_compartment_trapezoid_2tracks_v2(
    ...     coolp0="sperm.cool",
    ...     vec_file1="sperm_eigvec.bedgraph",
    ...     vec_file2="esc_eigvec.bedgraph",
    ...     region="chr1",
    ...     strata=50,
    ...     title="Sperm vs ESC compartment"
    ... )
    >>> fig.show()


**Source:** Line 1499 in [plot/hic.py](plot/hic.py#L1499)


---


### `plot_cool(coolp, title='', region='chr1', vmax=100, balance=False, ignore_diags=False, grid=True, donorm=False, return_data=False, **kwargs)`


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
    return_data: if True, return (fig, matrix, bins) instead of just fig.
        Default False.
Output:
    If return_data=False: plotly figure
    If return_data=True: (fig, matrix DataFrame, bins DataFrame)


**Source:** Line 238 in [plot/hic.py](plot/hic.py#L238)


---


### `plot_cool_track(coolp, track_files, region, title, balance=False, **args)`


Plot cooler with additional track files.
Input:
    coolp: path to cooler file.
    track_files: dict; key: track name, value: (path to track file, column name in track file)
    region: genome region to plot.
    title: name of the plot.
    balance: whether to load balanced cooler matrix.


**Source:** Line 395 in [plot/hic.py](plot/hic.py#L395)


---


### `plot_cool_with_tracks(coolp1: Union[str, Path], coolp2: Optional[Union[str, Path]]=None, tracks: Optional[List]=None, region: str='chr1', track_heights: Optional[List[float]]=None, heatmap_height: float=10, square_heatmap: bool=False, title: str='', balance: bool=True, vmax: Optional[float]=None, cmap: str='fall', ignore_diags: bool=False, adjust_brightness: bool=True, marks: Optional[List[Dict]]=None, **heatmap_kwargs)`

**Returns:** `go.Figure`


Plot Hi-C contact heatmap with optional tracks below.

This function creates a vertical layout with:
- Row 1: Hi-C contact matrix heatmap (single sample or tiling of two samples)
- Rows 2 to N: Additional tracks (BigWig, Gene, Eigen, IS, etc.)

When two cooler files are provided, the function uses tiling layout to
compare two contact matrices: lower triangle from coolp1, upper triangle
from coolp2.

Parameters
----------
coolp1 : str or Path
    Path to first cooler file.
coolp2 : str or Path, optional
    Path to second cooler file for tiling comparison. If None, only
    coolp1 is displayed. Default is None.
tracks : list, optional
    List of Track objects (BigWigTrack, GeneTrack, EigenTrack, etc.).
    Default is None (no tracks).
region : str, optional
    Genome region to plot, e.g., "chr1" or "chr1:1000000-5000000".
    Default is "chr1".
track_heights : list of float, optional
    Relative height ratios for each track. If None, all tracks have
    equal height of 2.0. Default is None.
heatmap_height : float, optional
    Relative height of heatmap row. Default is 10.
square_heatmap : bool, optional
    If True, constrain the heatmap x/y axes to equal data-unit scaling
    using Plotly scale anchors. This avoids row-height heuristics while
    preserving the same genomic x-axis domain for tracks. Default is False.
title : str, optional
    Plot title. Default is "".
balance : bool, optional
    Whether to use balanced cooler matrix. Default is True.
vmax : float, optional
    Maximum value for color scale. If None, auto-scaled. Default is None.
cmap : str, optional
    Colormap name. Default is "fall".
ignore_diags : bool, optional
    Whether to ignore diagonal values. Default is False.
adjust_brightness : bool, optional
    Whether to adjust brightness when tiling two matrices. When True,
    scales the brighter matrix to match the darker one. Default is True.
    Only used when coolp2 is provided.
marks : list of dict, optional
    List of mark dictionaries to highlight specific pixels. Each dict
    should contain:
    - chrom1, start1, end1: First genomic position
    - chrom2, start2, end2: Second genomic position
    - color: Optional, color of the mark (default: red)
    - shape: Optional, shape type ("circle", "rect") (default: circle)
    - width: Optional, line width (default: 2)
    - label: Optional, text label for the mark
    Default is None (no marks).
**heatmap_kwargs
    Additional keyword arguments passed to go.Heatmap.

Returns
-------
plotly.graph_objects.Figure
    Plotly figure with Hi-C heatmap and tracks.

Examples
--------
Single sample with tracks:

>>> from hic_basic.plot import BigWigTrack, GeneTrack
>>> from hic_basic.plot.hic import plot_cool_with_tracks
>>>
>>> fig = plot_cool_with_tracks(
...     coolp1="sample.cool",
...     tracks=[
...         BigWigTrack("H3K27ac.bw", bin_size=1000),
...         BigWigTrack("H3K4me3.bw", bin_size=1000),
...         GeneTrack("annotations.gtf")
...     ],
...     region="chr1:1000000-5000000",
...     title="Hi-C + Tracks"
... )
>>> fig.show()

Two samples tiling comparison with marks:

>>> fig = plot_cool_with_tracks(
...     coolp1="sample1.cool",
...     coolp2="sample2.cool",
...     marks=[
...         {"chrom1": "chr1", "start1": 100000, "end1": 110000,
...          "chrom2": "chr1", "start2": 200000, "end2": 210000,
...          "color": "red", "shape": "circle", "width": 2}
...     ],
...     region="chr1:0-500000",
...     title="Sample1 vs Sample2 (tiling)",
...     adjust_brightness=True
... )
>>> fig.show()

Notes
-----
- When coolp2 is provided, uses tiling_mat() to create a composite matrix
  with lower triangle from coolp1 and upper triangle from coolp2.
- All tracks are aligned to coolp1's bins.
- The adjust_brightness parameter helps visualize differences when two
  matrices have very different signal intensities.
- Marks are rendered as shapes (circles, rectangles) on top of the heatmap.


**Source:** Line 2171 in [plot/hic.py](plot/hic.py#L2171)


---


### `plot_cools(cools, regions, ncols=3, subplot_titles=None, sub_height=100, sub_width=100, margin=None, **kwargs)`


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


**Source:** Line 308 in [plot/hic.py](plot/hic.py#L308)


---


### `plot_saddle(file, title, vmin=10 ** (-1), vmax=10 ** 1, strength_square_size=10, **kwargs)`


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


**Source:** Line 941 in [plot/hic.py](plot/hic.py#L941)


---


### `plot_saddle_mpl(file, title, vmin=10 ** (-1), vmax=10 ** 1)`


No documentation available.


**Source:** Line 927 in [plot/hic.py](plot/hic.py#L927)


---


### `plot_tiling_compartment(coolps: List[Union[str, Path]], eigs_files: List[Union[str, Path]], region: str, title: str='', corr: bool=True, strip: bool=True, balance: bool=False, cmap: str='RdBu', donorm: bool=False, Enames: List[str]=None, **args)`

**Returns:** `go.Figure`


Plot two contact matrices in a tiling layout with eigenvector tracks.

This function creates a 2x2 subplot layout where:
- Top-left: Eigenvector track (vertical)
- Top-right: Combined heatmap (upper triangle from matrix 2, lower from matrix 1)
- Bottom-right: Eigenvector track (horizontal)

Args:
    coolps: List of two cooler file paths.
    eigs_files: List of two eigenvector file paths.
    region: Genomic region to plot (e.g., "chr1:1000000-5000000").
    title: Figure title.
    corr: Whether to plot correlation matrix. Default True.
    strip: Whether to strip consecutive zero rows/columns. Default True.
    balance: Whether to use balanced cooler matrix. Default False.
    cmap: Colormap name. Default "RdBu".
    donorm: Whether to use log normalization. Default False.
    Enames: Column names for eigenvector values. Default ["E1", "E1"].
    **args: Additional arguments passed to _plot_mat.

Returns:
    go.Figure: Plotly figure with tiling compartment plot.

Example:
    >>> fig = plot_tiling_compartment(
    ...     coolps=["sample1.cool", "sample2.cool"],
    ...     eigs_files=["eigs1.bedgraph", "eigs2.bedgraph"],
    ...     region="chr1:1000000-5000000",
    ...     title="Compartment Comparison"
    ... )
    >>> fig.show()


**Source:** Line 1673 in [plot/hic.py](plot/hic.py#L1673)


---


### `plot_trapezoid_with_tracks(coolp: Union[str, Path], vec_file: Union[str, Path], region: str, tracks: List, track_heights: Optional[List[float]]=None, trapezoid_height: Optional[float]=None, strata: int=20, title: str='', balance: bool=True, cmap: str='RdBu_r', quant: float=0.01, minmax: str='min', eigen_col: str='E1', fillna: bool=False, strip: bool=False, mask_eig_na: bool=True, total_height: int=800, width: int=700, auto_scale: bool=True, **heatmap_kwargs)`

**Returns:** `go.Figure`


Plot trapezoid Hi-C heatmap with arbitrary tracks below.

This function creates a vertical layout with:
- Row 1: Trapezoid compartment heatmap (rotated 45 degrees)
- Rows 2 to N: Additional tracks (BigWig, Gene, Eigen, IS, etc.)

Parameters
----------
coolp : str or Path
    Path to cooler file.
vec_file : str or Path
    Path to eigenvector file for compartment coloring.
region : str
    Genome region to plot, e.g., "chr1:1000000-5000000".
tracks : list
    List of Track objects (BigWigTrack, GeneTrack, EigenTrack, etc.).
track_heights : list of float, optional
    Relative height ratios for each track. If None, all tracks have
    equal height of 2.0.
trapezoid_height : float, optional
    Relative height of trapezoid row. If auto_scale=True, this parameter
    is ignored and height is calculated automatically. Set auto_scale=False
    to use a custom height. Default None (auto-calculated).
strata : int, optional
    Number of rows to keep in the trapezoid. Default 20.
title : str, optional
    Plot title. Default "".
balance : bool, optional
    Whether to use balanced cooler matrix. Default True.
cmap : str, optional
    Colormap name. Default "RdBu_r".
quant : float, optional
    Quantile for color scale clipping. Default 0.01.
minmax : str, optional
    How to compute color limits. Default "min".
eigen_col : str, optional
    Column name for eigenvector values. Default "E1".
fillna : bool, optional
    Whether to fill NaN values in the matrix. Default False.
strip : bool, optional
    Whether to strip consecutive zeros. Default False.
mask_eig_na : bool, optional
    Whether to mask NaN eigenvalues. Default True.
total_height : int, optional
    Total figure height in pixels. Default 800.
width : int, optional
    Figure width in pixels. Default 700.
auto_scale : bool, optional
    If True, automatically calculate trapezoid_height to maintain
    square pixel aspect ratio. Default True.
**heatmap_kwargs
    Additional keyword arguments passed to go.Heatmap.

Returns
-------
plotly.graph_objects.Figure
    Plotly figure with trapezoid heatmap and tracks.

Examples
--------
>>> from hic_basic.plot.track import BigWigTrack, GeneTrack
>>> from hic_basic.plot.hic import plot_trapezoid_with_tracks
>>>
>>> tracks = [
...     BigWigTrack("H3K27ac.bw", bin_size=1000, color="#1f77b4"),
...     BigWigTrack("H3K4me3.bw", bin_size=1000, color="#ff7f0e"),
...     GeneTrack("annotations.gtf", mode="all_transcripts")
... ]
>>> fig = plot_trapezoid_with_tracks(
...     coolp="sample.cool",
...     vec_file="eigvec.bedgraph",
...     region="chr1:1000000-5000000",
...     tracks=tracks,
...     track_heights=[2, 2, 3],  # gene track taller
...     title="Hi-C + Multi-omics"
... )
>>> fig.show()

Notes
-----
The function uses Track architecture for modular composition.
Track types supported:
- BigWigTrack: Continuous signal from bigWig files
- GeneTrack: Gene annotation from GTF/GFF3 files
- EigenTrack: Eigenvector compartment signal
- ISTrack: Insulation score track
- HeatmapTrack: Standard Hi-C heatmap

Auto-scaling behavior:
- When auto_scale=True, trapezoid_height is calculated to maintain
  square pixel aspect ratio based on the trapezoid matrix shape
- Formula accounts for strata, n_bins, figure dimensions, and track space


**Source:** Line 1980 in [plot/hic.py](plot/hic.py#L1980)


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


### `plot_cdps_mark(cdps, orig_annote, sample_col='sample_name', order_col='index_order', group_col='group', color_map=None, hline=False, **kwargs)`


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


## plot.plotly_template

**File:** `plot/plotly_template.py`

**Public functions:** 1


### `apply_template(fig, title=None, width=None, height=None)`


应用 hic_basic 模板到图形

Parameters
----------
fig : plotly.graph_objects.Figure
    Plotly 图形对象
title : str, optional
    标题（可选覆盖）
width : int, optional
    宽度（可选覆盖）
height : int, optional
    高度（可选覆盖）

Returns
-------
fig : plotly.graph_objects.Figure
    应用模板后的图形


**Source:** Line 102 in [plot/plotly_template.py](plot/plotly_template.py#L102)


---


## plot.render

**File:** `plot/render.py`

**Public functions:** 21


### `apply_rotation_commands(df, commands_str)`


No documentation available.


**Source:** Line 166 in [plot/render.py](plot/render.py#L166)


---


### `centelo_relpos(positions_2col, genome, dupref=False)`


No documentation available.


**Source:** Line 477 in [plot/render.py](plot/render.py#L477)


---


### `clip_b_pymol(_3dg, b_factor, png, clip=0, slab=2, cmap='magenta green, all, 0.005, 0.02', turn_cmd='', tmpdir=None, conda='pymol', **args)`


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


### `clip_centelo_pymol(_3dg_file, png, genome='mm10', clip=0, slab=2, tmpdir=None, cif_name=None, dupref=False, conda='pymol', turn_cmd='', **args)`


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


### `clip_single_centelo_pymol(_3dg_file, png, target_chroms=['chrX', 'chrY'], genome='mm10', clip=0, slab=2, tmpdir=None, cif_name=None, dupref=False, conda='pymol', turn_cmd='', **args)`


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


### `clip_single_territory_pymol(_3dg_file, png, target_chroms=['chrX', 'chrY'], clip=0, slab=2, tmpdir=None, conda=None, turn_cmd='', **args)`


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


### `clip_territory_pymol(_3dg, png, clip=0, slab=2, tmpdir=None, conda='pymol', turn_cmd='', **args)`


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


### `get_rotation_commands(target_vector, retvec=False)`


Calculate PyMOL turn commands to rotate the camera to face a given direction.

Input:
    target_vector: A list or numpy array representing the unit vector towards which the camera should be oriented.
    retvec: If True, return the rotation degrees along x, y, z axis as well.
Output:
    A list of string containing a sequence of 'turn' commands for PyMOL to adjust the camera's orientation accordingly.


**Source:** Line 99 in [plot/render.py](plot/render.py#L99)


---


### `highlight_surface_b_pymol(_3dg, b_factor, chain, png, cmap='magenta green, chain {}, 0.005, 0.02', turn_cmd='', tmpdir=None, conda='pymol', **args)`


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


### `render_clip_b_primary_view(_3dg_file, b_factor, outpng, view='lr', targets=None, **kwargs)`


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


### `render_clip_centelo_primary_view(_3dg_file, outpng, genome, clip, slab, tmpdir, view='lr', targets=None, **kwargs)`


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


### `render_clip_primary_view(_3dg_file, outpng, view='lr', targets=None, **kwargs)`


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


### `render_surface_centelo_primary_view(_3dg_file, outpng, genome, tmpdir, view='lr', targets=None, **kwargs)`


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


### `render_surface_primary_view(_3dg_file, outpng, tmpdir, view='lr', targets=None, **kwargs)`


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


### `surface_b_pymol(_3dg, b_factor, png, cmap='magenta green, all, 0.005, 0.02', tmpdir=None, conda='pymol', turn_cmd='', **args)`


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


### `surface_centelo_pymol(_3dg_file, png, genome='mm10', tmpdir=None, cif_name=None, dupref=False, conda='pymol', turn_cmd='', **args)`


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


**Source:** Line 258 in [plot/scanpy.py](plot/scanpy.py#L258)


---


### `plot_elbow(adata)`


No documentation available.


**Source:** Line 215 in [plot/scanpy.py](plot/scanpy.py#L215)


---


### `plot_pca(adata, color: str, PCs: Optional[Union[List[int], List[List[int]]]]=None, ncols: int=2, background: Optional[str]=None, **kwargs)`


Plot PCA results from AnnData object with flexible subplot configurations.

This function creates scatter plots of PCA components from single-cell data.
It can generate single plots or multiple subplots showing different combinations
of principal components.

Parameters
----------
adata : AnnData
    AnnData object containing PCA results in ``obsm['X_pca']``.
color : str
    Column name in ``adata.obs`` to use for coloring foreground points.
PCs : list of int or list of list of int, optional
    Specification of which principal components to plot:

    - ``None``: plots first two PCs (PC1 vs PC2).
    - ``List[int]``: plots all pairwise combinations of the specified PCs.
    - ``List[List[int]]``: plots specific pairs of PCs as subplots.
ncols : int, optional
    Number of columns for subplot layout when multiple plots are generated
    (default 2).
background : str, optional
    A :meth:`pandas.DataFrame.query` expression (same syntax as
    ``DataFrame.query``) that selects samples to treat as *background*.
    Background samples are rendered in light grey and drawn behind the
    foreground points.  Foreground samples (those **not** matching the
    expression) are colored normally according to *color*.  Defaults to
    ``None`` (all samples are foreground).
**kwargs
    Additional arguments passed to ``px.scatter`` for foreground traces.

Returns
-------
plotly.graph_objects.Figure
    Figure object containing the PCA plot(s).

Examples
--------
>>> # Plot first two PCs
>>> fig = plot_pca(adata, color='cell_type')

>>> # Plot first two PCs with a grey background group
>>> fig = plot_pca(adata, color='cell_type', background="batch == 'B'")

>>> # Plot all pairwise combinations of PCs 1-3
>>> fig = plot_pca(adata, color='cell_type', PCs=[0, 1, 2])

>>> # Plot specific PC pairs in subplots
>>> fig = plot_pca(adata, color='cell_type', PCs=[[0, 1], [0, 2], [1, 3]])


**Source:** Line 14 in [plot/scanpy.py](plot/scanpy.py#L14)


---


### `plot_tsne(adata, color)`


Plot scanpy tsne using plotly.


**Source:** Line 242 in [plot/scanpy.py](plot/scanpy.py#L242)


---


### `plot_umap(adata, color, key='X_umap', **kwargs)`


Plot scanpy umap using plotly.


**Source:** Line 224 in [plot/scanpy.py](plot/scanpy.py#L224)


---


## plot.scatter

**File:** `plot/scatter.py`

**Public functions:** 5


### `add_error_band(fig, df, orig_color=px.colors.qualitative.Plotly[0], debug=False, error_type='ci')`


Add error visualization to a figure.

Parameters
----------
fig : plotly.graph_objects.Figure
    Figure to add error to.
df : pandas.DataFrame
    DataFrame with x, y_upper, y_lower columns.
orig_color : str
    Color of original line, will use a grayer version for the ribbon.
debug : bool
    If True, show ribbon outline.
error_type : str
    "ci" for pooled CI ribbon, "sem" for varying-width ribbon from
    pre-computed SEM values.

Returns
-------
plotly.graph_objects.Figure
    Figure with error added.


**Source:** Line 145 in [plot/scatter.py](plot/scatter.py#L145)


---


### `line(df, x, y, color=None, color_discrete_map=None, with_error=False, frac=0.3, ci=0.95, debug=False, gray_bg=False, hover_data=None, hover_format=None, skip_smooth=False, error_type='ci', sem_col=None, **kwargs)`


Plot a single trend line with optional error bands and point coloring.

Creates a line plot for a single variable with optional confidence intervals.
Points can be colored according to a specified column in the dataframe.

Args:
    df (pd.DataFrame): Input curve data with one row per x value.
    x (str): Column name for x-axis values.
    y (str): Column name for y-axis values. Must already be one-to-one with x.
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
    hover_data (list, optional): Additional columns to show in hover tooltip.
    hover_format (dict, optional): Format specifiers for hover columns.
    skip_smooth (bool, optional): If True, bypass LOESS and plot raw y values.
                                  Defaults to False.
    error_type (str, optional): "ci" for pooled CI ribbon, "sem" for
                                ribbon from pre-computed SEM values at
                                each x. Defaults to "ci".
    sem_col (str, optional): Column name with pre-computed SEM values for
                             each x value. Required when error_type="sem".
    **kwargs: Additional keyword arguments passed to fig.update_layout().

Returns:
    plotly.graph_objs.Figure: Plotly figure object.

Notes:
    line() does not aggregate repeated observations. If your raw data has
    multiple rows per x value, aggregate it upstream so that each x maps to
    exactly one y, then pass SEM explicitly through sem_col if needed.

Examples:
    >>> line(curve_df, x='time', y='value', color='group')
    >>> line(curve_df, x='x', y='y', with_error=True, color_discrete_map={'A': 'red', 'B': 'blue'})
    >>> line(agg_df, x='bin', y='mean', with_error=True, error_type="sem", sem_col="sem")


**Source:** Line 316 in [plot/scatter.py](plot/scatter.py#L316)


---


### `lines(df, x, y, color_discrete_map=None, with_error=False, frac=0.3, ci=0.95, debug=False, gray_bg=False, skip_smooth=False, error_type='ci', sem_col=None, **kwargs)`


Plot multiple line trends with optional error visualization.

Args:
    df (pd.DataFrame): Input curve data with one row per x value.
    x (str): Column name for x-axis values.
    y (str or list): Column name(s) for y-axis values. Each y series must
        already be aligned one-to-one with x.
    color_discrete_map (dict, optional): Color mapping for each y column.
    with_error (bool): Whether to add error visualization.
    frac (float): Fraction of data for LOESS smoothing.
    ci (float): Confidence interval level (0-1).
    debug (bool): Enable debug mode for error bands.
    gray_bg (bool): Use gray color for error bands.
    skip_smooth (bool): If True, bypass LOESS and plot raw y values.
    error_type (str): "ci" for pooled CI ribbon, "sem" for ribbon from
        pre-computed SEM values at each x.
    sem_col (str or list[str]): Column name containing pre-computed SEM
        values for each x value, or a list aligned one-to-one with y when
        y is a list. Required when error_type="sem".
    **kwargs: Additional arguments passed to fig.update_layout().

Returns:
    plotly.graph_objects.Figure

Notes:
    lines() does not aggregate repeated observations. If your raw data has
    multiple rows per x value, compute summary statistics upstream and pass
    the aggregated curve plus an optional sem_col.
    When y is a list and error_type="sem", sem_col may be a single shared
    SEM column name or a list aligned one-to-one with y.


**Source:** Line 224 in [plot/scatter.py](plot/scatter.py#L224)


---


### `rgb_grayer(color, white_ratio=0.7, opacity=0.3)`


Convert a color to a lighter, more transparent version by blending with white
and reducing opacity.

Parameters:
    color: A list of three integers representing RGB values (0-255).

Returns:
    A string in "rgba(r, g, b, a)" format with adjusted color and opacity.


**Source:** Line 121 in [plot/scatter.py](plot/scatter.py#L121)


---


### `smooth_with_ci(df, x, y, frac=0.3, ci=0.95, skip_smooth=False, error_type='ci', sem_col=None)`


Apply LOESS smoothing to one-row-per-x curve data and estimate ribbons.

Parameters
----------
df : pandas.DataFrame
    Input DataFrame containing one y value for each x value. Repeated
    observations must be aggregated upstream before calling this function.
x : str
    Column name for independent variable (e.g., time points).
y : str
    Column name for dependent variable (e.g., aggregated mean values).
frac : float, optional
    Fraction of data used for each local regression (0 < frac <= 1).
    Default is 0.3.
ci : float, optional
    Confidence interval level (0 < ci < 1, e.g., 0.95 for 95% CI).
    Default is 0.95.
skip_smooth : bool, optional
    If True, bypass LOESS smoothing and use raw y values directly.
    Default is False.
error_type : str, optional
    Type of error estimation: "ci" for pooled residual CI, "sem" for
    per-point SEM. Default is "ci". When "sem", requires sem_col to be
    specified.
sem_col : str, optional
    Column name containing pre-computed SEM values for each x value.
    Required when error_type="sem". SEM is not computed internally.

Returns
-------
pandas.DataFrame
    DataFrame with columns:
    - 'x': Original x values (sorted)
    - 'y': Original y values
    - 'y_smooth': Smoothed y values (or raw y if skip_smooth=True)
    - 'y_upper': Upper bound of error band
    - 'y_lower': Lower bound of error band
    - 'y_sem': SEM values (if error_type="sem")

Notes
-----
This function assumes df already describes a single curve. It does not
group, bin, or aggregate repeated observations by x.


**Source:** Line 11 in [plot/scatter.py](plot/scatter.py#L11)


---


## plot.tech

**File:** `plot/tech.py`

**Public functions:** 1


### `tech_box(annote, grouping='condition', features=['defaults'], ncols=3, size=200, color_discrete_sequence=color_discrete_sequence, use_label=True, labels: dict=None, **kwargs)`

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

**Public functions:** 9


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


**Source:** Line 2307 in [plot/track.py](plot/track.py#L2307)


---


### `plot_arc_track(data: Union[str, Path, pd.DataFrame], region: str, title: Optional[str]=None, width: int=1000, height: int=260, output_path: Optional[Union[str, Path]]=None, **arc_kwargs)`

**Returns:** `go.Figure`


Plot BEDPE links as a standalone genomic arc track.

Parameters
----------
data : str, Path, or pandas.DataFrame
    BEDPE-like table accepted by :class:`ArcTrack`. The first six columns
    are chrom1, start1, end1, chrom2, start2, and end2. Optional score
    controls default arc width, and optional info supplies Plotly line
    attributes as a dict/JSON/Python-literal string.
region : str
    Genomic region to show, for example "chr1:0-5000000".
title : str, optional
    Figure title. If None, a title is generated from region.
width, height : int, optional
    Figure size in pixels. Defaults are 1000 by 260.
output_path : str or Path, optional
    PNG output path. If provided, the figure is written to this path.
**arc_kwargs
    Additional keyword arguments passed to :class:`ArcTrack`.

Returns
-------
plotly.graph_objects.Figure
    The generated arc track figure.


**Source:** Line 957 in [plot/track.py](plot/track.py#L957)


---


### `plot_gene_track(annotation_file: Union[str, Path], region: str, mode: str='all_transcripts', title: Optional[str]=None, height: int=400, width: int=1000, utr_height: float=0.3, exon_height: float=0.4, cds_height: float=0.6, row_spacing: float=0.2, show_labels: bool=True, label_font_size: int=10)`

**Returns:** `go.Figure`


Plot gene annotation track in genome browser style.

This function creates a Plotly figure showing gene structures including
exons, introns, CDS, and UTRs with strand direction indicators.

Parameters
----------
annotation_file : str or Path
    Path to GTF or GFF3 annotation file.
region : str
    Genomic region to display, e.g., "chr1:1000000-2000000".
mode : str, optional
    Display mode. Options:
    - "all_transcripts": Show all transcripts for each gene
    - "canonical": Show only canonical/longest transcript per gene
    - "gene_only": Show collapsed gene models (one row per gene)
    Default is "all_transcripts".
title : str, optional
    Plot title. If None, uses region as title.
height : int, optional
    Figure height in pixels. Default is 400.
width : int, optional
    Figure width in pixels. Default is 1000.
utr_height : float, optional
    Relative height of UTR boxes (0-1). Default is 0.30.
exon_height : float, optional
    Relative height of exon boxes (0-1). Default is 0.40.
cds_height : float, optional
    Relative height of CDS boxes (0-1). Default is 0.60.
row_spacing : float, optional
    Spacing between transcript rows. Default is 0.2.
show_labels : bool, optional
    Whether to show gene/transcript labels. Default is True.
label_font_size : int, optional
    Font size for labels. Default is 10.

Returns
-------
plotly.graph_objects.Figure
    Plotly figure with gene track visualization.

Examples
--------
>>> fig = plot_gene_track(
...     annotation_file="annotations.gtf",
...     region="chr1:1000000-2000000",
...     mode="all_transcripts",
...     title="Gene Track Example"
... )
>>> fig.show()

Notes
-----
TODO: Strand direction arrows are currently disabled pending improved
implementation. Planned improvements include:
- Fishbone-style arrows instead of solid triangles
- Proper direction rendering for both + and - strands
- Better visual integration with exon/intron structures


**Source:** Line 2979 in [plot/track.py](plot/track.py#L2979)


---


### `plot_gene_track_with_heatmap(annotation_file: Union[str, Path], region: str, heatmap_data: Optional[pd.DataFrame]=None, heatmap_title: Optional[str]=None, gene_track_height: int=300, heatmap_height: int=400, **kwargs)`

**Returns:** `go.Figure`


Plot gene track below a heatmap (e.g., Hi-C contact matrix).

Parameters
----------
annotation_file : str or Path
    Path to GTF or GFF3 annotation file.
region : str
    Genomic region to display.
heatmap_data : pd.DataFrame, optional
    Heatmap data (e.g., contact matrix). If None, only gene track is shown.
heatmap_title : str, optional
    Title for heatmap subplot.
gene_track_height : int, optional
    Height of gene track subplot. Default is 300.
heatmap_height : int, optional
    Height of heatmap subplot. Default is 400.
**kwargs
    Additional arguments passed to plot_gene_track.

Returns
-------
plotly.graph_objects.Figure
    Combined figure with heatmap and gene track.


**Source:** Line 3180 in [plot/track.py](plot/track.py#L3180)


---


### `plot_genome_heatmap(genome_data, chromosomes=None, **kwargs)`


Heatmap genome wide data. Sepaarate chromosomes by vertical lines and add chromosome names.
Input:
    genome_data: pd.DataFrame, index is multiindex with level 0 as chromosome names
        columns are samples. Will sort the index first.
    chromosomes: pd.DataFrame or Feature, the chromosome information used to sort the genome_data.
Output:
    fig: plotly.graph_objs._figure.Figure


**Source:** Line 2273 in [plot/track.py](plot/track.py#L2273)


---


### `plot_genome_phased_read_count_track(phased_read_counts, region=None, output_path: Optional[Union[str, Path]]=None, title: str='phased_read_count', count_cols: Tuple[str, ...]=('ref', 'alt', 'other'), colors: Optional[Dict[str, str]]=None, y_cap_quantile: float=0.995, width: int=1400, height: int=420)`


Plot a genome-wide phased read-count track.

This function is intended for internal phased-read count tables rather than
standard signal files such as bigWig. The input table must contain
chromosome/position coordinates and phased count columns.

Parameters
----------
phased_read_counts : pandas.DataFrame
    Table with chromosome, position, ref, alt, and other columns. The
    aliases chrom and pos are accepted. Optional start_pos and end_pos
    columns are used for interval hover text and region filtering.
region : str or tuple, optional
    Region to plot, either "chr:start-end" or (chrom, start, end). If None,
    all rows are shown with a genome-wide bin index when multiple
    chromosomes are present.
output_path : str or Path, optional
    PNG path. If provided, the figure is written to this path.
title : str, optional
    Figure title.
count_cols : tuple of str, optional
    Count columns to draw. Default is ("ref", "alt", "other").
colors : dict, optional
    Mapping from count column to trace color.
y_cap_quantile : float, optional
    Quantile used to cap the displayed y-axis. Default 0.995.
width, height : int, optional
    PNG export size in pixels when output_path is provided.

Returns
-------
plotly.graph_objects.Figure
    The generated read-count track figure.


**Source:** Line 2084 in [plot/track.py](plot/track.py#L2084)


---


### `plot_genome_scatter(genome_data, chromosomes=None, value_cols=None, shift_unit=0.05, shift_threshold=100, chorm_label_shift=0.05, y_extra_expansion=0.1, range_y=None, **kwargs)`


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


**Source:** Line 2222 in [plot/track.py](plot/track.py#L2222)


---


### `resolve_state_palette(states: List[str], palette: Union[str, Dict[str, str]]='auto')`

**Returns:** `Dict[str, str]`


Resolve state labels to colors using a built-in or custom palette.

Parameters
----------
states : list of str
    State labels in desired legend order.
palette : {"auto", "chromhmm_15", "chromhmm_18"} or dict, optional
    Named built-in palette or a custom state-to-color mapping. Missing
    labels fall back to Plotly's G10 qualitative sequence.

Returns
-------
dict
    Mapping from every input state label to a hex color.


**Source:** Line 432 in [plot/track.py](plot/track.py#L432)


---


### `sort_index_chrom(df, chromosomes)`


No documentation available.


**Source:** Line 2204 in [plot/track.py](plot/track.py#L2204)


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


**Source:** Line 405 in [plot/utils.py](plot/utils.py#L405)


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


**Source:** Line 447 in [plot/utils.py](plot/utils.py#L447)


---


### `expand_colors(colors, n_interpolate=1)`


在相邻颜色之间插入 n_interpolate 个新颜色，使颜色过渡自然。

参数:
    colors (List[tuple]): 输入颜色序列，每个颜色是 (r, g, b) 的 RGB 元组，范围 0~1。
    n_interpolate (int): 每对颜色之间插入的新颜色数量。

返回:
    List[tuple]: 扩展后的颜色序列，每个颜色为 (r, g, b) 的整数元组，范围 0~255。


**Source:** Line 478 in [plot/utils.py](plot/utils.py#L478)


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


**Source:** Line 558 in [plot/utils.py](plot/utils.py#L558)


---


### `hex_split(hex_color)`


解析十六进制颜色代码，并返回 R、G、B 和 A 分量。

:param hex_color: 十六进制颜色代码，可以是六位 (#RRGGBB) 或八位 (#RRGGBBAA)
:return: 一个元组 (R, G, B, A)，其中 A 是可选的透明度分量


**Source:** Line 377 in [plot/utils.py](plot/utils.py#L377)


---


### `hex_to_rgb(hex_color)`


将十六进制颜色字符串转换为 RGB 格式。

参数:
    hex_color (str): 十六进制颜色字符串，例如 "#4C78A8"

返回:
    tuple: 三元组 (r, g, b)，其中 r, g, b 分别为红、绿、蓝通道的值 (0-255)


**Source:** Line 343 in [plot/utils.py](plot/utils.py#L343)


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


**Source:** Line 416 in [plot/utils.py](plot/utils.py#L416)


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


**Source:** Line 265 in [plot/utils.py](plot/utils.py#L265)


---


### `pcolormesh_45deg(ax, matrix_c, start=0, resolution=1, *args, **kwargs)`


Helper to plot a matrix with 45 degree angle.


**Source:** Line 249 in [plot/utils.py](plot/utils.py#L249)


---


### `plot_color_sequence(colors, output_file=None)`


Plot a color sequence as a 1-row heatmap.

Args:
    colors (List[str] or List[tuple]): List of colors in hex or RGB format.
    output_file (str, optional): Path to save the plot as an image. If None, the plot is shown.


**Source:** Line 307 in [plot/utils.py](plot/utils.py#L307)


---


### `plotly_fig2array(fig)`


No documentation available.


**Source:** Line 548 in [plot/utils.py](plot/utils.py#L548)


---


### `rgb_to_hex(rgb_color)`


将 RGB 格式转换为十六进制颜色字符串。

参数:
    r (int): 红色通道的值 (0-255)
    g (int): 绿色通道的值 (0-255)
    b (int): 蓝色通道的值 (0-255)

返回:
    str: 十六进制颜色字符串，例如 "#4C78A8"


**Source:** Line 357 in [plot/utils.py](plot/utils.py#L357)


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


**Source:** Line 287 in [plot/utils.py](plot/utils.py#L287)


---


### `sub_genome_mat(mat, keep_regions)`


Keep only the regions in the keep_regions list.
Input:
    mat: index and columns are [chrom, start]
    keep_regions: list of regions to keep
        regions are tuples of (start, end)
Output:
    mat: matrix with only the regions in keep_regions


**Source:** Line 224 in [plot/utils.py](plot/utils.py#L224)


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


## sam

**File:** `sam.py`

**Public functions:** 7


### `convert_alignment_to_entry(alignment)`


Convert a pysam.AlignedSegment object to a list in the same format as a SAM entry.

Args:
    alignment (pysam.AlignedSegment): An alignment object from pysam.

Returns:
    list: A list representing the alignment in the same format as a SAM entry.


**Source:** Line 355 in [sam.py](sam.py#L355)


---


### `count_lines_in_file(file_path)`


Count the number of lines in a file efficiently.

Args:
    file_path (str): Path to the file.

Returns:
    int: Number of lines in the file.


**Source:** Line 387 in [sam.py](sam.py#L387)


---


### `count_reads_per_chrom(file_path: str, min_mapq: int=0, skip_unmapped: bool=True, skip_duplicates: bool=False, return_fraction: bool=False)`

**Returns:** `pd.Series`


Count reads per chromosome from a SAM or BAM file.

This function efficiently counts the number of reads mapped to each chromosome
in a SAM or BAM file, with optional filtering by mapping quality and read flags.

Args:
    file_path (str): Path to the SAM or BAM file.
    min_mapq (int): Minimum mapping quality threshold. Reads with MAPQ < min_mapq
                   will be excluded. Defaults to 0 (no filtering).
    skip_unmapped (bool): If True, skip unmapped reads. Defaults to True.
    skip_duplicates (bool): If True, skip PCR/optical duplicates. Defaults to False.
    return_fraction (bool): If True, return fractions instead of raw counts.
                           Defaults to False.

Returns:
    pandas.Series: A Series with chromosome names as index and read counts (or fractions)
                  as values, sorted by chromosome name.

Raises:
    FileNotFoundError: If the file does not exist.
    ValueError: If the file is not a valid SAM or BAM file.

Examples:
    Count reads per chromosome:

    >>> counts = count_reads_per_chrom("sample.bam")
    >>> print(counts)
    chr1     1000000
    chr2      900000
    chrX       30000
    chrY        3000
    dtype: int64

    Get fractions instead of counts:

    >>> fractions = count_reads_per_chrom("sample.bam", return_fraction=True)
    >>> print(fractions)
    chr1    0.5128
    chr2    0.4615
    chrX    0.0154
    chrY    0.0015
    dtype: float64

    Filter by mapping quality:

    >>> counts = count_reads_per_chrom("sample.bam", min_mapq=30)


**Source:** Line 401 in [sam.py](sam.py#L401)


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


**Source:** Line 283 in [sam.py](sam.py#L283)


---


### `get_file_type(file_path)`


Determine the file type by reading the file content.

Args:
    file_path (str): Path to the file to check.

Returns:
    str: The file type ('SAM', 'BAM', 'SAMgz', or 'unknown').


**Source:** Line 181 in [sam.py](sam.py#L181)


---


### `is_sam_content(text, max_lines=20)`


Check if text content matches SAM format characteristics.

Args:
    text (str): Text content to check
    max_lines (int): Maximum number of lines to check

Returns:
    bool: True if SAM format, False otherwise


**Source:** Line 248 in [sam.py](sam.py#L248)


---


### `parse_cigar(cigar: str)`

**Returns:** `List[Tuple[int, str]]`


Parse CIGAR string into operations and lengths.

Args:
    cigar: CIGAR string (e.g., "100M5D10I")

Returns:
    List of (length, operation) tuples (e.g., [(100, 'M'), (5, 'D'), (10, 'I')])


**Source:** Line 166 in [sam.py](sam.py#L166)


---


## scAB_embedding

**File:** `scAB_embedding.py`

**Public functions:** 8


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

**Public functions:** 10


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


**Source:** Line 455 in [structure/measure.py](structure/measure.py#L455)


---


### `angle_between_vectors(vector_a, vector_b)`


Calculate the angle between two 3D vectors.

Input:
    vector_a (np.array): First 3D vector.
    vector_b (np.array): Second 3D vector.
Output:
    angle_degrees (float): The angle between the two vectors in degrees.


**Source:** Line 558 in [structure/measure.py](structure/measure.py#L558)


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


**Source:** Line 414 in [structure/measure.py](structure/measure.py#L414)


---


### `append_pm_to_index(df)`


Append paternal and maternal information to the index.
Input:
    df: inner _3dg data structure (see parse_3dg in hires_utils)
Output:
    dat: data structure with paternal and maternal information
TODO: A more general way to determine paternal/maternal rather than name matching


**Source:** Line 527 in [structure/measure.py](structure/measure.py#L527)


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


**Source:** Line 477 in [structure/measure.py](structure/measure.py#L477)


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


**Source:** Line 341 in [structure/measure.py](structure/measure.py#L341)


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


**Source:** Line 366 in [structure/measure.py](structure/measure.py#L366)


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


**Source:** Line 591 in [structure/measure.py](structure/measure.py#L591)


---


### `pm_vector(_3dg)`


Compute the paternal/maternal vector.
Assume chromosome names are in the format of "chr1(pat)/chr1(mat)".
Input:
    _3dg: inner _3dg data structure (see parse_3dg in hires_utils)
Output:
    pm_vector: paternal centroid point to maternal centroid point vector (mat - pat)


**Source:** Line 541 in [structure/measure.py](structure/measure.py#L541)


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


### `consecutive_slice_pileup_bf(primary_views, targets, sample_table, features, step=2, upper=6, lower=-6, axis='ht', **args)`


No documentation available.


**Source:** Line 739 in [structure/pileup.py](structure/pileup.py#L739)


---


### `consecutive_slice_pileup_bf_parallel(primary_views, targets, sample_table, features, step=2, upper=6, lower=-6, axis='ht', prefix=None, threads=8, **args)`


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


### `pileup_bf_parallel_task(primary_views, targets, sample_table, features, step, grouping, sub, axis, prefix, key, *args)`


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

**Public functions:** 2


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


### `read_chromsizes(filepath_or, name_patterns=('^chr[0-9]+$', '^chr[XY]$', '^chrM$'), all_names=False, **kwargs)`


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


## visualize.annotations

**File:** `visualize/annotations.py`

**Public functions:** 2


### `annotate_structure(scene: StructureScene, annotations: AnnotationInput, *, layer: str='annotations', layer_label: Optional[str]=None, visible: bool=True, layer_order: int=0)`

**Returns:** `StructureScene`


Resolve and add typed annotations to an immutable structure scene.

Genome targets are resolved once in Python against the scene's stable
chromosome, particle, and polymer-bond order. The resulting scene embeds
the same resolved membership for interactive HTML and direct PNG output.
Calling the function again with the same layer id appends annotations to
that layer; annotation ids remain unique across the complete scene.

Parameters
----------
scene
    Existing validated structure scene.
annotations
    One annotation or a sequence of annotations. Genome-region targets may
    be unresolved; already resolved targets are retained and revalidated.
layer, layer_label, visible, layer_order
    Stable annotation-layer identity, display label, default visibility,
    and composition ordering key.

Returns
-------
StructureScene
    New scene containing resolved, immutable annotation membership and an
    annotation-aware checksum.


**Source:** Line 90 in [visualize/annotations.py](visualize/annotations.py#L90)


---


### `genomic_segment_annotation(region: str, *, id: str, label: Optional[str]=None, color: str='#ff0000', tube_width: float=8.0, outline_width: Optional[float]=None, outline_color: Optional[str]=None, material: str='unlit', opacity: float=1.0, occlusion: str='depth_tested', order: int=0)`

**Returns:** `StructureAnnotation`


Describe a polymer-segment annotation before scene resolution.

Parameters
----------
region
    Inclusive genome region such as ``chr1(pat):20000-30000``. A complete
    interval includes both endpoint particles and only real polymer bonds
    already present in the scene.
id, label
    Stable annotation identifier and optional display label.
color
    Requested semantic display color. It is preserved in ``unlit`` mode.
tube_width, outline_width, outline_color
    Independent annotation diameters and optional outline color in
    reference pixels at a 1400-pixel export width.
material
    ``"unlit"`` bypasses profile color transforms and lighting;
    ``"profile"`` opts into the active base-polymer profile.
opacity
    Annotation opacity in the interval ``(0, 1]``.
occlusion
    ``"depth_tested"`` for geometry-aware occlusion or the explicit
    ``"always_visible"`` overlay policy.
order
    Stable ordering key within the annotation layer.

Returns
-------
StructureAnnotation
    Unresolved typed annotation accepted by :func:`annotate_structure`.


**Source:** Line 26 in [visualize/annotations.py](visualize/annotations.py#L26)


---


## visualize.config

**File:** `visualize/config.py`

**Public functions:** 1


### `load_render_config(value: RenderConfigInput=None)`

**Returns:** `StructureRenderConfig`


Load and validate a structure render configuration.

``value`` may be an existing config, a mapping, a JSON string, a JSON file,
or ``None`` for defaults. Unknown fields and invalid values are rejected.


**Source:** Line 464 in [visualize/config.py](visualize/config.py#L464)


---


## visualize.html

**File:** `visualize/html.py`

**Public functions:** 2


### `view_structure(scene: StructureScene, *, title: Optional[str]=None, regions: Optional[Sequence[str]]=None, background: str='#151918', particle_size: float=1.0, bond_size: float=1.0, export_profile: str='molecular-studio-v1', render_config: RenderConfigInput=None)`

**Returns:** `StructureViewer`


Create a configurable offline browser viewer for a structure scene.

The returned object is lightweight and does not start a server or renderer.
Call ``write_html`` to create one portable HTML file. ``regions`` accepts
chromosome names or inclusive ranges such as ``chr1(mat):0-50Mb``; multiple
entries are displayed as a union. When supplied, ``render_config`` is the
authoritative region/profile request embedded for direct PNG rendering;
the older ``regions`` and ``export_profile`` arguments remain compatible.


**Source:** Line 122 in [visualize/html.py](visualize/html.py#L122)


---


### `write_structure_html(structure: StructureInput, output: Union[str, os.PathLike], *, color: ColorInput='chrom', palette: Optional[Mapping[Any, Any]]=None, cmap: str='viridis', vmin: Optional[float]=None, vmax: Optional[float]=None, missing_color: Any='#aeb7ba', polymer_bonds: bool=True, max_bond_gap: int=1000000, dataset_name: Optional[str]=None, metadata: Optional[Mapping[str, Any]]=None, title: Optional[str]=None, regions: Optional[Sequence[str]]=None, background: str='#151918', particle_size: float=1.0, bond_size: float=1.0, export_profile: str='molecular-studio-v1', render_config: RenderConfigInput=None, overwrite: bool=True)`

**Returns:** `Path`


Build a 3D genome scene and write it as one portable interactive HTML.

This convenience API performs no external process invocation. JavaScript,
CSS, typed arrays, colors, and polymer bonds are embedded in the output, so
the result can be downloaded and opened directly or served over Remote SSH.
``render_config`` embeds the same typed request accepted by
``render_structure_png``.

Returns
-------
pathlib.Path
    Path to the atomically written HTML artifact.


**Source:** Line 154 in [visualize/html.py](visualize/html.py#L154)


---


## visualize.render

**File:** `visualize/render.py`

**Public functions:** 1


### `render_structure_png(source: RenderSource, output: Union[str, os.PathLike], *, config: RenderConfigInput=None, backend: Optional[BrowserBackendConfig]=None, overwrite: bool=True)`

**Returns:** `RenderResult`


Render a 3dg structure, scene, viewer, or standalone HTML to PNG.

The function uses the repository-pinned Playwright runtime and a cached
headless Chromium. It does not start a server, open a GUI, use the network,
or download browser assets. ``config`` accepts a ``StructureRenderConfig``,
a strict mapping, JSON text, or a JSON file path. A ``StructureViewer``
supplies its embedded config when ``config`` is omitted.

Parameters
----------
source
    A 3dg path, native structure DataFrame, ``StructureScene``,
    ``StructureViewer``, or compatible standalone viewer HTML.
output
    Destination PNG path. The image is validated before atomic replacement.
config
    Reproducible region, camera, cross-section slab, layer, profile,
    colorbar, style, and raster settings.
backend
    Optional browser executable/module paths, timeout, viewport, and
    software-rendering policy. These execution details are kept separate
    from the scientific render configuration.
overwrite
    Whether an existing output may be replaced.

Returns
-------
RenderResult
    Output dimensions, file size, scene checksum, and backend/render report.

Notes
-----
Install the matching Chromium asset once with
``./run_in_env.sh --network internet playwright install --only-shell chromium``
when the backend reports that its executable is missing.


**Source:** Line 54 in [visualize/render.py](visualize/render.py#L54)


---


## visualize.structure

**File:** `visualize/structure.py`

**Public functions:** 2


### `structure_scene(structure: StructureInput, *, color: ColorInput='chrom', palette: Optional[Mapping[Any, Any]]=None, cmap: str='viridis', vmin: Optional[float]=None, vmax: Optional[float]=None, missing_color: Any='#aeb7ba', polymer_bonds: bool=True, max_bond_gap: int=1000000, dataset_name: Optional[str]=None, metadata: Optional[Mapping[str, Any]]=None)`

**Returns:** `StructureScene`


Build a renderer-independent scene from a 3dg path or native DataFrame.

Parameters
----------
structure
    A five-column ``.3dg`` path or a DataFrame. DataFrames may use a
    ``(chr, pos)`` MultiIndex with ``x, y, z`` columns, or flat
    ``chrom/chr, pos, x, y, z`` columns.
color
    ``"chrom"`` for chromosome territory colors, a DataFrame column name,
    one constant color specification, an ``N x 3`` RGB array, or a
    length-N sequence of color specifications. Numeric columns use
    ``cmap``; categorical columns use ``palette`` or a deterministic
    categorical palette. A string matching a column name takes precedence
    over interpreting that string as a color.
palette
    Optional category-to-color mapping. With ``color="chrom"``, keys are
    chromosome names.
cmap, vmin, vmax, missing_color
    Continuous color normalization and missing-value controls.
polymer_bonds
    Whether to connect adjacent genome-bin particles into polymer chains.
max_bond_gap
    Maximum genomic separation in base pairs for an adjacent-bin bond.
    Bonds never cross chromosome or haplotype boundaries.
dataset_name
    Display name stored in the scene. Defaults to the input file name.
metadata
    Additional JSON-compatible dataset metadata.

Returns
-------
StructureScene
    Validated positions, colors, genomic identity, and polymer bonds. The
    scene does not contain Three.js or other renderer-specific objects.


**Source:** Line 55 in [visualize/structure.py](visualize/structure.py#L55)


---


### `structure_scene_with_feature(structure: StructureInput, feature: FeatureInput, *, feature_name: str='feature', feature_value: Optional[str]=None, duplicate_haplotypes: bool=True, cmap: str='pymol-cpg', vmin: Optional[float]=None, vmax: Optional[float]=None, missing_color: Any='#aeb7ba', polymer_bonds: bool=True, max_bond_gap: int=1000000, dataset_name: Optional[str]=None, metadata: Optional[Mapping[str, Any]]=None)`

**Returns:** `StructureScene`


Build a scene colored by a genome-bin feature track.

Parameters
----------
structure
    A five-column 3dg path or native structure DataFrame.
feature
    A headerless ``chrom, pos, value`` TSV path, a DataFrame, or a
    two-level-indexed Series. DataFrames may use ``pos`` or ``start``.
feature_name
    Scientific label stored with the continuous color mapping and used by
    exported color bars.
feature_value
    Value column for a feature DataFrame. It may be omitted when exactly
    one non-coordinate column exists.
duplicate_haplotypes
    Match a haploid reference row to both ``(mat)`` and ``(pat)`` chains,
    while also normalizing optional ``chr`` prefixes.
cmap, vmin, vmax, missing_color
    Continuous feature color mapping. ``"pymol-cpg"`` reproduces the
    magenta-gray-green mapping used by the legacy PyMOL CpG workflow.
polymer_bonds, max_bond_gap
    Polymer-chain construction controls passed to :func:`structure_scene`.
dataset_name, metadata
    Optional scene identity and JSON-compatible metadata.

Returns
-------
StructureScene
    A renderer-independent scene with per-particle colors. Polymer bonds
    interpolate the colors of their endpoint particles.

Notes
-----
Missing feature rows use ``missing_color``. Duplicate feature keys are
rejected because their mapping to a particle would be ambiguous.


**Source:** Line 212 in [visualize/structure.py](visualize/structure.py#L212)


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


### `gen_adata(qc, cache_dir, rewrite=[], debug=False, **args)`


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


**Source:** Line 407 in [wet/afbb.py](wet/afbb.py#L407)


---


### `add_info(annote, rd)`


Add sample information from JSON files to an existing annotation DataFrame.

Args:
    annote: Existing annotation DataFrame
    rd: Directory containing JSON files with additional sample information

Returns:
    Combined DataFrame with original annotations and new sample information


**Source:** Line 279 in [wet/afbb.py](wet/afbb.py#L279)


---


### `add_mapping(annote, threads=16)`


Add DNA-library mapping and primary-mapping rate columns to annotation.

BAM paths are resolved via ``bubble_flow_touched()["hic_mapped"]``.  Samples
whose BAM file is missing are processed with ``mapping_rate(None)`` and
receive zero rates.

Parameters
----------
annote : pandas.DataFrame
    Annotation DataFrame indexed by sample name; must include a
    ``task_dirp`` column.
threads : int, optional
    Total worker threads (divided among parallel ``pysam.flagstat`` calls,
    default 16).

Returns
-------
pandas.DataFrame
    Original DataFrame with two new columns appended: ``"primary mapped"``
    and ``"mapped"``.


**Source:** Line 442 in [wet/afbb.py](wet/afbb.py#L442)


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


**Source:** Line 186 in [wet/afbb.py](wet/afbb.py#L186)


---


### `add_pairs_num(annote, threads, pairs_types=['pairs_c1', 'pairs_c12', 'pairs_c123', 'dip'])`


count contact number of pairs file in parallel
Input:
    annote: must have all pairs_type cols
Output:
    annote with extra cols, pairs_type + "_num" by default


**Source:** Line 87 in [wet/afbb.py](wet/afbb.py#L87)


---


### `add_umis(annote, umip)`


Adding umi stat to annote.
Accept a single umi count matrix path or an iterable (list/tuple) of umi count matrix paths.
Missing files will be skipped with a warning. If no valid files are found, assign 0 for both
"umis" and "genes" for all samples.

Input:
    umip: umi count matrix path or iterable of paths


**Source:** Line 293 in [wet/afbb.py](wet/afbb.py#L293)


---


### `check_RNA(task_dirp)`


Print RNA-seq QC statistics for a bubble_flow task directory.

Reads STAR alignment summary (``star_mapped/Log.final.out``) and
featureCounts assignment summary (``count_gene/gene_assigned.summary``)
from *task_dirp*, then prints and returns the combined statistics.

Parameters
----------
task_dirp : str
    Root directory of a bubble_flow task.


**Source:** Line 510 in [wet/afbb.py](wet/afbb.py#L510)


---


### `count_pairs(filename: str)`

**Returns:** `int`


Count the number of contact records in a 4DN pairs file.

Lines beginning with ``#`` are treated as header/comment lines and
excluded from the count.

Parameters
----------
filename : str
    Path to a gzip-compressed ``.pairs.gz`` file.

Returns
-------
int or None
    Number of contact records, or ``None`` if *filename* is not a string
    or the file does not exist.


**Source:** Line 58 in [wet/afbb.py](wet/afbb.py#L58)


---


### `get_info(log_d)`


Extract and combine information from all JSON files in a directory.

Processes JSON files in chronological order, merging their contents.
Each JSON file should contain sample information as key-value pairs.

Args:
    log_d: Directory containing JSON files with sample information

Returns:
    DataFrame where columns are sample IDs and rows are attributes


**Source:** Line 256 in [wet/afbb.py](wet/afbb.py#L256)


---


### `listdir_bytime(folder)`


List all files in a directory sorted by modification time (oldest first).

Args:
    folder: Path to directory to scan

Returns:
    List of full file paths sorted by modification time


**Source:** Line 239 in [wet/afbb.py](wet/afbb.py#L239)


---


### `mapping_rate(bam, threads=4)`


Compute the primary-mapped and overall mapped rates for a BAM file.

Parameters
----------
bam : str or None
    Path to the sorted BAM file. If ``None``, returns zeros.
threads : int, optional
    Number of threads for ``pysam.flagstat`` (default 4).

Returns
-------
dict
    Dictionary with keys ``"primary mapped"`` and ``"mapped"``, each
    holding the corresponding fraction of QC-passed reads.


**Source:** Line 107 in [wet/afbb.py](wet/afbb.py#L107)


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


**Source:** Line 151 in [wet/afbb.py](wet/afbb.py#L151)


---


### `pick_useful(annote, features=['basic'])`


Pick useful cols from meta dataframe.
Input:
    annote: meta dataframe
    features: list of features to pick, default words are "DNA" and "DNA_RNA"
Output:
    new df with useful cols.


**Source:** Line 629 in [wet/afbb.py](wet/afbb.py#L629)


---


### `real_file(i)`


Return the path if it points to an existing file, otherwise ``None``.

Parameters
----------
i : str or any
    A file path string to validate. Non-string values are returned as-is.

Returns
-------
str or None
    The original path if the file exists; ``None`` otherwise.


**Source:** Line 131 in [wet/afbb.py](wet/afbb.py#L131)


---


### `task_stat(task_dirp, ref=None, threads=32)`

**Returns:** `pd.DataFrame`


Get all infomation about this bubble_flow task in the form of a dataframe.
The meta dataframe is constructed by adding additional columns to "contacts_info.csv" generated by bubble_flow.
TODO:
    make decent dtypes for every column step by step.


**Source:** Line 529 in [wet/afbb.py](wet/afbb.py#L529)


---


### `task_stat_rna_wgs(task_dirp: str, ref: str)`

**Returns:** `pd.DataFrame`


Build a QC summary DataFrame for a combined RNA + WGS bubble_flow task.

Reads sample metadata from ``<task_dirp>/rd``, filters to a fixed set of
QC columns (mapping rates, read counts, phasing scores, etc.), appends UMI
and gene counts from the specified count matrix, and calculates derived
metrics via :func:`add_extra`.

Parameters
----------
task_dirp : str
    Root directory of the bubble_flow task.
ref : str
    Reference genome tag used to locate the count matrix at
    ``<task_dirp>/count_matrix_<ref>/counts.gene.tsv.gz``.

Returns
-------
pandas.DataFrame
    Per-sample QC DataFrame with columns for mapping rates, read counts,
    phasing metrics, UMIs, genes, and derived ratios.


**Source:** Line 564 in [wet/afbb.py](wet/afbb.py#L564)


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


**Source:** Line 216 in [wet/afbb.py](wet/afbb.py#L216)


---


### `zcount(filename: str, target: str)`

**Returns:** `int`


Count the number of lines matching a target string in a gzip-compressed file.

Parameters
----------
filename : str
    Path to the gzip-compressed file.
target : str
    String pattern to search for (passed to ``zgrep -c``).

Returns
-------
int
    Number of matching lines; 0 if none found or file is absent.


**Source:** Line 35 in [wet/afbb.py](wet/afbb.py#L35)


---


## wet.exp_record

**File:** `wet/exp_record.py`

**Public functions:** 6


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


### `drange(a, b, seats)`


No documentation available.


**Source:** Line 12 in [wet/exp_record.py](wet/exp_record.py#L12)


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


### `merge_meta(*dfs)`


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
    c_threshold: [low, high] for contacts, use None for unbounded (e.g., [None, 2000] means only upper limit)
    u_threshold: [low, high] for umis, use None for unbounded (e.g., [None, 2000] means only upper limit)
Output:
    new df with qc-passed rows


**Source:** Line 3 in [wet/qc.py](wet/qc.py#L3)


---


### `plot_qc(meta, ccol='contacts', ucol='umis', c_threshold=[150000, 1000000], u_threshold=[10000, 600000], **kwargs)`


No documentation available.


**Source:** Line 44 in [wet/qc.py](wet/qc.py#L44)


---


### `plot_qc(meta, ccol='contacts', ucol='umis', c_threshold=[150000, 1000000], u_threshold=[10000, 600000], **kwargs)`


No documentation available.


**Source:** Line 75 in [wet/qc.py](wet/qc.py#L75)


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