import concurrent.futures
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Union, List

sys.path.insert(0, "/share/home/ychi/dev/hires_utils")

import dask.dataframe as dd
import numpy as np
import pandas as pd
from hires_utils.hires_io import parse_3dg, parse_ref
from tqdm import tqdm
def spatial_binnify(left:tuple=(-70,-50,-40), right:tuple=(80, 60, 50), step:int=2, plot:bool=True)->pd.DataFrame:
    """
    Generate spatial grid in primary-coordinate(3D) space.
    Input:
        xyzf: first 3 col must be ["ht","dv","lr"], index must be ["chrom","start"]
        step: bin size
        plot: for plot compatibility, use left bound of each bin
    Output:
        3d-binned feature
        if plot, return df with index ["ht","dv","lr"]
        else, return df with index ["chrom","start"]
    """
    bins = {
        "ht" : pd.IntervalIndex.from_breaks(
            list(range(left[0],right[0],step))
        ),
        "dv" : pd.IntervalIndex.from_breaks(
            list(range(left[1],right[1],step))
        ),
        "lr" : pd.IntervalIndex.from_breaks(
            list(range(left[2],right[2],step))
        )
    }

    if not plot:
        grid = pd.DataFrame(
            index = pd.MultiIndex.from_product(
                bins.values(),
                names = ["ht","dv","lr"]
            )
        )
    else:
        grid = pd.DataFrame(
            index = pd.MultiIndex.from_product(
                [view.left for view in bins.values()],
                names = ["ht","dv","lr"]
            )
        )
    grid.sort_index(inplace=True)
    return grid
def sig_primary_coords(primary_views, sig, _3dg):
    """
    Transform xyz coordinates to primary coordinates.
    Input:
        primary_views: primary_views[sample]
        sig: 3 iterable, ht,dv, lr
        _3dg: xyz of sample
    Output:
        new _3dg; columns as head_tail coordinates, 
            dorsal_ventral coordinates,
            left_right coordinates; range 0-1
    """
    sig = dict(zip(
        ('head-tail', 'dorsal-ventral', 'left-right'),
        sig
    ))
    bases = pd.DataFrame(
        primary_views["bases"],
        columns = primary_views["name_of_vectors"]
    )
    extents = pd.Series(
        primary_views["extent"],
        index = primary_views["name_of_vectors"][:3]
    )
    o3dg = _3dg - bases["center"].values.reshape(-1,3)
    primary_coords = pd.concat(
        [(o3dg * (bases[name]*sig[name]).values.reshape(-1,3)).sum(axis=1) for 
            name in ["head-tail", "dorsal-ventral", "left-right"]],
        axis=1
    )
    primary_coords.columns = ["ht","dv","lr"]
    return primary_coords
# def pc_binagg_feature(xyzf, step, grouping=["ht","dv","lr"], agg=dict(density=("particle","sum")), plot=True):
#     """
#     Binnify feature and aggregate whithin in primary-coordinate(3D) space.
#     Input:
#         xyzf: first 3 col must be ["ht","dv","lr"]
#         step: bin size
#     Output:
#         3d-binned feature
#     """
#     bins = {
#         "ht" : pd.IntervalIndex.from_breaks(
#             list(range(-70,80,step))
#         ),
#         "dv" : pd.IntervalIndex.from_breaks(
#             list(range(-50,60,step))
#         ),
#         "lr" : pd.IntervalIndex.from_breaks(
#             list(range(-40,50,step))
#         )
#     }
#     data = pd.DataFrame({
#         view : pd.cut(
#                 xyzf[view], 
#                 bins=bins[view]
#                 ) 
#         for view in grouping
#     })
#     data = pd.concat(
#         [data, xyzf.iloc[:,3:]],
#         axis = 1
#     )
#     #return data
#     # binnify
#     bf = data.groupby(grouping).agg(**agg)
#     if plot:
#         bf.reset_index(inplace=True)
#         # using left bound(int) rather than Interval
#         for view in grouping:
#             bf[view] = bf[view].cat.rename_categories(
#                 bf[view].cat.categories.left
#             )
#         bf.set_index(grouping,inplace=True)
#     return bf
# def pileup_bf(
#     primary_views, targets, sample_table, 
#     features, pb_features=None, agg={"density":("particle","sum")},
#     sub=None, step=2,grouping=["ht","dv"]):
#     """
#     Pileup binned features.
#     Input:
#         primary_views: to get three bases of each sample; structured data.
#         targets: additional direction for each base. df.
#         sample_table: to get path of structure(3dg) file. df.
#         features: bulk features, df:=bin*feature
#         pb_features: dict of df, df:=bin*sample
#         agg: aggregation method; dict; find key in pb_features or columns of features.
#         sub: epression for primary views, used for slicing.
#     Output:
#         bfs
#     """
#     bfs = []
#     c = 0
#     for sample, row in targets.iterrows():
#         if row.abs().sum() < 3:
#             # abnormal alignment
#             continue
#         # signed primary coords, xyz df
#         spc = sig_primary_coords(
#             primary_views[sample],
#             row.values,
#             parse_3dg(sample_table.loc[sample,"gs"])
#         )
#         if sub is not None:
#             spc = spc.query(sub)
#         if pb_features is not None:
#             # get 
#             ad_pb_features = pd.concat(
#                 [pb_features[i][sample] for i in pb_features],
#                 axis=1,
#                 join="outer"
#             )
#             ad_pb_features.columns = [key for key in pb_features]
#         spc = pd.concat(
#             [spc, features] if pb_features is None else [spc, features, ad_pb_features],
#             axis = 1,
#             join = "inner" # the NA problem for spc, features may be better than ad_pb_features
#         )
#         bf = pc_binagg_feature(
#             spc,
#             step=step,
#             grouping=grouping,
#             agg=agg,
#             plot=True
#         )
#         # adding sample name
#         bf.columns = pd.MultiIndex.from_product(
#             [[sample],bf.columns]
#         )
#         bfs.append(bf)
#         if (c + 1) % 10 == 0:
#             print(c+1,"done")
#         c += 1
#     res = pd.concat(bfs, axis=1,join="outer")
#     return res
# --- pileup version 2 ---
def pc_bin_feature(xyzf:pd.DataFrame, step:int, plot:bool=True)->pd.DataFrame:
    """
    Binnify feature in primary-coordinate(3D) space.
    Input:
        xyzf: first 3 col must be ["ht","dv","lr"], index must be ["chrom","start"]
        step: bin size
    Output:
        3d-binned feature
        if plot, return df with index ["ht","dv","lr"]
        else, return df with index ["chrom","start"]
    """
    bins = {
        "ht" : pd.IntervalIndex.from_breaks(
            list(range(-70,80,step))
        ),
        "dv" : pd.IntervalIndex.from_breaks(
            list(range(-50,60,step))
        ),
        "lr" : pd.IntervalIndex.from_breaks(
            list(range(-40,50,step))
        )
    }
    # cut ht, dv, lr
    data = pd.DataFrame({
        view : pd.cut(
                xyzf[view], 
                bins=bins[view]
                ) 
        for view in ["ht","dv","lr"]
    })
    # append feature columns
    data = pd.concat(
        [data, xyzf.iloc[:,3:]],
        axis = 1
    )
    if data[["ht","dv","lr"]].isna().sum().sum() > 0:
        print("Warning: particles outside of the binning range.")
        data = data.dropna(subset=["ht","dv","lr"], how="any")
    if plot:
        data.reset_index(inplace=True,names=["chrom","start"])
        # using left bound(int) rather than Interval
        for view in ["ht","dv","lr"]:
            data[view] = data[view].cat.rename_categories(
                data[view].cat.categories.left
            )
        data.set_index(["ht","dv","lr"],inplace=True)
    data.sort_index(inplace=True)
    return data

def pc_binagg_feature(xyzf, step, grouping=["ht","dv","lr"],
                      agg=dict(density=("particle","sum")), min_particles={},
                      plot=True, global_min_particles=10):
    """
    Binnify feature and aggregate whithin in primary-coordinate(3D) space.
    Input:
        xyzf: first 3 col must be ["ht","dv","lr"]
        step: bin size
        min_particles: min particles in a bin to be valid, for each feature
        global_min_particles: default min particles, use this if feature not in min_particles
            set to None to disable all min_particles
    Output:
        3d-binned feature
    """
    # binnify
    bin_xyzf = pc_bin_feature(xyzf, step, plot=plot)
    # aggregation
    agg_bin_xyzf = bin_xyzf.groupby(grouping, observed=False).agg(**agg)
    valid_particles = bin_xyzf.groupby(grouping, observed=False).count() # all non NA particles, NA are created by missing feature ref
    physical_particles = bin_xyzf.fillna(0).groupby(grouping, observed=False).count() # all particles in the structure
    if global_min_particles is not None:
        for feature in agg:
            if feature in min_particles:
                feature_min_particle = min_particles[feature]
            else:
                feature_min_particle = global_min_particles
            agg_bin_xyzf.loc[valid_particles[agg[feature][0]] < feature_min_particle, feature] = pd.NA # orig col for aggregation
    return agg_bin_xyzf 
def pileup_bf(
    primary_views, targets, sample_table, 
    features, pb_features=None, agg={"density":("particle","sum")},
    sub=None, step=2, grouping=["ht","dv"], min_particles={}, global_min_particle=10
    ):
    """
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
    """
    target_samples = targets.index
    assert all(target_samples.isin(primary_views))
    assert all(target_samples.isin(sample_table.index))

    bfs = []
    c = 0
    for sample, row in targets.iterrows():
        if row.abs().sum() < 3:
            # abnormal alignment
            continue
        # signed primary coords, xyz df
        if "gs" in sample_table.columns:
            _3dg_col = "gs"
        elif "20k_g_struct1" in sample_table.columns:
            _3dg_col = "20k_g_struct1"
        else:
            raise ValueError("No 3dg column in sample_table")
        spc = sig_primary_coords(
            primary_views[sample],
            row.values,
            parse_3dg(sample_table.loc[sample,_3dg_col])
        )
        if sub is not None:
            spc = spc.query(sub)
        if pb_features is not None:
            # get 
            ad_pb_features = pd.concat(
                [pb_features[i][sample] for i in pb_features],
                axis=1,
                join="outer"
            )
            ad_pb_features.columns = [key for key in pb_features]
        spc = pd.concat(
            [spc, features] if pb_features is None else [spc, features, ad_pb_features],
            axis = 1,
            join = "outer"
        ).loc[spc.index]
        bf = pc_binagg_feature(
            spc,
            step=step,
            grouping=grouping,
            agg=agg,
            plot=True,
            min_particles=min_particles,
            global_min_particles=global_min_particle
        )
        # adding sample name
        bf.columns = pd.MultiIndex.from_product(
            [[sample],bf.columns]
        )
        bfs.append(bf)
        if (c + 1) % 10 == 0:
            print(c+1,"done")
        c += 1
    res = pd.concat(bfs, axis=1,join="outer")
    return res
def pileup_bf_mt_task(
    sample_name: str,
    primary_view: dict,
    directions,
    _3dg_file,
    features,
    cache_file,
    pb_features=None,
    pb_readers=None,
    agg=None,
    sub=None,
    step=2,
    grouping=["ht","dv","lr"],
    min_particles={},
    global_min_particle=10
):
    spc = sig_primary_coords(
        primary_view,
        directions,
        parse_3dg(_3dg_file)
    )
    if sub is not None:
        # get subset of structure
        spc = spc.query(sub)
    pb_readers = {} if pb_readers is None else pb_readers
    if pb_features is not None:
        ad_pb_features = pd.concat(
            {
                feature_name : \
                    pb_readers[feature_name](feature_file_dict[sample_name]) \
                        if feature_name in pb_readers else parse_ref(
                            feature_file_dict[sample_name],
                            value_name = feature_name
                            )[feature_name]
                for feature_name, feature_file_dict in pb_features.items()
            },
            axis = 1,
            join = "outer"
        )
    spc = pd.concat(
        [spc, features] if pb_features is None else [spc, features, ad_pb_features],
        axis = 1,
        join = "outer"
    ).loc[spc.index]
    bf = pc_binagg_feature(
        spc,
        step=step,
        grouping=grouping,
        agg=agg,
        plot=True,
        min_particles=min_particles,
        global_min_particles=global_min_particle
    )
    bf = bf.reset_index().assign(
        sample_name = sample_name
    )
    bf.to_parquet(
        cache_file
    )
    return cache_file
def echo_file(file_path):
    return file_path    
def pileup_bf_mt(
    sample_names: list,
    primary_views: dict,
    targets: pd.DataFrame,
    _3dg_files: Union[pd.DataFrame, List[str]],
    features: pd.DataFrame,
    outfile: Union[str, Path] = None,
    pb_features: Union[dict, pd.DataFrame] = None,
    pb_readers: dict = None,
    agg={"density":("particle","sum")},
    sub=None,
    step=2,
    grouping=["ht","dv","lr"],
    min_particles={},
    global_min_particle=10,
    cache_dir: Union[str, Path] = None,
    n_threads = 8
    ):
    """
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
    """
    # check outfile
    assert outfile is not None
    # check sample
    assert all([sample in primary_views for sample in sample_names])
    assert all([sample in targets.index for sample in sample_names])
    assert all([sample in _3dg_files.index for sample in sample_names]) \
        if isinstance(_3dg_files, pd.DataFrame) \
        else len(_3dg_files) == len(sample_names)
    # check _3dg_files
    if isinstance(_3dg_files, pd.DataFrame):
        assert ("gs" in _3dg_files.columns) or ("20k_g_struct1" in _3dg_files.columns)
    # check cache_dir
    assert cache_dir is not None
    # check pd_readers
    if pb_readers is not None:
        assert (pb_features is not None) and isinstance(pb_features, dict)

    gs_col = "gs" if "gs" in _3dg_files.columns else "20k_g_struct1"
    if isinstance(_3dg_files, pd.DataFrame):
        _3dg_files = _3dg_files[gs_col] # pd.Series
    else:
        _3dg_file = dict(zip(sample_names, _3dg_file))
    if isinstance(pb_features, pd.DataFrame):
        pb_features = pb_features.to_dict()
    if not Path(cache_dir).exists():
        os.makedirs(cache_dir, exist_ok=True)

    cache_files = []
    with ProcessPoolExecutor(n_threads) as executor:
        futures = []
        for i, sample_name in enumerate(sample_names):
            directions = targets.loc[sample_name] # pd.Series
            if directions.abs().sum() < 3:
                # abnormal alignment
                print("Skipping %s because of abnormal alignment" % sample_name)
                continue
            _3dg_file = _3dg_files[sample_name]
            if not Path(_3dg_file).exists():
                print("Skipping %s because of missing 3dg file" % sample_name)
                continue
            cache_file = Path(cache_dir) / ("%s.parquet" % sample_name)
            if cache_file.exists():
                #print("Cache file exists, loading from cache")
                future = executor.submit(
                    echo_file, # echo for aggregation
                    cache_file
                )
            else:
                future = executor.submit(
                    pileup_bf_mt_task,
                    sample_name,
                    primary_views[sample_name],
                    directions, 
                    _3dg_file,
                    features,
                    cache_file,
                    pb_features,
                    pb_readers,
                    agg=agg,
                    sub=sub,
                    step=step,
                    grouping=grouping,
                    min_particles=min_particles,
                    global_min_particle=global_min_particle,
                )
            futures.append(future)
        for future in tqdm(
            concurrent.futures.as_completed(futures),
            total=len(futures),
            desc="Processing cells"
            ):
            cache_files.append(future.result())
    # TODO: use queued multi-thread writing
    print("Combining results...")
    with pd.HDFStore(
        str(outfile),complevel=3,complib="zlib",mode="w"
    ) as store:
        for filep in tqdm(cache_files, desc="Writing to file"):
            df = pd.read_parquet(filep)
            store.put(
                "main", df, format="table", index=False,
                data_columns = ["lr","dv","ht","sample_name"],
                append=True,
                min_itemsize={"sample_name": 32}
            )
        print("Indexing...")
        store.create_table_index("main",columns=["lr","dv","ht"],optlevel=6)
    print("All done.")
    return outfile
def project_back(
    sample_name: str,
    primary_view: dict,
    directions:pd.Series,
    _3dg_file:pd.DataFrame,
    spatial_features:pd.DataFrame,
    cache_file,
    step=2,
):
    """
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
    """
    spc = sig_primary_coords(
        primary_view,
        directions,
        parse_3dg(_3dg_file)
    )
    # (ht, dv, lr) [chrom, start]
    bin_spc = pc_bin_feature(spc, step, plot=True)
    bin_spc = bin_spc.assign(
        sample_name = sample_name
    )

    bin_spc = pd.merge(
        bin_spc.reset_index(),
        spatial_features.reset_index(),
        on = ["ht","dv","lr"],
        how = "left"
    )

    bin_spc = bin_spc.drop(["ht","dv","lr"],axis=1)
    bin_spc = bin_spc.fillna(0)
    bin_spc.to_parquet(
        cache_file
    )
    return cache_file
def project_back_mt_task(
    sample_name: str,
    primary_view: dict,
    directions:pd.Series,
    _3dg_file:pd.DataFrame,
    spatial_features:pd.DataFrame,
    cache_file,
    step=2,
):
    """
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
    """
    try:
        project_back(
            sample_name,
            primary_view,
            directions,
            _3dg_file,
            spatial_features,
            cache_file,
            step=step
        )
        return cache_file
    except Exception as e:
        print("Error in %s: %s" % (sample_name, str(e)))
        return None
def project_back_mt(
    sample_names: list,
    primary_views: dict,
    targets: pd.DataFrame,
    _3dg_files: Union[pd.DataFrame, List[str]],
    spatial_features: pd.DataFrame,
    outfile: Union[str, Path] = None,
    step=2,
    cache_dir: Union[str, Path] = None,
    n_threads = 8
    ):
    """
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
    """
    # check outfile
    assert outfile is not None
    # check sample
    assert all([sample in primary_views for sample in sample_names])
    assert all([sample in targets.index for sample in sample_names])
    assert all([sample in _3dg_files.index for sample in sample_names]) \
        if isinstance(_3dg_files, pd.DataFrame) \
        else len(_3dg_files) == len(sample_names)
    # check _3dg_files
    if isinstance(_3dg_files, pd.DataFrame):
        assert ("gs" in _3dg_files.columns) or ("20k_g_struct1" in _3dg_files.columns)
    # check cache_dir
    assert cache_dir is not None

    gs_col = "gs" if "gs" in _3dg_files.columns else "20k_g_struct1"
    if isinstance(_3dg_files, pd.DataFrame):
        _3dg_files = _3dg_files[gs_col] # pd.Series
    else:
        _3dg_file = dict(zip(sample_names, _3dg_file))
    if not Path(cache_dir).exists():
        os.makedirs(cache_dir, exist_ok=True)

    cache_files = []
    with ProcessPoolExecutor(n_threads) as executor:
        futures = []
        for i, sample_name in enumerate(sample_names):
            directions = targets.loc[sample_name] # pd.Series
            if directions.abs().sum() < 3:
                # abnormal alignment
                print("Skipping %s because of abnormal alignment" % sample_name)
                continue
            _3dg_file = _3dg_files[sample_name]
            if not Path(_3dg_file).exists():
                print("Skipping %s because of missing 3dg file" % sample_name)
                continue
            cache_file = Path(cache_dir) / ("%s.parquet" % sample_name)
            if cache_file.exists():
                #print("Cache file exists, loading from cache")
                future = executor.submit(
                    echo_file, # echo for aggregation
                    cache_file
                )
            else:
                future = executor.submit(
                    project_back_mt_task,
                    sample_name,
                    primary_views[sample_name],
                    directions, 
                    _3dg_file,
                    spatial_features,
                    cache_file,
                    step=step,
                )
            futures.append(future)
        for future in tqdm(
            concurrent.futures.as_completed(futures),
            total=len(futures),
            desc="Processing cells"
            ):
            try:
                result = future.result()
                cache_files.append(result)
            except Exception as e:
                print("Error: %s" % str(e))
                continue
    # TODO: use queued multi-thread writing
    print("Combining results...")
    with pd.HDFStore(
        str(outfile),complevel=3,complib="zlib",mode="w"
    ) as store:
        for filep in tqdm(cache_files, desc="Writing to file"):
            df = pd.read_parquet(filep)
            store.put(
                "main", df, format="table", index=False,
                data_columns = ["chrom","start"],
                append=True,
                min_itemsize={"sample_name": 32, "chrom": 32}
            )
        print("Indexing...")
        store.create_table_index("main",columns=["chrom","start"],optlevel=6)
    print("All done.")
    return outfile
def consecutive_slice_pileup_bf(
    primary_views, targets, sample_table, 
    features, step=2, upper=6, lower=-6, axis="ht", **args
    ):
    slices = pd.IntervalIndex.from_tuples(
        list(zip(
            range(lower,upper,step),
            range(lower+step,upper+step,step)
            )),
        closed="left"
    )
    subs = {
        int(s.left) : '%s>=%d and %s<%d' % (axis,s.left,axis,s.right)
        for s in slices
        }
    grouping = [i for i in ["ht","dv","lr"] if i != axis]
    pileup_bfs = {
        left : pileup_bf(
            primary_views, targets, sample_table, 
            features,
            step=step,
            grouping=grouping,
            sub=sub,
            **args
        ) for left, sub in subs.items()
    }
    pileup_bfs = pd.concat(pileup_bfs, axis=0, join="outer")
    pileup_bfs.index.names = ["clip:"+axis, *grouping]
    return pileup_bfs
def pileup_bf_parallel_task(primary_views, targets, sample_table, features, step, grouping, sub, axis, prefix, key, *args):
    # This function would call the original pileup_bf function with the provided parameters.
    # Assuming pileup_bf is defined elsewhere and can accept these parameters.
    args = dict(args)
    if prefix is not None:
        assert key is not None
        cache_file = Path(prefix) / ("__%d__.parquet" % key)
        if os.path.exists(cache_file):
            print("Cache file exists, loading from cache")
            return cache_file
    res = pileup_bf(
        primary_views, targets, sample_table,
        features, step=step, grouping=grouping,
        sub=sub, **args
        )
    if prefix is not None:
        orig_index = res.index.names
        res = res.assign(
            **{axis: key}
        ).set_index(
            axis, append=True
            ).reorder_levels(
                [axis, *orig_index]
                )
        res = res.swaplevel(0,1,axis=1)
        res.to_parquet(cache_file)
        return cache_file
    else:
        return res

def consecutive_slice_pileup_bf_parallel(primary_views, targets, sample_table, features, step=2, upper=6, lower=-6, axis="ht", prefix=None, threads=8, **args):
    slices = pd.IntervalIndex.from_tuples([
        (left, left+step) for left in range(lower, upper, step)
    ], closed="left")
    
    subs = {
        int(s.left): '%s>=%d and %s<%d' % (axis, s.left, axis, s.right)
        for s in slices
    }
    
    grouping = [i for i in ["ht", "dv", "lr"] if i != axis]
    
    tasks = [(primary_views, targets, sample_table, features, step, grouping, sub, axis, prefix, key) + tuple(args.items()) for key, sub in subs.items()]
    key_of_tasks = [left for left, _ in subs.items()]
    
    if prefix is not None:
        os.makedirs(prefix, exist_ok=True)

    pileup_bfs = []
    with ProcessPoolExecutor(threads) as executor:
        futures = [executor.submit(pileup_bf_parallel_task, *task) for task in tasks]
        for future in futures:
            pileup_bfs.append(future.result())
    pileup_bfs = dict(zip(key_of_tasks, pileup_bfs))
    if prefix is not None:
        # just return the file paths
        # pileup_bfs = {
        #     key: pd.read_pickle(pileup_bfs[key]) for key in pileup_bfs
        # }
        return pileup_bfs
    pileup_bfs = pd.concat(pileup_bfs, axis=0, join="outer")
    pileup_bfs.index.names = ["clip:"+axis, *grouping]
    #pileup_bfs.index.names = [axis, *grouping]
    
    return pileup_bfs

# --- independent voxelize ---
def voxelize(_3dg, primary_view, target, step=2, plot=True, observed=True, lr90=False):
    xyz = sig_primary_coords(
        primary_view,
        target,
        _3dg
    )
    if lr90:
        flip_x = np.array([
            [-1, 0, 0],
            [0, 1, 0],
            [0, 0, 1]
        ])
        xyz = xyz.dot(flip_x.T)
        theta = np.radians(90)
        R_z = np.array([
            [np.cos(theta), np.sin(theta), 0],
            [-np.sin(theta), np.cos(theta), 0],
            [0, 0, 1]
        ])
        xyz = xyz.dot(R_z.T)
        #print(xyz)
        #print(type(xyz))
        xyz.columns = ["ht","dv","lr"]
    xyzd = xyz.assign(
        density = 1
    )
    bins = {
        "ht" : pd.IntervalIndex.from_breaks(
            list(range(-70,80,step))
        ),
        "dv" : pd.IntervalIndex.from_breaks(
            list(range(-50,60,step))
        ),
        "lr" : pd.IntervalIndex.from_breaks(
            list(range(-40,50,step))
        )
    }
    data = pd.DataFrame({
        view : pd.cut(
                xyzd[view], 
                bins=bins[view]
                ) 
        for view in ["ht","dv","lr"]
    })
    data = pd.concat(
        [data, xyzd.iloc[:,3]],
        axis = 1
    )
    if plot:
        data.reset_index(inplace=True,drop=True)
        # using left bound(int) rather than Interval
        for view in ["ht","dv","lr"]:
            data[view] = data[view].cat.rename_categories(
                data[view].cat.categories.left
            )
        data.set_index(["ht","dv","lr"],inplace=True)
    data.sort_index(inplace=True)
    voxel = data.groupby(
        level=[0,1,2],
        observed=observed
        )["density"].sum()
    return voxel
def mix_layer(x_voxel, axis="ht", near=1):
    """
    Mix the layer of the x_voxel with the near layer
    """
    mixed = x_voxel.copy()
    for i in range(1, near+1):
        mixed = mixed + x_voxel.shift({axis: i}, fill_value=0)
        mixed = mixed + x_voxel.shift({axis: -i}, fill_value=0)
    return mixed