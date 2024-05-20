import pandas as pd
import gzip
import os
from concurrent import futures

# All kinds of simple calculations on adata
def add_contact_describe(adata, inplace=True):
    """
    Get cell's contact decay statistics, defined in Nagano2017.
    Input:
        adata: AnnData object, must have uns["cdps"]
    Output:
        adata will extra obs cols
    """
    bounds = adata.uns["cdps"].columns.astype(int)
    short = adata.uns["cdps"].loc[:, (23_000 < bounds) & (bounds < 2_000_000)].sum(axis=1)
    mitotic = adata.uns["cdps"].loc[:, (2_000_000 < bounds) & (bounds < 12_000_000)].sum(axis=1)
    data = adata.obs.assign(short_r = short)
    data = data.assign(mitotic_r = mitotic)
    if inplace:
        adata.obs = data
        return None
    else:
        return data
def add_pmUMIs(adata, inplace=True):
    """
    Add cell's g1 UMI ratio.
    Input:
        adata: AnnData object, must have ["g1_UMIs","g2_UMIs"]
    Output:
        adata will extra obs cols
    """
    data = adata.obs.assign(pmUMIs = adata.obs["g1_UMIs"]/(adata.obs["g1_UMIs"] + adata.obs["g2_UMIs"]))
    if inplace:
        adata.obs = data
        return None
    else:
        return data
def chrom_cov(pairs: pd.DataFrame)-> pd.DataFrame:
    """
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
    """
    chr_pairs_count = pairs.groupby(
        ["chr1", "chr2"]
        ).size().rename(
            "count"
            ).reset_index()
    intra = chr_pairs_count.query('chr1 == chr2')[["chr1","count"]]
    intra.columns = ["chr", "intra"]
    inter = chr_pairs_count.query('chr1 != chr2')
    # add chrom1 inter for each chrom
    res = pd.merge(
        intra,
        inter.groupby("chr1")["count"].sum().rename("chrom1_inter"),
        left_on = ["chr"],
        right_index = True,
        how="left"
        )
    # add chrom2 inter for each chrom
    res = pd.merge(
        res,
        inter.groupby("chr2")["count"].sum().rename("chrom2_inter"),
        left_on = ["chr"],
        right_index = True,
        how="left"
        )
    res = res.fillna(0)
    res = res.set_index("chr")
    res = res.assign(
        tot = res["intra"] + res["chrom1_inter"] + res["chrom2_inter"],
    )
    return res[["intra","tot"]]
def pairs_coverage(pairs_fp, genome:str=None, binsize:int=1000000, flavor:str="hickit", sub:str=None)->pd.DataFrame:
    """
    Calculate coverage of pairs in bins.
    Input:
        pairs_fp: str, path to pairs file
        genome: str, genome name
        binsize: int, size of bins
        flavor: str, flavor of bins
        sub: only calculate coverage of a subset of genome
    Output:
        coverage: coverage of pairs in each bin
    """
    assert genome is not None, "Please specify genome"
    pairs = parse_pairs(pairs_fp)
    if sub is not None:
        pairs = pairs.query(sub)
    bins = GenomeIdeograph(genome).bins(binsize, bed=False, flavor=flavor)
    legs = pd.concat(
        [
            pairs[["chr1","pos1"]].rename(columns={"chr1":"chrom","pos1":"pos"}),
            pairs[["chr2","pos2"]].rename(columns={"chr2":"chrom","pos2":"pos"})
        ],
        axis=0, ignore_index=True
    )
    chunks = []
    for chrom, chunk in legs.groupby("chrom"):
        binned = pd.cut(
            chunk["pos"], bins = bins[chrom]
            ).value_counts()
        tidy_binned = pd.concat(
            [
                bins[chrom].to_series().rename("bin"),
                binned.rename("coverage")
                ],
            axis=1,
            join = "outer"
        ).loc[bins[chrom]].drop("bin", axis=1).fillna(0)
        tidy_binned.index = pd.MultiIndex.from_product(
            [[chrom], bins[chrom].left],
            names=["chrom","bin"]
        )
        chunks.append(tidy_binned)
    coverage = pd.concat(chunks, axis=0)
    return coverage
# cov = pairs_coverage(
#     "/shareb/ychi/repo/sperm_struct/formal_pipeline/pairs_c12/BJ8007.c12.pairs.gz",
#     genome = "mm10",
#     binsize = 20000,
#     )
def pairs_coverages(pairs_fps, sample_names, n_jobs=8, genome:str=None, binsize:int=1000000, flavor:str="hickit", sub=None)->pd.DataFrame:
    """
    Multiple pairs_coverage
    Input:
        pairs_fps: list of str, paths to pairs files
        sample_names: list of str, sample names
        genome, binsize, flavor: see pairs_coverage
    Output:
        n_bins * n_samples DataFrame
    """
    with ProcessPoolExecutor(n_jobs) as executor:
        futures = {
            executor.submit(pairs_coverage, pairs_fp, genome, binsize, flavor, sub): sample_name
            for pairs_fp, sample_name in zip(pairs_fps, sample_names)
        }
        results = []
        for future in tqdm(as_completed(futures), total=len(futures)):
            results.append(future.result())
        covs = pd.concat(dict(zip(sample_names, results)), axis=1,join="outer")
    return covs