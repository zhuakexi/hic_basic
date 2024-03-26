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