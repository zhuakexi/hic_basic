import csv
import gzip
import json
import os
from collections import namedtuple
from pathlib import Path

import numpy as np
import pandas as pd
from hires_utils.mmcif import chrom_rm_suffix, chrom_rm_prefix

from .hicio import get_ref_dir, load_json
from .utils import two_sets
ref_dir = Path(get_ref_dir())


### --- download data --- ###
import time

def download_geo(geo_id, outdir, method="ftp", timeout=30):
    """
    Download GEO data with enhanced messages, diagnostics, and timeout.
    
    Input:
        geo_id: GEO ID; string
        outdir: output directory; string
        method: download method; string; "ftp" or "https"
        timeout: timeout for downloads in seconds; int
    """
    print(f"Starting download for GEO ID: {geo_id} using method: {method}")
    os.makedirs(outdir, exist_ok=True)

    if method == "ftp":
        import ftplib
        try:
            ftp = ftplib.FTP("ftp.ncbi.nlm.nih.gov", timeout=timeout)
            ftp.login()
            ftp.cwd(f"geo/series/{geo_id[:-3]}nnn/{geo_id}/suppl/")
            files = ftp.nlst()
            print(f"Found {len(files)} files to download.")
            for file in files:
                print(f"Downloading {file}...")
                start_time = time.time()
                with open(os.path.join(outdir, file), "wb") as f:
                    ftp.retrbinary(f"RETR {file}", f.write)
                elapsed_time = time.time() - start_time
                print(f"Downloaded {file} in {elapsed_time:.2f} seconds.")
            ftp.quit()
        except ftplib.all_errors as e:
            print(f"FTP error: {e}")
    elif method == "https":
        import requests
        import re
        try:
            url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={geo_id}"
            print(f"Fetching download links from {url}...")
            r = requests.get(url, timeout=timeout)
            r.raise_for_status()
            links = re.findall(r'href=[\'"]?([^\'" >]+)', r.text)
            gz_links = [link for link in links if link.endswith(".gz")]
            print(f"Found {len(gz_links)} files to download.")
            for link in gz_links:
                print(f"Downloading {link}...")
                start_time = time.time()
                r = requests.get(link, timeout=timeout)
                r.raise_for_status()
                with open(os.path.join(outdir, os.path.basename(link)), "wb") as f:
                    f.write(r.content)
                elapsed_time = time.time() - start_time
                print(f"Downloaded {os.path.basename(link)} in {elapsed_time:.2f} seconds.")
        except requests.RequestException as e:
            print(f"HTTP error: {e}")
    else:
        raise ValueError("Method must be 'ftp' or 'https'")

    print(f"Download completed for GEO ID: {geo_id}. Files saved to {outdir}.")

### --- data manipulates --- ###


def dupref_annote(bins, ref):
    """
    Annotate diploid genome with haploid ref
    Input:
        bins: pd.DataFrame, 2-level index
        ref: pd.DataFrame, 2-level index
    Output:
        pd.DataFrame, 3dg with ref
    """
    bins.index.names = ["chrom","start"]
    ref.index.names = ["chrom","start"]

    bins = bins.reset_index()
    ref = ref.reset_index()
    bins = bins.assign(new_chrom=chrom_rm_suffix(bins["chrom"]))
    # for chr1<->1 inconsistency
    bins["new_chrom"] = chrom_rm_prefix(bins["new_chrom"])
    # strip chr, assume ref is haploid so don't need to strip suffix
    ref = ref.assign(new_chrom=chrom_rm_prefix(ref["chrom"])) 
    ref = ref.drop(columns=["chrom"])
    bins = pd.merge(
        bins,
        ref,
        on = ["new_chrom","start"],
        how = "left"
        ).drop(columns=["new_chrom"])
    bins = bins.set_index(["chrom","start"])
    return bins    
# --- basic genome features ---
def chromosomes(genome, order=False):
    """
    Get chromosome lengths.
    Input:
        genome: genome name or chromosome length table; string or pd.DataFrame;
            Now support hg19, hg19_dip, GRCh38, GRCh38_dip, mm10, mm10_dip
        order: if True, order chromosomes in a natural way (chr1, chr2 ...),
            if list order chromosomes by the list given;
    Return:
        pandas.DataFrame
    """
    files = {
        "hg19" : ref_dir / "hg19.len.tsv",
        "hg19_dip" : ref_dir / "hg19.dip.len.tsv",
        "GRCh38" : ref_dir / "GRCh38.len.tsv",
        "GRCh38_dip" : ref_dir / "GRCh38.dip.len.tsv",
        "mm10" : ref_dir / "mm10.len.tsv",
        "mm10_dip" : ref_dir / "mm10.dip.len.tsv",
    }
    if genome in files:
        data = pd.read_table(
            files[genome],
            index_col=0,
            names=["chrom", "length"]
            )
    else:
        try:
            data = pd.read_table(
                genome,
                index_col=0,
                names=["chrom", "length"]
                ) # input is a file path
        except FileNotFoundError:
            raise ValueError("ref: neither valid abbrevations nor valid reference file")
    if isinstance(order, bool):
        if order: # using default order
            data.index = pd.CategoricalIndex(data.index, categories=data.index, ordered=True)
    else: # provided order
        data.index = pd.CategoricalIndex(data.index, categories=order, ordered=True)
    return data
def fetch_centromeres(genome):
    """
    Get centromere regions.
    Returns:
        pandas.DataFrame
    """
    files = {
        "hg19_dip" : ref_dir / "hg19_centromeres.csv.gz",
        "mm10" : ref_dir / "mm10_gap.csv.gz"
    }
    if genome == "hg19_dip":
        # tediously adding the chr prefix
        # just read pre formatted file
        return pd.read_csv(files[genome])
    elif genome == "mm10":
        # mouse cytoband file does not have centromeric, gvar, and stalk regions
        # using gap files instead
        gap = pd.read_csv(files[genome])
        gap = gap.loc[gap["type"] == "centromere",["chrom","chromStart","chromEnd"]]
        gap = gap.reset_index(drop=True)
        gap.columns = ["chrom","start","end"]
        return gap
def fetch_cent_chromlen(genome):
    """
    Get centromere regions and chromosome lengths.
    Only work with mm10 for now.
    TODO: generalize to other genomes
    Input:
        fp (str): path to the gap file
    Returns:
        dict(chrom) of list[start, end, chrom_length]
    """
    # mouse cytoband file does not have centromeric, gvar, and stalk regions
    # using gap files instead
    cytoband_files = {
        "hg19" : ref_dir / "hg19_cytoBand.txt.gz",
        "GRCh38" : ref_dir / "GRCh38.cytoBandIdeo.txt.gz"# ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/CytoBandIdeo.txt.gz
    }
    gap_files = {
        "mm10" : ref_dir / "mm10_gap.csv.gz"
    }
    len_files = {
        "mm10" : ref_dir / "mm10.len.tsv"
    }
    DIP=False
    if genome.endswith("dip"):
        DIP=True
        base_genome, _ = genome.split("_")
    else:
        base_genome = genome
    # fetch base genome centromeres
    if base_genome == "mm10":
        # mouse cytoband file does not have centromeric, gvar, and stalk regions
        # using gap files instead
        with gzip.open(gap_files[base_genome],"rt") as f:
            fcsv = csv.reader(f)
            #header
            header = next(fcsv)
            # orig gap file is '#"bin"' here
            header[0] = "bin"
            # dtypes
            dtypes = [int, str, int, int, int, str, int, str, str]
            Bed = namedtuple("Bed", header)
            Cent = namedtuple("SRegion", ["chrom","start","end"])
            res = []
            for row in fcsv:
                row = Bed(*(convert(value) for convert, value in zip(dtypes, row)))
                if row.type == "centromere":
                    res.append(Cent(row.chrom, row.chromStart, row.chromEnd))
            res.append(Cent("chrY", 4050275, 4181802)) #see https://www.ncbi.nlm.nih.gov/grc/mouse/issues/MG-4414
        cent_chromlen = {row.chrom : [row.start, row.end] for row in res}
        # transform to dataframe
        cent_chromlen = pd.DataFrame(cent_chromlen).T
        cent_chromlen.columns = ["start","end"]
        cent_chromlen.index.name = "chrom"
    else:
        cytoband = pd.read_table(
            cytoband_files[base_genome],
            names=["chrom", "start", "end", "name", "gieStain"]
        )
        cent_chromlen = cytoband.query(
            'gieStain == "acen"'
            ).groupby(
                'chrom'
                ).aggregate(
                    {
                        "start": "min",
                        "end": "max"
                    }
                )
        #cent_chromlen = cent_chromlen.set_index("chrom")
    # additional treatment for diploid genomes
    if DIP:
        lengths = chromosomes(genome)
        lengths = lengths.assign(
            new_chrom = chrom_rm_suffix(lengths.index)
        )
        cent_chromlen = pd.merge(
            cent_chromlen,
            lengths,
            left_index=True,
            right_on="new_chrom",
            how="right"
        )
        cent_chromlen = cent_chromlen[["start", "end"]]
    else:
        pass
    # add length
    lengths = chromosomes(genome, order=True)
    cent_chromlen.index = pd.Categorical(
        cent_chromlen.index,
        categories=lengths.index, # length index is already ordered
        ordered=True
    )
    cent_chromlen = pd.concat(
        [cent_chromlen, lengths["length"].rename("chrom_length")],
        axis=1,
        join="inner"
    )
    cent_chromlen = cent_chromlen.sort_index()
    #print(cent_chromlen)
    return cent_chromlen
def id2name(genelist, genome):
    """
    Convert gene IDs to gene names.
    Input:
        genelist: list of gene IDs
        genome: genome name or gene ID to gene name mapping dataframe; string or pd.DataFrame
            gene ID as index when using a dataframe
    Output:
        gene names; list
    """
    id2name_files = {
        "mm10" : ref_dir / "mm10_id_name.csv.gz",
        "hg19" : ref_dir / "hg19_id_name.csv.gz"
    }
    if isinstance(genome, str):
        # using shipped id2name table if genome name is given
        ttable = pd.read_csv(id2name_files[genome], index_col=0)
    else:
        # using provided id2name table
        ttable = genome
    return ttable.loc[genelist, "gene_name"].tolist()
def name2id(genelist, genome, pickid=False):
    """
    Convert gene names to gene IDs.
    Input:
        genelist: list of gene names
        genome: genome name or gene ID to gene name mapping dataframe (name 2 id df will have non-unique index); string or pd.DataFrame
            gene id as index when using a dataframe
        pickid: if True, pick 1 gene ID for genes with multiple IDs, else return NA for genes with multiple IDs
    Output:
        gene IDs; list
    """
    id2name_files = {
        "mm10" : ref_dir / "mm10_id_name.csv.gz",
        "hg19" : ref_dir / "hg19_id_name.csv.gz"
    }
    if isinstance(genome, str):
        # using shipped name2id table if genome name is given
        ttable = pd.read_csv(id2name_files[genome], index_col=0)
    else:
        # using provided name2id table
        ttable = genome
    dup_genes = ttable["gene_name"].value_counts().loc[ttable["gene_name"].value_counts() > 1].index
    if set(genelist).intersection(dup_genes.index):
        # print warning if there are genes with multiple IDs
        print("Warning: there are genes with multiple IDs: %s" % ", ".join(set(genelist).intersection(dup_genes.index)))
    ttable = ttable.reset_index().drop_duplicates(subset="gene_name", keep="first").set_index("gene_name")
    if pickid:
        # pick 1 gene ID for genes with multiple IDs
        return ttable.loc[genelist, "gene_id"].tolist()
    else:
        # return NA for genes with multiple IDs
        ids = ttable.loc[genelist, "gene_id"].tolist()
        for i, g in enumerate(genelist):
            if g in dup_genes:
                ids[i] = np.nan
        return ids
def fetch_TSS(gnames, TSS, name_col="gene_name"):
    """
    Get the all TSS starts sites.
    Input:
        gnames: genename list
        TSS: genome name or TSS reference table; string or pd.DataFrame
        name_col: column name of gene name in TSS table; string
    Output:
        chromname, TSS id, TSS sites; pd.Series
    """
    TSS_files = {
        "mm10" : ref_dir / "mm10_TSS.csv.gz"
    }
    if isinstance(TSS, str):
        # using shipped TSS table if genome name is given
        TSS = pd.read_csv(
            TSS_files[TSS], 
            index_col=name_col
            )
    missing, _ = two_sets(gnames, TSS.index)
    if missing:
        print("Warning %s not in ref TSS table." % " ".join(list(missing)) )
        gnames = [gname for gname in gnames if gname in TSS.index]
    subdf = TSS.loc[gnames]
    # expand txStart string to long form dataframe
    # u_txStarts = subdf["txStart"]
    # u_txStarts = (i.split() for i in u_txStarts)
    # res = ([str(gname), str(chrom), i,int(txStart)] for gname, chrom, u_txStart in zip(subdf.index,subdf["seqname"],u_txStarts) 
    #        for i, txStart in enumerate(u_txStart))
    # res = pd.DataFrame(
    #     res,
    #     columns = ["gname","chrom","txStart_id","txStart"]
    # )
    return subdf.reset_index()


### --- field-specific data --- ###


# ZGA gene module
def Jichang2022_embryo_gene_module():
    real_path = os.path.join(get_ref_dir(), "Jichang2022_embryo_gene_module.csv.gz" )
    genes = pd.read_csv(real_path,index_col=0)
    minor_ZGA = genes["minor_ZGA"].dropna()
    major_ZGA = genes["major_ZGA"].dropna()
    totipotent = genes["totipotent"].dropna()
    pluripotent = genes["pluripotent"].dropna()
    maternal = genes["maternal"].dropna()
    return dict(zip(["minor_ZGA", "major_ZGA", "totipotent", "pluripotent", "maternal"], [set(i.values) for i in [minor_ZGA, major_ZGA, totipotent, pluripotent, maternal]]))
def XijinGe2017_embryo_gene_module():
    """
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
    """
    real_path = os.path.join(get_ref_dir(), "XijinGe2017.xlsx")
    sheet = pd.read_excel(real_path, sheet_name="S2 Gene clusters",dtype="string")
    sheet = sheet.set_index("Annotation")
    sheet = sheet.drop(["Clusters", "Genes"],axis=1)
    gene_sets = {}
    for i in sheet.index:
        gene_sets[i] = set(sheet.loc[i].dropna().astype("string"))
    return gene_sets
# cell cycle genes
def mouse_cell_cycle():
    """
    Cell cycle gene sets in mouse. Used in seurat or scanpy g1/s, g2/m score.
    Output:
        dict of list.
    """
    CC = {}
    CC["g2m"] = set(load_json(ref_dir / "mouse_g2m_genes.json"))
    CC["g1s"] = set(load_json(ref_dir / "mouse_s_genes.json"))
    return CC
def mouse_GO_cell_cycle():
    """
    Mouse cell cycle genes with term: `G2/M transition of mitotic cell cycle` and `G1/S transition of mitotic cell cycle`.
    Output:
        dict of list
    """
    mat = pd.read_csv(ref_dir / "GO_MM_mitotic_cell_cycle_genes.csv.gz", index_col=0)
    gene_sets_CC = {}
    gene_sets_CC["G2/M transition of mitotic cell cycle"] = set(mat.loc[mat["3"] == "GO:0000086","gene_symbol"].values)
    gene_sets_CC["G1/S transition of mitotic cell cycle"] = set(mat.loc[mat["3"] == "GO:0000082","gene_symbol"].values)
    return gene_sets_CC