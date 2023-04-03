import pandas as pd
import os
import json
from pathlib import Path
from .hicio import get_ref_dir, load_json
from .utils import two_sets
import gzip
from collections import namedtuple
import csv
ref_dir = Path(get_ref_dir())
# basic genome features
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
    gap_files = {
        "mm10" : ref_dir / "mm10_gap.csv.gz"
    }
    len_files = {
        "mm10" : ref_dir / "mm10.len.tsv"
    }
    if genome == "mm10":
        # mouse cytoband file does not have centromeric, gvar, and stalk regions
        # using gap files instead
        with gzip.open(gap_files[genome],"rt") as f:
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
    cent_chromlen = {row.chrom : [row.start, row.end] for row in res}
    if genome == "mm10":
        with open(len_files[genome],"r") as f:
            fcsv = csv.reader(f, delimiter="\t")
            #header
            #header = next(fcsv)
            header = ["chromosome", "length"]
            # dtypes
            dtypes = [str, int]
            Len = namedtuple("Len", header)
            res = []
            for row in fcsv:
                row = Len(*(convert(value) for convert, value in zip(dtypes, row)))
                res.append(row)
    for row in res:
        if row.chromosome in cent_chromlen:
            cent_chromlen[row.chromosome].append(row.length)
    #print(cent_chromlen)
    return cent_chromlen
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
def chromosomes(genome, order=False):
    """
    Get chromosome lengths.
    Returns:
        pandas.DataFrame
    """
    files = {
        "hg19" : ref_dir / "hg19.len.tsv",
        "hg19_dip" : ref_dir / "hg19.dip.len.tsv",
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
            print("ref: neither valid abbrevations nor valid reference file")
    if isinstance(order, bool):
        if order: # using default order
            data.index = data.index.as_ordered()
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