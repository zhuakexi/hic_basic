import pandas as pd
import os
import json
from pathlib import Path
from .io import get_ref_dir, load_json
ref_dir = Path(get_ref_dir())
# ZGA gene module
def Jichang2022_embryo_gene_module():
    real_path = os.path.join(get_ref_dir(), "Jichang2022_embryo_gene_module.csv.gz" )
    genes = pd.read_csv(real_path,index_col=0)
    minor_ZGA = genes["minor_ZGA"].dropna()
    major_ZGA = genes["major_ZGA"].dropna()
    totipotent = genes["totipotent"].dropna()
    pluripotent = genes["pluripotent"].dropna()
    maternal = genes["maternal"].dropna()
    return dict(zip(["minor_ZGA", "major_ZGA", "totipotent", "pluripotent", "maternal"], [i.values for i in [minor_ZGA, major_ZGA, totipotent, pluripotent, maternal]]))
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
        gene_sets[i] = list(sheet.loc[i].dropna().astype("string"))
    return gene_sets
# cell cycle genes
def mouse_cell_cycle():
    """
    Cell cycle gene sets in mouse. Used in seurat or scanpy g1/s, g2/m score.
    Output:
        dict of list.
    """
    CC = {}
    CC["g2m"] = load_json(ref_dir / "mouse_g2m_genes.json")
    CC["g1s"] = load_json(ref_dir / "mouse_s_genes.json")
    return CC
def mouse_GO_cell_cycle():
    """
    Mouse cell cycle genes with term: `G2/M transition of mitotic cell cycle` and `G1/S transition of mitotic cell cycle`.
    Output:
        dict of list
    """
    mat = pd.read_csv(ref_dir / "GO_MM_mitotic_cell_cycle_genes.csv.gz", index_col=0)
    gene_sets_CC = {}
    gene_sets_CC["G2/M transition of mitotic cell cycle"] = mat.loc[mat["3"] == "GO:0000086","gene_symbol"].values
    gene_sets_CC["G1/S transition of mitotic cell cycle"] = mat.loc[mat["3"] == "GO:0000082","gene_symbol"].values
    return gene_sets_CC