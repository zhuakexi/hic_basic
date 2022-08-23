import pandas as pd
import os
from .io import get_ref_dir
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