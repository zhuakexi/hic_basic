"""
Utility functions for threeD.
Using pure python because the way of pyMOL installation is not guaranteed.
"""
from collections import namedtuple
import os
import gzip
from pathlib import Path
import csv
# --- data ---
def get_ref_dir():
    """
    Return a path to src's ref
    """
    ref_dir = os.path.join(os.path.dirname(__file__), "ref")
    return os.path.join(ref_dir, "")
ref_dir = Path(get_ref_dir())
def fetch_centromeres(genome):
    """
    Get centromere regions.
    Only work with mm10 for now.
    TODO: generalize to other genomes
    Input:
        fp (str): path to the gap file
    Returns:
        dict(chrom) of list[start, end]
    """
    # mouse cytoband file does not have centromeric, gvar, and stalk regions
    # using gap files instead
    files = {
        "mm10" : "mm10_gap.csv.gz"
    }
    if genome == "mm10":
        # mouse cytoband file does not have centromeric, gvar, and stalk regions
        # using gap files instead
        with gzip.open(files[genome],"rt") as f:
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
    centromeres = {row.chrom : [row.start, row.end] for row in res}
    return centromeres