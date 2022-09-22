"""
DESCRIPTION:
Select particles in telomere flanking regions.

USAGE:
sele_telo obj genome flank

PARAMS:
obj (string or PyMOL selection)
    name of the PyMOL object/selection in which
    to find the telomere particles
genome (string)
    label of the genome used to find the telomere regions
flank (int)
    number of flanking particles to select, must be positive

RETURNS:
    a newly created selection with name `obj`_telo.
"""

from bisect import bisect_right, bisect_left
from collections import namedtuple
import gzip
import csv
from turtle import right
from pymol import cmd
def fetch_telomeres(genome):
    """
    Get centromere regions.
    Returns:
        dict(chrom) of list[left centromere, right centromere] of list[start, end]
    """
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
                if row.type == "telomere":
                    res.append(Cent(row.chrom, row.chromStart, row.chromEnd))
    telo = {}
    for row in res:
        if row.chrom not in telo:
            telo[row.chrom] = []
        telo[row.chrom].append((row.start, row.end))
    telo = {key : sorted(telo[key], key=lambda x : x[0]) for key in telo}
    return telo
def sele_telo(obj, genome, flank=1):
    flank = int(flank)
    SelName = obj + "_telo"

    # get all atoms (chrom, name)
    # cmd.iterate(selection str, expression str, space=)
    # expression string can use:
    #   1. default exposed varaibles
    #   2. stored.xxx
    #   3. custom variable space (space = myspace)
    myspace = {"bins":[]}
    cmd.iterate('(all)', "bins.append((chain, name))", space=myspace)
    
    # find atoms (chrom, name) in centromeres
    # tidy atoms
    atoms = {}
    for atom in myspace["bins"]:
        if atom[0] not in atoms:
            atoms[atom[0]] = []
        else:
            atoms[atom[0]].append(int(atom[1]))
    atoms = {chrom : sorted(atoms[chrom]) for chrom in atoms}
    # tidy centromeres reference
    telomeres = fetch_telomeres(genome)
    # find atoms in centromeres and select
    cmd.select(SelName, "None")
    for chrom in atoms:
        #print(chrom)
        #print(atoms.keys())
        #print(centromeres.keys())
        if chrom not in telomeres:
            continue
        left_end = bisect_right(atoms[chrom], telomeres[chrom][0][1])
        for atom in atoms[chrom][left_end:left_end+flank]: # select left telomere right flank
            cmd.select(SelName, SelName + " or (chain " + chrom + " and name " + str(atom) + ")")
        #print(telomeres[chrom])
        right_start = bisect_left(atoms[chrom], telomeres[chrom][1][0])
        for atom in atoms[chrom][right_start-flank:right_start]: # select right telomere left flank
            cmd.select(SelName, SelName + " or (chain " + chrom + " and name " + str(atom) + ")")
    return SelName
cmd.extend("sele_telo", sele_telo)