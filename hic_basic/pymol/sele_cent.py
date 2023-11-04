"""
DESCRIPTION:
Select particles in centromere regions.

USAGE:
sele_cent obj genome

PARAMS:
obj (string or PyMOL selection)
        name of the PyMOL object/selection in which
        to find the centromere particles
genome (string)
        name of the reference genome
flank (int)
        number of flanking particles to select, 
        positive for right flank, negative for left flank, 0 returns nothing,
        default 1 (right 1 flanking particle)
RETURNS:
a newly created selection with name `obj`_cent.
"""

from bisect import bisect_right, bisect_left
from collections import namedtuple
import gzip
import csv
from pymol import cmd
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
def sele_cent(obj, genome, flank=1):
    flank = int(flank)
    SelName = obj + "_cent"

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
        atoms[atom[0]].append(int(atom[1]))
    atoms = {chrom : sorted(atoms[chrom]) for chrom in atoms}
    # tidy centromeres reference
    centromeres = fetch_centromeres(genome)
    # find atoms in centromeres and select 
    cmd.select(SelName, "None")
    if flank > 0:
        for chrom in atoms:
            #print(chrom)
            #print(atoms.keys())
            #print(centromeres.keys())
            if chrom not in centromeres:
                continue
            start = bisect_right(atoms[chrom], centromeres[chrom][1])
            for atom in atoms[chrom][start:start+flank]:
                cmd.select(SelName, SelName + " or (chain " + chrom + " and name " + str(atom) + " and obj " + obj +")")
    elif flank < 0:
        for chrom in atoms:
            if chrom not in centromeres:
                continue
            end = bisect_left(atoms[chrom], centromeres[chrom][0])
            for atom in atoms[chrom][end+flank:end]:
                cmd.select(SelName, SelName + " or (chain " + chrom + " and name " + str(atom) + " and obj " + obj +")")
    else:
        print("flank must be positive or negative")
    # select right flanking according to chrom, name
    # select left flanking according to chrom, name
    return SelName
cmd.extend("sele_cent", sele_cent)