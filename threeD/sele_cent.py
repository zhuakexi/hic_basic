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
from pymol import cmd
from .utils import fetch_centromeres

def sele_cent(obj, genome, flank=1):
    SelName = obj + "_centromeres"

    # get all atoms (chrom, name)
    # cmd.iterate(selection str, expression str, space=)
    # expression string can use:
    #   1. default exposed varaibles
    #   2. stored.xxx
    #   3. custom variable space (space = myspace)
    myspace = {"bins":[]}
    cmd.iterate('(all)', bins.append((chain, name)), space=myspace)
    
    # find atoms (chrom, name) in centromeres
    # tidy atoms
    atoms = {}
    for atom in myspace["bins"]:
        if atom[0] not in atoms:
            atoms[atom] = []
        else:
            atoms[atom].append(int(atom[1]))
    atoms = {chrom : sorted(atoms[chrom]) for chrom in atoms}
    # tidy centromeres reference
    centromeres = fetch_centromeres(genome)
    centromeres = {row.chrom : [row.start, row.end] for row in centromeres}
    # find atoms in centromeres and select
    cmd.select(SelName, "None")
    if flank > 0:
        for chrom in atoms:
            start = bisect_right(atoms[chrom], centromeres[chrom][1])
            for atom in atoms[chrom][start:start+flank]:
                cmd.select(SelName, SelName + " or (chain " + chrom + " and name " + str(atom) + ")")
    elif flank < 0:
        for chrom in atoms:
            end = bisect_left(atoms[chrom], centromeres[chrom][0])
            for atom in atoms[chrom][end+flank:end]:
                cmd.select(SelName, SelName + " or (chain " + chrom + " and name " + str(atom) + ")")
    else:
        print("flank must be positive or negative")
    # select right flanking according to chrom, name
    # select left flanking according to chrom, name
    return SelName
cmd.extend("sele_cent", sele_cent)