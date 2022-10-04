"""
DESCRIPTION:
Color spectrum according to relative position to centromeres.

USAGE:
color_centelo, selection, genome

PARAMS:
selection: selection string
genome: predefined genome string
"""
from pymol import cmd
import gzip
from collections import namedtuple
import csv
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
        "mm10" : "mm10_gap.csv.gz"
    }
    len_files = {
        "mm10" : "mm10.len.csv"
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
            fcsv = csv.reader(f)
            #header
            header = next(fcsv)
            # dtypes
            dtypes = [str, int]
            Len = namedtuple("Len", header)
            res = []
            for row in fcsv:
                row = Len(*(convert(value) for convert, value in zip(dtypes, row)))
                res.append(row)
    for row in res:
        if row in cent_chromlen:
            cent_chromlen[row.chromosome].append(row.length)
    print(cent_chromlen)
    return cent_chromlen
def color_centelo(obj, genome):
    # get centromeres
    centromeres = fetch_cent_chromlen(genome)
    # get all atoms (chrom, name)
    # cmd.iterate(selection str, expression str, space=)
    # expression string can use:
    #   1. default exposed varaibles
    #   2. stored.xxx
    #   3. custom variable space (space = myspace)
    myspace = {"bins":[]}
    cmd.iterate('(all)', "bins.append((chain, name))", space=myspace)
    atoms = {}
    for atom in myspace["bins"]:
        if atom[0] not in atoms:
            atoms[atom[0]] = []
        atoms[atom[0]].append(int(atom[1]))
    atoms = {chrom : sorted(atoms[chrom]) for chrom in atoms}

    # find atoms in centromeres and select 
    for chrom in atoms:
        if chrom not in centromeres:
            continue
        for atom in atoms[chrom]:
            if atom < centromeres[chrom][0]:
                # left
                relpos = (centromeres[chrom][0] - atom) / (centromeres[chrom][0] - 0 )
            elif atom > centromeres[chrom][1]:
                # right
                relpos = (atom - centromeres[chrom][1]) / (centromeres[chrom][2] - centromeres[chrom][1])
            else:
                # centromere
                relpos = 0
            cmd.alter("%s and chain %s and atom %d" % (obj, chrom, atom), "b=%f" % relpos)
    cmd.spectrum("b", "blue_white_red", "%s" % obj, minimum=0, maximum=1)
cmd.extend("color_centelo", color_centelo)