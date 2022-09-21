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
RETURNS:
a newly created selection with name `obj`_cent.
"""

from pymol import cmd
from .utils import fetch_centromeres

def sele_cent(obj, genome):
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
    
    # select
    cmd.select(SelName, "None")
    # select right flanking according to chrom, name
    # select left flanking according to chrom, name
    
    return SelName
cmd.extend("sele_cent", sele_cent)