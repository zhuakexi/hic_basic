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

RETURNS:
a newly created selection with name `obj`_cent.
"""

from pymol import cmd

def sele_cent(obj, genome):
    SelName = obj + "_cent"
    pass
cmd.extend("sele_cent", sele_cent)