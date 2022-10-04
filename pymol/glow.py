"""
DESCRIPTION:
Color selections and transparent others.

USAGE:
glow color, obj, selection

PARAMS:
color: list of string
    color names
selection list of string
    name of the selections to color
"""
from pymol import cmd
def glow(colors, selections):
    colors = colors.split()
    selections = selections.split()
    # store copy of the selection
    for color, selection in zip(colors, selections):
        # color select object
        cmd.color(color, selection)
        cmd.show_as("spheres", selection)
        cmd.set("sphere_scale", 0.3, selection)
        # transparent others
    cmd.show_as("surface", "not " + "(" + " or ".join(selections) + ")")
    cmd.color("gray", "not " + "(" + " or ".join(selections) + ")")
    cmd.set("transparency", 0.4, "not " + "(" + " or ".join(selections) + ")")
cmd.extend("glow", glow)