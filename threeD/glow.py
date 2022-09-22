"""
DESCRIPTION:
Color select object and transparent others.

USAGE:
glow color, obj, selection

PARAMS:
color (string)
    color name
obj (string)
    name of the PyMOL object to store copy of the selection
selection (string)
    name of the selection to color
"""
from pymol import cmd
def glow(color, obj, selection):
    # store copy of the selection
    cmd.create(obj, selection)
    # color select object
    cmd.color(color, obj)
    cmd.show_as("spheres", obj)
    cmd.set("sphere_scale", 0.3, obj)
    # transparent others
    cmd.show_as("surface", "not " + obj)
    cmd.color("gray", "not " + obj)
    cmd.set("transparency", 0.4, "not " + obj)
cmd.extend("glow", glow)