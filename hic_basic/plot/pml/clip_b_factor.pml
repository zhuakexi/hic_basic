# pymol script to show clip b-factor
# Input:
#   $cif: .cif file
#   $png: output .png file

# --- profiles ---
set max_threads, 16
set connect_mode, 4

# --- load ---
load $cif, molecule

# --- shading ---
viewport 800, 800
clip slab, 1
clip move, $clip
set ray_shadows,0
set ambient, 1
set specular, off
set ray_opaque_background, off

# --- format molecule ---
as sticks, all
set_bond stick_radius, 0.25, all
bg_color white
spectrum b, $cmap

# --- rendering ---
ray
png $png, width=800, height=800

# --- quit ---
quit