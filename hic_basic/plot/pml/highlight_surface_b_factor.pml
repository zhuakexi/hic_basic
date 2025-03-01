# pymol script to show surface of molecules
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
set ray_shadows, 0

# --- format molecules ---
as sticks, all
set_bond stick_radius, 0.25, all
show sticks, chain $chain
spectrum b, $cmap

show sticks, not chain $chain
set transparency, 0.5, not chain $chain

bg_color white

# --- rotate camera ---
$turn_cmd

# --- rendering ---
ray
png $png, width=800, height=800

# --- quit ---
quit
