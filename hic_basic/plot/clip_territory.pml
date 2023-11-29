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
# set color for each chain
set_color red, [1.0, 0.0, 0.0]
set_color green, [0.0, 1.0, 0.0]
set_color blue, [0.0, 0.0, 1.0]
set_color yellow, [1.0, 1.0, 0.0]
set_color cyan, [0.0, 1.0, 1.0]
set_color magenta, [1.0, 0.0, 1.0]
set_color orange, [1.0, 0.65, 0.0]
set_color purple, [0.5, 0.0, 0.5]
set_color lime, [0.75, 1.0, 0.0]
set_color pink, [1.0, 0.75, 0.8]
set_color teal, [0.0, 0.5, 0.5]
set_color gold, [1.0, 0.84, 0.0]
set_color navy, [0.0, 0.0, 0.5]
set_color gray, [0.5, 0.5, 0.5]
set_color maroon, [0.5, 0.0, 0.0]
set_color olive, [0.5, 0.5, 0.0]
set_color aqua, [0.0, 1.0, 1.0]
set_color fuchsia, [1.0, 0.0, 1.0]
colors = ["red", "green", "blue", "yellow", "cyan", "magenta","orange", "purple", "lime", "pink", "teal", "gold","navy", "gray", "maroon", "olive", "aqua", "fuchsia"]
chains = []
iterate (all), chains.append(chain)
python
for i, chain in enumerate(list(set(chains))):
    color = colors[i % len(colors)]  # repeat when not enough colors
    cmd.color(color, f"chain {chain}")
    print(chain, color)
python end
bg_color white
# --- rendering ---
ray
png $png, width=800, height=800

# --- quit ---
quit