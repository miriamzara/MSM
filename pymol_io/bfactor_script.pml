reinitialize

load ../root/intermediate_outputs/bfactor/native_basin_rmsf_bfactor.pdb, hp35_native_rmsf



hide everything, all
show cartoon, hp35_native_rmsf


bg_color white
set cartoon_fancy_helices, 1
set cartoon_helix_radius, 0.45
set cartoon_transparency, 0.1
set dash_width, 1
set dash_gap, 0.25
set dash_radius, 0.08
set label_size, 18
set label_position, [1.5, 2, 1.5]



# Color by normalized RMSF stored in B-factor.
# Blue = rigid, white = intermediate, red = flexible.
spectrum b, blue_white_red, hp35_native_rmsf, minimum=0, maximum=100

# Putty cartoon emphasizes flexible regions by thickness.
#set cartoon_putty, on
#set cartoon_putty_scale_min, 0.35
#set cartoon_putty_scale_max, 1.8
#set cartoon_putty_transform, 7

# Highlight key hydrophobic-core residues / contacts.
# Adjust residue numbering if your PDB is 0-based vs 1-based.
#select key_core, hp35_native_rmsf and resi 6+17 and name CA
#show spheres, key_core
#set sphere_scale, 0.45, key_core
#color yelloworange, key_core

# Optional: mark other hydrophobic-core residues if useful.
# select hydrophobic_core, hp35_native_rmsf and resi 5+6+9+16+17 and name CA
# show spheres, hydrophobic_core
# set sphere_scale, 0.35, hydrophobic_core
# color orange, hydrophobic_core

# N- and C-termini
#select Nterm, hp35_native_rmsf and resi 0 and name CA
#select Cterm, hp35_native_rmsf and resi 34 and name CA

#show spheres, Nterm
#show spheres, Cterm
#set sphere_scale, 0.35, Nterm
#set sphere_scale, 0.35, Cterm
#color blue, Nterm
#color red, Cterm


select chain_labels, pdb2f4k and name CA and resi 1+35
label chain_labels, resi

# Visual cleanup
set ray_opaque_background, off
set antialias, 2
set ray_trace_mode, 1
set depth_cue, 0
set specular, 0.25
set shininess, 20

orient hp35_native_rmsf
zoom hp35_native_rmsf, 4

png out/native_basin_rmsf_pymol.png, width=2000, height=1600, dpi=300, ray=1