reinitialize

# -----------------------------
# Load native-basin ensemble
# -----------------------------
load ../root/movies/folded_core_consensus_frames_stride_100000/first_frame.pdb, ensemble
load_traj ../root/movies/folded_core_consensus_frames_stride_100000/traj.dcd, ensemble

# -----------------------------
# Align all ensemble states internally on the folded core
# -----------------------------

select aln_ensemble, ensemble and name CA and resi 6-10+14-19+22-29
count_atoms aln_ensemble

# Align all states to state 1 using the stable helical/core region
intra_fit aln_ensemble

# Split states after alignment
split_states ensemble, prefix=ens_
delete ensemble
delete aln_ensemble

# -----------------------------
# Display
# -----------------------------

hide everything, all
bg_color white

# Delete some frames if too cluttered
#delete ens_0001 or ens_0002 or ens_0003 or ens_0005 or ens_0006 or ens_0007 or ens_0008 or ens_0012

# All ensemble frames: thin, transparent cyan traces
show cartoon, ens_*
color gray80, ens_*
set cartoon_transparency, 0.55, ens_*
set cartoon_fancy_helices, 1
set cartoon_helix_radius, 0.2

# Representative frame: darker and less transparent

show cartoon, ens_0001
color gray80, ens_0001
set cartoon_transparency, 0.0, ens_0001
set cartoon_fancy_helices, 1
set cartoon_helix_radius, 0.2

# Highlight key 6--17 residues on representative frame
# numbered from 0, so i select 5--16

#select key_6_17, ens_* and resi 5+16
#color yelloworange, key_6_17



# Visual cleanup
set antialias, 2
set ray_opaque_background, off
set depth_cue, 0
set specular, 0.15
set shininess, 15
set orthoscopic, on

# Use a consistent view
orient ens_0001
zoom ens_0001, 4

# Fast output
png out/native_basin_ensemble_overlay_fast.png, width=2200, height=1600, dpi=300, ray=0

# Optional final high-quality render, only after the fast version looks right
png out/native_basin_ensemble_overlay_ray.png, width=2200, height=1600, dpi=300, ray=1