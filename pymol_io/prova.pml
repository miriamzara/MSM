# -----------------------------
# Load native-basin ensemble
# -----------------------------
load ../root/movies/tica_dihedrals_3state_macro0_frames_stride_1/first_frame.pdb, ensemble
load_traj ../root/movies/tica_dihedrals_3state_macro0_frames_stride_1/traj.dcd, ensemble

# -----------------------------
# Align all states on helix 2
# -----------------------------
# This fixes the segment before the Leu20-Pro21 loop,
# making differences around residues 19-22 easier to see.

select aln_ensemble, ensemble and name CA and resi 14-19
count_atoms aln_ensemble

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

# Full ensemble: subdued cartoons
show cartoon, ens_*
color gray80, ens_*
set cartoon_transparency, 0.70, ens_*
set cartoon_fancy_helices, 1
set cartoon_helix_radius, 0.18

# Representative frame
show cartoon, ens_0002
color gray50, ens_0002
set cartoon_transparency, 0.0, ens_0002

# -----------------------------
# Highlight the local loop
# -----------------------------
# Residues are numbered 0-34, so:
# Asn19 -> resi 18
# Leu20 -> resi 19
# Pro21 -> resi 20

select local_loop, ens_* and resi 17-22
show sticks, local_loop
set stick_radius, 0.13, local_loop

# Color the local backbone distinctly
color orange, ens_* and resi 18
color yelloworange, ens_* and resi 19
color magenta, ens_* and resi 20

# Keep side chains mainly on the representative frame
hide sticks, ens_* and resi 17-22 and not backbone
show sticks, ens_0002 and resi 17-22

# Explicitly emphasize backbone atoms around Leu20-Pro21
show spheres, ens_* and resi 19+20 and name CA
set sphere_scale, 0.25, ens_* and resi 19+20 and name CA

# Optional labels on representative frame
label ens_0002 and resi 18 and name CA, "Asn19"
label ens_0002 and resi 19 and name CA, "Leu20"
label ens_0002 and resi 20 and name CA, "Pro21"

set label_size, 18
set label_color, black
set label_outline_color, white

# -----------------------------
# Termini
# -----------------------------
select n_terminus_ref, ens_* and resi 0
select c_terminus_ref, ens_* and resi 34
color blue, n_terminus_ref
color red, c_terminus_ref

# -----------------------------
# Visual cleanup
# -----------------------------
set antialias, 2
set ray_opaque_background, off
set depth_cue, 0
set specular, 0.15
set shininess, 15
set orthoscopic, on

# Focus on the local loop rather than the whole protein
orient ens_0002 and resi 14-23
zoom ens_0002 and resi 14-23, 7

# Fast output
png out/tica_dihedrals_macro0_local_loop_fast.png, width=2200, height=1600, dpi=300, ray=0

# Optional high-quality render
#png out/tica_dihedrals_macro0_local_loop_ray.png, width=2200, height=1600, dpi=300, ray=1