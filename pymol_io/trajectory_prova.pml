reinitialize

# --- Load crystal/reference structure ---
load ../root/data/pdb2f4k_cut.ent, pdb2f4k
alter pdb2f4k, resi=str(int(resi)-41)
sort

# --- Load reconstructed trajectory ---
load ../root/data/backbone_reconstruction_first_frame.pdb, recon
load_traj ../root/data/backbone_reconstruction.dcd, recon

# --- Align trajectory object to crystal structure ---
align recon and name CA, pdb2f4k and name CA
intra_fit recon and name CA

hide everything, all
bg_color white

# --- Crystal structure in background ---
show cartoon, pdb2f4k
color gray70, pdb2f4k
set cartoon_transparency, 0.75, pdb2f4k

# Optional crystal aromatics / NLEs, very faint
select pdb_aromatics, pdb2f4k and resn PHE+TYR+TRP+HIS and not name N+C+O+CA
show sticks, pdb_aromatics
color gray50, pdb_aromatics
set stick_radius, 0.18, pdb_aromatics
set stick_transparency, 0.6, pdb_aromatics

select pdb_norleucine, pdb2f4k and resn NLE and not name N+C+O+CA
show sticks, pdb_norleucine
color black, pdb_norleucine
set stick_radius, 0.12, pdb_norleucine
set stick_transparency, 0.5, pdb_norleucine

# --- Moving trajectory on top ---
show cartoon, recon
color gray80, recon
set cartoon_fancy_helices, 1
set cartoon_helix_radius, 0.5
set cartoon_transparency, 0.1, recon

# Highlight aromatic positions in trajectory
select recon_aromatics, recon and resi 5+9+16 and name CA
show spheres, recon_aromatics
color orange, recon_aromatics
set sphere_scale, 0.45, recon_aromatics

# Termini
select Nterm, recon and resi 0 and name CA
select Cterm, recon and resi 34 and name CA

show spheres, Nterm or Cterm
set sphere_scale, 0.55, Nterm or Cterm
color blue, Nterm
color red, Cterm

# View
zoom recon or pdb2f4k
orient recon or pdb2f4k

set movie_fps, 1