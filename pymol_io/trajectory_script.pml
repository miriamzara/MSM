reinitialize

# topology first
load ../root/data/backbone_reconstruction_first_frame.pdb, recon

# load the trajectory into the existsing recon
load_traj ../root/data/backbone_reconstruction.dcd, recon

hide everything, all
bg_color white


# Reference frame
show cartoon, recon
set cartoon_fancy_helices, 1
set cartoon_helix_radius, 0.5
set cartoon_transparency, 0.1
#show spheres, recon and name CA
color gray80, recon
#set sphere_scale, 0.35, recon and name CA


# Highlight aromatic side chains
select aromatics, recon and resi 5+9+16 and name CA
show spheres, aromatics
color orange, aromatics
set sphere_scale, 0.45, aromatics



zoom recon
orient recon


select Nterm, recon and resi 0 and name CA
select Cterm, recon and resi 34 and name CA

color blue, Nterm
color red, Cterm

#label Nterm, "N"
#label Cterm, "C"

set label_size, 24
set label_color, black


set movie_fps, 10