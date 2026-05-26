

reinitialize

# topology first
load data/backbone_reconstruction_first_frame.pdb, recon

# load the trajectory into the existsing recon
load_traj data/backbone_reconstruction.dcd, recon

hide everything
show cartoon, recon
show spheres, recon and name CA
color gray80, recon

set sphere_scale, 0.35, recon and name CA
bg_color white

zoom recon
orient recon



select Nterm, recon and resi 0 and name CA
select Cterm, recon and resi 34 and name CA

color blue, Nterm
color red, Cterm

label Nterm, "N"
label Cterm, "C"

set label_size, 24
set label_color, black


set movie_fps, 10