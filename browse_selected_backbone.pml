reinitialize

load /Users/miriamzara/MSM/backbone_reconstructions/state0_pcca6/state0_pcca6.pdb, recon
load_traj /Users/miriamzara/MSM/backbone_reconstructions/state0_pcca6/state0_pcca6.dcd, recon

hide everything
show cartoon, recon
show spheres, recon and name CA
color gray80, recon

set sphere_scale, 0.35, recon and name CA
bg_color white

select Nterm, recon and resi 0 and name CA
select Cterm, recon and resi 34 and name CA

color blue, Nterm
color red, Cterm

label Nterm, "N"
label Cterm, "C"

set label_size, 24
set label_color, black

set movie_fps, 10
mset 1 x9999

frame 1
orient recon
zoom recon
