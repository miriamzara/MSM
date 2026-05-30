

reinitialize


load /Users/miriamzara/MSM/backbone_reconstructions/state5_pcca6/state5_pcca6.pdb, ens
load_traj /Users/miriamzara/MSM/backbone_reconstructions/state5_pcca6/state5_pcca6.dcd, ens

# after loading
remove not state 1+5+10+15+20+25+30+35+40

hide everything
intra_fit ens and name CA

# Show only first state as clean reference cartoon

create ref, ens, 1, 1
show cartoon, ref
color gray60, ref

# Show ensemble as CA trace only

show lines, ens and name CA
set all_states, on
color gray80, ens
set line_width, 1.5

# Highlight termini only on reference

select Nterm, ref and resi 0 and name CA
select Cterm, ref and resi 34 and name CA

show spheres, Nterm or Cterm
set sphere_scale, 0.45
color blue, Nterm
color red, Cterm
label Nterm, "N"
label Cterm, "C"
set label_size, 24
set label_color, black
bg_color white
orient ref
zoom ref