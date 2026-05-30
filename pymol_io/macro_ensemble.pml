reinitialize

load /Users/miriamzara/MSM/backbone_reconstructions/state0_pcca6/state0_pcca6.pdb, ens

hide everything
show cartoon, ens

# Align all models to state 1
intra_fit ens, 1

# Show all states simultaneously
set all_states, on

# Make cartoons thinner
set cartoon_tube_radius, 0.15
set cartoon_loop_radius, 0.25
set cartoon_oval_width, 0.15
set cartoon_oval_length, 0.8


# set background to white
bg_color white

# Ensemble appearance
color gray70, ens
set cartoon_transparency, 0.6, ens

# Reference structure (state 1)
create ref, ens, 1, 1
set cartoon_transparency, 0.0, ref
color gray30, ref
set cartoon_tube_radius, 0.25, ref
set cartoon_loop_radius, 0.25, ref

# N- and C-termini on reference only
select nterm, ref and resi 0
select cterm, ref and resi 34

color blue, nterm
color red, cterm

# Optional: highlight termini more clearly
show sticks, nterm or cterm

orient ref
zoom ref



# nicer display
set cartoon_transparency, 0.0, ens
set cartoon_transparency, 0.0, ens and state 1
set all_states, on
orient ens


png pcca6_state0.png, width=2000, height=1500, dpi=300, ray=1
quit