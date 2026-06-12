reinitialize

load ../root/data/pdb2f4k_cut.ent, pdb2f4k
alter pdb2f4k, resi=str(int(resi)-41)
sort

hide everything, all
show cartoon, pdb2f4k
color gray80, pdb2f4k

bg_color white
set cartoon_fancy_helices, 1
set cartoon_helix_radius, 0.5
set dash_width, 1
set dash_gap, 0.25
set dash_radius, 0.08
set label_size, 18
set label_position, [1.5, 2, 1.5]
set cartoon_transparency, 0.1

# Highlight aromatic side chains
select aromatics, pdb2f4k and resn PHE+TYR+TRP+HIS and not name N+C+O+CA
show sticks, aromatics
color orange, aromatics
set stick_radius, 0.28, aromatics


# Highlight mutated residues (Lys24Nle, Lys29Nle)
select norleucine, pdb2f4k and resn NLE and not name N+C+O+CA
show sticks, norleucine
color black, norleucine
set stick_radius, 0.1, norleucine



# Optional: show full aromatic residues slightly, including CA-CB connection
select aromatic_residues, pdb2f4k and resn PHE+TYR+TRP
show sticks, aromatic_residues and not name N+C+O
set stick_radius, 0.22, aromatic_residues and not name N+C+O


# Labels for aromatic and NLE residues
select special_labels, pdb2f4k and name CA and (resn PHE+TYR+TRP+HIS)
label special_labels, resi
set label_size, 18
set label_color, black


select chain_labels, pdb2f4k and name CA and resi 1+35
label chain_labels, resi