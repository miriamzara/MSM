# Expected usage:
# pymol trajectory_script.pml movie_name
#
# More robust usage:
# pymol -cq trajectory_script.pml -- movie_name




reinitialize

# ---------- Start Python for parsing filename argument ----------
python
import os
import sys
from pymol import cmd

if len(sys.argv) < 2:
    raise SystemExit("ERROR: missing movie folder argument, e.g. stride_100")

# Take the last command-line argument as the movie folder name
movie_name = sys.argv[-1]

base_path = os.path.join("..", "root", "movies", movie_name)

topology_path = os.path.join(base_path, "first_frame.pdb")
trajectory_path = os.path.join(base_path, "traj.dcd")

cmd.load(topology_path, "recon")
cmd.load_traj(trajectory_path, "recon")
python end

# ---------- End Python for parsing filename argument ----------


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


# Highlight PHI10 dihedral in this 0-based trajectory numbering:
# residue 10 in 1-based numbering corresponds to resi 9 here.
# PHI10 = C(8)-N(9)-CA(9)-C(9)
select phi10_atoms, recon and ((resi 8 and name C) or (resi 9 and name N+CA+C))
select phi10_backbone, recon and resi 8+9 and name N+CA+C+O
show sticks, phi10_backbone
color magenta, phi10_backbone
set stick_radius, 0.22, phi10_backbone
show spheres, phi10_atoms
color magenta, phi10_atoms
set sphere_scale, 0.38, phi10_atoms

distance phi10_bond_1, recon and resi 8 and name C, recon and resi 9 and name N
distance phi10_bond_2, recon and resi 9 and name N, recon and resi 9 and name CA
distance phi10_bond_3, recon and resi 9 and name CA, recon and resi 9 and name C
hide labels, phi10_bond_*
color magenta, phi10_bond_*
set dash_width, 2.5, phi10_bond_*
set dash_gap, 0.15, phi10_bond_*

dihedral phi10_dihedral, recon and resi 8 and name C, recon and resi 9 and name N, recon and resi 9 and name CA, recon and resi 9 and name C
color magenta, phi10_dihedral
set dash_width, 2.5, phi10_dihedral
hide labels, phi10_dihedral

label recon and resi 9 and name CA, "phi10"
set label_color, magenta, recon and resi 9 and name CA


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


set movie_fps, 5