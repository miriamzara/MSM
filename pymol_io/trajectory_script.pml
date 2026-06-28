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