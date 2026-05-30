#!/usr/bin/env python3

from pathlib import Path
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate a PyMOL script to browse a PDB+DCD trajectory frame by frame."
    )
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--dcd", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=Path("browse_backbone.pml"))
    parser.add_argument("--object-name", default="recon")
    parser.add_argument("--n-residues", type=int, default=35)
    parser.add_argument("--movie-fps", type=int, default=10)
    return parser.parse_args()


def main():
    args = parse_args()

    obj = args.object_name
    nterm_resi = 0
    cterm_resi = args.n_residues - 1

    pml = f"""reinitialize

load {args.pdb}, {obj}
load_traj {args.dcd}, {obj}

hide everything
show cartoon, {obj}
show spheres, {obj} and name CA
color gray80, {obj}

set sphere_scale, 0.35, {obj} and name CA
bg_color white

select Nterm, {obj} and resi {nterm_resi} and name CA
select Cterm, {obj} and resi {cterm_resi} and name CA

color blue, Nterm
color red, Cterm

label Nterm, "N"
label Cterm, "C"

set label_size, 24
set label_color, black

set movie_fps, {args.movie_fps}
mset 1 x9999

frame 1
orient {obj}
zoom {obj}
"""

    args.out.write_text(pml)
    print(f"Wrote PyMOL browser script: {args.out}")
    print(f"Open with: pymol {args.out}")


if __name__ == "__main__":
    main()