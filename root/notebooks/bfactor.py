import numpy as np
import mdtraj as md
import os
# -----------------------------
# Input files
# -----------------------------

topology_pdb = "../movies/folded_core_consensus_frames_stride_100/first_frame.pdb"
trajectory_dcd = "../movies/folded_core_consensus_frames_stride_100/traj.dcd"
#native_indices_file = "native_basin_frame_indices.txt"

os.makedirs("../intermediate_outputs/bfactor", exist_ok=True)
output_pdb = "../intermediate_outputs/bfactor/native_basin_rmsf_bfactor.pdb"
output_table = "../intermediate_outputs/bfactor/native_basin_rmsf_per_residue.csv"


# -----------------------------
# Load trajectory
# -----------------------------

native_traj = md.load_dcd(trajectory_dcd, top=topology_pdb)

#native_idx = np.loadtxt(native_indices_file, dtype=int)
#native_idx = native_idx[(native_idx >= 0) & (native_idx < traj.n_frames)]

#native_traj = traj[native_idx]

print("Total frames:", native_traj.n_frames)
print("Native-basin frames used:", native_traj.n_frames)


# -----------------------------
# Atom selections
# -----------------------------

top = native_traj.topology

# Since the trajectory is backbone-only, align using backbone atoms if present.
# For visualization and residue-wise RMSF, use C-alpha atoms.
align_atoms = top.select("backbone")
ca_atoms = top.select("name CA")

if len(ca_atoms) == 0:
    raise ValueError("No CA atoms found. Check atom names in the reconstructed topology.")

if len(align_atoms) == 0:
    print("No backbone selection found. Falling back to CA alignment.")
    align_atoms = ca_atoms

print("Alignment atoms:", len(align_atoms))
print("CA atoms:", len(ca_atoms))


# -----------------------------
# Align native-basin frames
# -----------------------------

# Use the first native-basin frame as reference.
# This removes rigid-body translation/rotation before computing RMSF.
native_traj.superpose(
    native_traj[0],
    atom_indices=align_atoms
)


# -----------------------------
# Compute C-alpha RMSF
# -----------------------------

ca_xyz = native_traj.xyz[:, ca_atoms, :]  # nm, shape = n_frames x n_CA x 3
mean_ca_xyz = ca_xyz.mean(axis=0)

rmsf_ca = np.sqrt(
    np.mean(
        np.sum((ca_xyz - mean_ca_xyz[None, :, :])**2, axis=2),
        axis=0
    )
)

# Normalize to 0--100 for PyMOL coloring
rmin = np.nanmin(rmsf_ca)
rmax = np.nanmax(rmsf_ca)

if rmax > rmin:
    rmsf_norm = 100.0 * (rmsf_ca - rmin) / (rmax - rmin)
else:
    rmsf_norm = np.zeros_like(rmsf_ca)


# -----------------------------
# Map CA RMSF to residues
# -----------------------------

residue_to_rmsf = {}

for atom_idx, rmsf_raw, rmsf_scaled in zip(ca_atoms, rmsf_ca, rmsf_norm):
    atom = top.atom(atom_idx)
    residue = atom.residue

    # MDTraj residue.index is zero-based.
    # residue.resSeq is the PDB residue number.
    residue_to_rmsf[residue.index] = {
        "resSeq": residue.resSeq,
        "resName": residue.name,
        "rmsf_nm": rmsf_raw,
        "rmsf_scaled": rmsf_scaled,
    }


# -----------------------------
# Save RMSF table
# -----------------------------

with open(output_table, "w") as f:
    f.write("residue_index,resSeq,resName,rmsf_nm,rmsf_scaled\n")

    for residue_index in sorted(residue_to_rmsf):
        d = residue_to_rmsf[residue_index]
        f.write(
            f"{residue_index},"
            f"{d['resSeq']},"
            f"{d['resName']},"
            f"{d['rmsf_nm']:.6f},"
            f"{d['rmsf_scaled']:.6f}\n"
        )

print("Saved:", output_table)


# -----------------------------
# Build representative structure
# -----------------------------

# Use the mean-aligned native-basin structure for a cleaner visualization.
representative = native_traj[0]
representative.xyz[0, ca_atoms, :] = mean_ca_xyz

# Save temporary PDB
tmp_pdb = "_tmp_native_mean.pdb"
representative.save_pdb(tmp_pdb)


# -----------------------------
# Rewrite B-factor column
# -----------------------------

# For each residue, assign the normalized CA RMSF to every atom in that residue.
# PDB B-factor column is columns 61--66.
with open(tmp_pdb, "r") as fin, open(output_pdb, "w") as fout:
    for line in fin:
        if line.startswith(("ATOM", "HETATM")):
            resseq = int(line[22:26])

            value = 0.0
            for d in residue_to_rmsf.values():
                if d["resSeq"] == resseq:
                    value = d["rmsf_scaled"]
                    break

            newline = line[:60] + f"{value:6.2f}" + line[66:]
            fout.write(newline)
        else:
            fout.write(line)

print("Saved:", output_pdb)
print("RMSF range, nm:", rmin, rmax)