#!/usr/bin/env python3
"""
Reconstruct an approximate protein backbone from phi/psi torsion angles.

Example
-------
python reconstruct_backbone.py \
    --data-folder ./data \
    --frames 0 1000 2000 3000 \
    --out-prefix selected_backbone \
    --out-folder ./data/selected_frames

By default, frame indices refer to rows of the full dihedral file before any
subsampling. If --frames is omitted, the script uses every --stride-th row.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import mdtraj as md
import numpy as np
from tqdm import tqdm


# ------------------------------------------------------------
# Geometry helpers
# ------------------------------------------------------------

def place_atom(
    a: np.ndarray,
    b: np.ndarray,
    c: np.ndarray,
    length: float,
    angle_deg: float,
    dihedral_deg: float,
) -> np.ndarray:
    """
    Place atom d given previous atoms a, b, c.

    Geometry:
        |c-d| = length
        angle(b-c-d) = angle_deg
        dihedral(a-b-c-d) = dihedral_deg

    Coordinates are in Angstrom.
    """
    theta = np.deg2rad(angle_deg)
    phi = np.deg2rad(dihedral_deg)

    bc = c - b
    bc /= np.linalg.norm(bc)

    n = np.cross(b - a, bc)
    n_norm = np.linalg.norm(n)
    if n_norm < 1e-12:
        raise ValueError("Cannot place atom: previous three atoms are nearly collinear.")
    n /= n_norm

    m = np.cross(n, bc)

    d = c + length * (
        -np.cos(theta) * bc
        + np.sin(theta) * (np.cos(phi) * m + np.sin(phi) * n)
    )
    return d


def place_oxygen(N: np.ndarray, CA: np.ndarray, C: np.ndarray) -> np.ndarray:
    """Approximate backbone carbonyl oxygen position."""
    C_O = 1.229
    CA_C_O_angle = 120.8
    return place_atom(N, CA, C, C_O, CA_C_O_angle, 180.0)


# ------------------------------------------------------------
# Backbone reconstruction
# ------------------------------------------------------------

def reconstruct_backbone(
    phi: np.ndarray,
    psi: np.ndarray,
    omega: np.ndarray | None = None,
    n_residues: int = 35,
) -> np.ndarray:
    """
    Reconstruct N, CA, C, O backbone atoms from torsions.

    Parameters
    ----------
    phi, psi
        Arrays with shape (n_frames, n_residues - 2), in degrees.
        Columns are assumed to correspond to residues 2, ..., n_residues - 1.
    omega
        Optional omega angles, same shape as phi/psi. If omitted, omega = 180 deg.
    n_residues
        Number of residues.

    Returns
    -------
    xyz
        Array with shape (n_frames, n_residues * 4, 3), in Angstrom.
        Atom order per residue: N, CA, C, O.
    """
    phi = np.asarray(phi, dtype=float)
    psi = np.asarray(psi, dtype=float)

    if phi.ndim == 1:
        phi = phi[None, :]
    if psi.ndim == 1:
        psi = psi[None, :]

    expected_cols = n_residues - 2
    if phi.shape != psi.shape:
        raise ValueError(f"phi and psi must have the same shape, got {phi.shape} and {psi.shape}.")
    if phi.shape[1] != expected_cols:
        raise ValueError(
            f"Expected {expected_cols} torsion columns for n_residues={n_residues}, "
            f"got {phi.shape[1]}."
        )

    n_frames = phi.shape[0]

    if omega is None:
        omega = 180.0 * np.ones_like(phi)
    else:
        omega = np.asarray(omega, dtype=float)
        if omega.ndim == 1:
            omega = omega[None, :]
        if omega.shape != phi.shape:
            raise ValueError(f"omega must have shape {phi.shape}, got {omega.shape}.")

    # Standard peptide geometry in Angstrom and degrees.
    N_CA = 1.458
    CA_C = 1.525
    C_N = 1.329

    C_N_CA = 121.7
    N_CA_C = 111.2
    CA_C_N = 116.2

    xyz = np.zeros((n_frames, n_residues * 4, 3), dtype=float)

    for t in tqdm(range(n_frames), desc="Reconstructing backbone"):
        atoms: list[np.ndarray] = []

        # First residue in a fixed reference frame.
        N0 = np.array([0.0, 0.0, 0.0])
        CA0 = np.array([N_CA, 0.0, 0.0])

        angle = np.deg2rad(N_CA_C)
        C0 = CA0 + CA_C * np.array([
            np.cos(np.pi - angle),
            np.sin(np.pi - angle),
            0.0,
        ])

        O0 = place_oxygen(N0, CA0, C0)
        atoms.extend([N0, CA0, C0, O0])

        N_prev = N0
        CA_prev = CA0
        C_prev = C0

        # Build residues 2, ..., n_residues.
        for i in range(1, n_residues):
            # Torsion columns correspond to residues 2, ..., n_residues-1.
            if 1 <= i <= n_residues - 2:
                col = i - 1
                phi_i = phi[t, col]
                psi_prev = psi[t, col - 1] if col > 0 else 180.0
                omega_prev = omega[t, col - 1] if col > 0 else 180.0
            else:
                phi_i = 180.0
                psi_prev = psi[t, -1] if i == n_residues - 1 else 180.0
                omega_prev = omega[t, -1] if i == n_residues - 1 else 180.0

            N_i = place_atom(
                N_prev, CA_prev, C_prev,
                length=C_N,
                angle_deg=CA_C_N,
                dihedral_deg=psi_prev,
            )

            CA_i = place_atom(
                CA_prev, C_prev, N_i,
                length=N_CA,
                angle_deg=C_N_CA,
                dihedral_deg=omega_prev,
            )

            C_i = place_atom(
                C_prev, N_i, CA_i,
                length=CA_C,
                angle_deg=N_CA_C,
                dihedral_deg=phi_i,
            )

            O_i = place_oxygen(N_i, CA_i, C_i)
            atoms.extend([N_i, CA_i, C_i, O_i])

            N_prev, CA_prev, C_prev = N_i, CA_i, C_i

        xyz[t] = np.array(atoms)

    return xyz


# ------------------------------------------------------------
# MDTraj topology and IO
# ------------------------------------------------------------

def build_backbone_topology(n_residues: int = 35) -> md.Topology:
    topology = md.Topology()
    chain = topology.add_chain()

    for i in range(n_residues):
        residue = topology.add_residue(f"RES{i + 1}", chain)

        N = topology.add_atom("N", md.element.nitrogen, residue)
        CA = topology.add_atom("CA", md.element.carbon, residue)
        C = topology.add_atom("C", md.element.carbon, residue)
        O = topology.add_atom("O", md.element.oxygen, residue)

        topology.add_bond(N, CA)
        topology.add_bond(CA, C)
        topology.add_bond(C, O)

        if i > 0:
            prev_C = list(topology.residue(i - 1).atoms)[2]
            topology.add_bond(prev_C, N)

    return topology


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Reconstruct approximate N/CA/C/O backbone coordinates from a dihedral time series."
    )
    parser.add_argument(
        "--data-folder",
        type=Path,
        default=Path("./data"),
        help="Folder containing the dihedral file. Default: ./data",
    )
    parser.add_argument(
        "--dihedral-file",
        type=str,
        default="hp35.dihs",
        help="Dihedral filename inside --data-folder, or an absolute/relative path. Default: hp35.dihs",
    )
    parser.add_argument(
        "--frames",
        type=int,
        nargs="*",
        default=None,
        help="Frame indices to select from the full dihedral time series. Example: --frames 0 1000 2000",
    )
    parser.add_argument(
        "--frames-file",
        type=Path,
        default=None,
        help="Optional text file containing frame indices, one per line or whitespace-separated.",
    )
    parser.add_argument(
        "--stride",
        type=int,
        default=1000,
        help="Used only when --frames and --frames-file are omitted. Default: 1000",
    )
    parser.add_argument(
        "--n-residues",
        type=int,
        default=35,
        help="Number of residues. Default: 35",
    )
    parser.add_argument(
        "--out-prefix",
        type=str,
        default="backbone_reconstruction",
        help="Output prefix. Files are written as <prefix>.pdb and <prefix>.dcd. Default: backbone_reconstruction",
    )
    parser.add_argument(
        "--out-folder",
        type=Path,
        default=None,
        help="Output folder. Default: same as --data-folder",
    )
    parser.add_argument(
        "--no-align",
        action="store_true",
        help="Disable CA Kabsch alignment to the first selected frame.",
    )
    return parser.parse_args()


def load_frame_indices(args: argparse.Namespace, n_total_frames: int) -> np.ndarray:
    if args.frames is not None and args.frames_file is not None:
        raise ValueError("Use either --frames or --frames-file, not both.")

    if args.frames_file is not None:
        indices = np.loadtxt(args.frames_file, dtype=int)
        indices = np.atleast_1d(indices)
    elif args.frames is not None:
        indices = np.array(args.frames, dtype=int)
    else:
        if args.stride <= 0:
            raise ValueError("--stride must be positive.")
        indices = np.arange(0, n_total_frames, args.stride, dtype=int)

    if indices.size == 0:
        raise ValueError("No frames selected.")
    if np.any(indices < 0) or np.any(indices >= n_total_frames):
        bad = indices[(indices < 0) | (indices >= n_total_frames)]
        raise IndexError(
            f"Frame indices out of range [0, {n_total_frames - 1}]: {bad[:10]}"
        )

    return indices


def main() -> None:
    args = parse_args()

    dihedral_path = Path(args.dihedral_file)
    if not dihedral_path.is_absolute():
        dihedral_path = args.data_folder / dihedral_path

    parent_folder = args.out_folder if args.out_folder is not None else args.data_folder
    parent_folder.mkdir(parents=True, exist_ok=True)

    out_folder = parent_folder / args.out_prefix
    out_folder.mkdir(parents=True, exist_ok=True)

    dihedral_timeseries = np.loadtxt(dihedral_path, delimiter=" ", dtype=float)
    if dihedral_timeseries.ndim != 2:
        raise ValueError(f"Expected a 2D dihedral array, got shape {dihedral_timeseries.shape}.")
    if dihedral_timeseries.shape[1] % 2 != 0:
        raise ValueError(
            f"Expected alternating phi/psi columns, but found an odd number of columns: "
            f"{dihedral_timeseries.shape[1]}."
        )

    frame_indices = load_frame_indices(args, dihedral_timeseries.shape[0])
    dihedral_movie = dihedral_timeseries[frame_indices]

    phi_series = dihedral_movie[:, 0::2]
    psi_series = dihedral_movie[:, 1::2]
    omega_series = 180.0 * np.ones_like(phi_series)

    print(f"Loaded dihedrals: {dihedral_timeseries.shape}")
    print(f"Selected frames: {frame_indices.size}")
    print(f"First selected indices: {frame_indices[:10]}")
    print(f"phi_series: {phi_series.shape}")
    print(f"psi_series: {psi_series.shape}")

    xyz_angstrom = reconstruct_backbone(
        phi_series,
        psi_series,
        omega_series,
        n_residues=args.n_residues,
    )

    topology = build_backbone_topology(args.n_residues)
    traj = md.Trajectory(
        xyz=xyz_angstrom / 10.0,  # Angstrom -> nm
        topology=topology,
    )

    if not args.no_align:
        ca_indices = traj.topology.select("name CA")
        traj.superpose(traj[0], atom_indices=ca_indices)

    pdb_path = out_folder / f"{args.out_prefix}.pdb"
    dcd_path = out_folder / f"{args.out_prefix}.dcd"

    traj.save_pdb(str(pdb_path))
    traj.save_dcd(str(dcd_path))

    # Also save the exact selected source frame indices for reproducibility.
    index_path = out_folder / f"{args.out_prefix}_selected_frames.txt"
    np.savetxt(index_path, frame_indices, fmt="%d")

    print(f"Wrote: {pdb_path}")
    print(f"Wrote: {dcd_path}")
    print(f"Wrote: {index_path}")


if __name__ == "__main__":
    main()
