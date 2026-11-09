# -*- coding: utf-8 -*-
"""Concatenate the reconstructed 2.8-ns prefix with the stored 197-ns trajectory.

The boundary is not continuous (CA RMSD 1.39 A, ligand RMSD 2.52 A).  The
output is a *reconstructed* 200-ns display/QC trajectory, not the original
continuous production trajectory.
"""
import os

import mdtraj as md
from mdtraj.formats import DCDTrajectoryFile

ROOT = r"D:\scratch\protein-drug\md\entrectinib_4DRE"
PREFIX = r"D:\scratch\protein-drug\md\entrectinib_4DRE_prefix\traj.dcd"
CURRENT = os.path.join(ROOT, "traj.dcd")
OUT = os.path.join(ROOT, "traj_200ns_reconstructed.dcd")

with DCDTrajectoryFile(OUT, "w") as w:
    total = 0
    for path in (PREFIX, CURRENT):
        with DCDTrajectoryFile(path) as r:
            while total < 2000:
                xyz, cell_lengths, cell_angles = r.read(n_frames=100)
                if xyz is None or len(xyz) == 0:
                    break
                take = min(len(xyz), 2000 - total)
                w.write(xyz[:take], cell_lengths=None if cell_lengths is None else cell_lengths[:take],
                        cell_angles=None if cell_angles is None else cell_angles[:take])
                total += take
                if take < len(xyz):
                    break

with DCDTrajectoryFile(OUT) as f:
    xyz = f.read(n_frames=1)[0]
print("saved", OUT, "first-frame atoms", xyz.shape[1])

with open(os.path.join(ROOT, "SPLICE_NOTES.md"), "w", encoding="utf-8") as fh:
    fh.write("""# Reconstructed 200 ns trajectory

- Prefix: `entrectinib_4DRE_prefix/traj.dcd`, 28 frames / 2.8 ns, regenerated
  from `positions.pkl` after the original 0-2.8 ns were lost.
- Stored production segment: `traj.dcd`, 1,972 frames / 197.2 ns.
- Concatenated output: `traj_200ns_reconstructed.dcd`, 2,000 frames.

**Boundary is not continuous**: CA RMSD 1.39 A and ligand RMSD 2.52 A between
the prefix endpoint and the first stored frame. Use this file for display/QC
only. Quantitative production analyses should use the stored 197.2-ns
trajectory, or explicitly flag frames 0-27 as reconstructed.
""")
print("wrote SPLICE_NOTES.md")
