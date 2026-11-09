# -*- coding: utf-8 -*-
"""统一 5-20 ns 窗口的接触面板（entrectinib / 566 / metformin）。"""
import csv
import os
import pickle
import shutil
import struct

import mdtraj as md
import numpy as np

RAW = r"D:\workspace\research\protein-drug\04-data\raw\pdb"
RAW_MAP = {
    "GSK625": r"D:\scratch\protein-drug\md\GSK625_5JFO\protein_raw.pdb",
    "entrectinib_4DTI": r"D:\workspace\research\protein-drug\04-data\raw\pdb\4DTI.pdb",
    "prednisolone_4DTI": r"D:\workspace\research\protein-drug\04-data\raw\pdb\4DTI.pdb",
}
DIRS = {
    "entrectinib": r"D:\scratch\protein-drug\md\entrectinib_4DRE",
    "566": r"D:\scratch\protein-drug\md\ligand566_4DRE",
    "metformin": r"D:\scratch\protein-drug\md\metformin_4DRE",
    "decoy": r"D:\scratch\protein-drug\md\decoyPQDR_4DRE",
    "GSK625": r"D:\scratch\protein-drug\md\GSK625_5JFO",
    "entrectinib_4DTI": r"D:\scratch\protein-drug\md\entrectinib_4DTI",
    "prednisolone_4DTI": r"D:\scratch\protein-drug\md\prednisolone_tebutate_4DTI",
}
OUT = r"D:\workspace\research\protein-drug\06-figures\contact-panels"


def patch_dcd(src, dst, natoms):
    shutil.copyfile(src, dst)
    data = open(dst, "rb").read()
    rec = 4 * natoms
    first = data.find(struct.pack("<i", rec))
    start = first + 8 + rec
    frame = 56 + 3 * (8 + rec)
    frames = 0
    while start + (frames + 1) * frame <= len(data):
        if struct.unpack("<i", data[start + frames * frame:start + frames * frame + 4])[0] != rec:
            break
        frames += 1
    with open(dst, "r+b") as fh:
        fh.seek(8)
        fh.write(struct.pack("<i", frames))
    return frames


def analyse(d, tag):
    system, topo = pickle.load(open(os.path.join(d, "system.pkl"), "rb"))
    natoms = topo.getNumAtoms()
    raw = []
    raw_path = RAW_MAP.get(tag, os.path.join(RAW, "4DRE.pdb"))
    for line in open(raw_path, encoding="utf-8", errors="replace"):
        if line.startswith("ATOM") and line[21] == "A":
            k = (int(line[22:26]), line[17:20].strip())
            if not raw or raw[-1] != k:
                raw.append(k)
    prot = [r for r in topo.residues() if r.name not in ("LIG", "HOH", "NA", "CL")]
    names = []
    for i, r in enumerate(prot):
        names.append("NAI" if r.name == "NAI" else f"{raw[i][1]}{raw[i][0]}")
    lig = [a.index for a in topo.atoms() if a.residue.name == "LIG"]
    ca = [a.index for a in topo.atoms() if a.residue.name != "NAI" and a.name == "CA"]

    tmp = os.path.join(d, "traj_window.dcd")
    n = patch_dcd(os.path.join(d, "traj.dcd"), tmp, natoms)
    start, end = 50, min(n, 201)  # 5-20 ns, 100 ps/frame
    with md.formats.DCDTrajectoryFile(tmp) as fh:
        res = fh.read()
        xyz = np.asarray(res[0], dtype=np.float32)[start:end]
        box = np.asarray(res[1], dtype=np.float32)[start:end] if res[1] is not None else None
    if box is not None and box.shape[0] == xyz.shape[0]:
        for f in range(xyz.shape[0]):
            pc = xyz[f][ca].mean(0)
            xyz[f][lig] -= np.round((xyz[f][lig].mean(0) - pc) / box[f]) * box[f]
    mtop = md.Topology.from_openmm(topo)
    traj = md.Trajectory(xyz, mtop)
    traj.superpose(traj, frame=0, atom_indices=ca)
    X = traj.xyz
    lx = X[:, lig, :]
    rows = []
    for i, r in enumerate(prot):
        ids = [a.index for a in r.atoms()]
        dd = lx[:, :, None, :] - X[:, ids, :][:, None, :, :]
        dmin = np.sqrt((dd * dd).sum(-1)).min(axis=(1, 2))
        rows.append((names[i], 100 * np.mean(dmin <= 4.0), float(dmin.mean()), float(dmin.min())))
    p = os.path.join(OUT, f"window5-20_{tag}.csv")
    with open(p, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["residue", "occ_5_20ns_pct", "mean_min_dist_A", "min_dist_A"])
        w.writerows(rows)
    rows.sort(key=lambda x: -x[1])
    print(f"== {tag}: {xyz.shape[0]} frames ==")
    for x in rows[:15]:
        print(f"  {x[0]:8s} {x[1]:5.0f}%  mean {x[2]:.2f}")


for tag, d in DIRS.items():
    analyse(d, tag)
