# -*- coding: utf-8 -*-
import os, pickle, shutil, struct
import mdtraj as md
import numpy as np

DIRS = {
    "entrectinib_4DTI": (r"D:\scratch\protein-drug\md\entrectinib_4DTI", 5),
    "prednisolone_4DTI": (r"D:\scratch\protein-drug\md\prednisolone_tebutate_4DTI", 1),
}
rows = []
for tag, (d, stride) in DIRS.items():
    system, topo = pickle.load(open(os.path.join(d, "system.pkl"), "rb"))
    natoms = topo.getNumAtoms()
    lig = [a.index for a in topo.atoms() if a.residue.name == "LIG"]
    ca = [a.index for a in topo.atoms() if a.residue.name != "NAI" and a.name == "CA"]
    src = os.path.join(d, "traj.dcd"); tmp = os.path.join(d, "traj_stats.dcd")
    shutil.copyfile(src, tmp)
    data = open(tmp, "rb").read(); rec = 4 * natoms
    first = data.find(struct.pack("<i", rec)); start = first + 8 + rec
    frame = 56 + 3 * (8 + rec); n = 0
    while start + (n + 1) * frame <= len(data):
        if struct.unpack("<i", data[start + n * frame:start + n * frame + 4])[0] != rec:
            break
        n += 1
    with open(tmp, "r+b") as fh:
        fh.seek(8); fh.write(struct.pack("<i", n))
    with md.formats.DCDTrajectoryFile(tmp) as fh:
        res = fh.read(stride=stride)
        xyz = np.asarray(res[0], dtype=np.float32)
        box = np.asarray(res[1], dtype=np.float32) if res[1] is not None else None
    if box is not None and box.shape[0] == xyz.shape[0]:
        for f in range(xyz.shape[0]):
            pc = xyz[f][ca].mean(0)
            xyz[f][lig] -= np.round((xyz[f][lig].mean(0) - pc) / box[f]) * box[f]
    mtop = md.Topology.from_openmm(topo)
    traj = md.Trajectory(xyz, mtop)
    traj.superpose(traj, frame=0, atom_indices=ca)
    lx = traj.xyz[:, lig, :]
    diff = lx - lx[0][None, :, :]
    r = np.sqrt((diff * diff).sum(2).mean(1))
    t = np.arange(len(r)) * stride * 0.1
    print(tag, "frames", len(r), "ns", round(t[-1], 1),
          "mean", round(float(r.mean()), 2), "max", round(float(r.max()), 2),
          "last10", round(float(r[t >= t[-1] - 10].mean()), 2))
    for ti, ri in zip(t, r):
        rows.append((tag, round(float(ti), 2), round(float(ri), 4)))

out = r"D:\workspace\research\protein-drug\06-figures\data\md_rmsd_extra.csv"
with open(out, "w", newline="") as fh:
    fh.write("system,time_ns,rmsd_A\n")
    for tag, ti, ri in rows:
        fh.write(f"{tag},{ti},{ri}\n")
print("saved", out)
