# -*- coding: utf-8 -*-
"""通用 MM-GBSA + IE：对任意 MD dir（system.pkl/traj.dcd/ligand.sdf/NAI.sdf）计算逐帧能量项与 IE。"""
import argparse
import csv
import json
import math
import os
import pickle
import shutil
import struct
import sys

sys.path.insert(0, r"D:\scratch\mmgbsa_entrectinib")

import MDAnalysis as mda
import numpy as np
import openmm
from openmm import VerletIntegrator, unit
from openmm.app import PDBFile
from openmm import Vec3

import compute_mmgbsa as cm


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


def write_pdb(d, outdir):
    system, topo = pickle.load(open(os.path.join(d, "system.pkl"), "rb"))
    pos = pickle.load(open(os.path.join(d, "positions.pkl"), "rb"))
    vecs = []
    for p in pos:
        v = []
        for x in p:
            v.append(float(x.value_in_unit(unit.nanometer)) if hasattr(x, "value_in_unit") else float(x))
        vecs.append(Vec3(*v))
    p = os.path.join(outdir, "complex_solvated.pdb")
    with open(p, "w", encoding="utf-8") as fh:
        PDBFile.writeFile(topo, vecs, fh)
    return topo.getNumAtoms()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--skip-ns", type=float, default=5.0)
    ap.add_argument("--end-ns", type=float, default=None)
    ap.add_argument("--stride", type=int, default=10)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    natoms = write_pdb(args.dir, args.out)
    dcd = os.path.join(args.out, "traj.dcd")
    n = patch_dcd(os.path.join(args.dir, "traj.dcd"), dcd, natoms)
    print("atoms", natoms, "frames", n, flush=True)

    lig_mol = cm.load_openff_molecule(os.path.join(args.dir, "ligand.sdf"))
    cof_mol = cm.load_openff_molecule(os.path.join(args.dir, "NAI.sdf"))
    complex_pdb, receptor_pdb, ligand_pdb = cm.strip_and_write_pdbs(
        os.path.join(args.out, "complex_solvated.pdb"), dcd,
        "resname LIG", "protein or resname NAI", __import__("pathlib").Path(args.out))
    cs, ct, _ = cm.build_implicit_system(str(complex_pdb), [lig_mol, cof_mol])
    rs, rt, _ = cm.build_implicit_system(str(receptor_pdb), [cof_mol])
    ls, lt, _ = cm.build_implicit_system(str(ligand_pdb), [lig_mol])

    def tag(system):
        out = {}
        for i, f in enumerate(system.getForces()):
            f.setForceGroup(i)
            out[i] = f.__class__.__name__
        return out

    nc, nr, nl = tag(cs), tag(rs), tag(ls)
    platform = openmm.Platform.getPlatformByName("CPU")
    ctx_c = openmm.Context(cs, VerletIntegrator(0.001 * unit.picosecond), platform)
    ctx_r = openmm.Context(rs, VerletIntegrator(0.001 * unit.picosecond), platform)
    ctx_l = openmm.Context(ls, VerletIntegrator(0.001 * unit.picosecond), platform)

    def terms(ctx, names, pos):
        ctx.setPositions(pos)
        out = {}
        for i, nm in names.items():
            out[nm] = out.get(nm, 0.0) + ctx.getState(getEnergy=True, groups={i}).getPotentialEnergy().value_in_unit(
                unit.kilocalorie_per_mole)
        return out

    u = mda.Universe(os.path.join(args.out, "complex_solvated.pdb"), dcd)
    solute = "not (resname HOH WAT TIP3 SOL NA CL Na+ Cl- K K+ NA+ CL-)"
    c_at = u.select_atoms(solute)
    r_at = u.select_atoms("protein or resname NAI")
    l_at = u.select_atoms("resname LIG")
    prot_ca = u.select_atoms("protein and name CA")
    lig = u.select_atoms("resname LIG")
    cof = u.select_atoms("resname NAI")

    start = int(args.skip_ns * 10)  # 100 ps/frame -> 10 frames/ns
    end = n if args.end_ns is None else min(n, int(args.end_ns * 10) + 1)
    rows = []
    for fi in range(start, end, args.stride):
        u.trajectory[fi]
        box = u.dimensions[:3] if u.dimensions is not None else None
        if box is not None and np.all(box > 0):
            pc = prot_ca.center_of_mass()
            for sel in (lig, cof):
                sel.positions -= np.round((sel.positions - pc) / box) * box
        ec = terms(ctx_c, nc, c_at.positions / 10.0)
        er = terms(ctx_r, nr, r_at.positions / 10.0)
        el = terms(ctx_l, nl, l_at.positions / 10.0)

        def agg(e):
            a = {}
            for k, v in e.items():
                if "Nonbonded" in k:
                    a["nonbonded"] = a.get("nonbonded", 0.0) + v
                elif "GB" in k:
                    a["GB"] = a.get("GB", 0.0) + v
                elif any(x in k for x in ("Bond", "Angle", "Torsion", "Improper")):
                    a["bonded"] = a.get("bonded", 0.0) + v
            return a

        ac, ar, al = agg(ec), agg(er), agg(el)
        de = ac.get("nonbonded", 0) - ar.get("nonbonded", 0) - al.get("nonbonded", 0)
        dg = ac.get("GB", 0) - ar.get("GB", 0) - al.get("GB", 0)
        db = ac.get("bonded", 0) - ar.get("bonded", 0) - al.get("bonded", 0)
        rows.append((fi, fi * 0.1, db, de, dg, db + de + dg))
        if len(rows) % 5 == 1:
            print(f"  frame {fi} ({fi*0.1:.1f} ns) dG {db+de+dg:.1f}", flush=True)

    with open(os.path.join(args.out, "mmgbsa_terms.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["frame", "time_ns", "dE_bonded", "dE_nonbonded", "dG_GB", "dG_total"])
        w.writerows(rows)

    e = np.array([r[3] for r in rows])
    g = np.array([r[5] for r in rows])
    kBT = 0.0019872041 * 300.0
    beta = 1.0 / kBT
    c = e - e.mean()
    mx = np.max(beta * c)
    ie = kBT * (mx + math.log(np.mean(np.exp(beta * c - mx))))
    out = {"n_frames": len(rows), "dE_nonbonded_mean": float(e.mean()), "dE_nonbonded_sd": float(e.std(ddof=1)),
           "dG_mean": float(g.mean()), "dG_sd": float(g.std(ddof=1)), "IE": float(ie),
           "dG_effective": float(g.mean() + ie)}
    json.dump(out, open(os.path.join(args.out, "mmgbsa_ie.json"), "w"), indent=2)
    print(json.dumps(out, indent=2), flush=True)


if __name__ == "__main__":
    main()
