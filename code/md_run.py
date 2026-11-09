# -*- coding: utf-8 -*-
"""OpenMM MD 生产（带 checkpoint / 断点续跑 / 100 ps 存帧）。

用法：
  python md_run.py --dir ... --ns 200
  python md_run.py --dir ... --ns 200 --resume   # 从 checkpoint 续跑
"""
import argparse
import json
import os
import shutil
import time

from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, Platform, Vec3, unit, app
from openmm.app import (CheckpointReporter, DCDReporter, ForceField, Modeller, PDBFile,
                        Simulation, StateDataReporter)
from openmm.app.element import Element
from openmmforcefields.generators import SystemGenerator
from openff.toolkit import Molecule
from openff.toolkit.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY

_orig = GLOBAL_TOOLKIT_REGISTRY.call


def _call(method, *a, **kw):
    if method == "assign_partial_charges" and kw.get("partial_charge_method") == "am1bcc":
        from openff.nagl_models import list_available_nagl_models
        kw["partial_charge_method"] = str(list_available_nagl_models()[-1])
    return _orig(method, *a, **kw)


GLOBAL_TOOLKIT_REGISTRY.call = _call


def load_mol(path, name):
    from rdkit import Chem
    mol = Chem.MolFromMolFile(path, removeHs=False)
    mol.SetProp("_Name", name)
    off = Molecule.from_rdkit(mol, allow_undefined_stereo=True)
    off.name = name
    return mol, off


def add_residue(modeller, sdf_path, resname):
    from rdkit import Chem
    mol = Chem.MolFromMolFile(sdf_path, removeHs=False)
    conf = mol.GetConformer()
    chain = modeller.topology.addChain()
    res = modeller.topology.addResidue(resname, chain)
    pos = list(modeller.positions)
    top_atoms = []
    for i, atom in enumerate(mol.GetAtoms()):
        el = Element.getByAtomicNumber(atom.GetAtomicNum())
        top_atoms.append(modeller.topology.addAtom(f"{el.symbol}{i + 1}", el, res))
        p = conf.GetAtomPosition(i)
        pos.append(Vec3(p.x / 10.0, p.y / 10.0, p.z / 10.0))
    for bond in mol.GetBonds():
        order = int(round(bond.GetBondTypeAsDouble()))
        order = order if order in (1, 2, 3) else 1
        modeller.topology.addBond(top_atoms[bond.GetBeginAtomIdx()],
                                  top_atoms[bond.GetEndAtomIdx()], order=order)
    modeller.positions = pos


def main():
    import pickle
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--ns", type=float, default=200.0)
    ap.add_argument("--equil-ps", type=float, default=100.0)
    ap.add_argument("--frame-ps", type=float, default=100.0)
    ap.add_argument("--padding-nm", type=float, default=1.0)
    ap.add_argument("--resume", action="store_true")
    args = ap.parse_args()
    d = args.dir
    ckpt = os.path.join(d, "checkpoint.chk")
    frame_steps = int(args.frame_ps / 0.002)
    resumed = args.resume and os.path.exists(ckpt)

    if resumed:
        # 续跑：加载 system/topology 与 checkpoint
        import pickle
        with open(os.path.join(d, "system.pkl"), "rb") as fh:
            system, topo = pickle.load(fh)
    else:
        prot = PDBFile(os.path.join(d, "protein_H.pdb"))
        nai_rd, nai = load_mol(os.path.join(d, "NAI.sdf"), "NAI")
        lig_rd, lig = load_mol(os.path.join(d, "ligand.sdf"), "LIG")
        sg = SystemGenerator(forcefields=["amber14-all.xml", "amber14/tip3pfb.xml"],
                             small_molecule_forcefield="openff-2.2.0",
                             molecules=[nai, lig],
                             periodic_forcefield_kwargs={"constraints": app.HBonds,
                                                         "nonbondedMethod": app.PME,
                                                         "nonbondedCutoff": 1.0 * unit.nanometer},
                             barostat=MonteCarloBarostat(1.0 * unit.bar, 300 * unit.kelvin, 25))
        sg.add_molecules([nai, lig])
        modeller = Modeller(prot.topology, prot.positions)
        modeller.addSolvent(ForceField("amber14-all.xml", "amber14/tip3pfb.xml"),
                            model="tip3p", padding=args.padding_nm * unit.nanometer,
                            ionicStrength=0.15 * unit.molar)
        add_residue(modeller, os.path.join(d, "NAI.sdf"), "NAI")
        add_residue(modeller, os.path.join(d, "ligand.sdf"), "LIG")
        system = sg.create_system(modeller.topology)
        topo = modeller.topology
        with open(os.path.join(d, "system.pkl"), "wb") as fh:
            pickle.dump((system, topo), fh)
        with open(os.path.join(d, "positions.pkl"), "wb") as fh:
            pickle.dump([tuple(p) for p in modeller.positions], fh)
        print("atoms", topo.getNumAtoms(), flush=True)

    integ = LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
    integ.setConstraintTolerance(1e-5)
    platform = Platform.getPlatformByName("OpenCL")
    props = {"Precision": "mixed"}
    sim = Simulation(topo, system, integ, platform, props)
    if args.resume and os.path.exists(ckpt):
        sim.loadCheckpoint(ckpt)
        print("resumed from checkpoint", flush=True)
    else:
        import pickle
        pos = pickle.load(open(os.path.join(d, "positions.pkl"), "rb"))
        sim.context.setPositions([Vec3(*p) for p in pos])
        sim.minimizeEnergy(maxIterations=1000)
        sim.step(int(args.equil_ps / 0.002))
        print("minimized+equil", flush=True)

    prod_steps = int(args.ns * 1e6 / 2)
    traj_path = os.path.join(d, "traj.dcd")
    if resumed and os.path.exists(traj_path):
        bak = os.path.join(d, "traj_backup_" + time.strftime("%Y%m%d_%H%M%S") + ".dcd")
        shutil.copy2(traj_path, bak)
        print("backup", bak, flush=True)
    sim.reporters.append(DCDReporter(traj_path, frame_steps, append=resumed))
    sim.reporters.append(StateDataReporter(os.path.join(d, "md.log"), frame_steps, step=True,
                                           potentialEnergy=True, temperature=True, density=True,
                                           speed=True, remainingTime=True, totalSteps=prod_steps,
                                           append=resumed))
    sim.reporters.append(CheckpointReporter(ckpt, frame_steps))
    t0 = time.time()
    sim.step(prod_steps)
    sim.saveCheckpoint(ckpt)
    dt = time.time() - t0
    out = dict(n_atoms=topo.getNumAtoms(), ns=args.ns, seconds=dt, ns_per_day=args.ns / dt * 86400)
    json.dump(out, open(os.path.join(d, "md_done.json"), "w"), indent=2)
    print(json.dumps(out), flush=True)


if __name__ == "__main__":
    main()
