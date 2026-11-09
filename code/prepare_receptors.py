"""Prepare InhA receptors for Vina docking (v2).

Fixes vs v1:
  * Superimposer direction was inverted. Bio.PDB's rotran maps MOVING onto FIXED, so to move the
    4U0J reference ligand into a target's crystal frame we must pass fixed=target, moving=4U0J.
    v2 verifies this by measuring the distance between the transformed ligand centroid and the
    target's own NADH centroid (should be a few Angstrom, not tens).
  * meeko refuses the NAD/NADH HETATM residue ("not in residue_templates"). Fallback:
    prepare the protein-only receptor with meeko, prepare NADH separately with obabel, then
    concatenate the two PDBQT files into one receptor file (standard workaround).
"""
import os
import subprocess

import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select, Superimposer

RAW = r"D:\workspace\research\protein-drug\04-data\raw\pdb"
PREP = r"D:\workspace\research\protein-drug\04-data\derived\receptor"
os.makedirs(PREP, exist_ok=True)
parser = PDBParser(QUIET=True)
COFACTORS = ("NAI", "NAD")
BOX = 22.0


class KeepProteinOnly(Select):
    def accept_residue(self, residue):
        return residue.get_parent().id == "A" and residue.id[0] == " "


def model_of(pdb_id):
    return parser.get_structure(pdb_id, os.path.join(RAW, f"{pdb_id}.pdb"))[0]


def ca_map(model, chain="A"):
    return {res.id[1]: res["CA"] for res in model[chain] if res.id[0] == " " and "CA" in res}


def het_atoms(model, resnames=COFACTORS, chain="A"):
    out = []
    for res in model[chain]:
        if res.get_resname() in resnames:
            out.extend(list(res))
    return out


def run(cmd):
    r = subprocess.run(cmd, capture_output=True, encoding="utf-8", errors="replace")
    return r.returncode, (r.stdout or ""), (r.stderr or "")


ref = model_of("4U0J")
ref_lig = [a for res in ref["A"] if res.get_resname() == "566" for a in res]
ref_lig_xyz = np.array([a.coord for a in ref_lig])
ref_ca = ca_map(ref)
print(f"reference ligand 566 from 4U0J: {len(ref_lig)} atoms")

summary = {}
for target in ("4DRE", "4DTI"):
    tgt = model_of(target)
    tgt_ca = ca_map(tgt)
    common = sorted(set(ref_ca) & set(tgt_ca))
    sup = Superimposer()
    # fixed = target (destination), moving = 4U0J (source) -> rotran maps 4U0J onto target
    sup.set_atoms([tgt_ca[i] for i in common], [ref_ca[i] for i in common])
    rot, tran = sup.rotran
    xyz = np.dot(ref_lig_xyz, rot) + tran
    center = xyz.mean(axis=0)

    cof = het_atoms(tgt)
    cof_xyz = np.array([a.coord for a in cof]) if cof else None
    d_ref = float(np.linalg.norm(center - cof_xyz.mean(axis=0))) if cof_xyz is not None else float("nan")
    print(f"\n{target}: CA RMSD {sup.rms:.3f} A over {len(common)} residues; "
          f"box center {center.round(2).tolist()}; cofactor atoms {len(cof)}; "
          f"box-to-cofactor-centroid distance {d_ref:.1f} A")

    # protein-only receptor
    prot = os.path.join(PREP, f"{target}_prot.pdb")
    io = PDBIO(); io.set_structure(tgt)
    io.save(prot, KeepProteinOnly())
    base = os.path.join(PREP, target)
    rc, out, err = run(["mk_prepare_receptor.exe", "--read_pdb", prot, "-o", base,
                        "-p", base + "_prot.pdbqt", "-v", base + "_vina_box.txt",
                        "--box_center", f"{center[0]:.3f}", f"{center[1]:.3f}", f"{center[2]:.3f}",
                        "--box_size", str(BOX), str(BOX), str(BOX), "--default_altloc", "A"])
    print(f"  meeko protein-only rc={rc}")

    # cofactor PDBQT via obabel
    cof_pdb = os.path.join(PREP, f"{target}_NADH.pdb")
    with open(cof_pdb, "w", encoding="utf-8") as fh:
        for a in cof:
            fh.write(f"HETATM{a.get_serial_number():5d} {a.get_name():<4s}{a.get_parent().get_resname():>4s} A"
                     f"{a.get_parent().id[1]:4d}    "
                     f"{a.coord[0]:8.3f}{a.coord[1]:8.3f}{a.coord[2]:8.3f}  1.00  0.00\n")
    cof_pdbqt = os.path.join(PREP, f"{target}_NADH.pdbqt")
    rc2, out2, err2 = run(["obabel", cof_pdb, "-O", cof_pdbqt, "-h", "--partialcharge", "gasteiger"])
    print(f"  obabel NADH rc={rc2} {os.path.getsize(cof_pdbqt) if os.path.exists(cof_pdbqt) else 0} B")

    # concatenate
    full = os.path.join(PREP, f"{target}_receptor_full.pdbqt")
    with open(full, "w", encoding="utf-8") as out_fh:
        for part in (base + "_prot.pdbqt", cof_pdbqt):
            if os.path.exists(part):
                with open(part, encoding="utf-8", errors="replace") as fh:
                    for line in fh:
                        if line.startswith(("ATOM", "HETATM")):
                            out_fh.write(line)
    n_atoms = sum(1 for line in open(full, encoding="utf-8", errors="replace") if line.startswith(("ATOM", "HETATM")))
    print(f"  merged receptor -> {os.path.basename(full)}  atoms={n_atoms}")
    summary[target] = dict(center=center.round(3).tolist(), rmsd=round(sup.rms, 3),
                           cofactor_atoms=len(cof), box_to_cofactor_A=round(d_ref, 2),
                           receptor_atoms=n_atoms)

print("\n=== summary ===")
for k, v in summary.items():
    print(k, v)
