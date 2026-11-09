# -*- coding: utf-8 -*-
"""构建带氢 MD 复合物：蛋白(PDBFixer/Modeller 加氢) + NAI(晶体坐标) + 配体(pose)。

用法：python md_build.py --rec 4DRE --pose-sdf ... --name entrectinib
输出（D:\\scratch\\protein-drug\\md\\<name>_<rec>\\）：
  complex_H.pdb  蛋白(含H)+NAI+LIG
  NAI.sdf        辅因子（含H）
  ligand.sdf     配体（含H，pose 坐标）
"""
import argparse
import os

import requests
from openmm import app, unit
from openmm.app import ForceField, Modeller, PDBFile
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")
PROJ = r"D:\workspace\research\protein-drug"
SCR = r"D:\scratch\protein-drug"
OUT = os.path.join(SCR, "md")
IDEAL = os.path.join(SCR, "NAI_ideal.sdf")


def nai_mol():
    if not os.path.exists(IDEAL):
        r = requests.get("https://files.rcsb.org/ligands/download/NAI_ideal.sdf", timeout=60)
        open(IDEAL, "wb").write(r.content)
    tpl = Chem.RemoveHs(Chem.MolFromMolFile(IDEAL, removeHs=False))
    lines = [ln for ln in open(os.path.join(PROJ, "04-data", "raw", "pdb", "4DRE.pdb"),
                               encoding="utf-8", errors="replace")
             if ln.startswith("HETATM") and ln[17:20].strip() == "NAI"]
    crystal = Chem.RemoveHs(Chem.MolFromPDBBlock("".join(lines) + "END\n", sanitize=False, removeHs=False))
    Chem.SanitizeMol(crystal)
    mol = AllChem.AssignBondOrdersFromTemplate(tpl, crystal)
    mol = Chem.AddHs(mol, addCoords=True)
    mol.SetProp("_Name", "NAI")
    return mol


def pdb_lines(mol, resname, chain, resseq, serial0=1):
    block = Chem.MolToPDBBlock(mol)
    out, serial = [], serial0
    for line in block.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            line = "HETATM" + f"{serial:5d}" + line[11:17] + f"{resname:>3}" + " " + chain + f"{resseq:4d}" + line[26:]
            out.append(line)
            serial += 1
    return out


def protein_with_h(rec):
    src = os.path.join(PROJ, "04-data", "raw", "pdb", f"{rec}.pdb")
    lines = [ln for ln in open(src, encoding="utf-8", errors="replace")
             if ln.startswith("ATOM") and ln[21] == "A"]
    tmp = os.path.join(SCR, "md", f"{rec}_prot.pdb")
    os.makedirs(os.path.dirname(tmp), exist_ok=True)
    open(tmp, "w", encoding="utf-8").write("\n".join(lines) + "\n")
    pdb = PDBFile(tmp)
    mod = Modeller(pdb.topology, pdb.positions)
    mod.addHydrogens(ForceField("amber14-all.xml"), pH=7.0)
    htmp = os.path.join(SCR, "md", f"{rec}_prot_H.pdb")
    with open(htmp, "w", encoding="utf-8") as fh:
        PDBFile.writeFile(mod.topology, mod.positions, fh)
    return [ln.rstrip("\n") for ln in open(htmp, encoding="utf-8", errors="replace")
            if ln.startswith(("ATOM", "HETATM"))]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rec", required=True)
    ap.add_argument("--pose-sdf", required=True)
    ap.add_argument("--name", required=True)
    args = ap.parse_args()

    lig = Chem.MolFromMolFile(args.pose_sdf, removeHs=False)
    lig.SetProp("_Name", "LIG")
    nai = nai_mol()

    outdir = os.path.join(OUT, f"{args.name}_{args.rec}")
    os.makedirs(outdir, exist_ok=True)
    Chem.MolToMolFile(nai, os.path.join(outdir, "NAI.sdf"))
    Chem.MolToMolFile(lig, os.path.join(outdir, "ligand.sdf"))

    prot = protein_with_h(args.rec)
    open(os.path.join(outdir, "protein_H.pdb"), "w", encoding="utf-8").write("\n".join(prot) + "\n")
    lines = prot + pdb_lines(nai, "NAI", "C", 1, 5000) + pdb_lines(lig, "LIG", "B", 1, 6000) + ["END"]
    p = os.path.join(outdir, "complex_H.pdb")
    open(p, "w", encoding="utf-8").write("\n".join(lines) + "\n")
    print("wrote", p)
    print("protein(H)", len(prot), "| NAI", nai.GetNumAtoms(), rdMolDescriptors.CalcMolFormula(nai),
          "| LIG", lig.GetNumAtoms(), rdMolDescriptors.CalcMolFormula(lig))


if __name__ == "__main__":
    main()
