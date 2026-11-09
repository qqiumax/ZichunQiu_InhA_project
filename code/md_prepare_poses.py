import os

from meeko import PDBQTMolecule, RDKitMolCreate
from rdkit import Chem

SCR = r"D:\scratch\protein-drug"
OUT = os.path.join(SCR, "md", "poses")
os.makedirs(OUT, exist_ok=True)

JOBS = [
    ("4DRE", "M_HAYYBYPASCDWEQ-UHFFFAOYSA-N", "entrectinib"),
    ("4DTI", "M_HAYYBYPASCDWEQ-UHFFFAOYSA-N", "entrectinib"),
    ("4DRE", "M_XZWYZXLIPXDOLR-UHFFFAOYSA-N", "metformin"),
    ("4DTI", "M_HUMXXHTVHHLNRO-KAJVQRHHSA-N", "prednisolone_tebutate"),
]

for rec, lig, name in JOBS:
    p = os.path.join(SCR, "confirm", f"{rec}_{lig}_s42_e16.pdbqt")
    if not os.path.exists(p):
        p = os.path.join(SCR, "confirm", f"{rec}_{lig}_s42_e8.pdbqt")
    if not os.path.exists(p):
        p = os.path.join(SCR, "docking", f"{rec}_{lig}.pdbqt")
    if not os.path.exists(p):
        print("MISSING", p)
        continue
    pm = PDBQTMolecule.from_file(p, skip_typing=True)
    mol = RDKitMolCreate.from_pdbqt_mol(pm)[0]
    out = os.path.join(OUT, f"{name}_{rec}.sdf")
    Chem.MolToMolFile(mol, out)
    print("ok", out, mol.GetNumAtoms())
