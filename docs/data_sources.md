# Data sources and accessions

## Protein structures (RCSB PDB)

| PDB | Role |
|---|---|
| 4DRE | Wild-type InhA, NADH-bound; primary receptor |
| 4DTI | S94A InhA, NADH-bound; sensitivity receptor |
| 4U0J | Wild-type InhA, NAD+-bound; cofactor-state sensitivity receptor |
| 5JFO | InhA–GSK625–NAD+ co-crystal; orthogonal positive control |
| 4R9S | InhA–NITD-916 structure; positive-control audit |
| 2NSD | InhA–NAD+–4PI structure; arylamide control |

## Ligand and activity data

- **ChEMBL target CHEMBL1849** (InhA): activity records used to build the 337
  actives. Filtering: `standard_flag = 1`, `standard_relation = '='`,
  `standard_units = 'nM'`, IC50 present; deduplicate by InChIKey skeleton;
  median pIC50 for replicates. See `code/clean_inha_activities.py`.
- **ChEMBL approved drugs**: `max_phase = 4`, retrieved 2026-11-01; 4,225
  records reduced to 2,147 unique small molecules after SMILES/element/salt/
  duplicate/MW filters. See `code/fetch_chembl_approved.py`.
- **TDC tox21** pool: source of property-matched decoys (`code/build_decoys.py`).
  Matching: MW ±10%, cLogP ±1.0, HBD ±1, HBA ±1, rotatable bonds ±1, formal
  charge; ECFP4 Tanimoto < 0.3 to every active; no recorded InhA activity.
- **BindingDB/ChEMBL**: 566 InhA IC50 = 10.7 µM (reported as exploratory).
- **GSK625**: 2 nM InhA inhibitor from the 5JFO co-crystal (orthogonal control).

## Software

AutoDock Vina 1.2.7; meeko 0.8.0; RDKit 2026.09.1; Open Babel 3.1.0;
OpenMM 8.6.1; OpenFF Sage 2.2.0; openmmforcefields 0.15.1; ParmEd; mdtraj
1.11.2; MDAnalysis 2.11.0. Full versions and roles are in
`environment/README.md`.
