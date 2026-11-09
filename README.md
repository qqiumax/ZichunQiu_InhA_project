# Calibrated computational triage of approved drugs against Mtb InhA

This repository accompanies the manuscript:

> **Calibrated computational triage of approved drugs against *Mycobacterium
> tuberculosis* InhA: docking, dynamics, and the limits of binding inference**
> Zichun Qiu, Shanghai Starriver Bilingual School

The study asks whether a decoy-calibrated docking protocol can rank approved
drugs against InhA, and what MD and end-point free energy add to that ranking.
It is a computational prioritization study, not an inhibition or antibacterial
claim.

## What is included

| Directory | Contents |
|---|---|
| `data/` | Derived tables: docking scores, enrichment metrics, tier/band assignments, contact occupancies, RMSD, MM/GBSA terms, bootstrap intervals |
| `code/` | Analysis and figure scripts used for the manuscript |
| `structures/` | Public PDB structures and the 5JFO protein file used for the orthogonal control |
| `environment/` | Python/RDKit/OpenMM environment notes and pinned versions |
| `docs/` | Reproducibility notes, trajectory manifest, validation report |
| `MANIFEST.csv` | SHA-256 checksums for packaged files |

## What is not included

Raw MD trajectories are **not** stored here: the seven production trajectories
total about 1.9 GB. Their paths, frame counts, lengths, and sizes are listed in
`docs/trajectory_manifest.csv`; they are available from the corresponding
author on reasonable request.

Public raw inputs (full ChEMBL JSON pages, TDC tox21 download, all candidate
PDB files) are not duplicated here. Accessions and retrieval routes are listed
in `docs/data_sources.md`.

## Quick start

The manuscript uses two Python environments:

- **Main analysis** (Python 3.12.10): RDKit, meeko, pandas, numpy, scipy,
  matplotlib.
- **MD/MM-GBSA** (`D:\scratch\mdenv` in the original project): OpenMM 8.6.1,
  OpenFF Sage 2.2.0, openmmforcefields 0.15.1, ParmEd, mdtraj 1.11.2,
  MDAnalysis 2.11.0.

Docking also requires the AutoDock Vina 1.2.7 executable and Open Babel 3.1.0.
See `environment/README.md` for the exact versions.

To regenerate the manuscript figures from the packaged derived tables:

```bash
python code/make_publication_figures.py
python code/make_md_unified_figures.py
python code/make_validation_figure.py
```

To recompute the tier/LE validation from the packaged docking tables:

```bash
python code/tier_and_property_validation.py
```

The MD/MM-GBSA scripts require the trajectories, which are not in this
repository; `docs/trajectory_manifest.csv` gives the expected directory names
and file lengths.

## Key results

- Enrichment: AUC 0.684 (4DRE), 0.713 (4DTI), 0.724 (4U0J); EF5 1.59, 1.65,
  1.41; permutation *p* = 0.001.
- Top descriptive band: prednisolone tebutate and entrectinib.
- Entrectinib/4DRE: 197 ns stored trajectory; 5/5 key contacts over 5–20 ns;
  effective MM-GBSA −43.5 kcal/mol versus −17.7 for the decoy.
- Validation: GSK625 (2 nM) reached no band; NITD-916 reached only the lower
  band; scaffold-held-out top-band recovery was 2.4%. The top band is therefore
  a prioritization, not a validated binder classifier.

## Data and code availability

Derived data and code are released with the manuscript. The repository DOI is
https://github.com/qqiumax/ZichunQiu_InhA_project. Raw MD trajectories are available from the
corresponding author on reasonable request.

## License

Code is released under the MIT License (`LICENSE`). Author-created data,
figures, and documentation are released under CC BY 4.0 (`LICENSE-DATA`).
Third-party PDB and ChEMBL material remains subject to its own terms.

## Citation

```bibtex
@article{Qiu2026InhA,
  author  = {Qiu, Zichun},
  title   = {Calibrated computational triage of approved drugs against
             Mycobacterium tuberculosis InhA: docking, dynamics, and the
             limits of binding inference},
  year    = {2026},
  note    = {Manuscript; https://github.com/qqiumax/ZichunQiu_InhA_project}
}
```
