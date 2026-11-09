# Derived data map

All files are CSV. They are sufficient to regenerate the main figures and the
reported summary values; raw trajectories are not included.

## Docking and calibration

| File | Contents |
|---|---|
| `dock_4DRE_all_joined.csv` | 4DRE scores for actives and decoys |
| `dock_4DTI_all_joined.csv` | 4DTI scores for actives and decoys |
| `dock_4U0J_all_joined.csv` | 4U0J sensitivity scores |
| `enrichment_summary.csv`, `enrichment_curves.csv` | AUC/EF/BEDROC summaries and curves |
| `dock_thresholds.csv`, `u0j_thresholds.csv` | Decoy/active percentile thresholds |
| `approved_hits.csv` | 2,147 approved drugs with scores, MW, cLogP, tier/band |
| `top_hits.csv` | Top-ranked approved drugs |
| `confirm_per_drug.csv`, `confirm_raw.csv` | Seed/exhaustiveness confirmation |
| `ddscore_per_ligand.csv`, `mutation_dds.csv` | 4DTI–4DRE mutation arm |
| `control_ex16_raw.csv`, `control_noise_raw.csv` | Control docking checks |

## MD, contacts, and MM/GBSA

| File | Contents |
|---|---|
| `md_rmsd.csv` | Per-frame ligand RMSD for the main systems |
| `md_rmsd_extra.csv` | Lead2/prednisolone RMSD |
| `md_rmsd_lead1_197.csv` | Full lead1 RMSD |
| `mmgbsa_summary.csv`, `mmgbsa_summary_all.csv` | MM/GBSA and IE summary |
| `mmgbsa_bootstrap.csv` | Block-bootstrap intervals |
| `mmgbsa_566_supp.csv` | Supplementary 566 values |
| `positive_control_audit.csv` | GSK625/NITD-916/4PI/566 audit |
| `tier_validation_leave_scaffold.csv` | Scaffold-held-out recovery |
| `property_stratified_enrichment.csv` | MW/LogP-stratified AUC/EF5 |
| `approved_ligand_efficiency.csv` | Heavy-atom ligand efficiency |

Contact occupancy CSVs are in the original project under
`06-figures/contact-panels/window5-20_*.csv`; the 566 panel is included in the
journal supplementary folder.
