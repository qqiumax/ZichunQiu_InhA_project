# Code map

| Script | Purpose |
|---|---|
| `prepare_receptors.py` | Receptor preparation and docking-box definition |
| `fetch_chembl_inha.py` | Retrieve InhA activity records from ChEMBL |
| `clean_inha_activities.py` | Build the 337-active calibration set |
| `build_decoys.py` | Property-matched decoy selection from TDC tox21 |
| `fetch_chembl_approved.py` | Build the 2,147 approved-drug library |
| `dock_batch.py` | AutoDock Vina docking |
| `analyze_approved.py` | Enrichment, tiers, and score summaries |
| `mutation_arm.py` | 4DTI S94A sensitivity analysis |
| `run_confirm.py` | Seed/exhaustiveness/4U0J confirmation |
| `tier_and_property_validation.py` | Scaffold-held-out validation, MW/LogP stratification, ligand efficiency |
| `make_validation_figure.py` | Validation figure (Fig. 4) |
| `contact_panel_window.py` | 5–20 ns contact occupancy |
| `md_prepare_poses.py`, `md_build.py`, `md_run.py` | MD system build and production |
| `mmgbsa_ie_generic.py` | MM/GBSA + interaction entropy |
| `final_md_stats.py` | Per-frame ligand RMSD |
| `make_publication_figures.py` | Main manuscript figures |
| `make_md_unified_figures.py` | Unified RMSD/contact/MM-GBSA figures |
| `render_inha_docking_3d.py`, `compose_docking_3d.py` | 3D graphical abstract |
| `splice_prefix.py` | Reconstructed 2.8-ns prefix (QC/display only) |

The scripts assume the original project directory layout. Paths to MD
trajectories can be replaced with local paths; the derived tables in
`../data/` allow the figure and validation scripts to run without trajectories.
