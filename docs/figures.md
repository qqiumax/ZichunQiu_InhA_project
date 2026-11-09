# Figure provenance

| Manuscript figure | File | Source data | Script |
|---|---|---|---|
| Graphical abstract | `fig_ga_docking_3d.pdf` | 4DRE docking pose + NADH | `render_inha_docking_3d.py`, `compose_docking_3d.py` |
| Fig. 1 | `fig1_workflow.pdf` | study design | `make_workflow_figure.py` |
| Fig. 2 | `fig1_retrospective_enrichment.pdf` | `enrichment_summary.csv`, `enrichment_curves.csv` | `make_publication_figures.py` |
| Fig. 3 | `fig2_docking_distributions.pdf` | `dock_distribution.csv`, `approved_hits.csv` | `make_publication_figures.py` |
| Fig. 4 | `validation_fig.pdf` | `positive_control_audit.csv`, `tier_validation_leave_scaffold.csv`, `property_stratified_enrichment.csv`, `approved_ligand_efficiency.csv` | `make_validation_figure.py` |
| Fig. 5 | `md_fig_rmsd_all.pdf` | `md_rmsd.csv`, `md_rmsd_extra.csv` | `make_md_unified_figures.py` |
| Fig. 6 | `md_fig_contacts_all.pdf` | `window5-20_*.csv` | `make_md_unified_figures.py` |
| Fig. 7 | `md_fig_mmgbsa_all.pdf` | `mmgbsa_summary_all.csv` | `make_md_unified_figures.py` |
| Fig. S1 | `fig6_docking_summary.pdf` | `approved_hits.csv`, `dock_thresholds.csv` | `make_publication_figures.py` |
| Fig. S2 | `fig7_confirmation.pdf` | `confirm_per_drug.csv`, `u0j_thresholds.csv` | `make_publication_figures.py` |
| Table S1 | `mmgbsa_566_supp.csv` | 566 MM/GBSA terms | `mmgbsa_ie_generic.py` |
