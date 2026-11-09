# Validation of the tier criteria — 2026-11-09

Purpose: test whether the original A/B/C tier criteria discriminate known
binders from decoys, and whether the approved-drug A tier is driven by
molecular size/lipophilicity.

All docking used the same protocol (Vina 1.2.7, 22 Å box,
exhaustiveness 8, seed 42, num_modes 1). Thresholds are the original
4DRE/4DTI decoy p1/p5 and active medians.

## 1. Positive-control audit

| Compound | 4DRE | 4DTI | Original tier | Note |
|---|---:|---:|---|---|
| GSK625 | −8.265 | −7.868 | **none** | 2 nM; native 5JFO/NAD⁺ |
| NITD-916 | −8.250 | −8.187 | **C** | direct InhA inhibitor |
| 4PI | −8.329 | −8.145 | **C** | 2NSD arylamide |
| 566 | −8.075 | −8.094 | **none** | weak, 10.7 µM, cross-docked |

**Result:** no known strong binder reached the top or middle band. GSK625 (2 nM) failed even
C because 4DTI was worse than the C threshold. The original criteria cannot
recover strong positive controls in the receptors actually used.

GSK693 could not be included because no public structure/SMILES was readily
retrievable within this audit.

## 2. Scaffold-held-out tier validation

Leave-one-Murcko-scaffold-out (337 actives; 1,000-fold-free):

| Band | Held-out actives recovered | Decoy false rate |
|---|---:|---:|
| top | **2.4%** | 0.89% |
| middle | **7.7%** | 1.78% |
| lower | **32.0%** | 12.76% |

top-band enrichment over decoys is only ~2.7×, and sensitivity is 2.4%.
The criteria are not usable as a binder classifier.

## 3. MW/LogP-stratified enrichment

Selected 4DRE bins (AUC / EF5):

| Stratum | n act / dec | AUC | EF5 |
|---|---:|---:|---:|
| MW 0–300 | 27 / 28 | 0.779 | 2.04 |
| MW 300–400 | 183 / 187 | 0.698 | 1.57 |
| MW 400–500 | 87 / 94 | 0.696 | 1.39 |
| MW ≥500 | 39 / 28 | 0.665 | 1.15 |
| LogP 4–10 | 156 / 155 | 0.640 | 1.25 |

Discrimination is weakest in the high-MW/high-LogP region where the top-band
approved drugs sit.

## 4. Ligand efficiency of the A tier

Heavy-atom ligand efficiency (LE = −score / heavy atoms), drug-like approved
molecules (HA ≥ 15, MW ≥ 200; n = 1,830):

| Compound | LE | Percentile |
|---|---:|---:|
| Prednisolone tebutate | 0.295 | 55th |
| Entrectinib | 0.236 | 23rd |

The two top-band drugs are not high-LE outliers; their raw Vina scores are
consistent with a size/lipophilicity advantage.

## Bottom line

The original A tier is **not validated** as a binder classifier:
strong known inhibitors are missed, scaffold-held-out recovery is low, and
the top-band drugs are not ligand-efficiency outliers. The paper must not
present top-band membership as evidence of binding.

## Viable fixes (in order)

1. Replace A/B/C with a **continuous, size-aware ranking** (score + LE +
   MW/LogP stratification) and report precision only if experimental labels
   exist.
2. Use **GSK625/NITD-916 as positive controls** to calibrate any post-docking
   criterion; if they fail, the criterion is invalid.
3. Add **scaffold-held-out recovery** as the primary criterion diagnostic.
4. If a binary tier is required, derive thresholds from **held-out** data and
   report sensitivity/false-positive rate, not just counts.
5. Experimental IC50 remains the only way to validate a binder claim.
