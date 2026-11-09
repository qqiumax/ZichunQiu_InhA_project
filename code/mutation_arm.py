"""Mutation arm: does S94A (4DTI) weaken predicted binding relative to wild type (4DRE)?

Pre-registered criterion (03-methods §5):
    dds = score(4DTI) - score(4DRE)        positive = worse binding in the mutant
    PASS  if median(dds over actives) >= +0.5 kcal/mol AND Wilcoxon signed-rank p < 0.05
    else  competitive explanation 2 stands (S94A does not block inhibitor binding directly)

Also reports the decoy arm as a control: a mutation effect that hits decoys equally is not a
drug-specific effect.
"""
import json
import os

import numpy as np
import pandas as pd
from scipy import stats

PROJ = r"D:\workspace\research\protein-drug"
RES = os.path.join(PROJ, "05-analysis", "results")          # metrics outputs
DOCK = os.path.join(PROJ, "04-data", "derived", "docking")  # run outputs (S4)

wt = pd.read_csv(os.path.join(DOCK, "dock_4DRE_all.csv"))[["lig_id", "score", "group"]]
mu = pd.read_csv(os.path.join(DOCK, "dock_4DTI_all.csv"))[["lig_id", "score"]]
df = wt.merge(mu, on="lig_id", suffixes=("_wt", "_mu")).dropna()
df["dds"] = df["score_mu"] - df["score_wt"]

act = df[df["group"] == "active"]["dds"]
dec = df[df["group"] == "decoy"]["dds"]

out = {
    "n_actives": int(len(act)), "n_decoys": int(len(dec)),
    "actives": {
        "median_dds": round(float(act.median()), 3),
        "mean_dds": round(float(act.mean()), 3),
        "iqr": [round(float(act.quantile(0.25)), 3), round(float(act.quantile(0.75)), 3)],
        "wilcoxon_p": round(float(stats.wilcoxon(act).pvalue), 4),
    },
    "decoys": {
        "median_dds": round(float(dec.median()), 3),
        "mean_dds": round(float(dec.mean()), 3),
        "wilcoxon_p": round(float(stats.wilcoxon(dec).pvalue), 4),
    },
    "actives_vs_decoys": {
        "mannwhitney_p": round(float(stats.mannwhitneyu(act, dec, alternative="two-sided").pvalue), 4),
        "median_difference": round(float(act.median() - dec.median()), 3),
    },
}

crit = (out["actives"]["median_dds"] >= 0.5) and (out["actives"]["wilcoxon_p"] < 0.05)
out["preregistered_criterion"] = {
    "rule": "median dds(actives) >= +0.5 kcal/mol AND Wilcoxon p < 0.05",
    "passed": bool(crit),
    "verdict": ("S94A weakens predicted binding of known inhibitors" if crit
                else "criterion not met -> competitive explanation 2 (S94A does not directly block inhibitors) stands"),
}

# how many actives / decoys got worse vs better in the mutant
out["direction"] = {
    "actives_worse_in_mutant": int((act > 0).sum()), "actives_better": int((act < 0).sum()),
    "decoys_worse_in_mutant": int((dec > 0).sum()), "decoys_better": int((dec < 0).sum()),
}

df.to_csv(os.path.join(RES, "ddscore_per_ligand.csv"), index=False, encoding="utf-8-sig")
with open(os.path.join(RES, "ddscore_summary.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=2)

print(json.dumps(out, ensure_ascii=False, indent=2))
