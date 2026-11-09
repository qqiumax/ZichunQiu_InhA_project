# -*- coding: utf-8 -*-
"""Scaffold-held-out tier validation + MW/LogP-stratified enrichment + LE."""
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from rdkit import Chem
from rdkit.Chem import Descriptors

ROOT = r"D:\workspace\research\protein-drug"
RES = os.path.join(ROOT, "05-analysis", "results")
OUT = os.path.join(ROOT, "06-figures", "data")

# ---- active/decoys scores ----
a4 = pd.read_csv(os.path.join(RES, "dock_4DRE_all_joined.csv"))
b4 = pd.read_csv(os.path.join(RES, "dock_4DTI_all_joined.csv"))
m = a4.merge(b4, on=["lig_id", "group"], suffixes=("_4d", "_4t"))
act = m[m.group == "active"][["lig_id", "score_4d", "score_4t", "pIC50_4d"]].copy()
dec = m[m.group == "decoy"][["lig_id", "score_4d", "score_4t"]].copy()

# ---- scaffolds ----
sc = pd.read_csv(os.path.join(ROOT, "04-data", "derived", "inha_ic50_clean.csv"))
sc["lig_id"] = "A_" + sc["inchikey_skel"]
act = act.merge(sc[["lig_id", "murcko_scaffold"]], on="lig_id", how="left")
act["murcko_scaffold"] = act["murcko_scaffold"].fillna("UNASSIGNED")

# ---- properties ----
props = pd.read_csv(os.path.join(ROOT, "04-data", "derived", "decoys_tox21.csv"))
aprop = props[["active_id", "a_mw", "a_clogp"]].drop_duplicates("active_id")
dprop = props[["decoy_id", "d_mw", "d_clogp"]].drop_duplicates("decoy_id")
act = act.merge(aprop, left_on="lig_id", right_on="active_id", how="left")
dec = dec.merge(dprop, left_on="lig_id", right_on="decoy_id", how="left")
act = act.rename(columns={"a_mw": "mw", "a_clogp": "clogp"})
dec = dec.rename(columns={"d_mw": "mw", "d_clogp": "clogp"})


def tier(s4, s4t, t4, t4t):
    if s4 <= t4["p1"] and s4t <= t4t["p1"]:
        return "A"
    if s4 <= t4["p5"] and s4t <= t4t["p5"]:
        return "B"
    if s4 <= t4["med"] and s4t <= t4t["med"]:
        return "C"
    return "-"


def thresholds(train, decoys):
    return ({"p1": np.percentile(decoys.score_4d, 1),
             "p5": np.percentile(decoys.score_4d, 5),
             "med": np.median(train.score_4d)},
            {"p1": np.percentile(decoys.score_4t, 1),
             "p5": np.percentile(decoys.score_4t, 5),
             "med": np.median(train.score_4t)})


# ---- leave-one-scaffold-out ----
rows = []
for scaf, test in act.groupby("murcko_scaffold"):
    if len(test) == 0:
        continue
    train = act[act.murcko_scaffold != scaf]
    if len(train) < 30:
        continue
    t4, t4t = thresholds(train, dec)
    rec = test.apply(lambda r: tier(r.score_4d, r.score_4t, t4, t4t), axis=1)
    rows.append({"scaffold": scaf, "n_test": len(test),
                 "A": int((rec == "A").sum()), "B": int((rec == "B").sum()),
                 "C": int((rec == "C").sum())})
v = pd.DataFrame(rows)
v["A_pct"] = 100 * v.A / v.n_test
v["B_pct"] = 100 * v.B / v.n_test
v["C_pct"] = 100 * v.C / v.n_test
v.to_csv(os.path.join(OUT, "tier_validation_leave_scaffold.csv"), index=False)
print("leave-scaffold-out actives", int(v.n_test.sum()),
      "A/B/C recovery %", round(v.A.sum() / v.n_test.sum() * 100, 1),
      round(v.B.sum() / v.n_test.sum() * 100, 1),
      round(v.C.sum() / v.n_test.sum() * 100, 1))

# ---- decoy false tier rates with full-data thresholds ----
t4, t4t = thresholds(act, dec)
drec = dec.apply(lambda r: tier(r.score_4d, r.score_4t, t4, t4t), axis=1)
print("decoy A/B/C false %", round((drec == "A").mean() * 100, 2),
      round((drec == "B").mean() * 100, 2), round((drec == "C").mean() * 100, 2))

# ---- MW / LogP stratified enrichment ----
def auc_ef(act_s, dec_s):
    y = np.r_[np.ones(len(act_s)), np.zeros(len(dec_s))]
    s = np.r_[act_s, dec_s]
    order = np.argsort(s)
    y = y[order]
    tp = np.cumsum(y); fp = np.cumsum(1 - y)
    tpr = tp / tp[-1]; fpr = fp / fp[-1]
    auc = np.trapezoid(tpr, fpr)
    k = max(1, int(round(0.05 * len(y))))
    ef5 = (y[:k].sum() / k) / (len(act_s) / len(y))
    return auc, ef5

mw_bins = [0, 300, 400, 500, 1e9]
lp_bins = [-10, 1, 2.5, 4, 10]
rows = []
for rec, a_s, d_s in [("4DRE", "score_4d", "score_4t")]:
    pass
for label, col in [("4DRE", "score_4d"), ("4DTI", "score_4t")]:
    for name, bins in [("MW", mw_bins), ("LogP", lp_bins)]:
        for i in range(len(bins) - 1):
            aa = act[(act.mw >= bins[i]) & (act.mw < bins[i + 1])] if name == "MW" else act[(act.clogp >= bins[i]) & (act.clogp < bins[i + 1])]
            dd = dec[(dec.mw >= bins[i]) & (dec.mw < bins[i + 1])] if name == "MW" else dec[(dec.clogp >= bins[i]) & (dec.clogp < bins[i + 1])]
            if len(aa) < 10 or len(dd) < 10:
                continue
            auc, ef = auc_ef(aa[col].values, dd[col].values)
            rows.append({"receptor": label, "stratifier": name,
                         "bin": f"{bins[i]}-{bins[i+1]}", "n_act": len(aa), "n_dec": len(dd),
                         "auc": round(float(auc), 3), "ef5": round(float(ef), 2)})
pd.DataFrame(rows).to_csv(os.path.join(OUT, "property_stratified_enrichment.csv"), index=False)
print(pd.DataFrame(rows).to_string(index=False))

# ---- ligand efficiency of approved drugs ----
h = pd.read_csv(os.path.join(RES, "approved_hits.csv"))
def heavy(s):
    m = Chem.MolFromSmiles(s) if isinstance(s, str) else None
    return m.GetNumHeavyAtoms() if m else np.nan
h["heavy_atoms"] = h["smiles_std"].map(heavy)
h["le"] = -h["score_wt"] / h["heavy_atoms"]
cols = ["pref_name", "tier", "score_wt", "score_mu", "mw", "alogp", "heavy_atoms", "le"]
h[cols].sort_values("le", ascending=False).to_csv(os.path.join(OUT, "approved_ligand_efficiency.csv"), index=False)
print("top LE:")
print(h[cols].dropna().sort_values("le", ascending=False).head(10).to_string(index=False))
print("A-tier:")
print(h[h.tier == "A"][cols].to_string(index=False))
print("Spearman score_wt vs MW", spearmanr(h.score_wt, h.mw).statistic)
