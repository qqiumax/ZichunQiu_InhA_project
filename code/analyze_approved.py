# -*- coding: utf-8 -*-
"""已批准药库筛选结果分析（判据见 03-methods/methods.md §10，2026-11-01 开跑前预注册）。

输入：04-data/derived/docking/dock_4DRE_approved.csv、dock_4DTI_approved.csv、derived/approved_drugs.csv
输出：05-analysis/results/approved_hits.csv、approved_summary.json
"""
import json
import os

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors
from scipy import stats

RDLogger.DisableLog("rdApp.*")
PROJ = r"D:\workspace\research\protein-drug"
DOCK = os.path.join(PROJ, "04-data", "derived", "docking")
RES = os.path.join(PROJ, "05-analysis", "results")

# 预注册阈值（methods §10，来自回顾性验证的诱饵/活性物分布）
TH = {
    "A": {"4DRE": -9.422, "4DTI": -9.205},   # 诱饵第 1 百分位
    "B": {"4DRE": -8.918, "4DTI": -8.759},   # 诱饵第 5 百分位
    "C": {"4DRE": -8.163, "4DTI": -8.024},   # 活性物中位数
}


def load(rec):
    p = os.path.join(DOCK, f"dock_{rec}_approved.csv")
    d = pd.read_csv(p)
    d = d.drop_duplicates(subset="lig_id", keep="last")
    return d.set_index("lig_id")


def main():
    wt = load("4DRE")
    mu = load("4DTI")
    meta = pd.read_csv(os.path.join(PROJ, "04-data", "derived", "approved_drugs.csv")).set_index("lig_id")

    m = meta.join(wt[["score"]].rename(columns={"score": "score_wt"}), how="left")
    m = m.join(mu[["score"]].rename(columns={"score": "score_mu"}), how="left")
    m["n_scored"] = m[["score_wt", "score_mu"]].notna().sum(axis=1)
    m["score_best"] = m[["score_wt", "score_mu"]].min(axis=1)
    m["score_worst"] = m[["score_wt", "score_mu"]].max(axis=1)
    m["dds"] = m["score_mu"] - m["score_wt"]
    m["tier"] = ""
    # 从弱到强赋值：A ⊂ B ⊂ C，最后写 A 才能保留最高层
    for t, th in (("C", TH["C"]), ("B", TH["B"]), ("A", TH["A"])):
        cond = (m["score_wt"] <= th["4DRE"]) & (m["score_mu"] <= th["4DTI"])
        m.loc[cond, "tier"] = t
    m.loc[m["tier"] == "", "tier"] = "-"
    m = m.sort_values(["tier", "score_best"])

    out = os.path.join(RES, "approved_hits.csv")
    m.reset_index().to_csv(out, index=False, encoding="utf-8-sig")

    scored = m[m["n_scored"] == 2]
    rep = dict(
        n_library=int(len(m)), n_scored=int((m["n_scored"] == 2).sum()),
        thresholds=TH,
        tier_counts={t: int((m["tier"] == t).sum()) for t in ("A", "B", "C")},
        frac_better_than_decoy_p1={rec: round(float((scored[f"score_{'wt' if rec=='4DRE' else 'mu'}"] <= TH["A"][rec]).mean()), 4)
                                   for rec in ("4DRE", "4DTI")},
        score_dist={rec: dict(median=round(float(scored[col].median()), 3),
                              p10=round(float(scored[col].quantile(.1)), 3),
                              p90=round(float(scored[col].quantile(.9)), 3))
                    for rec, col in (("4DRE", "score_wt"), ("4DTI", "score_mu"))},
        withdrawn_in_A=int((m["tier"] == "A").sum() and m.loc[m["tier"] == "A", "withdrawn"].sum()),
    )
    # 性质偏倚：库内打分 vs MW/LogP
    mw = scored["mw"].values
    lp = np.array([Crippen.MolLogP(Chem.MolFromSmiles(s)) if Chem.MolFromSmiles(s) else np.nan
                   for s in scored["smiles_std"]])
    ok = ~np.isnan(lp)
    rep["property_bias"] = dict(
        score_mw_rho=round(float(stats.spearmanr(scored["score_best"][ok], mw[ok]).statistic), 3),
        score_logp_rho=round(float(stats.spearmanr(scored["score_best"][ok], lp[ok]).statistic), 3))
    rep["top50"] = m.head(50)[["chembl_id", "pref_name", "tier", "score_wt", "score_mu", "mw", "first_approval"]].to_dict("records")
    json.dump(rep, open(os.path.join(RES, "approved_summary.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=2)
    print(json.dumps({k: rep[k] for k in ("n_library", "n_scored", "tier_counts", "frac_better_than_decoy_p1", "property_bias")},
                     ensure_ascii=False, indent=2))
    print("written:", out)


if __name__ == "__main__":
    main()
