# -*- coding: utf-8 -*-
"""A+B 候选的二次确认（methods §11.6）：
  * 两受体 ex16（seed 42）
  * 两受体 seed 43 / 44（ex8）
  * 第三受体 4U0J（ex8, seed 42）
稳定性判据（与 §11.6 一致）：ex16 与初筛差 ≤1.0；seed 差 ≤1.0；4U0J ≤ 其诱饵 5% 分位。

用法：python run_confirm.py [workers]
输出：results/confirm_per_drug.csv、results/confirm_summary.json
"""
import json
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd

PROJ = r"D:\workspace\research\protein-drug"
SCR = r"D:\scratch\protein-drug"
DOCK = os.path.join(PROJ, "04-data", "derived", "docking")
REC = os.path.join(PROJ, "04-data", "derived", "receptor")
RES = os.path.join(PROJ, "05-analysis", "results")
OUTDIR = os.path.join(SCR, "confirm")
BOX = 22.0
VINA = "vina"


def center(rec):
    txt = open(os.path.join(REC, f"{rec}_vina_box.txt"), encoding="utf-8").read()
    kv = dict(re.findall(r"(center_[xyz])\s*=\s*(-?[\d.]+)", txt))
    return kv["center_x"], kv["center_y"], kv["center_z"]


def run_one(job):
    rec, lig, seed, exh = job
    os.makedirs(OUTDIR, exist_ok=True)
    out = os.path.join(OUTDIR, f"{rec}_{lig}_s{seed}_e{exh}.pdbqt")
    ligp = os.path.join(SCR, "ligands", f"{lig}.pdbqt")
    cx, cy, cz = center(rec)
    cmd = [VINA, "--receptor", os.path.join(SCR, "receptor", f"{rec}.pdbqt"), "--ligand", ligp,
           "--center_x", cx, "--center_y", cy, "--center_z", cz,
           "--size_x", str(BOX), "--size_y", str(BOX), "--size_z", str(BOX),
           "--exhaustiveness", str(exh), "--num_modes", "1", "--seed", str(seed),
           "--cpu", "4", "--out", out]
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, encoding="utf-8", errors="replace")
    score = None
    if os.path.exists(out):
        for line in open(out, encoding="utf-8", errors="replace"):
            if line.startswith("REMARK VINA RESULT"):
                score = float(line.split()[3])
                break
    return dict(receptor=rec, lig_id=lig, seed=seed, exh=exh, score=score,
                seconds=round(time.time() - t0, 1),
                status="ok" if score is not None else f"fail_rc{r.returncode}")


def main():
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    hits = pd.read_csv(os.path.join(RES, "approved_hits.csv"))
    cand = hits[hits.tier.isin(["A", "B"])].copy()
    ligs = list(cand.lig_id)
    print(f"candidates: {len(ligs)} (A={int((cand.tier=='A').sum())}, B={int((cand.tier=='B').sum())})", flush=True)

    jobs = []
    for lig in ligs:
        for rec in ("4DRE", "4DTI"):
            jobs.append((rec, lig, 42, 16))
            jobs.append((rec, lig, 43, 8))
            jobs.append((rec, lig, 44, 8))
        jobs.append(("4U0J", lig, 42, 8))
    print(f"jobs: {len(jobs)}", flush=True)

    rows = []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(run_one, j) for j in jobs]
        for n, f in enumerate(as_completed(futs), 1):
            rows.append(f.result())
            if n % 20 == 0 or n == len(futs):
                el = time.time() - t0
                print(f"  {n}/{len(futs)}  {el/60:.1f} min  eta {(el/n*(len(futs)-n))/60:.1f} min", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(RES, "confirm_raw.csv"), index=False, encoding="utf-8-sig")

    # 4U0J 诱饵阈值
    d4 = pd.read_csv(os.path.join(DOCK, "dock_4U0J_all.csv"))
    dec = d4[d4.group == "decoy"].score.dropna()
    p1, p5 = float(dec.quantile(.01)), float(dec.quantile(.05))

    wide = df.pivot_table(index=["lig_id", "receptor"], columns=["seed", "exh"], values="score")
    summary = []
    for _, r in cand.iterrows():
        lig = r.lig_id
        out = dict(lig_id=lig, pref_name=r.pref_name, chembl_id=r.chembl_id, tier=r.tier)
        base = {"4DRE": r.score_wt, "4DTI": r.score_mu}
        for rec in ("4DRE", "4DTI"):
            key16 = (42, 16)
            try:
                v16 = wide.loc[(lig, rec), key16]
            except Exception:
                v16 = np.nan
            d16 = abs(float(v16) - float(base[rec])) if pd.notna(v16) else np.nan
            seeds = [wide.loc[(lig, rec), (s, 8)] for s in (43, 44) if (lig, rec) in wide.index]
            dseed = max([abs(float(v) - float(base[rec])) for v in seeds if pd.notna(v)]) if seeds else np.nan
            out[f"{rec}_ex16"] = round(float(v16), 3) if pd.notna(v16) else None
            out[f"{rec}_d16"] = round(d16, 3) if pd.notna(d16) else None
            out[f"{rec}_dseed"] = round(dseed, 3) if pd.notna(dseed) else None
        try:
            out["4U0J"] = round(float(wide.loc[(lig, "4U0J"), (42, 8)]), 3)
        except Exception:
            out["4U0J"] = None
        out["4U0J_ok"] = bool(out["4U0J"] is not None and out["4U0J"] <= p5)

        def ok(v):
            return v is not None and v <= 1.0

        out["stable"] = bool(ok(out.get("4DRE_d16")) and ok(out.get("4DTI_d16"))
                             and ok(out.get("4DRE_dseed")) and ok(out.get("4DTI_dseed"))
                             and out["4U0J_ok"])
        summary.append(out)
    S = pd.DataFrame(summary)
    S.to_csv(os.path.join(RES, "confirm_per_drug.csv"), index=False, encoding="utf-8-sig")
    rep = dict(n_candidates=len(cand), workers=workers, jobs=len(jobs),
               four_d_u0j_thresholds=dict(p1=round(p1, 3), p5=round(p5, 3)),
               stable_counts={t: int(((S.tier == t) & S.stable).sum()) for t in ("A", "B")},
               stable_A=S[S.tier == "A"].to_dict("records"),
               stable_B_n=int(((S.tier == "B") & S.stable).sum()))
    json.dump(rep, open(os.path.join(RES, "confirm_summary.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=2)
    print(json.dumps({k: rep[k] for k in ("n_candidates", "stable_counts", "four_d_u0j_thresholds")},
                     ensure_ascii=False, indent=2))
    print("written", os.path.join(RES, "confirm_per_drug.csv"))


if __name__ == "__main__":
    main()
