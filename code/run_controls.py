# -*- coding: utf-8 -*-
"""G3 收口对照实验（两个脚本合一）：
  ex16  —— 前 30 名（按 4DRE 最优打分）用 exhaustiveness=16 复跑，算 Kendall τ 与分数差
  noise —— 随机 50 配体 × 3 个 seed（42/43/44）× 2 受体，量单构象对接的重复噪声

用法：python run_controls.py ex16
      python run_controls.py noise
输出：05-analysis/results/control_ex16.json / control_noise.json
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
CTRL = os.path.join(SCR, "control")
BOX = 22.0
VINA = "vina"


def center(rec):
    txt = open(os.path.join(REC, f"{rec}_vina_box.txt"), encoding="utf-8").read()
    kv = dict(re.findall(r"(center_[xyz])\s*=\s*(-?[\d.]+)", txt))
    return kv["center_x"], kv["center_y"], kv["center_z"]


def run_vina(job):
    rec, lig_id, seed, exh, out_dir = job
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, f"{rec}_{lig_id}_s{seed}_e{exh}.pdbqt")
    lig = os.path.join(SCR, "ligands", f"{lig_id}.pdbqt")
    t0 = time.time()
    cmd = [VINA, "--receptor", os.path.join(SCR, "receptor", f"{rec}.pdbqt"),
           "--ligand", lig,
           "--center_x", center(rec)[0], "--center_y", center(rec)[1], "--center_z", center(rec)[2],
           "--size_x", str(BOX), "--size_y", str(BOX), "--size_z", str(BOX),
           "--exhaustiveness", str(exh), "--num_modes", "1", "--seed", str(seed),
           "--cpu", "4", "--out", out]
    r = subprocess.run(cmd, capture_output=True, encoding="utf-8", errors="replace")
    score = None
    if os.path.exists(out):
        for line in open(out, encoding="utf-8", errors="replace"):
            if line.startswith("REMARK VINA RESULT"):
                score = float(line.split()[3])
                break
    return dict(receptor=rec, lig_id=lig_id, seed=seed, exhaustiveness=exh,
                score=score, seconds=round(time.time() - t0, 1),
                status="ok" if score is not None else f"fail_rc{r.returncode}")


def dock_many(jobs, workers=6):
    rows = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(run_vina, j) for j in jobs]
        for n, f in enumerate(as_completed(futs), 1):
            rows.append(f.result())
            if n % 25 == 0 or n == len(futs):
                print(f"  {n}/{len(futs)}", flush=True)
    return pd.DataFrame(rows)


def ex16():
    D = pd.read_csv(os.path.join(DOCK, "dock_4DRE_all.csv"))
    top = D.nsmallest(30, "score")[["lig_id", "score"]].rename(columns={"score": "score_ex8"})
    jobs = [(rec, lid, 42, 16, os.path.join(CTRL, "ex16")) for rec in ("4DRE", "4DTI") for lid in top.lig_id]
    df = dock_many(jobs)
    df.to_csv(os.path.join(RES, "control_ex16_raw.csv"), index=False, encoding="utf-8-sig")
    out = {"n_ligands": int(len(top)), "exhaustiveness": [8, 16], "seed": 42, "by_receptor": {}}
    for rec in ("4DRE", "4DTI"):
        sub = df[df.receptor == rec].set_index("lig_id")
        base = pd.read_csv(os.path.join(DOCK, f"dock_{rec}_all.csv")).set_index("lig_id")["score"]
        common = [l for l in top.lig_id if l in sub.index and sub.loc[l, "status"] == "ok"]
        a = base.loc[common].values
        b = sub.loc[common, "score"].values
        tau = float(pd.Series(a).corr(pd.Series(b), method="kendall"))
        out["by_receptor"][rec] = dict(
            n=int(len(common)), kendall_tau=round(tau, 4),
            delta_mean=round(float(np.mean(b - a)), 3), delta_max_abs=round(float(np.max(np.abs(b - a))), 3),
            identical=int(np.sum(np.isclose(a, b, atol=1e-6))))
    json.dump(out, open(os.path.join(RES, "control_ex16.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print(json.dumps(out, ensure_ascii=False, indent=2))


def noise():
    rng = np.random.default_rng(42)
    ligs = list(pd.read_csv(os.path.join(DOCK, "dock_4DRE_all.csv")).lig_id)
    pick = list(rng.choice(ligs, size=50, replace=False))
    jobs = [(rec, lid, seed, 8, os.path.join(CTRL, "noise"))
            for rec in ("4DRE", "4DTI") for lid in pick for seed in (42, 43, 44)]
    df = dock_many(jobs)
    df.to_csv(os.path.join(RES, "control_noise_raw.csv"), index=False, encoding="utf-8-sig")
    out = {"n_ligands": len(pick), "seeds": [42, 43, 44], "exhaustiveness": 8, "by_receptor": {}}
    for rec in ("4DRE", "4DTI"):
        sub = df[(df.receptor == rec) & (df.status == "ok")]
        wide = sub.pivot_table(index="lig_id", columns="seed", values="score")
        sd = wide.std(axis=1, ddof=1)
        rng_ = wide.max(axis=1) - wide.min(axis=1)
        base = pd.read_csv(os.path.join(DOCK, f"dock_{rec}_all.csv")).set_index("lig_id")["score"]
        det = float(np.mean(np.isclose(wide[42].values, base.loc[wide.index].values, atol=1e-6)))
        out["by_receptor"][rec] = dict(
            n=int(len(wide)), sd_median=round(float(sd.median()), 3), sd_p90=round(float(sd.quantile(.9)), 3),
            range_median=round(float(rng_.median()), 3), range_max=round(float(rng_.max()), 3),
            seed42_reproduces_primary=round(det, 3))
    json.dump(out, open(os.path.join(RES, "control_noise.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print(json.dumps(out, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else ""
    if which == "ex16":
        ex16()
    elif which == "noise":
        noise()
    else:
        print("usage: python run_controls.py ex16|noise")
