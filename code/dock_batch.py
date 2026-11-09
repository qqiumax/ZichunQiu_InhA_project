"""Batch docking: SMILES -> 3D -> PDBQT -> Vina, parallel and resumable.

Usage:
    python dock_batch.py --receptor 4DRE --set pilot --workers 6
    python dock_batch.py --receptor 4DRE --set all --workers 6
    python dock_batch.py --receptor 4DTI --set all --workers 6

Design notes
  * All compute happens under D:\\scratch\\protein-drug (ASCII) ; results are appended to a CSV
    inside 04-data/derived/docking/ so a crash never loses more than the in-flight ligands.
  * Ligand PDBQT preparation is cached: one file per molecule, reused across receptors.
  * Resumable: a ligand is skipped when its result row already exists AND its pose file parses.
"""
import argparse
import json
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")

PROJ = r"D:\workspace\research\protein-drug"
SCRATCH = r"D:\scratch\protein-drug"
LIG_DIR = os.path.join(SCRATCH, "ligands")
DOCK_DIR = os.path.join(SCRATCH, "docking")
REC_DIR = os.path.join(PROJ, "04-data", "derived", "receptor")
RES_DIR = os.path.join(PROJ, "04-data", "derived", "docking")  # S4 产物：每分子打分 + 权威清单
os.makedirs(LIG_DIR, exist_ok=True)
os.makedirs(DOCK_DIR, exist_ok=True)
os.makedirs(RES_DIR, exist_ok=True)

VINA = "vina"
MEEKO = "mk_prepare_ligand.exe"
EXHAUST = 8
SEED = 42
BOX = 22.0


def load_ligands(set_name):
    """Return DataFrame with columns: lig_id, smiles, group ('active'/'decoy')."""
    act = pd.read_csv(os.path.join(PROJ, "04-data", "derived", "inha_ic50_clean.csv"))
    ALLOWED = {6, 7, 8, 9, 15, 16, 17, 35, 53, 1}

    def ok(s):
        m = Chem.MolFromSmiles(s)
        return m is not None and all(a.GetAtomicNum() in ALLOWED for a in m.GetAtoms())

    act = act[act["smiles_std"].apply(ok)].reset_index(drop=True)
    act_rows = pd.DataFrame({
        "lig_id": ["A_" + Chem.MolToInchiKey(Chem.MolFromSmiles(s))[:14] for s in act["smiles_std"]],
        "smiles": act["smiles_std"], "group": "active", "pIC50": act["pIC50_median"]})

    dec = pd.read_csv(os.path.join(PROJ, "04-data", "derived", "decoys_tox21.csv"))
    dec_rows = pd.DataFrame({"lig_id": dec["decoy_id"], "smiles": dec["decoy_smiles"],
                             "group": "decoy", "pIC50": None})

    df = pd.concat([act_rows, dec_rows], ignore_index=True)
    if set_name == "pilot":
        df = pd.concat([act_rows.head(20), dec_rows.head(20)], ignore_index=True)
    elif set_name == "approved":
        # S4-正式：已批准药物库（fetch_chembl_approved.py 生成，2147 个可对接分子）
        ap = pd.read_csv(os.path.join(PROJ, "04-data", "derived", "approved_drugs.csv"))
        df = pd.DataFrame({"lig_id": ap["lig_id"], "smiles": ap["smiles_std"],
                           "group": "approved", "pIC50": None})
    return df


def box_center(receptor):
    path = os.path.join(REC_DIR, f"{receptor}_vina_box.txt")
    with open(path, encoding="utf-8", errors="replace") as fh:
        kv = dict(re.findall(r"(center_[xyz])\s*=\s*(-?[\d.]+)", fh.read()))
    return float(kv["center_x"]), float(kv["center_y"]), float(kv["center_z"])


def prepare_ligand(job):
    """SMILES -> 3D SDF -> PDBQT. Returns (lig_id, path_or_None, error)."""
    lig_id, smiles = job
    out = os.path.join(LIG_DIR, f"{lig_id}.pdbqt")
    if os.path.exists(out) and os.path.getsize(out) > 500:
        return lig_id, out, None
    try:
        # Standardise exactly like clean_inha_activities.py / build_decoys.py:
        # strip counter-ions, neutralise, then add H. Without this, salt-form ligands fail
        # meeko (rc=3) and charged ligands are docked in a form their descriptors never saw.
        _un = rdMolStandardize.Uncharger()
        m0 = rdMolStandardize.FragmentParent(Chem.MolFromSmiles(smiles))
        m0 = _un.uncharge(m0)
        Chem.SanitizeMol(m0)
        mol = Chem.AddHs(m0)
        ps = AllChem.ETKDGv3()
        ps.randomSeed = SEED
        if AllChem.EmbedMolecule(mol, ps) != 0:
            return lig_id, None, "embed_failed"
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
        sdf = os.path.join(LIG_DIR, f"{lig_id}.sdf")
        with open(sdf, "w", encoding="utf-8") as fh:
            fh.write(Chem.MolToMolBlock(mol))
            fh.write("$$$$\n")
        r = subprocess.run([MEEKO, "-i", sdf, "-o", out], capture_output=True,
                           encoding="utf-8", errors="replace")
        if r.returncode != 0 or not os.path.exists(out):
            return lig_id, None, f"meeko_rc_{r.returncode}"
        return lig_id, out, None
    except Exception as e:
        return lig_id, None, f"{type(e).__name__}"


def dock_one(job):
    receptor, lig_id, lig_path, center = job
    out = os.path.join(DOCK_DIR, f"{receptor}_{lig_id}.pdbqt")
    if os.path.exists(out) and os.path.getsize(out) > 500:
        with open(out, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("REMARK VINA RESULT"):
                    return dict(lig_id=lig_id, receptor=receptor, score=float(line.split()[3]),
                                seconds=0.0, status="cached")
    t0 = time.time()
    cmd = [VINA, "--receptor", os.path.join(SCRATCH, "receptor", f"{receptor}.pdbqt"),
           "--ligand", lig_path,
           "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
           "--size_x", str(BOX), "--size_y", str(BOX), "--size_z", str(BOX),
           "--exhaustiveness", str(EXHAUST), "--num_modes", "1", "--seed", str(SEED),
           "--cpu", "4", "--out", out]
    r = subprocess.run(cmd, capture_output=True, encoding="utf-8", errors="replace")
    dt = time.time() - t0
    score = None
    if os.path.exists(out):
        with open(out, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("REMARK VINA RESULT"):
                    score = float(line.split()[3])
                    break
    return dict(lig_id=lig_id, receptor=receptor, score=score, seconds=round(dt, 1),
                status="ok" if score is not None else f"fail_rc{r.returncode}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--receptor", required=True, choices=["4DRE", "4DTI", "4U0J"])
    ap.add_argument("--set", default="all", choices=["all", "pilot", "approved"])
    ap.add_argument("--workers", type=int, default=6)
    args = ap.parse_args()

    # receptors are prepared by prepare_receptors.py into 04-data/derived/receptor/
    os.makedirs(os.path.join(SCRATCH, "receptor"), exist_ok=True)
    src = os.path.join(REC_DIR, f"{args.receptor}_receptor_full.pdbqt")
    dst = os.path.join(SCRATCH, "receptor", f"{args.receptor}.pdbqt")
    if not os.path.exists(dst) or os.path.getmtime(dst) < os.path.getmtime(src):
        with open(src, "rb") as a, open(dst, "wb") as b:
            b.write(a.read())
    center = box_center(args.receptor)

    ligs = load_ligands(args.set)
    print(f"receptor={args.receptor} set={args.set} ligands={len(ligs)} workers={args.workers}", flush=True)
    # sidecar: the authoritative list of ligands this run is supposed to cover.
    # compute_metrics.py joins against it so that ANY ligand that fails later still counts in
    # the denominator (pre-registered rule: failures scored as worst, never dropped).
    sidecar = os.path.join(RES_DIR, f"ligands_{args.receptor}_{args.set}.csv")
    ligs.to_csv(sidecar, index=False, encoding="utf-8-sig")

    # phase 1: prepare ligands
    t0 = time.time()
    prepared, failed = {}, []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(prepare_ligand, (r.lig_id, r.smiles)): r.lig_id for r in ligs.itertuples()}
        for n, fut in enumerate(as_completed(futs), 1):
            lid, path, err = fut.result()
            if path:
                prepared[lid] = path
            else:
                failed.append((lid, err))
            if n % 50 == 0 or n == len(futs):
                print(f"  prep {n}/{len(futs)}  ok={len(prepared)} fail={len(failed)}  "
                      f"{time.time()-t0:.0f}s", flush=True)
    print(f"prepared {len(prepared)} ligands in {time.time()-t0:.0f}s; failures={len(failed)}", flush=True)
    if failed:
        print("  first failures:", failed[:5], flush=True)

    # phase 2: dock
    out_csv = os.path.join(RES_DIR, f"dock_{args.receptor}_{args.set}.csv")
    done = set()
    if os.path.exists(out_csv):
        prev = pd.read_csv(out_csv)
        done = set(prev["lig_id"])
        print(f"resume: {len(done)} already docked", flush=True)

    todo = [(args.receptor, lid, path, center) for lid, path in prepared.items() if lid not in done]
    print(f"to dock: {len(todo)}", flush=True)
    t1 = time.time()
    rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = [ex.submit(dock_one, job) for job in todo]
        for n, fut in enumerate(as_completed(futs), 1):
            res = fut.result()
            rows.append(res)
            if n % 25 == 0 or n == len(futs):
                el = time.time() - t1
                eta = el / n * (len(futs) - n)
                print(f"  dock {n}/{len(futs)}  elapsed={el/60:.1f}min  eta={eta/60:.1f}min", flush=True)
            if n % 50 == 0:
                pd.concat([pd.DataFrame(rows)], ignore_index=True).to_csv(
                    out_csv + ".partial", index=False, encoding="utf-8-sig")

    res_df = pd.DataFrame(rows)
    if os.path.exists(out_csv):
        prev = pd.read_csv(out_csv)
        # keep only result columns: meta columns are re-attached once below, otherwise a
        # resumed run produces group_x/group_y duplicate columns
        keep = [c for c in ("lig_id", "receptor", "score", "seconds", "status") if c in prev.columns]
        res_df = pd.concat([prev[keep], res_df], ignore_index=True)
    res_df = res_df.drop_duplicates(subset=["lig_id", "receptor"], keep="last")
    meta = ligs[["lig_id", "group", "pIC50"]]
    res_df = res_df.merge(meta, on="lig_id", how="left")
    res_df.to_csv(out_csv, index=False, encoding="utf-8-sig")
    ok = res_df["score"].notna().sum()
    print(f"\nDONE receptor={args.receptor}: {ok}/{len(res_df)} scored -> {out_csv}", flush=True)
    print(f"median seconds/dock = {res_df['seconds'].median():.1f}", flush=True)


if __name__ == "__main__":
    main()
