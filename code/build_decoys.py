"""Build a property-matched decoy set for the InhA retrospective validation (v2).

Why v2: v1 computed descriptors and fingerprints on RAW pool SMILES while the actives had
already been standardised (FragmentParent + Uncharge). Four decoys were salts and failed
ligand preparation; eleven carried formal charges, so descriptors described a molecule that
would never be docked. v2 standardises everything identically, uses stable InChIKey-based IDs,
and tops the set up to exactly 1:1.

Pre-registered spec (G2, 03-methods/methods.md §2), unchanged:
  pool        : TDC tox21 (Harvard Dataverse id 4259612)
  match       : MW +-10%, cLogP +-1.0, HBD +-1, HBA +-1, rotB +-1, formal charge identical
  similarity  : ECFP4 (r=2, 2048 bit) Tanimoto < 0.3 against EVERY active
  ratio       : 1:1, unique decoy per active, sampling seed 42
  fallback    : one relaxation pass, flagged per decoy
"""
import gzip
import json
import os
import random
import urllib.request

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Crippen, Descriptors, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")

PROJ = r"D:\workspace\research\protein-drug"
RAW = os.path.join(PROJ, "04-data", "raw", "tdc")
DERIVED = os.path.join(PROJ, "04-data", "derived")
RESULTS = os.path.join(PROJ, "05-analysis", "results")
TOX21_PATH = os.path.join(RAW, "tox21.tab")
ALLOWED = {6, 7, 8, 9, 15, 16, 17, 35, 53, 1}
SEED = 42
UN = rdMolStandardize.Uncharger()
os.makedirs(RAW, exist_ok=True)
os.makedirs(DERIVED, exist_ok=True)
os.makedirs(RESULTS, exist_ok=True)


def standardise(m):
    m = rdMolStandardize.FragmentParent(m)
    m = UN.uncharge(m)
    Chem.SanitizeMol(m)          # RingInfo must be initialised before descriptor calls
    return m


def usable(smiles):
    m = Chem.MolFromSmiles(str(smiles))
    if m is None:
        return None
    try:
        m = standardise(m)
    except Exception:
        return None
    if m is None or m.GetNumAtoms() < 3:
        return None
    if not all(a.GetAtomicNum() in ALLOWED for a in m.GetAtoms()):
        return None
    return m


def props(m):
    return dict(mw=Descriptors.MolWt(m), clogp=Crippen.MolLogP(m),
                hbd=rdMolDescriptors.CalcNumHBD(m), hba=rdMolDescriptors.CalcNumHBA(m),
                rotb=rdMolDescriptors.CalcNumRotatableBonds(m), charge=Chem.GetFormalCharge(m))


def download_tox21():
    if os.path.exists(TOX21_PATH) and os.path.getsize(TOX21_PATH) > 10000:
        return
    url = "https://dataverse.harvard.edu/api/access/datafile/4259612"
    with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "fetch/1.0"}), timeout=120) as f:
        blob = f.read()
    if blob[:2] == b"\x1f\x8b":
        blob = gzip.decompress(blob)
    with open(TOX21_PATH, "wb") as fh:
        fh.write(blob)


download_tox21()

# ---------------- actives ----------------
act = pd.read_csv(os.path.join(DERIVED, "inha_ic50_clean.csv"))
act_mols, act_p, act_smi, act_ids, act_pic50 = [], [], [], [], []
for row in act.itertuples():
    m = usable(row.smiles_std)
    if m is None:
        continue
    act_mols.append(m)
    act_p.append(props(m))
    act_smi.append(Chem.MolToSmiles(m))
    act_ids.append(f"A_{Chem.MolToInchiKey(m)[:14]}")
    act_pic50.append(float(row.pIC50_median))
act_prop = pd.DataFrame(act_p)
act_fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) for m in act_mols]
print(f"actives usable: {len(act_mols)}")

# ---------------- pool ----------------
pool = pd.read_csv(TOX21_PATH, sep="\t")
if "X" in pool.columns:
    smiles_col = "X"
else:
    smiles_col = next(c for c in pool.columns
                      if sum(Chem.MolFromSmiles(str(v)) is not None
                             for v in pool[c].dropna().head(60)) >= 48)

cand_mols, cand_smi, cand_p = [], [], []
seen = set()
for s in pool[smiles_col].astype(str):
    m = usable(s)
    if m is None:
        continue
    smi = Chem.MolToSmiles(m)
    if smi in seen:
        continue
    seen.add(smi)
    cand_mols.append(m)
    cand_smi.append(smi)
    cand_p.append(props(m))
cand_prop = pd.DataFrame(cand_p)
cand_fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) for m in cand_mols]
print(f"pool usable (standardised, deduped): {len(cand_mols)}")

max_t = np.array([max(DataStructs.BulkTanimotoSimilarity(fp, act_fps)) for fp in cand_fps])
print(f"candidates with max Tanimoto < 0.3: {int((max_t < 0.3).sum())}")

# ---------------- matching ----------------
rng = random.Random(SEED)
used = np.zeros(len(cand_mols), dtype=bool)
TIERS = [dict(name="strict", mw=0.10, clogp=1.0, hb=1, rotb=1),
         dict(name="relaxed", mw=0.20, clogp=1.5, hb=2, rotb=2)]


def candidates_for(i, tier):
    ap = act_prop.iloc[i]
    mask = ((~used) & (max_t < 0.3)
            & (np.abs(cand_prop["mw"].values - ap["mw"]) <= tier["mw"] * max(ap["mw"], 1.0))
            & (np.abs(cand_prop["clogp"].values - ap["clogp"]) <= tier["clogp"])
            & (np.abs(cand_prop["hbd"].values - ap["hbd"]) <= tier["hb"])
            & (np.abs(cand_prop["hba"].values - ap["hba"]) <= tier["hb"])
            & (np.abs(cand_prop["rotb"].values - ap["rotb"]) <= tier["rotb"])
            & (cand_prop["charge"].values == ap["charge"]))
    return np.where(mask)[0]


matched = []          # list of dict(active_index, decoy_index, tier)
assigned_actives = set()

for tier in TIERS:
    for i in range(len(act_mols)):
        if i in assigned_actives:
            continue
        idxs = candidates_for(i, tier)
        if len(idxs) == 0:
            continue
        j = int(rng.choice(list(idxs)))
        used[j] = True
        assigned_actives.add(i)
        matched.append(dict(active_index=i, decoy_index=j, tier=tier["name"]))
    print(f"after {tier['name']}: {len(matched)} matched  (actives covered {len(assigned_actives)})")

# top-up to exact 1:1: extra decoys matched to already-covered actives (documented, flagged)
topup, guard = 0, 0
while len(matched) < len(act_mols) and guard < 200:
    guard += 1
    i = rng.choice(sorted(assigned_actives)) if assigned_actives else rng.randrange(len(act_mols))
    for tier in reversed(TIERS):
        idxs = candidates_for(i, tier)
        if len(idxs):
            j = int(rng.choice(list(idxs)))
            used[j] = True
            matched.append(dict(active_index=i, decoy_index=j, tier="topup"))
            topup += 1
            break
print(f"top-up decoys: {topup}  -> total {len(matched)}")

rows = []
seen_ids = {}

def unique_id(mol):
    """14-char InChIKey blocks ignore stereochemistry, so stereoisomers collide.
    Keep the first occurrence plain and suffix later ones deterministically."""
    base = f"D_{Chem.MolToInchiKey(mol)[:14]}"
    n = seen_ids.get(base, 0) + 1
    seen_ids[base] = n
    return base if n == 1 else f"{base}_{chr(96 + n)}"   # _b, _c, ...


for r in matched:
    i, j = r["active_index"], r["decoy_index"]
    rows.append(dict(
        decoy_id=unique_id(cand_mols[j]),
        active_id=act_ids[i], active_pIC50=act_pic50[i],
        active_smiles=act_smi[i], decoy_smiles=cand_smi[j],
        maxT=float(max_t[j]), tier=r["tier"],
        a_mw=act_prop.at[i, "mw"], a_clogp=act_prop.at[i, "clogp"], a_hbd=act_prop.at[i, "hbd"],
        a_hba=act_prop.at[i, "hba"], a_rotb=act_prop.at[i, "rotb"], a_charge=act_prop.at[i, "charge"],
        d_mw=cand_prop.at[j, "mw"], d_clogp=cand_prop.at[j, "clogp"], d_hbd=cand_prop.at[j, "hbd"],
        d_hba=cand_prop.at[j, "hba"], d_rotb=cand_prop.at[j, "rotb"], d_charge=cand_prop.at[j, "charge"],
    ))
dec = pd.DataFrame(rows)
dec.to_csv(os.path.join(DERIVED, "decoys_tox21.csv"), index=False, encoding="utf-8-sig")

report = dict(
    pool_source="Harvard Dataverse datafile 4259612 (TDC tox21)",
    standardisation="FragmentParent + Uncharge applied to pool AND actives before descriptors",
    pool_rows=int(len(pool)), pool_usable=int(len(cand_mols)),
    actives=int(len(act_mols)), decoys=int(len(dec)), ratio=f"{len(dec)}:{len(act_mols)}",
    tiers={t["name"]: int((dec["tier"] == t["name"]).sum()) for t in TIERS} | {"topup": int((dec["tier"] == "topup").sum())},
    maxT_median=round(float(dec["maxT"].median()), 3),
    maxT_max=round(float(dec["maxT"].max()), 3),
    delta_mean={k: round(float((dec[f"d_{k}"] - dec[f"a_{k}"]).abs().mean()), 3)
                for k in ("mw", "clogp", "hbd", "hba", "rotb")},
    decoys_multi_fragment=int(sum(len(Chem.GetMolFrags(Chem.MolFromSmiles(s))) > 1 for s in dec["decoy_smiles"])),
    decoys_with_charge=int(sum(any(a.GetFormalCharge() for a in Chem.MolFromSmiles(s).GetAtoms())
                               for s in dec["decoy_smiles"])),
    seed=SEED,
)
with open(os.path.join(RESULTS, "decoy_report.json"), "w", encoding="utf-8") as fh:
    json.dump(report, fh, ensure_ascii=False, indent=2)
print(json.dumps(report, ensure_ascii=False, indent=2))
