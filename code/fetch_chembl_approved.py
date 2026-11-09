# -*- coding: utf-8 -*-
"""S4-正式 · 已批准药物库采集（ChEMBL max_phase=4）。

- 原始分页 JSON 原样落盘到 04-data/raw/chembl/approved/（永不覆盖）
- 解析 + 标准化（FragmentParent → Uncharge → Sanitize）+ 元素白名单 → 04-data/derived/approved_drugs.csv
- 记录每个丢弃理由的条数，写 05-analysis/results/approved_library_report.json
- 不在这里抽子集；抽样规则见 03-methods/methods.md §6-A（预注册）

用法：python fetch_chembl_approved.py
"""
import hashlib
import json
import os
import time
import urllib.request

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")

PROJ = r"D:\workspace\research\protein-drug"
RAW = os.path.join(PROJ, "04-data", "raw", "chembl", "approved")
DER = os.path.join(PROJ, "04-data", "derived")
RES = os.path.join(PROJ, "05-analysis", "results")
BASE = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
PARAMS = "max_phase=4&limit=1000"
ALLOWED = {6, 7, 8, 9, 15, 16, 17, 35, 53, 1}


def fetch(offset, retries=6):
    url = f"{BASE}?{PARAMS}&offset={offset}"
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "research/1.0"})
            with urllib.request.urlopen(req, timeout=60) as f:
                return json.loads(f.read().decode("utf-8"))
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(min(30, 2 ** i))
    raise RuntimeError(f"fetch failed offset={offset}: {last}")


def standardize(smiles):
    if not smiles:
        return None, "no_smiles"
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None, "parse_fail"
    try:
        m = rdMolStandardize.FragmentParent(m)
        m = rdMolStandardize.Uncharger().uncharge(m)
        Chem.SanitizeMol(m)
    except Exception:  # noqa: BLE001
        return None, "standardize_fail"
    return m, None


def main():
    os.makedirs(RAW, exist_ok=True)
    os.makedirs(DER, exist_ok=True)
    os.makedirs(RES, exist_ok=True)

    pages, total = [], None
    offset = 0
    while True:
        path = os.path.join(RAW, f"page_offset{offset}.json")
        if os.path.exists(path) and os.path.getsize(path) > 1000:
            data = json.load(open(path, encoding="utf-8"))
            print(f"cached {path}", flush=True)
        else:
            data = fetch(offset)
            data["retrieved_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
            with open(path, "w", encoding="utf-8") as fh:
                json.dump(data, fh, ensure_ascii=False)
            print(f"fetched offset={offset} n={len(data['molecules'])}", flush=True)
            time.sleep(1.0)
        total = data["page_meta"]["total_count"]
        pages.extend(data["molecules"])
        offset += 1000
        if offset >= total:
            break

    rows, drop = [], {"no_smiles": 0, "parse_fail": 0, "standardize_fail": 0,
                      "not_small_molecule": 0, "element_filter": 0,
                      "duplicate_full_inchikey": 0, "mw_gt_1000": 0}
    seen = {}
    for m in pages:
        if m.get("molecule_type") != "Small molecule":
            drop["not_small_molecule"] += 1
            continue
        s = (m.get("molecule_structures") or {}).get("canonical_smiles")
        mol, err = standardize(s)
        if err:
            drop[err] += 1
            continue
        if not all(a.GetAtomicNum() in ALLOWED for a in mol.GetAtoms()):
            drop["element_filter"] += 1
            continue
        # 去重键用完整 InChIKey：盐型/水合物经 FragmentParent 后归一为同一分子；
        # 但立体异构体（对映体）是不同药物，不能按 14 位骨架合并。
        ik = Chem.MolToInchiKey(mol)
        if ik in seen:
            drop["duplicate_full_inchikey"] += 1
            continue
        mw = Descriptors.MolWt(mol)
        if mw > 1000:
            drop["mw_gt_1000"] += 1
            continue
        seen[ik] = m["molecule_chembl_id"]
        skel = ik[:14]
        props = m.get("molecule_properties") or {}
        rows.append(dict(
            lig_id=f"M_{ik}", inchikey=ik, group="approved",
            chembl_id=m["molecule_chembl_id"], pref_name=m.get("pref_name"),
            molecule_type=m.get("molecule_type"), max_phase=m.get("max_phase"),
            first_approval=m.get("first_approval"), withdrawn=bool(m.get("withdrawn_flag")),
            oral=m.get("oral"), parenteral=m.get("parenteral"),
            natural_product=m.get("natural_product"), polymer=m.get("polymer_flag"),
            inorganic=m.get("inorganic_flag"),
            full_mwt=props.get("full_mwt"), alogp=props.get("alogp"),
            smiles=Chem.MolToSmiles(mol), smiles_std=Chem.MolToSmiles(mol),
            mw=round(Descriptors.MolWt(mol), 2)))

    df = pd.DataFrame(rows).drop_duplicates(subset=["lig_id"])
    out = os.path.join(DER, "approved_drugs.csv")
    df.to_csv(out, index=False, encoding="utf-8-sig")

    rep = dict(source="ChEMBL REST molecule.json", query=PARAMS, total_count=total,
               pages=len(pages), retrieved_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
               kept=len(df), dropped=drop,
               n_withdrawn=int(df["withdrawn"].sum()),
               n_oral=int((df["oral"] == "Y").sum()) if "oral" in df else None,
               mw_max=float(df["mw"].max()) if len(df) else None,
               mw_median=float(df["mw"].median()) if len(df) else None,
               sha256=hashlib.sha256(open(out, "rb").read()).hexdigest())
    jp = os.path.join(RES, "approved_library_report.json")
    json.dump(rep, open(jp, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print(json.dumps(rep, ensure_ascii=False, indent=2))
    print("written:", out, jp)


if __name__ == "__main__":
    main()
