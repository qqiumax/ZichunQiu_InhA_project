"""Clean + deduplicate InhA activity records into a modelling table.

Input : 04-data/raw/chembl/CHEMBL1849/activities_all.json   (raw, untouched)
Output: 04-data/derived/inha_ic50_clean.csv     one row per unique molecule (median pIC50)
        04-data/derived/inha_long_clean.csv    one row per accepted activity record
        05-analysis/results/inha_cleaning_report.json summary + every dropped-record reason

Rules (fixed BEFORE looking at hit rates; recorded in the report):
  keep standard_type == IC50, standard_flag == 1, standard_relation == '=', units in {nM}
  pIC50 = 9 - log10(IC50_nM)
  drop unparsable SMILES; standardise = fragment_parent + uncharge, then canonical SMILES
  deduplicate per molecule: median pIC50; keep min/max/spread and n_measurements
"""
import json
import math
import os
import sys
from collections import Counter, defaultdict

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

RAW = r"D:\workspace\research\protein-drug\04-data\raw\chembl\CHEMBL1849\activities_all.json"
DERIVED = r"D:\workspace\research\protein-drug\04-data\derived"
RESULTS = r"D:\workspace\research\protein-drug\05-analysis\results"
os.makedirs(DERIVED, exist_ok=True)
os.makedirs(RESULTS, exist_ok=True)

with open(RAW, encoding="utf-8") as fh:
    payload = json.load(fh)
acts = payload["activities"]

UN = rdMolStandardize.Uncharger()
report = {"source": RAW, "retrieved_utc": payload.get("retrieved_utc"),
          "n_records_in": len(acts), "dropped": Counter(), "kept": 0}

rows = []
for a in acts:
    def drop(reason):
        report["dropped"][reason] += 1

    if a.get("standard_type") != "IC50":
        drop("type_not_IC50"); continue
    if str(a.get("standard_flag")) != "1":
        drop("flag_not_1"); continue
    if a.get("standard_relation") != "=":
        drop(f"relation_{a.get('standard_relation')}"); continue
    if a.get("standard_units") != "nM":
        drop(f"units_{a.get('standard_units')}"); continue
    try:
        ic50 = float(a.get("standard_value"))
    except (TypeError, ValueError):
        drop("value_not_numeric"); continue
    if ic50 <= 0:
        drop("value_nonpositive"); continue
    smi = a.get("canonical_smiles")
    if not smi:
        drop("no_smiles"); continue

    m = Chem.MolFromSmiles(smi)
    if m is None:
        drop("smiles_unparsable"); continue
    try:
        m2 = rdMolStandardize.FragmentParent(m)
        m2 = UN.uncharge(m2)
        std = Chem.MolToSmiles(m2)
        key = Chem.MolToInchiKey(m2)
    except Exception:
        drop("standardisation_failed"); continue
    if not std:
        drop("empty_after_standardisation"); continue

    rows.append({
        "molecule_chembl_id": a.get("molecule_chembl_id"),
        "canonical_smiles_raw": smi,
        "smiles_std": std,
        "inchikey": key,
        "ic50_nM": ic50,
        "pIC50": 9.0 - math.log10(ic50),
        "assay_chembl_id": a.get("assay_chembl_id"),
        "document_chembl_id": a.get("document_chembl_id"),
    })

report["kept"] = len(rows)
long_df = pd.DataFrame(rows)
long_df.to_csv(os.path.join(DERIVED, "inha_long_clean.csv"), index=False, encoding="utf-8-sig")

# ---- deduplicate per molecule (InChIKey full, then skeleton fallback) ----
long_df["inchikey_skel"] = long_df["inchikey"].str[:14]
grp = long_df.groupby("inchikey_skel")
agg = grp.agg(
    n_measurements=("pIC50", "size"),
    pIC50_median=("pIC50", "median"),
    pIC50_min=("pIC50", "min"),
    pIC50_max=("pIC50", "max"),
    n_assays=("assay_chembl_id", "nunique"),
    n_documents=("document_chembl_id", "nunique"),
    smiles_std=("smiles_std", "first"),
).reset_index()
agg["pIC50_spread"] = agg["pIC50_max"] - agg["pIC50_min"]

scaf = []
for s in agg["smiles_std"]:
    m = Chem.MolFromSmiles(s)
    try:
        scaf.append(MurckoScaffold.MurckoScaffoldSmiles(mol=m, includeChirality=False))
    except Exception:
        scaf.append("")
agg["murcko_scaffold"] = scaf

agg.to_csv(os.path.join(DERIVED, "inha_ic50_clean.csv"), index=False, encoding="utf-8-sig")

report["n_unique_molecules"] = int(len(agg))
report["n_unique_scaffolds"] = int(agg["murcko_scaffold"].nunique())
report["pIC50_median_mean"] = round(float(agg["pIC50_median"].mean()), 3)
report["pIC50_median_sd"] = round(float(agg["pIC50_median"].std()), 3)
report["pIC50_median_min"] = round(float(agg["pIC50_median"].min()), 3)
report["pIC50_median_max"] = round(float(agg["pIC50_median"].max()), 3)
report["n_with_multiple_measurements"] = int((agg["n_measurements"] > 1).sum())
report["median_spread_where_replicated"] = (
    round(float(agg.loc[agg["n_measurements"] > 1, "pIC50_spread"].median()), 3)
    if (agg["n_measurements"] > 1).any() else None)
report["largest_scaffold_cluster"] = int(agg["murcko_scaffold"].value_counts().iloc[0]) if len(agg) else 0
report["top_scaffold_share"] = (
    round(float(agg["murcko_scaffold"].value_counts().iloc[0] / len(agg)), 3) if len(agg) else None)

with open(os.path.join(RESULTS, "inha_cleaning_report.json"), "w", encoding="utf-8") as fh:
    json.dump(report, fh, ensure_ascii=False, indent=2)

print(json.dumps({k: (dict(v) if isinstance(v, Counter) else v) for k, v in report.items()},
                 ensure_ascii=False, indent=2))
print("\nwrote:", os.path.join(DERIVED, "inha_ic50_clean.csv"))
print(long_df[["pIC50"]].describe().to_string())
