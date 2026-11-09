"""Download ALL ChEMBL activity records for Mtb InhA (CHEMBL1849).

Provenance: ChEMBL REST v36/37 (https://www.ebi.ac.uk/chembl/api/data), target CHEMBL1849
(UniProt P9WGR1, Enoyl-[acyl-carrier-protein] reductase [NADH], Mtb, tax_id 1773).
Raw pages are saved verbatim; no filtering here (filtering happens in 05-analysis, so it stays auditable).
"""
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request

UA = {"User-Agent": "research-ingest/1.0"}
BASE = "https://www.ebi.ac.uk/chembl/api/data"
TID = "CHEMBL1849"
OUT = r"D:\workspace\research\protein-drug\04-data\raw\chembl\CHEMBL1849"
PAGE = 1000


def get(url, timeout=120, tries=10):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=timeout) as f:
                return json.loads(f.read().decode())
        except Exception as e:
            if i == tries - 1:
                raise
            print(f"    retry {i+1}/{tries-1} ({type(e).__name__})", flush=True)
            time.sleep(min(30, 4 + 3 * i))


os.makedirs(OUT, exist_ok=True)
total = None
offset = 0
page_no = 0
all_acts = []

while total is None or offset < total:
    q = urllib.parse.urlencode({"target_chembl_id": TID, "limit": PAGE, "offset": offset})
    url = f"{BASE}/activity.json?{q}"
    print(f"[page {page_no}] offset={offset} GET", flush=True)
    d = get(url)
    total = d.get("page_meta", {}).get("total_count", 0)
    acts = d.get("activities", [])
    print(f"[page {page_no}] got {len(acts)} (total_count={total})", flush=True)

    with open(os.path.join(OUT, f"activities_page_{page_no:02d}_offset{offset}.json"), "w", encoding="utf-8") as fh:
        json.dump(d, fh, ensure_ascii=False, indent=1)

    all_acts.extend(acts)
    page_no += 1
    offset += PAGE
    if not acts:
        break
    time.sleep(2)

with open(os.path.join(OUT, "activities_all.json"), "w", encoding="utf-8") as fh:
    json.dump({"target_chembl_id": TID, "api": BASE, "retrieved_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
               "n": len(all_acts), "page_meta_total": total, "activities": all_acts}, fh, ensure_ascii=False)

print(f"\nDONE: {len(all_acts)} activity records -> {OUT}", flush=True)
