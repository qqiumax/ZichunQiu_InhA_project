"""Environment self-check for this project.

Run this first after any restart, before trusting any downstream number.
Writes 00-admin/env-report.json and prints a human summary.

Design notes
  * Everything here is read-only except a tiny write test inside the ASCII scratch dir.
  * The non-ASCII path test is EXPECTED to fail (RDKit SDWriter); it is reported as
    "confirmed limitation" rather than an error, so we notice if it ever starts working.
"""
import json
import os
import platform
import subprocess
import sys
import time
import urllib.error
import urllib.request

PROJ = r"D:\workspace\research\protein-drug"
SCRATCH = r"D:\scratch\protein-drug"
OUT = os.path.join(PROJ, "00-admin", "env-report.json")
UA = {"User-Agent": "env-check/1.0"}
NET = [
    ("RCSB PDB", "https://files.rcsb.org/download/4DRE.pdb"),
    ("UniProt", "https://rest.uniprot.org/uniprotkb/P9WGR1.json"),
    ("ChEMBL REST", "https://www.ebi.ac.uk/chembl/api/data/target/CHEMBL1849.json"),
    ("Harvard Dataverse", "https://dataverse.harvard.edu/api/access/datafile/4259588"),
    ("PubChem", "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/aspirin/cids/JSON"),
]
PKGS = ["rdkit", "numpy", "pandas", "scipy", "sklearn", "Bio", "meeko", "torch", "openmm"]
CMDS = [
    ("vina", ["vina", "--version"]),
    ("obabel", ["obabel", "-V"]),
    ("mk_prepare_receptor", ["mk_prepare_receptor.exe", "--help"]),
    ("mk_prepare_ligand", ["mk_prepare_ligand.exe", "--help"]),
]

report = {"checked_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "os": platform.platform(),
          "python": sys.version.split()[0], "python_exe": sys.executable}

# --- packages -------------------------------------------------------------
import importlib
pk = {}
for m in PKGS:
    try:
        mod = importlib.import_module(m)
        pk[m] = getattr(mod, "__version__", "ok")
    except Exception as e:
        pk[m] = f"MISSING ({type(e).__name__})"
report["packages"] = pk

try:
    import torch
    report["cuda"] = {"available": bool(torch.cuda.is_available()),
                      "device": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None}
except Exception:
    report["cuda"] = {"available": None, "device": None}

# --- CLI tools ------------------------------------------------------------
tools = {}
for name, cmd in CMDS:
    try:
        r = subprocess.run(cmd, capture_output=True, encoding="utf-8", errors="replace", timeout=60)
        first = (r.stdout or r.stderr or "").strip().splitlines()
        tools[name] = first[0] if first else f"rc={r.returncode}"
    except Exception as e:
        tools[name] = f"ERROR {type(e).__name__}"
report["cli_tools"] = tools

# --- filesystem -----------------------------------------------------------
fs = {"scratch": SCRATCH, "scratch_exists": os.path.isdir(SCRATCH), "scratch_writable": False,
      "rdkit_write_ascii": None, "rdkit_write_nonascii": None}
os.makedirs(SCRATCH, exist_ok=True)
try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    RDLogger.DisableLog("rdApp.*")
    m = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.EmbedMolecule(m, randomSeed=42)
    block = Chem.MolToMolBlock(m)

    ascii_path = os.path.join(SCRATCH, "write_test.sdf")
    try:
        with open(ascii_path, "w", encoding="utf-8") as fh:
            fh.write(block + "$$$$\n")
        fs["rdkit_write_ascii"] = os.path.getsize(ascii_path) > 100
        os.remove(ascii_path)
    except Exception as e:
        fs["rdkit_write_ascii"] = f"ERR {type(e).__name__}"

    try:
        w = Chem.SDWriter(os.path.join(SCRATCH, "rdkit_sdw_test.sdf"))
        w.write(m)
        w.close()
        fs["scratch_writable"] = True
        os.remove(os.path.join(SCRATCH, "rdkit_sdw_test.sdf"))
    except Exception as e:
        fs["scratch_writable"] = f"ERR {type(e).__name__}"

    # non-ASCII path test: must use a REAL Chinese path. (The project tree is all-ASCII now,
    # so pointing this at 04-data/derived would trivially pass and prove nothing.)
    nonascii_dir = os.path.join(SCRATCH, "测试中文路径")
    os.makedirs(nonascii_dir, exist_ok=True)
    try:
        w = Chem.SDWriter(os.path.join(nonascii_dir, "x.sdf"))
        w.write(m)
        w.close()
        fs["rdkit_write_nonascii"] = "works (limitation no longer applies?)"
        os.remove(os.path.join(nonascii_dir, "x.sdf"))
    except Exception as e:
        fs["rdkit_write_nonascii"] = f"fails as expected ({type(e).__name__})"
    try:
        os.rmdir(nonascii_dir)
    except OSError:
        pass
except Exception as e:
    fs["note"] = f"rdkit probe failed: {type(e).__name__}"
report["filesystem"] = fs

# --- network --------------------------------------------------------------
net = {}
for name, url in NET:
    t0 = time.time()
    try:
        req = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=10) as f:
            net[name] = {"status": f.status, "ms": int((time.time() - t0) * 1000)}
    except urllib.error.HTTPError as e:
        net[name] = {"status": e.code, "ms": int((time.time() - t0) * 1000)}
    except Exception as e:
        net[name] = {"status": type(e).__name__, "ms": int((time.time() - t0) * 1000)}
report["network"] = net

try:
    req = urllib.request.Request("https://ifconfig.me/ip", headers=UA)
    with urllib.request.urlopen(req, timeout=10) as f:
        report["egress_ip"] = f.read(64).decode().strip()
except Exception as e:
    report["egress_ip"] = f"unknown ({type(e).__name__})"

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", encoding="utf-8") as fh:
    json.dump(report, fh, ensure_ascii=False, indent=2)

print(json.dumps(report, ensure_ascii=False, indent=2))
print("\nwritten:", OUT)
