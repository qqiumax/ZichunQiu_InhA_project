# -*- coding: utf-8 -*-
"""等药筛两个受体跑完 → 自动跑 analyze_approved.py → 写 screen_done.json。

用法：python run_screen_analysis.py（后台）
"""
import json
import os
import subprocess
import sys
import time

PROJ = r"D:\workspace\research\protein-drug"
DOCK = os.path.join(PROJ, "04-data", "derived", "docking")
RES = os.path.join(PROJ, "05-analysis", "results")
PY = sys.executable
LOG = os.path.join(RES, "screen_analysis.log")


def log(msg):
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(LOG, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def ready(rec):
    p = os.path.join(DOCK, f"dock_{rec}_approved.csv")
    if not os.path.exists(p):
        return False
    try:
        import pandas as pd
        d = pd.read_csv(p)
        return len(d) >= 2140 and d["score"].notna().sum() >= 2100
    except Exception:
        return False


def main():
    t0 = time.time()
    while time.time() - t0 < 12 * 3600:
        if ready("4DRE") and ready("4DTI"):
            log("both receptors ready")
            break
        time.sleep(60)
    else:
        log("TIMEOUT waiting for screening")
        json.dump(dict(status="timeout", at=time.strftime("%FT%TZ", time.gmtime())),
                  open(os.path.join(RES, "screen_done.json"), "w", encoding="utf-8"))
        return

    r = subprocess.run([PY, "-u", os.path.join(PROJ, "05-analysis", "analyze_approved.py")],
                       cwd=PROJ, capture_output=True, encoding="utf-8", errors="replace", timeout=3600)
    log(f"analyze rc={r.returncode} {(r.stdout or '')[-800:]}")
    if r.returncode == 0:
        json.dump(dict(status="done", at=time.strftime("%FT%TZ", time.gmtime())),
                  open(os.path.join(RES, "screen_done.json"), "w", encoding="utf-8"))
        log("DONE")


if __name__ == "__main__":
    main()
