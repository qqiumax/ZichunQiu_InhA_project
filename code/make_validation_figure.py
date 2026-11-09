# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FD = r"D:\workspace\research\protein-drug\06-figures\data"
OUT = r"D:\workspace\research\protein-drug\06-figures\main"
plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.linewidth": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "savefig.dpi": 600, "savefig.bbox": "tight"})
C = {"4DRE": "#0072B2", "4DTI": "#D55E00", "A": "#D55E00", "B": "#E69F00", "C": "#0072B2", "decoy": "#777777"}

fig, axes = plt.subplots(2, 2, figsize=(7.1, 5.0), layout="constrained")

# a positive controls
ax = axes[0, 0]
p = pd.read_csv(os.path.join(FD, "positive_control_audit.csv"))
y = np.arange(len(p))[::-1]
for i, row in p.iterrows():
    ax.scatter(row.score_4DRE, y[i] + 0.12, s=28, color=C["4DRE"], zorder=3)
    ax.scatter(row.score_4DTI, y[i] - 0.12, s=28, color=C["4DTI"], marker="s", zorder=3)
thr4 = [(-9.422, "p1"), (-8.918, "p5"), (-8.163, "median")]
for x, lab in thr4:
    ax.axvline(x, color="#999999", lw=0.7, ls=(0, (2, 2)))
    ax.text(x, len(p) - 0.35, lab, rotation=90, fontsize=5, ha="center", va="bottom", color="#555555")
ax.set_yticks(y)
label_map = {"-": "none", "C": "lower"}
ax.set_yticklabels([f"{c} ({label_map.get(t, t)})" for c, t in zip(p.compound, p.tier_original)])
ax.set(xlabel="Vina score (kcal/mol)", xlim=(-9.7, -7.6))
ax.scatter([], [], color=C["4DRE"], s=20, label="4DRE")
ax.scatter([], [], color=C["4DTI"], s=20, marker="s", label="4DTI")
ax.legend(loc="lower right", handlelength=0.8)
panel = ax.text(-0.22, 1.04, "a", transform=ax.transAxes, fontsize=9, fontweight="bold")

# b recovery
ax = axes[0, 1]
v = pd.read_csv(os.path.join(FD, "tier_validation_leave_scaffold.csv"))
n = v.n_test.sum()
rec = [100 * v.A.sum() / n, 100 * v.B.sum() / n, 100 * v.C.sum() / n]
false = [0.89, 1.78, 12.76]
x = np.arange(3); w = 0.35
ax.bar(x - w / 2, rec, w, color="#0072B2", label="held-out actives")
ax.bar(x + w / 2, false, w, color="#999999", label="decoys")
for i, (a, b) in enumerate(zip(rec, false)):
    ax.text(i - w / 2, a + 0.6, f"{a:.1f}%", ha="center", fontsize=5.4)
    ax.text(i + w / 2, b + 0.6, f"{b:.1f}%", ha="center", fontsize=5.4)
ax.set_xticks(x); ax.set_xticklabels(["Top", "Middle", "Lower"])
ax.set(ylabel="Rate (%)", ylim=(0, 36))
ax.legend(loc="upper left", handlelength=0.8)
ax.text(-0.18, 1.04, "b", transform=ax.transAxes, fontsize=9, fontweight="bold")

# c stratified AUC
ax = axes[1, 0]
s = pd.read_csv(os.path.join(FD, "property_stratified_enrichment.csv"))
s = s[s.stratifier == "MW"]
for rec in ["4DRE", "4DTI"]:
    d = s[s.receptor == rec]
    ax.plot(range(len(d)), d.auc, marker="o", ms=3, lw=1.0, color=C[rec], label=rec)
ax.set_xticks(range(len(d))); ax.set_xticklabels(["0-300", "300-400", "400-500", ">=500"], rotation=25, ha="right")
ax.set(ylabel="AUC", ylim=(0.60, 0.82))
ax.legend(loc="lower left", handlelength=0.8)
ax.text(-0.20, 1.04, "c", transform=ax.transAxes, fontsize=9, fontweight="bold")

# d LE
ax = axes[1, 1]
h = pd.read_csv(os.path.join(FD, "approved_ligand_efficiency.csv"))
h = h[(h.heavy_atoms >= 15) & (h.mw >= 200)].dropna(subset=["le"])
ax.hist(h["le"], bins=40, color="#BBBBBB", edgecolor="white", lw=0.3)
for name, col in [("PREDNISOLONE TEBUTATE", "#D55E00"), ("ENTRECTINIB", "#0072B2")]:
    v = h.loc[h.pref_name == name, "le"]
    if len(v):
        ax.axvline(float(v.iloc[0]), color=col, lw=1.0, ls="--")
        ax.text(float(v.iloc[0]), ax.get_ylim()[1] * 0.9, name.split()[0].title(),
                rotation=90, fontsize=5, color=col, ha="right", va="top")
ax.set(xlabel="Ligand efficiency (kcal mol$^{-1}$ per heavy atom)", ylabel="Approved drugs", xlim=(-0.6, 0.7))
ax.text(-0.18, 1.04, "d", transform=ax.transAxes, fontsize=9, fontweight="bold")

for ext in ("pdf", "png"):
    fig.savefig(os.path.join(OUT, f"validation_fig.{ext}"), dpi=600)
print("saved validation_fig")
