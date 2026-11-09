# -*- coding: utf-8 -*-
"""Unified MD figures: all RMSD, contacts, and MM-GBSA systems."""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd

FD = r"D:\workspace\research\protein-drug\06-figures\data"
CP = r"D:\workspace\research\protein-drug\06-figures\contact-panels"
OUT = r"D:\workspace\research\protein-drug\06-figures\main"

plt.rcParams.update({
    "font.family": "Arial", "font.size": 7,
    "axes.labelsize": 7.5, "axes.titlesize": 7.5,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
    "legend.fontsize": 6.3, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "pdf.fonttype": 42, "ps.fonttype": 42,
    "savefig.dpi": 600, "savefig.bbox": "tight",
})

COL = {
    "entrectinib": "#0072B2",
    "entrectinib_4DTI": "#D55E00",
    "prednisolone_4DTI": "#CC79A7",
    "GSK625": "#009E73",
    "decoy": "#882255",
    "metformin": "#56B4E9",
    "566": "#999999",
}
TITLE = {
    "entrectinib": "Entrectinib · 4DRE",
    "entrectinib_4DTI": "Entrectinib · 4DTI (S94A)",
    "prednisolone_4DTI": "Prednisolone tebutate · 4DTI",
    "GSK625": "GSK625 · 5JFO",
    "decoy": "Decoy · 4DRE",
    "metformin": "Metformin · 4DRE",
    "566": "566 · 4DRE (supp.)",
}


def panel(ax, letter, x=-0.16, y=1.10):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="top", ha="left")


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), dpi=600)
    plt.close(fig)
    print("saved", name)


# ---------------- RMSD ----------------
r1 = pd.read_csv(os.path.join(FD, "md_rmsd.csv"))
r2 = pd.read_csv(os.path.join(FD, "md_rmsd_extra.csv"))
r = pd.concat([r1, r2], ignore_index=True)

fig = plt.figure(figsize=(7.1, 6.0), layout="constrained")
gs = GridSpec(3, 3, figure=fig, height_ratios=[1.15, 1, 1])


def draw(ax, key, tmax, show_ylabel=False):
    d = r[r.system == key].sort_values("time_ns")
    dt = float(np.median(np.diff(d.time_ns)))
    win = max(2, int(round(2.0 / dt)))
    sm = d.rmsd_A.rolling(win, center=True, min_periods=1).mean()
    ax.axvspan(0, 5, color="#E6E6E6", alpha=0.7, lw=0)
    ax.plot(d.time_ns, d.rmsd_A, color=COL[key], lw=0.45, alpha=0.22)
    ax.plot(d.time_ns, sm, color=COL[key], lw=1.2)
    ax.plot([0, tmax], [3, 3], ls=(0, (2, 2)), color="#777777", lw=0.7)
    late = d[d.time_ns >= d.time_ns.max() - 10]
    ax.text(0.04, 0.91, f"last 10 ns: {late.rmsd_A.mean():.2f} Å",
            transform=ax.transAxes, fontsize=5.6, va="top")
    ax.set(xlim=(0, tmax), ylim=(0, 6))
    if show_ylabel:
        ax.set_ylabel("Ligand RMSD (Å)")


ax_top = fig.add_subplot(gs[0, :])
draw(ax_top, "entrectinib", 197, True)
ax_top.set_xlabel("Time (ns)")
ax_top.set_title(TITLE["entrectinib"] + " (WT, NADH)", pad=3)
panel(ax_top, "a", x=-0.06)

bottom = ["entrectinib_4DTI", "prednisolone_4DTI", "GSK625",
          "decoy", "metformin", "566"]
for i, key in enumerate(bottom):
    ax = fig.add_subplot(gs[1 + i // 3, i % 3])
    draw(ax, key, 20, i % 3 == 0)
    ax.set_xlabel("Time (ns)")
    ax.set_title(TITLE[key], pad=3)
    panel(ax, "bcdefg"[i], x=-0.18 if i % 3 == 0 else -0.12)
save(fig, "md_fig_rmsd_all")

# ---------------- Contacts ----------------
left = [("entrectinib", "Entrectinib"), ("decoy", "Decoy"),
        ("metformin", "Metformin"), ("566", "566")]
right = [("entrectinib_4DTI", "Entrectinib"), ("prednisolone_4DTI", "Prednisolone"),
         ("GSK625", "GSK625")]
dfs = {}
for key, _ in left + right:
    dfs[key] = pd.read_csv(os.path.join(CP, f"window5-20_{key}.csv")).set_index("residue")["occ_5_20ns_pct"]
allres = sorted(set().union(*[set(d.index) for d in dfs.values()]),
                key=lambda x: -max([dfs[k].get(x, 0) for k in dfs]))[:20]
keys = {"TYR158", "PHE149", "MET98", "GLY96", "ILE16"}

fig = plt.figure(figsize=(7.1, 4.9), layout="constrained")
gs = GridSpec(1, 3, width_ratios=[len(left), len(right), 0.10], wspace=0.05)
axl = fig.add_subplot(gs[0, 0])
axr = fig.add_subplot(gs[0, 1], sharey=axl)
cax = fig.add_subplot(gs[0, 2])

for ax, group, labels, letter, title in [
        (axl, left, [x[1] for x in left], "a", "4DRE / NADH"),
        (axr, right, [x[1] for x in right], "b", "4DTI or 5JFO")]:
    M = np.array([[dfs[k].get(res, 0) for k, _ in group] for res in allres])
    im = ax.imshow(M, aspect="auto", cmap="cividis", vmin=0, vmax=100)
    ax.set_xticks(range(len(group)))
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_title(title, pad=4)
    for i in range(len(allres)):
        for j in range(len(group)):
            ax.text(j, i, f"{M[i,j]:.0f}", ha="center", va="center", fontsize=4.8,
                    color="white" if M[i, j] < 45 else "black")
    ax.spines[:].set_visible(False)
    panel(ax, letter, x=-0.22 if letter == "a" else -0.10, y=1.07)

axl.set_yticks(range(len(allres)))
axl.set_yticklabels([f"$\\bf{{{res}}}$" if res in keys else res for res in allres], fontsize=6.0)
axr.tick_params(labelleft=False)
fig.colorbar(im, cax=cax, label="Occupancy (%)")
save(fig, "md_fig_contacts_all")

# ---------------- MM-GBSA ----------------
g = pd.read_csv(os.path.join(FD, "mmgbsa_summary_all.csv"))
order = ["entrectinib", "entrectinib_4DTI", "prednisolone_4DTI",
         "GSK625", "decoy", "metformin", "566"]
g = g.set_index("system").loc[order].reset_index()
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.1, 3.7), layout="constrained",
                               gridspec_kw={"width_ratios": [1.3, 1.0]})
y = np.arange(len(g))[::-1]
for i, row in g.iterrows():
    ax1.barh(y[i] + 0.17, row.dg_mmgbsa, height=0.30, color="#0072B2", lw=0)
    ax1.barh(y[i] - 0.17, row.ie, height=0.30, color="#E69F00", lw=0)
    ax1.plot(row.eff_dg, y[i], marker="D", ms=3.4, color="#111111", zorder=5)
    ax1.annotate(f"{row.eff_dg:.1f}", (row.eff_dg, y[i]), xytext=(0, 7),
                 textcoords="offset points", ha="center", fontsize=5.4,
                 bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.85))
ax1.plot([0, 0], [-0.6, len(g) - 0.4], color="#444444", lw=0.6)
ax1.set_yticks(y)
ax1.set_yticklabels([f"{s} ({r})" + ("*" if role == "supplementary" else "")
                     for s, r, role in zip(g.system, g.receptor, g.role)], fontsize=6.2)
ax1.set(xlabel="Energy (kcal/mol)", xlim=(-62, 30))
ax1.legend(handles=[
    plt.Line2D([], [], color="#0072B2", lw=5, label="ΔG$_{MM/GBSA}$"),
    plt.Line2D([], [], color="#E69F00", lw=5, label="IE (−TΔS)"),
    plt.Line2D([], [], marker="D", ls="", color="#111111", label="Effective ΔG"),
], loc="lower left", handlelength=1.0)
panel(ax1, "a")

rcol = {"4DRE": "#0072B2", "4DTI": "#D55E00", "5JFO": "#009E73"}
for _, row in g.iterrows():
    marker = "s" if row.receptor == "5JFO" else ("^" if row.role == "supplementary" else "o")
    ax2.scatter(row.mw, row.eff_dg, s=34, color=rcol.get(row.receptor, "#777777"),
                edgecolor="white" if row.role != "supplementary" else "#333333",
                linewidth=0.6, marker=marker, zorder=3)
    ax2.annotate(row.system, (row.mw, row.eff_dg), xytext=(4, 4), textcoords="offset points", fontsize=5.4)
ax2.axhline(0, color="#444444", lw=0.6)
ax2.set(xlabel="Molecular weight (Da)", ylabel="Effective ΔG (kcal/mol)", xlim=(80, 640), ylim=(-52, 2))
ax2.legend(handles=[
    plt.Line2D([], [], marker="o", ls="", color=rcol["4DRE"], label="4DRE"),
    plt.Line2D([], [], marker="o", ls="", color=rcol["4DTI"], label="4DTI (S94A)"),
    plt.Line2D([], [], marker="s", ls="", color=rcol["5JFO"], label="5JFO (orthogonal)"),
    plt.Line2D([], [], marker="^", ls="", mfc="none", mec="#333333", label="supplementary"),
], loc="lower right", handlelength=0.8)
panel(ax2, "b")
save(fig, "md_fig_mmgbsa_all")
