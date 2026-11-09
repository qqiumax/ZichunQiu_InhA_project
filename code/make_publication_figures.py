# -*- coding: utf-8 -*-
"""出版级图套件：富集、对接分布、MD RMSD、接触热图、MM-GBSA 分解。"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

FD = r"D:\workspace\research\protein-drug\06-figures\data"
CP = r"D:\workspace\research\protein-drug\06-figures\contact-panels"
OUT = r"D:\workspace\research\protein-drug\06-figures\main"
os.makedirs(OUT, exist_ok=True)

C = {
    "4DRE": "#0072B2", "4DTI": "#D55E00", "4U0J": "#009E73",
    "entrectinib": "#0072B2", "566": "#E69F00", "decoy": "#CC79A7",
    "metformin": "#56B4E9", "GSK625": "#009E73",
    "gray": "#8C8C8C", "light": "#D9D9D9",
}

plt.rcParams.update({
    "font.family": "Arial", "font.size": 7,
    "axes.labelsize": 7.5, "axes.titlesize": 8,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "pdf.fonttype": 42, "ps.fonttype": 42,
    "savefig.dpi": 600, "savefig.bbox": "tight",
})


def panel(ax, letter, x=-0.14, y=1.08):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="top", ha="left")


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), dpi=600)
    plt.close(fig)
    print("saved", name)


# ---------------- Fig 1: enrichment ----------------
cur = pd.read_csv(os.path.join(FD, "enrichment_curves.csv"))
summ = pd.read_csv(os.path.join(FD, "enrichment_summary.csv")).set_index("receptor")

fig = plt.figure(figsize=(7.1, 2.55), layout="constrained")
gs = GridSpec(1, 3, figure=fig, width_ratios=[1.0, 1.0, 1.35])
ax1, ax2 = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
sub = gs[0, 2].subgridspec(2, 1, hspace=0.12)
ax3a, ax3b = fig.add_subplot(sub[0]), fig.add_subplot(sub[1])

for rec in ["4DRE", "4DTI", "4U0J"]:
    d = cur[(cur.receptor == rec) & (cur.kind == "roc")]
    ax1.plot(d.x, d.y, color=C[rec], lw=1.1,
             label=f"{rec}  AUC {summ.loc[rec,'auc']:.3f}")
ax1.plot([0, 1], [0, 1], ls=(0, (2, 2)), color=C["gray"], lw=0.7)
ax1.set(xlim=(0, 1), ylim=(0, 1), xlabel="False positive rate", ylabel="True positive rate")
ax1.legend(loc="lower right", handlelength=1.2)
panel(ax1, "a")

for rec in ["4DRE", "4DTI", "4U0J"]:
    d = cur[(cur.receptor == rec) & (cur.kind == "ef")]
    ax2.plot(d.x, d.y, color=C[rec], lw=1.1)
ax2.plot([1, 20], [1, 1], ls=(0, (2, 2)), color=C["gray"], lw=0.7)
ax2.set(xlim=(1, 20), ylim=(0.8, 2.1), xlabel="Top-ranked (%)", ylabel="Enrichment factor")
panel(ax2, "b")

ypos = np.arange(len(summ))[::-1]
for i, rec in enumerate(summ.index):
    r = summ.loc[rec]
    ax3a.errorbar(r.auc, ypos[i], xerr=[[r.auc - r.auc_lo], [r.auc_hi - r.auc]],
                  fmt="o", ms=3.4, color=C["4DRE"], lw=1.0, capsize=2)
    ax3b.errorbar(r.ef5, ypos[i], xerr=[[r.ef5 - r.ef5_lo], [r.ef5_hi - r.ef5]],
                  fmt="s", ms=3.2, color=C["4DTI"], lw=1.0, capsize=2)
for ax in (ax3a, ax3b):
    ax.set_yticks(ypos)
    ax.set_yticklabels(summ.index)
    ax.set_ylim(-0.6, len(summ) - 0.4)
ax3a.set(xlim=(0.62, 0.78), xlabel="AUC", xticks=[0.65, 0.70, 0.75])
ax3a.tick_params(labelbottom=False)
ax3b.set(xlim=(1.1, 1.95), xlabel="EF5 (top 5%)", xticks=[1.2, 1.4, 1.6, 1.8])
ax3a.set_title("Metric (95% CI)", pad=2, fontsize=7)
panel(ax3a, "c", x=-0.45)
save(fig, "fig1_retrospective_enrichment")

# ---------------- Fig 2: docking score distributions ----------------
dist = pd.read_csv(os.path.join(FD, "dock_distribution.csv"))
thr = pd.read_csv(os.path.join(FD, "dock_thresholds.csv")).set_index("receptor")
hits = pd.read_csv(os.path.join(r"D:\workspace\research\protein-drug\05-analysis\results", "approved_hits.csv"))
try:
    from scipy.stats import gaussian_kde
except Exception:
    gaussian_kde = None

fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.7), layout="constrained", sharey=True)
for ax, rec in zip(axes, ["4DRE", "4DTI"]):
    d = dist[dist.receptor == rec]
    xg = np.linspace(-12.5, -4.0, 400)
    for grp, color, lab in [("decoy", C["gray"], "Decoys"), ("active", C["4DRE"], "Actives")]:
        v = d[d.group == grp].score.values
        if gaussian_kde is not None:
            yg = gaussian_kde(v, bw_method=0.25)(xg)
            ax.plot(xg, yg, color=color, lw=1.0, label=lab)
            ax.fill_between(xg, 0, yg, color=color, alpha=0.10, lw=0)
        else:
            ax.hist(v, bins=60, density=True, histtype="step", lw=0.9, color=color, label=lab)
    ap = d[d.group == "approved"].score.values
    if gaussian_kde is not None and len(ap) > 10:
        yg = gaussian_kde(ap, bw_method=0.30)(xg)
        ax.plot(xg, yg, color="#E69F00", lw=0.9, ls=(0, (4, 2)), label="Approved drugs")
        ax.fill_between(xg, 0, yg, color="#E69F00", alpha=0.06, lw=0)
    for name, key, ls in [("p1", "p1", (0, (1, 1.5))), ("p5", "p5", (0, (3, 1.5))),
                          ("median", "active_median", (0, (5, 2)))]:
        x = thr.loc[rec, key]
        ax.plot([x, x], [0, 0.52], ls=ls, color="#444444", lw=0.8)
        ax.text(x, 0.528, name, rotation=90, va="bottom", ha="center", fontsize=5.0, color="#444444")
    score_col = "score_wt" if rec == "4DRE" else "score_mu"
    s_e = hits.loc[hits.pref_name == "ENTRECTINIB", score_col].values
    s_p = hits.loc[hits.pref_name == "PREDNISOLONE TEBUTATE", score_col].values
    if len(s_e) and len(s_p):
        ax.text(0.98, 0.72, f"Top-band drugs\nentrectinib  {s_e[0]:.2f}\nprednisolone  {s_p[0]:.2f}",
                transform=ax.transAxes, ha="right", va="top", fontsize=5.3,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.85))
    ax.text(0.02, 0.96, rec, transform=ax.transAxes, fontsize=9, fontweight="bold", va="top")
    ax.set(xlabel="Vina score (kcal/mol)", xlim=(-12.5, -4.0), ylim=(0, 0.62))
axes[0].set_ylabel("Density")
axes[0].legend(loc="upper left", handlelength=1.0, ncol=1, bbox_to_anchor=(0.0, 0.90))
panel(axes[0], "a")
panel(axes[1], "b")
save(fig, "fig2_docking_distributions")

# ---------------- Fig 3: MD ligand RMSD ----------------
r = pd.read_csv(os.path.join(FD, "md_rmsd.csv"))
fig = plt.figure(figsize=(7.1, 4.6), layout="constrained")
gs = GridSpec(2, 3, figure=fig, height_ratios=[1.0, 1.0])

def draw_rmsd(ax, name, tmax=None, color=None):
    d = r[r.system == name]
    if tmax:
        d = d[d.time_ns <= tmax]
    dt = float(np.median(np.diff(d.time_ns)))
    win = max(2, int(round(2.0 / dt)))
    sm = d.rmsd_A.rolling(win, center=True, min_periods=1).mean()
    ax.axvspan(0, 5, color=C["light"], alpha=0.45, lw=0)
    ax.plot(d.time_ns, d.rmsd_A, color=color, lw=0.45, alpha=0.25)
    ax.plot(d.time_ns, sm, color=color, lw=1.2)
    ax.plot([d.time_ns.min(), d.time_ns.max()], [3, 3], ls=(0, (2, 2)), color=C["gray"], lw=0.7)
    late = d[d.time_ns >= d.time_ns.max() - 10]
    ax.text(0.03, 0.92, f"last 10 ns: {late.rmsd_A.mean():.2f} Å", transform=ax.transAxes,
            fontsize=5.8, va="top")
    ax.set(ylim=(0, 6), xlim=(0, tmax if tmax else d.time_ns.max() * 1.01))
    ax.spines[["top", "right"]].set_visible(False)

ax_top = fig.add_subplot(gs[0, :])
draw_rmsd(ax_top, "entrectinib", tmax=197, color=C["entrectinib"])
ax_top.set(xlabel="Time (ns)", ylabel="Ligand RMSD (Å)")
ax_top.set_title("Entrectinib · 4DRE (WT, NADH)", pad=3)
panel(ax_top, "a", x=-0.06)

for i, name in enumerate(["GSK625", "metformin", "decoy"]):
    ax = fig.add_subplot(gs[1, i])
    draw_rmsd(ax, name, tmax=20, color=C[name])
    ax.set(xlabel="Time (ns)", ylabel="Ligand RMSD (Å)" if i == 0 else "")
    ax.set_title({"GSK625": "GSK625 · 5JFO (orthogonal)", "metformin": "Metformin · 4DRE",
                  "decoy": "Decoy · 4DRE"}[name], pad=3)
    panel(ax, "bc"[i] if i < 2 else "d", x=-0.18 if i == 0 else -0.10)
save(fig, "fig3_md_ligand_rmsd")

# ---------------- Fig 4: contact heatmap ----------------
systems = ["entrectinib", "decoy", "metformin", "GSK625"]
labels = ["Entrectinib", "Decoy", "Metformin", "GSK625\n(5JFO)"]
dfs = {}
for s in systems:
    p = os.path.join(CP, f"window5-20_{s}.csv")
    dfs[s] = pd.read_csv(p).set_index("residue")["occ_5_20ns_pct"]
allres = sorted(set().union(*[set(d.index) for d in dfs.values()]),
                key=lambda x: -max([dfs[s].get(x, 0) for s in systems]))[:20]
M = np.array([[dfs[s].get(res, 0) for s in systems] for res in allres])
keys = {"TYR158", "PHE149", "MET98", "GLY96", "ILE16"}
fig, ax = plt.subplots(figsize=(7.1, 3.9), layout="constrained")
im = ax.imshow(M, aspect="auto", cmap="cividis", vmin=0, vmax=100)
ax.set_xticks(range(len(systems)))
ax.set_xticklabels(labels, rotation=25, ha="right")
ax.set_yticks(range(len(allres)))
ax.set_yticklabels([f"$\\bf{{{res}}}$" if res in keys else res for res in allres], fontsize=6.3)
for i in range(len(allres)):
    for j in range(len(systems)):
        ax.text(j, i, f"{M[i,j]:.0f}", ha="center", va="center", fontsize=5.2,
                color="white" if M[i, j] < 45 else "black")
ax.set_title("Contact occupancy, 5–20 ns (4 Å, heavy atoms)", pad=4)
ax.set_xlabel("")
fig.colorbar(im, ax=ax, label="Occupancy (%)", fraction=0.035, pad=0.02)
ax.spines[:].set_visible(False)
save(fig, "fig4_contact_heatmap")

# ---------------- Fig 5: MM-GBSA + IE ----------------
g = pd.read_csv(os.path.join(FD, "mmgbsa_summary.csv"))
g = g.sort_values("eff_dg")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.1, 3.0), layout="constrained",
                               gridspec_kw={"width_ratios": [1.25, 1.0]})
y = np.arange(len(g))[::-1]
for i, (_, row) in enumerate(g.iterrows()):
    ax1.barh(y[i] + 0.17, row.dg_mmgbsa, height=0.30, color=C["4DRE"], lw=0)
    ax1.barh(y[i] - 0.17, row.ie, height=0.30, color="#E69F00", lw=0)
    ax1.plot(row.eff_dg, y[i], marker="D", ms=3.6, color="#111111", zorder=5)
    ax1.annotate(f"{row.eff_dg:.1f}", (row.eff_dg, y[i]), xytext=(0, 7),
                 textcoords="offset points", ha="center", fontsize=5.6,
                 bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
ax1.plot([0, 0], [-0.6, len(g) - 0.4], color="#444444", lw=0.6)
ax1.set_yticks(y)
ax1.set_yticklabels([f"{r.system}" for _, r in g.iterrows()])
ax1.set(xlabel="Energy (kcal/mol)", xlim=(-62, 30))
ax1.legend(handles=[
    plt.Line2D([], [], color=C["4DRE"], lw=5, label="ΔG$_{MM/GBSA}$"),
    plt.Line2D([], [], color="#E69F00", lw=5, label="IE (−TΔS)"),
    plt.Line2D([], [], marker="D", ls="", color="#111111", label="Effective ΔG"),
], loc="lower left", handlelength=1.0, ncol=1)
panel(ax1, "a")

for _, row in g.iterrows():
    filled = row.status == "full"
    ax2.scatter(row.mw, row.eff_dg, s=34, color=C.get(row.system, C["gray"]),
                edgecolor="white" if filled else "#333333", linewidth=0.6,
                marker="o" if filled else "s", alpha=1.0, zorder=3,
                facecolors=C.get(row.system, C["gray"]) if filled else "none")
    ax2.annotate(row.system, (row.mw, row.eff_dg), xytext=(4, 4), textcoords="offset points", fontsize=5.8)
ax2.axhline(0, color="#444444", lw=0.6)
ax2.set(xlabel="Molecular weight (Da)", ylabel="Effective ΔG (kcal/mol)", xlim=(80, 640), ylim=(-52, 2))
ax2.legend(handles=[
    plt.Line2D([], [], marker="o", ls="", color=C["gray"], label="same system (4DRE/NADH)"),
    plt.Line2D([], [], marker="s", ls="", mfc="none", mec="#333333", label="orthogonal (5JFO/NAD+)"),
], loc="lower right", handlelength=0.9)
panel(ax2, "b")
save(fig, "fig5_mmgbsa_ie")

# ---------------- Fig 6: docking screen summary ----------------
try:
    from scipy.stats import spearmanr
except Exception:
    spearmanr = None

h = pd.read_csv(os.path.join(FD, "approved_hits_fig.csv"))
counts = pd.read_csv(os.path.join(FD, "tier_counts.csv")).set_index("tier")
mut = pd.read_csv(os.path.join(FD, "mutation_dds.csv"))
top = pd.read_csv(os.path.join(FD, "top_hits.csv"))
TIER_C = {"A": "#D55E00", "B": "#E69F00", "C": "#0072B2", "-": "#BDBDBD"}

fig, axes = plt.subplots(2, 2, figsize=(7.1, 5.2), layout="constrained")
ax = axes[0, 0]
order = ["A", "B", "C", "-"]
vals = counts.reindex(order)["n"].values
ax.barh(np.arange(len(order))[::-1], vals, color=[TIER_C[t] for t in order], height=0.62)
for i, v in enumerate(vals):
    ax.text(v * 1.15, len(order) - 1 - i, f"{int(v):,}", va="center", fontsize=6)
ax.set_yticks(np.arange(len(order))[::-1])
ax.set_yticklabels(["A", "B", "C", "no tier"])
ax.set(xscale="log", xlabel="Approved drugs (log scale)")
ax.set_xlim(1, 6000)
panel(ax, "a")

ax = axes[0, 1]
colors = h.tier.fillna("-").map(TIER_C)
ax.scatter(h.mw, h.score_wt, s=3.5, c=colors, alpha=0.35, linewidth=0)
aa = h[h.tier == "A"]
ax.scatter(aa.mw, aa.score_wt, s=22, c=TIER_C["A"], edgecolor="white", linewidth=0.6, zorder=5)
for _, r in aa.iterrows():
    if "PREDNISOLONE" in r.pref_name.upper():
        ax.annotate(r.pref_name.title(), (r.mw, r.score_wt), xytext=(-8, 10),
                    textcoords="offset points", fontsize=5.4, ha="right")
    else:
        ax.annotate(r.pref_name.title(), (r.mw, r.score_wt), xytext=(6, -10),
                    textcoords="offset points", fontsize=5.4)
if spearmanr is not None:
    rho = spearmanr(h.mw, h.score_wt, nan_policy="omit").statistic
    ax.text(0.03, 0.06, f"Spearman ρ = {rho:.2f}", transform=ax.transAxes, fontsize=5.8)
ax.legend(handles=[plt.Line2D([], [], marker="o", ls="", color=TIER_C[t],
                              label={"A": "A", "B": "B", "C": "C", "-": "no tier"}[t])
                   for t in ["A", "B", "C", "-"]],
          loc="upper right", fontsize=5.2, handlelength=0.6, labelspacing=0.25)
ax.set(xlabel="Molecular weight (Da)", ylabel="4DRE Vina score (kcal/mol)", xlim=(80, 650), ylim=(-12, -4))
panel(ax, "b")

ax = axes[1, 0]
data = [mut[mut.group == "active"].dds.values, mut[mut.group == "decoy"].dds.values,
        h.dds_calc.dropna().values]
bp = ax.boxplot(data, widths=0.5, patch_artist=True, showfliers=False,
                medianprops=dict(color="black", lw=0.8),
                boxprops=dict(lw=0.6), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
for patch, color in zip(bp["boxes"], ["#0072B2", "#999999", "#E69F00"]):
    patch.set_facecolor(color)
    patch.set_alpha(0.55)
rng = np.random.default_rng(7)
for i, d in enumerate(data):
    sub = d if len(d) < 250 else rng.choice(d, 250, replace=False)
    ax.scatter(np.full(len(sub), i + 1) + rng.normal(0, 0.06, len(sub)), sub,
               s=1.8, color="#333333", alpha=0.25, linewidth=0)
ax.axhline(0, color="#444444", lw=0.6)
ax.set_xticks([1, 2, 3])
ax.set_xticklabels(["Actives", "Decoys", "Approved"])
ax.set(ylabel="ΔΔscore (S94A − WT, kcal/mol)", ylim=(-4.5, 4.5))
panel(ax, "c")

ax = axes[1, 1]
t = top.sort_values("score_wt").head(16)
yy = np.arange(len(t))[::-1]
ax.scatter(t.score_wt, yy, s=26, c=[TIER_C.get(x, TIER_C["-"]) for x in t.tier],
           edgecolor="white", linewidth=0.5, zorder=3)
for y, v in zip(yy, t.score_wt):
    ax.text(v + 0.03, y, f"{v:.2f}", va="center", ha="left", fontsize=4.8)
ax.set_yticks(yy)
ax.set_yticklabels([n.title()[:22] for n in t.pref_name], fontsize=5.2)
ax.set(xlabel="4DRE Vina score (kcal/mol)", xlim=(-10.2, -8.5))
panel(ax, "d", x=-0.55)
save(fig, "fig6_docking_summary")

# ---------------- Fig 7: confirmation / sensitivity ----------------
c = pd.read_csv(os.path.join(FD, "confirmation_fig.csv"))
ut = pd.read_csv(os.path.join(FD, "u0j_thresholds.csv")).iloc[0]
fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.9), layout="constrained")
ax = axes[0]
for i, col in enumerate(["4DRE_dseed", "4DTI_dseed"]):
    v = c[col].values
    ax.scatter(np.full(len(v), i + 1) + rng.normal(0, 0.05, len(v)), v, s=8, color=C["4DRE"], alpha=0.55, linewidth=0)
    ax.boxplot([v], positions=[i + 1], widths=0.45, patch_artist=True, showfliers=False,
               boxprops=dict(facecolor="white", lw=0.6), medianprops=dict(color="black", lw=0.8),
               whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
    for _, row in c[c[col] > 0.2].iterrows():
        ax.annotate(row.pref_name.title(), (i + 1, row[col]), xytext=(6, 2),
                    textcoords="offset points", fontsize=4.6)
ax.set_xticks([1, 2])
ax.set_xticklabels(["4DRE", "4DTI"])
ax.set(ylabel="Δscore vs seed 42 (kcal/mol)", ylim=(0, 0.85))
panel(ax, "a")

ax = axes[1]
colors = c.tier.fillna("-").map(TIER_C)
for _, row in c.iterrows():
    stable = bool(row.stable)
    ax.scatter(row["4DRE_ex16"], row["4U0J"], s=22, color=TIER_C.get(row.tier, TIER_C["-"]),
               edgecolor="black" if not stable else "white", linewidth=0.6 if not stable else 0.4,
               marker="o" if stable else "X", zorder=3)
for _, row in c[c.tier == "A"].iterrows():
    if "PREDNISOLONE" in row.pref_name.upper():
        ax.annotate("Prednisolone tebutate", (row["4DRE_ex16"], row["4U0J"]),
                    xytext=(-10.20, -8.55), textcoords="data", fontsize=5.2, ha="left", va="top",
                    arrowprops=dict(arrowstyle="-", lw=0.5, color="#666666"))
    else:
        ax.annotate(row.pref_name.title(), (row["4DRE_ex16"], row["4U0J"]), xytext=(4, -7),
                    textcoords="offset points", fontsize=5.4)
ax.axhline(ut["p1"], color="#D55E00", lw=0.7, ls=(0, (2, 2)))
ax.axhline(ut["p5"], color="#E69F00", lw=0.7, ls=(0, (4, 2)))
ax.text(-10.15, ut["p1"] + 0.06, "4U0J p1", fontsize=5, color="#D55E00")
ax.text(-10.15, ut["p5"] + 0.06, "4U0J p5", fontsize=5, color="#E69F00")
ax.plot([-10.2, -8.6], [-10.2, -8.6], color="#999999", lw=0.5, ls=(0, (1, 2)))
ax.set(xlabel="4DRE score (kcal/mol)", ylabel="4U0J score (kcal/mol)", xlim=(-10.3, -8.5), ylim=(-12.2, -8.4))
ax.legend(handles=[
    plt.Line2D([], [], marker="o", ls="", mfc="white", mec="#333333", label="stable"),
    plt.Line2D([], [], marker="X", ls="", color="#333333", label="unstable"),
], loc="lower left", handlelength=0.8, fontsize=5.5)
panel(ax, "b")
save(fig, "fig7_confirmation")

print("all figures written to", OUT)
