# -*- coding: utf-8 -*-
"""
Paper A (MBE) — Figure 3 v2: Classification landscape and robustness (6 panels, 2x3)
===================================================================================
Design system per results/paperB/figure_spec_v1.md (Arial; GD #2166ac, RD #d6604d,
dual #e08214, neutral #bdbdbd; bold 11pt panel labels; PNG 300dpi + vector PDF).

Panel data sources:
  a | results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv
      (4,974 genes; gds_v7 x rds_v7 scatter, +-0.15 neutral band, hexbin density backdrop)
  b | same file: GDS component scores (gds_p_pct_resid / gds_lrt_pct_resid / gds_relax_pct /
      gds_selectome) by classification_v7 class (5 classes x 4 components boxplots)
  c | same file: RDS component scores (rds_doan / rds_tau / rds_nc / rds_brain) by class
  d | results/paper/sensitivity/sensitivity_summary.json threshold_grid 0.05-0.30
      (class sizes vs threshold; baseline 0.15; joint-1,000-perturbation agreement 94.6%)
  e | results/paper/sensitivity/per_gene_stability.csv
      (per-gene stability across 1,000 joint weight perturbations; GD median 1.000 / RD median 0.833)
  f | results/paper/revision_v7/pca_kmeans_crosstab.csv + kmeans_multistart.json
      (unsupervised PCA/k-means vs supervised classification; RD recovery 87.4% +- 5.1%)

Output: results/paper/figures_v2/Figure3_classification_robustness.png (300 dpi) + .pdf
"""
import os
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

if os.path.exists("/mnt/d"):
    ROOT = "/mnt/d/人类正选择基因项目"
else:
    ROOT = "D:/人类正选择基因项目"
V7 = os.path.join(ROOT, "results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv")
SENS = os.path.join(ROOT, "results/paper/sensitivity")
REV = os.path.join(ROOT, "results/paper/revision_v7")
OUTDIR = os.path.join(ROOT, "results/paper/figures_v2")
os.makedirs(OUTDIR, exist_ok=True)

C_GD, C_GD_RELAX, C_RD, C_DUAL, C_NEUT = "#2166ac", "#92c5de", "#d6604d", "#e08214", "#bdbdbd"
CLS_ORDER = ["gene-driven", "gene-driven (relaxed)", "dual-driven", "regulation-driven", "neutral"]
CLS_COLOR = {"gene-driven": C_GD, "gene-driven (relaxed)": C_GD_RELAX, "dual-driven": C_DUAL,
             "regulation-driven": C_RD, "neutral": C_NEUT}
CLS_SHORT = {"gene-driven": "GD", "gene-driven (relaxed)": "GD-rel.", "dual-driven": "dual",
             "regulation-driven": "RD", "neutral": "neutral"}

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "axes.linewidth": 0.7, "axes.edgecolor": "#333333",
    "xtick.direction": "out", "ytick.direction": "out",
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "axes.labelsize": 8.5,
    "axes.titlesize": 9.5, "legend.fontsize": 7.5,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})
FS_LABEL, FS_TITLE = 11, 9.5


def style_ax(ax, ygrid=False, xgrid=False):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if ygrid:
        ax.grid(axis="y", color="#e5e5e5", linewidth=0.6, zorder=0)
    if xgrid:
        ax.grid(axis="x", color="#e5e5e5", linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)


def panel_label(ax, letter, dx=-0.16, dy=1.04):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=FS_LABEL,
            fontweight="bold", va="top", ha="left")


# ---------------------------------------------------------------- load data
gc = pd.read_csv(V7)
with open(os.path.join(SENS, "sensitivity_summary.json"), encoding="utf-8") as f:
    sens = json.load(f)
stab = pd.read_csv(os.path.join(SENS, "per_gene_stability.csv"))
xtab = pd.read_csv(os.path.join(REV, "pca_kmeans_crosstab.csv"))
with open(os.path.join(REV, "kmeans_multistart.json"), encoding="utf-8") as f:
    km = json.load(f)

gds = gc["gds_v7"].to_numpy()
rds = gc["rds_v7"].to_numpy()
cls = gc["classification_v7"].to_numpy()

fig = plt.figure(figsize=(11, 8))
gs = fig.add_gridspec(2, 3, left=0.06, right=0.99, top=0.94, bottom=0.075,
                      wspace=0.34, hspace=0.5, width_ratios=[1.12, 1.06, 1.0])

# ================================================================ panel a
axa = fig.add_subplot(gs[0, 0])
axa.hexbin(gds, rds, gridsize=28, cmap="Greys", mincnt=1, alpha=0.45,
           linewidths=0.1, edgecolors="#dddddd", zorder=1)
xx = np.linspace(0, 1, 50)
axa.fill_between(xx, xx - 0.15, xx + 0.15, color="#f5d8d8", alpha=0.35, zorder=0)
axa.text(0.60, 0.475, "neutral band (|$\\Delta$| < 0.15)", rotation=45, fontsize=6.5,
         color="#999999", ha="center", va="center", zorder=2)
axa.plot(xx, xx + 0.15, color="#999999", lw=0.9, ls="--", zorder=1)
axa.plot(xx, xx - 0.15, color="#999999", lw=0.9, ls="--", zorder=1)
axa.plot(xx, xx, color="#666666", lw=0.7, ls=":", zorder=1)
for c in ["dual-driven", "regulation-driven", "gene-driven (relaxed)", "gene-driven"]:
    m = cls == c
    axa.scatter(gds[m], rds[m], s=8, alpha=0.75, c=CLS_COLOR[c],
                edgecolors="white", linewidths=0.15, zorder=3,
                label=f"{CLS_SHORT[c]} (n = {m.sum():,})")
axa.set_xlim(-0.02, 1.02)
axa.set_ylim(-0.02, 1.02)
axa.set_xlabel("GDS v7 (coding selection score)")
axa.set_ylabel("RDS v7 (regulatory score)")
axa.set_title("Classification landscape, 4,974 genes", loc="left", fontsize=FS_TITLE)
style_ax(axa)
panel_label(axa, "a", dx=-0.075)
handles_a, labels_a = axa.get_legend_handles_labels()
handles_a = [Patch(facecolor="#cccccc", edgecolor="#aaaaaa", label="neutral (n = 3,404)")] + handles_a
lega = axa.legend(handles=handles_a, loc="upper right", frameon=True, fontsize=6.6,
                  handletextpad=0.4, labelspacing=0.32, borderaxespad=0.1)
lega.get_frame().set_facecolor("white")
lega.get_frame().set_alpha(0.85)
lega.get_frame().set_edgecolor("none")
lega.set_zorder(5)

# ================================================================ panels b/c
def component_boxplot(ax, comps, comp_labels, title, letter, dx):
    ncomp = len(comps)
    bw = 0.15
    rng = np.random.default_rng(11)
    for ci, comp in enumerate(comps):
        for ki, c in enumerate(CLS_ORDER):
            vals = gc.loc[gc["classification_v7"] == c, comp].to_numpy()
            off = (ki - 2) * bw
            bp = ax.boxplot([vals], positions=[ci + off], widths=bw * 0.86,
                            patch_artist=True, showfliers=False,
                            medianprops=dict(color="black", linewidth=0.8),
                            whiskerprops=dict(linewidth=0.6), capprops=dict(linewidth=0.6),
                            boxprops=dict(linewidth=0.6))
            bp["boxes"][0].set_facecolor(CLS_COLOR[c])
            bp["boxes"][0].set_alpha(0.9 if c != "neutral" else 1.0)
            jit = rng.uniform(-bw * 0.25, bw * 0.25, size=min(len(vals), 120))
            vs = rng.choice(vals, size=min(len(vals), 120), replace=False)
            ax.scatter(np.full(len(vs), ci + off) + jit, vs, s=2.2,
                       color=CLS_COLOR[c], alpha=0.35, edgecolors="none", zorder=3)
    ax.set_xticks(range(ncomp))
    ax.set_xticklabels(comp_labels, fontsize=7.2)
    ax.set_ylim(-0.05, 1.05)
    ax.set_ylabel("Component score (0-1)")
    ax.set_title(title, loc="left", fontsize=FS_TITLE)
    style_ax(ax, ygrid=True)
    panel_label(ax, letter, dx=dx)


axb = fig.add_subplot(gs[0, 1])
component_boxplot(axb, ["gds_p_pct_resid", "gds_lrt_pct_resid", "gds_relax_pct", "gds_selectome"],
                  ["P-value\n(resid.)", "LRT\n(resid.)", "RELAX K", "Selectome"],
                  "GDS components by class", "b", dx=-0.135)

axc = fig.add_subplot(gs[0, 2])
component_boxplot(axc, ["rds_doan", "rds_tau", "rds_nc", "rds_brain"],
                  ["caMPRA", "tau", "nc-cons.", "brain tau"],
                  "RDS components by class", "c", dx=-0.135)

# shared legend for b/c (inside panel b, upper right over the near-zero Selectome group)
leg_bc = [Patch(facecolor=CLS_COLOR[c], label=CLS_SHORT[c]) for c in CLS_ORDER]
axb.legend(handles=leg_bc, loc="upper right", frameon=False, fontsize=6.4,
           handletextpad=0.35, labelspacing=0.28, borderaxespad=0.1, ncol=1,
           bbox_to_anchor=(1.0, 1.0))

# ================================================================ panel d
axd = fig.add_subplot(gs[1, 0])
tg = sens["threshold_grid"]
thr = [r["threshold"] for r in tg]
gd_n = [r["gene-driven"] + r["gene-driven (relaxed)"] for r in tg]
rd_n = [r["regulation-driven"] for r in tg]
neu_n = [r["neutral"] for r in tg]
axd.plot(thr, gd_n, "o-", color=C_GD, lw=1.4, ms=4.5, label="GD (incl. relaxed)")
axd.plot(thr, rd_n, "o-", color=C_RD, lw=1.4, ms=4.5, label="RD")
axd.plot(thr, neu_n, "o-", color=C_NEUT, lw=1.2, ms=4, label="neutral")
axd.axvline(0.15, color="#7f7f7f", ls="--", lw=0.8)
axd.text(0.153, 4360, "baseline $\\Delta$ = 0.15", fontsize=7, color="#555555", va="top", ha="left")
for t, g, r in zip(thr, gd_n, rd_n):
    axd.annotate(f"{r}", xy=(t, r), xytext=(-2, 7 if r < 250 else -11),
                 textcoords="offset points", fontsize=6.2, ha="right", color=C_RD)
    axd.annotate(f"{g}", xy=(t, g), xytext=(-2, 7), textcoords="offset points",
                 fontsize=6.2, ha="right", color=C_GD)
axd.set_xlabel("Classification threshold $\\Delta$")
axd.set_ylabel("Class size (genes)")
axd.set_title("Threshold sensitivity (0.05-0.30)", loc="left", fontsize=FS_TITLE)
axd.set_ylim(0, 4600)
axd.set_xlim(0.035, 0.31)
style_ax(axd, ygrid=True)
panel_label(axd, "d", dx=-0.075)
jp = sens["joint_perturbation_1000"]
axd.text(0.97, 0.60, f"1,000 joint weight perturbations:\n"
                     f"mean agreement with baseline = {jp['agreement_mean']*100:.1f}%",
         transform=axd.transAxes, ha="right", va="top", fontsize=7.2, color="#555555")
axd.legend(loc="upper left", frameon=False, fontsize=6.8, handletextpad=0.4,
           labelspacing=0.3, borderaxespad=0.1, bbox_to_anchor=(0.02, 1.0))

# ================================================================ panel e
axe = fig.add_subplot(gs[1, 1])
bins = np.linspace(0, 1, 41)
for c, col, lw in [("neutral", C_NEUT, 1.0), ("gene-driven", C_GD, 1.4),
                   ("regulation-driven", C_RD, 1.4)]:
    v = stab.loc[stab["classification_v7"] == c, "stability"].to_numpy()
    axe.hist(v, bins=bins, color=col, alpha=0.55, edgecolor="white", linewidth=0.3,
             label=f"{CLS_SHORT[c]} (n = {len(v):,})", zorder=2)
gd_med = stab.loc[stab["classification_v7"] == "gene-driven", "stability"].median()
rd_med = stab.loc[stab["classification_v7"] == "regulation-driven", "stability"].median()
axe.axvline(gd_med, color=C_GD, ls="--", lw=1.2, zorder=3)
axe.axvline(rd_med, color=C_RD, ls="--", lw=1.2, zorder=3)
axe.annotate(f"GD median\n{gd_med:.3f}", xy=(gd_med, 0.97), xycoords=("data", "axes fraction"),
             xytext=(-4, -2), textcoords="offset points", ha="right", va="top",
             fontsize=7, color=C_GD)
axe.annotate(f"RD median\n{rd_med:.3f}", xy=(rd_med, 0.72), xycoords=("data", "axes fraction"),
             xytext=(4, -2), textcoords="offset points", ha="left", va="top",
             fontsize=7, color=C_RD)
axe.set_xlabel("Per-gene stability across 1,000 perturbations")
axe.set_ylabel("Genes")
axe.set_title("Per-gene classification stability", loc="left", fontsize=FS_TITLE)
style_ax(axe, ygrid=True)
panel_label(axe, "e", dx=-0.135)
axe.legend(loc="upper left", frameon=False, fontsize=6.8, handletextpad=0.4,
           labelspacing=0.3, borderaxespad=0.1, bbox_to_anchor=(0.02, 1.0))

# ================================================================ panel f
axf = fig.add_subplot(gs[1, 2])
rows = xtab.set_index("row_0").loc[CLS_ORDER, ["0", "1", "2", "3"]]
M = rows.to_numpy()
Mnorm = M / M.sum(axis=1, keepdims=True) * 100
im = axf.imshow(Mnorm, cmap="Greys", vmin=0, vmax=100, aspect="auto")
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        dark = Mnorm[i, j] > 55
        axf.text(j, i, f"{M[i, j]}\n({Mnorm[i, j]:.0f}%)", ha="center", va="center",
                 fontsize=6.6, color="white" if dark else "#333333")
# highlight the RD row (best k-means cluster 0)
axf.add_patch(plt.Rectangle((-0.5, 2.5), 4, 1, fill=False, edgecolor=C_RD, lw=1.6))
axf.set_xticks(range(4))
axf.set_xticklabels([f"k-means {k}" for k in range(4)], fontsize=7.2)
axf.set_yticks(range(5))
axf.set_yticklabels([CLS_SHORT[c] for c in CLS_ORDER], fontsize=7.2)
axf.set_title("Unsupervised k-means vs supervised classes", loc="left", fontsize=FS_TITLE)
axf.set_xlabel("PCA + k-means (k = 4) cluster")
style_ax(axf)
for s in axf.spines.values():
    s.set_visible(False)
panel_label(axf, "f", dx=-0.135)
axf.text(0.985, -0.30, f"100 random seeds: RD recovery {km['rd_recovery_mean']*100:.1f}% "
                       f"$\\pm$ {km['rd_recovery_sd']*100:.1f}% (min {km['rd_recovery_min']*100:.0f}%)",
         transform=axf.transAxes, ha="right", va="top", fontsize=6.8, color="#555555",
         clip_on=False)

# ---------------------------------------------------------------- save
png = os.path.join(OUTDIR, "Figure3_classification_robustness.png")
pdf = os.path.join(OUTDIR, "Figure3_classification_robustness.pdf")
fig.savefig(png, dpi=300)
fig.savefig(pdf, bbox_inches="tight")
print("saved:", png)
print("saved:", pdf)

# ---------------------------------------------------------------- summary
print("\n=== Paper A Figure 3 v2 per-panel summary ===")
print("a | gene_classification_v7.csv | 4,974 genes | GD 1,214+52 / RD 293 / dual 11 / neutral 3,404")
print(f"b | GDS components | medians GD vs RD: " + "; ".join(
    f"{c.split('_',1)[1]} {gc.loc[gc['classification_v7']=='gene-driven',c].median():.2f}/"
    f"{gc.loc[gc['classification_v7']=='regulation-driven',c].median():.2f}"
    for c in ["gds_p_pct_resid", "gds_lrt_pct_resid", "gds_relax_pct", "gds_selectome"]))
print(f"c | RDS components | medians GD vs RD: " + "; ".join(
    f"{c.split('_',1)[1]} {gc.loc[gc['classification_v7']=='gene-driven',c].median():.2f}/"
    f"{gc.loc[gc['classification_v7']=='regulation-driven',c].median():.2f}"
    for c in ["rds_doan", "rds_tau", "rds_nc", "rds_brain"]))
print(f"d | threshold grid 0.05-0.30 | GD 1,382->918 / RD 772->55 | joint-1,000 agreement "
      f"{jp['agreement_mean']*100:.1f}% (SD {jp['agreement_sd']*100:.1f}%)")
print(f"e | per_gene_stability.csv | GD median {gd_med:.3f}, RD median {rd_med:.3f}, "
      f"stable core >=0.8: {sens['per_gene_stability']['stable_core_ge_0.8_fraction']*100:.1f}%")
print(f"f | pca_kmeans_crosstab | RD 259/293 in cluster 0 ({259/293*100:.0f}%) | "
      f"multistart recovery {km['rd_recovery_mean']*100:.1f}% +- {km['rd_recovery_sd']*100:.1f}%")
