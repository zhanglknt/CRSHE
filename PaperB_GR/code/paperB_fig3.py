# -*- coding: utf-8 -*-
"""
Paper B — Figure 3: Single-cell resolution: regulatory selection concentrates in neurons
6 panels, 2x3. Publication grade per results/paperB/figure_spec_v1.md.

Panel data sources (all under results/phase8_tissue_analysis/):
  a | hpa_wholebody_cellclass_enrichment.csv   (15 cell classes, LOO-brain RD OR bars + GD OR dots)
  b | hpa_wholebody_celltype_enrichment.csv    (154 cell types ranked by LOO-brain RD OR)
  c | hpa_wholebody_celltype_enrichment.csv + hpa_wholebody_summary.json
      (per-type RD OR neuronal n=9 vs non-neuronal n=145 across v7 / LOO-brain / LOO-tau; MWU p)
  d | hpa_single_nuclei_enrichment.csv + hpa_single_nuclei_summary.json
      (34 brain single-nuclei types, RD OR by brain cell class; class-level OR from json)
  e | hpa_wholebody_celltype_enrichment.csv    (all 9 neuronal types by LOO-brain RD OR, FDR stars)
  f | hpa_wholebody_cellclass_enrichment.csv   (15 classes, LOO-brain GD OR bars — null contrast to a)

Output: results/paperB/figures/Figure3_single_cell.png (300 dpi) + .pdf
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
from scipy import stats

# ---------------------------------------------------------------- paths
ROOT = os.environ.get("HSD_BASE") or (
    "/mnt/d/人类正选择基因项目" if os.path.exists("/mnt/d") else "D:/人类正选择基因项目")
DATA = os.path.join(ROOT, "results/phase8_tissue_analysis")
OUTDIR = os.path.join(ROOT, "results/paperB/figures")
os.makedirs(OUTDIR, exist_ok=True)

# ---------------------------------------------------------------- design system (figure_spec_v1.md)
C_GD = "#2166ac"        # gene-driven blue
C_GD_RELAX = "#92c5de"
C_RD = "#d6604d"        # regulation-driven red
C_NEUR = "#b2182b"      # neuronal highlight dark red
C_DUAL = "#e08214"
C_NEUT = "#bdbdbd"
C_RD_NS = "#fddbc7"     # RD n.s. desaturated
C_GD_NS = "#d1e5f0"     # GD n.s. desaturated
C_GRID = "#e5e5e5"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.linewidth": 0.7,
    "axes.edgecolor": "#333333",
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "axes.labelsize": 8.5,
    "axes.titlesize": 9.5,
    "legend.fontsize": 7.5,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

FS_LABEL = 11   # panel label a/b/c...
FS_TITLE = 9.5
FS_TICK = 7.5
FS_ANNOT = 7.5


def style_ax(ax, ygrid=False, xgrid=False):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if ygrid:
        ax.grid(axis="y", color=C_GRID, linewidth=0.6, zorder=0)
    if xgrid:
        ax.grid(axis="x", color=C_GRID, linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)


def panel_label(ax, letter, dx=-0.16, dy=1.04):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=FS_LABEL,
            fontweight="bold", va="top", ha="left")


def stars(p):
    if p < 1e-3:
        return "***"
    if p < 1e-2:
        return "**"
    if p < 0.05:
        return "*"
    return "n.s."


# ---------------------------------------------------------------- load data
cc = pd.read_csv(os.path.join(DATA, "hpa_wholebody_cellclass_enrichment.csv"))
ct = pd.read_csv(os.path.join(DATA, "hpa_wholebody_celltype_enrichment.csv"))
sn = pd.read_csv(os.path.join(DATA, "hpa_single_nuclei_enrichment.csv"))
with open(os.path.join(DATA, "hpa_wholebody_summary.json"), encoding="utf-8") as f:
    wb_sum = json.load(f)
with open(os.path.join(DATA, "hpa_single_nuclei_summary.json"), encoding="utf-8") as f:
    sn_sum = json.load(f)["hpa_single_nuclei"]
sn_class_or = {r["cell_class"]: r for r in sn_sum["cell_class_results"]}

# ---------------------------------------------------------------- figure layout
fig = plt.figure(figsize=(11, 8))
gs = fig.add_gridspec(2, 3, left=0.155, right=0.99, top=0.94, bottom=0.07,
                      wspace=0.50, hspace=0.52, width_ratios=[1.12, 1.0, 1.0])

# ================================================================ panel a
# 15 cell classes (whole body): LOO-brain RD OR bars + GD OR dots
axa = fig.add_subplot(gs[0, 0])
d = cc.sort_values("loo_brain_rd_OR", ascending=True).reset_index(drop=True)
y = np.arange(len(d))
colors_a = []
for _, r in d.iterrows():
    if r["cell_type_class"] == "neuronal cells":
        colors_a.append(C_NEUR)
    elif r["loo_brain_rd_fdr"] < 0.05:
        colors_a.append(C_RD)
    else:
        colors_a.append(C_RD_NS)
axa.barh(y, d["loo_brain_rd_OR"], height=0.68, color=colors_a, zorder=2)
axa.scatter(d["loo_brain_gd_OR"], y, s=22, color=C_GD, zorder=3,
            edgecolor="white", linewidth=0.5)
axa.axvline(1.0, color="#7f7f7f", linestyle="--", linewidth=0.8, zorder=1)
axa.set_yticks(y)
axa.set_yticklabels([f"{c} ({n})" for c, n in zip(d["cell_type_class"], d["n_types"])],
                    fontsize=7)
axa.set_xlim(0, 3.15)
axa.set_xlabel("Odds ratio of high-expression enrichment")
axa.set_title("15 cell classes, whole body (LOO-brain)", loc="left", fontsize=FS_TITLE)
style_ax(axa, xgrid=True)
panel_label(axa, "a", dx=-0.62)
leg_a = [Patch(facecolor=C_NEUR, label="RD neuronal"),
         Patch(facecolor=C_RD, label="RD FDR<0.05"),
         Patch(facecolor=C_RD_NS, label="RD n.s."),
         Line2D([0], [0], marker="o", color="none", markerfacecolor=C_GD,
                markeredgecolor="white", markersize=5.5, label="GD OR")]
axa.legend(handles=leg_a, loc="lower right", frameon=False, borderaxespad=0.1,
           handletextpad=0.35, labelspacing=0.32, fontsize=6.9)

# ================================================================ panel b
# all 154 cell types ranked by LOO-brain RD OR
axb = fig.add_subplot(gs[0, 1])
r = ct.sort_values("loo_brain_rd_OR", ascending=False).reset_index(drop=True)
rank = np.arange(1, len(r) + 1)
is_neu = r["cell_type_class"] == "neuronal cells"
is_sig = r["loo_brain_rd_fdr"] < 0.05
col_b = np.where(is_neu, C_NEUR, np.where(is_sig, C_RD, C_RD_NS))
axb.axhline(1.0, color="#7f7f7f", linestyle="--", linewidth=0.8, zorder=1)
axb.scatter(rank[~is_neu], r.loc[~is_neu, "loo_brain_rd_OR"], s=9,
            c=col_b[~is_neu.values], zorder=2, linewidth=0)
axb.scatter(rank[is_neu], r.loc[is_neu, "loo_brain_rd_OR"], s=16,
            c=C_NEUR, zorder=3, linewidth=0)
# top-5 labels: stacked list in the free upper-middle area (curve is below)
top5 = r.head(5)
axb.text(10, 3.52, "Top 5:", fontsize=7.2, fontweight="bold", color="#333333",
         ha="left", va="top")
for i, (_, row) in enumerate(top5.iterrows()):
    axb.text(10, 3.34 - i * 0.155, f"{row['cell_type']}  {row['loo_brain_rd_OR']:.2f}",
             fontsize=6.9, ha="left", va="top", color=C_NEUR)
axb.set_xlim(-4, len(r) + 4)
axb.set_ylim(0.5, 3.65)
axb.set_xlabel("Cell-type rank (n = 154)")
axb.set_ylabel("RD odds ratio (LOO-brain)")
axb.set_title("All 154 cell types ranked by RD OR", loc="left", fontsize=FS_TITLE)
style_ax(axb, ygrid=True)
panel_label(axb, "b", dx=-0.20)
leg_b = [Line2D([0], [0], marker="o", color="none", markerfacecolor=C_NEUR,
                markersize=5.5, label="Neuronal (n = 9)"),
         Line2D([0], [0], marker="o", color="none", markerfacecolor=C_RD,
                markersize=5, label="Non-neuronal, FDR<0.05"),
         Line2D([0], [0], marker="o", color="none", markerfacecolor=C_RD_NS,
                markersize=5, label="Non-neuronal, n.s.")]
axb.legend(handles=leg_b, loc="lower left", frameon=False, borderaxespad=0.1,
           handletextpad=0.3, labelspacing=0.35, bbox_to_anchor=(0.60, 0.40))

# ================================================================ panel c
# neuronal vs non-neuronal per-type RD OR across 5 classification variants
axc = fig.add_subplot(gs[0, 2])
variants = [("v7", "v7_rd_OR", "Full\n(v7)"),
            ("loo_brain", "loo_brain_rd_OR", "LOO-\nbrain"),
            ("loo_tau", "loo_tau_rd_OR", "LOO-\ntau"),
            ("loo_campra", "loo_campra_rd_OR", "LOO-\ncaMPRA"),
            ("loo_double", "loo_double_rd_OR", "Double-\nLOO")]
neu_mask = ct["cell_type_class"] == "neuronal cells"
bw = 0.34
rng = np.random.default_rng(7)
all_max, all_min = -np.inf, np.inf
for _, col, _ in variants:
    v = ct[col].to_numpy()
    all_max = max(all_max, np.nanmax(v))
    all_min = min(all_min, np.nanmin(v[v > 0]))
p_texts = []
for gi, (key, col, lab) in enumerate(variants):
    neu_v = ct.loc[neu_mask, col].to_numpy()
    non_v = ct.loc[~neu_mask, col].to_numpy()
    p2 = stats.mannwhitneyu(neu_v, non_v, alternative="two-sided").pvalue
    p_texts.append((gi, p2))
    for off, vals, fc in [(-bw / 2 - 0.02, non_v, "#d9d9d9"),
                          (bw / 2 + 0.02, neu_v, C_NEUR)]:
        bp = axc.boxplot([np.log2(vals)], positions=[gi + off], widths=bw * 0.78,
                         patch_artist=True, showfliers=False, medianprops=dict(color="black", linewidth=1),
                         whiskerprops=dict(linewidth=0.7), capprops=dict(linewidth=0.7),
                         boxprops=dict(linewidth=0.7))
        bp["boxes"][0].set_facecolor(fc)
        bp["boxes"][0].set_alpha(0.9 if fc == C_NEUR else 1.0)
        jitter = rng.uniform(-bw * 0.2, bw * 0.2, size=len(vals))
        axc.scatter(np.full(len(vals), gi + off) + jitter, np.log2(vals), s=4.5,
                    color=fc, edgecolor="none", alpha=0.55 if fc != C_NEUR else 0.9, zorder=3)
for gi, p2 in p_texts:
    ex = int(np.floor(np.log10(p2)))
    mant = p2 / 10.0 ** ex
    axc.text(gi, 3.82, f"{mant:.1f}$\\times$10$^{{{ex}}}$", ha="center", va="top",
             fontsize=6.3, color="#333333")
axc.text(-0.54, 3.82, "P:", ha="left", va="top", fontsize=6.3, color="#333333")
axc.axhline(0, color="#7f7f7f", linestyle="--", linewidth=0.8, zorder=1)
axc.set_xticks(range(len(variants)))
axc.set_xticklabels([v[2] for v in variants], fontsize=6.9)
yticks = [0, 1, 2, 3]
axc.set_yticks(yticks)
axc.set_yticklabels(["1", "2", "4", "8"])
axc.set_ylim(max(-1.35, np.log2(all_min) - 0.25), 3.9)
axc.set_xlim(-0.55, len(variants) - 0.45)
axc.set_ylabel("RD odds ratio (log$_2$ scale)")
axc.set_title("Neuronal vs non-neuronal RD OR", loc="left", fontsize=FS_TITLE)
style_ax(axc, ygrid=True)
panel_label(axc, "c", dx=-0.20)
leg_c = [Patch(facecolor=C_NEUR, label="Neuronal (n = 9)"),
         Patch(facecolor="#d9d9d9", label="Non-neuronal (n = 145)")]
axc.legend(handles=leg_c, loc="upper center", frameon=False, borderaxespad=0.1,
           handletextpad=0.4, labelspacing=0.35, ncol=2, fontsize=6.8,
           bbox_to_anchor=(0.5, 0.925), columnspacing=1.2)

# ================================================================ panel d
# brain single-nuclei detail: 34 types, RD OR colored by cell class
axd = fig.add_subplot(gs[1, 0])
class_colors = {
    "Neuronal cells": C_NEUR,
    "Glial cells": C_RD,
    "Ciliated cells": C_DUAL,
    "Endothelial and mural cells": C_GD_RELAX,
    "Mesenchymal cells": C_NEUT,
    "Blood and immune cells": C_GD_NS,
}
ds = sn.copy()
ds["class_rank"] = ds["cell_class"].map(
    {k: i for i, k in enumerate(["Neuronal cells", "Glial cells", "Ciliated cells",
                                 "Endothelial and mural cells", "Mesenchymal cells",
                                 "Blood and immune cells"])})
ds = ds.sort_values(["class_rank", "rd_OR"], ascending=[True, False]).reset_index(drop=True)
yd = np.arange(len(ds))[::-1]
for cls, sub in ds.groupby("cell_class", sort=False):
    yy = yd[sub.index.to_numpy()]
    axd.scatter(sub["rd_OR"], yy, s=20, color=class_colors[cls], zorder=3,
                edgecolor="white", linewidth=0.4)
axd.axvline(1.0, color="#7f7f7f", linestyle="--", linewidth=0.8, zorder=1)
axd.set_yticks(yd)
axd.set_yticklabels(ds["cell_type"], fontsize=6.2)
axd.set_xlim(0.4, 3.9)
axd.set_xlabel("RD odds ratio (LOO-brain)")
axd.set_title("Brain single-nuclei, 34 types", loc="left", fontsize=FS_TITLE)
style_ax(axd, xgrid=True)
panel_label(axd, "d", dx=-0.78)
leg_d = []
for cls in ["Neuronal cells", "Glial cells", "Ciliated cells",
            "Endothelial and mural cells", "Mesenchymal cells", "Blood and immune cells"]:
    rec = sn_class_or[cls]
    short = cls.replace(" cells", "").replace("Endothelial and mural", "Endothelial/mural")
    leg_d.append(Line2D([0], [0], marker="o", color="none", markerfacecolor=class_colors[cls],
                        markersize=5, label=f"{short} (class OR {rec['rd_OR']:.2f})"))
axd.legend(handles=leg_d, loc="upper left", frameon=True, borderaxespad=0.1,
           handletextpad=0.3, labelspacing=0.3, fontsize=6.6,
           facecolor="white", framealpha=0.9, edgecolor="none")

# ================================================================ panel e
# neuronal subtype breakdown: all 9 neuronal types (whole-body file)
axe = fig.add_subplot(gs[1, 1])
neu = (ct[ct["cell_type_class"] == "neuronal cells"]
       .sort_values("loo_brain_rd_OR", ascending=True).reset_index(drop=True))
ye = np.arange(len(neu))
cols_e = [C_NEUR if f < 0.05 else C_RD_NS for f in neu["loo_brain_rd_fdr"]]
axe.barh(ye, neu["loo_brain_rd_OR"], height=0.62, color=cols_e, zorder=2)
for i, (_, row) in enumerate(neu.iterrows()):
    axe.text(row["loo_brain_rd_OR"] + 0.05, i,
             f"{row['loo_brain_rd_OR']:.2f} {stars(row['loo_brain_rd_fdr'])}",
             va="center", ha="left", fontsize=7)
axe.axvline(1.0, color="#7f7f7f", linestyle="--", linewidth=0.8, zorder=1)
axe.set_yticks(ye)
axe.set_yticklabels(neu["cell_type"], fontsize=7)
axe.set_xlim(0, 4.1)
axe.set_xlabel("RD odds ratio (LOO-brain)")
axe.set_title("All nine neuronal cell types, whole body", loc="left", fontsize=FS_TITLE)
style_ax(axe, xgrid=True)
panel_label(axe, "e", dx=-0.42)
axe.text(0.985, 0.03, "*** FDR<0.001", transform=axe.transAxes, ha="right",
         va="bottom", fontsize=7, color="#555555")

# ================================================================ panel f
# GD has no cell-class preference: 15-class GD OR bars (null contrast to a)
axf = fig.add_subplot(gs[1, 2])
df_ = cc.sort_values("loo_brain_gd_OR", ascending=True).reset_index(drop=True)
yf = np.arange(len(df_))
cols_f = [C_GD if f < 0.05 else C_GD_NS for f in df_["loo_brain_gd_fdr"]]
axf.barh(yf, df_["loo_brain_gd_OR"], height=0.68, color=cols_f, zorder=2)
axf.scatter(df_["loo_brain_rd_OR"], yf, s=22, color=C_RD, zorder=3,
            edgecolor="white", linewidth=0.5)
axf.axvline(1.0, color="#7f7f7f", linestyle="--", linewidth=0.8, zorder=1)
axf.set_yticks(yf)
axf.set_yticklabels(df_["cell_type_class"], fontsize=7)
axf.set_xlim(0, 3.15)
axf.set_xlabel("Odds ratio of high-expression enrichment")
axf.set_title("GD flat across the same 15 classes", loc="left", fontsize=FS_TITLE)
style_ax(axf, xgrid=True)
panel_label(axf, "f", dx=-0.42)
leg_f = [Patch(facecolor=C_GD, label="GD FDR<0.05"),
         Patch(facecolor=C_GD_NS, label="GD n.s."),
         Line2D([0], [0], marker="o", color="none", markerfacecolor=C_RD,
                markeredgecolor="white", markersize=5.5, label="RD OR")]
axf.legend(handles=leg_f, loc="upper right", frameon=False, borderaxespad=0.1,
           handletextpad=0.35, labelspacing=0.32, fontsize=6.9)

# ---------------------------------------------------------------- save
png = os.path.join(OUTDIR, "Figure3_single_cell.png")
pdf = os.path.join(OUTDIR, "Figure3_single_cell.pdf")
fig.savefig(png, dpi=300)
fig.savefig(pdf, bbox_inches="tight")
print("saved:", png)
print("saved:", pdf)

# ---------------------------------------------------------------- verification summary
print("\n=== Figure 3 per-panel summary ===")
print(f"a | hpa_wholebody_cellclass_enrichment.csv | 15 classes | RD OR range "
      f"{cc['loo_brain_rd_OR'].min():.2f}-{cc['loo_brain_rd_OR'].max():.2f}; "
      f"neuronal RD OR {cc.loc[cc['cell_type_class']=='neuronal cells','loo_brain_rd_OR'].iloc[0]:.2f} "
      f"(FDR {cc.loc[cc['cell_type_class']=='neuronal cells','loo_brain_rd_fdr'].iloc[0]:.1e})")
print(f"b | hpa_wholebody_celltype_enrichment.csv | 154 types | top5: " +
      "; ".join(f"{r_['cell_type']} {r_['loo_brain_rd_OR']:.2f}" for _, r_ in r.head(5).iterrows()))
p2_all = {k: stats.mannwhitneyu(ct.loc[neu_mask, c], ct.loc[~neu_mask, c],
                                alternative="two-sided").pvalue for k, c, _ in variants}
print(f"c | neuronal n={neu_mask.sum()} vs non-neuronal n={(~neu_mask).sum()} | two-sided MWU p: " +
      ", ".join(f"{k} {p:.1e}" for k, p in p2_all.items()))
print(f"d | hpa_single_nuclei_enrichment.csv | {len(sn)} types, 6 classes | "
      f"class OR neuronal {sn_class_or['Neuronal cells']['rd_OR']:.2f} vs glial {sn_class_or['Glial cells']['rd_OR']:.2f}")
print(f"e | 9 neuronal types | RD OR {neu['loo_brain_rd_OR'].min():.2f}-{neu['loo_brain_rd_OR'].max():.2f}, "
      f"all FDR<0.05: {(neu['loo_brain_rd_fdr']<0.05).all()}")
print(f"f | GD OR range {df_['loo_brain_gd_OR'].min():.2f}-{df_['loo_brain_gd_OR'].max():.2f} "
      f"(flat, null contrast to panel a)")
