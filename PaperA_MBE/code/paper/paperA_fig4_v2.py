# -*- coding: utf-8 -*-
"""
Paper A (MBE) — Figure 4 v2: Non-circular validation of the regulation-driven class
====================================================================================
Design system per results/paperB/figure_spec_v1.md (Arial; GD #2166ac, RD #d6604d,
dual #e08214, neutral #bdbdbd; bold 11pt panel labels; PNG 300dpi + vector PDF).
Layout: 3 + 3 (top: LOO matrix, forest, evidence independence; bottom: purity-signal, HAR nesting, hCONDEL).

Panel data sources:
  a | results/paper/revision_v7/loo_full_matrix.csv
      (5 classification variants x {HAR, hCONDEL} enrichment OR heatmap, annotated
      with OR [95% CI] and Fisher p; row labels carry agreement with the full model)
  b | same file: forest plot of all 10 tests (5 variants x 2 evidence types),
      log-scale OR + 95% CI (CIs parsed from the "[lo, hi]" strings in the CSV)
  c | gene_classification_v7.csv + phase8_tissue_analysis/hcondel_gene_mapping.csv
      (HAR/hCONDEL gene-level overlap: 447 HAR-only / 29 both / 154 hCONDEL-only,
      chi2 phi; stratified RD enrichment per subset, Fisher exact)
  d | same file: purity-signal relation — RD class size vs hCONDEL OR (95% CI),
      points labelled by variant with agreement annotated
  e | results/phase8_tissue_analysis/campra_regulatory_summary.json
      (HAR-caMPRA nesting warning: RD all-HAR OR 4.79 (p = 4.2e-21) vs
      LOO-caMPRA active-HAR OR 1.84 (p = 0.012); GD 0.74 / 0.85 n.s.)
  f | results/phase8_tissue_analysis/hcondels_summary.json + effect_sizes_summary.json
      (hCONDEL main validation: RD OR 2.94 [1.92-4.51], p = 6.7e-06;
      GD OR 0.53 [0.33-0.80], enrichment n.s.)

Output: results/paper/figures_v2/Figure4_noncircular_validation.png (300 dpi) + .pdf
"""
import os
import json
import re
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
REV = os.path.join(ROOT, "results/paper/revision_v7")
TIS = os.path.join(ROOT, "results/phase8_tissue_analysis")
OUTDIR = os.path.join(ROOT, "results/paper/figures_v2")
os.makedirs(OUTDIR, exist_ok=True)

C_GD, C_RD, C_NEUT = "#2166ac", "#d6604d", "#bdbdbd"
C_HAR, C_HCONDEL = "#e08214", "#762a83"

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


def fmt_p(p):
    if p < 1e-3:
        m, e = f"{p:.1e}".split("e")
        return f"{m}×10$^{{{int(e)}}}$"
    return f"{p:.3f}".lstrip("0") if p < 0.01 else f"{p:.2f}"


def parse_ci(s):
    lo, hi = re.findall(r"[-\d.]+", s)
    return float(lo), float(hi)


# ---------------------------------------------------------------- load data
loo = pd.read_csv(os.path.join(REV, "loo_full_matrix.csv"))
with open(os.path.join(TIS, "campra_regulatory_summary.json"), encoding="utf-8") as f:
    campra = json.load(f)
with open(os.path.join(TIS, "hcondels_summary.json"), encoding="utf-8") as f:
    hcondel = json.load(f)
with open(os.path.join(TIS, "effect_sizes_summary.json"), encoding="utf-8") as f:
    eff = json.load(f)["odds_ratios"]

VARIANTS = ["full", "LOO-caMPRA", "LOO-tau", "LOO-nc", "LOO-brain"]
VLBL = {"full": "Full model", "LOO-caMPRA": "LOO-caMPRA", "LOO-tau": "LOO-tau",
        "LOO-nc": "LOO-nc-cons.", "LOO-brain": "LOO-brain-tau"}
loo = loo.set_index("variant").loc[VARIANTS].reset_index()
for ev in ["HAR", "hCONDEL"]:
    loo[[f"{ev}_lo", f"{ev}_hi"]] = loo[f"{ev}_CI"].apply(
        lambda s: pd.Series(parse_ci(s)))

fig = plt.figure(figsize=(11, 7.6))
gs = fig.add_gridspec(2, 3, left=0.065, right=0.985, top=0.93, bottom=0.085,
                      wspace=0.38, hspace=0.62,
                      width_ratios=[1.0, 1.0, 1.0], height_ratios=[1.0, 0.92])

# ================================================================ panel a
# LOO enrichment matrix heatmap (rows = variants, cols = evidence types)
axa = fig.add_subplot(gs[0, 0])
axb = fig.add_subplot(gs[0, 1])

M = loo[["HAR_OR", "hCONDEL_OR"]].to_numpy()
ima = axa.imshow(np.log2(M), cmap="YlOrRd", vmin=0, vmax=np.log2(10), aspect="auto")
for i in range(M.shape[0]):
    for j, ev in enumerate(["HAR", "hCONDEL"]):
        orv = M[i, j]
        lo_, hi_ = loo.loc[i, f"{ev}_lo"], loo.loc[i, f"{ev}_hi"]
        p_ = loo.loc[i, f"{ev}_p"]
        dark = np.log2(orv) > 1.9
        axa.text(j, i - 0.14, f"{orv:.2f}", ha="center", va="center",
                 fontsize=7.6, fontweight="bold",
                 color="white" if dark else "#333333")
        axa.text(j, i + 0.13, f"[{lo_:.2f}-{hi_:.2f}]", ha="center", va="center",
                 fontsize=6.0, color="white" if dark else "#555555")
        axa.text(j, i + 0.33, f"p={fmt_p(p_)}", ha="center", va="center",
                 fontsize=6.0, color="white" if dark else "#555555")
axa.set_xticks([0, 1])
axa.set_xticklabels(["HAR\nproximity", "hCONDEL\noverlap"], fontsize=7.4)
axa.set_yticks(range(5))
axa.set_yticklabels([f"{VLBL[v]}\n(agr. {a:.2f})" for v, a in
                     zip(loo["variant"], loo["agreement"])], fontsize=7.2)
axa.set_title("LOO enrichment matrix (OR, 95% CI, p)", loc="left", fontsize=FS_TITLE)
axa.tick_params(length=0)
for s in axa.spines.values():
    s.set_visible(False)
panel_label(axa, "a", dx=-0.30)
cb = fig.colorbar(ima, ax=axa, fraction=0.046, pad=0.03, ticks=[0, 1, 2, 3])
cb.ax.set_yticklabels(["1", "2", "4", "8"])
cb.set_label("Odds ratio", fontsize=7)
cb.ax.tick_params(labelsize=6.5)
cb.outline.set_visible(False)

# ================================================================ panel b
# Forest plot: all 10 tests (5 variants x HAR/hCONDEL)
rows = []
for i, v in enumerate(loo["variant"]):
    rows.append((v, "HAR", loo.loc[i, "HAR_OR"], loo.loc[i, "HAR_lo"],
                 loo.loc[i, "HAR_hi"], loo.loc[i, "HAR_p"]))
    rows.append((v, "hCONDEL", loo.loc[i, "hCONDEL_OR"], loo.loc[i, "hCONDEL_lo"],
                 loo.loc[i, "hCONDEL_hi"], loo.loc[i, "hCONDEL_p"]))
rows = rows[::-1]  # full model on top
ys = np.arange(len(rows)) * 1.0
for (v, ev, orv, lo_, hi_, p_), y in zip(rows, ys):
    col = C_HAR if ev == "HAR" else C_HCONDEL
    axb.plot([lo_, hi_], [y, y], color=col, lw=1.3, zorder=2)
    axb.plot([lo_, lo_], [y - 0.14, y + 0.14], color=col, lw=1.1, zorder=2)
    axb.plot([hi_, hi_], [y - 0.14, y + 0.14], color=col, lw=1.1, zorder=2)
    axb.scatter([orv], [y], s=26, color=col, zorder=3,
                marker="s" if ev == "HAR" else "o")
    axb.text(24, y, f"{orv:.2f} [{lo_:.2f}-{hi_:.2f}]", va="center", ha="left",
             fontsize=6.0, color="#444444")
axb.axvline(1.0, color="#7f7f7f", ls="--", lw=0.9, zorder=1)
axb.set_xscale("log")
axb.set_xlim(0.9, 110)
axb.set_xticks([1, 2, 4, 8, 16])
axb.set_xticklabels(["1", "2", "4", "8", "16"])
axb.set_yticks(ys)
axb.set_yticklabels([VLBL[v] if ev == "HAR" else "" for v, ev, *_ in rows],
                    fontsize=7.2)
# evidence-type legend (marker shape + color)
leg_b = [Line2D([0], [0], marker="s", color=C_HAR, lw=1.2, ms=5, label="HAR proximity"),
         Line2D([0], [0], marker="o", color=C_HCONDEL, lw=1.2, ms=5, label="hCONDEL overlap")]
axb.legend(handles=leg_b, loc="lower left", frameon=False, fontsize=6.6,
           handletextpad=0.4, labelspacing=0.3, borderaxespad=0.1, ncol=1,
           bbox_to_anchor=(0.33, 0.09))
axb.set_ylim(-0.7, len(rows) - 0.3)
axb.set_xlabel("Odds ratio (log scale), Fisher exact")
axb.set_title("All 10 LOO enrichment tests remain significant", loc="left",
              fontsize=FS_TITLE)
style_ax(axb, xgrid=True)
panel_label(axb, "b", dx=-0.135)
axb.text(0.985, -0.30, "All Holm-adjusted p < 0.05 across variants (Table 3)",
         transform=axb.transAxes, ha="right", va="top", fontsize=6.4, color="#555555")

# ================================================================ panel c
# Evidence independence: HAR vs hCONDEL overlap + stratified RD enrichment
from scipy.stats import fisher_exact as _fisher, chi2_contingency as _chi2

axc2 = fig.add_subplot(gs[0, 2])
_v7 = pd.read_csv(os.path.join(ROOT, "results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv"))
_hm = pd.read_csv(os.path.join(TIS, "hcondel_gene_mapping.csv"))
_v7["has_hcondel"] = _v7["gene_id"].isin(set(_hm["gene_id"]))
_v7["has_har"] = _v7["n_hars"] > 0
_tab = pd.crosstab(_v7["has_har"], _v7["has_hcondel"])
_chi2v, _phi_p, _, _ = _chi2(_tab)
_phi = float(np.sqrt(_chi2v / len(_v7)))

subsets = []
for nm, m in [("HAR only", _v7["has_har"] & ~_v7["has_hcondel"]),
              ("both", _v7["has_har"] & _v7["has_hcondel"]),
              ("hCONDEL only", ~_v7["has_har"] & _v7["has_hcondel"])]:
    sub = _v7[m]
    rd = int((sub["classification_v7"] == "regulation-driven").sum())
    rest = _v7[~m]
    rd_rest = int((rest["classification_v7"] == "regulation-driven").sum())
    orv, pv = _fisher([[rd, len(sub) - rd], [rd_rest, len(rest) - rd_rest]])
    subsets.append((nm, len(sub), rd, orv, pv))

ys2 = np.arange(3)[::-1]
cols2 = [C_HAR, "#7f7f7f", C_HCONDEL]
for (nm, n, rd, orv, pv), yv, col in zip(subsets, ys2, cols2):
    axc2.barh(yv, n, height=0.62, color=col, edgecolor="white", linewidth=0.5, zorder=3)
    axc2.text(n + 14, yv + 0.15, f"{n} genes; RD {rd} ({100 * rd / n:.1f}%)",
              va="center", ha="left", fontsize=6.6, color="#444444")
    axc2.text(n + 14, yv - 0.19, f"OR {orv:.2f}, p={fmt_p(pv)}",
              va="center", ha="left", fontsize=6.4, color=C_RD, fontweight="bold")
axc2.set_yticks(ys2)
axc2.set_yticklabels([s[0] for s in subsets], fontsize=7.4)
axc2.set_xlim(0, 1150)
axc2.set_xticks([0, 200, 400])
axc2.set_xticklabels(["0", "200", "400"])
axc2.set_xlabel("Genes in subset")
axc2.set_title("Evidence independence: HAR vs hCONDEL", loc="left", fontsize=FS_TITLE)
style_ax(axc2, xgrid=True)
panel_label(axc2, "c", dx=-0.30)
_n_hc = int(_tab[True].sum())
axc2.text(0.97, 0.60, f"Overlap: {subsets[1][1]}/{_n_hc} hCONDEL genes "
          f"({100 * subsets[1][1] / _n_hc:.1f}%);\nφ = {_phi:.3f}, χ² p = {fmt_p(_phi_p)}\n"
          "hCONDEL-only enrichment (bottom row)\nexcludes HAR re-detection",
          transform=axc2.transAxes, ha="right", va="top", fontsize=6.3, color="#555555")
print("panel c: phi=%.4f p=%.4g | subsets:" % (_phi, _phi_p),
      [(s[0], s[1], s[2], round(s[3], 2), "%.3g" % s[4]) for s in subsets])

# ================================================================ panel d
# Purity-signal relation: RD class size vs hCONDEL OR
axc = fig.add_subplot(gs[1, 0])
x = loo["regulation-driven"].to_numpy()
y = loo["hCONDEL_OR"].to_numpy()
ylo = loo["hCONDEL_lo"].to_numpy()
yhi = loo["hCONDEL_hi"].to_numpy()
agr = loo["agreement"].to_numpy()
for xi, yi, lo_, hi_, v, a in zip(x, y, ylo, yhi, loo["variant"], agr):
    axc.plot([xi, xi], [lo_, hi_], color=C_RD, lw=1.1, zorder=2, alpha=0.8)
    axc.scatter([xi], [yi], s=44, color=C_RD, edgecolors="white", linewidths=0.6,
                zorder=3)
    dxoff, dyoff, ha = 14, 0.10, "left"
    if v == "LOO-nc":
        dxoff, dyoff, ha = 16, -0.28, "left"
    elif v == "LOO-tau":
        dxoff, dyoff, ha = 14, 0.12, "left"
    elif v == "full":
        dxoff, dyoff, ha = 16, 0.10, "left"
    elif v == "LOO-brain":
        dxoff, dyoff, ha = 0, 0.16, "center"
    elif v == "LOO-caMPRA":
        dxoff, dyoff, ha = 0, -0.34, "center"
    axc.text(xi + dxoff, yi + dyoff, f"{VLBL[v]}\n(agr. {a:.2f})", fontsize=6.2,
             ha=ha, va="bottom" if dyoff > 0 else "top", color="#444444")
axc.axhline(1.0, color="#7f7f7f", ls="--", lw=0.9, zorder=1)
axc.set_xlabel("Regulation-driven class size (genes)")
axc.set_ylabel("hCONDEL enrichment OR")
axc.set_title("Purity-signal trade-off across variants", loc="left", fontsize=FS_TITLE)
axc.set_xlim(40, 1080)
axc.set_ylim(0.8, 7.4)
axc.set_xscale("log")
axc.set_xticks([100, 200, 400, 800])
axc.set_xticklabels(["100", "200", "400", "800"])
style_ax(axc, ygrid=True)
panel_label(axc, "d", dx=-0.22)

# ================================================================ panel e
# HAR-caMPRA nesting warning
axd = fig.add_subplot(gs[1, 1])
tests = ["All HARs\n(v7 classifier)", "Active HARs\n(LOO-caMPRA)"]
rd_or = [campra["rd_v7_all_hars"]["OR"], campra["rd_loo_campra"]["OR"]]
rd_p = [campra["rd_v7_all_hars"]["p"], campra["rd_loo_campra"]["p"]]
gd_or = [campra["gd_v7_all_hars"]["OR"], campra["gd_loo_campra"]["OR"]]
gd_p = [campra["gd_v7_all_hars"]["p"], campra["gd_loo_campra"]["p"]]
# CIs from effect_sizes where available
rd_ci = [(eff["caMPRA_activeHAR_RD_loo"]["ci95_low"], eff["caMPRA_activeHAR_RD_loo"]["ci95_high"])]
gd_ci = [(eff["caMPRA_activeHAR_GD_loo"]["ci95_low"], eff["caMPRA_activeHAR_GD_loo"]["ci95_high"])]
xx = np.arange(2)
w = 0.34
b1 = axd.bar(xx - w / 2, rd_or, width=w, color=C_RD, edgecolor="white",
             linewidth=0.5, zorder=3, label="RD")
b2 = axd.bar(xx + w / 2, gd_or, width=w, color=C_GD, edgecolor="white",
             linewidth=0.5, zorder=3, label="GD")
axd.errorbar([xx[1] - w / 2], [rd_or[1]],
             yerr=[[rd_or[1] - rd_ci[0][0]], [rd_ci[0][1] - rd_or[1]]],
             fmt="none", ecolor="#7f3b33", elinewidth=1.0, capsize=2.5, zorder=4)
axd.errorbar([xx[1] + w / 2], [gd_or[1]],
             yerr=[[gd_or[1] - gd_ci[0][0]], [gd_ci[0][1] - gd_or[1]]],
             fmt="none", ecolor="#3c5a78", elinewidth=1.0, capsize=2.5, zorder=4)
# group 0 (no CI): labels above bar top; group 1 (CI): labels above CI cap
axd.text(xx[0] - w / 2, rd_or[0] + 0.14, f"{rd_or[0]:.2f}", ha="center", fontsize=7.2,
         fontweight="bold", color=C_RD)
axd.text(xx[0] - w / 2, rd_or[0] + 0.62, f"p={fmt_p(rd_p[0])}", ha="center",
         fontsize=6.2, color="#555555")
axd.text(xx[1] - w / 2, rd_ci[0][1] + 0.12, f"{rd_or[1]:.2f}", ha="center", fontsize=7.2,
         fontweight="bold", color=C_RD)
axd.text(xx[1] - w / 2, rd_ci[0][1] + 0.58, f"p={fmt_p(rd_p[1])}", ha="center",
         fontsize=6.2, color="#555555")
axd.text(xx[0] + w / 2, gd_or[0] / 2, f"{gd_or[0]:.2f}", ha="center", va="center",
         fontsize=7.2, fontweight="bold", color="white")
axd.text(xx[0] + w / 2, gd_or[0] + 0.62, "n.s.", ha="center", fontsize=6.2, color="#555555")
axd.text(xx[1] + w / 2, gd_ci[0][1] + 0.12, f"{gd_or[1]:.2f}", ha="center", fontsize=7.2,
         fontweight="bold", color=C_GD)
axd.text(xx[1] + w / 2, gd_ci[0][1] + 0.58, "n.s.", ha="center", fontsize=6.2, color="#555555")
axd.axhline(1.0, color="#7f7f7f", ls="--", lw=0.9, zorder=1)
axd.set_xticks(xx)
axd.set_xticklabels(tests, fontsize=7.4)
axd.set_ylabel("Odds ratio (Fisher exact)")
axd.set_title("HAR-caMPRA nesting: conservative estimate", loc="left",
              fontsize=FS_TITLE)
axd.set_ylim(0, 6.3)
style_ax(axd, ygrid=True)
panel_label(axd, "e", dx=-0.20)
axd.legend(loc="upper right", frameon=False, fontsize=6.8, handletextpad=0.4,
           bbox_to_anchor=(1.0, 1.02))
axd.text(0.47, 0.82, "caMPRA-active elements are 100% nested\n"
         "within HARs; LOO-caMPRA removes this circularity",
         transform=axd.transAxes, ha="left", va="top", fontsize=6.4, color="#555555")

# ================================================================ panel f
# hCONDEL main validation
axe = fig.add_subplot(gs[1, 2])
rd_e = hcondel["rd_enrichment"]
gd_e = hcondel["gd_enrichment"]
# RD CI from loo_full_matrix (same source as panels a/b); GD CI: conditional MLE
full_row = loo.loc[loo["variant"] == "full"].iloc[0]
rd_lo, rd_hi = full_row["hCONDEL_lo"], full_row["hCONDEL_hi"]
gd_lo, gd_hi = eff["hCONDEL_GD_v7"]["ci95_low"], eff["hCONDEL_GD_v7"]["ci95_high"]
rows_e = [("RD", rd_e["OR"], rd_lo, rd_hi, rd_e["p"], C_RD,
           f"{hcondel['n_rd_with_hcondel']}/293 RD genes"),
          ("GD", gd_e["OR"], gd_lo, gd_hi, gd_e["p"], C_GD,
           f"{hcondel['n_gd_with_hcondel']}/1,214 GD genes")]
for k, (nm, orv, lo_, hi_, p_, col, cnt) in enumerate(rows_e):
    yv = 1 - k
    axe.plot([lo_, hi_], [yv, yv], color=col, lw=1.6, zorder=2)
    axe.plot([lo_, lo_], [yv - 0.12, yv + 0.12], color=col, lw=1.3, zorder=2)
    axe.plot([hi_, hi_], [yv - 0.12, yv + 0.12], color=col, lw=1.3, zorder=2)
    axe.scatter([orv], [yv], s=64, color=col, edgecolors="white", linewidths=0.7,
                zorder=3)
    sig = f"p={fmt_p(p_)}" if p_ < 0.05 else "n.s. (depleted)"
    tx, tha = (0.31, "left") if nm == "RD" else (5.7, "right")
    axe.text(tx, yv + 0.24, f"OR {orv:.2f} [{lo_:.2f}-{hi_:.2f}], {sig}",
             va="center", ha=tha, fontsize=7.0, color=col, fontweight="bold")
    axe.text(tx, yv - 0.26, f"hCONDEL overlap: {cnt}", va="center", ha=tha,
             fontsize=6.4, color="#555555")
axe.axvline(1.0, color="#7f7f7f", ls="--", lw=0.9, zorder=1)
axe.set_xscale("log")
axe.set_xlim(0.28, 6.0)
axe.set_xticks([0.5, 1, 2, 4])
axe.set_xticklabels(["0.5", "1", "2", "4"])
axe.set_yticks([1, 0])
axe.set_yticklabels(["RD", "GD"], fontsize=8)
axe.set_ylim(-0.55, 1.55)
axe.set_xlabel("hCONDEL enrichment OR (log scale)")
axe.set_title("hCONDEL main validation (v7)", loc="left", fontsize=FS_TITLE)
style_ax(axe, xgrid=True)
panel_label(axe, "f", dx=-0.20)
axe.text(0.98, 0.985, "183 genes overlap 583 hCONDELs;\nFisher exact, one-sided (greater)",
         transform=axe.transAxes, ha="right", va="top", fontsize=6.4, color="#555555")

out_png = os.path.join(OUTDIR, "Figure4_noncircular_validation.png")
out_pdf = os.path.join(OUTDIR, "Figure4_noncircular_validation.pdf")
fig.savefig(out_png, dpi=300)
fig.savefig(out_pdf)
print("saved:", out_png)
print("saved:", out_pdf)
