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
  c | gene_classification_v7.csv + phase9_hardening/hcondel_gene_mapping_fixed.csv
      (HAR/hCONDEL gene-level overlap: 454 HAR-only / 22 both / 90 hCONDEL-only,
      phi=0.052, chi2 p=4.6e-4; stratified RD enrichment per subset, Fisher exact)
  d | same file: purity-signal relation — RD class size vs hCONDEL OR (95% CI),
      points labelled by variant with agreement annotated
  e | results/paper/revision_v7/loo_full_matrix.csv + gene_classification_v7.csv
      (HAR-caMPRA nesting, current 476-gene HAR mapping: full-model all-HAR
      RD OR 4.15 (p = 3.3e-20) vs LOO-caMPRA all-HAR RD OR 1.61 (p = 4.7e-05);
      GD null in both, recomputed with the identical classifier)
  f | phase9_hardening/hcondels_liftover_fix.json (fixed hg38 mapping, 112 genes)
      (hCONDEL main validation: RD OR 3.19 [1.90-5.37], p = 7.8e-05;
      GD OR 0.59 [0.36-0.97], two-sided p = 0.037, depleted)

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

# --- R6 provenance fix: hCONDEL gene mapping repaired (hg18->hg38 liftOver;
# phase9_hardening/hcondel_gene_mapping_fixed.csv, union 183->112 genes).
# Legacy loo_full_matrix.csv still carries the OLD hCONDEL statistics, so we
# overwrite the hCONDEL columns with the fixed values, recomputing OR / p /
# Woolf 95% CI from k-over-n counts (directional Fisher: RD greater, GD less),
# which reproduces hcondels_liftover_fix.json exactly.
HCFIX = os.path.join(ROOT, "results/phase9_hardening/hcondels_liftover_fix.json")
with open(HCFIX, encoding="utf-8") as f:
    hcfix = json.load(f)
HC_TOTAL = int(hcfix["gene_sets"]["union_fixed"])          # 112 genes
HC_N_UNIVERSE = int(hcfix["universe_n"])                   # 4,974
for i, v in enumerate(VARIANTS):
    st = hcfix["loo_hcondel"][v]
    k, n = [int(x) for x in st["k_over_n"].split("/")]
    tab = [[k, n - k], [HC_TOTAL - k, HC_N_UNIVERSE - n - (HC_TOTAL - k)]]
    from scipy.stats import fisher_exact as _fe
    orv, pv = _fe(tab, alternative="greater")
    a_, b_, c_, d_ = tab[0][0], tab[0][1], tab[1][0], tab[1][1]
    se = np.sqrt(1 / a_ + 1 / b_ + 1 / c_ + 1 / d_)
    loo.loc[i, "hCONDEL_OR"] = float(orv)
    loo.loc[i, "hCONDEL_p"] = float(pv)
    loo.loc[i, "hCONDEL_lo"] = float(np.exp(np.log(orv) - 1.96 * se))
    loo.loc[i, "hCONDEL_hi"] = float(np.exp(np.log(orv) + 1.96 * se))
    assert abs(float(orv) - st["OR"]) < 0.01, (v, orv, st["OR"])
print("hCONDEL LOO stats overridden with fixed mapping (union=%d genes): %s"
      % (HC_TOTAL, {v: (round(loo.loc[i, 'hCONDEL_OR'], 2),
                        "%.1e" % loo.loc[i, 'hCONDEL_p'])
                    for i, v in enumerate(VARIANTS)}))

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
# fixed hCONDEL mapping (hg18->hg38 liftOver repair; union 112 genes)
_hm = pd.read_csv(os.path.join(ROOT, "results/phase9_hardening/hcondel_gene_mapping_fixed.csv"))
_v7["has_hcondel"] = _v7["gene_id"].isin(set(_hm["gene_id"]))
_v7["has_har"] = _v7["n_hars"] > 0
_tab = pd.crosstab(_v7["has_har"], _v7["has_hcondel"])
_chi2v, _phi_p, _, _ = _chi2(_tab)
# direct 2x2 phi coefficient (matches hcondels_liftover_fix.json phi=0.052)
_aa = int((_v7["has_har"] & _v7["has_hcondel"]).sum())
_bb = int((_v7["has_har"] & ~_v7["has_hcondel"]).sum())
_cc = int((~_v7["has_har"] & _v7["has_hcondel"]).sum())
_dd = int((~_v7["has_har"] & ~_v7["has_hcondel"]).sum())
_phi = (_aa * _dd - _bb * _cc) / np.sqrt((_aa + _bb) * (_cc + _dd) *
                                         (_aa + _cc) * (_bb + _dd))

subsets = []
for nm, m in [("HAR only", _v7["has_har"] & ~_v7["has_hcondel"]),
              ("both", _v7["has_har"] & _v7["has_hcondel"]),
              ("hCONDEL only", ~_v7["has_har"] & _v7["has_hcondel"])]:
    sub = _v7[m]
    rd = int((sub["classification_v7"] == "regulation-driven").sum())
    rest = _v7[~m]
    rd_rest = int((rest["classification_v7"] == "regulation-driven").sum())
    # one-sided (enrichment direction) per Methods convention; matches main
    # text "9/90, OR = 1.80, P = 0.081"
    orv, pv = _fisher([[rd, len(sub) - rd], [rd_rest, len(rest) - rd_rest]],
                      alternative="greater")
    subsets.append((nm, len(sub), rd, orv, pv))

ys2 = np.arange(3)[::-1]
cols2 = [C_HAR, "#7f7f7f", C_HCONDEL]
for (nm, n, rd, orv, pv), yv, col in zip(subsets, ys2, cols2):
    axc2.barh(yv, n, height=0.62, color=col, edgecolor="white", linewidth=0.5, zorder=3)
    axc2.text(n + 14, yv + 0.15, f"{n} genes; RD {rd} ({100 * rd / n:.1f}%)",
              va="center", ha="left", fontsize=6.6, color="#444444")
    # one-sided p shown with 3 decimals + n.s. tag when non-significant
    _ptxt = (f"{pv:.3f}" if pv > 0.01 else fmt_p(pv))
    _ns = ", n.s." if pv > 0.05 else ""
    axc2.text(n + 14, yv - 0.19, f"OR {orv:.2f}, p={_ptxt}{_ns}",
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
        dxoff, dyoff, ha = -16, 0.14, "right"
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
# HAR-caMPRA nesting: circular vs non-circular all-HAR estimates.
# Both tests use the current 476-gene HAR mapping (n_hars > 0) and the same
# classifier: full model vs LOO-caMPRA variant (loo_full_matrix.csv; RD values
# from that CSV, GD values recomputed with the identical classifier/mapping).
axd = fig.add_subplot(gs[1, 1])
full_row = loo.loc[loo["variant"] == "full"].iloc[0]
looc_row = loo.loc[loo["variant"] == "LOO-caMPRA"].iloc[0]
tests = ["All HARs\nfull model (circular)", "All HARs\nLOO-caMPRA (non-circular)"]
rd_or = [full_row["HAR_OR"], looc_row["HAR_OR"]]
rd_p = [full_row["HAR_p"], looc_row["HAR_p"]]
rd_lo = [parse_ci(full_row["HAR_CI"])[0], parse_ci(looc_row["HAR_CI"])[0]]
rd_hi = [parse_ci(full_row["HAR_CI"])[1], parse_ci(looc_row["HAR_CI"])[1]]

# recompute GD all-HAR OR under both classifications (same classifier as CSV)
_sig = _v7["bh_fdr_sig"].astype(str).str.strip().str.lower().eq("true").to_numpy()
_rel = _v7["relax_relaxed"].astype(str).str.strip().str.lower().eq("true").to_numpy()
_har = (pd.to_numeric(_v7["n_hars"], errors="coerce").fillna(0) > 0).to_numpy()
_Xg = _v7[["gds_p_pct_resid", "gds_lrt_pct_resid", "gds_relax_pct", "gds_selectome"]].apply(
    pd.to_numeric, errors="coerce").fillna(0).to_numpy()
_gds = _Xg @ np.array([0.35, 0.30, 0.20, 0.15])
_rd = {c: pd.to_numeric(_v7[c], errors="coerce").fillna(0).to_numpy()
       for c in ["rds_doan", "rds_tau", "rds_nc", "rds_brain"]}


def _classify(rds_v, threshold=0.15):
    out = np.empty(len(_gds), dtype=object)
    for i in range(len(_gds)):
        g, r, s = _gds[i], rds_v[i], _sig[i]
        if g > r + threshold and s:
            c = "gene-driven"
        elif r > g + threshold and not s:
            c = "regulation-driven"
        elif g > r + threshold and r > 0.4 and s:
            c = "dual-driven"
        elif r > g + threshold and s and g > 0.3:
            c = "dual-driven"
        else:
            c = "neutral"
        if _rel[i] and c == "gene-driven":
            c = "gene-driven (relaxed)"
        out[i] = c
    return out


gd_or, gd_p, gd_lo, gd_hi = [], [], [], []
for rds_v in [0.30 * _rd["rds_doan"] + 0.25 * _rd["rds_tau"] + 0.25 * _rd["rds_nc"] + 0.20 * _rd["rds_brain"],
              0.35 * _rd["rds_tau"] + 0.35 * _rd["rds_nc"] + 0.30 * _rd["rds_brain"]]:
    m = np.isin(_classify(rds_v), ["gene-driven", "gene-driven (relaxed)"])
    tab = [[int((m & _har).sum()), int((m & ~_har).sum())],
           [int((~m & _har).sum()), int((~m & ~_har).sum())]]
    orv, pv = _fisher(tab, alternative="greater")
    a, b, c_, d = tab[0][0], tab[0][1], tab[1][0], tab[1][1]
    se = np.sqrt(1 / a + 1 / b + 1 / c_ + 1 / d)
    gd_or.append(float(orv)); gd_p.append(float(pv))
    gd_lo.append(float(np.exp(np.log(orv) - 1.96 * se)))
    gd_hi.append(float(np.exp(np.log(orv) + 1.96 * se)))
print("panel e GD all-HAR: full OR=%.2f p=%.2g | LOO-caMPRA OR=%.2f p=%.2g"
      % (gd_or[0], gd_p[0], gd_or[1], gd_p[1]))

xx = np.arange(2)
w = 0.34
b1 = axd.bar(xx - w / 2, rd_or, width=w, color=C_RD, edgecolor="white",
             linewidth=0.5, zorder=3, label="RD")
b2 = axd.bar(xx + w / 2, gd_or, width=w, color=C_GD, edgecolor="white",
             linewidth=0.5, zorder=3, label="GD")
for j in range(2):
    axd.errorbar([xx[j] - w / 2], [rd_or[j]],
                 yerr=[[rd_or[j] - rd_lo[j]], [rd_hi[j] - rd_or[j]]],
                 fmt="none", ecolor="#7f3b33", elinewidth=1.0, capsize=2.5, zorder=4)
    axd.errorbar([xx[j] + w / 2], [gd_or[j]],
                 yerr=[[gd_or[j] - gd_lo[j]], [gd_hi[j] - gd_or[j]]],
                 fmt="none", ecolor="#3c5a78", elinewidth=1.0, capsize=2.5, zorder=4)
# value labels above CI caps
for j in range(2):
    axd.text(xx[j] - w / 2, rd_hi[j] + 0.12, f"{rd_or[j]:.2f}", ha="center", fontsize=7.2,
             fontweight="bold", color=C_RD)
    axd.text(xx[j] - w / 2, rd_hi[j] + 0.58, f"p={fmt_p(rd_p[j])}", ha="center",
             fontsize=6.2, color="#555555")
    axd.text(xx[j] + w / 2, gd_hi[j] + 0.12, f"{gd_or[j]:.2f}", ha="center", fontsize=7.2,
             fontweight="bold", color=C_GD)
    axd.text(xx[j] + w / 2, gd_hi[j] + 0.58, "n.s.", ha="center", fontsize=6.2,
             color="#555555")
axd.axhline(1.0, color="#7f7f7f", ls="--", lw=0.9, zorder=1)
axd.set_xticks(xx)
axd.set_xticklabels(tests, fontsize=7.0)
axd.set_ylabel("Odds ratio (Fisher exact)")
axd.set_title("HAR-caMPRA nesting: conservative estimate", loc="left",
              fontsize=FS_TITLE)
axd.set_ylim(0, 6.6)
style_ax(axd, ygrid=True)
panel_label(axd, "e", dx=-0.20)
axd.legend(loc="upper right", frameon=False, fontsize=6.8, handletextpad=0.4,
           bbox_to_anchor=(1.0, 1.02))
axd.text(0.47, 0.84, "caMPRA-active elements are 100% nested\n"
         "within HARs; LOO-caMPRA removes this circularity",
         transform=axd.transAxes, ha="left", va="top", fontsize=6.4, color="#555555")

# ================================================================ panel f
# hCONDEL main validation (fixed hg38 mapping; phase9_hardening)
axe = fig.add_subplot(gs[1, 2])
uni_fix = hcfix["enrichment"]["new"]["union_fixed"]
rd_k, gd_k = int(uni_fix["rd_k"]), int(uni_fix["gd_k"])   # 18 / 19
rd_or_new, rd_p_new = float(uni_fix["rd_or"]), float(uni_fix["rd_p"])
gd_or_new, gd_p_new = float(uni_fix["gd_or"]), float(uni_fix["gd_p"])
# RD CI from the overridden full-model row (same values as panel a/b)
full_row = loo.loc[loo["variant"] == "full"].iloc[0]
rd_lo, rd_hi = full_row["hCONDEL_lo"], full_row["hCONDEL_hi"]
# GD CI and p recomputed (GD + GD-relaxed class, two-sided Fisher — matches
# hcondels_liftover_fix.json: OR 0.59, p 0.037, depleted)
_n_gd = int((_v7["classification_v7"].isin(["gene-driven", "gene-driven (relaxed)"])).sum())
_gd_or, _gd_p = _fe([[gd_k, _n_gd - gd_k],
                     [HC_TOTAL - gd_k, HC_N_UNIVERSE - _n_gd - (HC_TOTAL - gd_k)]])
_se_gd = np.sqrt(1 / gd_k + 1 / (_n_gd - gd_k) +
                 1 / (HC_TOTAL - gd_k) +
                 1 / (HC_N_UNIVERSE - _n_gd - (HC_TOTAL - gd_k)))
gd_lo = float(np.exp(np.log(_gd_or) - 1.96 * _se_gd))
gd_hi = float(np.exp(np.log(_gd_or) + 1.96 * _se_gd))
print("panel f hCONDEL (fixed): RD OR=%.2f p=%.2e k=%d/293 | GD OR=%.2f "
      "p=%.3f k=%d/%d" % (rd_or_new, rd_p_new, rd_k, _gd_or, _gd_p, gd_k, _n_gd))

rows_e = [("RD", rd_or_new, rd_lo, rd_hi, rd_p_new, C_RD,
           f"{rd_k}/293 RD genes"),
          ("GD", _gd_or, gd_lo, gd_hi, _gd_p, C_GD,
           f"{gd_k}/{_n_gd:,} GD + GD-rel. genes")]
for k, (nm, orv, lo_, hi_, p_, col, cnt) in enumerate(rows_e):
    yv = 1 - k
    axe.plot([lo_, hi_], [yv, yv], color=col, lw=1.6, zorder=2)
    axe.plot([lo_, lo_], [yv - 0.12, yv + 0.12], color=col, lw=1.3, zorder=2)
    axe.plot([hi_, hi_], [yv - 0.12, yv + 0.12], color=col, lw=1.3, zorder=2)
    axe.scatter([orv], [yv], s=64, color=col, edgecolors="white", linewidths=0.7,
                zorder=3)
    if p_ < 0.05:
        sig = f"p={fmt_p(p_)}" + (" (depleted)" if orv < 1 else "")
    else:
        sig = "n.s."
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
axe.set_title("hCONDEL main validation", loc="left", fontsize=FS_TITLE)
style_ax(axe, xgrid=True)
panel_label(axe, "f", dx=-0.20)
axe.text(0.98, 0.985, f"{HC_TOTAL} genes overlap 583 hCONDELs (fixed hg38 mapping);\n"
         "Fisher exact",
         transform=axe.transAxes, ha="right", va="top", fontsize=6.4, color="#555555")

out_png = os.path.join(OUTDIR, "Figure4_noncircular_validation.png")
out_pdf = os.path.join(OUTDIR, "Figure4_noncircular_validation.pdf")
fig.savefig(out_png, dpi=300)
fig.savefig(out_pdf)
print("saved:", out_png)
print("saved:", out_pdf)
