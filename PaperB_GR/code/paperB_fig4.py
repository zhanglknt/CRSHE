# -*- coding: utf-8 -*-
"""
Paper B — Figure 4: Development, brain regions, independent validation
6 panels, 2x3. Publication grade per results/paperB/figure_spec_v1.md.

Panel data sources (all under results/phase8_tissue_analysis/):
  a | developmental_full_trajectory.csv   (GD/RD enrichment OR across 31 BrainSpan stages;
      25 pcw and 35 pcw have n_high_expr=0 -> shown as gaps. NOTE: spec said "mean expression"
      but the file holds enrichment OR per stage; OR trajectory plotted instead.)
  b | developmental_tau_comparison.csv    (per-gene developmental tau, LOO-brain classes;
      overall MWU p from developmental_analysis_summary.json)
  c | brainspan_region_enrichment.csv     (18 regions with computed enrichment; summary json
      lists 26 atlas regions, enrichment available for 18 -> data wins)
  d | prenatal_vs_postnatal.csv           (per-gene log2 prenatal/postnatal ratio by class)
  e | hcondels_summary.json               (RD OR=2.94 p=6.65e-06; GD OR=0.53 n.s.)
  f | campra_regulatory_summary.json      (active-HAR LOO RD OR=1.84 p=0.012; all-HAR v7 RD OR=4.79)

Output: results/paperB/figures/Figure4_development_validation.png (300 dpi) + .pdf
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
C_GD = "#2166ac"
C_GD_RELAX = "#92c5de"
C_RD = "#d6604d"
C_NEUR = "#b2182b"
C_DUAL = "#e08214"
C_NEUT = "#bdbdbd"
C_RD_NS = "#fddbc7"
C_GD_NS = "#d1e5f0"
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
FS_LABEL = 11
FS_TITLE = 9.5


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


def sci(p, digits=1):
    """format p value as x.x×10^-n (mathtext)"""
    if p >= 0.001:
        return f"{p:.{digits}f}".rstrip("0").rstrip(".")
    m, e = f"{p:.0e}".split("e")
    return f"{m}$\\times$10$^{{{int(e)}}}$"


# ---------------------------------------------------------------- load data
traj = pd.read_csv(os.path.join(DATA, "developmental_full_trajectory.csv"))
tau = pd.read_csv(os.path.join(DATA, "developmental_tau_comparison.csv"))
reg = pd.read_csv(os.path.join(DATA, "brainspan_region_enrichment.csv"))
pvn = pd.read_csv(os.path.join(DATA, "prenatal_vs_postnatal.csv"))
with open(os.path.join(DATA, "hcondels_summary.json"), encoding="utf-8") as f:
    hc = json.load(f)
with open(os.path.join(DATA, "campra_regulatory_summary.json"), encoding="utf-8") as f:
    ca = json.load(f)
with open(os.path.join(DATA, "developmental_analysis_summary.json"), encoding="utf-8") as f:
    dev_sum = json.load(f)

# ---------------------------------------------------------------- figure layout
fig = plt.figure(figsize=(11, 8))
gs = fig.add_gridspec(2, 3, left=0.075, right=0.99, top=0.93, bottom=0.115,
                      wspace=0.42, hspace=0.55, width_ratios=[1.25, 1.0, 1.0])

# ================================================================ panel a
# developmental enrichment trajectory across 31 BrainSpan stages
axa = fig.add_subplot(gs[0, 0])
t = traj.copy()
t.loc[t["n_high_expr"] == 0, ["gd_OR", "rd_OR"]] = np.nan
x = np.arange(len(t))
pre_end = 12  # index of '37 pcw' -> last prenatal stage
axa.axvspan(-0.6, pre_end + 0.5, color="#f2f2f2", zorder=0)
axa.axhline(1.0, color="#7f7f7f", linestyle="--", linewidth=0.8, zorder=1)
for col, pcol, c, lab in [("gd_OR", "gd_p", C_GD, "GD"), ("rd_OR", "rd_p", C_RD, "RD")]:
    axa.plot(x, t[col], color=c, linewidth=1.2, zorder=2)
    sig = t[pcol] < 0.05
    axa.scatter(x[sig], t.loc[sig, col], s=16, color=c, zorder=3, linewidth=0)
    axa.scatter(x[~sig & t[col].notna()], t.loc[~sig & t[col].notna(), col],
                s=14, facecolor="white", edgecolor=c, linewidth=0.8, zorder=3)
# mark the two no-data stages
for i in np.where(traj["n_high_expr"] == 0)[0]:
    axa.text(i, 0.865, "n.d.", ha="center", va="top", fontsize=6.2, color="#999999",
             rotation=90)
axa.text((pre_end - 0.5) / 2, 1.335, "prenatal", ha="center", fontsize=7.5, color="#555555")
axa.text((pre_end + 0.5 + len(t)) / 2, 1.335, "postnatal", ha="center", fontsize=7.5,
         color="#555555")
axa.set_xticks(x)
axa.set_xticklabels(t["stage"], rotation=90, fontsize=6.3)
axa.set_xlim(-0.6, len(t) - 0.4)
axa.set_ylim(0.85, 1.36)
axa.set_ylabel("Odds ratio of high-expression enrichment")
axa.set_xlabel("BrainSpan developmental stage", labelpad=2)
axa.set_title("Enrichment across 31 brain development stages", loc="left", fontsize=FS_TITLE)
style_ax(axa, ygrid=True)
panel_label(axa, "a", dx=-0.11)
leg_a = [Line2D([0], [0], color=C_GD, linewidth=1.2, marker="o", markersize=4.5, label="GD"),
         Line2D([0], [0], color=C_RD, linewidth=1.2, marker="o", markersize=4.5, label="RD"),
         Line2D([0], [0], marker="o", color="none", markerfacecolor="white",
                markeredgecolor="#555555", markersize=4.5, label="nominal P≥0.05 (open)")]
axa.legend(handles=leg_a, loc="lower left", frameon=False, borderaxespad=0.05,
           handletextpad=0.4, labelspacing=0.3, fontsize=6.8, ncol=1,
           bbox_to_anchor=(0.0, 0.02))

# ================================================================ panel b
# developmental tau: GD vs RD
axb = fig.add_subplot(gs[0, 1])
gd_tau = tau.loc[tau["classification_loo_brain"] == "gene-driven", "dev_tau"].to_numpy()
rd_tau = tau.loc[tau["classification_loo_brain"] == "regulation-driven", "dev_tau"].to_numpy()
p_tau = dev_sum["dev_tau_rd_vs_gd_p"]
vp = axb.violinplot([gd_tau, rd_tau], positions=[0, 1], widths=0.62,
                    showextrema=False)
for body, c in zip(vp["bodies"], [C_GD, C_RD]):
    body.set_facecolor(c)
    body.set_alpha(0.35)
    body.set_edgecolor("none")
bp = axb.boxplot([gd_tau, rd_tau], positions=[0, 1], widths=0.16, patch_artist=True,
                 showfliers=False, medianprops=dict(color="black", linewidth=1.1),
                 whiskerprops=dict(linewidth=0.7), capprops=dict(linewidth=0.7),
                 boxprops=dict(linewidth=0.8))
for patch, c in zip(bp["boxes"], [C_GD, C_RD]):
    patch.set_facecolor(c)
# Round-2 adjudication: dev-tau GD vs RD MWU one-sided P = 0.544 (summary json);
# two-sided P = 0.91 (verifier-manuscript confirmed = one-sided x 2 exactly).
axb.text(0.5, 1.065, "two-sided MWU P = 0.91 (n.s.)", ha="center", fontsize=7.8,
         color="#333333")
axb.set_xticks([0, 1])
axb.set_xticklabels([f"GD (n = {len(gd_tau)})\nmedian {np.median(gd_tau):.3f}",
                     f"RD (n = {len(rd_tau)})\nmedian {np.median(rd_tau):.3f}"], fontsize=8)
axb.set_ylim(0, 1.12)
axb.set_ylabel("Developmental tau (31 stages)")
axb.set_title("No GD–RD difference in stage specificity", loc="left", fontsize=FS_TITLE)
style_ax(axb, ygrid=True)
panel_label(axb, "b", dx=-0.20)

# ================================================================ panel c
# BrainSpan regions: GD vs RD enrichment (18 regions with data)
axc = fig.add_subplot(gs[0, 2])
rg = reg.sort_values("gd_OR", ascending=True).reset_index(drop=True)
yc = np.arange(len(rg))
for yi, (_, row) in zip(yc, rg.iterrows()):
    for off, orv, pv, c in [(0.19, row["gd_OR"], row["gd_p"], C_GD),
                            (-0.19, row["rd_OR"], row["rd_p"], C_RD)]:
        if pv < 0.05:
            axc.scatter(orv, yi + off, s=17, color=c, zorder=3, linewidth=0)
        else:
            axc.scatter(orv, yi + off, s=15, facecolor="white", edgecolor=c,
                        linewidth=0.8, zorder=3)
axc.axvline(1.0, color="#7f7f7f", linestyle="--", linewidth=0.8, zorder=1)
axc.set_yticks(yc)
axc.set_yticklabels(rg["region"], fontsize=7)
axc.set_xlim(0.96, 1.32)
axc.set_ylim(-0.7, len(rg) - 0.3)
axc.set_xlabel("Odds ratio (LOO-brain classification)")
axc.set_title("18 brain regions (LOO-brain)", loc="left", fontsize=FS_TITLE)
style_ax(axc, xgrid=True)
panel_label(axc, "c", dx=-0.20)
n_gd_nom = int((rg["gd_p"] < 0.05).sum())
n_rd_nom = int((rg["rd_p"] < 0.05).sum())
axc.text(0.0, -0.185, f"nominal P<0.05: GD {n_gd_nom}/18, RD {n_rd_nom}/18; no FDR<0.05",
         transform=axc.transAxes, ha="left", va="top", fontsize=6.8, color="#555555",
         clip_on=False)
leg_c = [Line2D([0], [0], marker="o", color="none", markerfacecolor=C_GD, markersize=5, label="GD"),
         Line2D([0], [0], marker="o", color="none", markerfacecolor=C_RD, markersize=5, label="RD"),
         Line2D([0], [0], marker="o", color="none", markerfacecolor="white",
                markeredgecolor="#555555", markersize=5, label="n.s. (open)")]
axc.legend(handles=leg_c, loc="lower left", frameon=True, borderaxespad=0.15,
           handletextpad=0.3, labelspacing=0.3, fontsize=6.8,
           facecolor="white", framealpha=0.95, edgecolor="none",
           bbox_to_anchor=(0.01, 0.01))

# ================================================================ panel d
# prenatal vs postnatal expression shift by class
axd = fig.add_subplot(gs[1, 0])
pvn2 = pvn[(pvn["prenatal_mean"] > 0) & (pvn["postnatal_mean"] > 0)].copy()
pvn2["l2r"] = np.log2(pvn2["prenatal_mean"] / pvn2["postnatal_mean"])
groups = [("gene-driven", "GD", C_GD), ("regulation-driven", "RD", C_RD),
          ("neutral", "Neutral", C_NEUT)]
arrs = [pvn2.loc[pvn2["classification_loo_brain"] == k, "l2r"].to_numpy() for k, _, _ in groups]
p_mwu = stats.mannwhitneyu(arrs[1], arrs[0]).pvalue
vp = axd.violinplot(arrs, positions=[0, 1, 2], widths=0.66, showextrema=False)
for body, (_, _, c) in zip(vp["bodies"], groups):
    body.set_facecolor(c)
    body.set_alpha(0.35)
    body.set_edgecolor("none")
bp = axd.boxplot(arrs, positions=[0, 1, 2], widths=0.15, patch_artist=True,
                 showfliers=False, medianprops=dict(color="black", linewidth=1.1),
                 whiskerprops=dict(linewidth=0.7), capprops=dict(linewidth=0.7),
                 boxprops=dict(linewidth=0.8))
for patch, (_, _, c) in zip(bp["boxes"], groups):
    patch.set_facecolor(c)
axd.axhline(0, color="#7f7f7f", linestyle="--", linewidth=0.8, zorder=1)
axd.set_xticks([0, 1, 2])
axd.set_xticklabels([f"{lab}\n(n = {len(a)})" for a, (_, lab, _) in zip(arrs, groups)],
                    fontsize=7.8)
axd.set_ylabel("log$_2$(prenatal / postnatal)")
axd.set_title("Prenatal-to-postnatal shift: no class effect", loc="left", fontsize=FS_TITLE)
# A8 (review round 1, adjudicated): two-sided MWU on log2(prenatal_mean/postnatal_mean),
# both means > 0, LOO-brain GD (excl. relaxed, n=1085) vs RD (n=605): P = 0.5673 (n.s.).
# NOTE: the 'ratio' column in prenatal_vs_postnatal.csv is corrupted (1.0 fill / inf);
# always recompute from prenatal_mean/postnatal_mean (as done above via pvn2["l2r"]).
axd.text(0.5, 0.965, "RD vs GD: two-sided MWU P = 0.57 (n.s.)", transform=axd.transAxes,
         ha="center", va="top", fontsize=7.8, color="#333333")
style_ax(axd, ygrid=True)
panel_label(axd, "d", dx=-0.22)

# ================================================================ panel e
# hCONDEL independent validation
axe = fig.add_subplot(gs[1, 1])
vals = [hc["gd_enrichment"]["OR"], hc["rd_enrichment"]["OR"]]
ps = [hc["gd_enrichment"]["p"], hc["rd_enrichment"]["p"]]
xe = [0, 1]
bars = axe.bar(xe, vals, width=0.52, color=[C_GD, C_RD], zorder=2)
for xi, v, p in zip(xe, vals, ps):
    st = "***" if p < 1e-3 else ("**" if p < 1e-2 else ("*" if p < 0.05 else "n.s."))
    if v > 2.0:  # tall bar: white label inside
        axe.text(xi, v - 0.14, f"OR = {v:.2f} {st}", ha="center", va="top",
                 fontsize=7.8, color="white", fontweight="bold")
    else:
        axe.text(xi, v + 0.09, f"OR = {v:.2f}\n{st}", ha="center", va="bottom",
                 fontsize=7.8, color="#333333")
axe.axhline(1.0, color="#7f7f7f", linestyle="--", linewidth=0.8, zorder=1)
axe.set_xticks(xe)
axe.set_xticklabels(["GD", "RD"], fontsize=8.5)
axe.set_ylim(0, 3.75)
axe.set_ylabel("Odds ratio")
axe.set_title("hCONDELs: independent RD enrichment", loc="left", fontsize=FS_TITLE)
axe.text(0.04, 0.955, f"RD P = {sci(hc['rd_enrichment']['p'])}",
         transform=axe.transAxes, ha="left", va="top", fontsize=7.2, color="#555555")
axe.text(0.04, 0.865, f"{hc['n_hcondels_parsed']} hCONDELs in {hc['n_genes_total']} genes; "
                      "non-circular", transform=axe.transAxes, ha="left", va="top",
         fontsize=7.2, color="#555555")
style_ax(axe, ygrid=True)
panel_label(axe, "e", dx=-0.20)

# ================================================================ panel f
# caMPRA active HARs (LOO) vs all HARs (v7)
axf = fig.add_subplot(gs[1, 2])
ctx = ["Active HARs\n(LOO, non-circular)", "All HARs\n(v7)"]
rd_v = [ca["rd_loo_campra"]["OR"], ca["rd_v7_all_hars"]["OR"]]
gd_v = [ca["gd_loo_campra"]["OR"], ca["gd_v7_all_hars"]["OR"]]
rd_p = [ca["rd_loo_campra"]["p"], ca["rd_v7_all_hars"]["p"]]
gd_p = [ca["gd_loo_campra"]["p"], ca["gd_v7_all_hars"]["p"]]
wbar = 0.34
xf = np.arange(2)
axf.bar(xf - wbar / 2 - 0.01, gd_v, width=wbar, color=C_GD, zorder=2, label="GD")
axf.bar(xf + wbar / 2 + 0.01, rd_v, width=wbar, color=C_RD, zorder=2, label="RD")
for xi, v, p in zip(xf - wbar / 2 - 0.01, gd_v, gd_p):
    st = "n.s." if p >= 0.05 else "***"
    axf.text(xi, v + 0.10, f"{v:.2f} {st}", ha="center", va="bottom", fontsize=7,
             color="#333333")
for xi, v, p in zip(xf + wbar / 2 + 0.01, rd_v, rd_p):
    st = "*" if 0.01 <= p < 0.05 else ("**" if p < 0.01 and p >= 1e-3 else "***")
    if p >= 0.05:
        st = "n.s."
    axf.text(xi, v + 0.10, f"{v:.2f} {st}", ha="center", va="bottom", fontsize=7,
             color="#333333")
axf.axhline(1.0, color="#7f7f7f", linestyle="--", linewidth=0.8, zorder=1)
axf.set_xticks(xf)
axf.set_xticklabels(ctx, fontsize=7.8)
axf.set_ylim(0, 5.6)
axf.set_ylabel("Odds ratio")
axf.set_title("caMPRA-active HARs near RD genes", loc="left", fontsize=FS_TITLE)
axf.legend(loc="upper left", frameon=False, borderaxespad=0.1, handletextpad=0.4,
           labelspacing=0.3, fontsize=7.2)
axf.text(0.03, 0.76, f"active-HAR RD: P = {ca['rd_loo_campra']['p']:.3f}\n"
                     f"all-HAR RD: P = {sci(ca['rd_v7_all_hars']['p'])}",
         transform=axf.transAxes, ha="left", va="top", fontsize=6.9, color="#555555")
style_ax(axf, ygrid=True)
panel_label(axf, "f", dx=-0.20)

# ---------------------------------------------------------------- save
png = os.path.join(OUTDIR, "Figure4_development_validation.png")
pdf = os.path.join(OUTDIR, "Figure4_development_validation.pdf")
fig.savefig(png, dpi=300)
fig.savefig(pdf, bbox_inches="tight")
print("saved:", png)
print("saved:", pdf)

# ---------------------------------------------------------------- verification summary
print("\n=== Figure 4 per-panel summary ===")
print(f"a | developmental_full_trajectory.csv | 31 stages (25/35 pcw no-data gaps) | "
      f"GD OR {traj['gd_OR'].min():.2f}-{traj['gd_OR'].max():.2f}, RD OR {traj['rd_OR'].min():.2f}-{traj['rd_OR'].max():.2f}; "
      f"GD nominal P<0.05 in {(traj['gd_p']<0.05).sum()}/31 stages, RD in {(traj['rd_p']<0.05).sum()}/31")
print(f"b | developmental_tau_comparison.csv | GD median tau {np.median(gd_tau):.3f} (n={len(gd_tau)}) "
      f"vs RD {np.median(rd_tau):.3f} (n={len(rd_tau)}); two-sided MWU P=0.91 "
      f"(one-sided summary json: {p_tau:.3f})")
print(f"c | brainspan_region_enrichment.csv | 18 regions (of 26 atlas regions) | "
      f"GD nominal {n_gd_nom}/18, RD nominal {n_rd_nom}/18, 0 FDR<0.05")
print(f"d | prenatal_vs_postnatal.csv | finite-ratio genes: GD {len(arrs[0])}, RD {len(arrs[1])}, "
      f"neutral {len(arrs[2])}; adjudicated two-sided MWU P=0.57 (fig4d_adjudication); "
      f"script recomputation (log2 of prenatal_mean/postnatal_mean, both>0): P={p_mwu:.4f}")
print(f"e | hcondels_summary.json | RD OR {hc['rd_enrichment']['OR']:.2f} (P={hc['rd_enrichment']['p']:.1e}) "
      f"vs GD OR {hc['gd_enrichment']['OR']:.2f} (P={hc['gd_enrichment']['p']:.2f})")
print(f"f | campra_regulatory_summary.json | active-HAR LOO: RD {ca['rd_loo_campra']['OR']:.2f} "
      f"(P={ca['rd_loo_campra']['p']:.3f}), GD {ca['gd_loo_campra']['OR']:.2f}; "
      f"all-HAR v7: RD {ca['rd_v7_all_hars']['OR']:.2f} (P={ca['rd_v7_all_hars']['p']:.1e})")
