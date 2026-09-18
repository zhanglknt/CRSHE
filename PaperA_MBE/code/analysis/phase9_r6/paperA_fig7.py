# -*- coding: utf-8 -*-
"""
Paper A (MBE) — Figure 7: An integrative answer to the coding-versus-regulatory
question (4 panels, 2x2)
==============================================================================
Design system per results/paperB/figure_spec_v1.md, identical to paperA_fig{1..6}
(Arial; GD #2166ac, RD #d6604d, dual #e08214, neutral #bdbdbd; bold 11pt panel
labels; PNG 300dpi + vector PDF).

Panel data sources (statistics recomputed on the fly where feasible):
  a | Empirical disjointness (Tier 4): gene_classification_v7.csv +
      phase9_hardening/hcondel_gene_mapping_fixed.csv (repaired hg18->hg38
      liftOver, 112 genes) — BUSTED coding-selection rate and RD+dual rate
      in genome / HAR-proximate / hCONDEL / combined subsets
      (complement-based Fisher p; coding rate n.s. in all three subsets)
  b | Functional specificity (Tier 1): tier1_universe_enrichment CSVs —
      top RD SynGO + GO_BP terms (-log10 p, q<0.05 filled); GD shown as
      summary annotation (0 FDR terms)
  c | External GWAS anchoring (Tier 2): tier2_gwas_enrichment.csv —
      forest plot of ORs with 95% CI for RD and GD, 8 traits, log axis
  d | Coverage sensitivity (Tier 3): tier3_campra_sensitivity.csv —
      RD vs n_pos curve with extrapolation, GD line, crossing star,
      HAR-restricted ceiling

Output: results/paper/figures_v2/Figure7_integrative_answer.png (300 dpi) + .pdf
Stdout prints per-panel key statistics for manuscript citation.
"""
import os
import textwrap

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

if os.path.exists("/mnt/d"):
    ROOT = "/mnt/d/人类正选择基因项目"
else:
    ROOT = "D:/人类正选择基因项目"
V7 = os.path.join(ROOT, "results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv")
HCONDEL = os.path.join(ROOT, "results/phase9_hardening/hcondel_gene_mapping_fixed.csv")
T1_SYN = os.path.join(ROOT, "results/phase9_hardening/tier1_universe_enrichment_regulation_driven_SynGO.csv")
T1_GO = os.path.join(ROOT, "results/phase9_hardening/tier1_universe_enrichment_regulation_driven_GO_BP.csv")
T2 = os.path.join(ROOT, "results/phase9_hardening/tier2_gwas_enrichment.csv")
T3 = os.path.join(ROOT, "results/phase9_hardening/tier3_campra_sensitivity.csv")
T3J = os.path.join(ROOT, "results/phase9_hardening/tier3_campra_sensitivity.json")
OUTDIR = os.path.join(ROOT, "results/paper/figures_v2")
os.makedirs(OUTDIR, exist_ok=True)

C_GD, C_GD_RELAX, C_RD, C_DUAL, C_NEUT = "#2166ac", "#92c5de", "#d6604d", "#e08214", "#bdbdbd"
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "axes.linewidth": 0.7, "axes.edgecolor": "#333333",
    "xtick.direction": "out", "ytick.direction": "out",
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "axes.labelsize": 8.5,
    "axes.titlesize": 9.5, "legend.fontsize": 7.5,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})
FS_LABEL, FS_TITLE = 11, 9.5
REPORT = []


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


# Taller layout + wider inter-column gap so that panel (b)'s full (untruncated)
# two-line term labels fit between the columns without touching panel (a).
fig = plt.figure(figsize=(7.2, 7.8))
gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.75,
                      left=0.15, right=0.97, top=0.93, bottom=0.10,
                      height_ratios=[1.25, 1.0])

# ============================================================ (a) disjointness
df = pd.read_csv(V7)
hc = pd.read_csv(HCONDEL)
hc_ids = set(hc["gene_id"].dropna().astype(str))
hc_syms = set(hc["gene_symbol"].dropna().astype(str)) - {""}
sym_valid = df["gene_symbol"].notna() & (df["gene_symbol"].astype(str).str.strip() != "")
df["is_hc"] = df["gene_id"].isin(hc_ids) | (sym_valid & df["gene_symbol"].isin(hc_syms))

N = len(df)
sig = df["bh_fdr_sig"].values.astype(bool)
rd_dual = df["classification_v7"].isin(["regulation-driven", "dual-driven"]).values

subsets = [("Genome-wide (N=4,974)", np.ones(N, dtype=bool)),
           ("HAR-prox. (n=476)", (df["n_hars"] > 0).values),
           ("hCONDEL-prox. (n=112)", df["is_hc"].values),
           ("Combined (n=566)", ((df["n_hars"] > 0) | df["is_hc"]).values)]
# fixed-mapping sanity: hCONDEL subset must be the 112-gene union, combined 566
assert int(df["is_hc"].sum()) == 112, f"hCONDEL n={int(df['is_hc'].sum())} != 112"
assert int(((df["n_hars"] > 0) | df["is_hc"]).sum()) == 566, "combined n != 566"

ax = fig.add_subplot(gs[0, 0])
xpos = np.arange(len(subsets))
w = 0.38
coding_rates, rd_rates, p_coding, p_rd = [], [], [], []
for name, mask in subsets:
    k = mask.sum()
    coding = sig & mask
    rdd = rd_dual & mask
    coding_rates.append(100 * coding.sum() / k)
    rd_rates.append(100 * rdd.sum() / k)
    # Fisher vs rest
    rest = ~mask
    tab_c = [[int((sig & mask).sum()), int((~sig & mask).sum())],
             [int((sig & rest).sum()), int((~sig & rest).sum())]]
    orr, p = stats.fisher_exact(tab_c)
    p_coding.append(p)
    tab_r = [[int((rd_dual & mask).sum()), int((~rd_dual & mask).sum())],
             [int((rd_dual & rest).sum()), int((~rd_dual & rest).sum())]]
    orr, p = stats.fisher_exact(tab_r)
    p_rd.append(p)
    REPORT.append(f"[a] {name.replace(chr(10), ' ')}: coding {coding.sum()}/{k} = {100*coding.sum()/k:.1f}% "
                  f"(p={p_coding[-1]:.2g}); RD+dual {rdd.sum()}/{k} = {100*rdd.sum()/k:.1f}% (p={p_rd[-1]:.2g})")

b1 = ax.bar(xpos - w / 2, coding_rates, w, color=C_GD, zorder=3, label="Coding selection (BUSTED FDR<0.05)")
b2 = ax.bar(xpos + w / 2, rd_rates, w, color=C_RD, zorder=3, label="Regulation-driven + dual")
for i, (v, p) in enumerate(zip(coding_rates, p_coding)):
    ax.text(xpos[i] - w / 2, v + 0.8, f"{v:.1f}%", ha="center", fontsize=6.8, color=C_GD, zorder=5)
    if i > 0:
        ax.text(xpos[i] - w / 2, 1.2, "n.s." if p > 0.05 else f"p={p:.2g}", ha="center",
                fontsize=6.5, color="white", rotation=90, va="bottom", zorder=5)
    else:
        ax.text(xpos[i] - w / 2, 1.2, "ref", ha="center", fontsize=6.5, color="white",
                rotation=90, va="bottom", zorder=5)
for i, (v, p) in enumerate(zip(rd_rates, p_rd)):
    ax.text(xpos[i] + w / 2 + 0.07, v + 0.8, f"{v:.1f}%", ha="center", fontsize=6.8, color=C_RD, zorder=5)
    if i > 0:
        ax.text(xpos[i] + w / 2, 0.6, "***" if p < 1e-20 else ("**" if p < 1e-3 else "*"),
                ha="center", fontsize=7.5, color="white", va="bottom", zorder=5)
ax.set_xticks(xpos)
ax.set_xticklabels([s[0] for s in subsets], fontsize=6.5,
                   rotation=25, ha="right", rotation_mode="anchor")
ax.set_ylabel("Rate (%)")
ax.set_ylim(0, 48)
ax.set_title("Empirical disjointness of the two axes", fontsize=FS_TITLE)
# legend inside upper-left of panel (a): away from panel (b)'s long y-labels
ax.legend(loc="upper left", frameon=False, fontsize=6.2,
          handlelength=1.2, labelspacing=0.3, borderaxespad=0.2)
style_ax(ax, ygrid=True)
panel_label(ax, "a", dx=-0.13)

# ============================================================ (b) Tier 1 terms
t1s = pd.read_csv(T1_SYN)
t1g = pd.read_csv(T1_GO)
# top terms: 5 SynGO + 5 GO_BP by p
top_s = t1s.nsmallest(5, "p").copy()
top_s["lib"] = "SynGO"
top_g = t1g.nsmallest(5, "p").copy()
top_g["lib"] = "GO BP"
top = pd.concat([top_s, top_g]).sort_values("p", ascending=True).reset_index(drop=True)

ax = fig.add_subplot(gs[0, 1])
yy = np.arange(len(top))[::-1]
cols = [C_RD if lib == "SynGO" else C_DUAL for lib in top["lib"]]
filled = top["p_bh"] < 0.05
for i, (y, neglog, f) in enumerate(zip(yy, -np.log10(top["p"]), filled)):
    ax.plot([0, neglog], [y, y], color=cols[i], lw=1.2, zorder=2)
    ax.plot(neglog, y, "o", ms=6, color=cols[i], markerfacecolor=cols[i] if f else "white",
            markeredgecolor=cols[i], zorder=3)
# GD reference: none FDR — dashed line at top GD GO_BP term p
gd_top_p = 2.99e-03  # top GD GO_BP term p (from tier1 GD CSV)
ax.axvline(-np.log10(gd_top_p), color=C_GD, ls="--", lw=1.2)
ax.set_ylim(-0.6, 10.6)
ax.text(-np.log10(gd_top_p) + 0.08, 10.5, "GD best term\n(0 FDR terms)",
        fontsize=6.5, color=C_GD, va="top")

# full, untruncated term names, wrapped to at most 2 lines of <=33 characters
labels = []
for _, r in top.iterrows():
    name = r["term"].split(" (GO:")[0]
    wrapped = textwrap.wrap(name, width=33)
    assert len(wrapped) <= 2, f"term needs >2 lines: {name}"
    labels.append("\n".join(wrapped))
ax.set_yticks(yy)
ax.set_yticklabels(labels, fontsize=6.2, linespacing=1.05)
ax.set_xlabel("$-\\log_{10}(P)$, hypergeometric\nvs 4,974-gene universe")
ax.set_title("Functional specificity (RD class)", fontsize=FS_TITLE)
ax.set_xlim(0, 7.2)
from matplotlib.lines import Line2D
ax.legend(handles=[
    Line2D([], [], marker="o", ls="-", color=C_RD, label="SynGO (11 FDR terms)", ms=5),
    Line2D([], [], marker="o", ls="-", color=C_DUAL, label="GO BP (18 FDR terms)", ms=5)],
    loc="lower right", frameon=False, fontsize=6.2,
    handlelength=1.2, labelspacing=0.3, borderaxespad=0.2)
style_ax(ax, xgrid=True)
panel_label(ax, "b", dx=-0.30)
for _, r in top.iterrows():
    REPORT.append(f"[b] {r['term']}: k={r['k_overlap']}/K={r['K_universe']} p={r['p']:.2e} q={r['p_bh']:.2e}")

# ============================================================ (c) GWAS forest
t2 = pd.read_csv(T2)
traits = ["Educational attainment", "Schizophrenia", "Major depression", "Bipolar disorder",
          "Cognitive performance", "Neuroticism", "ASD", "Intelligence"]
trait_labels = ["Educational\nattainment", "Schizophrenia", "Major\ndepression", "Bipolar\ndisorder",
                "Cognitive\nperformance", "Neuroticism", "ASD", "Intelligence"]
uni_n = 4974

ax = fig.add_subplot(gs[1, 0])
yy = np.arange(len(traits))[::-1]
for i, trait in enumerate(traits):
    for cls, color, dx in [("regulation_driven", C_RD, 0.0), ("gene_driven", C_GD, 0.0)]:
        row = t2[(t2["trait"] == trait) & (t2["class"] == cls)].iloc[0]
        k, n_cls, K = int(row["k_overlap"]), int(row["n_class"]), int(row["K_universe"])
        a, b = k, n_cls - k
        c, d = K - k, uni_n - n_cls - (K - k)
        se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d) if min(a, c) > 0 else 1.0
        orr = row["OR"]
        lo, hi = np.exp(np.log(orr) - 1.96 * se), np.exp(np.log(orr) + 1.96 * se)
        y = yy[i] + (0.18 if cls == "regulation_driven" else -0.18)
        ax.plot([lo, hi], [y, y], color=color, lw=1.4, zorder=2)
        ax.plot(orr, y, "s", ms=5, color=color, zorder=3)
        REPORT.append(f"[c] {trait} {cls}: k={k}/{n_cls} K={K} OR={orr:.2f} [{lo:.2f},{hi:.2f}] q={row['q_bh']:.3g}")
ax.axvline(1, color="#888888", lw=0.8, ls="-")
ax.set_yticks(yy)
ax.set_yticklabels(trait_labels, fontsize=6.8)
ax.set_xscale("log")
ax.set_xticks([0.25, 0.5, 1, 2, 4])
ax.set_xticklabels(["0.25", "0.5", "1", "2", "4"])
ax.set_xlabel("Odds ratio (95% CI), GWAS $P<5\\times10^{-8}$")
ax.set_title("External GWAS trait anchoring", fontsize=FS_TITLE)
ax.set_xlim(0.2, 6.5)
ax.legend(handles=[Line2D([], [], marker="s", ls="-", color=C_RD, label="regulation-driven (8/8 FDR-sig.)", ms=5),
                   Line2D([], [], marker="s", ls="-", color=C_GD, label="gene-driven (0/8; OR 0.56–1.03)", ms=5)],
          loc="lower right", frameon=False, fontsize=6.0,
          handlelength=1.2, labelspacing=0.3, borderaxespad=0.2)
style_ax(ax, xgrid=True)
panel_label(ax, "c", dx=-0.16)

# ============================================================ (d) Tier 3 curve
t3 = pd.read_csv(T3)
import json
with open(T3J, encoding="utf-8") as f:
    t3j = json.load(f)
slope = t3j["linear_fit"]["slope"]
intercept = t3j["linear_fit"]["intercept"]
n_cross = t3j["crossing"]["n_pos_star"]
gd_ref = t3j["gd_reference"]
n_har_ceiling = t3j["ceilings"]["doan_har_only"]["regulation-driven"]

ax = fig.add_subplot(gs[1, 1])
ax.errorbar(t3["n_pos"], t3["RD_mean"], yerr=t3["RD_sd"], fmt="o-", color=C_RD,
            capsize=2, lw=1.6, ms=4.5, label="RD count (100 reps)")
ax.axhline(gd_ref, color=C_GD, ls="--", lw=1.6, label=f"GD count = {gd_ref:,}")
xx = np.linspace(0, 3200, 100)
ax.plot(xx, slope * xx + intercept, color=C_RD, ls=":", lw=1.3, alpha=0.85)
ax.plot([n_cross], [gd_ref], marker="*", ms=15, color="#7B3294", zorder=5)
ax.annotate(f"crossing n* = {n_cross:,.0f}\n(59% of universe)",
            xy=(n_cross, gd_ref), xytext=(1650, 700),
            fontsize=6.8, color="#7B3294",
            arrowprops=dict(arrowstyle="->", color="#7B3294", lw=0.9))
# observed-range shading + HAR ceiling
ax.axvspan(0, 97, alpha=0.08, color="green")
ax.text(112, 130, "observed\n(97 = 2%)", ha="left", va="center", fontsize=6.3, color="green")
ax.plot([476], [n_har_ceiling], marker="D", ms=6, color=C_DUAL, zorder=5)
ax.annotate(f"HAR-restricted ceiling\n({n_har_ceiling} at 476 possible)",
            xy=(476, n_har_ceiling), xytext=(600, 1000), fontsize=6.8, color=C_DUAL,
            arrowprops=dict(arrowstyle="->", color=C_DUAL, lw=0.9))
ax.set_xlabel("caMPRA assay-positive genes\n(regulatory coverage)")
ax.set_ylabel("Regulation-driven count")
ax.set_title("Coverage sensitivity", fontsize=FS_TITLE)
ax.set_xlim(-60, 3300)
ax.set_ylim(0, 2050)
ax.legend(loc="upper left", frameon=False, fontsize=6.5)
style_ax(ax, ygrid=True)
panel_label(ax, "d", dx=-0.16)
REPORT.append(f"[d] slope={slope:.3f} intercept={intercept:.1f} crossing={n_cross:.0f} "
              f"HAR ceiling RD={n_har_ceiling}")

# ============================================================ save
fig.suptitle("", y=0.98)
fig.savefig(os.path.join(OUTDIR, "Figure7_integrative_answer.png"), dpi=300)
fig.savefig(os.path.join(OUTDIR, "Figure7_integrative_answer.pdf"))
print("\n".join(REPORT))
print("saved:", os.path.join(OUTDIR, "Figure7_integrative_answer.png"))
