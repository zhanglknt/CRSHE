# -*- coding: utf-8 -*-
"""
Paper B — Figure 2: Bulk-tissue level — coding selection ubiquitous,
regulatory flat (5 panels, layout 2+3). Genome Research design system.

Panel data sources
------------------
a,b | GTEx 30-tissue GD/RD enrichment:
      results/phase8_tissue_analysis/layer1_tissue_enrichment.csv
c   | tau distributions (LOO-tau, non-circular):
      results/phase8_tissue_analysis/data/phase8_master_table.csv (tau),
      results/phase8_tissue_analysis/layer2_statistical_tests.csv
d   | brain vs non-brain:
      results/phase8_tissue_analysis/layer1_brain_vs_nonbrain.csv
e   | per-tissue specificity:
      results/phase8_tissue_analysis/layer2_per_tissue_specificity.csv

Output: results/paperB/figures/Figure2_bulk_tissue.png (300 dpi) + .pdf
"""
import os
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

ROOT = os.environ.get("HSD_BASE") or (
    "/mnt/d/人类正选择基因项目" if os.path.exists("/mnt/d") else "D:/人类正选择基因项目")
P8 = os.path.join(ROOT, "results", "phase8_tissue_analysis")
OUTDIR = os.path.join(ROOT, "results", "paperB", "figures")
os.makedirs(OUTDIR, exist_ok=True)

GD = "#2166ac"; GDD = "#08306b"; RD = "#d6604d"; NHI = "#b2182b"
DUAL = "#e08214"; NEUT = "#bdbdbd"; GREY = "#737373"

plt.rcParams.update({
    "font.family": ["Arial", "DejaVu Sans"],
    "font.size": 8,
    "axes.titlesize": 9.5,
    "axes.labelsize": 8.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42,
})


def panel_label(ax, letter):
    ax.set_title(letter, loc="left", fontweight="bold", fontsize=11, pad=6)


# ---------------------------------------------------------------- data -----
# S10: unified LOO-tau variant (previous version mixed v7 / LOO-brain)
enr = pd.read_csv(os.path.join(P8, "layer1_tissue_enrichment_loo_tau.csv"))
bvn = pd.read_csv(os.path.join(P8, "layer1_brain_vs_nonbrain.csv"))
master = pd.read_csv(os.path.join(P8, "data", "phase8_master_table.csv"))
tests = pd.read_csv(os.path.join(P8, "layer2_statistical_tests.csv"))
spec = pd.read_csv(os.path.join(P8, "layer2_per_tissue_specificity.csv"))
variants = json.load(open(os.path.join(P8, "layer1_loo_variants_summary.json")))

enr = enr.sort_values("gd_OR", ascending=True).reset_index(drop=True)
assert len(enr) == 30
# S10: LOO-tau OR ranges — GD 1.19-1.41, RD 0.44-0.96; shared scale
XLIM = (0.2, 1.62)
summary = []

fig = plt.figure(figsize=(11, 7.4))
gs = gridspec.GridSpec(2, 6, figure=fig, height_ratios=[1.15, 1.0],
                       hspace=0.62, wspace=1.15,
                       left=0.10, right=0.975, top=0.94, bottom=0.08)

# ============================ panels a/b: GD & RD forests across 30 tissues
ypos = np.arange(len(enr))
brain_mask = enr["is_brain"].values

for col_idx, (prefix, color, color_brain, letter, title) in enumerate([
        ("gd", GD, GDD, "a",
         "GD enrichment, LOO-tau — FDR<0.05 in all 30 GTEx tissues"),
        ("rd", RD, NHI, "b",
         "RD enrichment, LOO-tau — FDR<0.05 in 0/30 tissues")]):
    ax = fig.add_subplot(gs[0, col_idx * 3:(col_idx + 1) * 3])
    panel_label(ax, letter)
    ax.set_title(title, loc="center", pad=6, fontsize=9)
    ors = enr[f"{prefix}_OR"].values
    cols = np.where(brain_mask, color_brain, color)
    ax.barh(ypos, ors, height=0.66, color=cols, edgecolor="none", zorder=2)
    ax.axvline(1.0, color=GREY, ls="--", lw=0.9, zorder=1)
    ax.set_xlim(*XLIM)
    ax.set_ylim(-0.8, len(enr) - 0.2)
    ax.set_yticks(ypos)
    if col_idx == 0:
        ax.set_yticklabels(enr["tissue"], fontsize=7.2)
        # A14: bold strictly limited to the two brain-related tissues
        for tick, is_b in zip(ax.get_yticklabels(), brain_mask):
            if is_b:
                tick.set_fontweight("bold")
                tick.set_color(color_brain)
    else:
        ax.set_yticklabels([])
    ax.set_xlabel("Enrichment odds ratio")
    ax.xaxis.set_major_locator(plt.MultipleLocator(0.4))
    ax.grid(axis="x", color="#e9e9e9", lw=0.6)
    ax.set_axisbelow(True)
    # OR value labels just right of each bar
    for y, v in zip(ypos, ors):
        ax.text(v + 0.022, y, f"{v:.2f}", va="center", ha="left",
                fontsize=6.2, color="#555555", zorder=4)

axa = fig.axes[0]
axa.text(0.985, 0.015, "bold = brain-related", transform=axa.transAxes,
         fontsize=6.8, color="#333333", ha="right", va="bottom",
         style="italic")
axb = fig.axes[1]
axb.text(0.97, 0.55, "no tissue significant;\nOR < 1 = depletion among\n"
         "highly expressed genes", transform=axb.transAxes, fontsize=7.2,
         color="#333333", ha="right", va="center", style="italic")

summary.append(f"a | GD OR {enr.gd_OR.min():.2f}-{enr.gd_OR.max():.2f}, "
               f"FDR<0.05 in {int(enr['gd_fdr_sig'].sum())}/30 "
               "(brain: Brain, Nerve highlighted)")
summary.append(f"b | RD OR {enr.rd_OR.min():.2f}-{enr.rd_OR.max():.2f}, "
               f"FDR<0.05 in {int(enr['rd_fdr_sig'].sum())}/30 "
               "(all depleted, OR<1)")

# ============================ panel c: tau distributions (LOO, non-circular)
axc = fig.add_subplot(gs[1, 0:2])
panel_label(axc, "c")
axc.set_title("Tissue specificity (tau), LOO-tau",
              loc="center", pad=6, fontsize=9)

gd_tau = master.loc[master["classification_loo_tau"] == "gene-driven",
                    "tau"].dropna().values
rd_tau = master.loc[master["classification_loo_tau"] == "regulation-driven",
                    "tau"].dropna().values
row = tests[tests["tau_source"] == "GTEx V11 (30 tissues)"].iloc[0]
# layer2_statistical_tests.csv was re-run two-sided (R1 fix): mw_p_loo and
# mw_p_circular are now already two-sided (0.412 / 3.24e-25) — use directly
p_loo = row["mw_p_loo"]
p_circ_2s = row["mw_p_circular"]

vp = axc.violinplot([gd_tau, rd_tau], positions=[1, 2], widths=0.72,
                    showmeans=False, showmedians=False, showextrema=False)
for body, c in zip(vp["bodies"], [GD, RD]):
    body.set_facecolor(c); body.set_alpha(0.75); body.set_edgecolor("none")
for i, (arr, c) in enumerate([(gd_tau, GD), (rd_tau, RD)], start=1):
    q25, med, q75 = np.percentile(arr, [25, 50, 75])
    axc.plot([i - 0.18, i + 0.18], [med, med], color="black", lw=1.4,
             zorder=4)
    axc.plot([i, i], [q25, q75], color="black", lw=2.4, alpha=0.55, zorder=3)
    rng = np.random.default_rng(7)
    jitter = rng.normal(i, 0.035, size=min(250, len(arr)))
    axc.scatter(jitter, arr[:len(jitter)], s=2.5, color=c, alpha=0.35,
                linewidths=0, zorder=2)
axc.set_xticks([1, 2])
axc.set_xticklabels([f"GD (n={len(gd_tau):,})", f"RD (n={len(rd_tau):,})"])
axc.set_ylabel("Tissue specificity tau (GTEx v11)")
axc.set_ylim(0, 1.18)
ymax = 1.06
axc.plot([1, 1, 2, 2], [ymax, ymax + 0.035, ymax + 0.035, ymax],
         color="#333333", lw=0.9)
axc.text(1.5, ymax + 0.055, f"two-sided MWU P = {p_loo:.2f} (n.s.)",
         ha="center", fontsize=8, color="#333333")
# cautionary inset: circular v7 contrast (bottom centre, below the violins)
axin = axc.inset_axes([0.295, 0.035, 0.41, 0.29])
axin.set_zorder(10)
meds = [row["gd_median_circular"], row["rd_median_circular"]]
axin.barh([1, 0], meds, height=0.52, color=[NEUT, NEUT],
          edgecolor="#999999", lw=0.6)
axin.set_yticks([1, 0]); axin.set_yticklabels(["GD", "RD"], fontsize=6.5)
axin.set_xlim(0, 1.12); axin.set_ylim(-0.55, 2.4)
axin.set_xticks([0, 0.5, 1.0])
for yi, v in zip([1, 0], meds):
    axin.text(v + 0.03, yi, f"{v:.2f}", fontsize=6.2, va="center",
              color=GREY)
# caption inside the inset (title text would overflow the inset width)
axin.text(0.56, 2.28, "v7 circular (artefact)", fontsize=5.9, color=GREY,
          ha="center", va="top")
axin.text(0.56, 1.72, f"P = {p_circ_2s:.1e} — not used", fontsize=5.9,
          color=GREY, ha="center", va="top")
axin.tick_params(axis="x", labelsize=6)
for s in ("top", "right"):
    axin.spines[s].set_visible(False)
axin.set_facecolor("#f7f7f7")
axin.patch.set_edgecolor("#aaaaaa")
axin.patch.set_linewidth(0.9)

summary.append(f"c | tau LOO-tau: GD med {np.median(gd_tau):.3f} (n={len(gd_tau)}) vs "
               f"RD med {np.median(rd_tau):.3f} (n={len(rd_tau)}), two-sided MWU "
               f"P={p_loo:.2f} n.s.; "
               f"circular v7 two-sided P={p_circ_2s:.1e} shown greyed")

# ============================ panel d: brain vs non-brain
axd = fig.add_subplot(gs[1, 2:4])
panel_label(axd, "d")
axd.set_title("Brain vs non-brain tissues, LOO-tau (mean OR)",
              loc="center", pad=6, fontsize=9)

# S2: use unified LOO-tau variant; no MWU p-value (n=2 granularity-limited)
gd_bvn = variants["loo_tau"]["brain_vs_nonbrain_gd"]
rd_bvn = variants["loo_tau"]["brain_vs_nonbrain_rd"]
groups = ["GD", "RD"]
brain_vals = [gd_bvn["brain_mean"], rd_bvn["brain_mean"]]
nonb_vals = [gd_bvn["nonbrain_mean"], rd_bvn["nonbrain_mean"]]
x = np.arange(2)
w = 0.34
b1 = axd.bar(x - w / 2, brain_vals, width=w, color=[GDD, NHI],
             label="Brain tissues (n=2)")
b2 = axd.bar(x + w / 2, nonb_vals, width=w, color=[GD, RD], alpha=0.55,
             label="Non-brain tissues (n=28)")
for xi, v in zip(x - w / 2, brain_vals):
    axd.text(xi, v + 0.03, f"{v:.2f}", ha="center", fontsize=7.6,
             color="#333333")
for xi, v in zip(x + w / 2, nonb_vals):
    axd.text(xi, v + 0.03, f"{v:.2f}", ha="center", fontsize=7.6,
             color="#333333")
axd.axhline(1.0, color=GREY, ls="--", lw=0.9)
axd.set_xticks(x); axd.set_xticklabels(groups, fontsize=9)
axd.set_ylabel("Mean enrichment OR")
axd.set_ylim(0, 2.3)
axd.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, 0.78),
           handlelength=1.2, borderaxespad=0.1, fontsize=7.2)
axd.text(0.5, 0.995,
         "direction only — n=2 brain tissues, test granularity limited "
         "(min two-sided P = 0.0046);\nn.s. under LOO-tau "
         f"(permutation p = {rd_bvn['perm_p_one_sided']:.3f}); "
         "not used as evidence",
         transform=axd.transAxes, ha="center", va="top", fontsize=6.9,
         color="#555555", style="italic")

summary.append(f"d | brain vs non-brain mean OR (LOO-tau): GD "
               f"{brain_vals[0]:.2f} vs {nonb_vals[0]:.2f}; RD "
               f"{brain_vals[1]:.2f} vs {nonb_vals[1]:.2f} "
               f"(perm p={rd_bvn['perm_p_one_sided']:.3f}; no MWU shown, "
               "n=2 granularity-limited)")

# ============================ panel e: per-tissue specificity scatter
axe = fig.add_subplot(gs[1, 4:6])
panel_label(axe, "e")
axe.set_title("Per-tissue specificity: GD vs RD genes",
              loc="center", pad=6, fontsize=9)

sig = spec["fdr_sig"].values
axe.scatter(spec.loc[~sig, "gd_mean_spec"], spec.loc[~sig, "rd_mean_spec"],
            s=26, color=NEUT, edgecolor="#8a8a8a", lw=0.5, zorder=3,
            label="n.s.")
# S5: Brain is the ONLY FDR-significant tissue (RD>GD) — positive evidence
brain_pt = spec[spec["tissue"] == "Brain"].iloc[0]
axe.scatter(spec.loc[sig, "gd_mean_spec"], spec.loc[sig, "rd_mean_spec"],
            s=90, color=RD, edgecolor=NHI, lw=1.4, zorder=5, marker="o",
            label=(f"Brain — RD > GD, FDR={brain_pt['fdr']:.3f}\n"
                   "(only FDR-sig tissue)"))
lims = (0.10, 0.56)
axe.plot(lims, lims, color=GREY, ls="--", lw=0.9, zorder=1)
axe.text(0.115, 0.40, "dashed: RD = GD", fontsize=6.8, color=GREY,
         ha="left", va="center", style="italic")
# S7: correct tissue names (no spinal cord / cerebellum in the 30 tissues);
# Nerve = tibial nerve, peripheral
offsets = {"Testis": (0.008, -0.026, "left"),
           "Nerve (periph.)": (0.012, 0.006, "left"),
           "Brain": (0.013, 0.008, "left"), "Ovary": (-0.012, 0.014, "right")}
for tissue, label in [("Testis", "Testis"), ("Nerve", "Nerve (periph.)"),
                      ("Brain", "Brain"), ("Ovary", "Ovary")]:
    r = spec[spec["tissue"] == tissue]
    if len(r) == 0:
        continue
    r = r.iloc[0]
    dx, dy, ha = offsets[label]
    weight = "bold" if tissue == "Brain" else "normal"
    axe.annotate(label, (r["gd_mean_spec"], r["rd_mean_spec"]),
                 xytext=(r["gd_mean_spec"] + dx, r["rd_mean_spec"] + dy),
                 fontsize=7.4, color="#222222", ha=ha, va="center",
                 fontweight=weight)
axe.set_xlim(*lims); axe.set_ylim(*lims)
axe.set_xlabel("GD genes: mean specificity in tissue")
axe.set_ylabel("RD genes: mean specificity in tissue")
axe.legend(frameon=False, loc="upper left", handletextpad=0.2, fontsize=7.0)

n_above = int((spec["rd_mean_spec"] > spec["gd_mean_spec"]).sum())
summary.append(f"e | S5: Brain is the ONLY FDR-sig tissue (RD>GD, "
               f"p={brain_pt['mw_p']:.1e}, FDR={brain_pt['fdr']:.3f}, "
               f"Cliff's d={brain_pt['cliffs_delta']:.3f}) — shown as "
               "positive evidence; labels Testis/Nerve(periph.)/Brain/Ovary")

# ----------------------------------------------------------------- save ----
fig.savefig(os.path.join(OUTDIR, "Figure2_bulk_tissue.png"), dpi=300)
fig.savefig(os.path.join(OUTDIR, "Figure2_bulk_tissue.pdf"),
            bbox_inches="tight")
plt.close(fig)

print("Figure 2 written to", OUTDIR)
for s in summary:
    print("  " + s)
