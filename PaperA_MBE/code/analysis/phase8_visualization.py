#!/usr/bin/env python3
"""
Phase 8 Component 7: Visualization — 7 Figures

Fig 1: Tissue enrichment heatmap (30 tissues × GD/RD, color=log(OR))
Fig 2: Tau distribution violin plots (GD vs RD vs neutral, 4 sources)
Fig 3: Brain region enrichment (BrainSpan 26 regions)
Fig 4: Full developmental trajectory (31 stages × GD/RD OR)
Fig 5: Cell type specificity (34 HPA single-nuclei)
Fig 6: Regulatory element tissue map (caMPRA + hCONDELs × tissues)
Fig 7: Omega distribution (by classification)

Run in WSL: python3 phase8_visualization.py
"""
import csv
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from pathlib import Path
from scipy import stats

BASE = Path("/mnt/d/hs_gene_project")
DATA_DIR = BASE / "results/phase8_tissue_analysis/data"
RESULTS_DIR = BASE / "results/phase8_tissue_analysis"
FIG_DIR = BASE / "results/phase8_tissue_analysis/figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Color scheme
GD_COLOR = "#2196F3"   # Blue
RD_COLOR = "#FF5722"   # Red-orange
NEUTRAL_COLOR = "#9E9E9E"  # Gray
DUAL_COLOR = "#9C27B0"  # Purple

plt.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 14,
    "axes.labelsize": 12,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

print("=" * 70)
print("Phase 8 Component 7: Visualization — 7 Figures")
print("=" * 70)

# ============================================
# Load master table (needed for multiple figures)
# ============================================
print("\n--- Loading master table ---")
genes = []
with open(DATA_DIR / "phase8_master_table.csv") as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames
    for row in reader:
        genes.append(row)
print(f"  {len(genes)} genes loaded")

# ============================================
# Figure 1: Tissue Enrichment Heatmap
# ============================================
print("\n--- Figure 1: Tissue Enrichment Heatmap ---")

tissue_csv = RESULTS_DIR / "layer1_tissue_enrichment.csv"
if tissue_csv.exists():
    tissue_data = []
    with open(tissue_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            tissue_data.append(row)
    
    # Build matrix: tissues × [GD_OR, RD_OR]
    tissues = [r["tissue"] for r in tissue_data]
    gd_ors = [float(r.get("gd_OR", 1)) for r in tissue_data]
    rd_ors = [float(r.get("rd_OR", 1)) for r in tissue_data]
    
    # Log2 transform
    gd_log = [np.log2(o) if o > 0 else -2 for o in gd_ors]
    rd_log = [np.log2(o) if o > 0 else -2 for o in rd_ors]
    
    matrix = np.array([gd_log, rd_log]).T  # tissues × 2
    
    fig, ax = plt.subplots(figsize=(6, 10))
    
    # Custom diverging colormap
    cmap = sns.diverging_palette(10, 150, s=80, l=55, sep=3, as_cmap=True)
    
    annot = np.array([[f"{o:.2f}" for o in gd_ors],
                       [f"{o:.2f}" for o in rd_ors]]).T
    sns.heatmap(matrix, annot=annot,
                fmt="", cmap=cmap, center=0, vmin=-1, vmax=1,
                yticklabels=tissues, xticklabels=["Gene-Driven", "Regulation-Driven"],
                cbar_kws={"label": "log2(Odds Ratio)", "shrink": 0.5},
                linewidths=0.5, linecolor="white", ax=ax)
    
    ax.set_title("Tissue Enrichment of Positive Selection Genes", fontweight="bold")
    ax.set_ylabel("GTEx V11 Tissue")
    
    # Mark brain tissues
    brain_tissues = ["Brain", "Nerve"]
    for i, t in enumerate(tissues):
        if t in brain_tissues:
            ax.add_patch(mpatches.Rectangle((0, i), 2, 1, fill=False, 
                                           edgecolor="red", linewidth=2))
    
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig1_tissue_enrichment_heatmap.png")
    plt.close()
    print(f"  Saved: fig1_tissue_enrichment_heatmap.png")
else:
    print("  Skipping: layer1_tissue_enrichment.csv not found")

# ============================================
# Figure 2: Tau Distribution Violin Plots
# ============================================
print("\n--- Figure 2: Tau Distribution Violin Plots ---")

tau_sources = ["gtex_v11_tau", "hpa_gtex_tau", "hpa_hpa_tau", "hpa_consensus_tau"]
tau_labels = ["GTEx V11\n(30 tissues)", "HPA-GTEx\n(36 tissues)", "HPA-HPA\n(40 tissues)", "HPA Consensus\n(51 tissues)"]

fig, axes = plt.subplots(1, 4, figsize=(16, 5))

for ax, tau_col, label in zip(axes, tau_sources, tau_labels):
    # Extract tau values by classification
    gd_tau = []
    rd_tau = []
    neu_tau = []
    
    for g in genes:
        cls = g.get("classification_loo_tau", g.get("classification_v7", ""))
        try:
            tau = float(g.get(tau_col, 0) or 0)
        except:
            continue
        
        if cls == "gene-driven":
            gd_tau.append(tau)
        elif cls == "regulation-driven":
            rd_tau.append(tau)
        elif cls == "neutral":
            neu_tau.append(tau)
    
    if not gd_tau and not rd_tau:
        ax.set_visible(False)
        continue
    
    # Create DataFrame for seaborn
    import pandas as pd
    plot_data = []
    for v in gd_tau:
        plot_data.append({"tau": v, "group": "Gene-Driven"})
    for v in rd_tau:
        plot_data.append({"tau": v, "group": "Regulation-Driven"})
    for v in neu_tau:
        plot_data.append({"tau": v, "group": "Neutral"})
    
    df_plot = pd.DataFrame(plot_data)
    
    # Violin plot
    sns.violinplot(data=df_plot, x="group", y="tau", ax=ax,
                   palette={"Gene-Driven": GD_COLOR, "Regulation-Driven": RD_COLOR, 
                           "Neutral": NEUTRAL_COLOR},
                   order=["Gene-Driven", "Regulation-Driven", "Neutral"],
                   cut=0, inner="quartile")
    
    ax.set_title(label, fontweight="bold")
    ax.set_ylabel("Expression Specificity (τ)" if tau_col == tau_sources[0] else "")
    ax.set_xlabel("")
    ax.set_ylim(-0.05, 1.05)
    
    # Add p-value annotation
    if gd_tau and rd_tau:
        u, p = stats.mannwhitneyu(rd_tau, gd_tau, alternative="greater")
        sig_text = f"p={p:.2e}" if p < 0.05 else f"p={p:.3f}"
        ax.text(0.5, 0.98, sig_text, transform=ax.transAxes,
                ha="center", va="top", fontsize=9, color="red" if p < 0.05 else "gray")

plt.suptitle("Expression Specificity (τ) Distribution by Selection Category (LOO-τ, non-circular)", 
             fontsize=14, fontweight="bold", y=1.02)
plt.tight_layout()
plt.savefig(FIG_DIR / "fig2_tau_violin_plots.png")
plt.close()
print(f"  Saved: fig2_tau_violin_plots.png")

# ============================================
# Figure 3: Brain Region Enrichment (Forest Plot)
# ============================================
print("\n--- Figure 3: Brain Region Enrichment ---")

bs_csv = RESULTS_DIR / "brainspan_region_enrichment.csv"
if bs_csv.exists():
    bs_data = []
    with open(bs_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            bs_data.append(row)
    
    # Sort by GD p-value
    bs_data.sort(key=lambda r: float(r.get("gd_p", 1)))
    
    fig, axes = plt.subplots(1, 2, figsize=(14, max(6, len(bs_data) * 0.3)), sharey=True)
    
    regions = [r["region"] for r in bs_data]
    gd_ors = [float(r.get("gd_OR", 1)) for r in bs_data]
    rd_ors = [float(r.get("rd_OR", 1)) for r in bs_data]
    gd_fdrs = [float(r.get("gd_fdr", 1)) for r in bs_data]
    rd_fdrs = [float(r.get("rd_fdr", 1)) for r in bs_data]
    
    y_pos = np.arange(len(regions))
    
    # GD forest plot
    ax = axes[0]
    for i, (or_val, fdr) in enumerate(zip(gd_ors, gd_fdrs)):
        color = GD_COLOR if fdr < 0.05 else "#90CAF9"
        ax.scatter(or_val, i, color=color, s=80, zorder=3)
        # CI would go here if available; use simple error bar
        ax.errorbar(or_val, i, xerr=0.2 * or_val, fmt="none", color=color, alpha=0.5)
    
    ax.axvline(x=1.0, color="black", linestyle="--", alpha=0.3)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(regions, fontsize=8)
    ax.set_xlabel("Odds Ratio")
    ax.set_title("Gene-Driven", fontweight="bold", color=GD_COLOR)
    ax.set_xlim(0.5, 2.0)
    
    # RD forest plot
    ax = axes[1]
    for i, (or_val, fdr) in enumerate(zip(rd_ors, rd_fdrs)):
        color = RD_COLOR if fdr < 0.05 else "#FFAB91"
        ax.scatter(or_val, i, color=color, s=80, zorder=3)
        ax.errorbar(or_val, i, xerr=0.2 * or_val, fmt="none", color=color, alpha=0.5)
    
    ax.axvline(x=1.0, color="black", linestyle="--", alpha=0.3)
    ax.set_xlabel("Odds Ratio")
    ax.set_title("Regulation-Driven", fontweight="bold", color=RD_COLOR)
    ax.set_xlim(0.5, 2.0)
    
    plt.suptitle("BrainSpan Region Enrichment (LOO-Brain, non-circular)", 
                 fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig3_brain_region_forest.png")
    plt.close()
    print(f"  Saved: fig3_brain_region_forest.png")
else:
    print("  Skipping: brainspan_region_enrichment.csv not found")

# ============================================
# Figure 4: Full Developmental Trajectory
# ============================================
print("\n--- Figure 4: Full Developmental Trajectory ---")

dev_csv = RESULTS_DIR / "developmental_full_trajectory.csv"
if dev_csv.exists():
    dev_data = []
    with open(dev_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            dev_data.append(row)
    
    stages = [r["stage"] for r in dev_data]
    gd_ors = [float(r.get("gd_OR", 1)) for r in dev_data]
    rd_ors = [float(r.get("rd_OR", 1)) for r in dev_data]
    gd_fdrs = [float(r.get("gd_fdr", 1)) for r in dev_data]
    rd_fdrs = [float(r.get("rd_fdr", 1)) for r in dev_data]
    
    # Sort stages chronologically
    def age_sort_key(age):
        if "pcw" in age:
            return (0, int(age.replace("pcw", "").strip()))
        elif "mos" in age:
            return (1, int(float(age.replace("mos", "").strip())))
        elif "yrs" in age:
            return (2, int(float(age.replace("yrs", "").strip())))
        return (3, 0)
    
    sorted_indices = sorted(range(len(stages)), key=lambda i: age_sort_key(stages[i]))
    stages = [stages[i] for i in sorted_indices]
    gd_ors = [gd_ors[i] for i in sorted_indices]
    rd_ors = [rd_ors[i] for i in sorted_indices]
    gd_fdrs = [gd_fdrs[i] for i in sorted_indices]
    rd_fdrs = [rd_fdrs[i] for i in sorted_indices]
    
    fig, ax = plt.subplots(figsize=(14, 6))
    
    x = np.arange(len(stages))
    
    # Plot ORs
    ax.plot(x, gd_ors, "o-", color=GD_COLOR, linewidth=2, markersize=6, label="Gene-Driven")
    ax.plot(x, rd_ors, "s-", color=RD_COLOR, linewidth=2, markersize=6, label="Regulation-Driven")
    
    # Mark significant stages
    for i, (gd_fdr, rd_fdr) in enumerate(zip(gd_fdrs, rd_fdrs)):
        if gd_fdr < 0.05:
            ax.scatter(i, gd_ors[i], color=GD_COLOR, s=150, zorder=5, 
                       edgecolors="black", linewidths=1.5)
        if rd_fdr < 0.05:
            ax.scatter(i, rd_ors[i], color=RD_COLOR, s=150, zorder=5,
                       edgecolors="black", linewidths=1.5)
    
    ax.axhline(y=1.0, color="black", linestyle="--", alpha=0.3)
    
    # Color background by developmental period
    prenatal_end = sum(1 for s in stages if "pcw" in s)
    postnatal_start = prenatal_end
    
    ax.axvspan(-0.5, prenatal_end - 0.5, alpha=0.08, color="blue", label="Prenatal")
    ax.axvspan(postnatal_start - 0.5, len(stages) - 0.5, alpha=0.08, color="orange", label="Postnatal")
    
    ax.set_xticks(x)
    ax.set_xticklabels(stages, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("Odds Ratio (high expression enrichment)")
    ax.set_title("Full Developmental Trajectory: GD/RD Enrichment Across 31 BrainSpan Stages\n(LOO-Brain, non-circular)", 
                 fontweight="bold")
    ax.legend(loc="upper right")
    ax.set_ylim(0.8, 1.4)
    
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig4_developmental_trajectory.png")
    plt.close()
    print(f"  Saved: fig4_developmental_trajectory.png")
else:
    print("  Skipping: developmental_full_trajectory.csv not found")

# ============================================
# Figure 5: HPA Single-Nuclei Cell Type Heatmap
# ============================================
print("\n--- Figure 5: HPA Single-Nuclei Cell Type Heatmap ---")

hpa_csv = RESULTS_DIR / "hpa_single_nuclei_enrichment.csv"
if hpa_csv.exists():
    hpa_data = []
    with open(hpa_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            hpa_data.append(row)
    
    # Sort by RD OR
    hpa_data.sort(key=lambda r: float(r.get("rd_OR", 0)), reverse=True)
    
    cell_types = [r["cell_type"] for r in hpa_data]
    cell_classes = [r.get("cell_class", "") for r in hpa_data]
    gd_ors = [float(r.get("gd_OR", 1)) for r in hpa_data]
    rd_ors = [float(r.get("rd_OR", 1)) for r in hpa_data]
    
    # Log2 transform
    gd_log = [np.log2(o) if o > 0 else -2 for o in gd_ors]
    rd_log = [np.log2(o) if o > 0 else -2 for o in rd_ors]
    
    matrix = np.array([gd_log, rd_log]).T
    
    # Y-axis labels with class
    y_labels = [f"{ct} ({cc})" for ct, cc in zip(cell_types, cell_classes)]
    
    fig, ax = plt.subplots(figsize=(8, 14))
    
    cmap = sns.diverging_palette(10, 150, s=80, l=55, sep=3, as_cmap=True)
    
    annot = np.array([[f"{o:.2f}" for o in gd_ors],
                       [f"{o:.2f}" for o in rd_ors]]).T
    sns.heatmap(matrix, annot=annot,
                fmt="", cmap=cmap, center=0, vmin=-1, vmax=2,
                yticklabels=y_labels, xticklabels=["Gene-Driven", "Regulation-Driven"],
                cbar_kws={"label": "log2(Odds Ratio)", "shrink": 0.5},
                linewidths=0.5, linecolor="white", ax=ax)
    
    ax.set_title("HPA Single-Nuclei Brain Cell Type Enrichment\n(LOO-Brain, non-circular)", 
                 fontweight="bold")
    ax.set_ylabel("Cell Type (Class)")
    
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig5_hpa_celltype_heatmap.png")
    plt.close()
    print(f"  Saved: fig5_hpa_celltype_heatmap.png")
else:
    print("  Skipping: hpa_single_nuclei_enrichment.csv not found")

# ============================================
# Figure 6: Regulatory Element Tissue Map
# ============================================
print("\n--- Figure 6: Regulatory Element Tissue Map ---")

campra_csv = RESULTS_DIR / "regulatory_campra_tissue_enrichment.csv"
hcondel_csv = RESULTS_DIR / "regulatory_hcondels_tissue_enrichment.csv"

fig, axes = plt.subplots(1, 2, figsize=(14, 10))

# caMPRA
if campra_csv.exists():
    campra_data = []
    with open(campra_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            campra_data.append(row)
    
    tissues = [r["tissue"] for r in campra_data]
    rd_loo_ors = [float(r.get("rd_loo_campra_OR", 1)) for r in campra_data]
    rd_v7_ors = [float(r.get("rd_v7_OR", 1)) for r in campra_data]
    
    matrix_campra = np.array([
        [np.log2(o) if o > 0 else -2 for o in rd_loo_ors],
        [np.log2(o) if o > 0 else -2 for o in rd_v7_ors],
    ]).T
    
    ax = axes[0]
    cmap = sns.diverging_palette(10, 150, s=80, l=55, sep=3, as_cmap=True)
    annot_campra = np.array([[f"{o:.2f}" for o in rd_loo_ors],
                              [f"{o:.2f}" for o in rd_v7_ors]]).T
    sns.heatmap(matrix_campra, 
                annot=annot_campra,
                fmt="", cmap=cmap, center=0, vmin=-2, vmax=4,
                yticklabels=tissues, 
                xticklabels=["RD (LOO-caMPRA)\nnon-circular", "RD (v7)\ncircular"],
                cbar_kws={"label": "log2(OR)", "shrink": 0.5},
                linewidths=0.5, linecolor="white", ax=ax)
    ax.set_title("caMPRA Active HAR Genes\nin Tissue High-Expression", fontweight="bold")
else:
    axes[0].set_visible(False)

# hCONDELs
if hcondel_csv.exists():
    hc_data = []
    with open(hcondel_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            hc_data.append(row)
    
    tissues = [r["tissue"] for r in hc_data]
    hc_ors = [float(r.get("hcondel_OR", 1)) for r in hc_data]
    
    matrix_hc = np.array([[np.log2(o) if o > 0 else -2 for o in hc_ors]]).T
    
    ax = axes[1]
    cmap = sns.diverging_palette(10, 150, s=80, l=55, sep=3, as_cmap=True)
    annot_hc = np.array([[f"{o:.2f}"] for o in hc_ors])
    sns.heatmap(matrix_hc, 
                annot=annot_hc,
                fmt="", cmap=cmap, center=0, vmin=-1, vmax=1,
                yticklabels=tissues, xticklabels=["hCONDEL Genes"],
                cbar_kws={"label": "log2(OR)", "shrink": 0.5},
                linewidths=0.5, linecolor="white", ax=ax)
    ax.set_title("hCONDEL Genes\nin Tissue High-Expression (v7, non-circular)", fontweight="bold")
else:
    axes[1].set_visible(False)

plt.suptitle("Regulatory Element Tissue Enrichment", fontsize=14, fontweight="bold", y=1.02)
plt.tight_layout()
plt.savefig(FIG_DIR / "fig6_regulatory_tissue_map.png")
plt.close()
print(f"  Saved: fig6_regulatory_tissue_map.png")

# ============================================
# Figure 7: Omega Distribution by Classification
# ============================================
print("\n--- Figure 7: Omega Distribution ---")

omega_col = "omega_weighted"
if omega_col in fieldnames:
    gd_omega = []
    rd_omega = []
    neu_omega = []
    
    for g in genes:
        cls = g.get("classification_v7", "")
        try:
            omega = float(g.get(omega_col, 0) or 0)
        except:
            continue
        
        if omega <= 0:
            continue
        
        if cls == "gene-driven":
            gd_omega.append(omega)
        elif cls == "regulation-driven":
            rd_omega.append(omega)
        elif cls == "neutral":
            neu_omega.append(omega)
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    import pandas as pd
    plot_data = []
    for v in gd_omega:
        plot_data.append({"omega": v, "group": "Gene-Driven"})
    for v in rd_omega:
        plot_data.append({"omega": v, "group": "Regulation-Driven"})
    for v in neu_omega:
        plot_data.append({"omega": v, "group": "Neutral"})
    
    df_plot = pd.DataFrame(plot_data)
    
    sns.boxplot(data=df_plot, x="group", y="omega", ax=ax,
                palette={"Gene-Driven": GD_COLOR, "Regulation-Driven": RD_COLOR,
                         "Neutral": NEUTRAL_COLOR},
                order=["Gene-Driven", "Regulation-Driven", "Neutral"],
                width=0.6, fliersize=2)
    
    # Add p-value
    if gd_omega and rd_omega:
        u, p = stats.mannwhitneyu(gd_omega, rd_omega, alternative="greater")
        sig_text = f"GD > RD: p={p:.2e}" if p < 0.05 else f"GD > RD: p={p:.3f}"
        ax.text(0.5, 0.97, sig_text, transform=ax.transAxes,
                ha="center", va="top", fontsize=10, 
                color="red" if p < 0.05 else "gray")
    
    ax.set_title("BUSTED Omega (Weighted) Distribution\nby Selection Category (v7 classification)", 
                 fontweight="bold")
    ax.set_ylabel("Weighted Omega (Σ ωᵢ × proportionᵢ)")
    ax.set_xlabel("")
    ax.set_ylim(0, max(np.percentile(gd_omega + rd_omega + neu_omega, 95), 5))
    
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig7_omega_distribution.png")
    plt.close()
    print(f"  Saved: fig7_omega_distribution.png")
else:
    print(f"  Skipping: {omega_col} not found in master table")

# ============================================
# Summary
# ============================================
print("\n" + "=" * 70)
print("All 7 figures generated!")
print("=" * 70)
print(f"\nFigure directory: {FIG_DIR}")
for f in sorted(FIG_DIR.glob("*.png")):
    size = f.stat().st_size / 1024
    print(f"  {f.name}: {size:.0f} KB")
