#!/usr/bin/env python3
"""
Phase 8 Component 2: Tissue-Specific Enrichment Analysis

For each GTEx V11 tissue:
- Define high-expression genes (top 25% or TPM > 1.0)
- Fisher's exact test: GD vs RD enrichment
- Also test GD vs background, RD vs background
- BH FDR correction + OR 95% CI

Also: Brain vs non-brain meta-comparison

Run in WSL: python3 phase8_layer1_tissue_enrichment.py
"""
import csv
import json
import os
import numpy as np
from scipy import stats
from collections import Counter
from pathlib import Path

BASE = Path("/mnt/d/hs_gene_project")
MASTER_CSV = BASE / "results/phase8_tissue_analysis/data/phase8_master_table.csv"
OUTPUT_DIR = BASE / "results/phase8_tissue_analysis"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 70)
print("Phase 8 Component 2: Tissue-Specific Enrichment Analysis")
print("=" * 70)

# ============================================
# Load master table
# ============================================
print("\n--- Loading master table ---")

genes = []
with open(MASTER_CSV) as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames
    for row in reader:
        genes.append(row)

print(f"  Loaded {len(genes)} genes, {len(fieldnames)} columns")

# Extract GTEx tissue columns
gtex_tissues = [h.replace("gtex_", "") for h in fieldnames if h.startswith("gtex_") and h != "gtex_v11_tau"]
print(f"  GTEx V11 tissues: {len(gtex_tissues)}")

# Brain tissues (GTEx)
BRAIN_TISSUES_GT = ["Brain", "Nerve"]
brain_tissue_cols = [t for t in gtex_tissues if any(b in t for b in BRAIN_TISSUES_GT)]
non_brain_tissue_cols = [t for t in gtex_tissues if t not in brain_tissue_cols]
print(f"  Brain tissues: {brain_tissue_cols}")
print(f"  Non-brain tissues: {len(non_brain_tissue_cols)}")

# ============================================
# Helper functions
# ============================================

def fisher_exact_test(group_genes, bg_genes, feature_genes, alternative="greater"):
    """Fisher's exact test for enrichment.
    
    group_genes: set of genes in the test group (e.g., GD)
    bg_genes: set of all genes (background)
    feature_genes: set of genes with the feature (e.g., high expression in tissue)
    """
    # 2x2 table
    #           Has feature | No feature
    # In group  |    a      |    b
    # Not group |    c      |    d
    
    a = len(group_genes & feature_genes)  # in group AND has feature
    b = len(group_genes - feature_genes)  # in group AND no feature
    c = len(feature_genes - group_genes)  # not in group AND has feature
    d = len(bg_genes - group_genes - feature_genes)  # not in group AND no feature
    
    table = [[a, b], [c, d]]
    odds, p = stats.fisher_exact(table, alternative=alternative)
    
    # 95% CI for OR
    if a > 0 and b > 0 and c > 0 and d > 0:
        ci_low = np.exp(np.log(odds) - 1.96 * np.sqrt(1/a + 1/b + 1/c + 1/d))
        ci_high = np.exp(np.log(odds) + 1.96 * np.sqrt(1/a + 1/b + 1/c + 1/d))
    else:
        ci_low, ci_high = None, None
    
    return {
        "OR": odds,
        "p": p,
        "ci_low": ci_low,
        "ci_high": ci_high,
        "a": a, "b": b, "c": c, "d": d,
        "n_group": len(group_genes),
        "n_feature": len(feature_genes),
        "pct_group_with_feature": 100 * a / max(len(group_genes), 1),
        "pct_bg_with_feature": 100 * (a + c) / max(len(bg_genes), 1),
    }

def bh_fdr(pvals):
    """Benjamini-Hochberg FDR correction."""
    pvals = np.array(pvals)
    n = len(pvals)
    sorted_idx = np.argsort(pvals)
    fdr = np.zeros(n)
    fdr[sorted_idx[-1]] = pvals[sorted_idx[-1]]
    for i in range(n - 2, -1, -1):
        rank = i + 1
        fdr[sorted_idx[i]] = min(pvals[sorted_idx[i]] * n / rank, fdr[sorted_idx[i + 1]])
    return fdr

# ============================================
# Define gene groups
# ============================================
print("\n--- Defining gene groups ---")

# Primary classification (classification_v7)
gd_genes = set(g["gene_id"] for g in genes if g["classification_v7"] == "gene-driven")
rd_genes = set(g["gene_id"] for g in genes if g["classification_v7"] == "regulation-driven")
dual_genes = set(g["gene_id"] for g in genes if g["classification_v7"] == "dual-driven")
gd_relaxed_genes = set(g["gene_id"] for g in genes if g["classification_v7"] == "gene-driven (relaxed)")
neutral_genes = set(g["gene_id"] for g in genes if g["classification_v7"] == "neutral")
all_genes = set(g["gene_id"] for g in genes)

# LOO-brain classification
gd_loo_brain = set(g["gene_id"] for g in genes if g.get("classification_loo_brain") == "gene-driven")
rd_loo_brain = set(g["gene_id"] for g in genes if g.get("classification_loo_brain") == "regulation-driven")

print(f"  GD (v7): {len(gd_genes)}")
print(f"  RD (v7): {len(rd_genes)}")
print(f"  Dual (v7): {len(dual_genes)}")
print(f"  GD-relaxed (v7): {len(gd_relaxed_genes)}")
print(f"  Neutral (v7): {len(neutral_genes)}")
print(f"  GD (LOO-brain): {len(gd_loo_brain)}")
print(f"  RD (LOO-brain): {len(rd_loo_brain)}")

# ============================================
# Layer 1: Tissue enrichment analysis
# ============================================
print("\n--- Layer 1: Tissue enrichment (GTEx V11, 30 tissues) ---")

results_layer1 = []

for tissue in gtex_tissues:
    col = f"gtex_{tissue}"
    
    # Get TPM values for this tissue
    tpm_values = np.array([float(g.get(col) or 0) for g in genes])
    
    # Define high-expression genes: top 25% OR TPM > 1.0, whichever gives more genes
    top25_threshold = np.percentile(tpm_values[tpm_values > 0], 75) if np.any(tpm_values > 0) else 0
    high_expr_tpm = set(g["gene_id"] for g, t in zip(genes, tpm_values) if t >= max(top25_threshold, 1.0))
    
    # Also define "expressed" genes (TPM > 0.1)
    expressed = set(g["gene_id"] for g, t in zip(genes, tpm_values) if t > 0.1)
    
    # Tests using classification_v7 (tissue TPM not in RDS, safe for all except brain tissues)
    # For brain tissues, use LOO-brain classification
    is_brain = tissue in brain_tissue_cols
    
    if is_brain:
        gd_set = gd_loo_brain
        rd_set = rd_loo_brain
        cls_label = "LOO-brain"
    else:
        gd_set = gd_genes
        rd_set = rd_genes
        cls_label = "v7"
    
    # Test 1: GD enrichment in high-expression genes
    res_gd = fisher_exact_test(gd_set, all_genes, high_expr_tpm)
    
    # Test 2: RD enrichment in high-expression genes
    res_rd = fisher_exact_test(rd_set, all_genes, high_expr_tpm)
    
    # Test 3: GD vs RD (GD in high-expr, RD as background)
    res_gd_vs_rd = fisher_exact_test(gd_set, gd_set | rd_set, high_expr_tpm)
    
    results_layer1.append({
        "tissue": tissue,
        "is_brain": is_brain,
        "classification_used": cls_label,
        "n_high_expr": len(high_expr_tpm),
        "n_expressed": len(expressed),
        "mean_tpm": np.mean(tpm_values),
        "median_tpm": np.median(tpm_values),
        
        "gd_OR": res_gd["OR"],
        "gd_p": res_gd["p"],
        "gd_ci_low": res_gd["ci_low"],
        "gd_ci_high": res_gd["ci_high"],
        "gd_n_high": res_gd["a"],
        "gd_pct": res_gd["pct_group_with_feature"],
        
        "rd_OR": res_rd["OR"],
        "rd_p": res_rd["p"],
        "rd_ci_low": res_rd["ci_low"],
        "rd_ci_high": res_rd["ci_high"],
        "rd_n_high": res_rd["a"],
        "rd_pct": res_rd["pct_group_with_feature"],
        
        "gd_vs_rd_OR": res_gd_vs_rd["OR"],
        "gd_vs_rd_p": res_gd_vs_rd["p"],
    })
    
    status = "  "
    if res_gd["p"] < 0.05:
        status += f"GD* "
    else:
        status += f"GD  "
    if res_rd["p"] < 0.05:
        status += f"RD* "
    else:
        status += f"RD  "
    
    print(f"  {tissue:25s} n_high={len(high_expr_tpm):4d} "
          f"GD:OR={res_gd['OR']:.2f},p={res_gd['p']:.2e} "
          f"RD:OR={res_rd['OR']:.2f},p={res_rd['p']:.2e} [{cls_label}]")

# BH FDR correction
gd_pvals = [r["gd_p"] for r in results_layer1]
rd_pvals = [r["rd_p"] for r in results_layer1]
gd_fdr = bh_fdr(gd_pvals)
rd_fdr = bh_fdr(rd_pvals)

for i, r in enumerate(results_layer1):
    r["gd_fdr"] = gd_fdr[i]
    r["rd_fdr"] = rd_fdr[i]
    r["gd_fdr_sig"] = gd_fdr[i] < 0.05
    r["rd_fdr_sig"] = rd_fdr[i] < 0.05

# Print FDR-significant results
print(f"\n  FDR-significant results:")
print(f"  {'Tissue':25s} {'GD OR':>8s} {'GD FDR':>10s} {'RD OR':>8s} {'RD FDR':>10s}")
for r in results_layer1:
    if r["gd_fdr_sig"] or r["rd_fdr_sig"]:
        gd_sig = "*" if r["gd_fdr_sig"] else " "
        rd_sig = "*" if r["rd_fdr_sig"] else " "
        print(f"  {r['tissue']:25s} {r['gd_OR']:8.2f}{gd_sig} {r['gd_fdr']:10.2e} "
              f"{r['rd_OR']:8.2f}{rd_sig} {r['rd_fdr']:10.2e}")

# Write Layer 1 output
layer1_csv = OUTPUT_DIR / "layer1_tissue_enrichment.csv"
layer1_fields = ["tissue", "is_brain", "classification_used", "n_high_expr", "n_expressed",
                 "mean_tpm", "median_tpm",
                 "gd_OR", "gd_p", "gd_fdr", "gd_fdr_sig", "gd_ci_low", "gd_ci_high",
                 "gd_n_high", "gd_pct",
                 "rd_OR", "rd_p", "rd_fdr", "rd_fdr_sig", "rd_ci_low", "rd_ci_high",
                 "rd_n_high", "rd_pct",
                 "gd_vs_rd_OR", "gd_vs_rd_p"]

with open(layer1_csv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=layer1_fields, extrasaction="ignore")
    writer.writeheader()
    for r in results_layer1:
        writer.writerow(r)

print(f"\n  Output: {layer1_csv}")

# ============================================
# Brain vs Non-Brain meta-comparison
# ============================================
print(f"\n--- Brain vs Non-Brain Meta-Comparison ---")

brain_results = [r for r in results_layer1 if r["is_brain"]]
non_brain_results = [r for r in results_layer1 if not r["is_brain"]]

# Mean OR
brain_gd_ors = [r["gd_OR"] for r in brain_results if not np.isnan(r["gd_OR"])]
brain_rd_ors = [r["rd_OR"] for r in brain_results if not np.isnan(r["rd_OR"])]
non_brain_gd_ors = [r["gd_OR"] for r in non_brain_results if not np.isnan(r["gd_OR"])]
non_brain_rd_ors = [r["rd_OR"] for r in non_brain_results if not np.isnan(r["rd_OR"])]

print(f"  Brain tissues ({len(brain_results)}):")
print(f"    Mean GD OR: {np.mean(brain_gd_ors):.3f}")
print(f"    Mean RD OR: {np.mean(brain_rd_ors):.3f}")
print(f"  Non-brain tissues ({len(non_brain_results)}):")
print(f"    Mean GD OR: {np.mean(non_brain_gd_ors):.3f}")
print(f"    Mean RD OR: {np.mean(non_brain_rd_ors):.3f}")

# Mann-Whitney U: brain vs non-brain OR distributions
if len(brain_gd_ors) > 1 and len(non_brain_gd_ors) > 1:
    u_gd, p_gd_brain = stats.mannwhitneyu(brain_gd_ors, non_brain_gd_ors, alternative="greater")
    print(f"  GD OR brain > non-brain: U={u_gd:.0f}, p={p_gd_brain:.4e}")

if len(brain_rd_ors) > 1 and len(non_brain_rd_ors) > 1:
    u_rd, p_rd_brain = stats.mannwhitneyu(brain_rd_ors, non_brain_rd_ors, alternative="greater")
    print(f"  RD OR brain > non-brain: U={u_rd:.0f}, p={p_rd_brain:.4e}")

# Stouffer's method for combining p-values
def stouffer_combine(pvals):
    """Stouffer's method: combine p-values using z-scores."""
    pvals = np.array(pvals)
    z_scores = stats.norm.ppf(1 - pvals)
    combined_z = np.sum(z_scores) / np.sqrt(len(pvals))
    combined_p = 1 - stats.norm.cdf(combined_z)
    return combined_p, combined_z

brain_gd_pvals = [r["gd_p"] for r in brain_results]
brain_rd_pvals = [r["rd_p"] for r in brain_results]
non_brain_gd_pvals = [r["gd_p"] for r in non_brain_results]
non_brain_rd_pvals = [r["rd_p"] for r in non_brain_results]

if len(brain_gd_pvals) > 1:
    st_gd, z_gd = stouffer_combine(brain_gd_pvals)
    print(f"  Stouffer GD (brain): z={z_gd:.2f}, p={st_gd:.2e}")

if len(brain_rd_pvals) > 1:
    st_rd, z_rd = stouffer_combine(brain_rd_pvals)
    print(f"  Stouffer RD (brain): z={z_rd:.2f}, p={st_rd:.2e}")

if len(non_brain_gd_pvals) > 1:
    st_gd_nb, z_gd_nb = stouffer_combine(non_brain_gd_pvals)
    print(f"  Stouffer GD (non-brain): z={z_gd_nb:.2f}, p={st_gd_nb:.2e}")

if len(non_brain_rd_pvals) > 1:
    st_rd_nb, z_rd_nb = stouffer_combine(non_brain_rd_pvals)
    print(f"  Stouffer RD (non-brain): z={z_rd_nb:.2f}, p={st_rd_nb:.2e}")

# Write brain vs non-brain summary
brain_meta_csv = OUTPUT_DIR / "layer1_brain_vs_nonbrain.csv"
with open(brain_meta_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["metric", "brain_mean", "brain_n", "non_brain_mean", "non_brain_n",
                     "mannwhitney_u", "mannwhitney_p", "stouffer_z_brain", "stouffer_p_brain",
                     "stouffer_z_nonbrain", "stouffer_p_nonbrain"])
    writer.writerow(["GD_OR", np.mean(brain_gd_ors), len(brain_gd_ors),
                     np.mean(non_brain_gd_ors), len(non_brain_gd_ors),
                     u_gd if 'u_gd' in dir() else "", p_gd_brain if 'p_gd_brain' in dir() else "",
                     z_gd if 'z_gd' in dir() else "", st_gd if 'st_gd' in dir() else "",
                     z_gd_nb if 'z_gd_nb' in dir() else "", st_gd_nb if 'st_gd_nb' in dir() else ""])
    writer.writerow(["RD_OR", np.mean(brain_rd_ors), len(brain_rd_ors),
                     np.mean(non_brain_rd_ors), len(non_brain_rd_ors),
                     u_rd if 'u_rd' in dir() else "", p_rd_brain if 'p_rd_brain' in dir() else "",
                     z_rd if 'z_rd' in dir() else "", st_rd if 'st_rd' in dir() else "",
                     z_rd_nb if 'z_rd_nb' in dir() else "", st_rd_nb if 'st_rd_nb' in dir() else ""])

print(f"\n  Output: {brain_meta_csv}")

# ============================================
# Summary JSON
# ============================================
summary = {
    "total_genes": len(genes),
    "n_tissues": len(gtex_tissues),
    "n_brain_tissues": len(brain_tissue_cols),
    "n_non_brain_tissues": len(non_brain_tissue_cols),
    "classification_v7": {
        "gene_driven": len(gd_genes),
        "regulation_driven": len(rd_genes),
        "dual_driven": len(dual_genes),
        "gene_driven_relaxed": len(gd_relaxed_genes),
        "neutral": len(neutral_genes),
    },
    "gd_fdr_sig_tissues": [r["tissue"] for r in results_layer1 if r["gd_fdr_sig"]],
    "rd_fdr_sig_tissues": [r["tissue"] for r in results_layer1 if r["rd_fdr_sig"]],
    "brain_vs_nonbrain": {
        "brain_gd_mean_OR": float(np.mean(brain_gd_ors)) if brain_gd_ors else None,
        "brain_rd_mean_OR": float(np.mean(brain_rd_ors)) if brain_rd_ors else None,
        "non_brain_gd_mean_OR": float(np.mean(non_brain_gd_ors)) if non_brain_gd_ors else None,
        "non_brain_rd_mean_OR": float(np.mean(non_brain_rd_ors)) if non_brain_rd_ors else None,
    },
}

summary_file = OUTPUT_DIR / "layer1_summary.json"
with open(summary_file, "w") as f:
    json.dump(summary, f, indent=2)

print(f"\n  Summary: {summary_file}")

print("\n" + "=" * 70)
print("Layer 1 Tissue Enrichment Analysis Complete!")
print("=" * 70)
