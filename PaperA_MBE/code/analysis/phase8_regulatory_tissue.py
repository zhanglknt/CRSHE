#!/usr/bin/env python3
"""
Phase 8 Component 5b: caMPRA Regulatory Element Tissue Mapping

Map active HARs (from Shin et al. 2024 caMPRA) to genes,
test tissue enrichment across GTEx V11 tissues.

Uses classification_loo_campra for non-circular testing
(caMPRA is in RDS, must remove it).

Also tests all 3,171 HARs using classification_v7 (HAR not in v7 RDS — safe).

Run in WSL: python3 phase8_regulatory_tissue.py
"""
import csv
import json
import numpy as np
from scipy import stats
from collections import defaultdict
from pathlib import Path

BASE = Path("/mnt/d/hs_gene_project")
ACTIVE_HARS_CSV = BASE / "results/doan2024_campra/active_hars.csv"
ALL_HARS_CSV = BASE / "results/doan2024_campra/hars_campra_merged.csv"
MASTER_CSV = BASE / "results/phase8_tissue_analysis/data/phase8_master_table.csv"
GENE_BED = BASE / "results/phase4_conservation/gene_body.bed"
OUTPUT_DIR = BASE / "results/phase8_tissue_analysis"

print("=" * 70)
print("Phase 8 Component 5b: caMPRA Regulatory Element Tissue Mapping")
print("=" * 70)

# ============================================
# Helper functions
# ============================================
def bh_fdr(pvals):
    pvals = np.array(pvals)
    n = len(pvals)
    sorted_idx = np.argsort(pvals)
    fdr = np.zeros(n)
    fdr[sorted_idx[-1]] = pvals[sorted_idx[-1]]
    for i in range(n - 2, -1, -1):
        rank = i + 1
        fdr[sorted_idx[i]] = min(pvals[sorted_idx[i]] * n / rank, fdr[sorted_idx[i + 1]])
    return fdr

def fisher_test(group_set, bg_set, feature_set):
    a = len(group_set & feature_set)
    b = len(group_set - feature_set)
    c = len(feature_set - group_set)
    d = len(bg_set - group_set - feature_set)
    if a + b == 0 or c + d == 0:
        return 1.0, 1.0, 0, 0
    table = [[a, b], [c, d]]
    odds, p = stats.fisher_exact(table, alternative="greater")
    return odds, p, a, len(group_set)

# ============================================
# Load master table
# ============================================
print("\n--- Loading master table ---")
genes = {}
with open(MASTER_CSV) as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames
    for row in reader:
        genes[row["gene_id"]] = row

all_genes = set(genes.keys())

# Non-circular: use LOO-caMPRA for caMPRA analysis
gd_campra = set(g for g, v in genes.items() if v.get("classification_loo_campra") == "gene-driven")
rd_campra = set(g for g, v in genes.items() if v.get("classification_loo_campra") == "regulation-driven")
neutral_campra = set(g for g, v in genes.items() if v.get("classification_loo_campra") == "neutral")

# Safe: use v7 for HAR analysis (HAR not in RDS)
gd_v7 = set(g for g, v in genes.items() if v["classification_v7"] == "gene-driven")
rd_v7 = set(g for g, v in genes.items() if v["classification_v7"] == "regulation-driven")
neutral_v7 = set(g for g, v in genes.items() if v["classification_v7"] == "neutral")

print(f"  LOO-caMPRA: GD={len(gd_campra)}, RD={len(rd_campra)}, Neutral={len(neutral_campra)}")
print(f"  v7: GD={len(gd_v7)}, RD={len(rd_v7)}, Neutral={len(neutral_v7)}")

# ============================================
# Load gene coordinates
# ============================================
print("\n--- Loading gene coordinates ---")
gene_coords = {}
gene_by_chrom = defaultdict(list)
with open(GENE_BED) as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) >= 4:
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            name = parts[3]
            gene_coords[name] = (chrom, start, end)
            gene_by_chrom[chrom].append((start, end, name))

for chrom in gene_by_chrom:
    gene_by_chrom[chrom].sort()

print(f"  {len(gene_coords)} genes with coordinates")

# ============================================
# Map active HARs to genes (±50kb)
# ============================================
print("\n--- Mapping active HARs to genes (±50kb) ---")

active_hars = []
with open(ACTIVE_HARS_CSV) as f:
    reader = csv.DictReader(f)
    for row in reader:
        active_hars.append(row)

print(f"  Active HARs: {len(active_hars)}")

EXTENSION = 50000
active_har_genes = set()

for har in active_hars:
    chrom = har.get("chr_hg38", "")
    try:
        start = int(float(har.get("start_hg38", 0)))
        end = int(float(har.get("end_hg38", 0)))
    except (ValueError, TypeError):
        continue
    
    if not chrom or not start or not end:
        continue
    
    chrom_genes = gene_by_chrom.get(chrom, [])
    for g_start, g_end, gene_id in chrom_genes:
        g_start_ext = max(0, g_start - EXTENSION)
        g_end_ext = g_end + EXTENSION
        if start < g_end_ext and end > g_start_ext:
            active_har_genes.add(gene_id)
            break

active_har_genes_in_set = active_har_genes & all_genes
print(f"  Genes with active HAR overlap: {len(active_har_genes)} (total), {len(active_har_genes_in_set)} (in set)")

# Also map ALL HARs (3,171) for comparison
print("\n--- Mapping all 3,171 HARs to genes (±50kb) ---")
all_har_genes = set()
with open(ALL_HARS_CSV) as f:
    reader = csv.DictReader(f)
    for row in reader:
        chrom = row.get("chr_hg38", "")
        try:
            start = int(float(row.get("start_hg38", 0)))
            end = int(float(row.get("end_hg38", 0)))
        except (ValueError, TypeError):
            continue
        
        if not chrom or not start or not end:
            continue
        
        chrom_genes = gene_by_chrom.get(chrom, [])
        for g_start, g_end, gene_id in chrom_genes:
            g_start_ext = max(0, g_start - EXTENSION)
            g_end_ext = g_end + EXTENSION
            if start < g_end_ext and end > g_start_ext:
                all_har_genes.add(gene_id)
                break

all_har_genes_in_set = all_har_genes & all_genes
print(f"  Genes with any HAR overlap: {len(all_har_genes)} (total), {len(all_har_genes_in_set)} (in set)")

# ============================================
# Enrichment of active HAR genes in GD vs RD
# ============================================
print("\n--- Active HAR Gene Enrichment ---")

# Using LOO-caMPRA (non-circular)
rd_active = len(active_har_genes_in_set & rd_campra)
gd_active = len(active_har_genes_in_set & gd_campra)
print(f"  LOO-caMPRA: RD genes with active HAR: {rd_active}, GD: {gd_active}")

odds_rd_campra, p_rd_campra, _, _ = fisher_test(rd_campra, all_genes, active_har_genes_in_set)
odds_gd_campra, p_gd_campra, _, _ = fisher_test(gd_campra, all_genes, active_har_genes_in_set)
print(f"    RD (LOO-caMPRA): OR={odds_rd_campra:.2f}, p={p_rd_campra:.2e}")
print(f"    GD (LOO-caMPRA): OR={odds_gd_campra:.2f}, p={p_gd_campra:.2e}")

# Using v7 (safe for HAR, not in RDS)
rd_v7_active = len(active_har_genes_in_set & rd_v7)
gd_v7_active = len(active_har_genes_in_set & gd_v7)
print(f"  v7: RD genes with active HAR: {rd_v7_active}, GD: {gd_v7_active}")

odds_rd_v7, p_rd_v7, _, _ = fisher_test(rd_v7, all_genes, active_har_genes_in_set)
odds_gd_v7, p_gd_v7, _, _ = fisher_test(gd_v7, all_genes, active_har_genes_in_set)
print(f"    RD (v7): OR={odds_rd_v7:.2f}, p={p_rd_v7:.2e}")
print(f"    GD (v7): OR={odds_gd_v7:.2f}, p={p_gd_v7:.2e}")

# All HARs (not just active)
odds_rd_all, p_rd_all, _, _ = fisher_test(rd_v7, all_genes, all_har_genes_in_set)
odds_gd_all, p_gd_all, _, _ = fisher_test(gd_v7, all_genes, all_har_genes_in_set)
print(f"\n  All HARs (3,171):")
print(f"    RD (v7): OR={odds_rd_all:.2f}, p={p_rd_all:.2e}")
print(f"    GD (v7): OR={odds_gd_all:.2f}, p={p_gd_all:.2e}")

# ============================================
# Tissue enrichment of active HAR genes
# ============================================
print("\n--- Tissue Enrichment of Active HAR Genes ---")

gtex_tissue_cols = [h for h in fieldnames if h.startswith("gtex_") and h != "gtex_v11_tau"]
BRAIN_TISSUES = ["Brain", "Nerve"]

campra_tissue_results = []

for tissue_col in gtex_tissue_cols:
    tissue_name = tissue_col.replace("gtex_", "")
    
    tpm_values = {}
    for gid, g in genes.items():
        try:
            tpm_values[gid] = float(g.get(tissue_col) or 0)
        except:
            tpm_values[gid] = 0
    
    vals = np.array(list(tpm_values.values()))
    threshold = np.percentile(vals[vals > 0], 75) if np.any(vals > 0) else 0
    high_expr = set(g for g, v in tpm_values.items() if v >= max(threshold, 1.0))
    
    # Non-circular: LOO-caMPRA for RD, also v7 for comparison
    odds_rd_loo, p_rd_loo, n_rd, _ = fisher_test(rd_campra, all_genes, high_expr & active_har_genes_in_set)
    odds_rd_v7, p_rd_v7_t, _, _ = fisher_test(rd_v7, all_genes, high_expr & active_har_genes_in_set)
    odds_gd_v7, p_gd_v7_t, _, _ = fisher_test(gd_v7, all_genes, high_expr & active_har_genes_in_set)
    
    campra_tissue_results.append({
        "tissue": tissue_name,
        "is_brain": tissue_name in BRAIN_TISSUES,
        "n_high_expr_har": len(high_expr & active_har_genes_in_set),
        "rd_loo_campra_OR": odds_rd_loo,
        "rd_loo_campra_p": p_rd_loo,
        "rd_v7_OR": odds_rd_v7,
        "rd_v7_p": p_rd_v7_t,
        "gd_v7_OR": odds_gd_v7,
        "gd_v7_p": p_gd_v7_t,
    })

# BH FDR
if campra_tissue_results:
    for key in ["rd_loo_campra_p", "rd_v7_p", "gd_v7_p"]:
        pvals = [r[key] for r in campra_tissue_results]
        fdrs = bh_fdr(pvals)
        fdr_key = key.replace("_p", "_fdr")
        for i, r in enumerate(campra_tissue_results):
            r[fdr_key] = fdrs[i]

print(f"\n  Active HAR gene tissue enrichment:")
print(f"  {'Tissue':25s} {'RD(LOO-cMPRA) OR':18s} {'FDR':10s} {'RD(v7) OR':10s} {'GD(v7) OR':10s}")
for r in sorted(campra_tissue_results, key=lambda x: x["rd_loo_campra_p"]):
    sig = "*" if r["rd_loo_campra_fdr"] < 0.05 else " "
    print(f"  {r['tissue']:25s} {r['rd_loo_campra_OR']:18.2f}{sig} {r['rd_loo_campra_fdr']:10.2e} "
          f"{r['rd_v7_OR']:10.2f} {r['gd_v7_OR']:10.2f}")

# ============================================
# Brain vs non-brain comparison
# ============================================
print("\n--- Brain vs Non-Brain Comparison ---")

brain_results = [r for r in campra_tissue_results if r["is_brain"]]
non_brain_results = [r for r in campra_tissue_results if not r["is_brain"]]

brain_ors = [r["rd_loo_campra_OR"] for r in brain_results if not np.isnan(r["rd_loo_campra_OR"])]
non_brain_ors = [r["rd_loo_campra_OR"] for r in non_brain_results if not np.isnan(r["rd_loo_campra_OR"])]

print(f"  Brain: mean OR={np.mean(brain_ors):.3f}")
print(f"  Non-brain: mean OR={np.mean(non_brain_ors):.3f}")

if len(brain_ors) > 1 and len(non_brain_ors) > 1:
    u, p = stats.mannwhitneyu(brain_ors, non_brain_ors, alternative="greater")
    print(f"  Brain > non-brain: U={u:.0f}, p={p:.4e}")

# ============================================
# Write outputs
# ============================================
print("\n--- Writing outputs ---")

campra_csv = OUTPUT_DIR / "regulatory_campra_tissue_enrichment.csv"
with open(campra_csv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(campra_tissue_results[0].keys()))
    writer.writeheader()
    for r in campra_tissue_results:
        writer.writerow(r)
print(f"  {campra_csv}")

# HAR gene mapping
har_gene_csv = OUTPUT_DIR / "har_gene_mapping.csv"
with open(har_gene_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["gene_id", "gene_symbol", "classification_v7", "classification_loo_campra", "has_active_har"])
    for gene_id in sorted(all_har_genes_in_set):
        sym = genes[gene_id].get("gene_symbol", "")
        cls_v7 = genes[gene_id]["classification_v7"]
        cls_campra = genes[gene_id].get("classification_loo_campra", "")
        has_active = "Y" if gene_id in active_har_genes_in_set else "N"
        writer.writerow([gene_id, sym, cls_v7, cls_campra, has_active])
print(f"  {har_gene_csv}")

# Summary
summary = {
    "n_active_hars": len(active_hars),
    "n_all_hars": 3171,
    "n_genes_active_har": len(active_har_genes_in_set),
    "n_genes_any_har": len(all_har_genes_in_set),
    "rd_loo_campra": {"OR": float(odds_rd_campra), "p": float(p_rd_campra)},
    "gd_loo_campra": {"OR": float(odds_gd_campra), "p": float(p_gd_campra)},
    "rd_v7_active": {"OR": float(odds_rd_v7), "p": float(p_rd_v7)},
    "gd_v7_active": {"OR": float(odds_gd_v7), "p": float(p_gd_v7)},
    "rd_v7_all_hars": {"OR": float(odds_rd_all), "p": float(p_rd_all)},
    "gd_v7_all_hars": {"OR": float(odds_gd_all), "p": float(p_gd_all)},
    "brain_mean_OR": float(np.mean(brain_ors)) if brain_ors else None,
    "non_brain_mean_OR": float(np.mean(non_brain_ors)) if non_brain_ors else None,
    "fdr_sig_tissues_rd_loo": [r["tissue"] for r in campra_tissue_results if r["rd_loo_campra_fdr"] < 0.05],
}

summary_file = OUTPUT_DIR / "campra_regulatory_summary.json"
with open(summary_file, "w") as f:
    json.dump(summary, f, indent=2)
print(f"  {summary_file}")

print("\n" + "=" * 70)
print("caMPRA Regulatory Tissue Mapping Complete!")
print("=" * 70)
