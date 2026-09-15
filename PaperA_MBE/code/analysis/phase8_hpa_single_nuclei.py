#!/usr/bin/env python3
"""
Phase 8 Component 4b: HPA Single-Nuclei Brain Cell Type Analysis

HPA data is in LONG format (Gene, Gene name, Cluster type, nCPM).
This script pivots to wide format and performs cell type enrichment.

Uses classification_loo_brain for non-circular brain analysis.

Run in WSL: python3 phase8_hpa_single_nuclei.py
"""
import csv
import json
import numpy as np
from scipy import stats
from collections import defaultdict
from pathlib import Path

BASE = Path("/mnt/d/hs_gene_project")
MASTER_CSV = BASE / "results/phase8_tissue_analysis/data/phase8_master_table.csv"
HPA_SN_TSV = BASE / "data/hpa_single_nuclei_brain/rna_single_nuclei_cluster_type.tsv"
HPA_CLUSTER_TYPES = BASE / "data/hpa_single_nuclei_brain/rna_single_nuclei_cluster_type_cluster_types.tsv"
OUTPUT_DIR = BASE / "results/phase8_tissue_analysis"

print("=" * 70)
print("Phase 8 Component 4b: HPA Single-Nuclei Brain Cell Type Analysis")
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
    for row in reader:
        genes[row["gene_id"]] = row

# Use LOO-brain classification for non-circular brain analysis
gd_genes = set(g for g, v in genes.items() if v.get("classification_loo_brain") == "gene-driven")
rd_genes = set(g for g, v in genes.items() if v.get("classification_loo_brain") == "regulation-driven")
neutral_genes = set(g for g, v in genes.items() if v.get("classification_loo_brain") == "neutral")
all_genes = set(genes.keys())

print(f"  GD (LOO-brain): {len(gd_genes)}")
print(f"  RD (LOO-brain): {len(rd_genes)}")
print(f"  Neutral (LOO-brain): {len(neutral_genes)}")

# ============================================
# Load cluster type metadata (neuronal/glial classification)
# ============================================
print("\n--- Loading cluster type metadata ---")
cluster_class = {}  # cluster_type -> cell_type_class
with open(HPA_CLUSTER_TYPES) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        ct = row.get("Cluster type", "")
        cc = row.get("Cell type class", "")
        if ct and cc:
            cluster_class[ct] = cc

print(f"  {len(cluster_class)} cluster types")
cell_classes = set(cluster_class.values())
print(f"  Cell type classes: {cell_classes}")

# ============================================
# Read HPA long-format TSV and pivot to wide
# ============================================
print("\n--- Reading HPA single-nuclei data (long format) ---")
print("  This may take a moment (665K rows)...")

# Accumulate: gene_id -> {cluster_type: [nCPM values]}
gene_ct_values = defaultdict(lambda: defaultdict(list))

with open(HPA_SN_TSV) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        gid = row.get("Gene", "")
        ct = row.get("Cluster type", "")
        ncpm_str = row.get("nCPM", "0")
        
        if not gid or not ct:
            continue
        
        # Remove Ensembl version if present
        eid = gid.split(".")[0] if "." in gid else gid
        
        try:
            ncpm = float(ncpm_str) if ncpm_str else 0.0
        except ValueError:
            ncpm = 0.0
        
        gene_ct_values[eid][ct].append(ncpm)

print(f"  Total genes in HPA: {len(gene_ct_values)}")

# Pivot to wide format: gene_id -> {cluster_type: mean nCPM}
all_cluster_types = sorted(set(ct for cts in gene_ct_values.values() for ct in cts.keys()))
print(f"  Cluster types found: {len(all_cluster_types)}")

gene_wide = {}  # gene_id -> {ct: mean_nCPM}
for eid, ct_vals in gene_ct_values.items():
    gene_wide[eid] = {ct: float(np.mean(vals)) for ct, vals in ct_vals.items()}

# Filter to our gene set
gene_wide_filtered = {g: v for g, v in gene_wide.items() if g in genes}
print(f"  Genes matched to our 4,974 set: {len(gene_wide_filtered)}")

# ============================================
# Cell Type Enrichment Analysis
# ============================================
print(f"\n--- Cell Type Enrichment Analysis ---")

celltype_results = []

for ct in all_cluster_types:
    # Get expression values for each gene
    gene_expr_ct = {}
    for gid in gene_wide_filtered:
        val = gene_wide_filtered[gid].get(ct, 0.0)
        gene_expr_ct[gid] = val
    
    if len(gene_expr_ct) < 100:
        continue
    
    # Define high-expression (top 25% among non-zero)
    vals = np.array(list(gene_expr_ct.values()))
    non_zero = vals[vals > 0]
    if len(non_zero) < 50:
        continue
    threshold = np.percentile(non_zero, 75)
    high_expr = set(g for g, v in gene_expr_ct.items() if v >= threshold and v > 0)
    
    if len(high_expr) < 10:
        continue
    
    # Fisher tests
    bg = set(gene_expr_ct.keys())
    gd_odds, gd_p, gd_n, _ = fisher_test(gd_genes & bg, bg, high_expr)
    rd_odds, rd_p, rd_n, _ = fisher_test(rd_genes & bg, bg, high_expr)
    
    cell_class = cluster_class.get(ct, "Unknown")
    
    celltype_results.append({
        "cell_type": ct,
        "cell_class": cell_class,
        "n_genes": len(gene_expr_ct),
        "n_high_expr": len(high_expr),
        "gd_OR": gd_odds,
        "gd_p": gd_p,
        "gd_n_high": gd_n,
        "rd_OR": rd_odds,
        "rd_p": rd_p,
        "rd_n_high": rd_n,
    })

# BH FDR
if celltype_results:
    gd_pvals = [r["gd_p"] for r in celltype_results]
    rd_pvals = [r["rd_p"] for r in celltype_results]
    gd_fdrs = bh_fdr(gd_pvals)
    rd_fdrs = bh_fdr(rd_pvals)
    for i, r in enumerate(celltype_results):
        r["gd_fdr"] = gd_fdrs[i]
        r["rd_fdr"] = rd_fdrs[i]

    print(f"\n  HPA Single-Nuclei Cell Types ({len(celltype_results)}):")
    print(f"  {'Cell Type':40s} {'Class':20s} {'GD OR':8s} {'GD FDR':10s} {'RD OR':8s} {'RD FDR':10s}")
    for r in sorted(celltype_results, key=lambda x: x["gd_p"]):
        gd_s = "*" if r["gd_fdr"] < 0.05 else " "
        rd_s = "*" if r["rd_fdr"] < 0.05 else " "
        print(f"  {r['cell_type']:40s} {r['cell_class']:20s} "
              f"{r['gd_OR']:8.2f}{gd_s} {r['gd_fdr']:10.2e} "
              f"{r['rd_OR']:8.2f}{rd_s} {r['rd_fdr']:10.2e}")

# Write HPA results
hpa_csv = OUTPUT_DIR / "hpa_single_nuclei_enrichment.csv"
with open(hpa_csv, "w", newline="") as f:
    if celltype_results:
        writer = csv.DictWriter(f, fieldnames=list(celltype_results[0].keys()))
        writer.writeheader()
        for r in celltype_results:
            writer.writerow(r)
print(f"\n  Output: {hpa_csv}")

# ============================================
# Neuronal vs Glial Comparison
# ============================================
print(f"\n--- Neuronal vs Glial Comparison ---")

# Group cell types by class
class_types = defaultdict(list)
for ct, cc in cluster_class.items():
    if ct in all_cluster_types:
        class_types[cc].append(ct)

print(f"  Cell type classes:")
for cc, cts in sorted(class_types.items()):
    print(f"    {cc}: {len(cts)} types")

# Compute per-gene mean expression for each cell class
class_results = []

for cc, ct_list in sorted(class_types.items()):
    if not ct_list:
        continue
    
    gene_expr_class = {}
    for gid, ct_vals in gene_wide_filtered.items():
        vals = [ct_vals.get(ct, 0.0) for ct in ct_list]
        gene_expr_class[gid] = float(np.mean(vals)) if vals else 0.0
    
    if len(gene_expr_class) < 100:
        continue
    
    vals = np.array(list(gene_expr_class.values()))
    non_zero = vals[vals > 0]
    if len(non_zero) < 50:
        continue
    threshold = np.percentile(non_zero, 75)
    high_expr = set(g for g, v in gene_expr_class.items() if v >= threshold and v > 0)
    
    bg = set(gene_expr_class.keys())
    gd_odds, gd_p, gd_n, _ = fisher_test(gd_genes & bg, bg, high_expr)
    rd_odds, rd_p, rd_n, _ = fisher_test(rd_genes & bg, bg, high_expr)
    
    class_results.append({
        "cell_class": cc,
        "n_types": len(ct_list),
        "n_genes": len(gene_expr_class),
        "n_high_expr": len(high_expr),
        "gd_OR": gd_odds,
        "gd_p": gd_p,
        "gd_n_high": gd_n,
        "rd_OR": rd_odds,
        "rd_p": rd_p,
        "rd_n_high": rd_n,
    })

print(f"\n  Cell class enrichment:")
print(f"  {'Class':30s} {'N types':8s} {'GD OR':8s} {'GD p':10s} {'RD OR':8s} {'RD p':10s}")
for r in class_results:
    print(f"  {r['cell_class']:30s} {r['n_types']:8d} "
          f"{r['gd_OR']:8.2f} {r['gd_p']:10.2e} "
          f"{r['rd_OR']:8.2f} {r['rd_p']:10.2e}")

# Write cell class results
class_csv = OUTPUT_DIR / "hpa_cell_class_enrichment.csv"
with open(class_csv, "w", newline="") as f:
    if class_results:
        writer = csv.DictWriter(f, fieldnames=list(class_results[0].keys()))
        writer.writeheader()
        for r in class_results:
            writer.writerow(r)
print(f"\n  Output: {class_csv}")

# ============================================
# Neuronal-specific vs Glial-specific gene comparison
# ============================================
print(f"\n--- Neuronal-specific vs Glial-specific Gene Expression ---")

neuronal_types = class_types.get("Neuronal cells", [])
glial_types = class_types.get("Glial cells", [])

if neuronal_types and glial_types:
    # Per-gene: neuronal specificity = mean(neuronal) / (mean(neuronal) + mean(glial))
    gene_neuronal_spec = {}
    for gid, ct_vals in gene_wide_filtered.items():
        n_vals = [ct_vals.get(ct, 0.0) for ct in neuronal_types]
        g_vals = [ct_vals.get(ct, 0.0) for ct in glial_types]
        n_mean = float(np.mean(n_vals)) if n_vals else 0
        g_mean = float(np.mean(g_vals)) if g_vals else 0
        total = n_mean + g_mean
        gene_neuronal_spec[gid] = (n_mean / total) if total > 0 else 0.5
    
    gd_spec = [gene_neuronal_spec[g] for g in gd_genes if g in gene_neuronal_spec]
    rd_spec = [gene_neuronal_spec[g] for g in rd_genes if g in gene_neuronal_spec]
    neutral_spec = [gene_neuronal_spec[g] for g in neutral_genes if g in gene_neuronal_spec]
    
    print(f"  Neuronal specificity score (0=glial, 1=neuronal):")
    print(f"    GD: median={np.median(gd_spec):.4f} (n={len(gd_spec)})")
    print(f"    RD: median={np.median(rd_spec):.4f} (n={len(rd_spec)})")
    print(f"    Neutral: median={np.median(neutral_spec):.4f} (n={len(neutral_spec)})")
    
    if gd_spec and rd_spec:
        u, p = stats.mannwhitneyu(rd_spec, gd_spec, alternative="greater")
        print(f"    RD > GD (neuronal): U={u:.0f}, p={p:.2e}")
        u2, p2 = stats.mannwhitneyu(gd_spec, rd_spec, alternative="greater")
        print(f"    GD > RD (neuronal): U={u2:.0f}, p={p2:.2e}")

# ============================================
# Summary
# ============================================
summary = {
    "hpa_single_nuclei": {
        "n_total_genes": len(gene_wide),
        "n_matched_genes": len(gene_wide_filtered),
        "n_cluster_types": len(all_cluster_types),
        "n_cell_classes": len(cell_classes),
        "gd_fdr_sig_types": [r["cell_type"] for r in celltype_results if r.get("gd_fdr", 1) < 0.05],
        "rd_fdr_sig_types": [r["cell_type"] for r in celltype_results if r.get("rd_fdr", 1) < 0.05],
        "cell_class_results": class_results,
    },
}

summary_file = OUTPUT_DIR / "hpa_single_nuclei_summary.json"
with open(summary_file, "w") as f:
    json.dump(summary, f, indent=2)

print(f"\n  Summary: {summary_file}")
print("\n" + "=" * 70)
print("HPA Single-Nuclei Analysis Complete!")
print("=" * 70)
