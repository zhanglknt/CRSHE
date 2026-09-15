#!/usr/bin/env python3
"""
Phase 8 Component 4: Brain Region Resolution

BrainSpan: 26 brain regions × 31 developmental stages
HPA single-nuclei: 34 brain cell types
GTEx V11: 13 brain sub-regions

Uses classification_loo_brain for non-circular brain analysis.

Run in WSL: python3 phase8_brain_resolution.py
"""
import csv
import json
import numpy as np
from scipy import stats
from collections import Counter, defaultdict
from pathlib import Path

BASE = Path("/mnt/d/hs_gene_project")
MASTER_CSV = BASE / "results/phase8_tissue_analysis/data/phase8_master_table.csv"
BRAINSPAN_EXPR = BASE / "data/brainspan/expression_matrix.csv"
BRAINSPAN_COLS = BASE / "data/brainspan/columns_metadata.csv"
BRAINSPAN_ROWS = BASE / "data/brainspan/rows_metadata.csv"
HPA_SN_TSV = BASE / "data/hpa_single_nuclei_brain/rna_single_nuclei_cluster_type.tsv"
OUTPUT_DIR = BASE / "results/phase8_tissue_analysis"

print("=" * 70)
print("Phase 8 Component 4: Brain Region Resolution")
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

# Define gene groups using LOO-brain
gd_genes = set(g for g, v in genes.items() if v.get("classification_loo_brain") == "gene-driven")
rd_genes = set(g for g, v in genes.items() if v.get("classification_loo_brain") == "regulation-driven")
neutral_genes = set(g for g, v in genes.items() if v.get("classification_loo_brain") == "neutral")
all_genes = set(genes.keys())

print(f"  GD (LOO-brain): {len(gd_genes)}")
print(f"  RD (LOO-brain): {len(rd_genes)}")
print(f"  Neutral (LOO-brain): {len(neutral_genes)}")

# ============================================
# Part 1: BrainSpan Analysis
# ============================================
print(f"\n{'='*50}")
print("Part 1: BrainSpan Developmental Expression Analysis")
print(f"{'='*50}")

# Load column metadata (samples)
print("\n  Loading BrainSpan metadata...")
columns = []
with open(BRAINSPAN_COLS) as f:
    reader = csv.DictReader(f)
    for row in reader:
        columns.append(row)
print(f"    {len(columns)} samples")

# Load row metadata (gene mapping)
print("  Loading BrainSpan gene mapping...")
row_map = {}  # row_idx -> ensembl_id
with open(BRAINSPAN_ROWS) as f:
    reader = csv.DictReader(f)
    for i, row in enumerate(reader):
        # Try to find Ensembl gene ID
        eid = row.get("ensembl_gene_id", row.get("gene_id", ""))
        if eid:
            # Remove version number if present
            eid = eid.split(".")[0]
            row_map[i] = eid

print(f"    {len(row_map)} genes mapped")

# Filter to our gene set
row_map_filtered = {k: v for k, v in row_map.items() if v in genes}
print(f"    {len(row_map_filtered)} genes overlap with our 4,974 set")

# Load expression matrix (RPKM)
print("  Loading BrainSpan expression matrix (183MB)...")
# Read header to get sample count
with open(BRAINSPAN_EXPR) as f:
    header = f.readline().strip().split(",")
    n_samples = len(header) - 1
    print(f"    {n_samples} samples in expression matrix")

# Read expression data for our genes only
gene_expr = {}  # ensembl_id -> [rpkm values across samples]
with open(BRAINSPAN_EXPR) as f:
    reader = csv.reader(f)
    header = next(reader)
    for i, row in enumerate(reader):
        if i in row_map_filtered:
            eid = row_map_filtered[i]
            try:
                expr_vals = [float(v) if v else 0.0 for v in row[1:]]
                gene_expr[eid] = expr_vals
            except ValueError:
                pass

print(f"    Expression data loaded for {len(gene_expr)} genes")

# Parse BrainSpan metadata
# Extract: brain region, developmental stage, age
brain_regions = set()
dev_stages = set()
sample_info = []

for i, col in enumerate(columns):
    region = col.get("structure_acronym", col.get("structure_name", ""))
    age = col.get("age", "")
    stage = col.get("developmental_stage", "")
    
    brain_regions.add(region)
    dev_stages.add(age)
    sample_info.append({"region": region, "age": age, "stage": stage})

print(f"    Brain regions: {len(brain_regions)}")
print(f"    Developmental stages: {len(dev_stages)}")

# Classify developmental periods
def get_dev_period(age_str):
    """Classify age string into developmental period."""
    if not age_str:
        return "unknown"
    age_str = str(age_str).lower()
    if "pcw" in age_str:
        pcw = int(age_str.replace("pcw", "").strip())
        if pcw <= 17:
            return "early_prenatal"
        else:
            return "late_prenatal"
    elif "mos" in age_str or "months" in age_str:
        return "early_postnatal"
    elif "yrs" in age_str or "years" in age_str:
        yrs = int(age_str.replace("yrs", "").replace("years", "").strip())
        if yrs <= 5:
            return "early_postnatal"
        else:
            return "late_postnatal"
    else:
        return "unknown"

# ============================================
# BrainSpan Region Enrichment
# ============================================
print(f"\n  --- BrainSpan Region Enrichment ---")

region_results = []

for region in sorted(brain_regions):
    if not region:
        continue
    
    # Get samples for this region
    region_sample_indices = [i for i, s in enumerate(sample_info) if s["region"] == region]
    
    if len(region_sample_indices) < 5:
        continue
    
    # Compute mean RPKM for each gene in this region
    gene_rpkm = {}
    for eid, expr_vals in gene_expr.items():
        region_vals = [expr_vals[i] for i in region_sample_indices if i < len(expr_vals)]
        if region_vals:
            gene_rpkm[eid] = np.mean(region_vals)
    
    if len(gene_rpkm) < 100:
        continue
    
    # Define high-expression genes (top 25%)
    rpkm_vals = list(gene_rpkm.values())
    threshold = np.percentile(rpkm_vals, 75)
    high_expr = set(g for g, v in gene_rpkm.items() if v >= threshold)
    
    # Fisher tests
    gd_odds, gd_p, gd_n, _ = fisher_test(gd_genes, all_genes, high_expr)
    rd_odds, rd_p, rd_n, _ = fisher_test(rd_genes, all_genes, high_expr)
    
    region_results.append({
        "region": region,
        "n_samples": len(region_sample_indices),
        "n_genes_expressed": len(gene_rpkm),
        "n_high_expr": len(high_expr),
        "gd_OR": gd_odds,
        "gd_p": gd_p,
        "gd_n_high": gd_n,
        "rd_OR": rd_odds,
        "rd_p": rd_p,
        "rd_n_high": rd_n,
    })

# BH FDR
if region_results:
    gd_pvals = [r["gd_p"] for r in region_results]
    rd_pvals = [r["rd_p"] for r in region_results]
    gd_fdrs = bh_fdr(gd_pvals)
    rd_fdrs = bh_fdr(rd_pvals)
    for i, r in enumerate(region_results):
        r["gd_fdr"] = gd_fdrs[i]
        r["rd_fdr"] = rd_fdrs[i]

    print(f"\n  BrainSpan regions ({len(region_results)}):")
    print(f"  {'Region':25s} {'N':5s} {'GD OR':8s} {'GD FDR':10s} {'RD OR':8s} {'RD FDR':10s}")
    for r in sorted(region_results, key=lambda x: x["gd_p"]):
        gd_s = "*" if r["gd_fdr"] < 0.05 else " "
        rd_s = "*" if r["rd_fdr"] < 0.05 else " "
        print(f"  {r['region']:25s} {r['n_samples']:5d} "
              f"{r['gd_OR']:8.2f}{gd_s} {r['gd_fdr']:10.2e} "
              f"{r['rd_OR']:8.2f}{rd_s} {r['rd_fdr']:10.2e}")

# Write BrainSpan region results
bs_region_csv = OUTPUT_DIR / "brainspan_region_enrichment.csv"
with open(bs_region_csv, "w", newline="") as f:
    if region_results:
        writer = csv.DictWriter(f, fieldnames=list(region_results[0].keys()))
        writer.writeheader()
        for r in region_results:
            writer.writerow(r)
print(f"\n  Output: {bs_region_csv}")

# ============================================
# Developmental Period Enrichment
# ============================================
print(f"\n  --- Developmental Period Enrichment ---")

# Group samples by developmental period
period_samples = defaultdict(list)
for i, s in enumerate(sample_info):
    period = get_dev_period(s["age"])
    period_samples[period].append(i)

print(f"  Developmental periods:")
for period, indices in sorted(period_samples.items()):
    print(f"    {period}: {len(indices)} samples")

period_results = []

for period, sample_indices in period_samples.items():
    if len(sample_indices) < 5 or period == "unknown":
        continue
    
    # Compute mean RPKM per gene in this period
    gene_rpkm = {}
    for eid, expr_vals in gene_expr.items():
        period_vals = [expr_vals[i] for i in sample_indices if i < len(expr_vals)]
        if period_vals:
            gene_rpkm[eid] = np.mean(period_vals)
    
    if len(gene_rpkm) < 100:
        continue
    
    rpkm_vals = list(gene_rpkm.values())
    threshold = np.percentile(rpkm_vals, 75)
    high_expr = set(g for g, v in gene_rpkm.items() if v >= threshold)
    
    gd_odds, gd_p, gd_n, _ = fisher_test(gd_genes, all_genes, high_expr)
    rd_odds, rd_p, rd_n, _ = fisher_test(rd_genes, all_genes, high_expr)
    
    period_results.append({
        "period": period,
        "n_samples": len(sample_indices),
        "n_genes": len(gene_rpkm),
        "n_high_expr": len(high_expr),
        "gd_OR": gd_odds,
        "gd_p": gd_p,
        "gd_n_high": gd_n,
        "rd_OR": rd_odds,
        "rd_p": rd_p,
        "rd_n_high": rd_n,
    })

print(f"\n  Developmental period results:")
for r in period_results:
    print(f"    {r['period']:20s} N={r['n_samples']:4d} "
          f"GD:OR={r['gd_OR']:.2f},p={r['gd_p']:.2e} "
          f"RD:OR={r['rd_OR']:.2f},p={r['rd_p']:.2e}")

# Prenatal vs Postnatal comparison
print(f"\n  --- Prenatal vs Postnatal ---")

# Compute per-gene prenatal and postnatal mean
prenatal_indices = period_samples.get("early_prenatal", []) + period_samples.get("late_prenatal", [])
postnatal_indices = period_samples.get("early_postnatal", []) + period_samples.get("late_postnatal", [])

gene_prenatal = {}
gene_postnatal = {}

for eid, expr_vals in gene_expr.items():
    if prenatal_indices:
        vals = [expr_vals[i] for i in prenatal_indices if i < len(expr_vals)]
        gene_prenatal[eid] = np.mean(vals) if vals else 0
    if postnatal_indices:
        vals = [expr_vals[i] for i in postnatal_indices if i < len(expr_vals)]
        gene_postnatal[eid] = np.mean(vals) if vals else 0

# Compute prenatal/postnatal ratio
gene_ratio = {}
for eid in gene_prenatal:
    post = gene_postnatal.get(eid, 0)
    if post > 0.1:
        gene_ratio[eid] = gene_prenatal[eid] / post
    elif gene_prenatal[eid] > 0.1:
        gene_ratio[eid] = float('inf')
    else:
        gene_ratio[eid] = 1.0

# Compare ratio between GD and RD
gd_ratios = [gene_ratio[g] for g in gd_genes if g in gene_ratio and np.isfinite(gene_ratio[g])]
rd_ratios = [gene_ratio[g] for g in rd_genes if g in gene_ratio and np.isfinite(gene_ratio[g])]

if gd_ratios and rd_ratios:
    u_pre, p_pre = stats.mannwhitneyu(rd_ratios, gd_ratios, alternative="greater")
    print(f"  Prenatal/Postnatal ratio RD > GD: U={u_pre:.0f}, p={p_pre:.2e}")
    print(f"    GD median ratio: {np.median(gd_ratios):.4f} (n={len(gd_ratios)})")
    print(f"    RD median ratio: {np.median(rd_ratios):.4f} (n={len(rd_ratios)})")

    # Also GD > RD
    u_gd, p_gd = stats.mannwhitneyu(gd_ratios, rd_ratios, alternative="greater")
    print(f"  Prenatal/Postnatal ratio GD > RD: U={u_gd:.0f}, p={p_gd:.2e}")

# Write developmental results
dev_csv = OUTPUT_DIR / "brainspan_developmental_stages.csv"
with open(dev_csv, "w", newline="") as f:
    if period_results:
        writer = csv.DictWriter(f, fieldnames=list(period_results[0].keys()))
        writer.writeheader()
        for r in period_results:
            writer.writerow(r)
print(f"\n  Output: {dev_csv}")

# ============================================
# Part 2: HPA Single-Nuclei Brain Cell Types
# ============================================
print(f"\n{'='*50}")
print("Part 2: HPA Single-Nuclei Brain Cell Types")
print(f"{'='*50}")

# Load HPA single-nuclei data
print("\n  Loading HPA single-nuclei data...")
hpa_sn = []
with open(HPA_SN_TSV) as f:
    reader = csv.DictReader(f, delimiter="\t")
    hpa_fields = reader.fieldnames
    for row in reader:
        hpa_sn.append(row)
print(f"    {len(hpa_sn)} rows, {len(hpa_fields)} columns")
print(f"    Columns: {hpa_fields[:5]}...")

# Identify gene ID column and cell type columns
gene_id_col = None
for candidate in ["Gene", "gene_id", "Gene name", "Ensembl"]:
    if candidate in hpa_fields:
        gene_id_col = candidate
        break

if not gene_id_col:
    gene_id_col = hpa_fields[0]
print(f"    Gene ID column: {gene_id_col}")

# Cell type columns (all numeric columns)
cell_type_cols = [c for c in hpa_fields if c != gene_id_col and "Brain" not in c and "Prediction" not in c]
print(f"    Cell types: {len(cell_type_cols)}")
if cell_type_cols:
    print(f"    First 5: {cell_type_cols[:5]}")

# Map HPA genes to our gene set
# HPA uses Ensembl IDs or gene symbols
hpa_gene_map = {}
unmapped = 0
for row in hpa_sn:
    gid = row.get(gene_id_col, "")
    if gid:
        # Try direct match (Ensembl ID)
        eid = gid.split(".")[0] if "." in gid else gid
        if eid in genes:
            hpa_gene_map[eid] = row
        elif gid.startswith("ENSG"):
            # Try without version
            if eid in genes:
                hpa_gene_map[eid] = row
            else:
                unmapped += 1
        else:
            # Try gene symbol match
            for g, v in genes.items():
                if v.get("gene_symbol") == gid:
                    hpa_gene_map[g] = row
                    break
            else:
                unmapped += 1

print(f"    Mapped: {len(hpa_gene_map)}, Unmapped: {unmapped}")

# Cell type enrichment analysis
celltype_results = []

for ct in cell_type_cols:
    # Get expression values for each gene
    gene_expr_ct = {}
    for gid, row in hpa_gene_map.items():
        try:
            val = float(row.get(ct, 0) or 0)
            gene_expr_ct[gid] = val
        except (ValueError, TypeError):
            gene_expr_ct[gid] = 0
    
    if len(gene_expr_ct) < 100:
        continue
    
    # Define high-expression (top 25%)
    vals = list(gene_expr_ct.values())
    threshold = np.percentile(vals, 75)
    high_expr = set(g for g, v in gene_expr_ct.items() if v >= threshold and v > 0)
    
    if len(high_expr) < 10:
        continue
    
    # Fisher tests
    gd_odds, gd_p, gd_n, _ = fisher_test(gd_genes & set(gene_expr_ct.keys()), 
                                          set(gene_expr_ct.keys()), high_expr)
    rd_odds, rd_p, rd_n, _ = fisher_test(rd_genes & set(gene_expr_ct.keys()),
                                          set(gene_expr_ct.keys()), high_expr)
    
    celltype_results.append({
        "cell_type": ct,
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
    print(f"  {'Cell Type':35s} {'GD OR':8s} {'GD FDR':10s} {'RD OR':8s} {'RD FDR':10s}")
    for r in sorted(celltype_results, key=lambda x: x["gd_p"])[:10]:
        gd_s = "*" if r["gd_fdr"] < 0.05 else " "
        rd_s = "*" if r["rd_fdr"] < 0.05 else " "
        print(f"  {r['cell_type']:35s} {r['gd_OR']:8.2f}{gd_s} {r['gd_fdr']:10.2e} "
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
# Neuronal vs Glial comparison
# ============================================
print(f"\n  --- Neuronal vs Glial Comparison ---")

neuronal_types = [ct for ct in cell_type_cols if any(k in ct.lower() for k in 
    ["excitatory", "inhibitory", "interneuron", "ca1", "ca4", "dentate", "spiny", "neuron", "corticothalamic", "intratelencephalic", "near-projecting", "mammillary", "midbrain-derived"])]
glial_types = [ct for ct in cell_type_cols if any(k in ct.lower() for k in 
    ["astrocyte", "oligodendrocyte", "opc", "bergmann", "glia", "macrophage", "microglia"])]

print(f"    Neuronal types: {len(neuronal_types)}")
print(f"    Glial types: {len(glial_types)}")

for group_name, group_types in [("Neuronal", neuronal_types), ("Glial", glial_types)]:
    if not group_types:
        continue
    
    # Compute mean expression across cell types in this group
    gene_expr_group = {}
    for gid, row in hpa_gene_map.items():
        vals = []
        for ct in group_types:
            try:
                vals.append(float(row.get(ct, 0) or 0))
            except:
                pass
        if vals:
            gene_expr_group[gid] = np.mean(vals)
    
    if len(gene_expr_group) < 100:
        continue
    
    vals = list(gene_expr_group.values())
    threshold = np.percentile(vals, 75)
    high_expr = set(g for g, v in gene_expr_group.items() if v >= threshold and v > 0)
    
    gd_odds, gd_p, gd_n, _ = fisher_test(gd_genes & set(gene_expr_group.keys()),
                                          set(gene_expr_group.keys()), high_expr)
    rd_odds, rd_p, rd_n, _ = fisher_test(rd_genes & set(gene_expr_group.keys()),
                                          set(gene_expr_group.keys()), high_expr)
    
    print(f"    {group_name}: GD OR={gd_odds:.2f} (p={gd_p:.2e}), RD OR={rd_odds:.2f} (p={rd_p:.2e})")

# ============================================
# Summary
# ============================================
summary = {
    "brainspan": {
        "n_samples": len(columns),
        "n_regions": len(brain_regions),
        "n_genes_matched": len(gene_expr),
        "gd_fdr_sig_regions": [r["region"] for r in region_results if r.get("gd_fdr", 1) < 0.05],
        "rd_fdr_sig_regions": [r["region"] for r in region_results if r.get("rd_fdr", 1) < 0.05],
        "dev_periods": {k: len(v) for k, v in period_samples.items()},
    },
    "hpa_single_nuclei": {
        "n_cell_types": len(cell_type_cols),
        "n_genes_mapped": len(hpa_gene_map),
        "gd_fdr_sig_types": [r["cell_type"] for r in celltype_results if r.get("gd_fdr", 1) < 0.05],
        "rd_fdr_sig_types": [r["cell_type"] for r in celltype_results if r.get("rd_fdr", 1) < 0.05],
    },
}

summary_file = OUTPUT_DIR / "brain_resolution_summary.json"
with open(summary_file, "w") as f:
    json.dump(summary, f, indent=2)

print(f"\n  Summary: {summary_file}")
print("\n" + "=" * 70)
print("Brain Region Resolution Analysis Complete!")
print("=" * 70)
