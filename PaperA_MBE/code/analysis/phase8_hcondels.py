#!/usr/bin/env python3
"""
Phase 8 Component 5a: hCONDELs Processing and Tissue Enrichment

Parse hCONDELs supplementary table (XLS with sub-header row),
map to genes using both coordinate overlap and gene name columns,
test tissue enrichment across GTEx V11 tissues.

Uses classification_v7 (hCONDELs NOT in RDS — non-circular).

Run in WSL: python3 phase8_hcondels.py
"""
import csv
import json
import re
import numpy as np
from scipy import stats
from collections import Counter, defaultdict
from pathlib import Path

BASE = Path("/mnt/d/hs_gene_project")
HCONDELS_XLS = BASE / "data/downloads/hcondels/hCONDELs_supplementary_table2.xls"
MASTER_CSV = BASE / "results/phase8_tissue_analysis/data/phase8_master_table.csv"
GENE_BED = BASE / "results/phase4_conservation/gene_body.bed"
OUTPUT_DIR = BASE / "results/phase8_tissue_analysis"

print("=" * 70)
print("Phase 8 Component 5a: hCONDELs Processing and Tissue Enrichment")
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

def parse_coord(coord_str):
    """Parse coordinate string like 'chr1:3551306-3551400'"""
    if not coord_str or not isinstance(coord_str, str):
        return None
    m = re.match(r'(chr\w+):(\d+)-(\d+)', coord_str)
    if m:
        return {"chrom": m.group(1), "start": int(m.group(2)), "end": int(m.group(3))}
    return None

# ============================================
# Load master table
# ============================================
print("\n--- Loading master table ---")
genes = {}
symbol_to_ensembl = {}
with open(MASTER_CSV) as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames
    for row in reader:
        genes[row["gene_id"]] = row
        sym = row.get("gene_symbol", "")
        if sym:
            symbol_to_ensembl[sym] = row["gene_id"]

all_genes = set(genes.keys())
gd_genes = set(g for g, v in genes.items() if v["classification_v7"] == "gene-driven")
rd_genes = set(g for g, v in genes.items() if v["classification_v7"] == "regulation-driven")
neutral_genes = set(g for g, v in genes.items() if v["classification_v7"] == "neutral")

print(f"  GD: {len(gd_genes)}, RD: {len(rd_genes)}, Neutral: {len(neutral_genes)}")
print(f"  Gene symbol map: {len(symbol_to_ensembl)}")

# ============================================
# Parse hCONDELs XLS (skip sub-header row)
# ============================================
print(f"\n--- Parsing hCONDELs XLS ---")

import pandas as pd
df = pd.read_excel(str(HCONDELS_XLS))

# First row is sub-header, skip it
df = df.iloc[1:].reset_index(drop=True)

# Rename columns for clarity
col_map = {
    "Human Conserved Sequence Deletion": "hcondel_name",
    "Unnamed: 1": "type",
    "Unnamed: 2": "pt2_size",
    "Unnamed: 3": "panTro2_coords",
    "Unnamed: 4": "hg18_size",
    "Unnamed: 5": "hg18_coords",
    "Unnamed: 6": "cons",
    "Upstream gene": "upstream_gene",
    "Unnamed: 8": "upstream_dist",
    "Within": "within_gene",
    "Downstream gene": "downstream_gene",
    "Unnamed: 11": "downstream_dist",
}
df = df.rename(columns=col_map)
print(f"  {len(df)} hCONDELs after removing sub-header")

# Parse coordinates
hcondels = []
for _, row in df.iterrows():
    name = str(row.get("hcondel_name", ""))
    htype = str(row.get("type", ""))
    
    # Try hg18 coordinates first (human genome)
    coord = parse_coord(str(row.get("hg18_coords", "")))
    if not coord:
        coord = parse_coord(str(row.get("panTro2_coords", "")))
    
    within = str(row.get("within_gene", "")).strip()
    if within == "nan" or not within:
        within = ""
    
    upstream = str(row.get("upstream_gene", "")).strip()
    if upstream == "nan" or not upstream:
        upstream = ""
    
    upstream_dist = row.get("upstream_dist", "")
    try:
        upstream_dist = int(float(upstream_dist)) if upstream_dist and str(upstream_dist) != "nan" else None
    except:
        upstream_dist = None
    
    downstream = str(row.get("downstream_gene", "")).strip()
    if downstream == "nan" or not downstream:
        downstream = ""
    
    downstream_dist = row.get("downstream_dist", "")
    try:
        downstream_dist = int(float(downstream_dist)) if downstream_dist and str(downstream_dist) != "nan" else None
    except:
        downstream_dist = None
    
    hcondels.append({
        "name": name,
        "type": htype,
        "coord": coord,
        "within": within,
        "upstream": upstream,
        "upstream_dist": upstream_dist,
        "downstream": downstream,
        "downstream_dist": downstream_dist,
    })

print(f"  Parsed {len(hcondels)} hCONDELs")
n_with_coord = sum(1 for h in hcondels if h["coord"])
n_with_within = sum(1 for h in hcondels if h["within"])
n_with_upstream = sum(1 for h in hcondels if h["upstream"])
n_with_downstream = sum(1 for h in hcondels if h["downstream"])
print(f"  With coordinates: {n_with_coord}")
print(f"  With within gene: {n_with_within}")
print(f"  With upstream gene: {n_with_upstream}")
print(f"  With downstream gene: {n_with_downstream}")

# ============================================
# Map hCONDELs to genes
# Strategy: 1) Within gene (direct overlap)
#           2) Upstream/Downstream within 50kb
#           3) Coordinate overlap with gene_body.bed (±50kb)
# ============================================
print(f"\n--- Mapping hCONDELs to genes ---")

# Method 1: Gene name mapping (Within, Upstream, Downstream)
hcondel_genes_by_name = set()

for hc in hcondels:
    # Within gene = direct overlap
    if hc["within"]:
        eid = symbol_to_ensembl.get(hc["within"])
        if eid:
            hcondel_genes_by_name.add(eid)
    
    # Upstream gene within 50kb
    if hc["upstream"] and hc["upstream_dist"] is not None:
        if hc["upstream_dist"] <= 50000:
            eid = symbol_to_ensembl.get(hc["upstream"])
            if eid:
                hcondel_genes_by_name.add(eid)
    
    # Downstream gene within 50kb
    if hc["downstream"] and hc["downstream_dist"] is not None:
        if hc["downstream_dist"] <= 50000:
            eid = symbol_to_ensembl.get(hc["downstream"])
            if eid:
                hcondel_genes_by_name.add(eid)

print(f"  Method 1 (gene name mapping, ±50kb): {len(hcondel_genes_by_name)} genes")
print(f"    In our gene set: {len(hcondel_genes_by_name & all_genes)}")

# Method 2: Coordinate overlap with gene_body.bed
print(f"\n  Method 2: Coordinate overlap with gene body (±50kb)...")
gene_coords = {}
with open(GENE_BED) as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) >= 4:
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            name = parts[3]
            gene_coords[name] = (chrom, start, end)

print(f"    Gene coordinates loaded: {len(gene_coords)}")

# Group genes by chromosome for efficiency
gene_by_chrom = defaultdict(list)
for gene_id, (chrom, start, end) in gene_coords.items():
    gene_by_chrom[chrom].append((start, end, gene_id))

# Sort by start position
for chrom in gene_by_chrom:
    gene_by_chrom[chrom].sort()

EXTENSION = 50000
hcondel_genes_by_coord = set()

for hc in hcondels:
    if not hc["coord"]:
        continue
    
    hc_chrom = hc["coord"]["chrom"]
    hc_start = hc["coord"]["start"]
    hc_end = hc["coord"]["end"]
    
    # Binary search for overlapping genes
    chrom_genes = gene_by_chrom.get(hc_chrom, [])
    if not chrom_genes:
        continue
    
    # Linear scan (could optimize with binary search, but 4974 genes is fine)
    for g_start, g_end, gene_id in chrom_genes:
        g_start_ext = max(0, g_start - EXTENSION)
        g_end_ext = g_end + EXTENSION
        if hc_start < g_end_ext and hc_end > g_start_ext:
            hcondel_genes_by_coord.add(gene_id)
            break

print(f"    Method 2 (coordinate overlap, ±50kb): {len(hcondel_genes_by_coord)} genes")
print(f"    In our gene set: {len(hcondel_genes_by_coord & all_genes)}")

# Combine both methods
hcondel_genes = hcondel_genes_by_name | hcondel_genes_by_coord
hcondel_genes = hcondel_genes & all_genes

print(f"\n  Combined (union): {len(hcondel_genes)} genes in our set")

# Classification of hCONDEL genes
hcondel_in_gd = hcondel_genes & gd_genes
hcondel_in_rd = hcondel_genes & rd_genes
hcondel_in_neutral = hcondel_genes & neutral_genes
print(f"  GD genes with hCONDEL: {len(hcondel_in_gd)}")
print(f"  RD genes with hCONDEL: {len(hcondel_in_rd)}")
print(f"  Neutral genes with hCONDEL: {len(hcondel_in_neutral)}")

# Fisher test: enrichment in hCONDEL genes
odds_rd, p_rd, _, _ = fisher_test(rd_genes, all_genes, hcondel_genes)
odds_gd, p_gd, _, _ = fisher_test(gd_genes, all_genes, hcondel_genes)
print(f"\n  hCONDEL enrichment:")
print(f"    RD: OR={odds_rd:.2f}, p={p_rd:.2e}")
print(f"    GD: OR={odds_gd:.2f}, p={p_gd:.2e}")

# ============================================
# Tissue enrichment of hCONDEL genes
# ============================================
print(f"\n--- Tissue enrichment of hCONDEL genes ---")

gtex_tissue_cols = [h for h in fieldnames if h.startswith("gtex_") and h != "gtex_v11_tau"]
BRAIN_TISSUES = ["Brain", "Nerve"]

hcondel_tissue_results = []

for tissue_col in gtex_tissue_cols:
    tissue_name = tissue_col.replace("gtex_", "")
    
    # Get TPM values
    tpm_values = {}
    for gid, g in genes.items():
        try:
            tpm_values[gid] = float(g.get(tissue_col) or 0)
        except:
            tpm_values[gid] = 0
    
    # Define high-expression genes
    vals = np.array(list(tpm_values.values()))
    threshold = np.percentile(vals[vals > 0], 75) if np.any(vals > 0) else 0
    high_expr = set(g for g, v in tpm_values.items() if v >= max(threshold, 1.0))
    
    # Fisher: hCONDEL genes enriched in high-expression of this tissue?
    odds, p, n_high, _ = fisher_test(hcondel_genes, all_genes, high_expr)
    
    hcondel_tissue_results.append({
        "tissue": tissue_name,
        "is_brain": tissue_name in BRAIN_TISSUES,
        "n_high_expr": len(high_expr),
        "hcondel_OR": odds,
        "hcondel_p": p,
        "hcondel_n_high": n_high,
    })

# BH FDR
if hcondel_tissue_results:
    pvals = [r["hcondel_p"] for r in hcondel_tissue_results]
    fdrs = bh_fdr(pvals)
    for i, r in enumerate(hcondel_tissue_results):
        r["hcondel_fdr"] = fdrs[i]
        r["hcondel_fdr_sig"] = fdrs[i] < 0.05

print(f"\n  hCONDEL tissue enrichment:")
print(f"  {'Tissue':25s} {'OR':8s} {'FDR':10s} {'Sig':4s}")
for r in sorted(hcondel_tissue_results, key=lambda x: x["hcondel_p"]):
    sig = "*" if r["hcondel_fdr_sig"] else " "
    print(f"  {r['tissue']:25s} {r['hcondel_OR']:8.2f} {r['hcondel_fdr']:10.2e} {sig:>4s}")

# ============================================
# Brain vs non-brain hCONDEL comparison
# ============================================
print(f"\n--- Brain vs non-brain hCONDEL enrichment ---")

brain_results = [r for r in hcondel_tissue_results if r["is_brain"]]
non_brain_results = [r for r in hcondel_tissue_results if not r["is_brain"]]

brain_ors = [r["hcondel_OR"] for r in brain_results if not np.isnan(r["hcondel_OR"])]
non_brain_ors = [r["hcondel_OR"] for r in non_brain_results if not np.isnan(r["hcondel_OR"])]

print(f"  Brain tissues: mean OR={np.mean(brain_ors):.3f}")
print(f"  Non-brain tissues: mean OR={np.mean(non_brain_ors):.3f}")

if len(brain_ors) > 1 and len(non_brain_ors) > 1:
    u, p = stats.mannwhitneyu(brain_ors, non_brain_ors, alternative="greater")
    print(f"  Brain > non-brain: U={u:.0f}, p={p:.4e}")

# ============================================
# Write outputs
# ============================================
print(f"\n--- Writing outputs ---")

hcondel_csv = OUTPUT_DIR / "regulatory_hcondels_tissue_enrichment.csv"
with open(hcondel_csv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(hcondel_tissue_results[0].keys()))
    writer.writeheader()
    for r in hcondel_tissue_results:
        writer.writerow(r)
print(f"  {hcondel_csv}")

# Also write hCONDEL gene list
hcondel_gene_csv = OUTPUT_DIR / "hcondel_gene_mapping.csv"
with open(hcondel_gene_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["gene_id", "gene_symbol", "classification_v7", "mapping_method"])
    for gene_id in sorted(hcondel_genes):
        sym = genes[gene_id].get("gene_symbol", "")
        cls = genes[gene_id]["classification_v7"]
        methods = []
        if gene_id in hcondel_genes_by_name:
            methods.append("name")
        if gene_id in hcondel_genes_by_coord:
            methods.append("coord")
        writer.writerow([gene_id, sym, cls, "+".join(methods)])
print(f"  {hcondel_gene_csv}")

# Summary
summary = {
    "n_hcondels_parsed": len(hcondels),
    "n_with_coord": n_with_coord,
    "n_with_within_gene": n_with_within,
    "n_genes_by_name": len(hcondel_genes_by_name & all_genes),
    "n_genes_by_coord": len(hcondel_genes_by_coord & all_genes),
    "n_genes_total": len(hcondel_genes),
    "n_gd_with_hcondel": len(hcondel_in_gd),
    "n_rd_with_hcondel": len(hcondel_in_rd),
    "n_neutral_with_hcondel": len(hcondel_in_neutral),
    "rd_enrichment": {"OR": float(odds_rd), "p": float(p_rd)},
    "gd_enrichment": {"OR": float(odds_gd), "p": float(p_gd)},
    "brain_mean_OR": float(np.mean(brain_ors)) if brain_ors else None,
    "non_brain_mean_OR": float(np.mean(non_brain_ors)) if non_brain_ors else None,
    "fdr_sig_tissues": [r["tissue"] for r in hcondel_tissue_results if r["hcondel_fdr_sig"]],
}

summary_file = OUTPUT_DIR / "hcondels_summary.json"
with open(summary_file, "w") as f:
    json.dump(summary, f, indent=2)
print(f"  {summary_file}")

print("\n" + "=" * 70)
print("hCONDELs Processing Complete!")
print("=" * 70)
