#!/usr/bin/env python3
"""
Phase 8: Build Master Table

Merges v7 classification + GTEx V11 tissue TPM + BUSTED v1 omega + 
HPA multi-source tau + RELAX + Conservation into a single master table.

Run in WSL: python3 phase8_build_master_table.py
"""
import csv
import os
import sys
import json
import numpy as np
from pathlib import Path

BASE = Path("/mnt/d/hs_gene_project")
CLASSIFICATION_CSV = BASE / "results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv"
GTEX_TPM_CSV = BASE / "results/gtex_v11_tau/gtex_v11_tissue_median_tpm.csv"
GTEX_TAU_CSV = BASE / "results/gtex_v11_tau/gtex_v11_tau.csv"
BUSTED_OMEGA_CSV = BASE / "results/phase8_tissue_analysis/busted_v1_omega.csv"
HPA_TAU_CSV = BASE / "results/phase7_gene_vs_regulation/phase7j_tissue_analysis/data/hpa_tau_values.csv"
OUTPUT_DIR = BASE / "results/phase8_tissue_analysis/data"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_CSV = OUTPUT_DIR / "phase8_master_table.csv"

print("=" * 70)
print("Phase 8: Build Master Table")
print("=" * 70)

# ============================================
# Step 1: Load v7 classification (base table)
# ============================================
print("\n--- Step 1: Loading v7 classification ---")

gene_data = {}
fieldnames = []
with open(CLASSIFICATION_CSV) as f:
    reader = csv.DictReader(f)
    fieldnames = list(reader.fieldnames)
    for row in reader:
        gid = row["gene_id"]
        gene_data[gid] = dict(row)

print(f"  Loaded {len(gene_data)} genes, {len(fieldnames)} columns")

# ============================================
# Step 2: Add GTEx V11 tissue TPM (30 tissues)
# ============================================
print("\n--- Step 2: Adding GTEx V11 tissue TPM ---")

gtex_tissues = []
gtex_count = 0
with open(GTEX_TPM_CSV) as f:
    reader = csv.DictReader(f)
    gtex_tissues = [h for h in reader.fieldnames if h != "gene_id"]
    for row in reader:
        gid = row["gene_id"]
        if gid in gene_data:
            for tissue in gtex_tissues:
                try:
                    gene_data[gid][f"gtex_{tissue}"] = float(row[tissue])
                except (ValueError, TypeError):
                    gene_data[gid][f"gtex_{tissue}"] = 0.0
            gtex_count += 1

# Add tau from GTEx V11
with open(GTEX_TAU_CSV) as f:
    reader = csv.DictReader(f)
    for row in reader:
        gid = row["gene_id"]
        if gid in gene_data:
            try:
                gene_data[gid]["gtex_v11_tau"] = float(row["tau"])
            except (ValueError, TypeError):
                gene_data[gid]["gtex_v11_tau"] = 0.0

print(f"  Added GTEx V11 TPM for {gtex_count} genes ({len(gtex_tissues)} tissues)")

# ============================================
# Step 3: Add BUSTED v1 omega
# ============================================
print("\n--- Step 3: Adding BUSTED v1 omega ---")

omega_count = 0
with open(BUSTED_OMEGA_CSV) as f:
    reader = csv.DictReader(f)
    for row in reader:
        gid = row["gene_id"]
        if gid in gene_data:
            gene_data[gid]["omega_max"] = float(row["omega_max"])
            gene_data[gid]["omega_weighted"] = float(row["omega_weighted"])
            gene_data[gid]["omega_positive"] = float(row["omega_positive"])
            gene_data[gid]["proportion_positive"] = float(row["proportion_positive"])
            gene_data[gid]["has_positive_selection"] = row["has_positive_selection"] == "True"
            gene_data[gid]["bg_omega_max"] = float(row["bg_omega_max"])
            gene_data[gid]["bg_omega_weighted"] = float(row["bg_omega_weighted"])
            omega_count += 1

print(f"  Added BUSTED v1 omega for {omega_count} genes")

# ============================================
# Step 4: Add HPA multi-source tau
# ============================================
print("\n--- Step 4: Adding HPA multi-source tau ---")

hpa_count = 0
with open(HPA_TAU_CSV) as f:
    reader = csv.DictReader(f)
    for row in reader:
        gid = row["ensembl_id"]
        if gid in gene_data:
            gene_data[gid]["hpa_tau_gtex"] = float(row.get("tau_gtex", 0)) if row.get("tau_gtex") else 0
            gene_data[gid]["hpa_tau_hpa"] = float(row.get("tau_hpa_own", 0)) if row.get("tau_hpa_own") else 0
            gene_data[gid]["hpa_tau_consensus"] = float(row.get("tau_consensus", 0)) if row.get("tau_consensus") else 0
            gene_data[gid]["hpa_max_tissue_gtex"] = row.get("max_tissue_gtex", "")
            gene_data[gid]["hpa_max_tissue_hpa"] = row.get("max_tissue_hpa_own", "")
            gene_data[gid]["hpa_max_tissue_consensus"] = row.get("max_tissue_consensus", "")
            gene_data[gid]["hpa_brain_mean_gtex"] = float(row.get("brain_mean_gtex", 0)) if row.get("brain_mean_gtex") else 0
            gene_data[gid]["hpa_non_brain_mean_gtex"] = float(row.get("non_brain_mean_gtex", 0)) if row.get("non_brain_mean_gtex") else 0
            gene_data[gid]["hpa_brain_mean_hpa"] = float(row.get("brain_mean_hpa_own", 0)) if row.get("brain_mean_hpa_own") else 0
            gene_data[gid]["hpa_non_brain_mean_hpa"] = float(row.get("non_brain_mean_hpa_own", 0)) if row.get("non_brain_mean_hpa_own") else 0
            gene_data[gid]["hpa_brain_mean_consensus"] = float(row.get("brain_mean_consensus", 0)) if row.get("brain_mean_consensus") else 0
            gene_data[gid]["hpa_non_brain_mean_consensus"] = float(row.get("non_brain_mean_consensus", 0)) if row.get("non_brain_mean_consensus") else 0
            hpa_count += 1

print(f"  Added HPA multi-source tau for {hpa_count} genes")

# ============================================
# Step 5: Build output fieldnames
# ============================================
print("\n--- Step 5: Building output fieldnames ---")

# Start with v7 fields
output_fields = list(fieldnames)

# Add GTEx tissue TPM columns
for tissue in gtex_tissues:
    output_fields.append(f"gtex_{tissue}")
output_fields.append("gtex_v11_tau")

# Add BUSTED v1 omega columns
for col in ["omega_max", "omega_weighted", "omega_positive", "proportion_positive",
            "has_positive_selection", "bg_omega_max", "bg_omega_weighted"]:
    if col not in output_fields:
        output_fields.append(col)

# Add HPA columns
for col in ["hpa_tau_gtex", "hpa_tau_hpa", "hpa_tau_consensus",
            "hpa_max_tissue_gtex", "hpa_max_tissue_hpa", "hpa_max_tissue_consensus",
            "hpa_brain_mean_gtex", "hpa_non_brain_mean_gtex",
            "hpa_brain_mean_hpa", "hpa_non_brain_mean_hpa",
            "hpa_brain_mean_consensus", "hpa_non_brain_mean_consensus"]:
    if col not in output_fields:
        output_fields.append(col)

print(f"  Total output columns: {len(output_fields)}")

# ============================================
# Step 6: Write master table
# ============================================
print("\n--- Step 6: Writing master table ---")

with open(OUTPUT_CSV, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=output_fields, extrasaction="ignore")
    writer.writeheader()
    for gid in sorted(gene_data.keys()):
        writer.writerow(gene_data[gid])

print(f"  Output: {OUTPUT_CSV}")
print(f"  Genes: {len(gene_data)}")
print(f"  Columns: {len(output_fields)}")

# ============================================
# Summary
# ============================================
print(f"\n--- Summary ---")
print(f"  Total genes: {len(gene_data)}")
print(f"  Total columns: {len(output_fields)}")
print(f"  GTEx V11 tissues: {len(gtex_tissues)}")
print(f"  With BUSTED v1 omega: {omega_count} ({100*omega_count/len(gene_data):.1f}%)")
print(f"  With HPA tau: {hpa_count} ({100*hpa_count/len(gene_data):.1f}%)")

# Classification counts
from collections import Counter
cls_counts = Counter(g["classification_v7"] for g in gene_data.values())
print(f"\n  Classification v7:")
for k, v in cls_counts.most_common():
    print(f"    {k}: {v}")

cls_loo_brain_counts = Counter(g.get("classification_loo_brain", "unknown") for g in gene_data.values())
print(f"\n  Classification LOO-brain:")
for k, v in cls_loo_brain_counts.most_common():
    print(f"    {k}: {v}")

# GTEx tissue list
print(f"\n  GTEx V11 tissues ({len(gtex_tissues)}):")
brain_tissues = ["Brain", "Nerve"]
brain_set = [t for t in gtex_tissues if any(b in t for b in brain_tissues)]
print(f"    Brain-related: {brain_set}")

print("\n" + "=" * 70)
print("Master Table Complete!")
print("=" * 70)
