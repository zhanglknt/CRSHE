#!/usr/bin/env python3
"""
Paper A (MBE) — P0: Supplementary Tables S1-S3
================================================
S1: Full gene classification table (4,974 genes × all component scores + stability)
S2: RELAX coverage details (1,701 genes with K/p/LRT + universe coverage flags)
S3: GDS/RDS scoring framework (formula + weights + provenance)

Output: results/paper/supplementary/  (CSV per table + single XLSX workbook)

Run (Windows): python.exe scripts/paperA_supplementary_tables.py
"""
import json
import os
import numpy as np
import pandas as pd

if os.path.exists("/mnt/d"):
    BASE = "/mnt/d/人类正选择基因项目"
else:
    BASE = "D:/人类正选择基因项目"

OUT_DIR = f"{BASE}/results/paper/supplementary"
os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 70)
print("Paper A — Supplementary Tables S1-S3")
print("=" * 70)

# ------------------------------------------------------------------
# S1: Full gene table (readable column names) + stability
# ------------------------------------------------------------------
print("\n--- S1: Full gene classification table ---")
v7 = pd.read_csv(f"{BASE}/results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv",
                 low_memory=False)
stab = pd.read_csv(f"{BASE}/results/paper/sensitivity/per_gene_stability.csv")

RENAME = {
    "gene_id": "Ensembl Gene ID",
    "gene_symbol": "Gene Symbol",
    "cds_length": "CDS Length (bp)",
    "classification_v7": "Classification (v7)",
    "busted_p": "BUSTED p-value",
    "busted_lrt": "BUSTED LRT",
    "bh_fdr": "BH FDR q-value",
    "bh_fdr_sig": "BH FDR < 0.05",
    "storey_q": "Storey q-value",
    "storey_sig": "Storey q < 0.05",
    "relax_K": "RELAX K",
    "relax_p": "RELAX p-value",
    "relax_fdr_sig": "RELAX FDR < 0.05",
    "has_relax_active": "RELAX Result Available",
    "relax_intensified": "RELAX Intensified (K>1, FDR sig)",
    "relax_relaxed": "RELAX Relaxed (K<1, FDR sig)",
    "has_selectome": "Selectome Positive Selection",
    "n_hars": "HARs within ±50kb",
    "has_doan_campra": "caMPRA Active HAR (Shin et al. 2024)",
    "tau": "GTEx tau (tissue specificity)",
    "brain_tau": "BrainSpan tau (brain specificity)",
    "phyloP_cds": "phyloP CDS Mean",
    "phyloP_gene": "phyloP Gene Body Mean",
    "phastCons_cds": "phastCons CDS Mean",
    "phastCons_gene": "phastCons Gene Body Mean",
    "nc_conservation": "Non-coding Conservation Index",
    "gds_p_pct_resid": "GDS Component: BUSTED p (CDS-residualized percentile)",
    "gds_lrt_pct_resid": "GDS Component: BUSTED LRT (CDS-residualized percentile)",
    "gds_relax_pct": "GDS Component: RELAX K Score",
    "gds_selectome": "GDS Component: Selectome",
    "gds_v7": "GDS Score (v7)",
    "rds_doan": "RDS Component: caMPRA",
    "rds_tau": "RDS Component: tau Percentile",
    "rds_nc": "RDS Component: Non-coding Acceleration Percentile",
    "rds_brain": "RDS Component: Brain tau Percentile",
    "rds_v7": "RDS Score (v7)",
    "classification_loo_campra": "Classification (LOO-caMPRA)",
    "classification_loo_tau": "Classification (LOO-tau)",
    "classification_loo_brain": "Classification (LOO-brain)",
    "stability": "Classification Stability (1000 perturbations)",
}
s1 = v7.rename(columns=RENAME).copy()
if "stability" not in s1.columns:
    s1 = s1.merge(stab[["gene_id", "stability"]].rename(
        columns={"gene_id": "Ensembl Gene ID", "stability": "Classification Stability (1000 perturbations)"}),
        on="Ensembl Gene ID", how="left")

# P0-D correction: the 'Selectome Positive Selection' column is overwritten with
# the corrected Selectome v6 NHX primate-ancestral flag (24 genes in the universe).
# The v9 classification and GDS component columns are intentionally left as
# published (the 3-gene sensitivity is quantified in the manuscript Limitations).
_sel_nhx = pd.read_csv(f"{BASE}/data/selectome/selectome_primate_positive_selection.tsv", sep="\t")
_nhx_ids = set(str(g).split(".")[0] for g in _sel_nhx["gene_id"].dropna())
_n_changed = int((s1["Selectome Positive Selection"].astype(bool)
                  != s1["Ensembl Gene ID"].map(lambda g: str(g).split(".")[0] in _nhx_ids)).sum())
s1["Selectome Positive Selection"] = s1["Ensembl Gene ID"].map(
    lambda g: str(g).split(".")[0] in _nhx_ids)
print(f"  S1 Selectome column corrected to NHX flag: {_nhx_ids and len(_nhx_ids)} source genes, "
      f"{_n_changed} rows changed")
print(f"  S1: {len(s1)} genes × {len(s1.columns)} columns")
s1.to_csv(f"{OUT_DIR}/TableS1_full_gene_classification.csv", index=False)

# ------------------------------------------------------------------
# S2: RELAX coverage details
# ------------------------------------------------------------------
print("\n--- S2: RELAX coverage details ---")
relax = pd.read_csv(f"{BASE}/results/phase4_hyphy/relax_results/relax_results_v3.csv")
universe = v7[["gene_id", "gene_symbol", "classification_v7"]].copy()

relax_d = relax.drop_duplicates("gene_id")
s2 = universe.merge(relax_d, on="gene_id", how="left")
s2["RELAX Run Available"] = s2["K"].notna()
s2["RELAX Default (K=1, p≥0.99)"] = s2["RELAX Run Available"] & (
    (s2["p_value"].fillna(1.0) >= 0.99) & (s2["K"].fillna(1.0) == 1.0))
s2["RELAX Coverage Category"] = np.where(
    ~s2["RELAX Run Available"], "No RELAX result (run failed/not attempted)",
    np.where(s2["RELAX Default (K=1, p≥0.99)"], "Default values (non-informative)",
             "Active result (K≠1 or p<0.99)"))

s2 = s2.rename(columns={
    "gene_id": "Ensembl Gene ID", "gene_symbol": "Gene Symbol",
    "classification_v7": "Classification (v7)",
    "K": "RELAX K", "LRT": "RELAX LRT", "p_value": "RELAX p-value",
    "fdr": "RELAX BH FDR", "fdr_significant": "RELAX FDR < 0.05"})
s2 = s2[["Ensembl Gene ID", "Gene Symbol", "Classification (v7)",
         "RELAX Run Available", "RELAX Coverage Category",
         "RELAX K", "RELAX LRT", "RELAX p-value", "RELAX BH FDR", "RELAX FDR < 0.05"]]
cov = s2["RELAX Coverage Category"].value_counts()
print(f"  S2: {len(s2)} genes; coverage: {dict(cov)}")
s2.to_csv(f"{OUT_DIR}/TableS2_RELAX_coverage.csv", index=False)

# ------------------------------------------------------------------
# S3: Scoring framework
# ------------------------------------------------------------------
print("\n--- S3: GDS/RDS scoring framework ---")
s3 = pd.DataFrame([
    # GDS
    {"Score": "GDS", "Component": "BUSTED p-value (CDS-length residualized, percentile)",
     "Weight": 0.35, "Data Source": "HyPhy BUSTED v3 (10-species alignment)",
     "Notes": "Residual of -log10(p) vs log10(CDS length); percentile rank"},
    {"Score": "GDS", "Component": "BUSTED LRT (CDS-length residualized, percentile)",
     "Weight": 0.30, "Data Source": "HyPhy BUSTED v3",
     "Notes": "Residual of LRT vs log10(CDS length)"},
    {"Score": "GDS", "Component": "RELAX K score", "Weight": 0.20,
     "Data Source": "HyPhy RELAX v3",
     "Notes": "Genes without active RELAX result assigned neutral 0.5"},
    {"Score": "GDS", "Component": "Selectome positive selection", "Weight": 0.15,
     "Data Source": "Selectome v6, primate ancestral branches (NHX parse)",
     "Notes": "Binary 0/1; see manuscript Limitations for flag-provenance sensitivity"},
    # RDS
    {"Score": "RDS", "Component": "caMPRA active HAR (Doan 2024)", "Weight": 0.30,
     "Data Source": "Shin et al. 2024 caMPRA (508 active HARs)",
     "Notes": "Binary: ≥1 active HAR within ±50kb"},
    {"Score": "RDS", "Component": "GTEx tau percentile", "Weight": 0.25,
     "Data Source": "GTEx v11 (30 tissues)",
     "Notes": "Tissue-specificity percentile rank"},
    {"Score": "RDS", "Component": "Non-coding acceleration percentile", "Weight": 0.25,
     "Data Source": "phyloP 100-way (UCSC)",
     "Notes": "Percentile of -(phyloP gene body − phyloP CDS)"},
    {"Score": "RDS", "Component": "Brain tau percentile", "Weight": 0.20,
     "Data Source": "BrainSpan (26 regions)",
     "Notes": "Brain-region-specificity percentile rank"},
])
s3.to_csv(f"{OUT_DIR}/TableS3_scoring_framework.csv", index=False)
print(f"  S3: {len(s3)} rows (4 GDS + 4 RDS components)")

# Classification rules (as a text block saved alongside)
rules = """Classification rules (Phase 7G v7; threshold = 0.15):
  gene-driven  : GDS > RDS + 0.15 AND BH FDR(BUSTED) < 0.05
  gene-driven (relaxed) : gene-driven AND RELAX K<1 with FDR<0.05 (relaxed constraint)
  regulation-driven : RDS > GDS + 0.15 AND BUSTED BH FDR >= 0.05
  dual-driven  : (GDS > RDS + 0.15 AND RDS > 0.4 AND FDR sig) OR (RDS > GDS + 0.15 AND FDR sig AND GDS > 0.3)
  neutral      : otherwise

Weight sensitivity (1000 joint perturbations, weights ~ U(0.5,1.5)×baseline, renormalized):
  mean agreement with baseline = 94.6% (SD 2.6%); GD median per-gene stability = 1.000;
  RD = 0.833; stable core (stability >= 0.8) = 87.9% of genes.
Leave-one-out variants (non-circular validation): LOO-caMPRA, LOO-tau, LOO-brain
  (weights of remaining components renormalized to 1.0)."""
with open(f"{OUT_DIR}/TableS3_classification_rules.txt", "w", encoding="utf-8") as f:
    f.write(rules)

# ------------------------------------------------------------------
# Combined XLSX workbook
# ------------------------------------------------------------------
print("\n--- Combined XLSX workbook ---")
xlsx_path = f"{OUT_DIR}/PaperA_Supplementary_Tables.xlsx"
with pd.ExcelWriter(xlsx_path, engine="openpyxl") as xw:
    s1.to_excel(xw, sheet_name="S1 Full gene table", index=False)
    s2.to_excel(xw, sheet_name="S2 RELAX coverage", index=False)
    s3.to_excel(xw, sheet_name="S3 Scoring framework", index=False)
    rules_df = pd.DataFrame({"Classification rules and sensitivity": rules.split("\n")})
    rules_df.to_excel(xw, sheet_name="S3b Rules", index=False)
print(f"  Written: {xlsx_path}")

# Summary JSON
summary = {
    "S1_genes": int(len(s1)), "S1_columns": int(len(s1.columns)),
    "S2_relax_coverage": {k: int(v) for k, v in cov.items()},
    "S3_components": {"GDS": {"p": 0.35, "LRT": 0.30, "RELAX_K": 0.20, "Selectome": 0.15},
                      "RDS": {"caMPRA": 0.30, "tau": 0.25, "nc_accel": 0.25, "brain_tau": 0.20}},
}
with open(f"{OUT_DIR}/supplementary_tables_summary.json", "w") as f:
    json.dump(summary, f, indent=2)

print("\n" + "=" * 70)
print("Supplementary tables complete")
print("=" * 70)
