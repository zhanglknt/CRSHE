# -*- coding: utf-8 -*-
"""Build Paper B (Genome Research) Excel tables:
- results/paperB/Main_Tables_PaperB.xlsx    (Table 1, Table 2)
- results/paperB/Supplemental_Tables_PaperB.xlsx (S1-S7)
Follows the Paper A convention: tables as Excel, titles + footnotes in-sheet.
"""
import csv
import json
import os
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment

BASE = os.environ.get("HSD_BASE") or "D:/人类正选择基因项目"
P8 = f"{BASE}/results/phase8_tissue_analysis"

TITLE_FONT = Font(bold=True, size=12)
HDR_FONT = Font(bold=True)


def write_sheet(ws, title, headers, rows, notes=None, num_fmt=None):
    ws["A1"] = title
    ws["A1"].font = TITLE_FONT
    r = 3
    for c, h in enumerate(headers, 1):
        cell = ws.cell(row=r, column=c, value=h)
        cell.font = HDR_FONT
        cell.alignment = Alignment(wrap_text=True, vertical="top")
    for row in rows:
        r += 1
        for c, v in enumerate(row, 1):
            cell = ws.cell(row=r, column=c, value=v)
            if num_fmt and c in num_fmt:
                cell.number_format = num_fmt[c]
    if notes:
        r += 2
        for n in notes:
            ws.cell(row=r, column=1, value=n)
            r += 1
    ws.freeze_panes = "A4"


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


# ================= Table 2 data (two-sided P and Cliff's delta from
# effect_sizes_summary.json; medians/ns/subgroups from
# loo_extended_neuronal_concentration.csv) =================
with open(f"{P8}/effect_sizes_summary.json", encoding="utf-8") as f:
    ES = json.load(f)
CD = ES["cliffs_delta_neuronal_vs_non"]
lo = {}
for row in read_csv(f"{P8}/loo_extended_neuronal_concentration.csv"):
    lo[(row["variant"], row["group"])] = row
# neuronal-class OR per variant from cellclass CSV
cls = {}
for row in read_csv(f"{P8}/hpa_wholebody_cellclass_enrichment.csv"):
    if row["cell_type_class"] == "neuronal cells":
        cls = row

VARMAP = [
    ("v7", "v7 (full model)"),
    ("loo_brain", "LOO-brain"),
    ("loo_tau", "LOO-tau"),
    ("loo_campra", "LOO-caMPRA"),
    ("loo_double", "double-LOO"),
]
P2 = {"v7": 3.8e-06, "loo_brain": 4.6e-06, "loo_tau": 5.3e-05,
      "loo_campra": 9.1e-06, "loo_double": 9.1e-06}
CLASSOR = {"v7": "v7_rd_OR", "loo_brain": "loo_brain_rd_OR",
           "loo_tau": "loo_tau_rd_OR", "loo_campra": "loo_campra_rd_OR",
           "loo_double": "loo_double_rd_OR"}

table2_rows = []
for key, label in VARMAP:
    n9 = lo[(key, "neuronal9")]
    b3 = lo[(key, "brain_neurons3")]
    r6 = lo[(key, "retinal_neurons6")]
    table2_rows.append([
        label, int(n9["gd_n"]), int(n9["rd_n"]),
        float(n9["group_median_rd_OR"]), float(n9["nonneuronal_median_rd_OR"]),
        P2[key], CD[key]["cliffs_delta"], int(n9["n_types_rd_fdr_sig_total"]),
        round(float(cls[CLASSOR[key]]), 2), "1 / 15",
        float(b3["group_median_rd_OR"]), float(r6["group_median_rd_OR"]),
    ])

# ================= Main tables =================
wb = Workbook()
write_sheet(
    wb.active,
    "Table 1. Data resources and key results summary.",
    ["Resource / analysis", "Scope", "Key result"],
    [
        ["Primate framework (companion paper)", "10 species, 4,974-gene universe",
         "GD 1,214; GD-relaxed 52; RD 293; dual 11; unclassified 3,404"],
        ["GTEx (2026-01-16 GCS release, GENCODE 47)", "30 tissues",
         "GD enriched 30/30 (OR 1.19-1.41); RD depleted in all (OR 0.44-0.96); LOO-tau classes"],
        ["Per-tissue specificity", "30 tissues",
         "Brain only FDR-significant tissue (RD > GD, P = 1.1e-03, FDR = 0.034, delta = 0.163)"],
        ["HPA single-cell atlas v25.1 (whole body)", "154 cell types, 15 classes",
         "Neuronal class OR 3.01 [2.52-3.60], P = 1.8e-33; rank 1 of 15 in all five variants"],
        ["HPA brain single-nuclei", "34 cell types (4,963 genes)",
         "RD FDR-significant in 33/34 types; neuronal class OR 3.25 [2.72-3.89]"],
        ["BrainSpan developmental", "26 regions (18 testable), 31 stages",
         "Dev. tau GD 0.59 vs RD 0.60 (P = 0.91); RD nominal in 10/18 regions, none FDR"],
        ["Prenatal vs postnatal", "1,085 GD / 605 RD genes",
         "P = 0.57 (two-sided MWU, log2 ratio); n.s."],
        ["hCONDELs (McLean et al. 2011)", "583 deletions, 183 genes",
         "RD OR 2.94 [1.85-4.55], P = 6.65e-06 (primary independent validation)"],
        ["caMPRA active HARs (Shin et al. 2024)", "508 active of 3,171",
         "RD (LOO-caMPRA) OR 1.84 [1.09-3.02], P = 0.012 (supporting evidence)"],
    ],
    notes=["Abbreviations: GD, gene-driven; RD, regulation-driven; OR, odds ratio; CI, confidence interval.",
           "Classification variants: v7 (full), LOO-brain, LOO-tau, LOO-caMPRA, double-LOO (see Methods)."])

write_sheet(
    wb.create_sheet("Table 2"),
    "Table 2. Neuronal concentration of regulatory selection across five classification variants.",
    ["Variant", "GD n", "RD n", "Neuronal median RD OR (n=9)",
     "Non-neuronal median (n=145)", "Two-sided MWU P", "Cliff's delta",
     "RD FDR<0.05 types (/154)", "Neuronal-class OR", "Class rank",
     "Brain-neuronal median (n=3)", "Retinal-neuronal median (n=6)"],
    table2_rows,
    notes=["Medians are per-cell-type RD enrichment odds ratios; class OR is the Fisher-exact odds ratio",
           "for the aggregated neuronal cell class (15-class analysis). All two-sided P values and",
           "Cliff's delta from effect_sizes_summary.json; per-variant P: v7 3.8e-06, LOO-brain 4.6e-06,",
           "LOO-tau 5.3e-05, LOO-caMPRA 9.1e-06, double-LOO 9.1e-06.",
           "double-LOO removes brain-specificity and tissue-tau components; remaining weights",
           "renormalized (caMPRA 0.545, nc-divergence 0.455). GD excludes gene-driven (relaxed)."])

out1 = f"{BASE}/results/paperB/Main_Tables_PaperB.xlsx"
wb.save(out1)
print("saved", out1)

# ================= Supplemental tables =================
wb2 = Workbook()
first = True


def sheet(name, title, csv_path, cols=None, notes=None):
    global first
    ws = wb2.active if first else wb2.create_sheet(name)
    first = False
    ws.title = name[:31]
    rows = read_csv(csv_path)
    headers = cols or list(rows[0].keys())
    data = [[r.get(h, "") for h in headers] for r in rows]
    for r in data:  # numeric conversion where possible
        for i, v in enumerate(r):
            try:
                r[i] = float(v)
            except (ValueError, TypeError):
                pass
    write_sheet(ws, title, headers, data, notes=notes)


sheet("S1a GTEx v7", "Supplemental Table S1a. GTEx 30-tissue enrichment, v7 classification.",
      f"{P8}/layer1_tissue_enrichment.csv",
      cols=["tissue", "is_brain", "gd_OR", "gd_p", "gd_fdr", "gd_fdr_sig",
            "rd_OR", "rd_p", "rd_fdr", "rd_fdr_sig", "gd_n_high", "rd_n_high"],
      notes=["v7 shown for comparison; primary analyses use the LOO-tau variant (S1b)."])
sheet("S1b GTEx LOO-tau", "Supplemental Table S1b. GTEx 30-tissue enrichment, LOO-tau classification.",
      f"{P8}/layer1_tissue_enrichment_loo_tau.csv",
      cols=["tissue", "is_brain", "gd_OR", "gd_p", "gd_fdr", "gd_fdr_sig",
            "rd_OR", "rd_p", "rd_fdr", "rd_fdr_sig", "gd_n_high", "rd_n_high"])
sheet("S2 cell types", "Supplemental Table S2. Whole-body 154 cell-type enrichment, five variants.",
      f"{P8}/hpa_wholebody_celltype_enrichment.csv",
      notes=["Per-type Fisher-exact ORs for GD and RD under v7, LOO-brain, LOO-tau, LOO-caMPRA,",
             "and double-LOO classifications; BH FDR per variant. Brain- and retinal-neuronal",
             "subgroup medians are summarized in main Table 2."])
sheet("S3 cell classes", "Supplemental Table S3. 15 cell-class enrichment (GD and RD) across variants.",
      f"{P8}/hpa_wholebody_cellclass_enrichment.csv",
      notes=["Class-level aggregation treats the cell class as the independent unit",
             "(primary single-cell evidence; see Methods)."])
sheet("S4 brain snRNA", "Supplemental Table S4. HPA brain single-nuclei, 34 cell-type enrichment.",
      f"{P8}/hpa_single_nuclei_enrichment.csv",
      notes=["4,963 genes matched; classification LOO-brain (see Methods)."])
sheet("S5a BrainSpan regions", "Supplemental Table S5a. BrainSpan region enrichment (18 testable regions).",
      f"{P8}/brainspan_region_enrichment.csv",
      notes=["RD nominally significant in 10/18 regions (A1C, DFC, IPC, ITC, M1C, MFC, OFC, S1C,",
             "STC, VFC; weakest IPC/STC); no region FDR-significant for either class."])
sheet("S5b BrainSpan periods", "Supplemental Table S5b. BrainSpan developmental-period enrichment.",
      f"{P8}/developmental_period_enrichment.csv")
sheet("S6a hCONDEL genes", "Supplemental Table S6a. hCONDEL-overlapping genes and classifications (v7).",
      f"{P8}/hcondel_gene_mapping.csv",
      notes=["183 genes: 27 RD, 27 GD, 125 unclassified, 2 dual (ASAP3, CNTN4),",
             "2 GD-relaxed (ST6GALNAC5, SNX27)."])
sheet("S6b HAR genes", "Supplemental Table S6b. HAR-proximate genes, caMPRA activity, classifications.",
      f"{P8}/har_gene_mapping.csv",
      notes=["has_active_har = Y indicates >=1 caMPRA-active HAR within +/-50 kb (Shin et al. 2024)."])
sheet("S7 threshold sens.", "Supplemental Table S7. Threshold sensitivity of single-cell enrichment.",
      f"{P8}/threshold_sensitivity_by_class.csv",
      notes=["High expression defined as non-zero nCPM >= 90th / 75th / 50th percentile (top-10/25/50%).",
             "The neuronal class ranks 1 of 15 under every threshold x variant combination;",
             "cell-type MWU significant at all thresholds (P < 1e-05)."])

out2 = f"{BASE}/results/paperB/Supplemental_Tables_PaperB.xlsx"
wb2.save(out2)
print("saved", out2)
