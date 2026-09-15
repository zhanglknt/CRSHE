# Figure and Table Code Map (v9)

Every figure and table in the manuscript has explicit generating code. All scripts are in
`code/paper/`; run from a shell with `HSD_BASE` set (see REPRODUCTION.md).
Figures are delivered as PDF (PNGs are also produced alongside); tables are delivered as Excel.

## Main figures (03_figures/, PDF)

| Figure | Package file | Generating script |
|---|---|---|
| Figure 1 (study design & classification framework, 5 panels) | `Figure1_analysis_pipeline.pdf` | `paperA_fig1_v2.py` |
| Figure 2 (coding-selection scan & score construction, 5 panels) | `Figure2_BUSTED_pvalue_distribution.pdf` | `paperA_fig2_v2.py` |
| Figure 3 (classification landscape & robustness, 6 panels) | `Figure3_GDS_RDS_scatter.pdf` | `paperA_fig3_v2.py` |
| Figure 4 (non-circular validation, 6 panels) | `Figure4_LOO_enrichment_matrix.pdf` | `paperA_fig4_v2.py` |
| Figure 5 (independent evidence & cross-study context, 6 panels) | `Figure5_independent_evidence.pdf` | `paperA_fig5.py` |
| Figure 6 (functional content of the classes, 6 panels) | `Figure6_class_functional_content.pdf` | `paperA_fig6.py` |

## Supplementary figures (06_supplementary_figures/, PDF)

| Figure | Package file | Generating script |
|---|---|---|
| Figure S1 (framework comparison) | `FigureS1_framework_comparison.pdf` | `paperA_v9_round4_figures.py` (Figure S1 section) |
| Figure S2 (threshold sensitivity) | `FigureS2_threshold_sensitivity.pdf` | `paperA_sensitivity_analysis.py` (Figures section) |
| Figure S3 (weight perturbation) | `FigureS3_weight_perturbation_stability.pdf` | `paperA_sensitivity_analysis.py` (Figures section) |
| Figure S4 (unsupervised PCA + k-means) | `FigureS4_unsupervised_PCA_kmeans.pdf` | `paperA_v9_round4_figures.py` (Figure S4 section) |

## Main tables (04_tables/, Excel)

| Table | Sheet in `Main_Tables_v9.xlsx` | Source data / generating script |
|---|---|---|
| Table 1 (classification summary) | `Table 1` | `paperA_tables_to_excel.py`; counts from `data/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv` |
| Table 2 (framework comparison) | `Table 2` | `paperA_tables_to_excel.py`; v5 vs v7 CSVs in `data/phase7_gene_vs_regulation/` |
| Table 3 (LOO matrix) | `Table 3` | `paperA_tables_to_excel.py`; `expected_outputs/revision_v7/loo_full_matrix.csv` (from `paperA_revision_v7_analysis.py`) |
| Table 4 (covariate sensitivity) | `Table 4` | `paperA_tables_to_excel.py`; `expected_outputs/revision_v7/covariate_sensitivity.csv` (from `paperA_v8_covariate_sensitivity.py`) |

## Supplementary tables (07_supplementary_tables/, Excel + CSV)

| Table | File / sheet | Generating script |
|---|---|---|
| Table S1 (full 4,974-gene classification) | `PaperA_Supplementary_Tables.xlsx` sheet `S1 Full gene table`; `TableS1_full_gene_classification.csv` | `paperA_supplementary_tables.py` |
| Table S2 (RELAX coverage) | `PaperA_Supplementary_Tables.xlsx` sheet `S2 RELAX coverage`; `TableS2_RELAX_coverage.csv` | `paperA_supplementary_tables.py` |
| Table S3 (scoring framework + rules) | `PaperA_Supplementary_Tables.xlsx` sheets `S3 Scoring framework` + `S3b Rules`; `TableS3_scoring_framework.csv` + rules txt | `paperA_supplementary_tables.py` |
| Table S2-1 (MEME site-level) | `Supplementary_Text_Tables_v9.xlsx` sheet `S2-1 MEME` | `paperA_tables_to_excel.py` |
| Table S3-1 (Tier 1 cross-referenced) | `Supplementary_Text_Tables_v9.xlsx` sheet `S3-1 Tier1` | `paperA_tables_to_excel.py` |
| Table S4-1 (top 10 gene-driven) | `Supplementary_Text_Tables_v9.xlsx` sheet `S4-1 Top GD` | `paperA_tables_to_excel.py` |
| Table S4-2 (top 10 regulation-driven) | `Supplementary_Text_Tables_v9.xlsx` sheet `S4-2 Top RD` | `paperA_tables_to_excel.py` |
| Table S4-3 (dual-driven) | `Supplementary_Text_Tables_v9.xlsx` sheet `S4-3 Dual` | `paperA_tables_to_excel.py` |
