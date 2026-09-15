# Reproduction Guide — Paper B (Genome Research)

## Environment

```bash
conda env create -f environment.yml   # or: pip install -r requirements.txt
conda activate human_selection
```

Submission-level analyses (figures, Excel tables) need Python only
(pandas, numpy, scipy, matplotlib, openpyxl). The three `phase8_*` analysis
scripts additionally document the original whole-workflow runs and require the
raw public datasets (HPA single-cell atlas, BrainSpan RNA-seq, GTEx v11) not
included in this package — see Data Availability in the manuscript.

## Layout & portability

All scripts resolve the project root via the `HSD_BASE` environment variable
(falling back to the original analysis machine paths). After extracting the
package:

```bash
export HSD_BASE=/path/to/GR_submission_PaperB   # Windows: set HSD_BASE=D:\path\to\...

# Stage the packaged data into the results/ layout the scripts expect:
mkdir -p $HSD_BASE/results/phase8_tissue_analysis/data
mkdir -p $HSD_BASE/results/phase7_gene_vs_regulation/phase7g_classification_v7
cp 07_data/*.csv 07_data/*.json $HSD_BASE/results/phase8_tissue_analysis/
cp 07_data/phase8_master_table.csv $HSD_BASE/results/phase8_tissue_analysis/data/
cp 07_data/gene_classification_v7.csv $HSD_BASE/results/phase7_gene_vs_regulation/phase7g_classification_v7/
```

`07_data/` contains every result table consumed by the figures and the Excel
builder (layer1/layer2 tissue analyses, HPA whole-body single-cell and brain
single-nuclei enrichments, BrainSpan developmental analyses, LOO neuronal
concentration, threshold sensitivity, hCONDEL/HAR mappings, effect-size
summaries) plus the two upstream inputs (the 4,974-gene v7 classification and
the phase-8 master table).

## Run order (submission-level reproduction)

| Step | Script | Input | Output |
|---|---|---|---|
| 1 | `paperB_tables_to_excel.py` | `07_data/` enrichment & LOO tables | `04_tables/Main_Tables_PaperB.xlsx`, `05_supplementary_tables/Supplemental_Tables_PaperB.xlsx` |
| 2 | `paperB_fig1.py` | v7 classification + master table + enrichment CSVs | `Figure1_study_design.png/pdf` (5 panels) |
| 3 | `paperB_fig2.py` | layer1/layer2 GTEx CSVs | `Figure2_bulk_tissue.png/pdf` (5 panels) |
| 4 | `paperB_fig3.py` | HPA whole-body + brain snRNA CSVs/JSONs | `Figure3_single_cell.png/pdf` (6 panels) |
| 5 | `paperB_fig4.py` | BrainSpan + hCONDEL + caMPRA CSVs/JSONs | `Figure4_development_validation.png/pdf` (5+ panels) |

Figures are written to `{HSD_BASE}/results/paperB/figures/` (PNG 300 dpi + PDF).

Original workflow scripts (kept for provenance, not runnable from this package
alone): `phase8_layer2_tau_analysis.py` (per-tissue specificity + tau tests),
`phase8_hpa_single_cell_wholebody.py` (154 cell-type enrichment from the HPA
tsv), `phase8_developmental_analysis.py` (BrainSpan trajectories and
prenatal/postnatal contrast). These require the raw public datasets; their
outputs are exactly the CSVs/JSONs provided in `07_data/`.

## Verification points

Regenerated figures/tables should reproduce these headline values
(independently re-derived during review; see manuscript Table 2):

- Bulk (LOO-tau): GD enriched 30/30 tissues (OR 1.19–1.41); RD depleted 0/30 (OR 0.44–0.96)
- Brain only FDR-significant per-tissue specificity (P = 1.1e-03, FDR = 0.034, Cliff's δ = 0.163)
- Neuronal class OR 3.01 [2.52–3.60], P = 1.8e-33; rank 1/15 in all five variants × three thresholds
- Five-variant two-sided MWU: 3.8e-06 / 4.6e-06 / 5.3e-05 / 9.1e-06 / 9.1e-06
- Brain snRNA neuronal class OR 3.25 [2.72–3.89]; developmental tau P = 0.91; prenatal/postnatal P = 0.57
- hCONDEL RD OR 2.94 [1.85–4.55], P = 6.65e-06; caMPRA (LOO) RD OR 1.84 [1.09–3.02], P = 0.012

## Cross-package dependency

The v7 classification (`07_data/gene_classification_v7.csv`) is the output of
the companion MBE submission (Paper A, `code/analysis/phase7g_v7_classification.py`);
within this package it is consumed as a fixed input.

## Note on `prenatal_vs_postnatal.csv`

`07_data/prenatal_vs_postnatal.csv` is the corrected version (2026-09-14):
the original `ratio` column contained a hidden <0.1-TPM fill rule, replaced by
pure prenatal/postnatal means; the authoritative prenatal/postnatal test
(two-sided, log2 ratio, both windows > 0, LOO-brain) gives P = 0.57.
