# Reproduction Guide

## Environment

```bash
conda env create -f environment.yml   # or: pip install -r requirements.txt
conda activate human_selection
```

External tools required only for the upstream pipeline (Phases 1-4):
HyPhy 2.5.100 (BUSTED/RELAX/MEME/GARD), PAML 4.10.10 (codeml), MAFFT 7.526, BLAST+ 2.17.
The submission-level analyses (classification robustness, tables, figures) need Python only.

## Layout

All scripts resolve the project root via the `HSD_BASE` environment variable
(falling back to the original analysis machine paths). After extracting the package:

```bash
export HSD_BASE=/path/to/MBE_submission_PaperA   # Windows: set HSD_BASE=D:\path\to\...
mv data results   # scripts expect input tables under {HSD_BASE}/results/...
```

The `code/paper/` and `code/analysis/phase7g_v7_classification.py` scripts are fully
portable (HSD_BASE-aware). The `code/pipeline/` scripts and most `code/analysis/phase8_*.py`
scripts document the original workflow and retain original absolute paths
(`/mnt/d/...`, WSL); they require the raw public datasets (Ensembl, GTEx v11, HPA,
BrainSpan, UCSC bigWigs) not included in this package (see Data Availability in the manuscript).

## Run order (submission-level reproduction)

| Step | Script | Input | Output |
|---|---|---|---|
| 1 | `analysis/phase7g_v7_classification.py` | BUSTED/RELAX/conservation/v5 tables | `results/.../gene_classification_v7.csv` (Table S1 basis) |
| 2 | `paper/paperA_sensitivity_analysis.py` | gene_classification_v7.csv | Figures S5-S6, `sensitivity_summary.json`, per-gene stability (feeds Table S1) |
| 3 | `paper/paperA_supplementary_tables.py` | v7 CSV + RELAX + stability | Tables S1-S3 + XLSX |
| 4 | `paper/paperA_figures.py` | v7 + v5 CSVs | Figures S1-S4 |
| 5 | `paper/paperA_revision_v7_analysis.py` | v7 CSV + hCONDEL mapping + stability | Table 3 (`loo_full_matrix.csv`), Figure S7-S8, `revision_v7_summary.json` |
| 6 | `paper/paperA_v8_covariate_sensitivity.py` | v7 CSV + human CDS FASTA + GTEx tau | Table 4 (`covariate_sensitivity.csv`), `kmeans_multistart.json` |

Expected outputs of steps 2, 5 and 6 are provided in `expected_outputs/` for verification
(random seed fixed at 42; results are exactly reproducible).

## Random seeds & determinism

`paperA_sensitivity_analysis.py` uses `numpy.random.default_rng(42)`; all 1,000
perturbation runs are fully deterministic.
