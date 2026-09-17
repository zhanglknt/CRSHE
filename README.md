# CRSHE — Coding Selection Dominates by Prevalence, Regulatory Selection by Specificity

Reproducibility repository for a pair of companion studies of human-lineage positive
selection across 10 primate species (4,974 one-to-one orthologous genes):

- **Paper A (submitted to *Molecular Biology and Evolution*)** — *"Coding Selection Dominates by
  Prevalence, Regulatory Selection by Specificity: A Comparability-Corrected 10-Primate
  Analysis."* Builds the gene-driven (GD) versus regulation-driven (RD)
  classification: 1,214 GD / 52 GD-relaxed / 293 RD / 11 dual-driven / 3,404 unclassified,
  with leave-one-out circularity control and hCONDEL/HAR validation.
  → see [`PaperA_MBE/`](PaperA_MBE/)

- **Paper B (submitted to *Genome Research*)** — *"Cell-type-resolution map of coding versus
  regulatory positive selection in human evolution."* Maps the two classes onto 30 GTEx
  bulk tissues, the Human Protein Atlas whole-body single-cell atlas (154 cell types),
  brain single-nuclei, and BrainSpan development. Headline: coding selection is
  expression-ubiquitous; regulatory selection is depleted in bulk but concentrated in
  neuronal cell types (class OR 3.01, rank 1/15 in all five classification variants).
  → see [`PaperB_GR/`](PaperB_GR/)

## Repository layout

| Directory | Content |
|---|---|
| `PaperA_MBE/code/` | Complete analysis code: `pipeline/` (Phases 1–4: CDS extraction, orthology, MAFFT, BUSTED/RELAX), `analysis/` (Phases 5–8: screening, v7 classification, tissue analyses), `paper/` (sensitivity analyses, tables, figures) |
| `PaperA_MBE/data/` | Input tables required to reproduce the classification and manuscript figures/tables |
| `PaperA_MBE/expected_outputs/` | Reference outputs for verification (random seed 42; exactly reproducible) |
| `PaperB_GR/code/` | Figure scripts, Excel table builder, and core phase-8 analyses (all `HSD_BASE`-aware) |
| `PaperB_GR/data/` | All result tables consumed by the Paper B figures/tables, plus upstream inputs (v7 classification, phase-8 master table) |

## Quick start

Each paper subdirectory contains its own reproduction guide
(`PaperA_MBE/code/REPRODUCTION.md`, `PaperB_GR/code/REPRODUCTION.md`) with
environment files (`environment.yml` / `requirements.txt`), data-staging
instructions, run order, and verification points. In brief:

```bash
conda env create -f PaperA_MBE/code/environment.yml && conda activate human_selection
export HSD_BASE=/path/to/CRSHE          # scripts resolve the root via HSD_BASE
# see the per-paper REPRODUCTION.md for staging and run order
```

## External tools and data

The submission-level analyses (classification, robustness, tables, figures) need
Python only. Re-running the upstream pipeline (Phases 1–4) additionally requires
HyPhy 2.5 (BUSTED/RELAX/MEME/GARD), PAML 4.10, MAFFT 7.5, and BLAST+, plus the
public source datasets (Ensembl, GTEx v11, Human Protein Atlas, BrainSpan, UCSC
conservation tracks, Shin et al. 2024 caMPRA) — see Data Availability in the
manuscripts.

## License

MIT (code). Data files retain the licenses of their upstream sources.

## Citation

Manuscripts under review. A Zenodo DOI for this repository will be provided
upon publication of the first of the two papers.
