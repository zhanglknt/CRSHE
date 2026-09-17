# Final Blind Verification — Paper A (MBE) Round-6 New Numbers

**Verifier:** verifier-final (independent) · **Date:** 2026-09-17
**Method:** all numbers recomputed from source data only (classification CSV components, raw Selectome NHX TSV, phase2 strict-ortholog table, phase8 hCONDEL mapping, tier3 sweep CSV, p2 GO-slim counts table). No analysis-JSON conclusions were imported; BUSTED BH-FDR was re-derived from raw p-values (0 mismatches vs stored column; 1,690/4,974 sig).

## A. Selectome correction — **PASS**
- Raw TSV: 96 rows → **94** unique gene_ids (2 rows lack gene_id); **24** fall in the 4,974 universe (symbol list matches).
- Class distribution of the 24: GD **6**/1,214, GDr **1**/52, dual **0**/11, RD **0**/293, neutral **17**/3,404.
- 2×2 = **[[7,1259],[17,3691]]**; Fisher OR = **1.2072** (95% CI **0.422–3.072**, conditional); one-sided greater P = **0.4127**; two-sided P = **0.6432**. All match (1.21, 0.42–3.07, 0.41, 0.64).
- Dual-method (Selectome ∧ BUSTED FDR<0.05): **9 genes** — CCDC122, CCHCR1, CD14, MRPS5, NUB1, PTER, RAD17, SLC38A9, TMEM171. Confirmed. (CD14/PTER have blank gene_symbol in the classification CSV; identity confirmed via ENSG ↔ Selectome gene_name.)

## B. Selectome sensitivity — **PASS**
- Stored flag 16 genes; new NHX flag 24; intersection **1** (C6orf118); symmetric difference **38**. ✓
- GDS = 0.35·p_pct + 0.30·lrt_pct + 0.20·relax_pct + 0.15·selectome verified from columns (max|Δ| 8e-7); classification re-derived 100% per-gene.
- With corrected flags: counts **1,215 / 52 / 293 / 11 / 3,403**; exactly **3 flips**: TXNDC16 GD→neutral (GDS 0.3893→0.2393), MRPS5 neutral→GD (0.1747→0.3247), PTER neutral→GD (0.3434→0.4934). All match.

## C. Tier 3 coverage sensitivity — **PASS**
- RD(n_pos) OLS on published sweep: slope **0.3270**, intercept **261.18**, r **0.99986**. ✓
- Crossing n* = (1,214−261.18)/0.3270 = **2,913.9 → 2,914 (58.58%)**. ✓
- GD fit slope **−0.1452**, intercept **1,279.99**; equal-coverage point **2,157.7 → 2,158 (43.38%)**. ✓
- Bootstrap (replicate-level, one draw per level per series): original seed 20260915 reproduces CI **[2,741, 3,146]** and **[2,019, 2,277]** exactly; independent seed gives [2,720, 3,131] and [2,048, 2,267] — consistent.
- Endpoint "measured" values ARE independently reproducible (deterministic classifier ceilings): doan=1 for all 4,974 → RD **1,906**; doan=1 for the 476 HAR-proximate → RD **410**. Line predictions: **416.8 ≈ 417** and **1,887.6 ≈ 1,888**. ✓

## D. Strict vs permissive — **PASS**
- From shared_genes_one2one_all10.csv (6,796 strict ids): strict **4,119** / permissive **855**. ✓
- BUSTED FDR<0.05: **1,345/4,119 = 32.65%** vs **345/855 = 40.35%**. ✓
- Fisher (2×2 [[1345,2774],[345,510]]): OR = **0.7167**, two-sided P = **2.129e-05**. Matches (0.72, 2.1e-05).

## E. Tier 4 coding-selection rates — **PASS** (one method-label caveat)
- HAR-proximate (n_hars>0): n **476**, sig **174**, rate **36.6%**; diff vs genome **+2.6pp**. ✓
- hCONDEL (id/symbol match to hcondel_gene_mapping.csv): n **183**, sig **56**, rate **30.6%**; diff **−3.4pp**. ✓
- Genome-wide **1,690/4,974 = 34.0%**; combined (HAR∪hCONDEL) **220/630 = 34.9%**, **+0.9pp**. ✓
- **Caveat:** the reported 95% CIs (−2.9~+8.3, −10.9~+5.0, −4.0~+6.1) reproduce only as the *difference of Wilson score bounds* (l₁−u₂, u₁−l₂) — the construction described in tier4_empirical_disjointness.json but there labeled "Newcombe hybrid score". The true Newcombe (1998) hybrid score (square-and-add) interval is narrower: **[−1.8, +7.2], [−9.7, +3.8], [−2.9, +5.0]**. All variants include 0, so the "indistinguishable" conclusion is unaffected. Recommend renaming the method or citing the true values.

## F. GO-slim GD-class q — **PASS** (rounding note)
- 40 tests (4 classes × 10 buckets); one-sided Fisher p recomputed from k_overlap/n_class/K_universe — matches stored p to 1e-16; BH re-derived.
- gene_driven minimum q = **0.4583** (Metabolism & proteostasis) — exact match.
- Note: manuscript wording "all q ≥ 0.46" is literally false (0.4583 < 0.46); it holds only after rounding. Suggest "minimum q = 0.46 (rounded)" or "all q > 0.45".

## G. LOO counts — **PASS**
- Recomputed all five variants from component columns with documented weight renormalizations (full 0.30/0.25/0.25/0.20; LOO-caMPRA 0/0.35/0.35/0.30; LOO-tau 0.40/0/0.33/0.27; LOO-nc 0.40/0.33/0/0.27; LOO-brain 0.375/0.3125/0.3125/0).
- RD counts: full **293**, LOO-caMPRA **800**, LOO-tau **148**, LOO-nc **110**, LOO-brain **613**. All match; every other class count in loo_full_matrix.csv also reproduced exactly (e.g., LOO-nc 1,369/71/110/8/3,416).

## Cross-cutting finding (submission package)
**MBE_submission_PaperA/01_main_text/manuscript_english_v9.md still contains the OLD, superseded Selectome numbers** ("Eight genes … Fisher OR = 2.94, P = 0.030"; "all q ≥ 0.38"), whereas the corrected numbers are in **results/paper/manuscript_english_v9.md**. The pending package rebuild (task #89) must use the corrected manuscript or the retracted numbers will ship.

## Verdict
**7/7 items PASS** at the claimed precision. Two minor wording/method-label issues (E CI label, F "≥0.46" rounding) and one packaging risk (stale submission copy) — none affects any scientific conclusion.

*Verification scripts and raw outputs: `results/phase9_hardening/verification/vf_*.py`, `vf_*_out.txt`.*
