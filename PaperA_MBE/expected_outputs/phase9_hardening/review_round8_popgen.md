# Review — Round 8 (Blind), Population & Quantitative Genetics

**Manuscript:** Coding Selection Dominates by Prevalence, Regulatory Selection by Specificity: A Comparability-Corrected 10-Primate Analysis (v9)
**Reviewer expertise:** GWAS architecture, polygenic traits, gene-mapping conventions, confounding control, power of selection tests
**Score:** 6.5 / 10
**Recommendation:** Major Revision

---

## Summary

The authors compare coding versus regulatory modes of human-specific positive selection across 4,974 orthologs in 10 primates, combining BUSTED/RELAX coding tests with a composite regulatory score (RDS), leave-one-out (LOO) validation against HARs and hCONDELs, universe-background functional enrichment, and GWAS Catalog anchoring across 8 neuropsychiatric/cognitive traits plus 4 negative controls. The manuscript is exceptionally self-critical: the RELAX run-list confound, the HAR–caMPRA nesting, the coverage asymmetry (2–10% vs. 100%), and — most importantly — the gene-length confounding of the GWAS anchoring are all disclosed in the main text rather than buried. The length-matched and logistic-corrected SynGO/GO-BP enrichment (7/11 SynGO, 18/18 GO BP surviving) and the hCONDEL LOO validation are the strongest results. My main concerns are (1) the GWAS anchoring, which by the authors' own analysis is a long-gene/mappability phenomenon yet still carries the "specificity" narrative in the title, abstract and Fig. 7c; (2) an asymmetry in length handling — GDS is residualized on CDS length at the score level, RDS is never residualized on gene span, so the classification itself embeds the very confound that downstream corrections then chase; and (3) several loose ends around the 48-test one-sided BH family, the undocumented joint-LOO weight scheme, and a numeric inconsistency in the MDD depletion q-value range. None of these is fatal, but the first two require analysis-level work, not rewording.

---

## Major comments

### M1. The GWAS anchoring is a gene-length artifact by the authors' own showing — it should not carry the "specificity" claim, and the residual EA/SCZ signal needs a length- and LD-aware gene-set test

The unadjusted headline (all 8 neuropsychiatric/cognitive traits enriched in RD, OR 1.9–4.0, GD null) appears in the abstract, the Results ("External disease and trait anchoring"), Fig. 7c, and the concluding remarks. The same manuscript then demonstrates that:

- two of four negative controls (height OR ≈ 2.2, T2D OR ≈ 2.1) are enriched in RD at magnitudes overlapping the neuro traits (intelligence OR = 1.95 is *smaller* than height's 2.21);
- logistic adjustment on log₁₀ gene span abolishes the enrichment for **all 12 traits** (all adjusted q > 0.44) **and** abolishes the neuro-versus-control family difference (adjusted pooled OR 1.23 vs. 1.06, P = 0.24);
- only coarse length-tertile CMH and alternative mapping conventions retain educational attainment and schizophrenia.

Given this, the trait anchoring cannot support "Regulatory Selection by Specificity" (title) as trait evidence. The abstract does disclose the abolition (good), but the first sentence of the GWAS claim ("enriched for all eight … with gene-driven genes null throughout") is what a reader — and a citation — will take away, and Fig. 7c visually presents the unadjusted forest plot as the anchoring result with the one-sided 48-test q as canonical. The honest framing, which the Discussion mostly reaches, is: *the anchoring is a long-gene mappability phenomenon; the length-robust specificity evidence is the synaptic gene-set enrichment.* I ask that:

1. Fig. 7c either shows the length-adjusted estimates alongside the unadjusted ones (paired forest plot), or the panel is explicitly titled "unadjusted";
2. the abstract sentence be reordered so the abolition is not a "however" clause after the claim — e.g., state up front that unadjusted enrichment is universal but length-explained;
3. the residual EA/SCZ survival under CMH and alternative mappings be tested with a method that models gene length and LD jointly — MAGMA gene-set analysis (which conditions gene-level association on gene size and LD) or stratified-LDSC with an RD-linked annotation. Fisher tests on MAPPED_GENE tokenization are the weakest design available for exactly the bias the paper diagnoses. If MAGMA/S-LDSC is out of scope, the EA/SCZ residual should be presented as "unresolved", not "retained".

Relatedly, the pooled-OR comparison (neuro 3.00 vs. control 2.02; z = 3.49, P = 4.8 × 10⁻⁴) treats the 12 trait tests as independent. They are not: EA–cognitive performance rg ≈ 1, SCZ–bipolar rg ≈ 0.7, and the RD/GD sets are identical across traits. The fixed-effect pooling CI is therefore anti-conservative; the trait-label permutation (P = 0.009, 1,000 permutations) partially rescues it, but with 12 labels the granularity is coarse and the control family has only 4 members (Fisher 8/8 vs. 2/4, P = 0.091 — underpowered, as the authors note). The claim "strongest and most consistent" currently rests substantially on 2 of 4 controls happening to be null; with, say, bone density or eosinophil count as controls the pattern could look different. Control-trait selection should be justified a priori (why these four?), and ideally expanded.

### M2. Asymmetric length handling: GDS is score-level residualized on CDS length; RDS is never residualized on gene span

The framework residualizes the GDS substitution components on log₁₀(CDS length) — correctly removing length-inflated BUSTED significance at the *score* level. The RDS side receives no analogous treatment, even though the dominant RD component (non-coding conservation divergence, largest contributor for 52% of RD genes) correlates ρ = 0.47 with gene span, producing the 4.7× span skew that then propagates into *every* downstream validation: GWAS mappability (M1), ±50 kb HAR/hCONDEL proximity windows (the manuscript itself shows the hCONDEL enrichment is confined to the long-gene tertile, OR = 2.06), and the SynGO/GO enrichment that requires post-hoc length-matched rescue.

Post-hoc covariate adjustment inside each enrichment test is not equivalent to a length-comparable classification: the 7/11 SynGO survival is conditional on a class definition that is itself span-laden. A span-residualized RDS sensitivity classification (residualize the nc-divergence component — or the composite RDS — on log₁₀ span, reclassify, and re-run the LOO HAR/hCONDEL matrix and the SynGO/GO enrichment) would settle whether the RD class and its "length-robust synaptic core" survive design-level length control. If they do, the paper is substantially strengthened; if they do not, the specificity claim needs re-scoping. Either way this is the single most informative analysis currently missing, and it is cheap relative to what has already been done.

A secondary note: GDS residualization removes the *mean* length trend of BUSTED statistics but not the *power* differential (short genes have lower LRT power regardless of residualization), so even the coding side retains a length-dependent detection variance. This deserves one sentence in Limitations.

### M3. The 48-test one-sided BH family and the joint-LOO variant: four loose ends

(i) **One-sided canonical status is post hoc.** The manuscript candidly states there is no preregistration. One-sided Fisher tests in the enrichment direction are conventional, and I verified that the weakest canonical call (intelligence, one-sided q = 0.0432) is not an artifact of the 48-family choice — under a neuro-only 8×4 family the same p = 0.0090 would give q ≈ 0.036. Still, the family definition deserves a sentence noting that including the 16 mostly-null control tests *raises* BH critical values for the neuro tests, i.e., the unified family is *less* conservative for the headline claims than a neuro-only family, not more. As written ("headline and control traits together") a reader may assume the opposite.

(ii) **The joint LOO-brain-tau weight scheme is undocumented.** Methods and Supplementary Text S4 give the renormalized weights for the four single-component removals but not for the joint removal (presumably caMPRA 0.545 / nc 0.455 after renormalization). This is a reproducibility gap in the variant that now carries the expression-confounding rebuttal (6/8 traits, 7/74 SynGO). Please state the joint weights explicitly, and report the overlap between the joint-RD 403 genes and the full-model 293 — the mechanism of the enlargement (renormalization upweights the span-correlated nc-divergence component, which is computable for every gene) should be stated, since it predicts — and the span data confirm (joint-RD median 105.8 kb ≈ 4× unclassified) — that the joint variant is *more*, not less, length-driven. The joint variant answers the expression-circularity question but tightens the length confound; this trade-off is currently invisible.

(iii) **Numeric inconsistency in the MDD GD depletion.** Main text: "depletion for major depression … not under LOO reclassification (two-sided q = 0.11–0.15)". Table S5 gives LOO-brain 0.1426, LOO-tau 0.1478, and joint LOO-brain-tau 0.0811. The stated range matches neither the two single LOOs (0.14–0.15) nor all three variants (0.08–0.15). Correct the range and state which variants it spans.

(iv) **Table 3 omits the joint variant.** The LOO matrix reports full + four single removals for HAR/hCONDEL; the joint LOO-brain-tau classification (n = 403) is used for GWAS and SynGO but its HAR/hCONDEL enrichment is never shown. Add the row for completeness, or state why it is excluded.

---

## Minor comments

1. Methods, GWAS section: "trait sets span from 103 (regulation-driven ∩ educational attainment) to 822 universe genes" mixes two different quantities (an overlap count and a trait-set size) into one range; the smallest trait set is actually Crohn disease (95). Rewrite.
2. The length-correction power asymmetry is not discussed: logistic adjustment of term membership on class will retain large GO BP terms far more readily than small SynGO terms (≥3 universe genes). The 18/18 GO BP versus 7/11 SynGO survival contrast could partly reflect term size/power rather than biology; a brief power note (or per-term n) is warranted before calling the SynGO losses substantive.
3. Fig. 7c legend: "all eight traits are FDR-significant for regulation-driven under the unified 48-test one-sided family" — please mark in the figure itself that intelligence is marginal under the two-sided reference (q = 0.061), since the figure is the object most likely to be reused out of context.
4. The headline GD prevalence (1,214; 24.4%) inherits the --srv No detection inflation: the srv sensitivity is applied to the detection rate (34% → ~28%) but never propagated through the classification. A one-line bound on how the GD count/prevalence ratio would shift under the srv point estimate would close this.
5. Negative-control trait choice: height and T2D are enriched, Crohn and LDL are not. Note for readers that height GWAS genes are themselves enriched for long genes and for enhancer-linked biology, so height is arguably a *positive control for the length/mappability confound* rather than a pure negative control — this interpretation actually strengthens the M1 diagnosis and is worth stating.
6. The permutation floor (1/1,000) is stated for gene-set permutations; please also state it for the trait-label permutation (P = 0.009 is fine, but with 12 labels the number of distinct 8-vs-4 partitions is limited — report the permutation design, not just the seed).
7. Abstract, "parity would require ~59% assay coverage" — consider adding "(55–63%)"; the CI exists in the main text and the bare 59% overstates precision of a linear extrapolation that the Limitations themselves flag as an upper bound.
8. Supplementary Text S4 states the LOO anchoring family as "8 traits × 4 classes" in one sentence and "12-trait × 4-class" for controls in the next; align the wording so the 48-test unified family is described identically everywhere.

---

## Verification appendix (independent checks performed)

- **S5 CSV structure:** 192 data rows = 4 variants × 12 traits × 4 classes; header fields as described in Methods. ✓
- **Hand recomputation of ORs from counts (full variant):** SCZ/RD: (36×4460)/(257×221) = 2.827 vs. CSV 2.8269 ✓; Height/RD: 2.2146 vs. 2.2146 ✓; MDD/GD_strict: 0.5632 vs. 0.5632 ✓.
- **8/8 full-variant neuro enrichment (one-sided q < 0.05):** confirmed for all 8 traits; weakest = intelligence q = 0.0432 (two-sided q = 0.0610, matches text "0.061"). ✓
- **LOO-brain 6/8** (ASD q = 0.189, intelligence q = 0.060 fail) ✓; **LOO-tau 5/8** (ASD, intelligence, bipolar fail) ✓; **joint LOO-brain-tau 6/8, n = 403** (ASD, intelligence fail) ✓.
- **Negative controls:** height and T2D FDR-significant in RD under all four variants; Crohn and LDL never significant. ✓ (Note: joint-variant height OR = 2.46 slightly exceeds the main-text "1.9–2.3" range, which is accurate only for full/LOO-brain/LOO-tau.)
- **8/8 vs. 2/4 Fisher P = 0.091:** reproduced (hypergeometric one-sided = 45/495 = 0.0909). ✓
- **EA set size:** 822/4,974 = 16.53% ✓; EA overlap 103 ✓.
- **Inconsistency found:** MDD GD_strict LOO two-sided q range (see M3-iii): text 0.11–0.15 vs. CSV 0.081–0.148.
- **Table 1/2/4 counts** internally consistent with the classification summary in the abstract and Fig. 1c.

---

## Score justification

6.5/10. The methodological candor is exemplary and the core comparative-genomics machinery (LOO validation, coverage quantification, hCONDEL independence, disjointness analysis) is sound and well documented. The deduction reflects: a title/abstract-level claim that rests on an anchoring analysis the paper itself refutes after length adjustment (M1); a design-level length asymmetry between GDS and RDS that post-hoc corrections cannot fully substitute for (M2); and a set of fixable but real documentation/numeric gaps (M3). If M1 and M2 are addressed analytically — MAGMA or S-LDSC for the residual trait claims, and a span-residualized RDS sensitivity classification — this becomes a strong paper whose central message (prevalence vs. specificity, conditional on evidence coverage) is genuinely useful to the field.
