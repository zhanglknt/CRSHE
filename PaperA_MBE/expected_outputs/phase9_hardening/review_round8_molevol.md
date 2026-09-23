# Review — Round 8 (Molecular Evolution Methods)

**Manuscript:** *Coding Selection Dominates by Prevalence, Regulatory Selection by Specificity: A Comparability-Corrected 10-Primate Analysis* (v9)
**Reviewer expertise:** molecular evolution methodology (dN/dS models, HyPhy BUSTED/RELAX/GARD/MEME, alignment quality, CpG/gBGC confounding, multiple testing)
**Materials reviewed:** main text, Supplementary Text v9 (S1–S6), Main Tables v9, figure legends (Fig. 1–7), and the machine-readable result files referenced in the text (`srv_primary_estimate.json`, `alignment_rerun_sensitivity.json`, `neutral_sim_calibration.json`, `neutral_sim_parametric.json`, `gard_troubleshoot.json`, `topology_sensitivity.json`, `stratified_rerun_sensitivity.json`, `cpg_masking_stats.summary.json`).

---

## Overall assessment

This is an unusually self-critical BUSTED/RELAX scan of 4,974 one-to-one primate orthologs with a comparability-corrected gene-driven vs. regulation-driven classification. From the molecular-evolution side, the paper does many things right that most scans omit: a stratified 400-gene sensitivity battery with a baseline replication (0.51 pp, r = 0.97 — internally consistent), CpG masking (−1.3 pp), an alternative topology (−1.3 pp), a *completed* --srv Yes re-run (400/400, removing the earlier timeout bias), two independent null calibrations (evolver fixed-parameter and per-gene parametric bootstrap), honest disclosure of the RELAX run-list confounding, the K = 0 boundary artifact, the Selectome identifier bug with sensitivity re-classification, and the bounded-support/chi-bar-square anatomy of the BUSTED p-value distribution. The Limitations section (16 items) is among the most candid I have seen.

That said, the two new quantitative centerpieces of this revision — the 28.0% srv-adjusted point estimate and the +14.8 pp codon-aware realignment — are both presented more confidently than the underlying data allow, and the recombination (GARD) thread remains unresolved in a way that interacts with the paper's own ILS discussion. The problems are addressable; none requires redoing the production scan. I recommend **Major Revision**.

**Score: 6/10 — Major Revision.**

---

## Major issues (in order of severity)

### M1. The codon-aware realignment (+14.8 pp) is over-interpreted, and the interpretation is asymmetric

The text (Discussion, "Interpreting the 34% rate"; Limitations item 2) concludes that because protein-guided realignment with gap-column filtering *raises* detection from 49.4% to 64.1%, "alignment noise attenuated power rather than inflating false positives." The underlying data do not support this directional reading:

- **The two pipelines barely agree at the gene level.** Per `alignment_rerun_sensitivity.json`: Jaccard between significant sets = 0.538; 96 genes gained vs. 38 lost; 19.6% of baseline-significant genes lost; and — most telling — the Spearman correlation of p-values between the two pipelines on the same 393 genes is only **0.385**, versus 0.97 for the baseline replication. If realignment merely "removed noise" from the same underlying signal, p-values should remain highly correlated. A correlation of 0.39 means the realigned pipeline is measuring something substantially different, not the same thing more precisely.
- **59.1% of alignment columns were removed on average (median 61.5%, p95 81.5%).** After filtering, BUSTED is testing a different, much smaller, and non-random site set (median 685.5 retained columns). Gap-rich columns are exactly where alignment ambiguity, frameshifts, mis-assembly, and paralogous mapping concentrate; deleting them and retaining cleaner, more variable blocks changes both the signal and the effective null. A detection-rate increase after removing 60% of (mostly conserved or ambiguous) columns is equally consistent with (i) power recovery, (ii) enrichment of the retained set for high-substitution blocks including gBGC tracts and residual misalignment, and (iii) simple loss of the conserved sites that dilute the LRT. The data presented cannot distinguish these.
- **The realignment QC itself flags the production alignments as artifact-rich:** 68,680 in-frame stop codons masked across 378/400 genes, and 254 genes with non-codon-multiple sequence lengths. This is a startling rate of frameshift-level content in production alignments. It cuts both ways: it could mean the *production* 34% is partly misalignment-driven (high dN/dS from misaligned blocks), which is precisely the false-positive direction the authors claim to have excluded. At minimum this observation belongs in the main text, not only in a JSON.
- **Asymmetric adoption of perturbation results.** When a perturbation lowers the rate (srv, −8.7 pp) it is treated as correcting an upward bias and used to revise the headline estimate; when a perturbation raises the rate (+14.8 pp) it is treated as revealing attenuated power. Whichever direction supports "the signal is real" is adopted. The defensible joint conclusion is instead: **the detection rate is sensitive to alignment pipeline at the ±15 pp level and gene-level calls are fragile (Jaccard ≈ 0.54), so any point estimate carries large methodological uncertainty.**
- The newly significant 96 genes are not further characterized (RELAX K? Selectome? conservation?). And the realigned set was not itself run under --srv Yes, so the two largest perturbations were never crossed.

*Required:* (i) soften the causal claim to "consistent with power attenuation but not distinguishable from site-set enrichment"; (ii) report the Jaccard and p-correlation (0.39) in the main text alongside the +14.8 pp; (iii) report the stop-codon/frameshift QC numbers; (iv) ideally run the srv × realignment cross (or state why not); (v) acknowledge that production used nucleotide-mode MAFFT for coding sequence, and justify.

### M2. The 28.0% srv-extrapolated point estimate uses a biased estimator for a stratified sample and is reported without uncertainty

The extrapolation multiplies the universe rate (34.0%) by the *relative* rate change (40.7/49.4 = 0.825) measured on the stratified 400-gene sample. But the sample is 50% FDR-significant by design, while the universe is 34%, and the srv effect is strongly asymmetric across strata (per `srv_primary_estimate.json`): among baseline-significant genes, 44/194 (22.7%) flip to non-significant; among baseline-non-significant genes, 10/199 (5.0%) flip *to* significant. A stratification-weighted estimator — applying each stratum's flip rate to the corresponding universe stratum — gives 34.0 × (1 − 0.227) + 66.0 × 0.050 ≈ **29.6%**, not 28.0%. The proportional-ratio estimator is only correct if the relative srv effect is homogeneous across strata, which the authors' own flip counts contradict. The 1.6 pp difference is not large, but the manuscript presents "≈28.0% (≈1,394 genes)" as *the* point estimate in both Results and Discussion with no confidence interval and no statement of the estimator assumption. Sampling error alone (n = 400, ~54 status changes) puts several pp of uncertainty on this number; a rough binomial-level CI would be on the order of ±3 pp.

More importantly, **the srv sensitivity is never propagated into the classification that constitutes the paper's headline.** GDS requires BUSTED FDR significance; the flip rates above imply that, if the srv result generalizes, roughly one in five gene-driven calls (and a similar fraction of the regulation-driven class boundary, which requires BUSTED *non*-significance) is srv-sensitive. The prevalence claim "24.4% vs. 5.9%" — the paper's title claim — is built entirely on srv = No. At minimum, re-run the classification on the 400-gene subset under srv = Yes and report class-level stability; or state explicitly that the prevalence numbers are srv = No quantities with an estimated ~17% relative deflation of the coding side.

*Required:* report the stratified estimator (or both), a CI, and the classification-level srv sensitivity.

### M3. Recombination/ILS remains an unbounded confound; the dismissal argument is partly self-contradictory, and the decisive control was not run

The GARD narrative has improved (the standard-mode non-convergence is convincingly diagnosed as a configuration pathology — exhaustive candidate-breakpoint search on ~100%-variable, gap-rich alignments vs. a 600 s timeout; `gard_troubleshoot.json` reproduction tests are solid). But the substantive conclusion is weak:

- All 50/50 top genes return "breakpoints" with a median Δc-AIC of 1,255 — enormous model-improvement values — which are then declared alignment artifacts because they co-localize with alignment blocks and show pathological partition branch lengths (>100 subst/site). That diagnosis is plausible, but the obvious control exists and was not run: **GARD on the codon-aware cleaned alignments** (which the authors themselves produced for M1). If breakpoints vanish after block cleaning, the artifact interpretation is confirmed; if they persist on clean alignments, recombination or gene-tree discordance is implicated. As it stands, "non-informative" is asserted, not demonstrated.
- Supplementary Text S2 argues the primary inferences are uncompromised because "(2) the 10-species primate phylogeny is well established and concordant across independent studies." This is a non-sequitur that directly contradicts the manuscript's own Discussion, which (correctly) notes ~30% of the genome carries discordant genealogies in the hominine radiation and invokes SPILS as a dN/dS-inflating mechanism. Species-tree concordance says nothing about within-gene topology heterogeneity. Reason (2) should be deleted or rewritten.
- The single alternative-topology re-run (human sister to gorilla) tests one alternative *species* tree; it does not address gene-tree discordance, which is the mechanism the authors themselves raise. The current framing ("bounding the practical impact") overstates what that analysis bounds.

### M4. "BH validity is exact" is contradicted by the authors' own simulation, and FPR/FDR terminology is misused

Limitations item 2 states that because the boundary atom is the expected chi-bar-square signature "leaving p-values uniform on (0, 0.5), ... the validity of the Benjamini–Hochberg procedure is exact rather than approximate." But `neutral_sim_calibration.json` reports a KS test of the simulated null p-values against U(0, 0.5) with **p = 3.6 × 10⁻⁷ — uniformity is rejected** (in the conservative direction: mean p below 0.5 is 0.31 vs. 0.25; nominal-5% FPR is 1.8%). The empirical evidence therefore supports "conservative," not "exact" — and the manuscript says "conservative" in one paragraph and "exact" in another. Additionally:

- The chi-bar-square mixture is an asymptotic approximation; combined with the KS rejection, "exact" should be replaced with "valid and mildly conservative in simulation."
- "Empirical BH false-discovery rate ... 0 of 400" and "empirical FDR of 0.26% (1/388)" are misnomers: on all-null data with 0–1 total rejections, what is estimated is a **per-test false-positive rate** (with an upper 95% bound around ~0.75–0.9% at n ≈ 400), not an FDR among the 1,690 production rejections. The direction of the conclusion (the 34% is not a bounded-support artifact) is fine; the terminology and implied precision are not.
- The parametric bootstrap's 12/400 failed runs (388 denominator) are not mentioned anywhere in the text.
- One genuine limitation of the parametric bootstrap deserves a sentence: branch lengths fitted under the null partially absorb real human-branch signal into the simulated null (direction: inflates simulated FPR, i.e., conservative for the FDR claim — but it should be stated). The disclosed limitations (no indels, single-ω, no SRV) are otherwise appropriately candid.

### M5. MEME/site-level and ω-distribution framing need tightening

- Fig. 2e reports a median gene-wide ML ω of 1.96 with 62.1% of genes above 1. For a 10-primate tree with a short human foreground, a majority of genes having gene-wide ω > 1 will surprise every MBE reader and prima facie suggests alignment/branch-length pathology rather than biology. The parenthetical caveat (gene-wide ω averages across sites) does not explain why the *median* is ~2. This needs a real explanation (e.g., foreground-branch weight in the gene-wide estimate, alignment artifact contribution per M1's stop-codon findings) or readers will discount the entire scan.
- MEME completed only 8/20 genes, selected by gene-wide significance, with long alignments timing out — the 7/8 site-support rate is correctly disclaimed, but given M1 (production alignments contain frameshift-level content), MEME on stop-masked alignments would be more convincing; at minimum note that the MEME alignments were the production ones.

---

## Minor issues

1. Abstract/Results report 604/4,974 = "12.1%"; Discussion and Concluding remarks say "12.2%". 604/4,974 = 12.14% — pick one.
2. Table 2 lists current unclassified as 68.5%; Table 1 and the correct rounding of 3,404/4,974 = 68.4%. 
3. The four sensitivity re-runs use different comparison bases (393 common for srv, 391 for CpG, 371 for topology, 377/397 group results; topology lost 23/400 to timeouts, never mentioned in text). A small summary table of n per run and comparison basis would prevent reader confusion; the topology timeout attrition (5.8%) should be disclosed like the srv one was.
4. Production used HyPhy 2.5.100, re-runs 2.5.101 (disclosed in JSON only). State it in Supplementary Text S1.
5. The Storey π₀ analysis "adjusted for the bounded distribution" is described only by its outcome (π₀ = 1.0, same 1,690 genes). One or two sentences on the estimator (which λ / how the bounded support enters) are needed for reproducibility.
6. The 1,866-gene BUSTED attrition is acknowledged as length/convergence-related, and the timeouts were scaled by alignment length — so the universe (and hence the 34% rate, and every downstream class) is conditional on tractable alignments. One sentence acknowledging this selection on the universe itself is warranted, especially since longer alignments time out more (cf. MEME, GARD).
7. Limitation 8's strict-subset robustness paragraph is good; please also report the permissive-call fraction among the 96 realignment-gained genes (M1) — if enriched, it further undermines the "power not artifacts" reading.
8. Limitation 12 discloses linearity/class-migration caveats for the coverage extrapolation — good. The same candor is missing for the srv extrapolation (M2).
9. Fig. 4f vs. Table 3 use different gene-driven denominators (with/without relaxed class); both are defensible but the reader has to discover the difference — cross-reference them.
10. The evolver null used κ = 2 while fitted κ median is 3.02; the parametric bootstrap supersedes it, but say explicitly that the evolver simulation is retained as the fixed-parameter complement (the current text implies they are two equivalent calibrations).

---

## Specific revision suggestions (actionable)

1. Rewrite the realignment paragraph: report Jaccard 0.54, Spearman 0.39, 59% column removal, and stop-codon QC; replace the causal claim with a bounded statement; run or decline-with-justification the srv × realignment cross.
2. Replace the srv point estimate with a range: report both the ratio extrapolation (28.0%) and the stratification-weighted estimator (~29.6%), with a binomial/bootstrap CI; add the classification-stability re-run under srv = Yes on the 400-gene subset.
3. Run GARD (Faster mode) on the codon-aware cleaned alignments for the same 50 genes, or explicitly justify why block co-localization is sufficient evidence; delete S2 reason (2).
4. Replace "exact" with "valid and mildly conservative"; relabel "empirical FDR" as per-test FPR with an upper bound; disclose the 388/400 denominator; add the branch-length-absorption caveat.
5. Add the minor-issue fixes (12.1 vs 12.2, 68.4 vs 68.5, comparison-base table, HyPhy version, Storey details, universe-attrition sentence).

## Internal-consistency spot checks (all passed unless noted)

- 1,690/4,974 = 34.0%; class counts sum to 4,974; 1,214/293/52/11 percentages check.
- RELAX chain 1,660 → 728 → 124 → 117 (−7 degenerate K = 0) → 52 + 65 consistent; 1,610 active + 50 default = 1,660.
- srv: 49.36% → 40.71% on 393 common, −17.5% relative, ×34.0% = 28.0%, ≈1,394 — all match `srv_primary_estimate.json` (estimator caveat per M2).
- Realignment: +14.76 pp, 80.4% retention, 59.1% mean column removal — match JSON.
- Baseline replication 0.51 pp, r = 0.97, 98.5% — match JSON.
- Topology −1.35 pp (46.6%), 92.1% retained — match JSON.
- Null calibrations: evolver 45.8% point mass / 1.8% nominal / 0 rejections; parametric 73.7% / 1.29% / 1 of 388 — match JSONs.
- HAR 81/293 vs. 395/4,681 → OR 4.15; combined subset 476 + 112 − 22 = 566; GD hCONDEL OR 0.59 (0.36–0.97) — consistent.
- Bounded-support stress test 1,690 → 1,639 (−3.0%) internally consistent with BH thresholds ≤ 0.05 ≪ 0.5.

---

**Recommendation: Major Revision (6/10).** The sensitivity infrastructure is genuinely strong and the candor is exemplary; what is needed is not more computation so much as more symmetric interpretation of the computation already done, two bounded-estimator fixes, and one additional GARD control run.
