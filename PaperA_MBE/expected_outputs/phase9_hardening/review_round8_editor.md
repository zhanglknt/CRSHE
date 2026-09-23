# Round 8 — Associate Editor Assessment (MBE)

**Manuscript:** "Coding Selection Dominates by Prevalence, Regulatory Selection by Specificity: A Comparability-Corrected 10-Primate Analysis"
**Files reviewed:** manuscript_english_v9.md (full), Supplementary_Text_v9.md, Main_Tables_v9.md, cover_letter_v1.md, Figures 1–7 (figures_v2 PNG)
**Role:** Associate Editor, first-read evaluation. This is a fresh read; no prior rounds seen.
**Date:** 2026-09-23

---

## 1. Overall recommendation

**Score: 8/10. Decision: Minor Revision (send out / accept-in-principle trajectory).**

This is a mature, unusually self-critical submission. The two title claims are each supported by a distinct, internally validated evidence chain, and the manuscript has been engineered to survive exactly the objections a methods-literate MBE referee would raise. I would send this out with confidence; if reviews mirror the manuscript's own candor, I expect a path to acceptance after one revision. The issues below are editorial/narrative-level, not fatal, and none requires new computation — except that several require firm textual resolution before the paper is read by hostile referees.

## 2. Are the two title pillars supported?

**Pillar 1 — "Coding Selection Dominates by Prevalence":** Supported, with the qualifier honestly built into the claim itself. 1,214 GD (24.4%) vs. 293 RD (5.9%), and the coverage-sensitivity analysis (crossing at ~59% of the universe; HAR-restricted ceiling 410) correctly converts the prevalence statement from a naive class-size comparison into a conditional one ("dominates *unless* regulatory assays reach ~59% coverage"). The title's "by Prevalence" framing is therefore accurate and defensible. The residual vulnerability is that "prevalence" of coding selection itself rests on BUSTED FDR hits whose 34.0% rate the authors themselves call an upper bound; the --srv Yes point estimate (28.0%) and the RELAX-confirmed rate (12.1%) are both still ≫ 5.9%, so the qualitative claim survives every recalibration on offer. Good.

**Pillar 2 — "Regulatory Selection by Specificity":** Supported, but only after the manuscript's own corrections are read carefully. The strong version (8/8 neuropsychiatric/cognitive GWAS traits enriched in RD, 0/8 in GD) is length-fragile — the authors say so themselves (all adjusted q > 0.44). The durable specificity claim is the synaptic gene-set enrichment (7/11 SynGO + 18/18 GO BP surviving length correction) plus the hCONDEL primary validation. An editor reading the title against the abstract might expect the GWAS anchoring to be the "specificity" evidence and will then discover it is largely a long-gene effect. The abstract does disclose this ("length adjustment abolishes the trait enrichment"), so this is honest — but it sits one clause away from the 8/8 claim, and a careless reader (or a press office) will quote the 8/8. This is my main narrative concern (see §4, Q1).

## 3. Narrative honesty assessment

The limitations section is exemplary in coverage: RELAX run-list confounding, BUSTED upper-bound status, bounded-support BH calibration, neutral simulation, parametric bootstrap, gBGC, enhancer turnover, ILS/SPILS, MEME/GARD non-convergence, weight-choice non-preregistration. Disclosure is full without being self-immolating — each limitation is paired with a quantified bound on its impact. Two exceptions where disclosure tilts toward self-defense:

1. **The GWAS trait anchoring is presented as a headline result (Fig. 7c, "8/8") while its length-sensitivity is disclosed in prose.** The figure itself carries no length-correction caveat in panel (c); a reader skimming figures gets the strong version. Panel (b) does carry the "GD 0 FDR terms" comparison. Consider a one-line caveat inside the Fig. 7c legend ("trait enrichment is largely a long-gene effect; see text") — the legend is currently silent on this.
2. **Abstract sentence structure on HAR/hCONDEL.** "…remain enriched for human accelerated regions (odds ratio 1.6–9.4) and human-specific conserved deletions (2.1–4.7) under every classification variant (the primary independent validation)" — the parenthetical grammatically attaches to both evidence types, but the text is emphatic that only hCONDEL is the *primary independent* validation and that the HAR enrichment is partially nested with the caMPRA component. A referee will catch this and read it as overselling in the one place overselling is most visible.

## 4. Numeric consistency spot-check (abstract / text / tables / figures / cover letter)

Checked: 1,690/4,974 = 34.0% ✓; 1,214 (24.4%) ✓; 293 (5.9%) ✓; 52 (1.0%) ✓; 11 dual ✓; 3,404 (68.4%) ✓; srv Yes 28.0% ≈ 1,394 ✓ (40.7/49.4 × 34.0 = 28.0) ✓; crossing 2,914 (59%) ✓ vs. Fig. 7d ✓; equal-size point 2,158 (43%) ✓; HAR OR range 1.61–9.39 matches Table 3 ✓; hCONDEL 2.13–4.67 ✓; LOO matrix counts (1,020/800, 1,339/148, 1,369/110, 1,099/613) match Table 3 and Fig. 1e ✓; RELAX chain 1,701 → 1,660 → 728 → 124 → 117 internally consistent (728−124=604) ✓; 1,960/4,974 = 39.4% (p=0.5 atom) ✓; 592/4,974 = 11.9% (p=0) ✓; combined subset 476+112−22 = 566 ✓; counterfactual 1,906 vs. residual GD 543 consistent ✓; CpG/topology −1.3pp each consistent between text and cover letter ✓.

**Discrepancies found (all minor, all fixable in minutes):**

1. **12.1% vs. 12.2%** — 604/4,974 = 12.14%. The Abstract says 12.1%; the Discussion and Concluding Remarks say "approximately 12.2%". Pick one (12.1% is the correct rounding). This is the only flat inconsistency I found, and it sits in a headline sentence about the "true positive-selection rate".
2. **"+14.8 points" vs. displayed 49.4% → 64.1%** — the displayed values differ by 14.7 points. Presumably underlying values (≈49.36% → 64.12%) justify 14.8, but a referee doing mental arithmetic on the printed numbers will flag it. Show one more decimal or adjust the delta.
3. **Table 2 unclassified 68.5% vs. Table 1 68.4%** — 3,404/4,974 = 68.44%. Table 1 is right.
4. **Table 4 baseline row GD = 1,216 vs. 1,214 elsewhere** — disclosed in the table footnote (independent reimplementation, 0.04%), so acceptable, but expect a referee query anyway; consider bolding the footnote reference in the table title.
5. **"Discoveries" format mismatch (editorial, not numeric).** The cover letter targets MBE's *Discoveries* category. Main text is ~10,100 words with 7 figures and 4 tables — far beyond any short-format envelope and long even for a full MBE research article. Either retarget as a full Research Article (my recommendation; the supplementary apparatus supports it) or cut by half. Submitting this length under Discoveries invites a desk return.

## 5. MBE fit

Strong. The paper is about a genuinely MBE question (King–Wilson coding-vs-regulatory debate) answered with MBE-native methods (BUSTED/RELAX/MEME on a primate phylogeny), and its durable contribution is methodological practice (LOO circularity control, coverage auditing) of direct interest to anyone running comparative selection scans — exactly the MBE readership. The suggested reviewers (Kosakovsky Pond, Scally, Fraser, Wertheim) are appropriate and cover methods, primate genomics, and regulatory evolution. The cover letter is effective: it leads with the validation rather than the class sizes, and pre-emptively frames the 34% upper-bound issue — a sophisticated move that disarms the most likely hostile first read.

## 6. Submission-readiness gaps (must fix before anything leaves the building)

1. **Funding section is a placeholder** ("[Funding sources and grant numbers to be added.]").
2. **AI-use disclosure placeholder** left in Author contributions ("[If AI-based tools were used…]"). MBE/OUP requires an explicit statement either way.
3. **Zenodo DOI not yet assigned** — acceptable at submission but must be resolved by revision; confirm the GitHub repo is public at submission.
4. Word count / article-type decision (§4.5 above).

## 7. Minor editorial notes

- The phrase "regulation-driven genes are enriched for synaptic gene sets (11 SynGO terms at FDR < 0.05)" in the abstract, followed later by "7 of the 11 … survive explicit gene-length correction", is fine, but the abstract would be stronger if the length-robust subset (7/11, 18/18) were the lead number — that is the durable result.
- Discussion's srv/codon-aware passage stacks five rates (49.4, 40.7, 28.0, 64.1, 34.0) in one paragraph. Consider a small table or a figure panel; this is the densest numeric passage in the paper and the one most likely to generate confused referee arithmetic (cf. §4.2).
- Fig. 2a annotation "max p with FDR<0.05 = 0.02" is informative; keep.
- The hCONDEL-only subset attenuation (OR = 1.80, P = 0.081) is handled honestly with a power calculation; good practice.
- Table 3's note on span tertiles (hCONDEL enrichment confined to the long-gene stratum, OR = 2.06) is the right control and should be cited where the GWAS length correction is discussed, since it shows the hCONDEL validation behaves differently from the GWAS anchoring under length stratification — a genuine strength worth one sentence in the Discussion.

## 8. Top 3 issues for the authors

1. **Reconcile the "specificity" pillar's presentation with its own length correction.** The 8/8 GWAS anchoring (unadjusted) currently sits in the abstract, Fig. 7c, and the concluding remarks as co-equal with the length-robust SynGO signal, but the authors' own analysis shows the trait enrichment is largely a long-gene effect (all adjusted q > 0.44). The fix is presentational, not analytical: lead everywhere with the length-robust gene-set result, demote the 8/8 to "unadjusted, externally defined but length-sensitive", and add a length caveat to the Fig. 7c legend.
2. **Fix the flat numeric inconsistencies** (12.1% vs. 12.2%; +14.8 vs. displayed 14.7pp; Table 2 68.5% vs. 68.4%) and clarify the abstract parenthetical that currently implies both HAR and hCONDEL are "the primary independent validation".
3. **Resolve submission mechanics before review:** article type (this is a full Research Article at ~10,100 words, not a Discoveries piece), funding statement, AI-use disclosure, Zenodo DOI.

## 9. Editor's summary

After eight rounds of simulated scrutiny this reads like a paper that has already survived real review. The science is conditionally framed where it must be, the validation hierarchy (hCONDEL primary, HAR partially nested, caMPRA LOO not counted) is stated with unusual precision, and the remaining problems are of the kind an editor can see being fixed in one revision cycle: presentation hierarchy of the specificity claim, three numeric nits, and submission mechanics. **8/10, Minor Revision.**
