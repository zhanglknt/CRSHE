# Blind Verification Report — Paper A R6 Revision

**验证人**: general-purpose-30（独立统计验证员）
**日期**: 2026-09-18
**对象**: `results/paper/manuscript_english_v9.md`（含补充文本、主表 xlsx、cover letter、投稿包）
**方法**: 不信任任何中间报告，全部从原始数据独立重算：
- `results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv`（4,974 基因分类 + 分量分数）
- `results/phase9_hardening/hcondel_gene_mapping_fixed.csv`、`loo_variant_per_gene.csv`、`tier4_empirical_disjointness.csv`、`tier2_anchoring_master_table.csv`、`gwas_confound_calibration.json`、`bh_bounded_support_check.json`、`neutral_sim_pvalues.csv`、`stratified_rerun_pvalues.csv`、`topology_sensitivity.json`
- `results/phase4_hyphy/final_analysis_v3/busted_results_v3.csv`（生产 BUSTED p 值）
- **原始外源数据**: `data/gwas/gwas-catalog-download-associations-v1.0-full.tsv`（1,150,105 行全量重扫）、`data/downloads/gtex_v11/gencode.v47.basic.annotation.gtf.gz`（基因跨度重算）
- 复算脚本与全部输出: `scripts/_blind_verify_r6/`（outputs/ 下 A/B/B2/C/D 五份计算日志）

统计口径均按稿件 Methods 声明执行（单侧富集 Fisher、Wald log CI、Newcombe hybrid score CI、Mann–Whitney、CMH、逆方差合并 logOR、Li & Ji Meff、置换种子 20260917 精确复现）。

---

## 总体判定

**62 项检查：56 项 PASS，6 项 FAIL（均为轻微数字/口径问题，无一影响任何结论的方向或显著性）。**

核心三链（hCONDEL 修复链、GWAS 长度校准链、BUSTED 校准链）的**所有 headline 数字均从原始数据精确复现**（多数到小数点后 3 位）。6 项 FAIL 全部是范围端点、口径混用或旧文件残留的表述级问题，逐项列出如下，供编辑层决定是否修订。

---

## 链 1：hCONDEL（17 项：15 PASS / 2 FAIL）

| # | 声明（稿件位置） | 独立复算值 | 判定 |
|---|---|---|---|
| 1.1 | hCONDEL 基因 112、覆盖 2.3%（L86、L124） | 112/4,974 = 2.25% | **PASS** |
| 1.2 | RD 富集 OR=3.19, P=7.8e-5, CI 1.90–5.37（L76、Fig4f、Table 3） | k=18/293, OR=3.1940, p=7.809e-5, CI (1.901, 5.366) | **PASS**（精确） |
| 1.3 | GD（含 relaxed, n=1,266）双侧耗竭 OR=0.59, CI 0.36–0.97, P=0.037（Fig4f） | k=19/1,266, OR=0.5923, p=0.0371, CI (0.360, 0.974) | **PASS**（见注 a） |
| 1.4 | 5 个 LOO 变体 hCONDEL OR：3.19/2.13/4.60/4.67/2.56 + CI（Table 3） | 3.194/2.132/4.598/4.669/2.560；CI 全部吻合（如 LOO-nc 4.67 (2.368, 9.206)） | **PASS**（精确） |
| 1.5 | 摘要 "deletions (2.1–4.7)"；L74 hCONDEL p 范围 2.3e-5–5.3e-4 | LOO OR 范围 2.13–4.67；p 范围 2.27e-5–5.34e-4 | **PASS** |
| 1.6 | L74 HAR p 范围 "3.3 × 10⁻²⁰ to 4.7 × 10⁻⁵" | 五变体 HAR p 实际范围 **3.59 × 10⁻²⁵（LOO-nc）– 4.74 × 10⁻⁵**；3.3e-20 仅为 full 模型值，未覆盖 LOO-tau (3.42e-21) 与 LOO-nc (3.59e-25) | **FAIL**（轻微，方向保守：实际 p 更小，"全部通过 Holm" 的结论反而更强） |
| 1.7 | HAR∩hCONDEL 22、19.6%、22/476=4.6%、φ=0.052、χ² p=4.6e-4（L76、Fig4c） | both=22, 19.6%, 4.6%；φ=0.0520（未校正 χ²=13.43）；Yates 校正 χ² p=4.60e-4 | **PASS**（见注 b） |
| 1.8 | Fig 4c 三分片：HAR-only 72/454 OR=3.67 P=4.3e-16；双重 9/22 OR=11.38 P=1.9e-6；hCONDEL-only 9/90 OR=1.80 P=0.081 | 72/454 OR=3.6664 p=4.280e-16；9/22 OR=11.379 p=1.904e-6；9/90 OR=1.800 p=0.0814 | **PASS**（精确） |
| 1.9 | Tier 4 三行 complement 率 + Newcombe CI + Fisher p（L124） | HAR 36.55% vs 33.70%, +2.9pp (−1.6, +7.5), p=0.222；hCONDEL 30.36% vs 34.06%, −3.7pp (−11.6, +5.4), p=0.480；combined 35.51% vs 33.78%, +1.7pp (−2.4, +6.0), p=0.423 — 三行全部逐位吻合 | **PASS**（精确） |
| 1.10 | Table 3 HAR 列 4.15/1.61/6.39/9.39/2.78 + CI | 4.146/1.614/6.386/9.388/2.775；CI 全吻合 | **PASS**（精确） |
| 1.11 | Fig 7a RD 率 17.0–19.1%、P 3.8e-5–1.7e-25 | 19.12/16.96/17.67%；p 3.79e-5–1.67e-25 | **PASS** |
| 1.12 | L86：41/59=69.5%、基线 293/3,284=8.9%、caMPRA 97 (1.95%)、HAR 476 (9.6%) | 41/59、293/3,284、97 (1.95%)、476 (9.57%) 全部吻合 | **PASS** |
| 1.13 | L124 caMPRA 耗竭 8.2% vs 25.8%, OR=0.26, P=2.7e-5 | 8/97=8.25% vs 1,258/4,877=25.79%, OR=0.259, p=2.70e-5（GD 含 relaxed 口径） | **PASS** |
| 1.14 | Methods 双路线：70/112、Jaccard 0.625、name OR=2.90 P=0.0023、coord OR=3.58 P=2.1e-5；583→581 映射、2 失败 | 70/112、0.6250、2.897/0.00226、3.582/2.140e-5；liftover 581/583、失败 2 个 | **PASS**（精确） |
| 1.15 | L74 caMPRA LOO 检验 OR=1.84, CI 1.09–3.02, P=0.012, 校正 P=0.059 | 该数值仅在**旧 89 基因 caMPRA 集**下成立（LOO-caMPRA RD=800, k=23, K=89 → OR=1.8424, p=0.0117；与 phase8 `campra_regulatory_summary.json` 一致）；用**生产管线 97 基因集**复算同一检验：k=23/800, **OR=1.64, p=0.032**（Holm×5≈0.16） | **FAIL**（轻微：数字来自被取代的 phase8 口径；定性结论"校正后不显著"在两种口径下均成立，且新口径下更不显著） |
| 1.16 | Table 3 注 "all survive Holm (p<0.001)" | 最大原始 p=5.34e-4；Holm 校正后最大 5.34e-4 | **PASS** |
| 1.17 | Fig 4a "agreement 0.84–1.00"、Fig 4d "110–800 基因、OR>2.1" | LOO 一致性 0.836–0.928（1.00 为 full 自身）；RD 110–800、hCONDEL OR 最小 2.13 | **PASS**（口径说明） |

**注 a**：Fig 4f 图例写 "gene-driven genes (OR = 0.59...)"，该值对应 **GD 含 relaxed (19/1,266)**；strict GD (17/1,214) 为 OR=0.55, p=0.019。Table 3 注明 GD 不含 relaxed，两处 "gene-driven" 口径不同，建议图例注明。
**注 b**：φ=0.052 来自未校正 χ²，而 χ² p=4.6e-4 来自 Yates 校正 χ²（各自均可复现，属约定混用）。

---

## 链 2：GWAS 长度校准（14 项：12 PASS / 2 FAIL）

全部从原始 GWAS Catalog TSV（1,150,105 行）+ GENCODE v47 GTF 独立重建 12 个性状基因集（12/12 与主表逐一吻合）后重算。

| # | 声明 | 独立复算值 | 判定 |
|---|---|---|---|
| 2.1 | RD 中位长度 133.6 kb vs 28.6 kb、4.7×、Wilcoxon P=1.5e-71（L126、摘要、Discussion） | 133,590 vs 28,648, ratio 4.663, MWU p=1.453e-71 | **PASS**（精确；见注 c） |
| 2.2 | nc 分量 Spearman ρ=0.47，其他分量 \|ρ\|≤0.09 | rds_nc ρ=0.4719；doan 0.0861 / tau −0.0691 / brain 0.0648 | **PASS**（精确） |
| 2.3 | logistic 校正后 12 性状全部 q>0.44；EA 2.99→1.22、Height 2.22→1.06、SCZ 2.83→1.27 | 全部 12 个性状 OR_adj/p 逐一精确复现（EA 1.219, Height 1.063, SCZ 1.272）；min q=0.4479 | **PASS**（精确） |
| 2.4 | CMH：EA 1.87, P=3.9e-6；SCZ 1.95, P=8.2e-4 | 1.8653/3.886e-6；1.9486/8.154e-4（从三分量表独立重算 CMH 统计量） | **PASS**（精确） |
| 2.5 | 组间：pooled 3.00 (2.60–3.47) vs 2.02 (1.70–2.39)，z=3.49, P=4.8e-4；置换 P=0.009；校正后 1.23 vs 1.06, P=0.24 | 3.003/2.015, z=3.492, p=4.786e-4；置换（种子 20260917 精确复现）p=0.0090；1.234/1.063, p=0.238 | **PASS**（精确，含置换逐位置复现） |
| 2.6 | author-reported：EA OR=2.99, P=0.015；SCZ OR=1.81, P=0.007 | 从原始 TSV REPORTED GENE(S) 重建：EA 2.991/0.0151；SCZ 1.814/0.0072 | **PASS**（精确） |
| 2.7 | Jaccard 均值 0.148、Meff=3（Li & Ji 2005） | 0.1482；特征值 >1 截断和 3.20 → Li&Ji 装箱 = 3 | **PASS** |
| 2.8 | RD 8 神经性状 OR 1.9–4.0 全 FDR 显著；GD 0.56–1.03 均不显著（L126） | 神经 RD OR 1.947–4.005；神经 GD OR 0.563–1.033，q 全 >0.046 | **PASS**（精确） |
| 2.9 | EA 103/293；822 (16.5%)；trait set 103–822（L50、L126） | k=103/293；822/4,974=16.53%；103–822 | **PASS** |
| 2.10 | LOO：LOO-brain 6/8、LOO-tau 5/8 显著；GD 全 null；Height/T2D 三分类下均 FDR 显著 | 6/8（缺 INT、ASD）、5/8（缺 INT、ASD、BP）✓；GD 三分类全 null ✓；Height q=1.5e-9/7.5e-12/1.4e-5、T2D q=1.5e-5/3.8e-5/1.1e-3 均 <0.05 ✓ | **PASS**（显著性部分） |
| 2.11 | 同句 "height (OR = 2.2–2.3) 和 T2D (OR = 2.1–2.3) ... under all three classifications" | 主表实际：Height **2.215 / 1.878 / 2.335**；T2D **2.148 / 1.770 / 2.336** — LOO-brain 值（1.88、1.77）落在所述范围之外 | **FAIL**（轻微：OR 区间端点错；FDR 显著性声明本身成立） |
| 2.12 | "GWAS-mapped genes are longer than unmapped genes for every trait (median 1.6–3.5-fold)" | 12 性状实际中位比 **1.29（LDL）– 4.12（认知表现）**；方向（全部 >1）成立 | **FAIL**（轻微：范围端点错，应为 1.3–4.1） |
| 2.13 | EA-RD 主导分量 59% vs 48%、χ²=12.3, P=0.006、brain 1/103（L70、L126） | 61/103=59.2% vs 92/190=48.4%；4×2 表 χ²=12.33, df=3, p=0.0063；brain 1/103 | **PASS**（精确） |
| 2.14 | 摘要 "all adjusted q > 0.44"、"adjusted OR 1.02–1.64" | min q=0.4479 ✓；OR 范围 1.02–1.64 覆盖 11 个富集性状（Crohn 校正 OR=0.33，为耗竭方向，在范围外） | **PASS**（附口径说明） |

**注 c**："versus 28.6 kb for the remainder of the universe" 的对照集实为**同时剔除 GD（含 relaxed）**后的 3,415 基因（neutral+dual）；纯非 RD 集（4,681）中位为 28.8 kb。数值可复现，措辞略欠精确。

---

## 链 3：BUSTED 校准（11 项：10 PASS / 1 FAIL）

| # | 声明 | 独立复算值 | 判定 |
|---|---|---|---|
| 3.1 | p 恰 = 0.5 的基因 1,960、无 >0.5（Limitations 2） | 1,960 精确；0 个 >0.5；39.4% at 0.5、11.9% at 0 | **PASS**（精确） |
| 3.2 | 压缩校正 1,690→1,639（−3.0%） | BH(α=0.05)=1,690；p*=2p 后 BH=1,639（−3.02%） | **PASS**（精确） |
| 3.3 | 模拟零分布：45.8% at p=0.5、无 >0.5、1.8% p<0.05、经验 FDR 0/400 | 183/400=45.75%、0、7/400=1.75%、BH 拒绝 0/400 | **PASS**（精确） |
| 3.4 | 真实数据 p<0.05 占 35.8% | 1,782/4,974 = 35.83% | **PASS** |
| 3.5 | baseline 49.4% vs 生产 49.9%、r=0.97、98.5% 一致 | 194/393=49.36% vs 0.4987、r=0.9728、98.47%；差 0.509pp（四舍五入 0.5pp） | **PASS**（见注 d） |
| 3.6 | CpG 47.8%、−1.3pp、93.8% 保留 | 共同 391 基因：47.83% vs 49.10%（−1.28pp）、180/192=93.75% | **PASS**（见注 e） |
| 3.7 | 拓扑 46.6%、−1.3pp、92.1% 保留 | 共同 371 基因：46.63% vs 47.98%（−1.35pp）、164/178=92.13% | **PASS**（精确） |
| 3.8 | srv：42.6%→34.5%（289 可比基因）、−8.0pp、74.8% 保留、110/400 超时、289 基因 | 289 共同基因：42.56%→**34.60%**（100/289）、−7.96pp、92/123=74.80%、110 超时 ✓ | **FAIL**（轻微：34.5 应为 34.6；34.5% 是全部 290 个 srv 基因的率 100/290=34.48%，与 "on the 289 comparable genes" 的口径不符。−8.0pp、74.8%、110 超时均正确） |
| 3.9 | "relative 19% reduction that extrapolates to ~27% universe-wide"；"27–34%" | (42.56−34.60)/42.56=18.7%≈19%；34.0%×(1−0.187)=27.6%≈27 — 外推逻辑**自洽成立**（备选保留率逻辑 34.0%×0.748=25.4%，稿件未采用且已明确声明其外推逻辑） | **PASS** |
| 3.10 | Limitations (2) 全部数字 | −3.0%、45.8%、1.8%、0/400、35.8%、−1.3pp×2、−8.0pp、74.8%、110/400 均吻合（除 3.8 的 34.5 一处） | **PASS** |
| 3.11 | Cover letter：1,690→1,639（≤3%）、0% FDR、CpG/拓扑各 −1.3pp、27–34% | 全部与复算一致 | **PASS** |

**注 d**："within 0.5 percentage points" 实为 0.509pp，四舍五入为 0.5，属边界通过。
**注 e**：Discussion 括号 "(47.8% vs. 49.4%)" 混用分母（47.8% 为共同 391 基因口径，49.4% 为 baseline 全部 393 基因口径）；一致口径应为 47.8 vs 49.1 或 48.1 vs 49.4。−1.3pp 与 93.8% 本身正确。

---

## 链 4：文档一致性（6 项：5 PASS / 1 FAIL）

| # | 声明 | 复核结果 | 判定 |
|---|---|---|---|
| 4.1 | 摘要所有 GWAS/hCONDEL 数字与正文一致 | 1,690 (34.0%)、604 (12.2%)、52、1,214/293/11、HAR 1.6–9.4、hCONDEL 2.1–4.7、30.4–36.6% vs 33.7–34.1%、11 SynGO、OR 1.9–4.0、q>0.44、4.7×、2–10%、293–1,906、24.4% vs 5.9% — 全部与正文及复算一致 | **PASS** |
| 4.2 | Methods provenance：3,171/3,169/3,161/3,168 | `doan2024_campra/doan_summary.json`：total_hars_doan=3,171、liftover_success=3,169、har_coord_overlap=3,161、total_hars_ours=3,168 | **PASS** |
| 4.3 | 583/510、508/3,171=16.1%、97/89 | 583（hcondels_liftover_fix total）✓、581/2 失败 ✓、508 ✓、16.1% ✓（510 为 McLean 原文外部数字，无法本地核验，已注明出处）；**"89 of the 97 within the 4,974-gene universe" 与生产数据矛盾**：v7 分类表中 has_doan_campra=True 的 **97 个基因全部位于 4,974 universe 内且全部 HAR-proximate**（97/4,974=1.95%，正是 L86 与 Fig 1d 所用口径）；"89" 溯源自被取代的 phase8 旧汇总（campra_regulatory_summary.json，其 HAR 基因集为 364 个，亦与现行 476 不符） | **FAIL**（该从句与 L86/L104/tier3 及 v7 数据自相矛盾） |
| 4.4 | 新引文 5 条存在且正文引用 | Li & Ji 2005、Mendes & Hahn 2016、Scally & Durbin 2012、Scally et al. 2012、Vanderpool et al. 2020 — 引文列表与正文（Methods L50、Discussion L172）均存在 | **PASS** |
| 4.5 | 全文无残留旧值（183、0.0417、OR=2.94、30.6%、n=630、15.8%、"strongest and most consistent"、"Supplementary Table S2"→S5） | 对稿件/补充文本/cover letter 三份文件正则检索：**零命中**；GWAS 主表引用正确指向 Supplementary Table S5（L50）。唯一近似项：Limitations (13) 保留一处连字符用法 "strongest-and-most-consistent"（"we therefore frame the trait link as strongest-and-most-consistent rather than exclusive"），为刻意保留的定性表述，与被删除的旧断言不同 | **PASS** |
| 4.6 | Table 3 表注与数据一致性 | GD/RD 列 1,214/293、1,020/800、1,339/148、1,369/110、1,099/613 与 loo_variant_per_gene.csv 逐一吻合 | **PASS** |

---

## 链 5：图表/投稿包一致性（3 项：3 PASS）

| # | 声明 | 复核结果 | 判定 |
|---|---|---|---|
| 5.1 | Main_Tables_v9.xlsx（results/paper/ 与 MBE_submission_PaperA/04_tables/ 两处）Table 3 hCONDEL 列 = 3.19/2.13/4.60/4.67/2.56 + CI | 两处文件完全相同：3.19 (1.90–5.37)、2.13 (1.41–3.24)、4.60 (2.52–8.40)、4.67 (2.37–9.21)、2.56 (1.66–3.94) — 与独立复算逐位一致 | **PASS** |
| 5.2 | PaperA_Supplementary_Tables.xlsx 有 "S5 GWAS anchoring" sheet 且 108 数据行 | 两处副本均有该 sheet；表头后数据行 = **108**（12 性状 × 3 分类 × 3 类），首行 EA RD OR=2.827 与主表/复算一致 | **PASS** |
| 5.3 | MBE_submission_PaperA_single.pdf 存在且 ≥50 页 | 存在（项目根目录），**53 页** | **PASS** |

---

## FAIL 项汇总（6 项，均轻微）

1. **L74 HAR p 值范围**：写 "3.3 × 10⁻²⁰ to 4.7 × 10⁻⁵"；实际五变体范围为 3.6 × 10⁻²⁵（LOO-nc）– 4.7 × 10⁻⁵。下端点只取了 full 模型。hCONDEL 范围正确。
2. **Methods provenance "89 of the 97 within the 4,974-gene universe"**：生产数据中 97 个 caMPRA 基因全部在 universe 内（与 L86 "97 (1.95%)"、L104、Fig 1d、tier3 一致）；"89" 来自被取代的 phase8 汇总（其 HAR 集 364 个 vs 现行 476）。
3. **L126 Height/T2D OR 范围**："height (OR = 2.2–2.3) / T2D (OR = 2.1–2.3) under all three classifications"；主表实际 Height 2.215/1.878/2.335、T2D 2.148/1.770/2.336，LOO-brain 值（1.88、1.77）在范围外。三分类下 FDR 显著的声明本身成立。
4. **L126 性状长度比范围**："median 1.6–3.5-fold"；实际 1.29（LDL）–4.12（认知表现）。方向（每性状 >1）成立。
5. **Discussion srv 率 "34.5%"**：289 可比基因口径下应为 34.6%（100/289）；34.5% 是 290 个 srv 基因的全组率。−8.0pp/74.8%/110 超时均正确。
6. **L74 caMPRA LOO 检验 OR=1.84/P=0.012/校正 P=0.059**：仅在旧 89 基因 caMPRA 集下复现（k=23, K=89, RD=800 → OR=1.8424）；生产 97 基因集下同一检验为 OR=1.64, P=0.032（Holm ≈0.16）。定性结论（校正后不显著、不作独立验证）在两种口径下均成立。

**附带口径备注（不计 FAIL）**：Fig 4f "gene-driven OR=0.59" 实为 GD 含 relaxed（strict 为 0.55）；CpG 括号 "47.8% vs 49.4%" 混用分母；"remainder of the universe" 实际对照集同时剔除了 GD；φ（未校正 χ²）与 χ² p（Yates）约定混用；"adjusted OR 1.02–1.64" 不含 Crohn 的 0.33（耗竭方向）；baseline 复现差 0.509pp（四舍五入 0.5）。

## 结论

R6 修订的三条关键证据链——hCONDEL 坐标修复链、GWAS 基因长度混杂校准链、BUSTED 有界支撑/敏感性校准链——**全部从原始数据（含 115 万行 GWAS Catalog 与 GENCODE GTF 全量重扫）独立复现，headline 数字无一错误**。所发现的 6 处 FAIL 均为端点/口径/旧文件残留的表述级问题：3 处（#1、#3、#5）是范围端点或分母口径偏差，2 处（#2、#6）源自同一被取代的 phase8 caMPRA 旧分析，1 处（#4）是范围端点。全部不改变任何结论的方向、显著性或论证结构；其中 #1、#6 的修正方向反而强化稿件声明。建议在下一轮文字修订中一并核销上述 6 处及附带口径备注。
