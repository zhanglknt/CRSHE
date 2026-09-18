# R6 Popgen Re-Check — 复核报告（Reviewer: Population & Quantitative Genetics）

**复核对象**：response_round6.md 中对 review_round6_popgen.md 各项意见的处置，以及修订后 manuscript_english_v9.md / Supplementary_Text_v9.md / gwas_confound_calibration.json / tier2_anchoring_master_table.csv（Supplementary Table S5）。

**复核结论**：P0-3 与 P0-4 属实质修复且处置诚实；P1-2 / P1-3 已修复；P1-1 大体修复但存在选择性报告的残留问题；另发现一处正文数字与数据不符。

---

## 逐条判定

### P0-3（GWAS 基因长度/映射混杂）— **CLOSED（附 1 项 P1 残留，见下）**

数据核验（gwas_confound_calibration.json）：

1. **混杂存在性确认**：RD 基因中位 span 133,590 bp vs 其余 28,6648 bp（4.66×，Wilcoxon p=1.5e-71）——与正文 "133.6 kb versus 28.6 kb (4.7-fold; Wilcoxon P = 1.5 × 10⁻⁷¹)" 一致。分量归因正确：nc 分量与 span 的 Spearman ρ=0.472，其余分量 |ρ|≤0.086（caMPRA 0.086、tau −0.069、brain 0.065）——与正文 "ρ = 0.47 ... versus |ρ| ≤ 0.09" 一致。我 R6 评审中预测的机制（nc 分量与基因跨度相关 → 长基因更易入 RD 类 + 更易被位置映射）被数据完全证实。
2. **Logistic 模型设定正确**：RD 状态 ~ trait 基因状态 + log10(span)（GENCODE v47 TSS–TES，4,974 基因全覆盖、0 缺失），每性状分别拟合，len_logOR≈2.19–2.23（长度项本身强显著）。12 性状调整后 OR 1.02–1.64、全部 q>0.44（BIP 最高 1.64、p=0.075 未达显著）——与正文 "adjusted OR 1.02–1.64, all q > 0.44; educational attainment 2.99 → 1.22, height 2.22 → 1.06, schizophrenia 2.83 → 1.27" 逐数吻合。正文如实报告"abolishes the enrichment for all 12 traits"，未做任何粉饰。
3. **CMH 残余混杂声明诚实**：正文明确限定 CMH "removes only between-stratum confounding"（仅去除层间混杂、不处理层内残余），这是准确且克制的表述；EA（common OR 1.87, p=3.9e-6）与 SCZ（1.95, p=8.2e-4）的数据与 JSON 的 P0_3c_tertiles.cmh 完全一致。分层数据本身也显示混杂并未在层内消失：EA 在第 3 长度层内 OR=1.99（p=6.7e-6），但该层富集未做长度内进一步校正——正文用"Coarser analyses bound the effect"的措辞定位这两项结果（作为界而非净效应），定位正确。
4. **文本落实**：Abstract 重写句（"4.7× the genomic median length, and gene-length correction abolishes the trait enrichment (all adjusted q > 0.44), leaving the synaptic gene-set enrichment as the length-robust functional signal"）、Results "Third, gene length is a structural confounder" 整段、Discussion "A conditional answer"（"largely a long-gene effect"）、Limitations (15)、Concluding remarks（"gene-length-sensitive"）——五处全部到位，口径一致。
5. **Supplementary Table S5 已打包**：PaperA_Supplementary_Tables.xlsx 含 "S5 GWAS anchoring" sheet（111 行 = 表头 + 108 数据行，与 Methods 声明的 12 性状 × 3 分类 × 3 类 = 108 行一致）；抽查 CSV 行（SCZ RD：n=257/k=36/OR=2.827/q=6.7e-6 等）与 JSON 及正文吻合。

**残留问题（升级为新的 P1，选择性报告）**：
(a) CMH 粗分层下 **height（common OR 1.33, p=0.029）与 T2D（1.42, p=0.028）同样保持名义显著**（JSON P0_3c_tertiles.cmh），且对 CMH 报告的 5 个性状做 BH 后 q≈0.035 亦存活——正文只写 EA/SCZ 存活，未提阴性对照同样存活。同理 locus-dedup 约定下 height 反而更显著（OR=3.06, p=1.2e-4），response_round6.md 中"intelligence 与 height 失去显著"仅对 author-reported 映射成立（height author-reported p=0.16），对 locus-dedup 不成立。当前写法给读者留下"粗分层下仅 EA/SCZ 独存"的偏乐观印象。要求：在 Results 该句补充 height/T2D 在 CMH 与 locus-dedup 下的表现，或改写为"educational attainment and schizophrenia survive across all three coarser conventions, whereas the negative controls survive only some"（若属实）并附完整数字。

### P0-4（"strongest and most consistent" 声明）— **CLOSED**

1. **正式组间检验已补**：未校正神经族 pooled OR 3.00（2.60–3.47）vs 对照族 2.02（1.70–2.39），z=3.49、p=4.8e-4，性状标签置换 p=0.009（1,000 次）——与 JSON P0_4a/P0_4b 一致；置换经验 p（0.009）与解析 p（4.8e-4）同向且均显著，检验设定（log-OR 加权合并、组间 z）合理。
2. **一致性比较如实报告**：8/8 vs 2/4，Fisher p=0.091，正文如实写"directional consistency"差异未达显著；LOO 变体下 6/8、5/8 亦保留原文。旧措辞 "strongest and most consistent" 已全文清除，替换为"the neuropsychiatric family shows stronger enrichment than the control family ... and higher directional consistency"＋长度校正后组间差异消失（1.23 vs 1.06, p=0.24）——两个层面（未校正更强 / 长度解释）均如实呈现，无 cherry-picking。
3. Abstract 与 Concluding remarks 同步降级为 "directionally robust though gene-length-sensitive and not exclusive"。

### P1-1（MAPPED_GENE 分词映射敏感性）— **CLOSED WITH MINOR RESIDUE**

两种替代约定已跑：locus-dedup（±1 Mb 聚类去重）与 author-reported。EA（OR 2.99→3.16 / 2.99）与 SCZ（→3.95 / 1.81）均保持显著，正文 "author-reported-gene mapping retains both (EA OR=2.99, P=0.015; SCZ OR=1.81, P=0.007)" 与 JSON 一致。残留：height 在 locus-dedup 下更显著一事未报（并入上方 P1-3 残留问题一起处理即可）。

### P1-2（trait 重叠与有效检验数）— **CLOSED**

Jaccard 矩阵均值 0.148（范围 0.069–0.310，最大对 ASD–intelligence 0.31）、Li & Ji (2005) Meff=3、Kaiser=2，已写入 Methods（"effective number of tests 3 by Li and Ji 2005"）并新增引文；Statistical conventions 保留 PRDS 句。处理到位。

### P1-3（阴性对照功效与 CI 完整性）— **CLOSED**

功效表（P1_3_power_table / S5）覆盖 12 性状：n_universe_mapped_genes、k_RD_overlap、expected、OR、95% CI、CI_width_log、Fisher p、BH q 齐全。Crohn n=95、k=3、OR=0.52、log CI 宽 2.31；LDL n=256、k=22、OR=1.54（0.98–2.43）、log CI 宽 0.91——LDL 的 CI 上界 2.43 排除了大幅富集，null 判定可信；Crohn 反向点估计（OR<1）亦支持 null 不是纯功效问题。null 可信的结论成立。

### 额外发现（新 P2，数字不符）

正文 Results（GWAS anchoring 段）："GWAS-mapped genes are longer than unmapped genes for every trait (median 1.6–3.5-fold)"——与 JSON 不符：实际范围 **1.29–4.12**（LDL 1.29、intelligence 1.62、Crohn 1.80、ASD 1.88、height 2.26、T2D 2.28、SCZ 2.48、MD 2.72、neuroticism 2.93、EA 2.78、BIP 3.52、cognitive performance 4.12）。"every trait" 的方向性正确，但区间数字两端均错。要求改为 "(median 1.3–4.1-fold)"。

---

## 汇总

| 原意见 | 判定 | 备注 |
|---|---|---|
| P0-3 长度/映射混杂 | **CLOSED** | 模型设定正确、数据-文本逐数吻合、五处诚实化到位；CMH 残余混杂声明诚实 |
| P0-4 组间检验 | **CLOSED** | 检验合理、置换与解析互证、旧措辞清除、两个层面如实呈现 |
| P1-1 映射敏感性 | **CLOSED（小残留）** | height locus-dedup 显著未报 |
| P1-2 Jaccard/Meff | **CLOSED** | Meff=3 已入 Methods |
| P1-3 功效表 | **CLOSED** | null 可信（LDL CI 上界 2.43） |
| 新 P1（选择性报告） | **仍需修改** | CMH/locus-dedup 下 height、T2D 亦存活，正文未报 |
| 新 P2（数字不符） | **仍需修改** | "1.6–3.5-fold" 应为 "1.3–4.1-fold" |

**审稿人建议**：两处新问题均为局部文字修补（不需重跑任何分析），修复后我对本稿的推荐由 Major Revision 上调至 **Minor Revision**。作者对不利结果的披露姿态（"abolishes the enrichment for all 12 traits" 直接进摘要）值得肯定。
