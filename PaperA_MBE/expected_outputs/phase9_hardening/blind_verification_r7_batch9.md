# R7 (batch 9A–9F) 盲交叉验证报告

验证员：general-purpose-30（独立统计盲验证，不信任中间报告，全部从原始数据重算）
日期：2026-09-23
对象：`results/paper/manuscript_english_v9.md`（306 行版）+ 4 个新计算 JSON/CSV + S5 重建 + 补充 S6 + Table 3 脚注
方法：每项给出我的独立复算值 vs 稿件值。原始数据来源：`family_census_r7.csv`、`loo_joint_gwas_matrix.csv`、`TableS5_anchoring_master.csv`、`length_correction_syngo(_terms).csv/json`、`hcondel_sensitivity_r7.json`、`relax_derivation_r7.json`、`loo_joint_brain_tau_summary.json`、`stratified_rerun_pvalues.csv`、`gene_classification_v7.csv`、`hcondel_gene_mapping_fixed.csv`、`busted_all_results.csv`。

---

## 1. SynGO 长度校正整合 — PASS

| 稿件声明 | 我的独立核算 | 判定 |
|---|---|---|
| 7/11 SynGO 词条同时通过两种长度校正 | `length_correction_syngo_terms.csv`：SynGO 11 行，survives_both=True 7 行 / False 4 行 | PASS |
| 9/11 单独通过 logistic | 同表 p_bh_logistic<0.05 恰为 9/11（0.0038/0.0038/0.0070/0.0070/0.0263/0.0383/0.0353/0.0383/0.0407 通过；0.0642、0.0760 未过） | PASS |
| 18/18 GO BP 全部通过 | 同表 GO_BP 18 行 survives_both=True 18 行 | PASS |
| LOO-brain/tau SynGO 7/5，GO BP 15/1；联合 LOO 403 基因 → 7 SynGO、GD 0 | `length_correction_syngo.json` 各键值一致；正文与讨论（"7, 5, and 7"）一致 | PASS |

## 2. hCONDEL 敏感性整合（以 JSON 为准）— PASS

| 稿件值 (L76) | JSON / 我的复算 | 判定 |
|---|---|---|
| 名称∩坐标交集 n=70，OR=3.40，P=6.8×10⁻⁴ | 3.4039 / 6.77×10⁻⁴（k=12/293） | PASS |
| 510 条经序列验证子集 → 97 基因，OR=3.28，P=1.4×10⁻⁴ | liftover 509/510；并集 97；3.2803 / 1.43×10⁻⁴（k=16/293） | PASS |
| 仅 hCONDEL 亚组 9/90，事后功效 93% | posthoc_power_mc=0.9329；n_for_80pct=57 | PASS |
| 跨度三分位 T3 OR=2.06，P=0.012 | T3 17/73，OR=2.055，P=0.0124（T1 0/15、T2 1/24 梯度一致） | PASS |
| nc 耦合 0.667 vs 0.496，MWU P=5.2×10⁻⁶；top decile 23.2%/9.7%；LOO-nc OR=4.67 | 0.667471/0.4962805，P=5.20×10⁻⁶；4.669 | PASS |

## 3. 统一 48 家族 + S5 重建 — PASS

- `family_census_r7.csv` 144 行；Intelligence RD：two-sided q48=**0.061000**（稿件 "q = 0.061"），one-sided q48=**0.04322**（稿件口径 0.0432），sig_unified_48(two-sided)=False、one-sided=True。PASS
- 单侧 48 家族下 8/8 神经性状 RD 全显著（q 7.7e-6–0.043）；双侧下 7/8（intelligence 临界）。PASS
- `TableS5_anchoring_master.csv` 两副本均 **192 行** = 4 变体（full/LOO-brain/LOO-tau/LOO-brain-tau）× 12 性状 × 4 类；首行数值与 census 完全一致；xlsx S5 表两副本均 194 行（192 数据 + 表头 + 注释行）。PASS
- "Meff" 全文 0 命中；无残留 32/36-test 家族表述（6 种正则扫描均 0 命中）。PASS

## 4. MDD 修正 — PASS（附口径备注 N1）

- 稿件 "OR = 0.56, two-sided q = 0.035"：census 全分类 GD-strict OR=0.5631，p_two=0.00977，q(S5-36)=**0.035167** ✓。
- "LOO q = 0.11–0.15"：LOO-brain q48=0.1426/q36=0.1070；LOO-tau q48=0.1478/q36=0.1109 → 区间 0.107–0.148 覆盖 "0.11–0.15" ✓；联合 LOO q48=0.0811 同样不显著，与 "not under LOO reclassification" 定性一致。
- bipolar "non-significant depletion"：q36=**0.5432** ✓（数值未印出，无矛盾）。
- 独立 Fisher 复算（k=24, n=154, K=1214, N=4974）：OR=0.56315 与 census 完全一致。

## 5. RELAX 推导链 — PASS

`relax_derivation_r7.json`：1,660 valid → 728 FDR-sig → 124 K<1 / 604 K>1；117 relaxed = 52 GD-relaxed + 65 neutral；7 个 K=0 退化 → 5 neutral + 2 GD。稿件 L64 链条与此一致。

## 6. 联合 LOO — PASS

- 权重 doan 0.5455 / nc 0.4545；RD **403**、GD 1207、GDr 49、dual 21（`loo_joint_brain_tau_summary.json`）✓。
- GWAS 矩阵 48 行：RD 行 SCZ q=2.4e-6、EA 5.4e-16、MDD 0.0228、Bip 3.4e-4、Neuro 0.0070、CogPerf 4.0e-4 显著；ASD p=0.0586、Intelligence p=0.0928 不显著 → **6/8**，与稿件 "retains six of eight (ASD and intelligence not significant)" 完全一致。独立 Fisher 复算 MDD：k=23,n=154,K=403,N=4974 → OR=2.0514 与矩阵一致。
- LOO-brain 6/8、LOO-tau 5/8 我从 census 逐行复核一致；讨论区 "six, five, and six" 三处一致。
- 对照性状区间（R6 FAIL #3 修复确认）：Height 2.215/1.878/2.335 → "1.9–2.3" ✓；T2D 2.148/1.770/2.336 → "1.8–2.3" ✓。

## 7. batch 9A 文字组 — PASS（附 FAIL-F1）

- chi-bar-square 零分布签名段落在位；BH 下界论证完整。PASS
- Tier 4 不校正理由（三个预设假设、零为所声明、FWER  inflated type-II）+ 六检验家族 "q < 0.02" + GD hCONDEL 消耗 "nominal and would not survive"（two-sided P=0.037）均在位且与我的 R6 复算一致。PASS
- GDS 有效权重披露："combined nominal weight (0.65) effectively represents a single evidence type" 在位。PASS
- φ = 0.052（未校正 χ²）、Yates P = 4.6×10⁻⁴、22 observed vs **10.7 expected**（112×476/4974=10.72，22/10.72=2.05 ≈ "≈2×"）✓。
- EA 4×2 χ²=12.33，df = 3，P=0.006 ✓（我 R6 复算 12.33/3/0.0063）；"59% vs 48%" ✓。
- Shao：4/18、期望 6.1（18×1690/4974=6.116）、OR=0.55（我算 0.5541）、one-sided P=0.21（我算 depletion p=0.2135）✓。
- seeds 20260915 / 20260923（各 1–2 命中）/ 20260917（trait permutation）均在位 ✓。
- 508 N2A 措辞在位（两处上下文完整）✓ —— 但见 F1。

## 8. 摘要 — 数字全部 PASS；2 处备注

逐项复核：1,690 (34.0% = 33.98%)、604、52、1,214/293/11、HAR LOO OR 1.6–9.4（我 R6 算 1.614–9.388，含联合 2.775 在界内）、hCONDEL 2.1–4.7（2.132–4.669）、11 SynGO、7/11 + all 18、八性状、4.7×、q>0.44、2–10%（97/4974=1.95%，476/4974=9.57%）、~59%、24.4% vs 5.9%（1214/4974=24.41%，293/4974=5.89%）——全部一致。
- **FAIL-F2**："604 (12.2%)"（摘要与正文 "approximately 12.2% of analyzed genes (604/4,974)"）：604/4974 = **12.143%**，规范舍入为 12.1%，差 0.06pp。正文同句已给确切分数 604/4,974，故读者可自纠，但百分比应改 12.1%。
- **备注 N4**：摘要词数我用三种分词口径计得 258–264 词（空白切分 262、去 markdown 261、含字母数字 token 258），与 "248 词" 规格不符。词数不印在稿件中，属元信息问题，供参考。

## 9. S6 专节 — PASS（附备注 N3）

- 583 目录、510 序列验证、581/583 liftOver 成功（2 失败点名 hCONDEL.449/577）、无 panTro2 回退 ✓。
- 金标准 44 基因：旧映射 7 直接 + 32 错位 vs 修复后 40 直接 + 2 延伸 + 2 错位（"40/2/2 vs 7/32"）✓。
- 名称 80 / 坐标 102 / 一致 70 / 并集 112；Jaccard 0.625 = 70/112 精确 ✓；name-only OR=2.90、coord-only OR=3.58、union 3.19 与我 R6 复算一致。
- 三个基因清单与 `hcondel_gene_mapping_fixed.csv` 在符号级**逐一精确匹配**（intersection 70/70、name-only 10/10、coord-only 23 个有名称符号全部一致）。
- **备注 N3**：CSV 坐标单路径唯一符号实为 24 个（其中 1 个空白符号 = 无 HGNC 名的 Ensembl 基因），行级 32 条 Ensembl 记录；并集 112 是 Ensembl 级。S6 "Coordinate-route only (n = 23)" 计的是有名符号，70+10+23=103≠112。文中 "some symbols map to multiple Ensembl gene IDs" 部分覆盖此口径，但复现者可能困惑，建议加注 "(23 named symbols; 24 genes including one without an HGNC symbol; 32 Ensembl records)"。

## 10. Table 3 脚注 — PASS

`Main_Tables_v9.md` 第 34 行脚注：full 133.6 / LOO-caMPRA 82.0 / LOO-tau 160.2 / LOO-brain 95.4 / 联合 105.8 kb —— 与 `hcondel_sensitivity_r7.json` 中位数（133,590 / 82,017.5 / 160,195 / 95,420 / 105,793 bp）逐一相符 ✓。

## 11. 一致性扫描 — PASS（除 F1/F2）

- 无残留 36/32-test 家族、无 Meff、Fig 7c 图例 "unified 48-test one-sided family (OR 1.9–4.0)"（实际 1.95–4.00）✓。
- srv 超时偏差复核：以生产 bh_fdr_sig（1,690 集）计，400 样本 200/400=**50.0%**（分层设计如此），110 个 srv 超时基因中 76 个 = **69.1%** → 稿件 "69% versus 50%" ✓（Fisher OR=2.99，P=3.6×10⁻⁶，"disproportionately production-significant" 成立）。"~26–28%" 与我 R6 两种外推（26.8%/27.6%）一致 ✓。
- **FAIL-F1**：Methods HAR 段 "our 508/3,171 = 16.1%"：508/3171 = **16.02%**（应为 16.0%）；16.1% 对应分母 3,161（508/3161=16.07%）。同一段 Methods caMPRA 小节 "3,161 lifted; 508 (16.1%)" 自洽，惟 "508/3,171 = 16.1%" 这一等式内部不一致，差 0.1pp。
- **备注 N1**：MDD "two-sided q = 0.035" 是 S5-36 家族值；统一 48 家族下为 0.047（同样 <0.05，定性不变）。该句未标注家族，而上下文正强调统一 48 家族，建议标 "(36-trait family)" 或改用 0.047。
- **备注 N2**：Discussion "of the 1,690 BUSTED-significant genes, 36.4% show FDR-significant intensified selection and 7.0% relaxed constraint" —— 36.4%/7.0% 的分母是 1,660 个有效 RELAX（604/1660=36.39%，117/1660=7.05%），非 1,690（35.7%/6.9%）。Results 节已明确 1,660，此句口径偏松。
- **备注 N5**：LOO-brain 下 intelligence "marginal"（two-sided q48=0.107）措辞偏松但可接受。

## R6 遗留 FAIL 修复确认（顺带复核）

L74 HAR p 区间 "3.6×10⁻²⁵ to 4.7×10⁻⁵" ✓；caMPRA LOO "OR = 1.64, CI 1.02–2.63, P = 0.032, corrected P = 0.16" ✓（我算 1.640/0.0315）；srv "42.6% to 34.6%" ✓；Height/T2D 区间 ✓（见 §6）；S6 已取代旧 183 基因坐标映射的表述 ✓。

---

## 总评

**判定：通过（2 处小数值 FAIL + 5 处口径备注）。** R7 batch 9A–9F 的全部头条数字——SynGO 7/11+18/18、hCONDEL 四项敏感性、统一 48 家族（intelligence q=0.061/0.0432）、MDD 修正、RELAX 链、联合 LOO 403/6-5-6、S5 192 行双副本、S6 映射链、Table 3 脚注五个中位数、srv 69% vs 50%——均经我从原始数据独立复算确认，无一实质性错误。两处 FAIL 均为 ≤0.1pp 的舍入/分母口径问题（508/3,171=16.0% 而非 16.1%；604/4,974=12.1% 而非 12.2%），不影响任何结论，建议投稿前顺手修正。备注 N1–N5 为可选的口径标注改进。
