# R7 batch 9G 盲验证报告（WSL 侧四块结果数字核销）

验证员：general-purpose-30（独立统计盲验证；以 JSON 为权威，不信任转述）
日期：2026-09-23
源 JSON：`results/phase9_hardening/{srv_primary_estimate, alignment_rerun_sensitivity, neutral_sim_parametric, gard_troubleshoot}.json`
被验证文件：`manuscript_english_v9.md`、`cover_letter_v1.md`、`supplementary/Supplementary_Text_v9.md`

---

## 1. srv Yes 完成版主估计 — **FAIL（口径框架混用，0.3pp）**

| 稿件表述（Discussion + Limitations (2)，两处同文） | JSON / 我的复算 | 判定 |
|---|---|---|
| "all runs finished / all 400 runs completed in the final re-run" | n_srv_valid=400, n_missing_after_rerun=0 | PASS |
| "49.4% ... on the 393 comparable genes" | rate_srvNo_on_common=0.493639 = 194/393 | PASS |
| "**41.0%** on the 393 comparable genes" | on-393 实为 **40.7%**（rate_srvYes_on_common=0.407125 = 160/393）；41.0% 是 rate_srvYes_full400=0.41（164/400，全 400 框架） | **FAIL** |
| "**−8.4 points**" | 共同基因框架精确差 = 49.364−40.712 = **8.65pp → −8.7**；−8.4 来自混合框架 49.364−41.0=8.36 | **FAIL（随上）** |
| "a 17.5% relative reduction" | 1−0.407125/0.493639 = 0.17526 → 17.5%（共同框架比率） | PASS（但与同句 41.0% 不自洽：49.4%×(1−0.175)=40.7%，非 41.0%） |
| "extrapolating to ≈28.0% universe-wide" | 0.339767×(0.407125/0.493639)=0.280220 → 28.0%；universe_rate_bounds=[0.28022, 0.28022] | PASS |
| "roughly 1,394 of 4,974 genes" | 0.280220×4974 = 1393.8 → 1,394（point_estimate_genes_of_4974=1394） | PASS |
| "superseding the earlier −8.0-point estimate from a 289-gene timeout-biased subset" | r6_partial_estimate: n_common=289, 0.426→0.345 (−8.1pp≈−8.0)；超时偏差方向（稀释显著性）与 R6 我的独立复核一致 | PASS |
| cover letter "point estimate at 28.0%" | 同上 | PASS |

**修正建议**（两处：Discussion L142 段、Limitations (2)）：将 "41.0%" 改为 "40.7%"、"−8.4 points" 改为 "−8.7 points"（共同 393 基因框架，与 17.5% 和 28.0% 外推自洽）；或保留 41.0% 但改写框架，如 "--srv Yes gives 41.0% across all 400 completed runs (40.7% versus 49.4% on the 393 comparable genes, −8.7 points, a 17.5% relative reduction)"。涉及 JSON 字段：`rate_srvYes_on_common` vs `rate_srvYes_full400`。

## 2. codon-aware 重比对 — PASS

| 稿件表述 | JSON / 复算 | 判定 |
|---|---|---|
| "49.4% to 64.1% (+14.8 points on the same 393 genes)" | rate_other_common=0.641221（64.1%）、rate_baseline_common=0.493639；delta_pp=14.7583 → +14.8 成立（精确值非 64.1−49.4=14.7） | PASS |
| "80.4% of baseline-significant genes retained" | retention_of_baseline_sig=0.804124（156/194） | PASS |
| "59.1% of columns removed on average" | frac_removed_mean=0.591175 | PASS |
| Methods (v)："protein-based MAFFT followed by back-translation, removing columns with >50% gaps" | method.realignment 一致（>50% gap 列剔除、端部修剪、终止子 NNN 屏蔽） | PASS |
| cover letter "moves detection upward to 64.1%" | 同上 | PASS |
| "alignment noise attenuated power rather than inflating false positives" | 方向与 switch_stratification（各 tercile gained>lost）一致，属合理解读 | PASS |

## 3. 参数化 bootstrap — PASS

| 稿件表述 | JSON / 复算 | 判定 |
|---|---|---|
| "empirical FDR of 0.26% (1/388)" | empirical_BH_FDR_alpha0.05=0.002577；n_BH_rejections_0.05=1；n_busted_results=388；1/388=0.258% | PASS |
| "73.7% of bootstrap p-values at the 0.5 boundary" | p_eq_0.5_point_mass=0.737113（hist_deciles 286/388=0.737，总和 388 闭合） | PASS |
| "median ω = 0.54, κ = 3.02, F3x4 frequencies, fitted branch lengths" | omega_null_mean_median=0.53751 → 0.54；kappa_median=3.01589 → 3.02；F3x4 与拟合枝长见 description | PASS |
| 三条局限："omits indels, single-ω null approximation, does not model SRV" | limitations[] 三条逐字对应 | PASS |
| （对照）"empirical null FDR is 0/400 under the evolver null"；Methods "45.8% at exactly 0.5" | comparison_r6_fixed_param: FDR=0/400、mass=45.8%，与 R6 链一致 | PASS |

## 4. GARD — PASS（附 2 条口径备注）

| 稿件表述 | JSON / 复算 | 判定 |
|---|---|---|
| Methods L32 + Limitations (3)：standard 模式不收敛诊断为配置病理（穷举候选断点 vs 600 s 上限、100% 可变位点高 gap 比对），非数据问题 | diagnosis.root_cause 一致；生产 50/50 零字节失败；--rv None 无效；900 s 复现仍超时 | PASS |
| "Faster-mode re-runs completed all 50 alignments" | rerun_results: n_valid=50/50（3600 s 超时、JSON>50B 判据） | PASS |
| SI S2："breakpoint at site 200 with Δc-AIC ≈ 1,519 between competing models in replicate runs" | ENSG00000177119：Faster rc=2 → site 200、Δc-AIC 1519.8；rc=4 → site 200、1519.3（复现一致） | PASS |
| "breakpoints co-localize with alignment-block artifacts（>100 substitutions per site）" | caveat：分区树 >100 subst/site、拓扑打乱，判定为比对块伪影 | PASS |

- **备注 A（SI 口径）**：SI S2 称 GARD 以 "nucleotide mode with 3–4 rate classes" 运行；JSON 生产命令为 `--type codon --rate-classes 4`（复现测试为 rc=2 与 rc=4）。建议改为 "codon mode" 且速率类写 "4 (2 in one reproduction configuration)" 或 "2–4"。不影响结论。
- **备注 B（JSON 内部，未入稿）**：`delta_cAIC_median`=1255.1387 与 `time_elapsed_median_s`=593 均为 50 个值的**上中位数**（第 26 顺位，恰好都是 ENSG00000064763）；偶数样本常规中位数应为 1210.89 / 592。该两值未出现在任何稿件文本中，无需行动，仅供 JSON 维护者知悉。

## 5. 一致性扫描 — PASS

| 禁用旧口径 | MS | CL | SI |
|---|---|---|---|
| "27–34%" / "27-34" | 0 | 0 | 0 |
| "42.6%" | 0 | 0 | 0 |
| "34.6%" | 0 | 0 | 0 |
| "289 comparable" | 0 | 0 | 0 |
| "110 of 400" | 0 | 0 | 0 |

说明：Limitations (2) 中 "superseding the earlier −8.0-point estimate from a 289-gene timeout-biased subset" 是对被取代历史估计的明确标注，属有意保留的溯源表述，不视为残留。

## 6. 算术自检 — PASS（第 4 项除外，见 FAIL-1）

| 自检式 | 结果 | 判定 |
|---|---|---|
| 0.34×(1−0.175) ≈ 0.2805 → 28.0% | 0.2805；精确链 0.339767×0.824742=0.280220 → 28.0% | PASS |
| 0.2802×4974 ≈ 1,394 | 1393.8 → 1,394 | PASS |
| 64.1−49.4=14.7 vs 稿件 +14.8 | JSON delta_pp=14.7583（精确值 64.122−49.364），按精确值 +14.8 成立 | PASS |
| 49.4−41.0=8.4 | 仅混合框架成立；共同 393 框架精确差 8.6514 → 8.7 | **FAIL-1** |

---

## 总评

**判定：1 项 FAIL（srv 主估计的 41.0%/−8.4pp 与 "on the 393 comparable genes" 框架混用，幅度 0.3pp）+ 2 条口径备注，其余全部 PASS。** codon-aware 重比对、参数化 bootstrap、GARD 三块与 JSON 完全一致（含 SI S2 的 site 200 / Δc-AIC≈1,519 复现性描述）；旧口径残留扫描三文件均为零；外推算术闭合。修正 FAIL-1 只需在 Discussion 与 Limitations (2) 两处把 41.0%→40.7%、−8.4→−8.7（或改写框架说明），即可使 49.4%→x%、相对降幅 17.5%、外推 28.0% 三者完全自洽。
