# R6 修订复核 — molevol-methods 审稿人

**复核范围**: 本人在 `review_round6_molevol_methods.md` 中提出的 2 项 P0、5 项 P1、7 项 P2。证据核查基于 `response_round6.md`、修订后 `manuscript_english_v9.md`、`Supplementary_Text_v9.md`，以及 `results/phase9_hardening/` 下四个新 JSON（neutral_sim_calibration、bh_bounded_support_check、stratified_rerun_sensitivity、topology_sensitivity）与更新后的 tier4_empirical_disjointness.csv。未修改任何稿件文件。

**总体结论**: 两项 P0 均已实质性解决（证据链从"披露"升级为"测量"），P1 全部关闭，P2 六项关闭、一项（Zenodo DOI）待投稿前完成。发现 **2 处小的文本残留错误**需要修正（均不改变结论），详见末尾。

---

## P0 复核

### P0-1（有界双峰 p → BH 校准）: **CLOSED**

三段式证据链完整且方法正确：

1. **解析论证**（bh_bounded_support_check.json）：截断机制下 BH 拒绝集数学不变（临界值 i/n×0.05 ≤ 0.05 << 0.5，0.5 处点质量不可能跨越阈值；<0.5 的 p 值相对次序不变）——推导正确，1,690→1,690。压缩机制（最不利情形）下 p\*=2p 校正仅 1,690→1,639（−3.0%）。两种机制的效应方向与幅度均被量化，正是我 P0-1 要求的"要么模拟、要么有界支撑校正"，两者都做了。
2. **中性模拟**（neutral_sim_calibration.json）设计审查——抽样与参数均合理：evolver M0 ω=1（最激进的可信零假设；真实零基因多受纯化选择 ω<1，只会更保守，因此该零下的保守性结论只会被加强）、κ=2、F3x4 来自全部 4,974 条真实比对、逐基因树、**真实缺失数据掩码**（关键优点：缺失模式本身不造成假阳性）、同一 BUSTEateD 管线与 300s 超时。机制判定为截断（45.8% 点质量、0 个 p>0.5）——与真实数据 39.4% 点质量形态一致；名义尾 1.75% < 5% → 保守；经验 BH FDR = 0/400。真实数据 35.8% p<0.05 vs 零假设 1.75% → 小 p 富集为真实信号/模型误设，非有界支撑伪影。
3. 稿件文本（Limitations (2)、Results L60、Discussion L142）如实整合，且 Results 仍保留 Storey 检查的"internal consistency only"限定——措辞分寸恰当。

残留小问题（不阻塞关闭）：模拟 400 个基因的超时/完成数未在 JSON 中报告（从比例合计推断全部完成）；建议在补充材料补一句。

### P0-2（分层重跑测量 34% 上限）: **CLOSED（附 1 处必须修正的事实性错误）**

设计审查：400 基因按 显著性 × 长度三分位 × LRT 三分位 分层（200+200），四组对照（baseline/CpG/topo/srv）。优点：
- baseline 复跑（49.4% vs 49.9%，p 相关 0.97，逐基因一致 98.5%）证明环境可复现性，是可信的内部对照；
- 正确处理了 BH 家族大小依赖（400 基因族内 BH + embedded 全宇宙族诊断双口径）；
- 比对质量描述统计已交付（我 P0-2 的子要求）：500 基因样本，平均成对一致性 86.3%（p5 74.9%）、每序列 gap 比例中位数 55.8%、每比对 9–10 物种。gap 比例偏高但已如实披露；
- CpG 屏蔽定义明确（含跨密码子边界 CpG，全物种掩码，平均 8.7% 密码子）；
- srv 为唯一实质敏感项（共同 289 基因上 −8.0pp，相对 −19%，外推 ~27%），稿件表述改为"27–34% 取决于 srv 建模"（Methods L32 / Results L60 / Discussion L142 / Limitations (2) 四处一致）——这正是我要求的"把上限变成测量"。

**必须修正的事实性错误（小，但方向反了）**: Discussion L142 与 Limitations (2) 称 srv 组 "110 of 400 runs timed out, **biased toward smaller genes**"。我从 stratified_rerun_sensitivity.json 复算：srv 完成子集的生产检出率为 42.8%（124/290），而全样本为 50%（200/400）——**超时基因中有 76/110（69%）是生产显著基因，即超时富集于显著基因**。鉴于稿件自己报告显著性与 CDS 长度正相关（r=0.18），超时更可能偏向**更长**基因（baseline 组 7 个超时中 4 个显著 ≈ 57%，与样本比例一致，说明超时偏倚是 srv 计算量所致）。该 JSON 全库检索亦无 "smaller" 任何佐证。此括注方向错误，且遗漏了更重要的可迁移性警告：−8.0pp 是在显著基因被稀释的子集上估计的。修正建议：改为 "110 of 400 srv runs timed out; the timed-out set is enriched for production-significant genes (69% vs. 50%), so the −8.0-point estimate derives from a significance-depleted subset"，并可补一句绝对法外推（34%−8pp ≈ 26%）与相对法（~27.6%）一致，说明 27–34% 区间对该偏倚稳健。

## P1 复核

| 项 | 判定 | 证据 |
|---|---|---|
| P1-1 子集含于参照 + CSV 列误标 | **CLOSED（附 1 处残留）** | CSV 已增 complement_* 列、真 Newcombe 列（本人手算 HAR vs 补集 Newcombe = (−1.6%, +7.5%)，与 CSV (−0.0158, 0.0748) 逐位吻合），旧列正确重标为 wilson_bound_subtraction；Results L124 改 complement 口径并明示理由；Fig 7a 图例同步。**残留**：Methods L52 仍写 "compared with the genome-wide rate (Fisher exact)"，与 L124 矛盾，需一处改词 |
| P1-2 RELAX K>1 方向歧义 | **CLOSED** | L64 完整限定（K 缩放含纯化分量、适应性解读依赖与 BUSTED 的合取、操作性定义非金标准）；结论段 12.2% 同样操作化 |
| P1-3 权重 a priori 不可验证 | **CLOSED** | L36 改为 "fixed ... before any downstream enrichment was inspected; there is no preregistration by which this could be independently verified, so weight-choice sensitivity is quantified directly"——正是我要求的表述；dual GDS>0.3 的事后性已明示 |
| P1-4 specificity 弱化 | **CLOSED（超出要求）** | 合并组的长度混杂分析（logistic 校正后 12 性状全不显著、CMH 分层 EA/SCZ 存活、组间检验 z=3.49/置换 P=0.009、长度校正后 1.23 vs 1.06 也如实报告）比我的要求更进一步；结论句 "unadjusted ... directionally robust though gene-length-sensitive and not exclusive" 分寸恰当 |
| P1-5 MEME 选择偏倚 | **CLOSED** | L32 明示最长比对最易超时、7/8 不可读作总体率 |

## P2 复核

| 项 | 判定 |
|---|---|
| P2-1 Selectome 单双侧 | CLOSED（L80: one-sided P = 0.41 (two-sided 0.64)） |
| P2-2 "ten strongest terms" 指代 | CLOSED（Fig 7b 图例明确 five SynGO + five GO BP） |
| P2-3 RELAX run-list 简化 | CLOSED（L32 简化，S1 保留核算） |
| P2-4 ±50 kb 窗口敏感性 | CLOSED（Limitation (10)：±25/±100 kb 单调稀释，全部 P ≤ 3e-4） |
| P2-5 Fig 2e omega 注 | CLOSED（备择模型 ML 估计 + "不单独暗示 episodic selection" 限定） |
| P2-6 GTEx 日期措辞 | CLOSED（正文 L28 与 S1 L21 逐字一致） |
| P2-7 Zenodo DOI | **STILL OPEN**（L188 仍为 "DOI to be assigned upon release"；随交付链重建完成即可，投稿前必须填入实际 DOI） |

## 顺带核查（与我领域相关的合并组修复）

- 替代拓扑（(human,gorilla) 姐妹群，389/400 完成，−1.3pp，92.1% 保留）已正确并入 ILS 段（L172），SPILS 机制引文（Mendes & Hahn 2016）恰当，"实际影响小但建议 gene-tree-aware 重分析"的分寸正确。
- hCONDEL 坐标修复（183→112）传播一致：Fig 7a 图例 n=112/566、Tier4 CSV、Limitation (10) 的窗口 OR 序列（3.95/3.09/2.36）与新基数自洽。

## 复核结论

**两项 P0 关闭，稿件在我专长维度（dN/dS 检测统计基础、比对/拓扑/CpG 混杂、FDR 校准）已达到 MBE 可发表标准。** 剩余事项：2 处文本残留修正（srv 超时方向括注；Methods L52 "genome-wide" 改 "complement"）+ Zenodo DOI 落地。修正后无需再审。

*复核人: r6-molevol-methods（2026-09-18）*
