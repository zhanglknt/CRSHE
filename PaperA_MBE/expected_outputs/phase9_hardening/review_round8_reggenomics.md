# Round 8 盲审报告 — 调控基因组学视角（MBE）

**Reviewer:** r8-reggenomics（人类调控基因组学：HARs/hCONDELs、增强子进化、MPRA/caMPRA 功能验证、脑发育调控）
**稿件:** manuscript_english_v9.md + Supplementary_Text_v9.md（含 S6）+ Main_Tables_v9.md + Figures 4/7（figures_v2）
**评审性质:** 第 8 轮模拟盲审，此前未见过此稿，只审不改。

---

## 1. 总评与判定

本文试图定量回答 King–Wilson 问题（编码 vs 调控何者主导人类特异性正选择），其方法论自觉程度在我审过的同类稿件中属上乘：BUSTED/RELAX 分解强化选择与松弛约束、CDS 长度残差化、完整的 LOO 验证矩阵、hCONDEL hg18→hg38 liftOver 双路由修复（S6，含 44 基因 within-gene 金标准审计）、caMPRA-HAR 嵌套的自我拆穿与降级、以及覆盖不对称的定量外推（293→1,906 反事实 + 59% 交叉点）。所有关键阴性结果（hCONDEL-only 子集不显著、caMPRA LOO 校正后不显著、GWAS 长度校正后全灭）均如实上报。这是值得明确肯定的。

但挑剔地看，作为"主要独立验证"的 hCONDEL 富集在三个地方仍未闭环：(i) 全文最强的循环性控制变体（联合 LOO-brain-tau，RD 293→403）的 HAR/hCONDEL 富集从未被报告，而 Figure 4 把单组件 LOO-brain 的数据误标为 "LOO-brain-tau"，构成硬性图-文矛盾；(ii) 用 93% post-hoc power 论证"衰减是信息性的而非功效伪影"是统计学上无效的 observed-power 论证，且其目标效应（全集 OR=3.19）本身被双支持基因抬高，存在自指；(iii) 剔除 HAR 后最干净的 hCONDEL-only 子集（9/90，OR=1.80，P=0.081）并不显著——真正独立于 HAR 的验证仅落在 9 个基因上，"primary independent validation"的措辞强度超出子集分解所支持的水平。另有一处概念张力：标题的 "Dominates by Prevalence" 是一个被作者自己证明 evidence-limited 的类大小比较。

以上问题全部可以在不重算主流程的前提下修复（补一行表、一次正式的异质性检验、改图、改写两段话）。

**判定：Minor Revision**（若联合 LOO-brain-tau 变体下 hCONDEL 富集无法存活，或异质性检验显示双支持 vs hCONDEL-only 显著不同质，则重新评估为 Major Revision）
**总分：6.5/10**

## 2. 分项评分

| 维度 | 分数 | 说明 |
|---|---|---|
| 新颖性 | 7 | 问题经典；新颖性在 comparability-corrected 框架、LOO 矩阵与覆盖量化，而非概念或数据突破 |
| 方法严谨性 | 7 | liftOver 修复链、金标准审计、双路由敏感性均出色；扣在 M1（联合变体缺口+图误标）与 M4（长度混杂仅到 tertile 粒度） |
| 统计分析 | 6 | 家族划分、Holm/BH、Newcombe CI、单侧约定都规范；扣在 M2 的 observed-power 误用与缺失的异质性检验 |
| 结果可信度 | 6 | 富集信号本身可信；"primary independent validation"与标题 "Dominates" 的强度超出证据 |
| 写作与透明度 | 8 | 罕见的自我批评式写作；限定语嵌套过深但几乎所有我想挑的问题作者都已自曝 |
| **总体** | **6.5** | 修后可达 7.5–8 |

## 3. Major comments

**M1. Figure 4 把 LOO-brain 误标为 "LOO-brain-tau"，且真正的联合变体（RD=403）的 HAR/hCONDEL 富集在全文缺失。**
Figure 4a 最后一行标注 "LOO-brain-tau (agr. 0.90)"，数值 HAR OR=2.78 [2.21–3.48]、hCONDEL OR=2.56 [1.66–3.94]；Figure 4b/4d 同标。但这些数值与 Table 3 的 **LOO-brain**（单组件移除）行逐一相同，S4 也明确单组件 LOO-brain 给出 RD=613（Fig 4d 中该点类别规模 ~613），而正文第 126 行与 Methods 给出的联合 LOO-brain-tau（权重 caMPRA 0.5455/nc 0.4545）扩大为 **RD=403**。也就是说：要么图把单组件变体误标成了联合变体（最可能），要么图与 Table 3/S4 用了不同数据却数值完全相同（不可能）。这是硬性图-文矛盾，必须修正。
更重要的是由此暴露的缺口：Table 3 的 LOO 矩阵只有 full + 四个单组件变体，**联合 LOO-brain-tau 变体的 HAR/hCONDEL 富集在全文任何位置都没有报告**——而它恰恰是表达组件（合计 45% RDS 权重）全移除后的最强循环性控制变体，其 SynGO（7 个 FDR 术语）与 GWAS（6/8 性状）结果在正文被正面引用，却没有对应的 HAR/hCONDEL 行。我还注意到一个必须一并说明的概念问题：联合变体下 RDS 仅剩 caMPRA（HAR 嵌套）与 nc 组件（phyloP 耦合），因此该变体对 HAR 检验更循环、对 hCONDEL 检验更耦合——作者需要明确说明该变体在两个验证轴上的可解释性边界，而不是简单地不报告。
**要求：** (a) 修正 Figure 4a/4b/4d 的变体标签；(b) 在 Table 3 增补联合 LOO-brain-tau 行（HAR 与 hCONDEL OR/CI/P），或在正文明确论证该变体下两个检验均不可解释并说明理由。

**M2. "93% post-hoc power"论证在统计上无效，且目标效应自指；应以正式的异质性检验替代。**
正文第 76 行：hCONDEL-only 子集 9/90、OR=1.80、P=0.081（n.s.），作者随即以"93% post-hoc power to detect the full-set effect at this sample size"断言"the attenuation is informative rather than a power artifact"。两处问题：(i) observed power 是 P 值的一对一单调变换，用一个不显著的结果去算它"本来有多大功效"不能提供任何超出 P 值本身的信息（Hoenig & Heisey 2001 的经典批评）；(ii) 更微妙的是，功效计算的目标效应取的是全集 OR=3.19，而全集估计本身被双支持基因（9/22，OR=11.38）显著抬高——用被双支持基因抬高的效应量去论证 hCONDEL-only 子集"本该检测到它"，在逻辑上是自指的。统计上正确的问题是：**双支持子集与 hCONDEL-only 子集的富集是否显著不同质**（Breslow–Day 或对数 OR 之差的 Z 检验；粗略估算 log(11.38)−log(1.80)≈1.85，SE≈0.8 量级，很可能 P≈0.02–0.05，即可检验）。若不同质显著，则诚实的结论是"hCONDEL 验证的强度部分依赖与 HAR 共定位的基因座"，而非"衰减是信息性的"。
**要求：** 删除或重写 post-hoc power 论证，代之以子集间异质性的正式检验，并据实调整"informative"一句的措辞。

**M3. "Primary independent validation"的强度超出子集分解所支持的水平；最干净的 HAR-独立子集不显著。**
hCONDEL 全集富集 OR=3.19 仅由 18 个 RD∩hCONDEL 基因驱动；子集分解（Fig 4c，披露得很诚实）显示信号强烈集中于双支持基因（9/22，OR=11.38，P=1.9e-6），而剔除全部 HAR 后的 hCONDEL-only 子集（9/90，OR=1.80，P=0.081）不显著。此外两证据类型的重叠本身显著高于独立期望（22 vs 10.7，≈2×，Yates P=4.6e-4）——作者用了"largely non-overlapping"这一准确但偏弱的措辞，却没有正面承认这个超额重叠意味着双支持基因座很可能就是被两种方法共同标记的**最长、最深保守的基因座**（与 T3-only 富集互为印证）。这把我引向一个审稿人必然会提的问题：所谓"独立于 HAR 的验证"，目前经验上只落在 9 个 hCONDEL-only RD 基因上，且该子集不显著。我不认为这推翻了验证（方向一致、LOO-nc 下 OR=4.67、双路由与 510 验证子集敏感性均稳健），但"primary independent validation"的措辞需要与子集证据对齐。
另：Figure 4c 的面板标注"hCONDEL-only enrichment (bottom row) excludes HAR re-detection"是图内过度声明——一个不显著的 OR=1.80 不能"排除"任何事，正文措辞谨慎而图注不谨慎，请统一。
**要求：** (a) 正文与 Fig 4 图注/面板标注中把验证强度表述降为"full-set enrichment + 多路由敏感性稳健，HAR-独立子集方向一致但功效不足"；(b) 修正 Fig 4c 的 "excludes HAR re-detection" 标注。

**M4. hCONDEL 富集的长度混杂控制只到 tertile 粒度；建议与 SynGO 分析同规格的长度匹配置换。**
hCONDEL-proximate 基因中位跨度 113.8 kb（vs 全基因组 31.0 kb），RD 类自身中位 133.6 kb；±50 kb 窗口的命中概率随基因长度近似线性上升。作者用跨度 tertile 分层回应：富集局限于长基因层 T3（OR=2.06，P=0.012），两个短层为 null，并解读为"等长基因内 RD 过量仍存在"。方向对，但 T3 内部跨度可从 ~50 kb 到 >2 Mb，tertile 是相当粗的匹配；而同一篇稿件在 SynGO 富集上已经用了 1,000 次 decile-matched 长度零分布 + logistic 校正这一更高规格。**验证主轴（hCONDEL）的长度校正规格反而低于功能富集轴**，这不一致。
**要求：** 对 hCONDEL（及 HAR）富集做与 SynGO 同规格的长度匹配置换检验（对 RD 类做 decile-stratified span-matched null，或在 logistic 框架内以 log₁₀ span 为协变量检验 class×hCONDEL 关联），至少对全集 OR=3.19 与 T3 层 OR=2.06 两个关键数字各补一个长度校正后的 P 值。低成本，闭环 M3。

## 4. Minor comments

1. **S4 与主文对 GWAS 校正家族的描述不一致：** 主文 Methods 与 Results 用统一 48 检验家族（12 性状 × 4 类，含阴性对照），S4 的 LOO 段落却写"BH across the 8 traits × 4 classes family"。请统一为实际采用的家族。
2. **Fig 7d 只画了一个交叉点定义：** 正文给出两个——固定 GD=1,214 的交叉（n*=2,914，59%）与双计数联合建模的等规模点（2,158，43%）；图中仅有 ~2,914 的星标。图注应说明采用的是哪一定义，或两点并标。
3. **联合 LOO-brain-tau 后 RD 293→403 扩大的组成未分解：** 新增的 110 个基因来自哪里（原 unclassified 还是原 GD）？其 SynGO/GWAS 贡献是否集中于新成员（即信号是否由类别膨胀带入）？一句话或一张小表即可，否则"扩类后 7 个 SynGO 术语仍存活"难以排除是新成员的被动带入。
4. **标题 "Dominates by Prevalence" 与自身 evidence-limited 框架的张力：** 摘要与正文已充分限定（59% 覆盖才达 parity），但标题是断言式；且需指出 59% 交叉场景要求非 HAR 元件进入 assay，而作者自己以"活性非人类特异"为由排除了 Shin et al. 的 CNE/variable 层——即交叉场景在作者自定的保守证据定义下不可达，真正可达的 HAR-restricted 天花板是 410（远低于 GD 的 1,214）。这使 prevalence 结论在现有 assay 设计下其实相当稳，但论证链中这两处的关系（交叉点 vs 天花板）值得在 Discussion 用一小段理顺，而不是散落在两处。
5. **Fig 4a 标签 "LOO-nc-cons." 与 Table 3 "LOO-nc" 命名不统一**（同一变体两个名字），请统一。
6. **Fig 7c 图例 "regulation-driven (8/8 FDR-sig.)" 未注明是单侧家族：** 主文说明 intelligence 在双侧参考下 q=0.061（marginal），图内未标注；建议在图注加一句单侧/双侧口径，避免该图被单独引用时过度解读。
7. **hCONDEL 双支持基因（n=22）本身值得一张小表：** 这 22 个基因是全文验证链条上信息密度最高的集合（含 CNTN4、DAB1、ESRRG 等），列出其类别归属（9 RD / 13 其他）与跨度中位数，能让读者直接评估 M3/M4 的关切。S6 已给全部基因清单，补这个子集几乎零成本。
8. **S6 基因清单计数说明清晰（103 symbols vs 112 Ensembl records），值得肯定；** 但请在 S6 同时给出 510 验证子集 → 97 基因这一映射的机器可读路径（目前只有 112 并集的路径）。
9. **正文第 74 行 hCONDEL P 值范围 "2.3e-5 to 5.3e-4" 与 Fig 4a 核对无误；** HAR 范围 "3.6e-25 to 4.7e-5" 无误。Abstract 的 "HAR 1.6–9.4 / hCONDEL 2.1–4.7" 与 Table 3 一致。除 M1 外图-文数字一致性整体良好。

## 5. 已核对一致的关键数字（抽样）

- hCONDEL 全集：18/293 RD 重叠，OR=3.19 [1.90–5.37]，P=7.8e-5（正文、Table 3、Fig 4f 一致，2×2 表复算成立）。
- 敏感性链：交集 70 基因 OR=3.40（P=6.8e-4）、510 验证子集 97 基因 OR=3.28（P=1.4e-4）、双路由 name OR=2.90 / coord OR=3.58、T3 OR=2.06（P=0.012）、nc 耦合 0.667 vs 0.496（MWU P=5.2e-6）、LOO-nc OR=4.67——正文、S6、Table 3 三方一致。
- Fig 4c 子集划分自洽：HAR-only 454 / hCONDEL-only 90 / 双支持 22，合计 566（Fig 7a combined n=566 ✓）；RD 计数 72+9+9=90，加 dual 与 Fig 7a 的 19.1%/17.0%/17.7% 兼容。
- Fig 7a：genome-wide RD+dual 6.1%（304/4,974 ✓）；编码侧 36.6%/30.4%/35.5% 与正文 CI 一致。
- Fig 7d：观测 97=2%、天花板 410@476、交叉 ~2,914（59%）与正文一致（两点定义见 Minor 2）。
- LOO 类规模：293/800/148/110/613（S4 复算声明）与 Table 3、Fig 4a/d 的类别规模轴一致。

## 6. 什么改动能让我提高分数

完成 M1（修图 + 补联合变体行）、M2（异质性检验替代 observed-power）、M4（hCONDEL 长度匹配置换）三项计算/图形修复，并按 M3 对齐"primary independent validation"的措辞——这四项完成后方法严谨性与可信度均可升至 8，总体 7.5–8。若 M1(b) 显示联合变体下 hCONDEL 富集不存活或 M2 异质性检验显著且作者拒绝下调措辞，则维持 Major Revision 级别的保留意见。
