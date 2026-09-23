# Review Round 7 — Reviewer: Statistical Genomics & Reproducibility Audit (MBE)

**Manuscript:** Coding Selection Dominates by Prevalence, Regulatory Selection by Specificity: A Comparability-Corrected 10-Primate Analysis (v9)

**Reviewer expertise:** multiple testing / test-family design, circularity & information leakage, effect size & power, sensitivity-analysis design, reproducibility.

**声明：** 本轮为全新盲审，未参考既往任何轮次意见。文中行号指 `manuscript_english_v9.md`；"Supp" 指 `Supplementary_Text_v9.md`；表号指 `Main_Tables_v9.md`。

---

## 1. 总评与推荐决定

**一句话总评：** 这是一篇在循环性审计与覆盖不对称量化上接近领域范本的稿件，算术与内部一致性经我抽查全部通过；但"长度校正杀死 GWAS 锚定后，突触基因集富集成为唯一幸存的长度稳健信号"这一新的结论支柱本身从未接受任何长度控制检验，且全文的检验家族划分存在互相重叠、可被择优报告的结构——在这两点解决之前，标题级结论的证据链不闭合。

**推荐决定：Major Revision**（不需新数据，所有要求的分析均可用现有管线在数日内完成；但有一处结论级声明目前是无证据的断言，且检验家族结构需要系统性重整，故不宜定为 Minor）。

---

## 2. 评分表（1–10 整数）

| 维度 | 评分 | 简评 |
|---|---|---|
| 新颖性 | 8 | 组件均成熟（BUSTED/RELAX/LOO/富集），但"可比性校正 + LOO 全矩阵 + 覆盖反事实"的整合审计框架对领域有示范价值；conditional answer 的提法诚实且恰当。 |
| 方法严谨性 | 8 | 分类与验证侧的防循环设计（证据排除 + LOO）是全文最强部分；扣分在 GDS 权重共线性与 LOO 仅做单分量移除（见 M3/M5）。 |
| 统计分析 | 6 | 单点算术我复核无误（Selectome OR=1.21、φ=0.052/0.434、HAR OR=4.15、Shao OR=0.55、Holm 0.032×5=0.16 等均自洽）；但家族划分重叠可择优（M2）、核心"长度稳健"声明未检验（M1）、若干非零声明未校正（M4）。 |
| 结果可信度 | 7 | LOO-hCONDEL 链与覆盖敏感性令人信服；GWAS 锚定已被作者自己的长度校正证伪并诚实披露；剩余可信度缺口集中在 M1 与 M3。 |
| 写作清晰度 | 9 | Statistical conventions 一节（L54）堪称模板；局限披露异常完整。扣分仅在两处数字歧义（117 vs 52；OR 方向）。 |
| **总体** | **7** | 修好后可达 8–9；当前被 M1/M2 压在 7。 |

---

## 3. Major comments

### M1. "长度稳健的突触富集"是全文新的结论支柱，却从未接受任何长度控制检验——这是当前证据链上最大的洞

稿件已经（值得赞扬地）证明：RD 基因长度为中位 133.6 kb，是其余基因（28.6 kb）的 4.7 倍（L126：Wilcoxon P = 1.5 × 10⁻⁴⁵）；RDS 的非编码保守分歧分量与基因跨度 Spearman ρ = 0.47，即**分类规则本身按构造偏向长基因**；GWAS 锚定在 log₁₀(span) 校正下全线失守（"adjusted OR 0.33–1.64, all q > 0.44"，L126）。作者据此把结论退守到"leaving the synaptic gene-set enrichment as the length-robust functional signal"（Abstract，L12；L166 同义重复）。

但 SynGO/GO-BP 的超几何富集（L114）以 4,974 基因为背景，**不含任何长度项**。而突触/神经发育基因恰是基因组中最长的一类基因（neurexin、CNTN、cadherin 家族均为 megabase 级）——这恰是 GWAS 侧混杂的同一机制。逻辑链条是：nc 分量 ∝ 长度 → RD 类富集长基因 → 突触基因是长基因 → RD 富集 SynGO。此链条的每一环稿件自己都提供了数字。L126 末尾的辩护（"the synaptic gene-set enrichment … is a property of class membership rather than locus mapping"）只回应了映射约定，没有回应长度。**因此"length-robust"目前是一个断言，不是一个结果。**

**要求（任一即可，成本均低）：**
(a) 从 unclassified 基因中按 log₁₀(span) 分层抽样构建长度匹配背景（如 10 次重抽），重跑 11 个 SynGO 与 18 个 GO-BP 术语的超几何检验，报告匹配后仍显著的术语数；或
(b) 对每个幸存术语做 logistic 回归：term membership ~ class + log₁₀(span)，报告 class 的调整后 OR 与 q；或
(c) 长度三分层内的分层超几何/CMH 检验（与 GWAS 侧已做的 CMH 完全同构，管线可直接复用）。
若富集在长度控制后存活——哪怕是缩减后的子集——这将真正成为全文最硬的生物学结论；若失守，Abstract 与 Concluding remarks（L180）中"carries the specific functional signal"的表述必须再次降级。无论结果如何，这一分析是结论成立的前提，不是可选的稳健性装饰。

### M2. 检验家族互相重叠且规模不一，存在"择优家族"空间；S5 的 3 类结构与主文 4 类家族直接矛盾

主文 L54 声明的家族划分本身是好的，但拼起来看有问题：

1. **Tier 2 主家族 = 8 traits × 4 classes = 32 检验（L50）；阴性对照家族 = 12 traits × 4 classes = 48 检验"per variant"（L54）。** 8 个头条神经性状同时出现在两个家族里，且在两家族中会得到不同的 BH q 值；头条声明（"all FDR < 0.05 across 32 tests"，L126）使用的是较小的 32 家族。q 值对家族规模单调，先验地把"自己关心的假设"划进小家族、把"对照"扩成大家族，正是 BH 使用中最常见的宽松化手法。我接受两组假设在科学上可区分，但同一批 p 值不应在两个家族中各校正一次再择小者报告——应声明一个主家族（建议 48，即全部 12 个性状一体校正，或在 Supp 中同时给出两种 q 并说明头条数字来自哪一个）。
2. **LOO 变体的家族又不同**：Supp L43 用"8 traits × 4 classes"（Tier 2 LOO）与"12-trait × 4-class"（阴性对照 LOO），即同一套检验在 full/LOO-brain/LOO-tau 三种分类下各形成一个家族，三个家族间不再校正。变体间分类一致率 0.84–1.00（Fig 4a 图注），检验高度相关，Holm/BH 在其间均有效，但"每个变体自成一家族"的切法应明确写出并论证。
3. **硬性矛盾：Table S5 被描述为"12 traits × 3 classification variants × 3 gene classes; 108 rows"（L50）**——12×3×3=108 自洽；但主文与 Supp 的所有家族声明都是 **4 classes**（8×4=32、12×4=48）。3 类与 4 类不可能同时成立：究竟 GD(relaxed) 与 dual 是否进入了 GWAS 家族？若 S5 只含 3 类而家族按 4 类校正，读者无法复算任何 q 值。
4. **Meff = 3（Li & Ji 2005）被引用却从未被使用**（L50）。且用 trait 基因集的 mean Jaccard = 0.148 推 Meff 不是 Li & Ji 的原程序——该法要求检验统计量（或标记）的相关矩阵特征值，Jaccard 重叠只是代理。若 Meff 不进入任何校正，请删除或降为描述；若要保留，请给出 12 性状检验相关矩阵与特征值分解（Supp 一张表即可）。

**要求：** 提供一张"检验普查表"（建议放 Supp）：全文每一个 p 值、所属家族、家族规模、校正方法、单/双侧。这同时解决 M4 与审稿人/读者复算 q 值的全部困难。

### M3. LOO 只做单分量移除；表达轴（tau + brain = RDS 权重的 45%）的联合移除缺失，而这恰是衰减信号指向的方向

LOO 设计的论证逻辑是"移除可能与验证证据共享信息的分量"（L40, L74）。但 tau（GTEx 30 组织特异性）与 BrainSpan 脑区特异性是同一信息轴的两个拷贝：SynGO 术语按突触功能定义，天然富集脑表达基因；神经精神 GWAS 基因同样脑表达富集（L126 已承认）。单分量移除时另一个表达分量仍在分数里——LOO-brain 下 tau 还在、LOO-tau 下 brain 还在，循环通道只被关了一半。稿件自己的数字已经显示衰减梯度：GWAS 8/8 → 6/8（LOO-brain）→ 5/8（LOO-tau）；SynGO 11 → 7 → 5（L114, L126）。**按此梯度，联合移除是最可能改变结论的那个变体，而它恰恰是唯一没做的。**

**要求：** 增加 LOO-expression 变体（caMPRA/nc 两分量按既有的圆整权重方案重整），重跑 Tier 1（SynGO/GO-BP）与 Tier 2（12 traits）全套检验，纳入 Table 3 的 LOO 矩阵与 Holm 家族。若结论存活，循环性论证才算闭合；若大幅衰减，Discussion L158 的"partial circularity … bounded but not eliminated"需要相应改写。

另外一处应明示的嵌套：hCONDEL 按定义位于**高保守非编码区**（McLean et al. 2011 的筛选逻辑），而 nc 分歧分量正是奖励基因体非编码保守——即 hCONDEL 验证与 nc 分量存在与 HAR–caMPRA 同构的潜在嵌套。LOO-nc 下 hCONDEL OR = 4.67（Table 3）实际上已经回答了这个问题，但正文只讨论了 HAR–caMPRA 嵌套（L92），未讨论 hCONDEL–nc 嵌套。请在 L76 或 L150 处显式指出 LOO-nc 行就是对这一嵌套的检验。

### M4. "不校正"声明的范围越界：Tier 4 的不校正约定不能覆盖同节中的非零声明

L54 为 Tier 4 子集检验不校正给出了可接受的辩护（预设假设、零为所关心声明、校正会膨胀 II 型错误）。但同一分析块（L124）里混入了多个**方向为非零**的声明，它们不享受该辩护：

- Spearman ρ = −0.06, P = 7.5 × 10⁻⁶（"two axes largely disjoint"的证据之一）；
- caMPRA 阳性基因的 coding-driven 耗竭 OR = 0.26, P = 2.7 × 10⁻⁵；
- GD 类 hCONDEL 耗竭 OR = 0.59, two-sided P = 0.037（Fig 4f）——P = 0.037 在任何家族规模 ≥2 的校正下都不存活，而它被画进了主图；
- Fig 7a 中 RD 在三个子集的富集（P = 3.8 × 10⁻⁵ to 1.7 × 10⁻²⁵）。

此外三个子集（HAR / hCONDEL / combined = 两者之并）天然不独立，即便不校正也应声明相关性。这些 p 值大多强到校正后仍显著，所以结论大概率不变——但"claim of interest is the null"的豁免条款必须只适用于零声明，非零声明要么并入声明的家族、要么明确标注 nominal。

### M5. GDS 的名义权重与实际信息量脱节：同一 LRT 的两个统计量占 0.65 权重，Selectome 分量近乎常数

GDS = BUSTED 显著性（0.35）+ LRT（0.30）+ RELAX（0.20）+ Selectome（0.15）（L36）。问题：

1. 前两个分量是**同一次 LRT 检验的两个函数**（显著性由 LRT 的 p 值决定），残差化只去除了长度，不去除两者间近乎确定性的相关；GDS 的有效权重结构因此是"一个 BUSTED 检验 ≈ 0.65 + RELAX 0.20 + Selectome 0.15"，与 L36 所述"balance substitution-based and expression/functional evidence"的设计意图不符——GDS 内部根本没有第二种独立证据可平衡。
2. Selectome 分量在 4,974 基因中仅 24 个非零（0.5%），对 99.5% 的基因是常数 0，名义权重 0.15 实际近乎不产生区分度。
3. 权重敏感性分析（L44, L96；±50% 单权重、1,000 次联合扰动、seed 42）扰动的是**名义**权重；在分量强相关下，名义权重的扰动范围与有效权重的扰动范围并不对应。

我并不要求重做分类框架（这不现实也不必要），但要求：(a) 在 Methods 或 Limitations 中显式承认 BUSTED-p 与 LRT 分量的共线性与 Selectome 分量的退化性，给出"有效权重"的一句话刻画（例如各分量与 GDS 的相关或 Shapley 式分解——L70 的 dominant-component 分解已有雏形，GD 类 67% 由 p 分量主导、30% 由 LRT 主导，恰恰印证了这一点）；(b) 把"weights were fixed to balance … evidence"（L36）的措辞软化为名义权重。

---

## 4. Minor comments

1. **可复现性——置换种子缺失。** L44 给出了权重扰动的 seed 42，但 L48 的 1,000 次 gene-label 置换与 L52 的 downsampling（100 replicates × 11 levels）未给种子。请在 Statistical conventions 一并列出。
2. **置换深度不足以支撑所引显著性。** 置换 floor = 1/1,000（L54，无伪计数，可），但 L114 称最强 10 个术语"retain P ≤ 0.002 under 1,000 permutations"——对超几何 q 低至 2.1 × 10⁻⁴ 的术语，1,000 次置换无法分辨 p < 10⁻³，置换验证对这些术语只提供了"p ≤ 0.002"的弱确认。要么对 top 术语加大置换次数（≥10⁵），要么把措辞改为"permutation confirms P ≤ 0.002, consistent with but not resolving the parametric tail"。
3. **χ² 自由度未报。** L70 与 Supp L41："59% vs 48%; χ² = 12.3, P = 0.006"。df = 1 时 χ² = 12.3 对应 P ≈ 4.5 × 10⁻⁴；P = 0.006 对应 df = 3（4 分量 × 2 组）。引文只给了单分量的 59/48 对比，读者会按 df = 1 复算而得到不一致的 p。请写明检验构造与 df。
4. **RELAX 计数链有两处歧义（L64）。** (a) 117 个 FDR 显著 relaxed，但只有 52 个被重分类为"gene-driven (relaxed)"——其余 65 个的去向（推测：不满足 GDS 边际而未入 GD 类）未说明；(b) "seven FDR-significant genes with degenerate K = 0 … are excluded"——这 7 个是在 117 之内还是之外？两个数字各一句话即可澄清。
5. **Limitation (8) 的 OR 方向与算术（L176）。** "40.4% versus 32.7% of strict-call genes; OR = 0.72"——0.72 是 strict:permissive 方向（0.327/0.673 ÷ 0.404/0.596 ≈ 0.717），与句子的叙述顺序相反；"22.0% versus 9.2%; OR = 1.52"无法由所述 2×2 复算（所述比例给出 OR ≈ 2.8），请说明该 OR 的对照组到底是什么。
6. **Shao 检验的单双侧未声明（L102）。** "OR = 0.55, P = 0.21"——方向是耗竭，按 L54 的约定（greater 单侧）该 p 无法解释；是耗竭方向单侧还是双侧？Table S5 若已有 p_two_sided 列，请在正文指明所用列。
7. **LOO 变体间一致率 0.84–1.00 的分母被 unclassified 多数主导。** LOO-caMPRA 的 RD 从 293 膨胀到 800、LOO-brain 到 613（Table 3），整体一致率仍有 0.84，说明该指标对 RD 类成员的稳定性不敏感。请在 Fig 4a 图注或 Supp 补每类一致率（尤其 RD 的 per-class agreement），否则"classification agreement 0.84–1.00"会让读者高估 LOO 变体间的相似性，从而高估"every variant survives"的独立复制含义。
8. **Table 3 与 Fig 4f 的 GD 分母不一致。** Table 3 注明 GD 排除 relaxed（n = 1,214），Fig 4f 的耗竭检验包含 relaxed（n = 1,266）。两处各有道理，但请在 Table 3 脚注交叉指明。
9. **"4.7× genomic median"（Abstract L12）与 133.6 vs 28.6 kb（L126）的对照集不同**——后者是"neither class"基因（n = 3,415）而非全基因组中位数。数字接近，措辞请统一。
10. **attrition 与长度机制的潜在交互。** Supp L7 承认 BUSTED 失败率随比对长度上升（6,840 → 4,974，27% 丢失），即宇宙本身偏向较短基因；而 RD 类又偏向长基因（M1）。宇宙的长度截断如何影响 RD 类组成与各类富集的分母，值得在 Limitations 加一句讨论（不要求重跑）。
11. **中性模拟的检验力定位（L176）。** 400 个 ω = 1 模拟校准了 size（1.8% < 5%，BH FDR 0/400），但没有 ω > 1 的尖峰，不能校准混合模型下的 FDR；0/400 的 95% 上界约 0.75%，请给出区间而非点值。另：1.8% 的保守性部分可能来自 evolver 与真实比对的差异（如二核苷酸上下文），宜加半句限定。
12. **k-means 回收率的事后选择偏倚（L90）。** "best cluster captured 87.4% of 293 RD"是 4 个簇中择优；纯度 ~10% 意味着该簇约 2,500 基因。请补 Adjusted Rand Index 或 normalized mutual information（对簇选择不敏感），或报告每个监督类在各簇的完整分布（Fig 3f 已有雏形，正文化即可）。
13. **hCONDEL-only 子集的检验力（L76）。** 9/90, OR = 1.80, P = 0.081 被正确解读为功效不足；请补一句该样本量下的最小可检 OR（MDE at 80% power），把"非显著"量化而非定性带过。
14. **阴性对照家族只有 4 个性状（L126）。** Fisher P = 0.091（8/8 vs 2/4）已诚实标注受限于家族规模；建议未来扩展对照家族（如 autoimmune、metabolic、anthropometric 各 2–3 个）以提高该比较的鉴别力——可在 Discussion 作为 future work 一句带过。
15. **编辑性占位符。** Data availability 的 Zenodo DOI 仍为"to be assigned"（L188）；Author contributions 保留 AI 工具披露的条件句占位（L194）——投稿前必须落实。

---

## 5. 什么改动能让我提高分数

**完成 M1 的长度匹配/长度校正 SynGO 重检并如实报告结果（无论存活与否），加上 M2 的"全文检验普查表"——仅这两项，总体分即从 7 升至 9，推荐转为 Accept。**

---

## 附：本轮抽查通过的内部一致性（供编辑参考）

以下为我独立复算无误的项目：Selectome 2×2（7/1,266 vs 17/3,708 → OR = 1.21 ✓）；hCONDEL–HAR 重叠 φ = 0.052 ✓ 与 caMPRA–HAR 嵌套 φ = 0.434 ✓；HAR 富集 81/293 vs 395/4,681 → OR = 4.15 ✓；Shao 条件分析 4/18 vs 期望 6.1 → OR = 0.55 ✓；Holm 校正 0.032 × 5 = 0.16 ✓；分类计数 1,214/293/11 与百分比 24.4/5.9/0.2% ✓；RELAX 604/1,660 = 36.4% ✓；覆盖外推线性模型（截距 ≈ 261、slope 0.33）对 2,914 交叉点、476 → 417、4,974 → ≈1,900 三个锚点均自洽 ✓；GO BP 名义期望 2,507 × 0.05 ≈ 125 ✓；LOO 五变体 × 2 证据的 Holm 家族即便合并为 10 检验，最大 p（5.3 × 10⁻⁴ × 10）仍 < 0.05 ✓（LOO 结论对家族切法稳健）。稿件的算术纪律罕见地好；本评审的分歧全部在设计与披露层面，不在计算层面。
