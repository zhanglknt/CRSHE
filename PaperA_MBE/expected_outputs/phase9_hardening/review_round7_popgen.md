# Review Round 7 — Reviewer: Population & Quantitative Genetics (MBE)

**Manuscript:** Coding Selection Dominates by Prevalence, Regulatory Selection by Specificity: A Comparability-Corrected 10-Primate Analysis (v9)

**Reviewer expertise:** GWAS architecture, polygenic traits, gene-mapping conventions, confounding control, power of selection tests.

---

## 1. 总评与推荐决定

**一句话总评：** 这是一篇自我纠偏力度罕见、统计披露近乎范本的稿件——作者已主动承认 GWAS anchoring 的幅度由基因长度驱动并降级了"specificity"叙事——但**全文唯一留存的功能特异性卖点（"the synaptic gene-set enrichment as the length-robust functional signal"）本身从未接受过作者自己建立的基因长度检验**，加之 HAR/hCONDEL 邻近富集同样未做长度匹配零假设，使标题级结论目前建立在一个未检验的断言上。

**推荐决定：Major Revision**（所有要求的分析均可用现有数据完成，无需新实验或新数据；若长度匹配检验支持作者断言，则可迅速转为 Minor）。

## 2. 评分表（1–10 整数）

| 维度 | 评分 | 简评 |
|---|---|---|
| 新颖性 | 7 | 组件均非首创，但 comparability-correction + 覆盖不对称量化 + 对己不利结果的完整披露构成方法学示范 |
| 方法严谨性 | 7 | 分类/验证侧接近示范级；长度混杂处理层次齐全（logistic/CMH/双映射）但止步于 GWAS，未延伸到其余下游富集 |
| 统计分析 | 7 | 多重检验族定义、单边约定、置换校准均规范；但存在族定义前后不一致与一处与 S5 表矛盾的"non-significant depletion"表述 |
| 结果可信度 | 6 | prevalence 与 disjointness 结论可信；"specificity"的可信核心（SynGO）未做长度匹配零检验，目前是不可证伪状态 |
| 写作清晰度 | 8 | 信息密度极高而结构清晰；Limitations 与 Statistical conventions 堪称范本 |
| 总体 | 6 | 距离可发表只差一个关键零假设检验与若干一致性修复 |

---

## 3. Major comments

**M1. "Length-robust functional signal" 是断言而非结果：SynGO 富集从未做过基因长度校正。**
摘要（第 12 行）与 Discussion（第 166 行）均写道：长度校正消灭了全部 12 个性状的 GWAS 富集后，"leaving the synaptic gene-set enrichment as the length-robust functional signal" / "The length-robust core of the regulatory mode's specificity is the synaptic gene-set enrichment"。但通读全文，SynGO/GO-BP 富集（第 114 行；Fig. 7b）只做过两类稳健性检验：universe 背景超几何检验与 LOO-brain/LOO-tau 重分类——**两者都没有移除长度轴**。稿件自己的数据恰恰预示该轴是主要混杂：(i) RD 基因中位跨度 133.6 kb vs 28.6 kb（4.7×，第 126 行）；(ii) 驱动 52% RD 成员的非编码保守分歧分量与基因跨度 Spearman ρ = 0.47（第 126 行）——该分量在 LOO-brain/LOO-tau 下仍然保留，故"LOO 下 SynGO 仍有 7/5 个术语显著"完全不能排除长度机制；(iii) 突触基因（neurexin/cadherin/CNTN 家族等，恰好也是 RD 榜首 CLSTN2、CDH10、DAB1 的同类）是全基因组最长的基因类别之一。换言之，SynGO 富集与 GWAS 富集暴露于**完全相同**的混杂结构，作者却只对后者做了校正、对前者直接冠以 "length-robust"。这是当前稿件最严重的证据缺口，因为标题的 "by Specificity" 与摘要的最终落点都压在这一断言上。**要求**：(a) 对每个 FDR 显著的 SynGO/GO-BP 术语，在长度匹配零假设下重检（如从 universe 中按 log₁₀ span 分层抽取与 RD 等大的随机基因集 1,000 次，报告富集的置换 p 与效应量衰减）；或 (b) 以 logistic 回归（术语成员 ~ RD 状态 + log₁₀ span）报告校正后 OR；(c) 若衰减后仍显著，"length-robust" 成立并可保留标题；若不再显著，标题与摘要必须相应降级。

**M2. 主要独立验证（hCONDEL/HAR 邻近富集）同样未控制基因长度的机械捕获效应。**
第 76 行将 hCONDEL 富集（OR = 3.19, P = 7.8 × 10⁻⁵）定为 "primary independent validation"，Fig. 4f 同。但 HAR/hCONDEL 映射采用 ±50 kb 邻近窗口（Limitation 10）：基因越长，其 ±50 kb 邻域覆盖的基因组越大，被 HAR/hCONDEL 命中的先验概率越高。RD 基因恰长 4.7×，因此即便两类元件与 RD 生物学毫无关系，也会机械性产生 OR > 1 的富集。LOO 矩阵（Table 3）只能排除评分分量与被验证据之间的循环，不能排除这一几何混杂；LOO-nc 变体（移除与长度相关的 nc 分量后 HAR OR 反而升至 9.39）提示信号不全是长度，但也可能由保留的 caMPRA 分量（与 HAR 完全嵌套，第 92 行自认 100% nesting）解释。**要求**：(a) 报告 Table 3 五个变体各自 RD 类的中位基因跨度，使读者可判断各变体长度混杂的相对大小；(b) 对 hCONDEL（核心验证）与 LOO-caMPRA all-HAR（清洁估计）做长度分层或长度匹配零检验（长度三分层 CMH 或分层置换即可，与 GWAS 侧已做者同法）；(c) 若 hCONDEL 富集在长度匹配后仍显著，本稿的证据链将真正闭环——我预期这是最可能出现的结果，但必须先做出来。

**M3. 长度混杂的处理层次齐全，但对"幸存信号"的解释框架需要统一：连续校正全灭 vs 粗分层部分幸存，应以后者更弱而非更强来解读。**
第 126 行报告了三层校正：logistic log-length（12/12 性状全灭，adjusted OR 0.33–1.64，全 q > 0.44；族间差异 p = 0.24）、CMH 长度三分层（EA 1.87、SCZ 1.95 幸存，height/T2D 边缘）、双映射约定（选择性幸存）。作者正确地注明 CMH "removes only between-stratum confounding"。但随后 "the anchoring is therefore directionally robust across mapping conventions and coarse stratification" 的措辞，以及 Discussion（第 166 行）"educational attainment and schizophrenia survive coarse length stratification and alternative mapping conventions"，仍把粗校正下的幸存者当作正面证据陈述。逻辑上应反过来：当连续协变量校正使所有效应归零（且调整后点估计 11/12 性状仍在 1.02–1.64 区间、方向未反转），三分层内的幸存者最 parsimonious 的解释是**层内残余混杂**——三分层的层内跨度仍相差数倍，而 ρ = 0.47 的剂量关系在层内继续起作用。**要求**：(a) 在正文明确写出这一解释层级（连续校正 > 分层校正 > 未校正），并将 EA/SCZ 的 CMH 幸存表述为 "consistent with residual within-stratum confounding or a modest residual effect, which the present design cannot distinguish"；(b) 顺手可做增强分析：长度十分层 CMH 或限制性立方样条 logistic，观察 EA/SCZ 是否随分层细化而进一步衰减——若衰减，即证实残余混杂解释。

**M4. 映射约定敏感性做了，但"author-reported 幸存"恰是循环性最强之处，需要点破。**
第 126 行：author-reported 映射下仅 EA（OR 2.99）与 SCZ（OR 1.81）幸存。作者将此列为稳健性证据。然而神经精神 GWAS 的作者报告基因绝大多数是通过脑 eQTL 共定位、脑表达优先或脑功能注释提名的——即 author-reported 约定系统性地把"脑表达轴"重新引入基因集，而 RDS 的 45% 权重是表达分量（第 126 行第一点自认）。于是"author-reported 下 EA/SCZ 幸存"与"RD 含表达分量"之间存在结构性共谋，其证据价值低于 locus-dedup 约定（位置映射，与表达注释无关）。**要求**：在正文点破这一不对称，说明 author-reported 约定的幸存不能作为与表达轴无关的证据；locus-dedup 约定下幸存的 EA/SCZ/intelligence 反而是该约定自身也含位置-长度偏倚的未校正结果，二者均不构成独立于长度与表达轴的锚定。

**M5. 正文与 S5 表存在一处数字级矛盾：major depression 的 GD 耗竭并非 "non-significant"。**
第 126 行写道 GD 类 "none significant, with non-significant depletion for major depression and bipolar disorder"。但 Supplementary Table S5（full × Major depression × gene_driven_strict 行）给出 OR = 0.563（95% CI 0.36–0.87），p_two_sided = 0.0098，q_bh = 0.035——按 S5 自己声明的 "conservative reference"（双侧）该耗竭是 FDR 显著的（bipolar 的耗竭 p = 0.226，确不显著）。我理解正文采用单边富集方向约定（第 54 行），在此约定下 GD 行 q = 1.0；但 "non-significant depletion" 这一短语在双侧尺度上是错的，且一个有 FDR 显著保护的效应（GD 对 MDD 耗竭近半）本身是信息——它提示两类基因在性状架构上的对立比正文所述更强。**要求**：统一口径——要么正文改为 "gene-driven genes show no enrichment for any trait; major depression shows a significant depletion under the two-sided reference (OR = 0.56, q = 0.035)"，要么在 S5 中同时给出单边 q 使两处的显著性标签可比。

**M6. 多重检验族定义在三处互相矛盾，需统一并写清实际使用的族。**
第 50 行（Methods）称 "Benjamini–Hochberg corrected across all 32 trait × class tests"；第 54 行称 "Tier 2: all 32 trait × class tests; … GWAS negative controls: all 48 trait × class tests per variant"；而 S5 表头称 "q_bh is Benjamini-Hochberg across the 36 tests within each classification variant"（12 性状 × 3 类 = 36）。32、48、36 三个族大小互不相同，读者无法判断 S5 中 q 值的实际分母，也无法复算（我用 36 检验重算 full 变体可复现表内 q，说明实际执行的是 36/变体；那么第 50/54 行的 "32" 与 "48" 就是过时文本）。此外作者已计算有效检验数为 3（Li & Ji，第 50 行）却未使用——既然 BH 在正相关下仍有效但保守，至少应说明为何保留名义族。**要求**：统一 Methods 与 S5 表头的族定义；一句话说明不采用 effective-tests 校正的理由。

**M7. 锚定分析还能支撑什么？——我的裁决与建议的最终表述。**
回答"该锚定分析还能支撑摘要中的任何 claim 吗"：能支撑的有且仅有——"RD 基因比 GD 基因更频繁地落在神经精神/认知 GWAS 位点中（未校正 OR 1.9–4.0，方向在全部映射约定与 LOO 变体下一致），且该关联的幅度主要由基因长度解释（校正后全 null，族间差异 p = 0.24）"。它**不能**支撑：神经特异性（height/T2D 同样富集，长度校正前族间差异在 1,000 次置换下 p = 0.009 但校正后消失）；历史选择靶点（作者已自认 present-day variation caveat）；以及任何独立于表达轴的锚定（M4）。摘要第 12 行目前的表述（"directionally robust but gene-length sensitive"）基本诚实，但最后半句把落点交给了未检验的 SynGO 断言（M1）。因此：**标题的 "by Specificity" 当前由 M1/M2 两项未做的检验支撑**。这也是我把决定定在 Major 而非 Minor 的唯一原因——这不是润色问题，而是结论的承重墙尚未验收。

---

## 4. Minor comments

**m1.** 阴性对照的 null 部分按功效排序：Crohn（95 基因）与 LDL（256）是 4 个对照中最小的两个 trait set，height（2,191，占 universe 44%）与 T2D（561）是最大的两个（S5 n_trait_universe 列）。"2/4 对照富集"与"按功效排序的检出"无法区分；Crohn 的 null（OR 0.52，CI 0.16–1.64）尤其可能功效不足。建议在讨论 2/4 模式时明确注明这一点。

**m2.** height 的 trait set 覆盖 universe 的 44%（2,191/4,974）。当"阴性对照"的映射基因接近半数基因组时，Fisher 检验实际测量的是"基因是否容易被 GWAS Catalog 映射"，这正是长度检测器——作者事实上也如此使用它，但建议在 Methods 一句话点明，避免读者误读为生物学对照。

**m3.** 8 个神经性状并非 8 次独立观察：EA/intelligence/cognitive performance 很大程度共享 SSGAC 系研究，SCZ/BIP/MDD/neuroticism 高度遗传相关。作者已报告 mean Jaccard 0.148 与有效检验数 3，但建议在 S5 或补充文本中给出 8×8 的逐对重叠矩阵（或至少最大逐对 Jaccard），使 "8/8" 的折扣可核查。

**m4.** 族间 z = 3.49 的 pooled-OR 比较对 8 个相关性状的对数 OR 做了加权合并，性状间协方差未入模型，SE 偏乐观；作者同时给的 trait-label 置换 p = 0.009 才是可辩护的那个。建议正文把置换 p 作为主结果、z 检验降为参考（或注明 z 假设独立性）。

**m5.** S5 表头将长度校准结果指向 "gwas_confound_calibration.json"。JSON 不是期刊补充材料的适当载体；请把 logistic 调整后 OR/CI/q、CMH 三分层结果、双映射结果作为 S5 的额外表或新增 S6 sheet 落入 xlsx。

**m6.** 单边检验约定（第 54 行 "one-sided in the enrichment direction"）使所有 GD 行的 q = 1.0，掩盖了 MDD 的显著耗竭（见 M5）。建议在 S5 增加 depletion 方向的 q 列，或直接统一为双侧（作者已把双侧称为 "conservative reference"，那就让它成为主口径）。

**m7.** MAPPED_GENE 分词约定对含连字符/撤回符号（如 "C4orf48"-style 旧符号、antisense、"GENE1 - GENE2" 区间记法）的处理规则未说明；多基因位点逐一计数对基因密集区（MHC 等）的膨胀未量化。建议在补充方法中给出分词正则与清洗规则，并报告去除 MHC 区（chr6:25–34 Mb）后核心性状 OR 的变化（一行即可）。

**m8.** GWAS Catalog 全量 115 万行混合了不同年代、样本量与祖源的 study；12 个性状均未说明祖源构成（推测以 EUR 为主）。选择事件早于出非洲，以 EUR 为主的当下变异架构做锚定的外推限制值得在 Limitation (9) 中加半句。

**m9.** 摘要中 "odds ratios 1.9–4.0" 与 "gene-length correction abolishes the trait enrichment" 同句并置，信息上正确但修辞上前者会先声夺人；建议把未校正 OR 区间后紧跟 "(height and type 2 diabetes show similar enrichment)"——其实作者已这样做，只是在摘要中的顺序仍可再压实一句，把 "directionally robust, magnitude length-driven" 的八字结论直接写进摘要末句。

**m10.** 第 50 行 "trait sets span from 103 (regulation-driven ∩ educational attainment) to 822"——103 是交集数、822 是 trait set 大小，两个不同量纲被并列在一个 "span from … to …" 中，请改写。

**m11.** Limitation (8) 报告 permissive-call 基因 coding-selection 率更高（40.4% vs 32.7%）且 strict 子集结论不变，处理得当；但 strict 子集分析同样未做长度校正的 GWAS 重检（intelligence 在 strict 子集已掉到 P = 0.050），建议与 M1–M3 的修订合并报告。

---

## 5. 什么改动能让我提高分数

**一句话：用作者自己已有的数据做出长度匹配（或长度协变量）零假设下的 SynGO 与 hCONDEL/HAR 富集重检——若幸存则在正文坐实 "length-robust"、若衰减则相应降级标题与摘要——同时修复 M5/M6 两处与 S5 表的数字级矛盾，我的总体评分将从 6 升至 8，并转为 Minor Revision 推荐。**
