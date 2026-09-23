# Review Round 8 — Reviewer: Statistical Genomics & Reproducibility Audit (MBE)

**Manuscript:** Coding Selection Dominates by Prevalence, Regulatory Selection by Specificity: A Comparability-Corrected 10-Primate Analysis (v9)

**Reviewer expertise:** multiple testing / test-family architecture, circularity & information leakage, effect size & power, extrapolation validity, reproducibility.

**声明：** 本轮为全新盲审，未参考既往任何轮次意见。行号指 `manuscript_english_v9.md`；"Supp" 指 `Supplementary_Text_v9.md`；表号指 `Main_Tables_v9.md` 与 `TableS5_anchoring_master.csv`。所有注明"复算"的数字均由我从原始计数/表中独立重算。

---

## 1. 总评与推荐决定

**一句话总评：** 这是我近年审到的在多重检验纪律与自我证伪上最严格的稿件之一——我对 Table S5 全部 192 行的 BH q 值做了逐行复算（四个变体 × 48 检验，单双侧两族），与文件值误差在机器精度内；OR、φ、p 值抽查近二十处全部吻合；两处经验 null 校准设计得当。问题集中在三处**推断的呈现方式**而非计算错误：(i) srv 外推的 28.0%/1,394 这一进入 Results 层面的点估计没有任何不确定度量化，且其相对尺度移植假设未被声明或给出替代口径；(ii) hCONDEL-only 子集的"93% post-hoc power"论证框架不当（数字本身我可复现，但有严格正确的直接检验可用）；(iii) 非单调的 LOO 衰减链 8/8→6/8→5/8→6/8 未获解释，类大小（功效）与分量移除在链中混杂。另有一簇小的数字/表述错误需要勘正。均无需新数据。

**推荐决定：Minor Revision**（全部意见可由文字修订 + 已有管线的分钟级重算解决；无任何一条动摇类级结论）。

---

## 2. 评分表

| 维度 | 评分 | 简评 |
|---|---|---|
| 新颖性 | 8 | 单组件均成熟，但"有界支撑 BH 论证 + 双 null 校准 + LOO 全矩阵 + 覆盖反事实"的组合审计对比较基因组学有方法示范价值。 |
| 方法严谨性 | 8 | 分类/验证分离、LOO 逐分量移除、长度混杂的自我拆解均属范本级；扣分在 GDS 名义权重（0.65 实为同一 LRT 的两个单调变换）与外推链。 |
| 统计分析 | 7 | 计算精度满分（见 §5 复算记录）；扣分在 srv 外推无 CI、post-hoc power 框架、非单调衰减链未解释、Tier 4 校正口径前后不一。 |
| 结果可信度 | 8 | 类级富集信号经多重独立路线收敛；作者对自己 GWAS 锚定的长度校正证伪诚实且彻底。个位数 RD/dual 归属的权重依赖性已如实披露。 |
| 写作清晰度 | 7 | Statistical conventions（L54）是教科书式一节；但多处括号内并置不同参照系的统计量（见 m4/m7），两处"parity"数字（59% vs 43%）易混淆。 |
| **总体** | **7.5** | 修好后 8.5–9；当前被 §3 的三项呈现级问题压住。 |

---

## 3. Major comments（三条，均为呈现/推断框架问题，非计算错误）

### M1. srv 外推的 28.0%（≈1,394/4,974）是一个没有不确定度的点估计，且移植假设未声明

L60/L142/L176：400 基因分层样本（200 显著 + 200 非显著，按比对长度与 LRT 三分位平衡）在 393 个可比基因上测得 49.4%→40.7%（−8.7pp，相对 −17.5%），再按**相对口径**外推到全宇宙 34.0%×(1−0.175)≈28.0%。我复算：34/194 = 17.53%，0.33977×(1−0.1753) = 28.02%，0.280×4,974 ≈ 1,394——算术全部吻合。问题在推断层：

1. **无置信区间。** 393 个基因的二项波动本身即 ±4.9pp（40.7% ± 1.96×2.5%），传导后 28.0% 的区间大约为 25–31%。该数字出现在 Results 首节并支撑"the rate should be interpreted as an upper bound... point estimate at 28.0%"的表述，必须带区间（对 393 个基因做 bootstrap 或对翻转计数做二项区间即可）。
2. **移植口径未声明、未给替代。** 相对口径隐含"翻转率在各显著性层内恒定"——等价于翻转全部发生在显著层（34/194），这在抽样设计下恰好可传输，是合理的；但**绝对口径**（−8.7pp 直接减）给出 25.3%。两个口径相差 2.7pp，读者有权看到这个括号。更规范的做法是说明分层抽样比例并给 IPW 估计——作者已有每层抽样分数，成本为零。
3. **双向翻转未披露。** CpG 与拓扑重跑给了 retention（93.8%/92.1%），唯独 srv 只有净值。若 --srv Yes 下有非显著基因翻正，净 −8.7pp 背后的总搅动更大，外推解释会变。请报告 srv 的 retention 与反向翻转数。
4. **相关的不对称：** codon-aware 重比对 +14.8pp（64.1%）被解释为"alignment noise attenuated power rather than inflating false positives"（L142），但该管线变体**没有匹配的中性模拟校准**（evolver/bootstrap 均基于生产管线的 gap mask）。去掉 59.1% 的列后 LRT 分布改变，方向性解释（偏向"产能管线可信"的一侧）缺一条腿。要么补一句明确的不确定性声明，要么对该变体跑同一 evolver null。

### M2. "93% post-hoc power"论证框架不当——且存在一个严格更强、成本为零的直接检验

L76：hCONDEL-only 子集（n=90）富集不显著（9/90, OR=1.80, P=0.081），作者以"93% post-hoc power to detect the full-set effect at this sample size"论证"attenuation is informative rather than a power artifact — the enrichment concentrates in genes supported by both evidence types"。

先说公道话：这个 93% **本身可复现**。我用 Monte Carlo 精确 Fisher（单侧，α=0.05，备择取全集 OR=3.19 对应的 p₁=0.166 vs p₀=0.059，n=90 vs 4,884）得到 power = 0.930，与稿件一致；且以预设效应量（而非观测效应量）做功效计算规避了 Hoenig–Heisey 陷阱的经典形式。

但逻辑上该句做的推断是"衰减为真"（即子集 OR 不同于全集 OR），而功效计算回答的是另一个问题（若全集效应为真，子集能否检出）。支持其结论的正确工具是**直接比较两个子集**：双重支持 9/22 (41%) vs hCONDEL-only 9/90 (10%)，我复算 Fisher OR = 6.23，单侧 P = 1.4 × 10⁻³（双侧同值）——**显著**，且这比功效论证更强、更直接。请用此检验替换（或并列）post-hoc power 句式；顺带把 Fig. 4c 图注中"directionally consistent but individually non-significant"升级为这一直接比较。这也顺带修复了"attenuation is informative"目前依赖读者接受功效论证的软肋。

### M3. 非单调衰减链 8/8→6/8→5/8→6/8 未获解释：分量移除与功效在链中混杂

我从 Table S5 逐格复核了这条链（单侧 48 族 q<0.05 的 RD-神经性状计数）：full 8/8；LOO-brain 6/8（ASD q=0.19、intelligence q=0.060 出局）；LOO-tau 5/8（再失 bipolar q=0.43）；joint LOO-brain-tau 6/8（bipolar 以 q=3.4×10⁻⁴ **回归**）。稿件只报告计数，未指出链条非单调，更未解释。

机制其实清楚：LOO-tau 把 RD 类缩到 148 基因，joint 反而扩到 403——**类大小（即功效）在变体间剧烈变化（110–800，Table 3），衰减链把"移除了哪个分量"与"检验功效变了多少"混在一起**。bipolar 在移除更多分量后回归，正是功效回升（148→403）而非证据增强。这不是错误，但作为依赖该链做"部分循环性量化"立论（L158：SynGO 7/5/7、GWAS 6/5/6）的稿件，必须显式声明非单调性及其功效来源，否则读者会把 5/8→6/8 误读为证据的非单调。**建议：**(i) 一句话声明链条非单调且由类大小驱动；(ii) 可选但强烈建议——把 joint 类下采样至 148 基因（100 次重复）重跑 8 个性状，给出去除功效混杂后的纯组成效应。成本：分钟级。

---

## 4. Minor comments

**m1. Tier 4 校正口径在 Methods 与 Results 之间矛盾。** L54（Statistical conventions）："The subset tests in the empirical-disjointness analysis (Tier 4) are reported **uncorrected**"；L124 却报告"both this depletion and the correlation survive Benjamini–Hochberg correction across the six tests of this paragraph (q < 0.02)"。实际情况（我梳理）是：三个 null-claim 子集检验不校正（合理，claim 是 null），三个正 claim（ρ=−0.06、caMPRA 耗竭、GD-hCONDEL 耗竭）在六检验族内 BH——设计本身自洽且偏保守，但 Methods 句子按字面是错的，请改写为明确区分两个亚族。

**m2. 一簇数字/表述勘误（均不影响结论）：**
- (a) L126："not under LOO reclassification (two-sided q = **0.11–0.15**)"。Table S5 实测：LOO-brain 0.143、LOO-tau 0.148、joint 0.081——三个变体的正确范围是 **0.081–0.148**；即便只取单分量变体也是 0.14–0.15。"0.11–0.15"两头都对不上。
- (b) L50："trait sets span from **103** (regulation-driven ∩ educational attainment) to 822"。103 是 RD∩EA 的**重叠数**而非性状集大小；Table S5 中最小性状集是 Crohn 的 95。句子把两个量混为一谈。
- (c) L126："author-reported-gene mapping retains only educational attainment (OR = **2.99**, P = 0.015)"——OR 与主约定（2.987）在三位有效数字上完全相同，而 P 从 8.7×10⁻¹⁶ 变到 0.015。映射约定改变几乎必然改变 k/n 从而改变 OR；完全相同的 OR 伴随五个数量级的 P 移动，疑似从主约定行复制粘贴。请核对源表。
- (d) 12.1% vs 12.2%：Abstract 与 L134 同述 604/4,974——前者 12.1%（正确），L134 与 Concluding remarks 作 12.2%（604/4,974 = 12.14%）。统一为 12.1%。
- (e) Limitation (8)："higher permissive-call fraction than the regulation-driven class (**22.0% versus 9.2%**; OR = **1.52**, P = 5.8×10⁻⁷)"。我复算：OR = 1.52（P≈3.4×10⁻⁷）对应的是 GD vs **全部非 GD**；而 22.0% vs 9.2% 是 GD vs RD，其 OR = **2.78**（P≈7×10⁻⁷）。同一括号内并置两个不同参照系的统计量，请拆开或统一参照。
- (f) L86："Among genes with any functional regulatory evidence, 69.5% (41/59)"——59 的定义（BUSTED 非显著 ∩ caMPRA 阳性）全文未出现，读者无法重建；且"any functional regulatory evidence"按字面应含 HAR/hCONDEL，与实际口径（仅 caMPRA）不符。请给出明确定义。
- (g) Limitation (8)："intelligence falling to a marginal OR = 1.69 (P = **0.050**)"——边界值报告请给三位小数（0.050 恰好压线的写法在审稿中必然被问）。
- (h) Fig. 3d 图注："joint weight-**and-threshold** perturbation (1,000 replicates)" vs Methods L44："1,000 joint **weight** perturbations"。二者是否含阈值联合扰动需统一表述。

**m3. "empirical FDR of 0.26% (1/388)"术语不精确。** 核查 `neutral_sim_parametric.json`：400 个模拟中 388 个返回结果（12 个失败未在文中披露），BH 在 α=0.05 下拒绝 1 个。在纯 null 下，该拒绝即假发现，FDP=1；1/388 是**null 下的 BH 拒绝率**（每检验），不是 FDR。按此率折算，真实数据期望假拒绝 ≈ 4,974×0.0026 ≈ 13，隐含 FDR ≈ 13/1,690 ≈ 0.8%——结论（保守）不变，但术语与 12/388 的损耗请写准。同理 L176 evolver 句"the empirical BH false-discovery rate ... is 0 of 400"宜作"empirical rejection rate 0/400"。

**m4. 48 检验族的内部结构值得一句话。** 族内 4 个类中 GD_strict ⊂ GD_all（共享 96% 成员）、dual 仅 11 基因（其 OR 的 CI 跨 0.08–46，近乎无信息），因此 48 的名义族长大于有效检验数——方向保守、且已正确援引 BY2001 处理正相关，但类侧嵌套与性状侧相关不同，请补一句承认。另外 Table S5 的 192 = 4 变体 × 48 是**按变体分别成族**的：跨变体不重校正（作为敏感性分析合理），但这一架构选择目前只能由读者从 q 值反推，请在 Methods 写明。

**m5. 种子披露接近完整，缺两处。** 已披露：20260915（Tier 1 置换）、20260923（长度匹配 null）、20260917（性状标签置换）、42（权重扰动）。缺：覆盖下采样（11 水平 × 100 重复，Fig. 7d）与跨越点 bootstrap CI（2,741–3,146）的种子；k-means 的 100 个种子值（或生成规则）。另：置换 P 下限 1/1,000 无伪计数已声明，很好。

**m6. 两个"parity"数字需要消歧。** Abstract 与 L162 的 ~59%（2,914，GD 计数固定在 1,214）与 L128 的 43%（2,158，两类均随证据增长重估）回答不同问题，但都用"match/equal size"措辞。建议在首次出现处用半句话点明口径差异，否则 59% vs 43% 会被当作内部矛盾。

**m7. 条件率比较的分母。** L86 "69.5% vs 8.9%" 已正确标注"largely reflects scoring, not biology"——很好；但 8.9%（293/3,284）是 RD 的构造性基率，两数之比没有推断内容，建议删去比较只留描述，或明确标注此为机械结果。

---

## 5. 抽查复算记录（全部通过，除 §4 已列勘误）

| 项目 | 稿件值 | 复算值 | 结果 |
|---|---|---|---|
| Table S5 全表 BH q（4 变体 × 48 检验 × 单/双侧，768 个 q 值） | 文件值 | 逐行重算，max\|Δq\| < 1×10⁻¹⁵ | **PASS（精确）** |
| 衰减链计数 | 8/8→6/8→5/8→6/8 | 从 S5 q 值逐格重数 | **PASS** |
| OR（SCZ 2.827、EA 2.987、bipolar 4.005、height 2.215 等） | 表值 | 2×2 重算至 4 位小数 | **PASS** |
| hCONDEL 全集富集 | OR=3.19, P=7.8×10⁻⁵ | k=18：OR=3.194, Fisher P=7.81×10⁻⁵ | **PASS（精确）** |
| post-hoc power | 93% | MC 精确 Fisher（单侧）= 0.930 | **PASS（可复现；框架见 M2）** |
| 双重支持 vs hCONDEL-only | 未报告 | OR=6.23, P=1.4×10⁻³ | 建议补报（M2） |
| caMPRA 编码耗竭 | OR=0.26, P=2.7×10⁻⁵ | OR=0.259, Fisher 双侧 P=2.70×10⁻⁵ | **PASS（精确）** |
| HAR∩hCONDEL 重叠 | φ=0.052, Yates P=4.6×10⁻⁴ | φ=0.0520, Yates P=4.60×10⁻⁴ | **PASS（精确）** |
| 期望重叠 10.7 | 476×112/4,974 = 10.72 | **PASS** |
| srv 外推 | 49.4→40.7 (−17.5%) → 28.0% ≈1,394 | 34/194=17.53%；0.33977×0.8247=28.02%；×4,974≈1,394 | **PASS（算术；推断见 M1）** |
| BH 压力测试 | 1,690→1,639 (−3.0%) | 51/1,690 = 3.02% | **PASS** |
| Selectome 零结果 | OR=1.21, P=0.41 | OR=1.207, Fisher 单侧 P=0.413 | **PASS** |
| Shao 条件分析 | OR=0.55, P=0.21；期望 6.1 | OR=0.554, P=0.214；18×0.3398=6.12 | **PASS** |
| 池化 OR 家族对比 | z=3.49, P=4.8×10⁻⁴ | z=3.47, P=5.2×10⁻⁴（CI 舍入差） | **PASS** |
| 性状族 Fisher | 8/8 vs 2/4, P=0.091 | P=0.0909 | **PASS（精确）** |
| LOO 权重重整化 | 四组归一化方案 | 逐组验证和为 1 且比例正确 | **PASS** |
| Holm（caMPRA LOO） | 0.032 → 0.16 | ×5 | **PASS** |
| 分类计数 | 1,214+52+293+11+3,404 | =4,974 | **PASS** |
| RELAX 链 | 728=604+124；124−7=117；604/1,660=36.4% | 复算一致 | **PASS** |
| Table 2 迁移 | −406 (−25.1%)；+168 (+134.4%) | 复算一致 | **PASS** |
| 摘要词数 | 236–240（三口径 <250） | 我的机械分词 251，扣除十进制/千分位数字的分词膨胀后 ≈239 | **PASS** |
| 覆盖率 | 97 (1.95%) / 476 (9.6%) / 112 (2.3%) / RELAX 33.4% | 97/4,974=1.950%；476→9.57%；112→2.25%；1,660→33.37% | **PASS** |
| Spearman ρ=−0.06, P=7.5×10⁻⁶ | — | 无结近似 P≈2.3×10⁻⁵；考虑 p=0/0.5 两处原子造成的大量结，方向与量级合理，无法在无原始数据下精确复核 | 基本可信 |

**未能复算（需作者澄清，非指控）：** author-reported 映射的 EA OR=2.99（m2c）；"59"分母（m2f）。

---

## 6. 优点记录（供编辑参考）

1. 有界支撑 + chi-bar-square 的 BH 论证正确且完整：0.5 原子恰是 χ̄² null 的理论签名，BH 临界值 ≤0.05 远离边界故精确有效；×2 压力测试（−3.0%）与双 null 校准（evolver 0/400；参数化 bootstrap 拒绝率 0.26%）构成三角验证。这是本文方法上最出彩的一笔。
2. 单侧 greater 约定在全文（含 Table S5 双侧参照列、GD 耗竭走双侧）执行一致，且在 Methods 显式声明——罕见地干净。
3. 对自己 GWAS 锚定的长度混杂拆解（logistic + CMH + 双映射约定）是负结果报告的范本；synaptic 基因集侧补做了长度匹配 null（1,000 次十分位抽样）与 logistic 校正的双重检验，堵住了 GWAS 侧暴露的同一机制。
4. LOO 设计把循环性从"讨论中的承认"变成"矩阵中的数字"，hCONDEL 作为主独立验证的遴选逻辑（发现口径独立、未被人类细胞系测定、与 HAR φ=0.052）论证严密。
5. GDS 名义权重 0.65 实为同一 LRT 的两个单调变换——作者自己点破（L36），这种自我披露应被鼓励而非惩罚。

---

## 7. 结论

**7.5/10，Minor Revision。** 核心统计架构（48 族单侧 BH、LOO-Holm、Tier 1 宇宙背景、双 null 校准）经逐行复算无一计算错误；三条 major 意见全部属于推断呈现层（外推无 CI、功效论证替代方案、衰减链混杂声明），加上一簇可逐条勘正的小数字错误。作者团队展现出罕见的自我审计文化，预计一轮小修即可达标。
