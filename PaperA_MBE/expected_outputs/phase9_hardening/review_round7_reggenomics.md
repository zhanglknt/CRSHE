# Round 7 盲审报告 — 调控基因组学视角（MBE）

**Reviewer:** r7-reggenomics（人类调控基因组学：HARs/hCONDELs、增强子进化、MPRA/caMPRA、脑发育调控）
**稿件:** manuscript_english_v9.md + Supplementary_Text_v9.md + Main_Tables_v9.md + Figures 4/7（figures_v2）

---

## 1. 总评与决定

这是一篇在循环性控制和证据透明度上远超领域常规水准的工作：hCONDEL 的 hg18→hg38 liftOver 修复链完整披露、caMPRA 嵌套问题被作者自己拆穿并降级处理、全部关键阴性结果（hCONDEL-only 子集不显著、GWAS 长度校正后全灭）均如实上报——但"以 hCONDEL 为主要独立验证"的论证仍有一个未被正面处理的共享保守性 ascertainment 轴，使该验证的"独立性"比我期望的弱一档。

**推荐决定：Minor Revision**（若下文 Major 1 的交集敏感性检验推翻了 hCONDEL 富集，则升级为 Major Revision）

## 2. 评分表（1–10）

| 维度 | 分数 | 说明 |
|---|---|---|
| 新颖性 | 7 | 问题经典（King–Wilson），新颖性主要在 comparability-corrected 框架与 LOO 验证矩阵，而非数据或概念突破 |
| 方法严谨性 | 8 | liftOver 双路由核对、LOO 全矩阵、窗口敏感性、协变量扩展均属上乘；扣分项见 Major 1–3 |
| 统计分析 | 8 | 单侧约定、Holm/BH 家族划分、Yates、Newcombe CI、功效披露意识都好；hCONDEL-only 缺功效声明 |
| 结果可信度 | 7 | 富集信号本身可信；"独立验证"的措辞强度略超证据实际强度 |
| 写作清晰度 | 7 | 极度诚实但句子过长、限定语嵌套过深，关键 caveat 埋在 16 条 Limitations 里 |
| **总体** | **7** | 达到 MBE 水准，修后可达 8 |

## 3. Major comments

**M1. hCONDEL 验证的"独立性"回避了与 nc-divergence 组件共享的保守性 ascertainment 轴。** 稿件把 hCONDEL 的独立性论证建立在两点上：与 HAR 低重叠（22/112，φ=0.052；第 76 行）和"未在人类细胞系中 assay"（第 76、150 行）。但 hCONDEL 的定义本身是"在黑猩猩及其他哺乳动物中深度保守、在人类中缺失的元件"（McLean et al. 2011），而 RDS 的最大权重组件之一 nc conservation divergence（phyloP gene − phyloP CDS，权重 0.25，且是 52% RD 基因的 dominant 组件，第 168 行）同样构建于 100-vertebrate phyloP 深度保守性轨道上。这意味着：非编码保守性高的基因座天然拥有更多可缺失的保守元件，hCONDEL-overlap 概率随之机械性升高——RDS 的 nc 组件与 hCONDEL 验证共享同一个 ascertainment 轴。稿件在 Limitations (7) 只说 "hCONDELs partially mitigate" 细胞系问题，从未讨论这一耦合。LOO-nc 变体（hCONDEL OR = 4.67，95% CI 2.37–9.21；Table 3）在统计上确实缓解了这一担忧——nc 组件移除后富集反而更强——但这是我认为必须在正文第 76 段的"独立性"论证中明确写出的一句话，否则审稿人很容易指控作者选择性定义了"独立"的两个轴而回避第三个。**要求：** 在第 76 行段落和 Discussion 第 150 行处加入对 conservation-ascertainment 耦合的讨论，并明确引用 LOO-nc 结果作为回应。

**M2. 112 基因集是 name 路由（80）与 coord 路由（102）的并集，缺少交集（n=70）敏感性检验。** 第 28 行披露两路由仅在 70/112 基因上一致（Jaccard 0.625），且各自单独均显著（name OR = 2.90, P = 0.0023；coord OR = 3.58, P = 2.1 × 10⁻⁵）。但主分析用的是并集（80 + 102 − 70 = 112）。并集是三种选择中最宽松的一种——它纳入了 42 个仅由单一路由支持的基因。Jaccard 0.625 并不高，考虑到 hCONDEL 的靶基因归属本来就比 HAR 更不确定（McLean 自己的经典例子是远距离增强子缺失），单路由基因的错误归属风险不可忽略。**要求：** 报告交集（n = 70）下的 RD 富集 OR 与 P 值；若交集仍显著（预期会，因为并集 OR = 3.19 落在两单路由之间），一句话即可，但这能把 M1 的论证闭环。同时请明确说明主分析采用并集的理由。

**M3. hCONDEL 目录用的是全部 583 个位点，而非 McLean 原文 sequence-validated 的 510 个。** 第 28 行写明"583 hCONDELs from McLean et al. (2011; their Supplementary Table 2; 510 of the 583 were sequence-validated in that study)"，Fig. 4f 标注也是 "112 genes overlap 583 hCONDELs"。也就是说 73 个未经验证的低置信缺失进入了验证集。对一篇以 hCONDEL 为"主要独立验证"的稿件，验证集纯度应当从紧。**要求：** 要么改用 510 个验证位点重跑富集（预期 112 基因集变化很小），要么明确论证纳入 73 个未验证位点不影响结论。这是低成本高收益的修复。

**M4. "主要独立验证"的强度表述超出了子集分解所支持的水平。** Fig. 4c 的三子集分解（双支持 9/22 OR = 11.38；HAR-only 72/454 OR = 3.67；hCONDEL-only 9/90 OR = 1.80, P = 0.081 n.s.）被诚实披露，值得肯定。但请注意这组数字的含义：hCONDEL 全集富集（OR = 3.19）中信号最强的是与 HAR 重叠的那 22 个基因；剔除 HAR 后，纯 hCONDEL 子集在 n = 90 下不显著。作者用"assay 独立性 + 全集富集强度"来支撑主验证地位（第 76 行），这在逻辑上成立，但在经验上，纯 hCONDEL 子集功效极低（基线 RD 率 ~9–10%，n = 90 时即使真实 OR = 3，Fisher 单侧功效也只有 ~0.4–0.5）。**要求：** 在第 76 行或 Fig. 4c 图注中给出 hCONDEL-only 子集检验的功效量化（post-hoc power 或达到 80% 功效所需样本量），把"方向一致但不显著"从一句修辞变成一个可评估的统计陈述。

**M5. caMPRA 证据链的处理总体诚实，但"508 active HARs"的活性定义与 Shin et al. 原文口径的差异需要更显眼的位置。** 第 28 行括号内说明 16.1%（508/3,171，任一 assay day 元件级活性）与 Shin 报告的 14.2%（N2A 聚合口径）不同——这是合理的再分析，但它埋在方法长句的括号里。此外 caMPRA 活性是在一种细胞背景（N2A）中测得的 enhancer 活性，既非体内、也非发育时间序列；Limitations (7) 承认"activity ≠ human-specific regulatory change"，但正文第 104 行"experimentally validated active HARs"的措辞偏强。**要求：** 正文首次出现 508 处即标注口径差异与细胞背景限制；"experimentally validated"建议改为"caMPRA-active (N2A)"。

## 4. Minor comments

1. **第 28 行 liftOver 程序细节不足：** 请说明所用 chain 文件（hg18→hg38 直链还是经 hg19 中转）与 pyliftover 版本；2 个无 chain 位点的处理方式已说明（剔除），很好。
2. **第 28 行 HAR 目录异质性：** Girskis 3,171 合并目录整合了 7 个发现研究，各研究的物种比对深度与显著性阈值不同；建议在 Supplementary 中给出按发现研究分层的 caMPRA 活性率或 HAR-proximate 富集，确认信号不集中于某一旧目录（如 Pollard 2006 的 HAR1–HAR49 系列）。
3. **±50 kb 基因归属窗口（HAR 与 hCONDEL 皆然）：** Limitations (10) 的 25/50/100 kb 敏感性（OR 单调稀释但仍显著）令人安心；但窗口映射本质是 crude proxy，已知相当比例 HAR/hCONDEL 经 3D 染色质接触跳过最近基因。建议在未来工作一句中提及 TAD/Hi-C 或 fetal-brain eQTL 归属作为改进方向（不必本次补做）。
4. **第 76 行 φ = 0.052 的解读：** 期望重叠（476 × 112/4,974 ≈ 10.7）下观察到 22，实为 ~2 倍富集（Yates P = 4.6 × 10⁻⁴ 也确证非随机）。"largely non-overlapping"的表述准确，但建议补一句"重叠本身显著高于随机期望"，避免读者把 φ 小误读为"两证据类型统计独立"。
5. **Fig. 4a/4b 与 Table 3 的 RD 类别规模问题：** LOO 变体间 RD 类别规模从 110（LOO-nc）到 800（LOO-caMPRA）相差 7 倍，OR 在类别定义剧烈变化下仍 >2 是稳健性证据，但 Fig. 4d 的 "purity–signal trade-off" 框架暗示作者已意识到 OR 与类别纯度不可兼得；建议在正文用一句话点明 LOO 矩阵解读的是"信号存活"而非"同一类别"。
6. **Fig. 7a 红条标签：** 图例 "Regulation-driven + dual" 与子集内比率（17.0–19.1%）的检验在正文中（第 124 行）只给了编码选择一侧的 CI/P；红条侧的富集检验 P 值（3.8 × 10⁻⁵ 至 1.7 × 10⁻²⁵，Fig. 7 图注）建议同时在正文给出。
7. **BrainSpan 组件（RDS 权重 0.20）：** tau 跨 26 个脑区计算，但脑区间表达高度相关，tau 对相关性结构敏感；且 BrainSpan 混合了产前/产后样本。建议在 Supplementary 一句话说明脑区清单与发育阶段处理。
8. **第 92 行 "100% nesting" 的表述：** "all 97 caMPRA-positive genes are also HAR-proximate ... because caMPRA-active elements are themselves HARs"——逻辑上是定义性嵌套而非经验发现，建议措辞改为 "by design"，以免被误读为独立观察。
9. **Fig. 4f 与第 54 行统计约定的一致性：** GD 侧耗竭检验用了双侧 P = 0.037（恰当，因为耗竭方向非预设），而第 54 行声明"all enrichment tests one-sided in the enrichment direction"——建议在图注中明确该检验为双侧（目前只标注了 "two-sided P"，正文第 293 行图注里已有，很好，但请在 Statistical conventions 处加一句例外说明）。
10. **Supplementary Text 中缺少 hCONDEL 修复链的专节：** liftOver 修复的全部细节（双路由、Jaccard、被取代的错误版本）只存在于主文第 28 行一句话与"revision history"的引用中；建议在 Supplementary Text 增加 S6 节给出 80/102/70 三集合的完整基因清单与逐位点 liftOver 状态表。这是可复现性的关键一环，目前读者无法独立核对 112 基因集。
11. **Fig. 7d 线性外推：** 两端点验证（476 → 预测 417/实测 410；4,974 → 预测 1,888/实测 1,906）做得很好；但 Limitations (12) 已承认"ignores class migration"——建议在图注中重复这一 caveat，因为该图最容易被单独引用。
12. **第 168 行 caMPRA 另两类元件（CNE 1,188/5,155、variable 405/1,702）的讨论有价值：** 建议给出一句定量的"若纳入 CNE 活性层，assay-positive 覆盖将从 2% 升至约 X%"，让读者直观感受 HAR-restricted 选择的保守程度（可放 Supplementary）。

## 5. 一句话：什么改动能让我提高分数

在正文正面承认 hCONDEL 与 nc-divergence 组件共享保守性 ascertainment 轴（M1），并补做 hCONDEL 交集敏感性（M2）与 510 验证位点重跑（M3）——这三项完成后，方法严谨性与结果可信度均可升至 9，总体 8。
