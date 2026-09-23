# Round 7 评审：MBE 副编辑（Associate Editor）编辑评估

**稿件**：Coding Selection Dominates by Prevalence, Regulatory Selection by Specificity: A Comparability-Corrected 10-Primate Analysis（manuscript_english_v9.md + Supplementary Text v9 + Main Tables v9 + cover letter v1 + Figures 1–7）
**评审人角色**：MBE Associate Editor（第 7 轮，此前未见过此稿）
**评估范围**：只审不改。期刊契合度 / 标题-摘要-结论一致性 / 诚实化后叙事力度 / 图表质量 / 投稿成熟度 / 是否送审

---

## 1. 编辑决定建议：**Major Revision（建议送审，预期一审结论为大修）**

这是一篇期刊契合度很高、方法学自觉程度远超同类稿件的比较基因组学研究，值得进入正式外审。King–Wilson 问题是 MBE 读者的核心议题；稿件对 BUSTED/RELAX 的联用、循环性消除、证据覆盖不对称的定量化（59% 交叉点、HAR 受限上限 410）构成真实的方法学贡献，且 6 轮内部修复（hCONDEL 坐标 liftOver、GWAS 长度校准、BUSTED srv/中性模拟校准链）后，文本的诚实度罕见地高——几乎每个已知弱点都被作者自己量化并写明。但正因为它现在是结论式标题，编辑必须按"标题主张是否被最终版本的证据支撑"来把关：**长度校正后，"Regulatory Selection by Specificity" 这一半标题的实证基础已从"8/8 神经精神 GWAS 性状富集"收缩为"突触基因集（SynGO/GO-BP）富集"单一支柱，而这一支柱尚未做与 GWAS 锚定同等级别的基因长度敏感性检验**。RD 基因长度是基因组中位数的 4.7 倍，而突触/神经元基因恰是全基因组最长的基因类群之一；既然长度校正能把 OR 1.9–4.0 的性状富集全部打到 q > 0.44，就不能排除同样的机制部分或全部解释 SynGO 富集。LOO-brain/LOO-tau 解决的是表达组分循环性，不是长度混杂。在这一检验补上（或标题相应降级）之前，我不能以当前标题放行。加之摘要超 MBE 上限（304 vs 250 词）、Funding/AI 披露/Zenodo DOI 三处占位符未填，属大修而非小修。预计一轮聚焦修订即可，无需多轮。

## 2. 评分表（1–10 整数）

| 维度 | 分数 | 简评 |
|---|---|---|
| 期刊契合度 | 8 | 灵长类正选择扫描 + BUSTED/RELAX 方法学 + 调控演化，正中 MBE 读者群；方法审计框架对该领域有直接可迁移性 |
| 新颖性 | 7 | 科学问题经典（King–Wilson），新意在于"可比性校正 + 覆盖度定量化"框架而非概念本身；59% 交叉点与 disjointness 检验是可被引用的增量 |
| 方法严谨性 | 7 | 敏感性/校准链条异常完整（LOO 全矩阵、1000 次权重扰动、协变量、中性模拟、srv、拓扑、映射惯例）；扣分点：无比对过滤、GARD 全线不收敛、复合评分权重本质上是启发式的、SynGO 富集未做长度校正 |
| 叙事与写作 | 7 | 逻辑骨架清晰（prevalence → disjointness → anchoring → coverage），"两条战线"模型立得住；但 caveat 密度过高，长句+破折号堆叠，Discussion 近半篇幅是限制陈述，读者负担重 |
| 图表质量 | 7 | 7 张多面板图风格统一、信息量大、标注专业（Fig 1/4/7 尤其好）；扣分点：Fig 6 富集 term 标签被截断（"Regulation Of Lymphocyte Prolifera..."）、Fig 2d 拟合线渲染越界怪异、Fig 3d 与正文口径需一句话澄清（见 §3） |
| 投稿成熟度 | 6 | 距可投还差一轮聚焦修订：一个关键补分析 + 三处占位符 + 摘要压缩；不需要结构性重做 |
| **总体** | **7** | 有真实贡献、诚实、成熟度高但未完工；标题成色问题是唯一硬约束 |

## 3. 给作者的三条最重要建议（按优先级）

**（1）为"specificity"主张补齐长度敏感性检验——这是结论式标题当前唯一的实证支柱，必须先加固。**
对 RD 类的 SynGO（11 terms）与 GO-BP（18 terms）宇宙背景超几何富集，做与 GWAS 锚定同等级的基因长度校正：长度匹配的对照基因集置换（按 log₁₀ gene span 分层抽样背景）、长度协变量逻辑回归、或长度三分位分层超几何，三取其二即可。若 SynGO 富集在长度控制后仍然 FDR 显著，标题原样成立，且在正文中明示"SynGO 是 length-robust"这一claim从此有据；若显著性消失或大幅衰减，则必须降级标题与摘要的 specificity 措辞（如改为 "carries the detectable functional signal" 之类非对称性更弱的表述），并把"突触信号是否独立于长基因效应"列为开放问题。无论结果如何，这一节的加入都会让稿件更强。

**（2）补齐投稿就绪项并把摘要压到 MBE 上限内。**
三处占位符必须在投稿前填实：Funding（"[Funding sources and grant numbers to be added.]"）、作者贡献段的 AI 使用披露括号、Data availability 的 Zenodo DOI（"to be assigned upon release"——MBE 要求录用时有可解析的存档 DOI，投稿时至少给出预留 DOI）。摘要当前 304 词，超 MBE 约 250 词上限；建议删掉八个 GWAS 性状的逐一 OR 罗列（反正下一句就被长度校正撤回），用一句"trait-level anchoring was directionally consistent but abolished by gene-length correction"概括，把省下的篇幅给 59% 交叉点——这是全文最有引用价值的数字。

**（3）收敛叙事密度，统一图-文口径。**
Discussion 目前几乎每段以限制陈述收尾（coverage transparency / circularity / gBGC / ILS / cross-study），科学上可敬、阅读上疲惫。建议：将 Limitations 已有的 16 条保留，但把正文中与其重复的次级 caveat（如 Shao 交叉研究细节、caMPRA 三元素类讨论）迁入补充材料；把"two largely disjoint fronts"模型提前为 Discussion 第一段的主题句，让 caveat 服务于模型而非淹没模型。图文口径：正文阈值敏感性报的是严格 GD（1,382→918），Fig 3d 画的是 GD incl. relaxed（1,556→942），二者可推出一致但读者会对不上号——在 Fig 3d 图注或正文加一句口径说明即可。Fig 6 的 term 标签截断请改用换行或缩写字典重绘；Fig 2d 残差化前后面板的拟合线越界渲染需修正。

## 4. 标题与摘要的具体改进意见

**标题。** "Coding Selection Dominates by Prevalence" 一半是可辩护的：prevalence 比较已显式条件化于覆盖度（59% 交叉点、HAR 受限上限 410 均在正文量化），且 RELAX 确认集（12.2%）与 GD 复合分类（24.4%）两种口径下方向一致。"Regulatory Selection by Specificity" 一半当前成色不足：原文所依据的 GWAS 特异性已被作者自己的长度校正撤销（12 个性状全部 q > 0.44），且阴性对照（身高、T2D）同幅度富集说明该锚定本就不具神经特异性；现存特异性证据 = 突触基因集富集 + hCONDEL 独立验证，前者未做长度控制、后者是"类存在性"验证而非"功能特异性"验证。两个选项：（a）完成建议（1）的长度校正后保留现标题；（b）立即 softened，例如 "Coding Selection Dominates by Prevalence; Regulatory Selection Carries the Synaptic Signal: A Comparability-Corrected 10-Primate Analysis"——把 specificity 替换为已被数据直接支撑的具体内容（synaptic signal），即使长度校正结果不利也无需再改题。

**摘要。** 除压缩到 ≤250 词外，三处具体修改：
1. "enriched for all eight neuropsychiatric and cognitive GWAS traits (odds ratios 1.9–4.0) — and, to a similar degree, for height and type 2 diabetes" 这句先扬后抑的结构占用约 60 词，建议整句替换为对长度敏感性的直接陈述，并点明"the length-robust core is the synaptic gene-set enrichment"（该句已存在，保留即可）。
2. 末句 "trait anchoring directionally robust but gene-length sensitive" 与标题的 "Specificity" 之间存在张力，编辑上建议摘要结论句与标题措辞严格对齐——标题改了，这句就顺了。
3. "full-coverage counterfactual: 293–1,906" 建议改写为更具信息量的形式（如 "regulatory class size scales linearly with assay coverage and would match the coding class only at ~59% genome coverage"），把 59% 这个可被媒体和综述引用的数字直接放进摘要。

**Cover letter（顺带）。** 写得专业且诚实（主动前置 34% 上界问题），建议审稿人名单合适；唯一提醒：Pollard 与 Noonan 是被大量复用数据（HAR/hCONDEL 目录）的原作者，部分期刊会视为潜在利益冲突，建议各备一名替代人选。

---

## 附：编辑核查中发现的具体事实点（供作者与审稿人参考，非独立审稿意见）

- 摘要词数 304（MBE 上限 ~250）；全文约 11,400 词（含参考文献与图注），正文主体约 7,500–8,000 词，对 MBE Discoveries 可接受。
- Fig 3d 数值（GD incl. relaxed: 1,556→942）与正文（strict GD: 1,382→918）可相互推出（差值为 relaxed 类在各 margin 下的大小：174→52→24），非错误，但需口径说明。
- 摘要 HAR LOO OR 范围 1.6–9.4、hCONDEL 2.1–4.7 与 Table 3（1.61–9.39 / 2.13–4.67）一致；分类计数（1,214/52/293/11/3,404）在摘要、Table 1、Fig 1c、Fig 3a 间一致。
- Limitations 第 (15) 条与摘要、正文 §"A conditional answer" 关于长度校正的表述三处一致（all adjusted q > 0.44；EA 与 SCZ 仅在粗分层/替代映射下保留）。
- 补充文本 S4 的 LOO 重算（RD = 293/800/148/110/613）与 Table 3 完全一致，复现性声明可信。
