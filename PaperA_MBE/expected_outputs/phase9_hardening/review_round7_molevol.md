# Review — Round 7 (Molecular Evolution Methods, MBE)

**Manuscript:** Coding Selection Dominates by Prevalence, Regulatory Selection by Specificity: A Comparability-Corrected 10-Primate Analysis (v9)

## 1. 总评与推荐决定

This is an unusually self-critical comparative-genomics study whose BUSTED+RELAX design and circularity/coverage audits are methodologically valuable, but the headline 34% detection rate rests on an alignment and model configuration (non-codon-aware alignment, no filtering, --srv No) that the authors' own sensitivity analyses identify as the least defensible choice, and the calibration experiments cannot detect the dominant false-positive mechanisms. **Major Revision.**

## 2. 评分表（1–10）

| 维度 | 分数 |
|---|---|
| 新颖性 | 7 |
| 方法严谨性 | 5 |
| 统计分析 | 7 |
| 结果可信度 | 6 |
| 写作清晰度 | 8 |
| 总体 | 6 |

## 3. Major comments

**M1. 比对流程不足以支撑 dN/dS 全基因组扫描的 headline 数字。** Methods（manuscript line 28）只说 "MAFFT alignment and PAML-format QC"，正文 Limitations（line 176）承认 "no alignment filtering (GUIDANCE/trimAl) was applied"。对密码子 dN/dS 分析而言，更根本的问题是：MAFFT 默认是核苷酸比对，不保证密码子框完整；基于核苷酸的 CDS 比对会在近缘灵长类中引入少量 but 系统性的移码 gap 与错配密码子列，而这些恰恰是基因级 LRT 的主要假阳性来源之一（这是 branch-site / gene-wide 检验文献中的标准结论）。领域惯例是 PRANK-codon、MACSE 或 translatorX 之类的 codon-aware 比对加列质量过滤。作者在 Limitations 中把 "no alignment filtering" 列为使 34% 成为上界的理由之一，但从未实测其影响（CpG/拓扑/srv 三个敏感性都测了，唯独比对质量这个文献中最主要的假阳性来源没有测）。**必须**：至少在与现有 400 基因分层重跑相同的样本上，用 codon-aware 重比对 + GUIDANCE/trimAl 过滤重跑 BUSTED，给出百分点变化，与 CpG/srv 并列报告。

**M2. 中性模拟设计无法校准最主要的假阳性机制，0/400 的说服力被高估。** evolver 模拟（line 32、line 176：ω = 1、F3x4、κ = 2、真实缺失模式）有两个结构性盲区：(i) 若模拟未引入 indel（正文未说明引入），模拟序列经 MAFFT 比对几乎无错误，因此该实验**在设计上不可能**检测 M1 中的比对错误假阳性——"empirical FDR 0/400" 只证明流程对理想比对是保守的，而这几乎没有人怀疑；(ii) 零模型取全基因 ω = 1 与现实不符：真实零假设应含强烈的纯化选择成分（BUSTED 的 null 正是带 ω < 1 组分、ω ≥ 1 组分权重为 0 的混合分布）。正确的校准是对每个基因在其**拟合的 BUSTED null 参数**下参数化模拟（parametric bootstrap），这同时能检验 srv 未建模（M3）造成的 I 类错误膨胀。另请明确说明模拟是否含同义替换率变异（srv）；若不含，则该模拟对 srv-No 生产运行最关键的失效模式是沉默的。此外，400 个模拟中 0 个 BH 发现给出的 FDR 上界相当宽松，建议扩大规模或明确报告置信上界。

**M3. 生产运行采用 --srv No，而作者自己的数据表明 srv 是唯一实质性敏感轴——headline 应建立在更可信的设置上。** HyPhy 团队对大规模扫描的建议是建模 srv；作者的分层重跑显示 --srv Yes 使检出率从 42.6% 降到 34.6%（−8.0pp，line 142），且 110/400 超时基因 69% 为生产显著（偏长基因），即估计来自显著性被稀释的子集。"绝对/相对外推 ~26–28%" 依赖对非随机缺失的 289 个基因的外推，脆弱性高。我理解全宇宙 srv Yes 重跑成本高，但当前论文的叙事是 "34% 是主结果、27% 是敏感性的下限"——更诚实的框架应反过来：srv Yes 估计是主估计，srv No 是上界。至少在摘要与 Discussion 的 headline 表述（"27–34% 取决于 srv 建模"）中明确哪一端对应哪个模型假设，并对外推的置信度给出区间而非点值。

**M4. GARD 50/50 全部不收敛更像实现问题而非数据性质，重组混杂目前完全无界。** Supplementary Text S2 报告 GARD 在 50 个最显著基因上 "all 50 alignments" 数值收敛失败。10 物种、数百至数千密码子的灵长类比对上 GARD 全军覆没并不典型，强烈提示运行配置问题（rate classes、初始树、序列数 vs. 核苷酸模式组合）。作者已诚实地将其报告为限制，但后果是：整篇论文（BUSTED、MEME）都假设无重组，而没有任何一个基因得到过成功的重组筛查——S2 中 "simulation studies suggest partial robustness" 的辩护没有引文支撑且作者自己也承认 "cannot verify this property in our own data"。**必须**：排查失败原因（尝试 codon 模式、减少 rate classes、分批或更新 HyPhy 版本），至少让一部分基因通过 GARD；若确实全部失败，应如实说明所用 HyPhy/GARD 版本与完整错误信息，并考虑用其他断点检测方法交叉验证。

**M5. RELAX 松弛约束基因计数内部矛盾：117 vs 52。** Results（line 64）先写 "117 (7.0%) FDR-significant relaxed constraint (0 < K < 1; seven FDR-significant genes with degenerate K = 0 boundary estimates … are excluded)"，紧接着写 "The 52 genes that are BUSTED-positive with FDR-significant 0 < K < 1 are reclassified"。RELAX 只跑在 BUSTED 显著基因上，因此 "BUSTED-positive 且 FDR-显著 0<K<1" 的基因数按前文应为 117−7=110，而非 52；后文（Table 1、摘要）又统一用 52。若 52 = 110 中再通过 gene-driven margin（GDS > RDS + 0.15）者，必须明确写出这一推导链；否则是数字错误。读者目前无法复算。

**M6. "压缩零分布、双倍 p 值" 检验的统计论证不正确（结论侥幸不受影响）。** Limitations（line 176）称 "under the more conservative interpretation that the null is compressed to U(0, 0.5) and p-values should be doubled, the significant count changes from 1,690 to 1,639"。但在 BUSTED 的 chi-bar-square 零分布（0.5χ²₀ + 0.5χ²₁）下，LRT = 0 处的 0.5 原子恰好补偿了连续部分缺失的质量：对任意 t < 0.5，P(p ≤ t) = t，即 p 值在 (0, 0.5) 上**精确随机均匀**，不存在需要 "双倍" 的压缩。BH 在临界值 ≤ 0.05 时的有效性是精确的而非近似的，这反而是比 "1,690→1,639" 更强的论证。请修正这段推理（保留数值检验无妨，但应表述为过度保守的稳健性检查，而非 "the more conservative interpretation"）。这段文字目前会让懂混合卡方的审稿人对全文的统计论证质量产生不成比例的怀疑。

## 4. Minor comments

1. "BUSTED v3"（line 32）引用的是 Murrell et al. 2015（原始 BUSTED）；v3 实现应对应更新的 HyPhy 版本引文（如 Weaver et al. 2018 / HyPhy 2.5），并请报告具体 HyPhy 版本号。
2. 模拟参数 κ = 2（line 32）缺乏依据；灵长类编码区 κ 通常更高（~3–5）。F3x4 频率从所有比对合并估计 vs. 逐基因估计的选择也应说明——逐基因 F3x4 更接近生产运行的拟合条件。
3. srv 重跑 "400 中 110 超时" 后应剩 290，正文写 "289 comparable genes"——差的 1 个基因去向请说明。
4. 592 个 p = 0 精确值（Fig. 2d、S1(ii)）实为数值下溢（< 1e-300 下限）；建议正文统一表述为 p < 1e-300 而非 "p = 0"，后者在统计上不规范。
5. Fig. 2a 的 p 值分布在 (0, 0.5) 中段几乎为空、两端堆积——这种极端双峰更像功率二分（长基因全显著、短基因全不显著）而非信号/噪声混合。正文（line 60）提到双峰违反均匀零假设，但未讨论其与 CDS 长度功率梯度的关系；考虑到 GDS 已做长度残差化，这一讨论是必要的。
6. Storey π₀ = 1.0 分析（line 60）作者自己承认是 "internal consistency check, not independent evidence"——既然不提供信息，建议删除或移入补充，避免给读者留下 "做了独立 FDR 校准" 的印象。
7. MEME 部分（line 32、S2）请报告 12 个超时基因与 8 个完成基因的比对长度分布对比，量化 "coverage biased toward strong gene-wide signals" 的程度；目前只有定性表述。
8. 前景枝仅为人类末端枝（newick 中 human#1，S1）意味着 "human-specific" 实际涵盖人–黑猩猩/倭黑猩猩分歧后约 6 My 的全部事件；建议在 Results 开头一句点明，避免读者误读为近端（如 <1 My）事件。
9. GC3 协变量仅用人类 CDS 计算（line 42）——对 gBGC 这恰是正确做法（gBGC 是谱系特异的），但请明说这一理由，否则读者会质疑为何不用祖先重建或跨物种均值。
10. 摘要中 "empirical FDR 0/400" 与正文一致，但摘要未提模拟为 ω = 1 理想比对零假设（见 M2）——建议在摘要或 headline 数字旁加半句限定。
11. Table 4 基线行 "99.96% 一致" 的脚注（percentile-convention edge cases）值得赞赏，但 0.04% 对应约 2 个基因，建议直接说明是哪些基因及方向。
12. Fig. 7d 中外推交叉点给了两个版本（单一线性外推 ~2,914 / 双模型交叉 2,158），正文（line 128）均有交代，但图注（line 305）只提 ~2,900/59%；请使图注与正文一致或在图注中说明两版定义。

## 5. 什么改动能让我提高分数

在全宇宙（或至少与现有分层样本相同的更大样本）上以 codon-aware 重比对 + 列过滤 + --srv Yes 重跑 BUSTED，并用含 indel、按逐基因拟合 null 参数的参数化模拟重新校准 FDR——把 "27–34%" 的区间收敛为一个建立在最可信模型假设上的单一数字。
