# R6 审稿回复 — 决策与修改日志（滚动更新）

来源：`review_round6_merged.md`（4位新专家：primate-genomics 7.5 / popgen 7.0 / molevol-methods 7.5 / editor 8.0，均分 7.5）

## 盲验证 + 最终关闭（2026-09-15 深夜）

- **盲验证员（fresh agent，独立重算）**：62 项检查 **56 PASS / 6 FAIL**；三链 headline 全部精确复现（含 115 万行 GWAS TSV + GENCODE GTF 全量重扫、种子 20260917 置换逐位置复现）。6 项 FAIL 全部为表述级（无一影响结论）：
  - #1 HAR p 范围 3.3e-20→**3.6e-25**（下端点漏了 LOO-nc，实际更强）✅ batch 7
  - #2 "89 of the 97"→**all 97** 在 universe 内（89 溯自被取代的 phase8 口径）✅ batch 7
  - #3 Height/T2D OR 范围 2.2–2.3/2.1–2.3→**1.9–2.3/1.8–2.3**（LOO-brain 值在范围外）✅ batch 7
  - #4 长度比 1.6–3.5→1.3–4.1 ✅（batch 6 已修，验证员读的是旧版）
  - #5 srv 34.5%→**34.6%**（289 可比基因口径；34.5 为 290 全组率）✅ batch 7
  - #6 caMPRA LOO 检验 1.84/0.012/0.059 为旧 89 基因集口径→生产 97 基因集 **1.64 (CI 1.02–2.63) / 0.032 / Holm 0.16**（定性结论两口径一致）✅ batch 7（自算复现 k=23/97, OR=1.640, p=0.0315）
- **交付链最终重建**：docx×4 → 包同步 → 单 PDF 54 页（11/11 探针过）→ zip 6.15MB/180 文件/CRC OK → CRSHE 提交 **c9fab1a** 推送
- 附带修复：loo_full_matrix.csv（Table 3 源文件）hCONDEL 列三处刷新（results/paper + 包内 + CRSHE），消除与 Table 3 xlsx 的口径不一致
- 验证报告入 CRSHE：blind_verification_r6.md + 两份 recheck

## 审稿人复核（2026-09-15）

- **r6-molevol-methods**：P0-1/P0-2 均 **CLOSED**，"稿件在 dN/dS 检测统计基础维度已达 MBE 可发表标准，修正 2 处文本残留后无需再审"。两处残留已修（batch 6）：①srv 超时偏倚方向写反（非"偏小基因"——复算超时基因 69% 为生产显著、即偏**长**基因；−8.0pp 来自显著稀释子集已披露，绝对/相对外推一致 ~26–28%）②Methods L52 "genome-wide rate"→"respective complement sets"。小建议已采纳：模拟句补 "all 400 completed"。
- **r6-popgen**：P0-3/P0-4 + P1 全部 **CLOSED**，复核意见 "Minor revision → accept"。新增 2 个改进项已修（batch 6）：③CMH 与映射敏感性完整披露（CMH 下 height 1.33/T2D 1.42 亦名义显著、intelligence 不显著；locus-dedup 下 intelligence OR=4.09 与 height OR=3.06 显著、author-reported 仅 EA/SCZ）④长度比范围 1.6–3.5→1.3–4.1 倍。
- 复核文件：review_round6_molevol_recheck.md / review_round6_popgen_recheck.md
- 盲验证员（独立重算三链）：在跑
- 遗留：Zenodo DOI 占位符（P2-7，交付时替换）

## P0 项处置

| # | 问题 | 处置 | 状态 |
|---|------|------|------|
| P0-1 | 34% 上界需有界双峰 p→BH 模拟校准 | fix-busted-wsl **全部完成**：任务1 解析（截断 BH 不变 / 压缩 −3.0%）+ 任务2 中性模拟（400 个 ω=1 evolver + 真实缺失模式 + 同一管线）：机制判定=截断（45.8% 点质量，无 p>0.5）；p<0.05 仅 1.75%（<5% 名义，偏保守）；经验 BH FDR=0/400（α=0.05/0.10/0.25 全 0）；真实数据 35.8% p<0.05 vs 零假设 1.75% —— 小 p 富集非伪影。已写入 Limitations (2) + cover letter。环境基线复跑一致性：检出率 49.4% vs 49.9%，p 相关 0.97，逐基因一致 98.5% | ✅ 已修复 |
| P0-2 | 34% 需分层重跑（--srv Yes / CpG mask） | fix-busted-wsl **完成**（948 分钟）：baseline 复跑一致性极佳（49.4% vs 49.9%，r=0.97，98.5%）；CpG 屏蔽 −1.3pp（47.8%，保留 93.8%）；替代拓扑 −1.3pp（46.6%，保留 92.1%）；**--srv Yes 唯一实质敏感**：−8.0pp（42.6%→34.5%，289 共同基因，110 超时偏小基因），相对降 19% → 宇宙外推 ~27%。表述改为 "27–34% 取决于 srv 建模"（batch 5 六处：Methods 敏感性设计 / Interpreting 34% 测量化 / L60 / Limitations(2) / ILS 拓扑证据句 / cover letter） | ✅ 已修复 |
| P0-3 | GWAS 基因长度/映射混杂 | fix-gwas 完成：**审稿人质疑成立** — RD 基因中位 133.6kb vs 28.6kb（4.7×，p=1.5e-71），nc 分量驱动（ρ=0.47）；logistic 校正后 12 性状全部失去显著（OR 1.02–1.64，q>0.44）；CMH 粗分层 EA（OR=1.87, p=3.9e-6）与 SCZ（1.95, p=8.2e-4）存活；作者报告基因映射 EA p=0.015、SCZ p=0.007 存活。文本已集成（batch 3，15 处） | ✅ 已修复 |
| P0-4 | "strongest and most consistent" 缺组间 log-OR 检验 | fix-gwas 完成：未校正神经族 3.00 (2.60–3.47) vs 对照族 2.02 (1.70–2.39)，z=3.49 p=4.8e-4，置换 p=0.009；长度校正后 1.23 vs 1.06 (p=0.24) — "stronger" 亦被长度解释。一致性 8/8 vs 2/4（Fisher p=0.091，功效有限）。文本已集成 | ✅ 已修复 |
| P0-5 | HAR/hCONDEL 来源审计（3,171 vs Pollard 202；583 vs McLean 510） | fix-provenance 完成：3,171=Girskis 2021 七研究合并目录（Shin 2024 复用 mmc2，hg19→hg38 liftOver 3,169；HARsv2 原生 hg38 3,168 驱动分类管线）；583=McLean 完整目录（510=序列验证子集）。Methods provenance 段已替换（E9），Introduction 加 Girskis 目录注（E11） | ✅ 已修复（正文层） |
| P0-6 | ILS 结构性混杂（~25–30%） | fix-provenance（Discussion 段 + 引文，batch 3 E8）✅ + fix-busted-wsl 替代拓扑实测（−1.3pp，92.1% 保留，实际影响小）已并入 ILS 段（batch 5 B5） | ✅ 已修复 |
| P0-新 | hCONDEL hg18/hg38 坐标错配（修复中发现） | ✅ 已闭环（见上方 batch 4 记录）：183→112，全链路数字+图+表三方一致 | ✅ 已修复 |

## P1 项处置（fix-gwas 完成）

| 项 | 结果 |
|----|------|
| P1-1 映射敏感性 | EA/SCZ 两种替代约定下均显著（作者报告 EA OR=2.99 p=0.015、SCZ OR=1.81 p=0.007）；intelligence（n=70）与 height 失去显著 — 部分支持特异性 |
| P1-2 Jaccard/Meff | 均值 0.148、最大 ASD-intelligence 0.31；Li&Ji 有效检验数=3（Kaiser=2），8 性状 BH 偏保守 |
| P1-3 功效表 | 含 CI 宽度（Crohn log 宽 2.31、LDL 0.91，null 可信） |
| P1-4 anchoring 主表 | tier2_anchoring_master_table.csv（108 行：12 性状×3 分类×RD/GD_strict/GD_all，n/OR/CI/p/q）→ 待打包为 Supplementary Table S2 xlsx |

## P1/P2 项处置（文本层已完成）

| 项 | 内容 | 处置 |
|----|------|------|
| 段落压缩（editor） | Results "Evidence-availability asymmetry" 与 Fig 7d coverage 段重叠 | L86 压缩：删除与 L128 重复的反事实推演句，改为前向引用 "The coverage-sensitivity analysis below (Fig. 7d) quantifies this dependence" ✅ |
| GTEx 措辞（molevol P2-6） | 正文 L28 与 S1 L21 日期措辞不一致 | 统一为 "open-access release 2026-01-16; updated bulk RNA-seq files of 2026-05-19, which separate the laser-capture microdissection pilot samples" ✅ |
| Cover letter 34%（editor） | 预防性说明 | cover_letter_v1.md 加入上界框架句（不预支分层重跑结果）：34% 为上界，类级结论在 strict 子集与更严阈值下不变 ✅ |
| Fig 6 降级（editor，可选） | 降为补充图 | **婉拒**：Fig 6 承载功能内容叙事（5 处正文引用），MBE 允许 7 张主图；压缩诉求已通过段落合并实现；降级将触发 Fig 7→6 级联重编号（10 处引用）+ 补充材料交叉引用重写，风险收益不成比例 |

## 前 batch 已落地（v9 文本 18 处）

- Batch 1（13 处方法学措辞）：树澄清（枝长为起始值，BUSTED 重优化）/ RELAX 运行清单简化 / MEME 选择偏倚句 / 无预注册声明 / PRDS 句 + Benjamini & Yekutieli 2001 / K>1 操作性连用 / "ten strongest terms" / Selectome 94→24 保守低召回 / 单侧 P=0.41（双侧 0.64）/ 取样理由 / Fig 2e omega 注 / 12.2% 操作性限定
- Batch 2（5 处数字整合）：complement-based Tier4 比较（正文+摘要+Fig 7a 图例）/ Limitation 8 strict-only / Limitation 10 窗口敏感性
- Batch 3（15 处，fix-gwas + fix-provenance 集成）：摘要 GWAS 句+结论句诚实化 / Methods 校准方法段（logistic+CMH+双映射+Jaccard+主表 S2）/ Results 组间检验替代 "strongest and most consistent" / Results 第三混杂（长度）整段 / Discussion "in unadjusted tests"+长度段 / Limitations (15) 长度混杂 (16) segdup / ILS Discussion 段 / Methods provenance 段（Girskis 目录+liftOver 口径+583/510 调和+97/89 口径+14.2% vs 16.1%）/ Introduction Girskis 目录注 / 新引文 5 条（Li&Ji 2005, Mendes&Hahn 2016, Scally&Durbin 2012, Scally et al. 2012, Vanderpool 2020）
- 残留清理 2 处：结论段 "unadjusted odds ratios + directionally robust though gene-length-sensitive"、S1 同步长度校准句
- fix-subsets 已完成：Tier4 complement 修复（HAR +2.85pp P=0.222）、strict-only 4,119 全结论稳健、窗口敏感性 25/50/100kb 单调稀释

## 待办（等 worker）

1. ~~fix-provenance hCONDEL liftOver 重跑~~ ✅ 完成：183→112 基因（金标准 within 恢复 7→40/44，旧坐标法仅 7/103 保留证实伪影）；RD OR 2.94→**3.19**（p=7.8e-5）更干净；GD 显著耗竭 OR=0.59 p=0.037；LOO OR 区间 2.3–3.9→**2.1–4.7**；φ 0.0417→0.052；Tier4 全部 ns 结论不变（hCONDEL −3.7pp p=0.48 / 合并 +1.7pp p=0.42）；名称法与坐标法 Jaccard 0.625 双法各自显著。正文 batch 4（12 处）已落地。**诚实化代价**：hCONDEL-only 子集 9/90 OR=1.80 P=0.081 不再显著（旧 17/154 p=0.013），Fig 4c 与 L76 已改写为"方向一致但不单独显著，非冗余论证主要依赖 assay 独立性+全集富集"
2. fix-busted-wsl → BH 校准数字 / 分层重跑 / 替代拓扑 → Limitations (2) 更新
3. ~~anchoring 主表~~ 待打包 TableS2 xlsx（数据已就绪 tier2_anchoring_master_table.csv）
4. ~~fig-fix 已派~~ ✅ 图件全部完成（fig-lead-a + fig-fix 双 worker 汇合互证）：Fig 4（a/b/d 动态读 JSON 重算；c 单侧 greater 口径修正 p=0.081 ✓ doubly p=1.9e-6 ✓；f GD 双侧 0.037 ✓）+ Fig 7a（n=112/566 + assert）→ results/paper/figures_v2/（17:43 最终版）；Table 3 主表 md + xlsx 两处（results/paper/ 与包内 04_tables/）均已核验新值（3.19/2.13/4.60/4.67/2.56 + k/n + 单侧 p）✓。**图-文-表三方一致性已核验**
5. 交付链重建 + 交叉验证
6. **用户决策点：标题** — "Regulatory Selection by Specificity" 的 GWAS 支柱被长度混杂削弱后，specificity 主要落在 SynGO 突触富集（每分类变体下 RD 独有）。选项：A. 维持标题（默认）；B. "...Functional Specificity"；C. "...Targeted Function"。待用户定夺。


## Batch 8 (2026-09-15): 盲验证 6 条附带口径备注修复

对象：blind_verification_r6.md 附带口径备注（不计 FAIL，但按用户指令一并核销）。每条均先从 phase9_hardening JSON/CSV 独立核数再修改，共 7 处替换（φ 双处）：

| # | 备注 | 独立核数 | 修改 |
|---|------|---------|------|
| 1 | Fig 4f GD OR=0.59 实为含 relaxed | strict GD k=17/1,214 OR=0.5479；GD+relaxed k=19/1,266 OR=0.5923（hcondel_gene_mapping_fixed.csv × v7 分类，gene_id 键） | 图例改为 "including relaxed-constraint genes (n = 1,266; OR = 0.59 …; the strict gene-driven class alone gives OR = 0.55)" |
| 2 | CpG 括号混用分母 | 共同 391 基因：47.83% vs 49.10%（Δ=−1.28pp）；保留 180/192=93.8% | "(47.8% vs. 49.1% across the 391 genes common to both runs; 93.8% …)" |
| 3 | remainder of the universe 措辞 | 对照集=3,415（neutral 3,404 + dual 11，剔 GD±relaxed）；非 RD 全集 4,681 中位 28.8kb | 改为 "for genes in neither driven class (n = 3,415, unclassified plus dual-driven…)"；L50 Fisher 口径（RD vs 4,681 含 GD）经复算（EA OR 2.987 ✓）无需改 |
| 4 | φ 未校正 vs χ² p Yates 混用 | phi=0.052（未校正 χ²），chi2_p=4.60e-4（Yates） | 正文 L76 + 图例 L293 均标注 "φ = 0.052 from the uncorrected χ²; Yates-corrected (p/P) = 4.6 × 10⁻⁴" |
| 5 | adj OR 1.02–1.64 不含 Crohn 0.33 | 12 性状调整 OR：Crohn 0.331 唯一耗竭方向，其余 11 个 1.022–1.641 | 改为 "adjusted OR 0.33–1.64 … Crohn disease is the sole depletion-direction estimate at 0.33, the remaining 11 traits span 1.02–1.64" |
| 6 | baseline 复现差 0.509pp | 0.49873−0.49364=0.509pp | "within 0.5" → "within 0.51 percentage points" |

脚本：.workbuddy/tmp_r6fix8.py（7 处 count==1 断言全过）。补充：Supplementary_Text_v9.md 与 Main_Tables_v9.md 扫描确认无同类表述，无需改。


### Batch 8 交付状态 (2026-09-15 晚)
- 复核: general-purpose-30 聚焦复核 7/7 PASS（blind_verification_r6_batch8_recheck.md），无新错误；采纳其微调建议第 3 处改 "neither the gene-driven nor the regulation-driven class"（dual-driven 严格属两个 driven class）
- 交付链: 4 docx 重建→包同步→单 PDF 54 页（18/18 探针含 8 条 batch-8 新探针）→zip 6.15MB/180 文件 CRC OK→CRSHE 提交并推送 c9fab1a..5adf311
- 注: PDF 长短语探针因换行断裂需用短形（'0.51'、'regulation-driven class (n = 3,415'）
