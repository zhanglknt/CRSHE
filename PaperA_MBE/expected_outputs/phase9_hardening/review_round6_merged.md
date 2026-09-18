# MBE 模拟审稿 Round 6 合并报告（v9 终版，2026-09-17）

**稿件**：manuscript_english_v9.md（R5 全部 P0/P1 修复 + 盲验证 7/7 PASS 之后的版本）
**全新四人专家团**（与 R5 四位无重叠）：灵长类比较基因组学 / 群体数量遗传学 / 分子演化方法学 / MBE 副编辑视角

## 一、总评分

| 审稿人 | 总分 | 推荐 | 分项（N/R/S/C） |
|--------|------|------|-----------------|
| 灵长类比较基因组学 | 7.5 | Minor Revision | 7.5 / 7 / 7 / 8.5 |
| 群体数量遗传学 | 7.0 | Major Revision | 7 / 7 / 7 / 8 |
| 分子演化方法学 | 7.5 | Major Revision | —（详见单份） |
| MBE 副编辑 | 8.0 | **Send to review**（预期 minor-to-moderate 后接收） | 契合 8.5 / 叙事 8 / 图表 8 / 包完整 7 |
| **均值** | **7.5** | 2 Major + 1 Minor + 1 送审 | |

**六轮轨迹**：4.8（3/3 Major）→ 6.7（3/3 Minor）→ 7.2（3/3 Accept）→ 7.3（3/3 Accept）→ 7.5（4/4 Minor）→ **7.5（新团更严：2 Major + 1 Minor + 1 送审）**

单份报告：`review_round6_{primate_genomics,popgen,molevol_methods,editor}.md`

## 二、P0 合并清单（去重后 6 项）

### 编码侧检测校准（2 项，molevol）
- **P0-1 有界双峰 p 值 → BH 偏差方向未定**：39.4% p 值恰在 0.5；若零分布为 U(0,0.5)，BH 约保守 2 倍，Storey π₀=1.0 是退化检查。修法：中性模拟校准或有界支撑校正
- **P0-2 34% 上限只披露未测量**：要求 300–500 基因分层子集上 --srv Yes / GUIDANCE 过滤 / CpG 屏蔽三组对照重跑；GARD 全不收敛本身是比对质量红旗，稿件未追问

### GWAS 锚定（2 项，popgen）
- **P0-3 富集未控制基因长度/映射混杂**：MAPPED_GENE 多按距离指派、长基因更易被映射，RDS 的 nc 分量与基因跨度相关；coding 侧做了 CDS 残差化而 GWAS 侧零处理。修法：logistic 回归（RD ~ GWAS + log10 基因长度）或分层 Fisher 重检 12 性状；height 富集可能即此通用机制的表现
- **P0-4 "strongest and most consistent" 无正式组间检验**：intelligence OR 1.95 < height 2.2–2.3；LOO 下神经性状 6/8、5/8 而阴性对照三分类全稳。修法：组间 log OR 比较，否则降级表述（含 Abstract）

### 数据与谱系（2 项，primate）
- **P0-5 HAR/hCONDEL 集合来源未注明且与引文不符**：Pollard 202 vs 文中 3,171；McLean 510 vs 583——hCONDEL 是 "primary independent validation"，必须交代确切来源/基因组版本/lift-over
- **P0-6 ILS 结构性混杂零讨论**：人–黑猩猩–大猩猩节点约 25–30% 基因座谱系不一致，对人类末端支前景的 BUSTED 是结构性混杂。修法：Discussion 定量讨论 + 高 ILS 区域敏感性

## 三、P1 合并（去重后 13 项）

| # | 问题 | 来源 |
|---|------|------|
| 1 | MAPPED_GENE 分词敏感性无量化（locus 去重 + 最近基因两套重跑） | popgen |
| 2 | 8 性状遗传相关，"all eight" 修辞高估（trait set 重叠未报） | popgen |
| 3 | Crohn/LDL null 的功效与 CI 未报告 | popgen |
| 4 | 缺 GWAS anchoring 主表（n/OR/CI/q × LOO 变体） | popgen |
| 5 | variant-density 替代假说（现在变异密度 vs 历史选择）未讨论 | popgen |
| 6 | 分支长度 "calibrated tree" 表述 | primate |
| 7 | 无外群（marmoset/squirrel monkey 作为新世界猴外群的处理） | primate |
| 8 | Selectome 94→24 映射损耗解释 | primate |
| 9 | 子集含于参照的独立性（HAR 子集 ⊂ universe 的检验独立性） | molevol |
| 10 | RELAX K>1 "intensified" 方向歧义 | molevol |
| 11 | 权重 a priori 声明不可验证 | molevol |
| 12 | height/T2D 负对照对 specificity 表述的约束 | molevol |
| 13 | 编辑可选：Fig. 6 降补充图、三处占位符、两处覆盖率论述合并 | editor |

**已核销（误报）**：molevol P1-1 "tier4 CSV 与稿件不符" —— 已验证三处副本（results/包/CRSHE）逐字节一致，CSV 为比例单位（0.0258=2.58pp），稿件为百分点（+2.6），数值精确对应；系单位误读，非文件错误（可选改进：CSV 列名标注单位）。

## 四、共识判读（四方汇聚）

**共识优点**：
- 循环性卫生"教科书级"（LOO 全矩阵 / caMPRA 校正后不显著主动撤销 / Selectome 零结果如实报告）——molevol 与 primate 独立同词
- 自我审计诚实度罕见（12 条量化 Limitations、Shao 分母纠错、阴性对照自曝）
- 可复现性六轮最佳：molevol 独立复算 8 组数字全部精确一致
- hCONDEL 独立验证设计获 popgen 点名为最大亮点
- 编辑确认 **MBE 为最优定位**（GBE/PLOS Gen 受众匹配度低；Nat Comms 条件式结论桌拒中高风险）

**汇聚攻击点（三方以上汇聚）**：
1. **标题 "Dominates by Prevalence" 继承未校准的检测层**（molevol P0-1/2 + primate strict/permissive + editor 风险项）——若分层重跑检出率缩水，标题与摘要的流行度表述需改
2. **GWAS 锚定是标题级卖点却是最弱证据链**（popgen 双 P0 + molevol P1-12 + editor 风险项）——特异性叙事面临"RD 对一切充分 empowered 的复杂性状普遍富集"的退化风险
3. 标题-结论张力（全覆盖反事实 293→1,906 可逆转）为外审最可能攻击点（editor + R5 已有）

## 五、与 R5 对照

| | R5（v9 首审） | R6（v9 修复后，新团） |
|---|---|---|
| 均分 | 7.5 | 7.5 |
| 判定 | 4/4 Minor | 2 Major + 1 Minor + 1 送审 |
| P0 性质 | 修复执行类（LOO 未重跑、Selectome 作废、口径矛盾） | **分析深化类**（混杂校正、模拟校准、来源追溯、ILS）——无需推翻既有结果，但要求新的对照分析 |
| 数字可信度 | R3 抽查全过 + 盲验证 7/7 | molevol 复算 8 组全过（唯一不符项为误报） |

**解读**：R6 分数持平但判定更严，符合"每轮换更挑剔视角"的预期；P0 全部是"补对照/补校准"类而非"结果错误"类——既有结论框架未被推翻，但投稿前需评估是否执行（P0-3/P0-4 为 GWAS 侧可本地快速执行；P0-1/P0-2 需 WSL 重跑 BUSTED 子集；P0-5 为文献考证；P0-6 为讨论+敏感性）。
