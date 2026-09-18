# Batch 8 Recheck — 七处口径修改交叉验证

**验证人**: general-purpose-30（R6 盲验证员）
**对象**: `results/paper/manuscript_english_v9.md`（batch 8 修改后版本）
**方法**: 逐处正则核验新文本存在/旧文本消失 + 与本验证员 R6 独立复算值（`scripts/_blind_verify_r6/outputs/`）逐位对照 + 上下文语法与相邻数字一致性检查。只读验证。

## 逐项结果

| # | 修改 | 独立复算对照 | 判定 |
|---|---|---|---|
| 1 | baseline 复现差 "within 0.51 percentage points" | 我方复算：production-on-same-genes 0.49873 − baseline 0.49364 = **0.509pp → 0.51** ✓；旧文 "within 0.5 percentage points" 0 残留；相邻 "p-value correlation 0.97, 98.5% per-gene concordance" 与复算（r=0.9728, 98.47%）一致 | **PASS** |
| 2 | CpG 分母 "(47.8% vs. 49.1% across the 391 genes common to both runs; 93.8% …retained)" | 我方复算：共同 391 基因 cpg 187/391=**47.83%**、baseline 192/391=**49.10%**（差 1.28pp→−1.3pp ✓）、保留 180/192=**93.75%→93.8%** ✓；"49.4%" 全文 0 残留；口径统一后括号内自洽 | **PASS** |
| 3 | 长度对照集 "for genes in neither driven class (n = 3,415, unclassified plus dual-driven; 4.7-fold; Wilcoxon P = 1.5 × 10⁻⁷¹)" | 我方复算（与分析代码定义一致）：对照集 = neutral 3,404 + dual 11 = **3,415** ✓，中位 **28,648→28.6 kb** ✓，ratio **4.663→4.7** ✓，MWU **1.453e-71→1.5 × 10⁻⁷¹** ✓；旧文 "for the remainder of the universe" 在该处 0 残留（余下 2 处 "remainder of the universe" 属 Methods Fisher 检验描述与 Fig 5e Selectome 句，与长度对照无关，合法保留） | **PASS**（附措辞注：dual-driven 严格说同时属于两个 driven class，"neither driven class" 对其略欠精确，但括号内 "unclassified plus dual-driven" 已显式定义集合构成，n=3,415 无歧义） |
| 4 | 调整 OR 范围 "(adjusted OR 0.33–1.64, all q > 0.44; Crohn disease is the sole depletion-direction estimate at 0.33, the remaining 11 traits span 1.02–1.64; …)" | 我方复算：Crohn OR_adj=**0.331→0.33** ✓（唯一 <1 者 ✓）；其余 11 性状 **1.022（Intelligence）–1.641（Bipolar）→1.02–1.64** ✓；EA **1.219→1.22**、height **1.063→1.06**、SCZ **1.272→1.27** ✓；min q=0.4479>0.44 ✓；单独的旧文 "adjusted OR 1.02–1.64" 0 残留（新文内 "span 1.02–1.64" 为有意保留） | **PASS** |
| 5 | 正文 φ 标注 "φ = 0.052 from the uncorrected χ², Yates-corrected p = 4.6 × 10⁻⁴" | 我方复算：未校正 χ²=13.43 → **φ=0.0520→0.052** ✓；Yates 校正 χ²=12.27 → **p=4.60e-4** ✓；旧式 "φ = 0.052, χ² p = 4.6 × 10⁻⁴" 0 残留；同句 19.6%、22/476=4.6%、9/90 OR=1.80 P=0.081 均与复算一致 | **PASS** |
| 6 | Fig 4c 图例 "φ = 0.052 from the uncorrected χ²; Yates-corrected P = 4.6 × 10⁻⁴" | 与第 5 项同源数字，均与复算一致；图例内 22/112 (19.6%)、9/22 OR=11.38 P=1.9e-6、72/454 OR=3.67 P=4.3e-16 亦与复算一致 | **PASS** |
| 7 | Fig 4f 图例 "gene-driven genes including relaxed-constraint genes (n = 1,266; OR = 0.59, 95% CI 0.36–0.97, two-sided P = 0.037; the strict gene-driven class alone gives OR = 0.55)" | 我方复算：GD+relaxed k=19/1,266 **OR=0.5923→0.59**、CI (0.360, 0.974)→**0.36–0.97**、双侧 **p=0.0371→0.037** ✓；strict GD k=17/1,214 **OR=0.5479→0.55** ✓；旧文无口径限定的版本 0 残留 | **PASS** |

## 附带核查

- **未见新错误**：7 处修改句语法通顺、无重复表述；相邻数字（133.6 kb、ρ=0.47、|ρ|≤0.09、OR 3.19/CI 1.90–5.37/P 7.8e-5、1.23 vs 1.06 P=0.24 等）均与 R6 独立复算一致。
- "OR = 0.55" 全文另有 3 处命中，均为 Shao et al. 交叉研究句（OR=0.55, P=0.21）的既有内容，与本修改无关、数值巧合但语境独立，不构成矛盾。
- 顺带确认：同段落 "median 1.3–4.1-fold"（此前 R6 报告 FAIL #4）已同步修复，与复算（1.29–4.12）一致。
- "Yates" 一词全文 3 处：2 处为本修改（正文+图例），1 处为引文作者名（Frankish… Yates A），合法。

## 结论

**7/7 全部 PASS：新文本均在位、旧文本均已清除、全部数字与本验证员的独立复算逐位一致，且未引入任何新错误。**（第 3 项附一条可选的措辞微调建议：dual-driven 严格而言同属两个 driven class，如需绝对精确可将 "neither driven class" 改为 "neither the gene-driven nor the regulation-driven class"，但现有括号定义已无歧义。）
