# 盲验证报告 R8：GARD 清洁比对对照数字独立复算

**验证员**：verify-gard-r8（独立盲验证，不信任任何中间报告）
**日期**：2026-09-24
**唯一输入**：
- `results/phase9_hardening/gard_clean_control_r8.json`（per_gene 50 条，逐字段重算）
- `results/phase7_gene_vs_regulation/phase7g_classification_v7/gene_classification_v7.csv`（4,974 行）
**被验文本**：`results/paper/manuscript_english_v9.md`（L32/L106/L142/L176）、`results/paper/supplementary/Supplementary_Text_v9.md`（S2）
**方法**：从 per_gene 逐条记录独立重算全部统计量（不复用 JSON 顶层 summary；summary 仅作交叉核对，10 项全部吻合），再与稿件逐字比对。

---

## 基础复算结果（供各声明引用）

| 量 | 复算值 |
|---|---|
| verdict 分布 | diff_site 35 / same_site_healthy 13 / same_site_patho 2 / vanished 0 |
| 旧病态树（old_pathological_tree flag） | 44（与 old_max_bl>100 完全一致） |
| 清洁后病态树（new_rc4_max_bl>100） | 8（其中 2 棵由健康转为病态：ENSG00000092529、ENSG00000166866） |
| 旧病态→新健康 | 38 / 44 = 86.36% |
| median old maxBL（全 50） | **2,244.145** |
| median old maxBL（仅 44 病态子集） | **2,599.735** |
| median new rc4 maxBL（全 50） | **1.57** |
| median new rc4 maxBL（仅 38 修复子集） | **1.435** |
| median new rc4 maxBL（13 同位点健康子集） | **1.58** |
| 13 子集 Δc-AIC（new_rc4_delta） | min 30.7 / max 2,717.7 / median 216.7 |
| 35 个 diff_site 移位幅度（\|new_rc4_bp−old_bp\|/clean_kept_cols，已独立重算字段值，一致） | min 0.0985（ENSG00000183621，66/670）与 0.0994（ENSG00000137727，78/785）；max 0.7992（ENSG00000168961，195/244）。33/35 落在 [0.10, 0.80]，其余 2 个为 9.85%/9.94% |
| rc=2 | 49/50 有断点（ENSG00000151461 为 null）；全 49 中 43 个与 rc4 位点 ±50 codon 一致；**13 子集 13/13 一致**（且 13 个 rc4、rc2 位点均在旧位点 ±50 codon 内） |
| 13 子集 classification_v7（CSV 独立核对） | 13/13 全为 "gene-driven"（strict）；50 筛查集 = 49 strict + 1 relaxed；CSV 全宇宙 gene-driven(strict) = **1,214**；13/1,214 = **1.071%** |
| JSON 顶层 summary 交叉核对 | 10/10 项与独立复算完全一致 |

---

## 逐条声明核对

### 声明 1（Methods，L32）
> "no breakpoint vanished outright, but 70% shifted position and 38 of 44 pathological trees were repaired — while 13 of 50 genes retain a same-site breakpoint on healthy trees"

| 子项 | 稿件 | 复算 | 判定 |
|---|---|---|---|
| vanished | 0 | 0 | PASS |
| shifted | 70% | 35/50 = 70.0% | PASS |
| repaired | 38 of 44 | 38/44 | PASS |
| same-site healthy | 13 of 50 | 13/50 | PASS |

**判定：PASS**

### 声明 2（Results，L106 括号句）
> "70% shifted position, median maximum branch length 2,244 → 1.57 substitutions per site — while 13 of 50 genes retain a same-site breakpoint"

| 子项 | 稿件 | 复算 | 判定 |
|---|---|---|---|
| shifted | 70% | 70.0% | PASS |
| median maxBL 2,244 → 1.57 | 全筛查集旧→新中位数（箭头式表述，未绑定子集） | 全 50：2,244.145 → 1.57 | PASS |
| same-site | 13 of 50 | 13/50 | PASS |

**判定：PASS**

### 声明 3（Limitations (3)，L176）
> "38 of 44 pathological partition trees **(median maximum branch length 2,244 substitutions per site)** became healthy on clean alignments **(median 1.57)**, and 35 of 50 original breakpoint positions (70%) shifted by 10–80% of alignment length"
> "13 of 50 genes (26%) retain a breakpoint at the same site ... (Δc-AIC 30.7–2,718, median 217; rc = 2 concordant in 13 of 13)"
> "13 of 1,214 gene-driven genes, ~1%"

| 子项 | 稿件 | 复算 | 判定 |
|---|---|---|---|
| 38 of 44 | 38/44 | 38/44 | PASS |
| **44 病态树中位数 2,244** | 括号语法上修饰 "44 pathological partition trees" | **44 棵病态树子集中位数 = 2,599.7**；2,244.145 是全 50 树（含 6 棵健康）的中位数 | **FAIL** |
| **修复后 (median 1.57)** | 括号语法上修饰 38 棵修复树 | **38 棵修复树子集中位数 = 1.435**；1.57 是全 50 树清洁后中位数 | **FAIL** |
| shifted 70%，幅度 10–80% | 35/50，范围 10–80% | 70.0%；实测 9.85–79.9%（两个端点值 9.85%/9.94% 与 79.9% 四舍五入即为 10%/80%） | PASS（舍入边界，见备注 a） |
| 13 of 50 (26%) | 26% | 26.0% | PASS |
| Δc-AIC 30.7–2,718, median 217 | — | 30.7 / 2,717.7 / 216.7 | PASS |
| rc=2 concordant 13 of 13 | — | 13/13（±50 codon） | PASS |
| 13 of 1,214，~1% | — | CSV strict GD = 1,214；13/1,214 = 1.071% | PASS |

**判定：FAIL（中位数挂载对象错误）**。2,244 与 1.57 本身是真数字（全 50 树前/后中位数，Results L106 与补充文本箭头式表述均正确），但在 Limitations 句中括号直接修饰 "44 pathological partition trees" 与隐含的 38 棵修复树，按此自然语法解读两个子集中位数应为 **2,599.7 → 1.435**。修复建议（任选其一）：
1. 改为与 Results/补充文本一致的箭头式："38 of 44 pathological partition trees became healthy on clean alignments (screen-wide median maximum branch length 2,244 → 1.57 substitutions per site)"；或
2. 改用子集真实中位数："(median 2,600) became healthy (median 1.44)"。

### 声明 4（Discussion，L142 插入句）
> "13 of the 50 screened genes — all in the gene-driven class — retain a same-site breakpoint ... this affects ~1% of the gene-driven class"

| 子项 | 稿件 | 复算 | 判定 |
|---|---|---|---|
| 13 of 50 | — | 13/50 | PASS |
| all gene-driven | — | CSV 核对 13/13 全为 gene-driven(strict) | PASS |
| ~1% of gene-driven class | — | 13/1,214 = 1.071% | PASS |

**判定：PASS**

### 声明 5（补充文本 S2，第 (i) 层）
> "38 of the 44 pathological partition trees (86%) became healthy (median maximum branch length 2,244 → 1.57; only 8 of 50 exceed 100 substitutions per site after cleaning), 35 of 50 original breakpoint positions (70%) shifted by 10–80%"

| 子项 | 稿件 | 复算 | 判定 |
|---|---|---|---|
| 38/44 (86%) | 86% | 86.36% | PASS |
| 2,244 → 1.57 | 箭头式全筛查集表述 | 2,244.145 → 1.57（全 50） | PASS |
| 8 of 50 exceed 100 | — | 恰好 8 个 new_rc4_max_bl>100 | PASS |
| 70% shifted by 10–80% | — | 70.0%；实测 9.85–79.9% | PASS（同备注 a） |

**判定：PASS**

### 声明 6（补充文本 S2，第 (ii) 层）
> "No breakpoint vanished outright (0 of 50), 13 of 50 genes (26%) retain a breakpoint at the same site (±50 codons) ... healthy partition trees (median maximum branch length 1.58), strong model support (Δc-AIC 30.7–2,718, median 217), rc = 2 concordance in 13 of 13 cases; 2 further genes retain a same-site call with still-pathological trees"

| 子项 | 稿件 | 复算 | 判定 |
|---|---|---|---|
| 0 of 50 vanished | — | 0 | PASS |
| 13/50 (26%) same site ±50 codons | — | 13/50 = 26.0%；13 个 rc4 位点距旧位点 0–48 codon，全在 ±50 内 | PASS |
| median maxBL 1.58 | 13 子集 | 13 子集中位数 = 1.58 | PASS |
| Δc-AIC 30.7–2,718, median 217 | — | 30.7 / 2,717.7 / 216.7 | PASS |
| rc=2 concordance 13/13 | — | 13/13（rc2 位点均在 rc4 ±50 codon 内） | PASS |
| 2 further same-site patho | — | 2（ENSG00000134899 maxBL 2,603.66；ENSG00000176658 maxBL 103.72） | PASS |

**判定：PASS**

### 声明 7（补充文本 S2，第 (ii)/(iii) 层）
> "all 13 fall within the gene-driven class (13 of 1,214, ~1.1%)"；"(iii) The persistent subset is ~1.1% of the gene-driven class"

| 子项 | 稿件 | 复算 | 判定 |
|---|---|---|---|
| all 13 gene-driven | — | 13/13 strict GD（CSV 独立核对） | PASS |
| 13 of 1,214 | — | CSV 全宇宙 strict GD = 1,214 | PASS |
| ~1.1% | — | 13/1,214 = 1.071% ≈ 1.1% | PASS |

**判定：PASS**

---

## 备注

a. **移位幅度 "10–80%"**：35 个 diff_site 基因实测移位范围为 9.85%–79.9%（最小两个 9.85%、9.94%）。两个端点四舍五入后即为 10% 与 80%，判 PASS，但知会此舍入边界。
b. **rc=2 全筛查集一致性**：49 个完成 rc=2 的基因中 43 个与 rc4 位点 ±50 codon 一致（6 个不一致：ENSG00000012822/000092529/00139131/00139567/00166866/00177119）。稿件未对全筛查集 rc2 一致性作声明，仅声明 13 子集 13/13（属实），无需修改。
c. **分母 1,214 vs 1,266**：JSON 的 genome_wide_baseline 把 "gene-driven (relaxed)" 计入得 1,266；稿件用 strict GD = 1,214，与 CSV 独立计数一致，且 13 子集全为 strict GD，分母选择正确。
d. **2 棵逆向转化**（旧健康→新病态：ENSG00000092529、ENSG00000166866）使"旧病态→新健康 = 38"与"新病态总数 = 8"自洽（44−6=38；6 旧病态未修复 + 2 新增 = 8）。稿件未提逆向转化，但所声明数字均与此自洽，无矛盾。

---

## 总体判定

**6 / 7 条声明 PASS；声明 3（Limitations (3)，L176）FAIL** —— "(median maximum branch length 2,244 ...)" 与 "(median 1.57)" 两处括号在语法上绑定 44 病态树/38 修复树子集，但数字是全 50 树的中位数；子集真实中位数为 2,599.7 与 1.435。

**总体：FAIL（不可直接交付）**。仅需对 Limitations (3) 该句做一处文字修复（改为箭头式全筛查集表述或改用子集中位数 2,600 → 1.44），其余全部数字经得起独立复算。修复后无需重跑数据，复核该句即可转为 PASS。

---

## 修复复核（2026-09-24，team-lead 按方案 1 修复后）

**落盘核实**（不信转述，直接 grep `manuscript_english_v9.md`）：
- 新文本 "screen-wide median maximum branch length 2,244 → 1.57 substitutions per site" 在 L176 出现 **1 次**；旧文本 "partition trees (median maximum branch length ... (median 1.57))" **0 处残留**。
- 修复后全句逐字提取："...38 of 44 pathological partition trees became healthy on clean alignments (screen-wide median maximum branch length 2,244 → 1.57 substitutions per site), and 35 of 50 original breakpoint positions (70%) shifted by 10–80% of alignment length..."，与 team-lead 转述完全一致。

**数字复核**：
- "screen-wide median maximum branch length 2,244 → 1.57" 为全筛查集箭头式表述：全 50 旧树中位数 2,244.145 → 全 50 清洁后中位数 1.57，与独立复算一致，PASS。
- 表述与 Results L106（"median maximum branch length 2,244 → 1.57 substitutions per site"）及 Supplementary S2（"median maximum branch length 2,244 → 1.57"）完全一致，三处口径统一。
- 句内其余子项（38 of 44、35 of 50、70%、10–80%）前轮已判 PASS，未受改动。

**声明 3 复核判定：PASS。**

## 最终判定

**7 / 7 条声明全部 PASS。总体：PASS（可交付）。**

GARD 清洁比对对照写入稿件（主文 Methods/Results/Discussion/Limitations + 补充文本 S2 三层叙事）的全部数字均经得起从 `gard_clean_control_r8.json` per_gene 记录与 `gene_classification_v7.csv` 的独立重算。
