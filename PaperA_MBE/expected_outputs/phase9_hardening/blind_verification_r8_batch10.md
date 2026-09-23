# R8 batch 10 全量盲验证报告（总核销）

验证员：general-purpose-30（独立盲验证）
日期：2026-09-24
方法：不信任任何转述，全部数字从 9 个权威源 JSON（`results/phase9_hardening/`）与原始文件独立重算；四文件（manuscript_english_v9.md、Main_Tables_v9.md、cover_letter_v1.md、supplementary/Supplementary_Text_v9.md）逐短语扫描；两张图 Read 看图核对。
范围：batch 10A（12 处）+ 10B（5 处）+ 10C（33 处）+ 摘要压缩，8 项清单逐项核销。

## 总评

**8/8 项 PASS。** 未发现任何 FAIL。两处观察性备注（不要求修改），一处 JSON 内部口径问题（不影响稿件）。

---

## 1. srv 新主口径 — PASS（全项）

| 检查点 | 稿件 | JSON 复算 | 结论 |
|---|---|---|---|
| 分层加权主估计 | 29.2% [26.4–32.1] | stratified=0.291922, CI [0.26442, 0.32093] → 29.2% [26.4, 32.1] | PASS |
| 绝对数 ≈1,450 | ≈1,450 | 0.291922 × 4974 = 1452.0 | PASS |
| 下夹口径仅作下限 | 28.0% / 25.3% 标注为 lower bound | ratio=0.280220, absdiff=0.253245 | PASS |
| flip 表 | 44/194 flip-out vs 10/199 flip-in | 0.22680 / 0.05025 | PASS |
| retention | 77.3% | 150/194 = 0.77320 | PASS |
| GD 层 srv 率 | 24.3% | 35/144 = 0.24306 | PASS |
| srvYes GD vs RD | 19.6 [17.6–21.7] vs 6.2 [5.9–6.7], factor 3.2 | 19.620 [17.559, 21.713] / 6.210 [5.891, 6.705], ratio = 3.1595 → 3.2 | PASS |
| 旧 1,394 零残留 | — | "1,394" 四文件 0 命中 | PASS |
| 9G FAIL-1 修复 | "49.4% to 40.7% on the 393 comparable genes" | rate_srvYes_on_common = 160/393 = 0.407125 → 40.7%（共同框架正确） | PASS（修复确认） |
| 算术自检 | — | 0.339767×(151/196) + 0.660233×(9/197) = 0.29192 闭合；151/196=0.770408、9/197=0.045685 与 JSON 一致 | PASS |

## 2. GWAS 降级 — PASS（全项）

| 检查点 | 稿件 | JSON 复算 | 结论 |
|---|---|---|---|
| 主口径 0/48 | 0/48 survive length-matched permutation | gwas_length_matched_perm_r8.json: full 0/48 | PASS |
| LOO-brain EA | q=0.048（唯一幸存格） | LOO-brain 1/48, EA q=0.04795 → 0.048 | PASS |
| span×CDS 双残差 | 0/12（×2 变体） | 0/12 | PASS |
| 最弱名义格 | Bipolar 0.017 / EA 0.020 | 0.01698 / 0.01998 | PASS |
| pooled z | Jaccard z=2.84, p=0.004；overlap z=2.06, p=0.039 | z=2.84494 p=0.0044419；z=2.06217 p=0.03919 | PASS |
| 148 下采样 | 中位 6/8 | rd_downsample_decay_r8.json: median 6.0（P(≤5)=0.48） | PASS |
| LOO-tau 位置 | 5/8 居 48 分位 | 5/8 与下采样分布中位邻域一致 | PASS |
| "family-level residual signal" 六处一致 | 摘要 @1795 / integrative answer @47322 / Limitations @73880 / 结语 @96111 前 / 图例区 | 稿件 5 处 + response_round6.md 1 处 = 6 处，措辞一致 | PASS |

## 3. 验证链重新定位 — PASS（全项）

| 检查点 | 稿件 | JSON 复算 | 结论 |
|---|---|---|---|
| HAR 长度置换 | p=0.001 | hcondel_length_perm_r8.json: 0.000999（×2 处一致） | PASS |
| hCONDEL 置换链 | 0.22 / 0.69 / 0.032 / 0.035 / 0.076 / 0.126 | 0.22178 / 0.68831 / 0.03197 / 0.03497 / 0.07592（span-resid）/ 0.126 | PASS |
| LOO-nc 限定形式 | OR=4.67, p=0.032 | 独立 Fisher 复算 4.6686 / 0.03197 | PASS |
| dual-support 限定形式 | OR=11.38, p=0.035 | 11.3792 / 0.03497 | PASS |
| 异质性 | OR=6.23, p=1.4e-3 | 独立复算 6.2308 / 0.0014354 | PASS |
| joint 主验证 | HAR 3.86 [3.02–4.95]；hCONDEL 3.05 [1.90–4.88] | 独立复算 3.8641 (p=7.15e-23) / 3.0481 (p=2.24e-5) | PASS |
| joint 基数 | 403 / 1,256 / 21 | 403/1256/21 | PASS |
| 跨验证交集 | 215 / 188 / 78 | 215/188/78 | PASS |
| 跨度 | 105.8 / 72.8 kb | 105793 / 72839 bp → 105.8 / 72.8 kb | PASS |
| "primary independent validation" 残留 | — | 四文件 0 命中 | PASS |
| "93% post-hoc power" 残留 | — | "post-hoc power" 与 "93% power" 均 0 命中 | PASS |

## 4. P0-2 变体 B — PASS（全项）

| 检查点 | 稿件 | JSON 复算（rds_spanresid_r8.json） | 结论 |
|---|---|---|---|
| RD 数 | 237 | 237（主残差瑕疵版 1675、变体 C 595 均未入稿） | PASS |
| span ratio | 4.2×→2.8× | v7 4.2114 → variant B 2.7708（87894/31721） | PASS |
| Spearman ρ | −0.008 | −0.0077715 → −0.008 | PASS |
| HAR | 3.52, p=3.5e-13 | 独立复算 3.5210 / 3.49e-13 | PASS |
| hCONDEL | 2.97, p=7.8e-4 | 独立复算 2.9718 / 7.78e-4 | PASS |
| SynGO | 7/7 | 7/7 | PASS |
| GO_BP | 12/18 | 12/18 | PASS |
| 置换敏感性 | 0.126 / 0.076 / 0.001 | rds_spanresid_perm_r8.json: 0.12587 / 0.07592 / 0.000999 | PASS |
| 瑕疵中间版残留 | — | "1,675" 与 RD=223 四文件 0 命中 | PASS |

## 5. 残留扫描 — PASS

16 个禁用串在 manuscript / Main_Tables / cover letter / SI 全部 0 命中：
"12.2%"、"68.5%"、"BH exact"、"empirical FDR"、"did not converge"、"phylogeny is well established"、"primary independent validation"、"post-hoc power"、"93% power"、"1,675"、"1,394"、"re-detection"、"223"、"49.4% to 41.0%"（旧 9G 错版）、"0.543"（旧 bipolar q）、"258 词" 旧摘要特征。
L106 GARD 新措辞与 SI S6 名称/坐标新口径均已落稿（batch-9 备注 N3 建议已采纳）。

## 6. 图-文-表三方一致 — PASS

- **Figure 4**（Read 看图）：panel a/b 末行 = LOO-brain-tau (joint)，HAR 3.86 [3.02–4.95]、hCONDEL 3.05 [1.90–4.88]，与正文/Table 3 完全一致；panel c 无 "excludes HAR re-detection" 字样；panel f RD 3.19 [1.90–5.37]、GD 0.59 [0.36–0.97] p=0.04 与 census 一致。
- **Table 3 joint 行**：1,256† / 403 / 3.86 (3.02–4.95) / 3.05 (1.90–4.88) ✓；脚注 215/188/78、105.8 kb、72.8 kb ✓。
- **Figure 7c**（Read 看图）：红色长度警示 "RD 8/8 FDR-sig. — but uncorrected for gene length (see text); 0/48 survive length-matched permutation" ✓。

## 7. 摘要 — PASS

词数三口径：249 / 248 / 243，均 ≤250。摘要内数字（29.2% [26.4–32.1]、≈1,450、0/48、3.86/3.05、12.1%）与正文及 JSON 自洽；604/4974=12.143% → "12.1%" 已正确落稿（batch-9 FAIL F2 修复确认）。

## 8. 算术自检 — PASS

- 分层估计量：0.339767×(151/196) + 0.660233×(9/197) = 0.29192 ✓
- 0.291922 × 4974 = 1452.0 ✓
- 19.6203 / 6.2100 = 3.159 → "3.2-fold" ✓
- 44/194=0.2268、10/199=0.0503、150/194=0.7732、35/144=0.2431 全部与 JSON 逐位一致 ✓
- 4 个 Fisher 独立复算（joint HAR/hCONDEL、异质性、变体 B HAR/hCONDEL）全部与 JSON 逐位一致 ✓

---

## 观察性备注（非 FAIL，不要求修改）

- **O1**：变体 B 的 "217/293 = 74.1% 保留率"（JSON: retained 217 of 293）数字与 JSON 一致但稿件未写明该分数表述。属取舍问题，非错误。
- **O2**：`audit_r8.json` 记录 tier4 JSON 仍含旧 n=183 数（CSV 已为 n=112）——JSON 内部滞后，不影响稿件任何数字。

## 结论

R8 batch 10 全部修复落稿**核销通过**：8/8 项 PASS、0 FAIL。此前 batch-9 的 2 处 FAIL（16.1%、12.2%）与 batch-9G 的 1 处 FAIL（srv 框架混用）均已确认修复。稿件数字与 9 个源 JSON 完全一致，可进入下一轮。
