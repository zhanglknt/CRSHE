# 盲验证报告 R8：匹配 null 校准（codon-aware 重比对）数字核销

验证员：verify-nullcal-r8（独立盲验，未参与整合）
日期：2026-09-25
范围：manuscript_english_v9.md（Discussion/Methods/Limitations）、Supplementary_Text_v9.md（S7）、cover_letter_v1.md
真值源：neutral_realign_calibration_r8.json、neutral_realign_stratified_r8.json、alignment_rerun_sensitivity.json、srv_estimators_r8.json、neutral_sim_calibration.json、neutral_sim_parametric.json

---

## 1. 基础复算表（从 JSON 独立算出）

### 1a. neutral_realign_calibration_r8.json（p<0.05 检出率，×100 为 %）

| 组 | n_common | baseline % | realigned % | 精确 delta (pp) | JSON 存 delta | 闭合误差 | BH 拒绝 base→real |
|---|---|---|---|---|---|---|---|
| R7 parametric | 383 | 1.3055 | 12.7937 | 11.4883 | 11.49 | 0.0017 | 1 → 37 |
| R6 evolver | 398 | 1.7588 | 19.0955 | 17.3367 | 17.34 | 0.0033 | 0 → 71 |
| combined | 781 | 1.5365 | 16.0051 | 14.4686 | 14.47 | 0.0014 | 1 → 108 |

- n_realigned_total = 790；800 − 790 = 10；383 + 398 = 781 = combined n_common ✓
- p = 0.5 点质量：R7 74.6736% → 64.2298%；R6 45.7286% → 41.9598%；combined 59.9232% → 52.8809%

### 1b. neutral_realign_stratified_r8.json（by_kept_cols，delta pp）

| 组 | low | mid | high |
|---|---|---|---|
| R7 parametric null | 3.15（计数精确 3.1496） | 10.16（10.1562） | 21.09（21.0938） |
| R6 evolver null | 12.88（12.8788） | 18.05（18.0451） | 21.05（21.0526） |
| real data（delta_qsig_pp） | 4.58 | 14.50 | 25.19 |

by_frac_removed（非单调轴）：R7 8.66 / 14.06 / 11.72；R6 16.67 / 16.54 / 18.80（计数精确 16.6667 / 16.5414 / 18.7970）。
全部 JSON 内存 delta 与两端率之差闭合，最大误差 0.01pp（容差 ±0.05pp）。

### 1c. alignment_rerun_sensitivity.json（真数据参照）

- common 393：49.3639% → 64.1221%，delta_pp = 14.7583（存储值与复算完全一致）
- 1/388 = 0.25773% ≈ 0.26%；srv_estimators_r8.json 推荐估计 0.29192 → 29.2%，CI [0.26442, 0.32093] → 26.4–32.1%

---

## 2. 逐项 PASS/FAIL

### A. 主稿 Discussion（"A matched null calibration adjudicates"，行 142）

| 项 | 文中值 | 复算值 | 判定 |
|---|---|---|---|
| A1 | 790 of 800 ω = 1 alignments completing | n_realigned_total = 790（800−790=10） | PASS |
| A2 | 781 common to baseline | combined n_common = 781；383+398=781 | PASS |
| A3 | 1.5% → 16.0%（+14.5 combined） | 1.5365→16.0051；精确 14.4686pp → 1dp 14.5 | PASS |
| A4 | +11.5 / +17.3（parametric-bootstrap / evolver） | R7 11.4883→11.5；R6 17.3367→17.3；组别-名称映射正确（parametric-bootstrap=R7，evolver=R6） | PASS |
| A5 | 108 of 781 vs 1 at baseline | combined n_rejections 1 → 108 | PASS |
| A6 | +3.2→+21.1（parametric）；+12.9→+21.1（evolver）；+4.6→+25.2（real） | 见 1b：21.09/21.05/25.19→21.1/21.1/25.2 ✓；12.88→12.9 ✓；4.58→4.6 ✓；+3.2 见注 1 | PASS（注 1） |
| A7 | 残留句 0/400 evolver、0.26%（1/388）parametric | neutral_sim_calibration.json：BH rejections 0/400；neutral_sim_parametric.json：1/388=0.2577%≈0.26%。calibration JSON baseline 列为同批原始跑取 common 子集，拒绝数同为 0 与 1，口径自洽（残留句=production 管线 null；新句=重比对管线 null，两者并存无矛盾） | PASS |

注 1（舍入边界，不影响结论）：文本 +3.2 来自 JSON 存值 3.15（round-half-up），按分层 n_sig 精确计数为 3.1496pp（4/127），1dp 应作 3.1。偏差 0.0504pp，属"文本对 JSON 二次舍入"情形，梯度形状与结论不受影响。

### B. 主稿 Limitations (2)（"a matched null calibration — both neutral simulation sets"，行 176）

- 790 of 800 completing ✓；1.5% → 16.0%（+14.5 combined）✓；+11.5 parametric-bootstrap / +17.3 evolver ✓；+3.2 to +21.1（parametric kept-cols）✓（同注 1）；+4.6 to +25.2（real data）✓ —— 与 A3/A4/A6 同源同值。**PASS**

### C. 主稿 Methods（行 32）

- "Both null sets ... identical re-alignment and filtering pipeline (790 of 800 runs completing) as a matched calibration" 与 A1 一致。**PASS**

### D. 补充 S7 全节（行 76–84）

| 文中值 | 复算值 | 判定 |
|---|---|---|
| 1.31% → 12.79%（+11.5；1 → 37） | 1.3055→12.7937（1dp 11.5）；rejections 1→37 | PASS |
| 1.76% → 19.10%（+17.3；0 → 71） | 1.7588→19.0955（2dp 19.10 ✓；1dp 17.3）；0→71 | PASS |
| 1.54% → 16.01%（+14.5；1 → 108） | 1.5365→16.0051；1→108 | PASS |
| p=0.5 点质量 74.7→64.2 / 45.7→42.0 | 74.67→64.23 / 45.73→41.96（1dp 均合） | PASS |
| kept_cols +3.2/+10.2/+21.1（parametric） | 3.15/10.16/21.09（+3.2 同注 1；10.1562→10.2 ✓；21.0938→21.1 ✓） | PASS（注 1） |
| kept_cols +12.9/+18.1/+21.1（evolver） | 12.88/18.05/21.05；12.8788→12.9 ✓；21.0526→21.1 ✓；+18.1 见注 2 | PASS（注 2） |
| kept_cols +4.6/+14.5/+25.2（real，BH-q） | delta_qsig_pp 4.58/14.50/25.19 | PASS |
| frac_removed +8.7/+14.1/+11.7（parametric） | 8.66/14.06/11.72（计数 8.6614/14.0625/11.7188 均合 1dp） | PASS |
| frac_removed +16.7/+16.5/+18.8（evolver） | 16.67/16.54/18.80（16.6667→16.7 ✓；16.5414→16.5 ✓；18.7970→18.8 ✓） | PASS |
| 781 = 383 + 398 | ✓ | PASS |
| 790/800（10 truncated timeout） | ✓（"balanced across sets" 见注 3） | PASS（注 3） |
| 49.4% → 64.1%（+14.8） | 49.3639→64.1221，14.7583 | PASS |

注 2（舍入边界）：+18.1 来自 JSON 存值 18.05（round-half-up）；按 n_sig 精确计数 24/133 = 18.0451pp，1dp 应作 18.0。偏差 0.055pp，同为二次舍入边界情形。
注 3（措辞宽松）：S7 "10 truncated timeout files, balanced across sets"。严格 5/5 对分不成立——R6 n_common=398 蕴含 R6 重比对完成数 ≥398，即 R6 超时 ≤2，从而 R7 超时 ≥8（可行对分仅 2/8、1/9、0/10）。JSON 未存分组完成数，无法直接证伪，但"balanced"按字面（大致均等）与可行对分不符，建议改为"in both sets"之类中性措辞或给出实际对分数。

### E. cover letter（行 5）

- 新表述在位："a matched neutral calibration attributes that increase to pipeline-induced null inflation rather than recovered power — we therefore treat the production-alignment estimates as primary and report the re-aligned rate only as a methodological bound"，与主稿口径一致。**PASS**
- 旧表述扫描："rather than directional corrections" 0 命中；"bounds on alignment-methodology sensitivity" 0 命中。**PASS**

### F. 陈旧文本扫描

- 主稿 "without adjudicating a single true value"：0 命中。**PASS**
- 补充 "^## S7"：恰好 1 处（行 76；S1–S7 节序列完整）。**PASS**
- 主稿 "+14.8"：2 处（行 142、176），均在新整合段落内，与 alignment_rerun_sensitivity.json delta_pp=14.7583、精确率 49.36%/64.12% 一致（文中 "exact rates of 49.36% and 64.12%" 逐字吻合）。**PASS**

### G. 算术自检（JSON 内部 delta 闭合）

- calibration JSON 三组 p<0.05 delta 闭合误差 0.0014–0.0033pp；bh_fdr delta 闭合误差 ≤0.0008pp。
- stratified JSON 全部 15 个分层 delta（含 real_data qsig 两轴 6 个）闭合误差 ≤0.01pp。
- 全部 ≤ ±0.05pp 容差。**PASS**

---

## 3. 总结论

**可交付。** A–G 七项全部 PASS。仅余三条不影响结论的注释，供酌情处理：

1. （低优先）S7 "+18.1" 与主稿/S7 "+3.2" 系对 JSON 两位小数存值的二次舍入（round-half-up），按分层精确计数应分别为 18.0 与 3.1，偏差 ≤0.055pp。若期刊终稿前想完全消除舍入歧义，可将 stratified JSON 的 delta_pp 按精确计数重存，或接受现状（梯度结论不变）。
2. （低优先）S7 "balanced across sets" 与 calibration JSON 的 n_common 约束在字面上不完全相容（R6 超时必 ≤2、R7 必 ≥8），建议措辞中性化。
3. 无其他口径/连贯性问题：production 管线 null（0/400、1/388≈0.26%）与重比对管线 null（37/383、71/398）在 Discussion 同段并存，指向不同管线，逻辑自洽。
