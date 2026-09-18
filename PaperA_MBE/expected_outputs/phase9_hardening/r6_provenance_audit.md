# R6 数据溯源审计报告（P0-5 核心支撑 + ILS/segdup/primate-P1 支撑材料）

**审计人**: fix-provenance（比较基因组学数据溯源）  
**日期**: 2026-09-17  
**对象**: manuscript_english_v9.md（Paper A）及其数据管线  
**性质**: 仅审计与考证，未改动稿件

---

## 摘要（一句话结论）

稿件 "HARs (3,171)" 的真实来源是 **Girskis et al. 2021 (Neuron) 合并的 3,171 条 HAR 目录**（整合 Pollard 2006 等 7 项发现研究的重叠合并结果），被 **Shin et al. 2024 (Cell Genomics) caMPRA 复用**，我们取其 mmc2 坐标表（hg19）经 pyliftover 转 hg38；"hCONDELs (583)" 来自 **McLean et al. 2011 补充表 2 的完整目录**（共 583 条，其中 510 条经序列验证——稿件数字本身正确，但须补调和句）。**审计中发现一个严重的方法学问题**：hCONDEL 坐标法基因映射存在 hg18/hg38 基因组版本错配，需重算（名称法复算后结论存活：RD OR=2.90, p=0.0023）。

---

## 任务 1【P0-5 核心】HAR/hCONDEL 集合溯源

### 1.1 HAR "3,171" 的完整证据链

**文献考证（已联网核实）**：

1. **Pollard et al. 2006**（*PLoS Genet* 2:e168）只报道 202 个 HAR——审稿人指出的问题属实。稿件 Introduction 将其作为 HAR 概念引文是恰当的，但它**不是** 3,171 的来源。
2. **3,171 的直接来源是 Girskis et al. 2021**（Girskis KM, Stergachis AB, DeGennaro EM, Doan RN, ..., Walsh CA. *Rewiring of human neurodevelopmental gene regulatory programs by human accelerated regions*. **Neuron 109(20):3239–3251.e7, 2021**, 10.1016/j.neuron.2021.08.005）。原文（Results）："We then applied the caMPRA approach to all **3,171** reported HARs (Bird et al., 2007; Bush and Lahn, 2008; Gittelman et al., 2015; Hubisz and Pollard, 2014; Lindblad-Toh et al., 2011; Pollard et al., 2006; Prabhakar et al., 2008)"——即 **7 项独立发现研究合并去重后的目录**。该文成功捕获 3,132 个并建立 **HARhub UCSC track hub（hg38）**，数据存 GEO **GSE180714**。
3. **Shin et al. 2024**（*Cell Genom* 4(8):100609, 10.1016/j.xgen.2024.100609）的 STAR Methods："The set of 3171 HARs examined in this study were selected from a number of different studies that identified HARs separately... Identified HARs that overlap were merged."——**同一目录的复用**。其 mmc2 即该目录的坐标表（hg19）+ mmc3 caMPRA 活性结果。
4. **稿件现状**：引文列表有 Shin 2024、Pollard 2006、Prabhakar 2006、McLean 2011，**缺 Girskis 2021**——这是引文缺口。

**本地文件内证（data/ 与 results/）**：

| 文件 | 内容 | 数量 | 坐标版本 | 用途 |
|---|---|---|---|---|
| `data/downloads/doan_2024/mmc2.xlsx`（HARs sheet） | Shin 2024 HAR 坐标表 | 3,171 | hg19 | 上游来源 |
| `results/doan2024_campra/hars_campra_merged.csv` | mmc2 + pyliftover + caMPRA 合并 | 3,171 行 | hg19→hg38（成功 3,169） | phase8 全 HAR 基因映射（稿件 "3,171" 出处） |
| `results/doan2024_campra/active_hars.csv` | caMPRA 活性 HAR | 508 行 | hg38 | 稿件 "508 elements" 出处 |
| `data/human/hgTables.txt` | UCSC Table Browser 导出（ID 形如 `HARsv2_0001`，推断为 HARhub hg38 导出） | **3,168** | hg38（原生） | phase7e/f HAR 比对与分类用 `n_hars`（"476 HAR-proximate genes" 出处） |
| `results/phase8_tissue_analysis/har_gene_mapping.csv` 等 | 输出 | — | — | 476（HARsv2 集）与 364（mmc2 集）两套并存 |

**两套 HAR 集合的一致性**：`doan_summary.json` 显示坐标重叠 3,161（mmc2 liftOver 后 vs 本地 HARsv2）——即 98%+ 为同一目录，但**稿件内并存两套计数**：3,171（mmc2，phase8）与 3,168（HARsv2，分类管线），且 "476 HAR-proximate" 来自后者、"3,161 lifted" 与 `doan_summary.json` 的 liftover_success=3,169 不完全一致。**"HARsv2 共 3,167" 的记忆线索不确——实测 3,168 行**。建议稿件统一表述并给出调和句。

**caMPRA 活性数字核实**：508（active_any）/226（active_both），稿件 "508 (16.1%)" 正确（508/3,171=16.02%）。注意 Shin 2024 正文报道 HAR 活性率 14.2%（N2A 细胞），与我们 16.1% 的口径（D1/D3 任一天、p.adj<0.05、log2>0，element 级聚合）不同——Methods 已写明标准即可，但建议加一句与原文 14.2% 的差异说明。

### 1.2 hCONDEL "583" 考证（问题解决）

- McLean et al. 2011（*Nature* 471:216–219）**共鉴定 583 个 hCONDEL**；正文中 "510" 指其中经 Trace Archive 计算验证的子集（510/583 = 87.5%；39 个另经 PCR 验证）。Nature 官网补充材料说明明确写明 "Supplementary Table 2 ... providing access to **all 583** hCONDELs"。
- 我们使用的 `data/downloads/hcondels/hCONDELs_supplementary_table2.xls`（874 KB，OLE2 .xls）实测：584 数据行，去 1 行子表头后 **583 条唯一记录**（hCONDEL.1–583；interg 410 / intron 165 / exon 8），坐标为 panTro2 与 hg18（Build 36.1）双套。**稿件 583 正确**，问题只是未写明"583=完整目录、510=验证子集"的调和句。

### 1.3 ⚠️ 严重发现：hCONDEL 坐标法映射的基因组版本错配

`phase8_hcondels.py`（第 123–126 行）解析 hg18 坐标后**未做 liftOver**，直接与 `results/phase4_conservation/gene_body.bed`（**GRCh38/hg38**，GENCODE v47）做 ±50 kb 区间重叠（Method 2，产出 121 个"coord"基因；与名称法 80 个合并为稿件中的 183 个 hCONDEL 基因）。

**量化验证**（本次审计新增，脚本见 C:\temp\prov\hc_check.py）：以 McLean 自己注释的 "Within" 基因为金标准，44 个 Within 基因落在 4,974 universe 内的 hCONDEL 中，hg18 坐标直接与 hg38 基因体重叠仅 **7/44 (16%)**；**32/44 (73%) 完全落错位置**（连 ±50 kb 都无法恢复；如 hCONDEL.1/TP73：hg18 chr1:3.617Mb，GRCh38 中 TP73 位于 chr1:3.66Mb，尚属好的情形；chr1 大偏移 >1 Mb 的座位比比皆是）。

**影响**：稿件中 183 个 hCONDEL 基因里约 103 个为"仅坐标法"命中，其中相当部分是**随机错配基因**。受连锁影响的稿件数字：183 基因、29 个 HAR∩hCONDEL 双证据基因、φ=0.042、摘要中 LOO hCONDEL OR 2.3–3.9 区间、Tier 4 hCONDEL 子集编码选择率（30.6–36.6%）、coverage 3.7%、Fig 4c/4f、Fig 7 相关 panel。

**名称法复算（build-safe，McLean 自带 within/upstream/downstream ≤50 kb 注释，80 基因）**：
- RD: 12/293, **OR = 2.90, p = 0.0023**（对比发表的 OR=2.94, p=6.7×10⁻⁶）
- GD: 11/1214, OR = 0.49, p = 0.994

**结论：定性结论（RD 富集 hCONDEL、作为主要独立验证）存活，但 p 值与全部连锁数字必须重算**。建议修复路径（二选一，推荐 A）：
- A. 对 583 条 hg18 坐标做 pyliftover hg18→hg38（管线已有 pyliftover 依赖）后重跑 Method 2，并保留名称法交叉验证；
- B. 退守名称法（80 基因）并在 S1 说明保守性。

### 1.4 基因组版本与坐标处理汇总

| 数据 | 原始版本 | 处理 | 目标版本 |
|---|---|---|---|
| HAR（mmc2） | hg19 | pyliftover（3,169/3,171 成功；3,161 与 HARsv2 坐标一致） | hg38 |
| HAR（HARsv2/hgTables.txt） | hg38（原生） | 无 | hg38 |
| hCONDEL | hg18/panTro2 | **未 liftOver（bug）** | 与 hg38 基因体直接比对 |
| 基因体 BED | GRCh38（GENCODE v47） | — | hg38 |
| 映射方法 | 自定义 Python 区间重叠，基因体 ±50 kb，非最近基因 | — | — |

### 1.5 可直接写入 Methods 的 provenance 段落（英文）

> **Regulatory element provenance.** Human accelerated regions were taken from the merged catalog of 3,171 HARs compiled by Girskis et al. (2021), which integrates and de-overlaps calls from seven discovery studies (Pollard et al. 2006; Bird et al. 2007; Prabhakar et al. 2008; Bush and Lahn 2008; Lindblad-Toh et al. 2011; Hubisz and Pollard 2014; Gittelman et al. 2015). Coordinates and caMPRA activity calls were obtained from Shin et al. (2024; hg19 coordinates in their Supplementary Table 2) and lifted to GRCh38/hg38 using pyliftover (3,169 of 3,171 loci mapped; loci with unparseable or unmapped coordinates were excluded, leaving 3,161 elements concordant with an independent hg38 export of the same catalog, N = 3,168). Elements were assigned to GENCODE v47 gene bodies within ±50 kb. Human-specific conserved deletions comprised the complete catalog of 583 hCONDELs from McLean et al. (2011; their Supplementary Table 2; 510 of the 583 were sequence-validated in that study), with gene assignment via the within/upstream/downstream annotations of McLean et al. (≤50 kb) and ±50 kb coordinate overlap after hg18→hg38 lift-over. Note that the frequently cited figure of 510 hCONDELs refers to the sequence-validated subset, whereas we analyzed the full catalog.

（若采用名称法方案，末句改为 "with gene assignment via the within/upstream/downstream annotations of McLean et al. (≤50 kb), which are internally consistent with the hg18 coordinates of that study; coordinate-based mapping after hg18→hg38 lift-over yielded concordant results."）

### 1.6 稿件需修改处清单（共 7 处）

1. **Methods 第 28 行 "Regulatory evidence:" 句**：替换为上文 1.5 provenance 段落（现句未交代 3,171/583 的确切出处）。
2. **引文列表补 Girskis et al. 2021**（现缺失；这是 3,171 目录的版权来源）。Shin 2024 引文信息已核实无误（Cell Genom 4:100609）。
3. **Introduction（第 18 行）**：Pollard 2006/Prabhakar 2006 作为概念引文可保留，但若读者会误以为它们是计数来源，可加 "(see also Girskis et al. 2021 for a merged catalog)"。
4. **hCONDEL 调和句**：加入 "583 (complete catalog; 510 sequence-validated)" 表述（可并入 1.5 段落）。
5. **⚠️ hCONDEL 全部连锁数字重算**（liftOver 修复后）：183/29/φ/OR 2.3–3.9（摘要）/Tier4 子集率/coverage/Fig 4c、4f/Fig 7。名称法下 RD OR=2.90, p=0.0023。
6. **HAR 计数一致性**：统一 "3,171（mmc2 口径）/3,169 成功 liftOver/3,161 进入映射" 与 "3,168（HARsv2 口径，476 基因）" 的表述；修正 "3,161 coordinates were lifted" 与 liftover_success=3,169 的矛盾（3,169 为 liftOver 成功数，3,161 为通过下游坐标校验数，须写清）。
7. **活性 HAR 基因数口径**：campra 摘要 n_genes_active_har=89（universe 内）vs 稿件 "97 genes within ±50 kb"（全集口径）——二选一并写明口径；另建议注明与 Shin 原文 14.2%（N2A）的口径差异。

---

## 任务 2【P0-6 支撑】ILS 定量材料

### 2.1 权威数字（均已联网核实）

| 事实 | 数值 | 来源 |
|---|---|---|
| 人–黑猩猩–大猩猩 ILS | **约 30% 的基因组**中，大猩猩比人–黑猩猩彼此更接近任一方；**编码基因周边 ILS 更少见**（重要缓解点！） | Scally et al. 2012, *Nature* 483:169–175（大猩猩基因组论文） |
| 人–黑猩猩–倭黑猩猩 | >3% 的人类基因组更接近倭黑猩猩或黑猩猩之一 | Prüfer et al. 2012, *Nature* 486:527–531 |
| 框架性综述 | 短内支 + 大 Ne → 基因树不一致；直接从基因组数据定谱系 | Scally & Durbin 2012, *Nat Rev Genet* 13:745–754, 10.1038/nrg3185 |
| 机制（SPILS） | 不一致基因树上发生的替换在固定物种树上须解释为多次替换 → 系统性拉长/缩短特定末端支；**固定物种树做正选择检验会抬高假阳性** | Mendes & Hahn 2016, *Syst Biol* 65:711–721, 10.1093/sysbio/syw018 |
| branch-site 对拓扑敏感 | 基因树选择直接影响正选择位点的推断 | Diekmann & Pereira-Leal 2016, *Evol Bioinform* 11(Suppl 2):11–17, 10.4137/EBO.S30902 |
| 灵长类普遍不一致 | 多个快速辐射 + 古老种间渗入；基因树不一致普遍 | Vanderpool et al. 2020, *PLoS Biol* 18:e3000954 |
| 现代方法学综述（可选） | 不一致存在下的比较方法 | Hibbins & Hahn 2023, *PNAS* 120:e2220389120 |
| 实证：物种树 vs 基因树选择检验结果差异（可选） | BUSTED/aBSREL/M1a-M2a 结果随树选择系统性改变 | Chang et al. 2024, *GBE*（murine rodents，PMC10491188） |

### 2.2 我们的补充树内支长度（`results/species_tree/species_tree_paml_foreground_fixed.nwk`）

- 人–(黑猩猩,倭黑猩猩) 内支：**0.003** 替换/位点
- ((人,黑猩猩科),大猩猩) 内支：**0.003**
- 人末端支：0.006；大猩猩末端支：0.012
- 与团队记忆的 "内支长 0.003–0.006" 一致（另一个 0.004 为 (人猿科,猴科) 内支，0.010 为更深层内支）。

### 2.3 Discussion 段落草稿（英文，~150 词，可直接粘贴后微调）

> A further caveat concerns incomplete lineage sorting (ILS). Speciation intervals in the hominine radiation were short relative to ancestral effective population sizes, so that roughly 30% of the genome carries a genealogy discordant with the species tree, with gorilla closer to human or chimpanzee than the latter are to each other (Scally et al. 2012; Scally and Durbin 2012). Our calibrated tree shows the same signature: the internal branches separating human from (chimpanzee, bonobo) and from gorilla are each only 0.003 substitutions per site. Under gene-tree discordance, substitutions that occurred on discordant genealogies must be resolved as multiple changes on a fixed species tree, systematically distorting terminal branch lengths (SPILS; Mendes and Hahn 2016) and inflating apparent lineage-specific dN/dS. Because BUSTED conditions on a single species-tree topology (Murrell et al. 2015), part of the human-foreground signal may reflect this artifact, although ILS is markedly rarer around coding genes (Scally et al. 2012), and pervasive primate discordance (Vanderpool et al. 2020) argues for gene-tree-aware reanalysis as future work.

### 2.4 引文清单（稿件格式）

- Scally A, Dutheil JY, Hillier LW, et al. 2012. Insights into hominid evolution from the gorilla genome sequence. *Nature*. 483:169–175. 10.1038/nature10842.
- Scally A, Durbin R. 2012. Direct determination of hominid phylogeny from genomic data. *Nat Rev Genet*. 13:745–754. 10.1038/nrg3185.
- Mendes FK, Hahn MW. 2016. Gene tree discordance causes apparent substitution rate variation. *Syst Biol*. 65:711–721. 10.1093/sysbio/syw018.
- Diekmann Y, Pereira-Leal JB. 2016. Gene tree affects inference of sites under selection by the branch-site test of positive selection. *Evol Bioinform*. 11(Suppl 2):11–17. 10.4137/EBO.S30902.
- Vanderpool D, Minh BQ, Lanfear R, et al. 2020. Primate phylogenomics uncovers multiple rapid radiations and ancient interspecific introgression. *PLoS Biol*. 18:e3000954. 10.1371/journal.pbio.3000954.
- （可选）Prüfer K, Munch K, Hellmann I, et al. 2012. The bonobo genome compared with the chimpanzee and human genomes. *Nature*. 486:527–531. 10.1038/nature11128.
- （可选）Hibbins MS, Hahn MW. 2023. Phylogenomic comparative methods: accurate evolutionary inferences in the presence of gene tree discordance. *PNAS*. 120:e2220389120.

---

## 任务 3【primate P1-6】区段重复（segmental duplication）过滤检查

**检查结果：未做任何 segdup 过滤。** 全部 `scripts/` 下 grep `segdup|segmental|duplication filter` 零命中；phase2 严格一对一 orthology 只做了 "orthologues, not paralogs" 的 Ensembl Compara 过滤与 reciprocal BLAST；重复屏蔽（RepeatMasker/soft-mask）也未应用于 CDS 提取（CDS 本身基本无重复，影响主要在旁系同源与多比对区）。本地无 UCSC segdup track 数据（downloads 目录无 genomicSuperDups），故本次无法直接量化 universe 基因与 segdup 的重叠率。

**建议写入 Methods/Limitations 的句子（英文）**：

> No explicit segmental-duplication filtering was applied to the gene universe. Orthology relied on Ensembl Compara one-to-one calls plus reciprocal BLAST, which excludes most paralogous loci but does not remove genes embedded in recent segmental duplications, where mis-assembly or paralogous mapping can inflate apparent lineage-specific substitution rates. A post hoc intersection of the 4,974 gene bodies with the UCSC hg38 segmental-duplication track (genomicSuperDups) is provided in Supplementary Text SX [若执行] / remains a limitation of the present analysis [若不执行].

**快速可执行方案**（供团队决定）：下载 UCSC hg38 `genomicSuperDups.txt.gz`（~数 MB），与 `gene_body.bed` 做区间交集，报告 4,974 universe 中 segdup 重叠基因数与分类分布（GD/RD/neutral），30 分钟内可完成。

---

## 任务 4【primate P1-4】2023 灵长类系统发育基因组学对话素材

**准确引文（已核实）**：

- Kuderna LFK, Gao H, Janiak MC, et al. 2023. **A global catalog of whole-genome diversity from 233 primate species**. *Science*. 380(6648):906–913. 10.1126/science.abn7829.
- Shao Y, Zhou L, Li F, et al. 2023. **Phylogenomic analyses provide insights into primate evolution**. *Science*. 380:eabn6919. 10.1126/science.abn6919.（稿件已引，为 82 PSGs 来源）
- Zoonomia：Christmas MJ, Kaplow IM, Genereux DP, et al. (Zoonomia Consortium). 2023. **Evolutionary constraint and innovation across hundreds of placental mammals**. *Science*. 380(6643):eabn3943. 10.1126/science.abn3943.

**各自能/不能提供什么（弹药）**：

| 资源 | 能提供 | 不能提供 |
|---|---|---|
| Kuderna 2023（233 种） | 覆盖 86% 属、全部 16 科的高覆盖 WGS；种内多态性（π）与气候/社会性关联；物种级核 DNA 系统树与分化时间；突变率的种间差异；"曾被认为是人类特异的错义突变在其它灵长类中反复出现"（直接削弱"人类特异替换"叙事，对我们 10 物种前景检验是重要外部检验） | 无 MPRA/caMPRA 式功能验证；单基因 dN/dS 前景检验（其框架是谱系基因组学而非逐基因选择检验）；无细胞类型分辨率表达 |
| Shao 2023（50 种，82 PSGs） | 50 物种谱系的前景选择检验（更深分类单元、更多内支平均化 ILS）；MPRA 功能数据（我们已用于交叉验证） | 人末端支特异性（其检验针对各支而非仅人类）；universe 与我们 4,974 只有 18 个可分析交集（稿件已如实报告） |
| Zoonomia 2023（240 哺乳动物） | 跨哺乳类的约束注释（3.32 亿受约束碱基，80% 在外显子外）；人类基因组约束基线；可作非编码元件功能的进化注释层 | 灵长类内部分辨率有限（深哺乳类比对稀释近期灵长类信号）；无逐基因人类支选择检验；无实验功能验证 |

**可用的对话角度**：(i) 我们的 10 物种人末端支检验与 Shao/Kuderna 的多支检验互补——前者问"人支上发生了什么"，后者受 ILS 影响更小但特异性更弱；(ii) Kuderna 的"人类特异错义突变反复出现"直接支持我们对 34% BUSTED 率作为上界的保守解读；(iii) Zoonomia 约束层可与我们的 phyloP/phastCons 组件互为印证。

---

## 附录：审计方法与可复现性

- 本地核查：`doan2024_parse.py`、`phase8_regulatory_tissue.py`、`phase8_hcondels.py`、`extract_conservation_scores.py`、`phase5_human_specific_screening.py`、`doan_summary.json`、`campra_regulatory_summary.json`、`hcondels_summary.json`、`hcondel_gene_mapping.csv`、`species_tree_paml_foreground_fixed.nwk`、`data/human/hgTables.txt`、hCONDEL xls 原始解析。
- 新增量化脚本（本次审计）：`C:\temp\prov\hc_check.py`（hg18/hg38 错配率）、`C:\temp\prov\hc_nameonly.py`（名称法 RD/GD 富集复算）。
- 项目内临时文件已清理；C:\temp\prov 保留脚本与中间输出供复核。

## 后记（2026-09-17）：修复方案 A 已执行

§1.3 的修复已由 `scripts/phase9_hcondels_liftover_fix.py` 完成（不覆盖旧脚本）：583 条 hg18 坐标 pyliftover 转 hg38，成功 581/583（失败：hCONDEL.449、hCONDEL.577；panTro2→hg38 链不可用）。金标准验证：within 基因直接恢复 40/44（修复前 7/44），错位 2/44（修复前 32/44）。旧坐标法 103 个基因仅 7 个被修复后方法保留——证实旧坐标集合基本为伪影。全部新旧数字对照见 `results/phase9_hardening/hcondels_liftover_fix_comparison.csv` 与 `hcondels_liftover_fix.json`；`tier4_empirical_disjointness.csv` 已更新（备份 `.bak_pre_hcondel_fix`）；替代映射表 `hcondel_gene_mapping_fixed.csv`。核心变化：hCONDEL 基因 183→112（名称 80 / 修复坐标 102 / 交集 70）；RD 富集 OR 2.94→**3.19**（p 6.7e-6→7.8e-5，仍显著）；HAR∩hCONDEL 29→22（φ 0.042→0.052，χ² p 0.0049→4.6e-4）；LOO hCONDEL OR 区间 2.27–3.91→**2.13–4.67**（p 区间 2.5e-6–6.4e-4→2.3e-5–5.3e-4，五个变体全显著）；Tier4 hCONDEL 子集 BUSTED 率 30.6%→30.4%（vs 补集 34.1%，差异仍不显著，结论不变）。所有定性结论存活。
