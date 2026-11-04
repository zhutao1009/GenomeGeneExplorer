# 核心算法与坐标约定

本文档说明 Genome Gene Explorer 里**必须小心处理**的部分：坐标系、链方向、区段划分、参考匹配、比对判定、压缩包角色识别、数据集库、性能设计、彩色文档生成。
改动 `gge/extract.py`、`gge/validate.py`、`gge/annotation.py`、`gge/bundle.py`、`gge/library.py`、`gge/seqdoc.py` 之前请先读完对应章节。

| 章节 | 内容 | 对应模块 |
| --- | --- | --- |
| 1–2 | 坐标转换、链方向 | `coords.py` |
| 3 | 注释解析与自检 | `annotation.py` |
| 4 | 八类序列提取 | `extract.py` |
| 5–6 | 参考匹配、比对与判定 | `validate.py` |
| **7** | **压缩包内文件角色识别（打分表、mRNA/CDS 长度吻合、BGZF、数据集清单）** | `bundle.py` |
| **8** | **数据集库表结构、判重、删除归属** | `library.py` |
| 9 | 性能设计（含彩色文档的体积模型） | — |
| **10** | **彩色序列文档：最小 OOXML 包、行模型、Word 验收方法** | `seqdoc.py` |
| 11 | 容易踩的坑 | — |

---

## 1. 坐标系：只允许一次转换

| 表示 | 用在哪里 | 约定 |
| --- | --- | --- |
| **1-based 闭区间** `[s, e]` | 注释原始值、`Segment`、界面显示、导出表格、FASTA 标题 | 与 GFF/GTF 文件完全一致，不做任何“顺手 +1/−1” |
| **0-based 半开区间** `[start0, end0)` | Python 切片、`Block.regions` 的块内偏移、`local_offsets()` 返回值 | 只用于取子串 `seq[start0:end0]` |
| **基因组坐标** | 只出现在 `Segment`、`Block.genomic_interval`、`Region.genomic_start/end` | 永远是 1-based 闭区间 |

唯一允许做换算的地方是 `gge/coords.py`：

```python
to_0based_closed(s, e) -> (s - 1, e)      # 1-based 闭 → 0-based 半开（长度不变）
to_1based_closed(s0, e0) -> (s0 + 1, e0)  # 反向
```

**铁律**：一段序列的长度恒等于 `e - s + 1`（1-based）或 `end0 - start0`（0-based）。任何一处长度对不上，就是这里被改了两次。分屏序列里的 `Region.start0/end0` 是**块内偏移，不是基因组坐标**——把它们当基因组坐标用，是这类软件最常见的错误来源（剪接后的坐标本来就与基因组坐标无关）。

`Block` 用 `coord_mode` 区分三种坐标语义：

| `coord_mode` | 含义 | 面板 |
| --- | --- | --- |
| `genomic` | 序列是基因组连续片段，`coord_start` 是该块第一个碱基的基因组坐标 | 基因组连续序列、上游、下游 |
| `spliced` | 序列是拼接产物，坐标是**转录本内**位置（1 起） | 拼接 mRNA、CDS |
| `per_segment` | 每段各自从 1 开始 | 外显子、内含子、UTR |

`coord_reverse` 只在负链的基因组块上为 `True`（序列是反向互补后的，基因组坐标沿序列方向递减），用于把块内偏移映射回基因组坐标：

```python
Block.region_at(offset) -> 基因组坐标
```

---

## 2. 链方向：先倒序，再反向互补

这是第二容易出错的地方。负链转录本按转录顺序读，等价于：

1. 把所有外显子/外显子片段**按基因组坐标降序**排列（`order_for_transcript()`）；
2. 顺序拼接每个片段的**反向互补**（`revcomp()`，IUPAC 感知：`R↔Y`、`K↔M`、`B↔V`、`D↔H`、`N↔N`、大小写保留）。

**只做 1 或只做 2 都会得到错误序列**，而且结果长度完全正确、GC 含量也近似，很难从表象发现。`tests/test_extract.py` 用人工可读的小序列专门锁住这两步（负链基因 `g2` 的 mRNA 与手工拼接结果逐碱基相等）。

上下游同样按**转录方向**定义，不是按坐标方向：

| 链 | 上游区间 | 下游区间 |
| --- | --- | --- |
| `+` | `[tx_start − up, tx_start − 1]` | `[tx_end + 1, tx_end + dn]` |
| `−` | `[tx_end + 1, tx_end + up]` | `[tx_start − dn, tx_start − 1]` |

再各自按转录方向取反向互补。区间先 `clip_interval()` 夹到 `[1, 序列长度]`，夹掉的部分在面板里写明“实际只能取到 N bp（已到序列边界）”，不会悄悄给你一段长度不对的序列。

---

## 3. 注释解析（`gge/annotation.py`）

六个阶段，每阶段结束都做自检：

1. **`detect_format()`**：看有没有 `##gff-version`、属性里是 `ID=` 还是 `gene_id "..."`，决定方言；
2. **`classify()`**：按第 3 列 + 属性特征判定 gene / mRNA / exon / CDS / 其它；
3. **建转录本**：`mRNA|transcript|mRNA_parent` 等类型建 `Transcript`，`Parent` 可指向多个（一个 exon 共享给多个转录本的情况是合法的，**不复制、不合并**）；
4. **挂子特征**：exon/CDS 按 `Parent` 挂到对应转录本，**按转录本各自的链方向排序**（同一份 exon 挂到正链和负链两个转录本，顺序不同）；
5. **补全**：没有 exon 但有 CDS 的转录本，用 CDS 反推外显子（`exons_derived=True`，界面标注“外显子由 CDS 推导”）；没有 gene 行的转录本，合成一个 gene；
6. **一致性与索引**：统计、建 `tx_index`/`name_index`、收集警告。

容错与警告（**报告而不是修正**）：

| 情况 | 处理 |
| --- | --- |
| 坐标 `start > end` | 交换并记警告（示例数据 `g7`），保留原始行号便于回溯 |
| `Parent` 指向不存在的转录本 | 记“N 个子特征的 Parent 在注释中不存在”，**不**偷偷挂到 gene 上 |
| 外显子相互重叠 | 保留原样并警告；内含子计算时重叠处取 `None`（不产生负数长度区间） |
| 列数 < 8 / 无法解析 | 跳过该行并计入 `stats["skipped"]`，在文件行提示 |
| 同一转录本 exon 重复 | 去重后警告 |
| seqid 在 FASTA 中不存在 | 提取时抛 `SeqIDMismatchError`，列出两边名字并给“是不是 chr1 / 1 / NC_xxx 的差异”提示 |

`AnnotationModel.stats` 里 `genes / transcripts / exons / cds / seqids / skipped / warnings` 是界面文件行显示的那串数字，也是“解析结果”可核对之处。

---

## 4. 八类序列的提取（`gge/extract.py`）

`extract_gene(provider, gene, tx, upstream, downstream) -> GeneExtraction` 一次性构造全部面板，**互不共享可变状态**（多转录本切换时不会串味，这一点有专门测试）。

| 面板 | 构造方式 | 着色区域 |
| --- | --- | --- |
| 基因组连续序列 | `[min(tx_start, up_s, dn_s), max(tx_end, up_e, dn_e)]` 单块，负链整体 revcomp | 上游/外显子/内含子/CDS/UTR/下游，按**真实基因组坐标**投影 |
| 上游 / 下游 | 各自的连续区间，单块 | 上游或下游单色 |
| 外显子序列 | 每 exon 一块，按转录顺序排列，块内编号 `外显子 i/n` | 外显子色 + 覆盖其上的 CDS/UTR |
| 内含子序列 | 相邻 exon 之间的空隙，负链同样降序编号 | 内含子色（负链片段已 revcomp） |
| 拼接 mRNA | 所有 exon 顺序拼接（负链倒序 + revcomp） | 拼接线处 CDS/UTR 颜色延续 |
| CDS | `tx.ordered_cds()` 顺序拼接 | CDS 单色，附 phase 信息 |
| UTR | 5′UTR = CDS 起点前的 exon 片段，3′UTR = CDS 终点后的片段 | 紫 / 青 |

**UTR 是推导值，不是注释内容**：NCBI/RefSeq 的 GFF 只有 gene/mRNA/exon/CDS，没有 UTR 行（草鱼注释里 UTR feature 数为 0），Ensembl 有时会额外给 UTR feature，但本程序统一按「外显子 − CDS」计算，不读注释里的 UTR。判据以**转录本坐标**为准：`first_cds_pos` 之前的 exon 碱基是 5′UTR，`last_cds_pos` 之后的为 3′UTR；夹在首尾 CDS 之间、但未被 CDS 覆盖的空隙**不标 UTR**（那说明注释里 CDS 不连续，属异常，校验会提示）。

**界面上必须如实标注来源**，否则用户会以为注释里写了 UTR：面板标题固定为 `UTR 序列（推导）`，标题栏徽标显示 `注释无 UTR` 或 `注释有 UTR 行`（`Panel.badge`，`GeneExtraction.utr_in_annotation` 决定），每段序列标题为 `5′ UTR（推导）`，块说明里写出依据（`首个 CDS 起点 x − 1 为止` / `末个 CDS 终点 y + 1 起算`）。`export.py` 的摘要与 `seqdoc.py` 的长度行也带同样的标记，导出的 FASTA 不会被误当成注释内容。

**区域着色优先级**（`REGION_PRIORITY`，数值大者覆盖小者）：

```
CDS(6) > UTR(5) > 外显子(4) > 内含子(3) > 上游(2) > 下游/侧翼(1)
```

实现上是逐碱基写一个 `bytearray` 类别码，再 `kinds_to_runs()` 压缩成连续区间给绘制用。逐碱基数组的代价是每个碱基 1 字节（3800 bp 的窗口 = 3.8 KB），换来的是“任意重叠都自然取到正确颜色”，不必写优先级判断的补丁逻辑。

**基因组块里的区域是投影出来的**：先算每个 exon 在块内的偏移，再把 CDS/UTR 的基因组区间与 exon 区间求交（`intersect()`），交出来的部分才有颜色。所以不会出现“整条 exon 被当成 CDS”这种错误显示——只有注释里真的标了 CDS 的部分才是橙色。

**CDS phase 只报告，不裁剪**。GFF 的 `phase` 描述的是“该 CDS 片段首个完整密码子之前还有几个碱基属于上一个密码子”，它**不意味**要丢掉那些碱基。按 phase 删碱基会直接让 CDS 长度算错。所以：

- 面板标题显示 `phase=0/1/2`；
- 校验里报告 `CDS 长度 % 3`、起始/终止密码子；
- 只有当注释本身不自洽（长度非 3 的倍数）时给出警告，**绝不修改序列**。

单外显子转录本没有内含子：面板返回 `Panel(is_empty=True, empty_text="该转录本仅 1 个外显子，无内含子")`，不抛异常（这是“正常结果”而不是“错误”）。

---

## 5. 参考序列匹配（`gge/validate.py`）

**绝不按文件顺序取第一条**，也绝不“第一条能对上就用”。流程：

1. `load_reference()` 读 FASTA，每条记录的 `token` = 标题第一个空白分隔字段（`>g1.1 cdna ...` → `g1.1`），并保留完整标题；
2. `score_candidates()` 对**所有**记录打分：

| 条件 | 加分 |
| --- | --- |
| transcript ID 与 token 或标题任一字段完全一致 | +100 |
| gene ID 完全一致 | +85 |
| 标题包含 transcript ID | +60 |
| 标题包含 gene ID | +50 |
| 标题包含 gene symbol（长度 ≥ 3） | +45 |
| 长度与期望长度完全相同 | +40 |
| 长度差异 ≤ max(3, 1%) | +12 |

3. 决策：

| 情况 | 状态 |
| --- | --- |
| 无候选（得分为 0） | `未匹配到参考`，附上参考文件前 3 条标题供核对 |
| 最高分有并列（分数差 < 1e-6） | `需人工确认`，列出全部候选 |
| 最高分 < 40 | `需人工确认`（得分过低不可靠） |
| 否则 | `matched`，使用最高分那条 |

4. 用户在下拉框里手动指定后，直接用指定记录，状态说明改为“按用户手动指定的对应关系进行校验”。

因此“两条参考都叫 g1.1”这种数据问题会停成 `需人工确认` 而不是随机挑一条——示例数据里专门放了这样一对等分候选。

### 5.1 蛋白参考的特殊处理（RefSeq / Ensembl 标题里没有转录本 ID）

RefSeq 的 `protein.faa` 标题长这样：

```text
>XP_051755185.1 zinc finger protein 750-like [Ctenopharyngodon idella]
```

只有**蛋白 ID**，既没有转录本 ID 也没有基因 ID。而上面那张打分表全靠转录本/基因 ID，
再加上早期版本给蛋白传 `expected_len=0`（注释写着"蛋白没有 transcript 长度可比较"），
等于把长度这条唯一的非 ID 证据也关了——结果就是：**蛋白参考文件永远"未匹配到参考"**，
用户只能从 5 万多条候选里手动挑一条。这是实测草鱼数据时才暴露出来的。

现在两条路一起走：

1. **从已匹配的参考 CDS 标题里取 `protein_id`**（确定性证据，不是猜）：

   ```text
   lcl|NC_067220.1_cds_XP_051755185.1_1 [gene=LOC127515509] [protein=…] [protein_id=XP_051755185.1] [gbkey=CDS]
   ```

   CDS 那一路已经逐碱基验证过（`identical`），它指向的蛋白就是唯一正确的对照。
   `protein_id_from_cds()` 用正则取 `protein_id=`，再确认该 ID 在蛋白文件里存在。
   **只在唯一命中时采用**：标题里出现两个都存在的蛋白 ID 就返回空字符串，
   退回打分 + 人工确认——宁可让你选，也不替你猜。
2. **取不到就靠长度打分**：`expected_len` 改为翻译结果去掉终止符的长度
   （参考蛋白 FASTA 一般不写终止符，所以 519 对 520 属于正常差异）。
   长度一致的记录会拿到 +40 分浮出来；如果有多条同长度并列，仍然停在
   `需人工确认` 并把候选列出来。

注意**不要**用"序列完全一致"来反推匹配关系：那会让校验变成自证——凡是能匹配上的
必然一致，差异永远查不出来。匹配只能依据 ID 与长度这类外部证据。

---

## 6. 比对与判定

### 6.1 比对

```python
PairwiseAligner: mode=global, match=+2, mismatch=-3,
                 open_gap=-6, extend_gap=-1, end_gap=-0.5
```

全局比对（Needleman–Wunsch）保证 indel 被**定位**出来，而不是按位置硬比。逐列扫比对结果：

| 列形态 | 记为什么 | 事件坐标 |
| --- | --- | --- |
| ref 有碱基 / obs 有碱基 / 相同 | match | — |
| ref 有 / obs 有 / 不同 | mismatch | `ref_pos`, `obs_pos`（1-based） |
| ref 是 gap | **insertion**（obs 相对 ref 多出碱基） | `obs_pos + 1` |
| obs 是 gap | **deletion**（ref 相对 obs 多出碱基） | `ref_pos + 1` |

注意方向容易搞反：**参考序列里少了一个碱基，在比对结果里表现为“观测序列多出一个碱基”**（insertion）。示例数据的 `demo_mrna.fa` 里那条 1 bp 缺失的参考，报告里就是 insertion，这是正确的，不是 bug。

`identity = matches / aligned_columns × 100`（分母是比对列数，含 gap 列）。`identical` 当且仅当 mismatch = insertion = deletion = 0。

### 6.2 超长序列

`max_bases = 60_000`。超过则走 `_coarse_compare()`：找公共前缀 + 公共后缀，中段只比较长度差（多出即 insertion/deletion，等长则数 mismatch），并在报告里标 `coarse=True` 说明“差异定位精度下降”。这样 200 kb 的转录本也不会把界面卡死。

### 6.3 状态

| 状态 | 触发条件 |
| --- | --- |
| `完全一致` | 比对无任何 mismatch/indel |
| `存在差异` | 有差异，且能给出差异位点 |
| `未匹配到参考` | 参考文件没有可用记录 / 该转录本无 CDS |
| `需人工确认` | 候选并列或得分过低 |
| `未运行` | 未加载该类型的参考文件 |

除此之外还有一组**与参考无关**的编码区检查（`observe_cds_checks`）：长度是否为 3 的倍数、起始密码子（含替代起始）、终止密码子、内部终止密码子数量、翻译产物前 40 aa。这些即使没有参考序列也会跑。

差异存在时，报告里会给出“排查方向”（注释与参考是否同一组装版本、外显子边界是否偏移、链方向、Transcript ID 是否指向不同转录本、可变剪接），**措辞明确是可能性而非结论**，程序不会断言原因。

---

## 7. 压缩包内文件角色的识别（`gge/bundle.py`）

### 7.1 出发点：文件名不可信

用户给的包里常见 `assembly_v1.fna`、`sequences_a.fna`、`sequences_b.fna`、`all.fa` 这种名字，
同一个物种不同来源的命名习惯也完全不同。**任何基于文件名的判断都只能是加分项，不能是决定项。**
所以识别的骨架是：先嗅探格式 → 再统计内容 → 最后**与注释做交叉验证**，交叉验证才是决定性的。

### 7.2 只读头部，不解包

扫描阶段绝不把文件解到磁盘（几十 GB 的包也不能先解再判）：

| 常量 | 值 | 用途 |
| --- | --- | --- |
| `_HEAD_BYTES` | 256 KB | 每个成员只读头部这么多字节，用于判断"是不是文本 / 是不是 FASTA / 列数够不够 8 列" |
| `_FASTA_SAMPLE_BYTES` | 32 MB | FASTA 统计最多读这么多字节 |
| `_FASTA_SAMPLE_RECORDS` | 5000 | 以及最多读这么多条记录 |

超过上限时 `FastaStats.complete = False`，界面与报告里的数字加「以上」「（采样）」字样——
**宁可说"至少这么多"，不假装是精确值**。注释的统计用 `_summarize_lines()` 流式扫（`tx_exon_len` / `tx_cds_len`
两个 dict 记每转录本的 exon 和 / CDS 和），几百万行也只占几十 MB，不构建完整模型。

两个采样上限是 `max_bytes: int = 0` 这种"函数内回填"的参数，**不是**默认参数值——
默认值在 `def` 时就绑定了，测试无法 monkeypatch（踩过）。

### 7.3 打分表（`_score_candidate`）

证据累加，最后取最高分角色；**没有任何一条证据是"一票决定"**：

| 证据 | 加分 | 条件 |
| --- | --- | --- |
| 蛋白字母构成 | 蛋白 **+70** | ACGTUN 占比 < 0.80；一旦判定为蛋白，核酸三个角色的分数直接清零（否则 `pep.faa` 会因为"长度能被 3 整除"抢 CDS） |
| 文件名关键词 | 每个命中角色 **+30**（`_NAME_WEIGHT`） | 只加分。ASCII 关键词按词边界匹配（`(^\|[^a-z])w([^a-z]\|$)`），中文关键词按子串匹配 |
| 序列名命中 `seqid` | 基因组 **+55** | FASTA 标题的 token 出现在注释的 seqid 集合里 |
| 序列名命中转录本 ID | mRNA **+26**、CDS **+26** | 两者都加，靠下面的长度吻合再分开 |
| 记录数 ≈ 转录本数 | mRNA **+20**、CDS **+20** | `0.7 ≤ n_records / n_transcripts ≤ 1.4` |
| 记录数远少于转录本数 | 基因组 **+22** | 比值 < 0.15 且记录数 ≤ 1000 |
| 长度吻合 | mRNA / CDS / 蛋白 **+45 × fit** | fit 见 7.4，要求 fit ≥ 0.5 才计分 |
| 长度能被 3 整除 | CDS **+12** | ≥ 90% 的记录长度是 3 的倍数，且记录数 ≥ 3 |
| 长度多数不能被 3 整除 | mRNA **+8** | ≤ 40%，同样要求记录数 ≥ 3（成熟 mRNA 含 UTR，通常不是 3 的倍数） |
| 无注释时的体积启发 | 基因组 **+18** | 文件 ≥ 20 MB 而记录数 ≤ 1000 |

判定门槛与置信度：

```text
最高分 < 25                    → 角色置空（"未识别"），证据写「证据不足，无法自动判定，请手动指定角色」
confidence = clamp(score/100, 0.2, 0.99)     # 只是给人看的相对把握，不是概率
ambiguous = (最高分 - 次高分) < 12 且 最高分 < 100
```

`ambiguous` 不会改变结论，只是往证据里追加一句「与「参考 CDS」得分接近，建议确认」。
**唯一的例外**是"整个包里只有一个核酸 FASTA、又没有任何注释可交叉验证"：
此时所有交叉证据都不存在，程序按基因组处理，但 `confidence` 固定 0.3，并强制附一条
「只能凭…当作基因组，请确认」的警告——**猜可以，但必须让你知道这是猜的**。

### 7.4 mRNA 与 CDS 靠长度吻合分开

两者都是核酸、都带转录本 ID，字母构成完全一样。真正的区别只有一个：**长度和谁对得上**。

```python
tol_abs = max(3.0, 0.02 * expect)          # 相对 2%、绝对至少 3 bp（容忍参考版本的小差异）
ok = abs(expect - length) <= tol_abs
```

- mRNA 记录的长度 ≈ 该转录本**所有外显子长度之和**（`tx_exon_len`）；
- CDS 记录的长度 ≈ 该转录本**所有 CDS 段长度之和**（`tx_cds_len`）；
- 蛋白记录的长度 ≈ `CDS/3`（去掉终止密码子后再容差 1 aa），所以 `_length_fit(..., allow_third=True)`。

标题 token 到转录本 ID 的映射交给 `_lookup_token()`：先试原样，再试 `|` / `:` 前的部分，
再剥掉 `_cds` `.cds` `-cds` `_mrna` `_transcript` `_protein` `_pep` 这类常见后缀，最后做一次大小写不敏感回退。

**至少 2 条记录比对得上才算证据**（`total < 2` 直接返回 0）——小基因组只给 3 条参考序列时，
1 条巧合的概率太高。同理"能被 3 整除"要求 ≥ 3 条记录，否则 `GATTACA` 这种单条记录就能骗到分数
（这两条都是实测踩出来的：早期版本的门槛设成 5 / 10，导致 3 条记录的 CDS 文件被判"未识别"）。

### 7.5 一个角色只能有一个文件

`_assign_unique()` 保留同角色里分数最高的那个，**把其余的降级为"未识别"并在证据里写明被谁压下去**，
再由 `_demoted_warnings()` 汇总进 `scan.warnings`。绝不静默挑一个：
"你给了两个基因组"是用户必须知道的信息。

### 7.6 压缩格式与 BGZF

- 容器：`zipfile.is_zipfile()` / `tarfile.is_tarfile()` 按**内容**判断，不看扩展名；`tar` 的成员大小取 `TarInfo.size`（不是 `file_size`）；
- 成员内部：`.gz` `.bz2` `.xz` 透明解开（嗅探与统计都在解压流上做）；
- **BGZF 例外**：gzip 头的 extra field 里带 `BC` 子字段，说明这是 bgzip 生成的可随机访问压缩文件。
  pyfaidx 能直接读 `.fa.gz`（BGZF），所以**不能**再解一层；而普通 gzip 的 `.fa.gz` pyfaidx 无法建索引，
  必须解成 `.fa`。判错的后果是要么白占一倍磁盘，要么 `open_fasta` 直接打不开——两个方向都有测试锁住。

### 7.7 解包的安全边界

| 风险 | 处理 |
| --- | --- |
| zip-slip / 目录穿越 | `safe_filename()` 拍平成 basename，`../` `..\` 与路径分隔符一律丢弃 |
| Windows 非法字符 | `< > : " \| ? *` 与控制字符替换为 `_` |
| DOS 保留名 | `CON` `PRN` `AUX` `NUL` `COM1..9` `LPT1..9` → 加前缀 `_`（`CON.fa` → `_CON.fa`） |
| 重名 | `_unique_name()` 追加 `_2` `_3` |
| 磁盘不够 | `need = 未压缩总大小 × 0.35 + 128 MB`，不够就在解包前报错，不留半成品 |
| 中途取消 / 失败 | `tasks.import_bundle_from_archive()` 用 `try/except` 包住整个解包过程，异常时 `shutil.rmtree(dest)` |
| 重复导入 | 解包成功后按来源路径删掉旧记录与旧目录，再登记新记录（同一包不会堆副本） |

---

### 7.8 数据集清单：识别能力的另一半（导出侧）

"按内容猜角色"是给**别人**的包用的；自己导出的包不需要猜。`write_dataset_archive()`
把当前这套文件打成 zip，并在根目录放一份 `gge_dataset.json`：

```json
{"format": "gge-dataset", "version": 1, "name": "…", "note": "…",
 "created_at": "2026-10-29T09:12:44", "generator": "Genome Gene Explorer 1.1.0",
 "stats": {"n_genes": 31801, "n_transcripts": 70058, "n_seqids": 31, "total_bases": 893220000},
 "files": [{"role": "genome", "member": "data/xxx.fna", "basename": "xxx.fna", "size": 905000000}]}
```

导入时 `_read_manifest_from_archive()` 会先找它，`_apply_manifest()` 按清单直接落地角色：

| 现状 | 处理 |
| --- | --- |
| 清单存在且 `format` 正确 | **先把所有候选的 `role` 清零**，再按清单指定；`score=300`、`confidence=0.99`、`ambiguous=False`，`evidence` 写明"按清单指定" |
| 清单里某个 `member` 在包里找不到 | 进 `ignored` 并说明"可能被改名或删掉"，不要静默少一个角色 |
| 同一文件被清单写成两个角色 | 只保留第一个，另一个标"已用作 X，此处忽略" |
| 清单是坏的（截断 / 手改错） | `json` 解析失败 → 返回 `None` → 完全退回按内容打分，包照样能打开 |
| 清单有 `stats.n_genes` | 跳过注释全量统计（`_summarize_lines`），直接构造 `AnnotationSummary`。524 MB 注释的这一趟是 60–70 s，跳过它才是"打开快"的关键 |
| 清单 `version` 高于本程序 | 目前只比对 `format`；未知版本仍按字段名尽力读取（缺字段当 0 / 空），不抛错 |

"清零再指定"这一步容易漏：如果只给命中的候选赋值，那些只靠打分得到的角色会残留在
其它 `FileCandidate` 上，于是"清单里 3 个文件、界面却显示 5 个角色"。测试
`test_manifest_wins_over_content_guessing` 钉的就是这个。

---

## 8. 数据集库（`gge/library.py`）

### 8.1 表结构

只用标准库 `sqlite3`，不引入新依赖，也不需要管理员权限。

```sql
meta(key, value)                                    -- schema_version
datasets(id, name, source, fmt, n_genes, n_transcripts,
         n_seqids, total_bases, note, created_at, updated_at, last_used)
files(id, dataset_id → datasets(id) ON DELETE CASCADE,
      role, member, path, size, UNIQUE(dataset_id, role))
```

- `PRAGMA foreign_keys = ON` 每个连接都要设（SQLite 默认关），否则 `CASCADE` 不生效，删记录会留孤儿行；
- `UNIQUE(dataset_id, role)` 让"一个数据集里一个角色只能有一个文件"在数据库层面也成为约束；
- 库里只存**元数据与路径**，序列仍在解包目录里。这样打开软件不必读几十 GB 的基因组，数据也能放在移动硬盘上；
- 删掉数据库不会丢数据；删掉数据目录则会在列表里显示"文件缺失"（而不是报错崩溃）。

### 8.2 排序与"最近使用"

```sql
ORDER BY COALESCE(NULLIF(last_used, ''), created_at) DESC, id DESC
```

时间戳是 `datetime.now().isoformat(timespec="microseconds")`（**微秒**，不是秒）：
同一秒内连续导入两个数据集时，秒级时间戳会让顺序变得不确定；末尾再补 `id DESC` 兜底。

### 8.3 按来源判重

`get_by_source(source)` 用来源压缩包的**绝对路径**（`os.path.abspath`）匹配。导入前先查，命中就替换。
这样"我把包重新打了一遍、放在同一个位置"不会在库里留下两份看起来一样的数据。

注意它的边界：把压缩包**改名或换目录**后再导入，会被当成一套新数据（路径不同了）。
这是有意的取舍——不做内容哈希，是因为对大包算校验和要整读一遍，代价远高于收益；
真要合并，请在库对话框里「移除」旧的那套。

### 8.4 删除时只删自己的文件

```python
def _owned_dir(self, file_path: str) -> str:
    root = os.path.abspath(store_root())                 # .../library/datasets
    target = os.path.abspath(os.path.dirname(file_path))
    if os.path.commonpath([root, target]).lower() == root.lower() and target != root:
        return target                                    # 允许递归删除
    return ""                                            # 用户自己的目录：只删记录，不动文件
```

- `commonpath` 会正确处理 `..`，比字符串 `startswith` 安全（后者会把 `datasets_evil/` 当成子目录）；
- 不同盘符会抛 `ValueError`，必须捕获（Windows 上 C: 与 D: 就是这种情形）；
- `target != root` 防止把整个存储根目录删掉；
- 界面上的确认框据此给出两种不同措辞：「将删除库内 5 个文件」 vs 「只从列表移除，文件保留在原处」。

---

## 9. 性能设计

| 关注点 | 做法 |
| --- | --- |
| 大基因组 | `pyfaidx` 随机访问，只读需要的区间；`.fai` 建立失败时小文件退回内存模式（`PlainFastaProvider`，< 64 MB），大文件直接报错并说明原因 |
| 索引并发 | `PyfaidxProvider` 内部加 `RLock`；索引缺失自动重建一次 |
| 界面不冻结 | 所有 IO/解析/提取/校验都在 `QThread`（`Job`）里跑，信号回主线程；连续请求进队列，一次只跑一个 |
| 可取消 | 任务函数接收 `cancel_check` 回调，长循环里轮询；取消后抛 `CancelledError`，界面回到可用状态 |
| 绘制 | 画布按视口虚拟化：只绘制可见行，不构建每碱基的 QWidget/QGraphicsItem；颜色用 run-length 压缩后的区间批量填充 |
| 缓存 | 同一（基因 + 转录本 + 上下游长度 + 参考文件版本）的提取结果缓存，参数不变时切回不重新计算 |
| 内存 | 3800 bp 窗口的逐碱基类别数组约 3.8 KB；面板只保留 `Panel` 数据，不保留 QPixmap |
| 扫描压缩包 | 只读成员头部（256 KB）与采样（32 MB / 5000 条），几 GB 的包也不落盘、不整解；实在要读全量时才解包 |
| 解包 | 1 MiB 分块流式复制 + 进度 + 可取消，不在内存里攒整个文件 |
| 注释解析 | 逐行流式读入 → 分类 → 建层级。属性列先判 `'"' in text` 再决定要不要跑 GTF 正则（见 §10-13） |
| 注释内存 | **167 万个 exon/CDS 行不保留属性字典**（只留 `tx_key`/`gene_key`），`seqid`/`ftype`/`strand` 用 `sys.intern` 共用；片段只保留 `strand`；挂载完成后 `del child_rows` |
| 注释进度 | 按「已读字符数 / 文件字节数」换算到 15→70%，之后 72/85/92 三个里程碑；不再固定报 40% |

### 9.1 实测基准（NCBI 草鱼 HZGC01，Windows 11 / 16 GB）

| 阶段 | 数据量 | 耗时 | 峰值内存 |
| --- | --- | --- | --- |
| `load_genome`（已有 `.fai`） | 863 MiB / 31 条 | 0.2 s | 几十 MB |
| `load_annotation_file` | 524 MB / 178 万行 | **≈ 70 s**（优化前 214 s） | **1.26 GB**（优化前 3.2 GB） |
| `load_reference_file` mRNA | 259.8 Mb / 67,808 条 | 2.8 s | 随记录数线性 |
| `load_reference_file` CDS | 134.7 Mb / 59,323 条 | 1.8 s | 同上 |
| `load_reference_file` 蛋白 | 44.8 Mb / 59,198 条 | 1.1 s | 同上 |
| `extract_and_validate`（单基因 3 路参考） | — | 4.2 s | 忽略不计 |

**注释解析是唯一的重活**，且只做一次：之后换基因只重新提取那几 kb，不再碰文件。

优化分三步做的，每一步都用真实数据量过：

| 改动 | 耗时 | 峰值内存 |
| --- | --- | --- |
| 原始实现 | 214 s | 3.2 GB |
| 属性列先判引号再跑 GTF 正则 + 值先判 `%` 再 unquote | 71 s | 3.2 GB |
| 片段不再复制整字典（`_segment_attrs`） | 60 s | 2.8 GB |
| **子特征行不再保留属性字典 + 字符串 interning** | ≈ 70 s | **1.26 GB** |

最后一步把 CPU 换回了一点（≈ +10 s），但这是必要的取舍：15 GB 内存、剩 3 GB 的机器上，
2.8 GB 的峰值会直接把 GUI 拖进换页——实测窗口里注释加载 4 分钟都没结束，看起来像死机。
1.26 GB 之后同样的操作 70 s 完成。

再往下要改成分两遍流式扫描（第二遍只读子特征、不建行对象），预计还能省几百 MB，
但那属于架构改动，收益与回归风险不成比例。

---

### 9.2 彩色文档的体积模型

`document.xml` 的大小基本由 **run 数** 决定（每个 run 约 200 B：`rFonts` + `sz` + `color` + `shd`）。
run 数又由"着色方式"决定，这是选默认值时最该算清楚的一件事：

| 着色方式 | run 数 ≈ | 实测每碱基开销 | 50 kb 面板的 docx |
| --- | --- | --- | --- |
| 按区域（默认） | 行数 ×（结构变化次数 + 分组空格数） | ≈ 4 B/bp | 几百 KB |
| 区域背景 + 碱基字母色 | ≈ 碱基数 × 0.9 | ≈ 210 B/bp | 8–10 MB |
| 只按碱基字母色 | ≈ 碱基数 × 0.9 | ≈ 200 B/bp | 8–10 MB |

所以默认是"按区域"：它既与界面所见一致，又把"颜色变化次数"从 O(碱基) 降到 O(结构)。
想加碱基配色也行，但要先把体积算给用户看（`estimate()` 用的就是上表的系数，再按 zip
压缩 4 倍折算），并在 20 万 bp 以上给出提醒。

合并规则（`BlockPainter.runs`）：

1. 先用 region 把行切成 `(s, e, kind)`，缺口补 `other`（区间查找用游标，只前进不回退）；
2. 段内再按 10 bp 分组切开，组间插一个**不带底纹**的空格（否则色块连成一条，看不清分组）；
3. 同色相邻的碱基合并成一个 run（`mode=region` 时整段就 1 个 run）；
4. 每行坐标来自 `Block.first_coord(line, bpl)` —— 与界面左栏共用同一份逻辑，不另算一遍。

### 9.3 怎么确认 Word 真的认了底纹

不能靠"看起来像"，可以让 **Word 自己说**：用 COM 打开导出的 docx，能算出页数就说明没有触发
"需要修复"，再另存为 HTML（`wdFormatHTML = 8`），检查那份 HTML 里是否出现
`background:#bbf7d0` 与 `font-family:"Consolas"`——出现了就说明 `<w:shd>` 被解析成了字符底纹、
`<w:rFonts>` 被解析成了字体。本机实测（Office 16）：无修复提示、页数正常、底纹与字体都在。

> 沙箱里跑 COM 会报 `0x80080005 服务器运行失败`，必须脱离沙箱执行；脚本见 §10 末尾。

### 9.4 启动自动恢复（`auto_restore`）

**问题**：老版本启动时只把上次的 5 个路径填回左侧列表（状态写"上次使用（未加载）"），
数据并没有加载。用户看到文件名都在，以为已经就绪，一搜索只得到"还没有加载注释文件"——
这是个**假就绪**状态，比空列表更误导。

**做法**：`__init__` 末尾 `QTimer.singleShot(120, self._auto_restore_files)`，按键序
`genome → annotation → mrna → cds → protein` 依次 `load_role()`，全部走原有的排队机制。

几个必须有的守卫：

| 守卫 | 原因 |
| --- | --- |
| `auto_restore` 关 或 `existing_path("genome_path")` 为空 | 没有可恢复的，或用户明确要手动 |
| `provider/model` 已存在 或 `_job/_queue` 非空 | 命令行 `--genome`、`--import-bundle`、示例数据已经在加载，不能抢队列 |
| `120 ms` 延时而不是 `0` | 窗口还没画完就开始加载，进度条与状态栏的更新会"压"在首帧里，用户看不到任何反馈 |
| 恢复期间 `_restoring=True` 时不弹"注释有提醒" | 开机就被对话框挡住是最招人烦的；提醒改写到状态栏与「数据校验」面板 |

**代价与取舍**：内存模型决定了注释每次启动都要重解析一遍（524 MB 的 GFF 约 60–75 s）。
这段时间是后台任务，界面可正常操作、随时可取消，但确实"打开就要等"。所以开关做成可关
（`文件 → 启动时自动加载上次的文件`，值存 `auto_restore`），给超大数据的用户留出口。
被移动/删除的文件由 `Settings.existing_path()` 过滤掉，不会让启动卡在报错对话框上。

> 顺手补的一件事：旧版的 `load_annotation_file` / `load_reference_file` **没有** `cancel_check`
> 参数，`Job` 按签名判断是否注入取消钩子，于是"取消"对这两个任务完全无效——500 MB 的 GFF
> 只能等它读完（实测点了取消 15 s 还在跑）。现在 `read_rows` 每 5 万行、外层挂载循环每 20 万行
> 各查一次取消标记，`read_fasta_records` 每 2 万行同样处理；实测 500 MB GFF 按下取消后
> **1.0 s** 结束（`tests/test_annotation.py` 里用 12 万行的合成 GFF 锁住了这个行为）。

**自检必须隔离**：`--self-test` 以前直接 `MainWindow()`，用的是用户真实的 `settings.json`
与 `GGE_LIBRARY_DIR`，结果每跑一次自检就往用户库里塞一套「演示数据集」，还把用户上次的真实
路径覆盖成示例数据（用户再打开就看到一堆演示文件）。现在自检改用 `<outdir>/selftest_settings.json`，
并把 `GGE_LIBRARY_DIR` 指到 `<outdir>/_library`。判据：跑完自检后真实库的 `datasets/` 目录数量不变。

**顺带恢复上次的基因**：`_save_settings` 记 `last_gene`；注释加载完、genome 也在时
`_auto_extract_if_ready()` 会检查 `last_gene in model.genes`（**必须查存在性**，否则同名不同版
的注释会直接命中 `select_gene()` 里的 `_error` 模态框，把启动吊死），命中就 `select_gene()`。

---

## 10. 彩色序列文档的生成（`gge/seqdoc.py`）

### 10.1 为什么自己拼 OOXML

需求是"把界面上那张彩色序列表贴进 Word"。纯文本天生带不了底色；`python-docx` 只是把 OOXML
包了一层，为几个色块再拖一个依赖进 PyInstaller 不划算。Word 真正需要的最小包只有三份部件：

| 部件 | 作用 |
| --- | --- |
| `[Content_Types].xml` | 声明 `word/document.xml` 是主文档 |
| `_rels/.rels` | 包级关系，指向正文 |
| `word/document.xml` | 正文（段落 + run） |

文字级底色用 `<w:shd w:val="clear" w:color="auto" w:fill="RRGGBB"/>`，写在 `w:rPr` 里——
就是 Word 工具栏"字符底纹"存盘时写的东西（Word 2007+ 支持）。页面尺寸/页边距放
`w:sectPr`；没有它 Word 会按默认 Letter 排版，看起来"怎么跑到右边去了"。

### 10.2 一份行模型，三个渲染器

`document_lines()` 把一次提取摊平成 `DocLine` 列表：`title` / `info` / `legend` / `heading` /
`subhead` / `note` / `line` / `spacer`。序列行只存"坐标字符串 + 若干 `(文本, 底色, 字色)`"。
docx 与 html 两个渲染器都只消费这个列表，所以**颜色、坐标、10 bp 分组位置必然一致**；
剪贴板富文本（`clipboard_html`）用的是同一个渲染器的 `full_page=False` 形态。

### 10.3 每个 run 自带字体

一度想让 `w:pPr/w:rPr` 承载默认字体，好把 `document.xml` 缩小一半（一条 60 bp 的行有约
70 个 run）。规范上 `w:pPr/w:rPr` 是"段落标记"的字符格式，Word 通常会继承，但不同阅读器
（WPS、LibreOffice、各种预览器）继承行为不一致；而序列行一旦丢字体就整体错版，省下的字节
不值得冒这个风险。**现在每个 run 都写 `w:rFonts` + `w:sz`**，只省掉 `w:szCs`/`w:cs` 这类
对 CJK/复杂文种才起作用的项。

### 10.4 HTML 侧的两个细节

- 用 `background-color` 而不是简写 `background`：Word 粘贴 HTML 时只认前者；
- 剪贴板里不放 `<!DOCTYPE>`：那会让部分宿主（WPS、某些富文本框）把它当成"要打开一个文件"。

### 10.4.1 坐标列是开关，不是排版常量

`document_lines(..., include_coord=False)` 时序列行的 `coord` 留空，渲染器据此**整段跳过**坐标
run 和它后面那两个占位空格：

```python
runs = ""
if line.coord:                       # 空串就不写 run，也不写那两个空格
    runs += _srun(line.coord, color="#9CA3AF", size_pt=font_pt)
    runs += _srun("  ", size_pt=font_pt)
```

为什么不做成"空字符串也照样写 run"：docx 里 `<w:t></w:t>` 仍是合法 run，一个 200 行的文档会白多
200 个空 run（约 30 KB），而且 HTML 侧会留下 `<span class="c"></span>` 这种谁都看不懂的碎片。
**默认关闭**，界面只在导出对话框里给一个「每行左侧显示坐标列」的勾选（选择记在
`settings.json` 的 `doc_include_coord`）。判据也简单：坐标用浅灰 `#9CA3AF`，所以测试里
"`w:val="9CA3AF"` 在不在 `document.xml` 里"就是"有没有导出坐标列"的可靠信号（正文信息里的
`chr1:1,001-1,800` 是另一回事，不能拿数字去 grep）。

关掉坐标后**段落数与带坐标版本完全一致**（只是每段少了前缀），这一点在
`tests/test_seqdoc.py::test_docx_without_coord` 里逐段比对过——防止哪天"顺手"把序列行也一起改了。

### 10.5 用 Word 做验收（PowerShell，脱离沙箱跑）

```powershell
$word = New-Object -ComObject Word.Application
$word.Visible = $false
$word.DisplayAlerts = 0
$doc = $word.Documents.Open("D:\tmp\out.docx", $false, $true)   # ReadOnly
"pages=" + $doc.ComputeStatistics(2)          # wdStatisticPages
$doc.SaveAs2("D:\tmp\out_word.html", 8)       # wdFormatHTML
$doc.Close(0)
$word.Quit()
```

随后在这份 HTML 里查 `background:#bbf7d0` 与 `Consolas`：前者证明字符底纹生效，后者证明等宽字体生效。

---

## 11. 容易踩的坑（写给未来的维护者）

1. **不要**在界面层做坐标换算。界面只消费 `Block.regions`（块内偏移）和 `Block.genomic_interval`（基因组坐标）。
2. **不要**在负链上先 revcomp 再排序，或只做其中一步。
3. **不要**用 `Segment.start/end` 直接 `seq[start:end]`——那是 1-based 闭区间，切片会错一位。
4. **不要**为了“让长度好看”去裁 CDS 或补 N。
5. **不要**把 `Block.regions` 的偏移当作剪接坐标以外的含义：在 `per_segment` 模式下每个块都从 1 开始。
6. **不要**合并共享的外显子对象：一个 exon 可以同时属于正链与负链转录本，排序方向不同，缓存键必须包含 `tx_id`。
7. 面板初始滚动位置由 `_scroll_to_interesting()` 决定，且滚动值是行高的整数倍——否则会出现半行被切顶的观感问题。
8. **不要**把扫描用的采样上限写成 `def f(max_bytes=_SAMPLE)`——默认值在 `def` 时就绑定了，测试没法改。用 `max_bytes: int = 0` 再在函数体里回填。
9. **不要**对 BGZF（bgzip 压缩的 FASTA）再解一层 gzip：pyfaidx 本来能直接读它；反过来普通 gzip 的 `.fa.gz` 必须解成 `.fa`，否则 `open_fasta` 建不了索引。判据是 gzip extra field 里的 `BC` 子字段。
10. **不要**用字符串 `startswith` 判断"文件在不在库目录里"：`datasets_evil/` 也满足。用 `os.path.commonpath`（并捕获跨盘符的 `ValueError`）。
11. 测试里 `QApplication` 必须**会话级共用一个**（`tests/conftest.py` 的 `qapp`）：创建第二个会直接 `RuntimeError: Please destroy the QApplication singleton`。窗口销毁前还要 `MainWindow.shutdown()` 等后台线程退出，否则 Qt 会崩掉整个进程（实测会把 pytest 打成 SIGTERM）。
12. **不要**用 `except Exception: pass` 兜底再降级——真实事故：pyfaidx 0.9 把 `Fasta.index` 改名成 `records`，包装层读 `.index` 抛 `AttributeError`，被这行吞掉后只剩内存模式的"文件过大，请安装 pyfaidx"，一个 900 MB 的基因组直接打不开，而 pyfaidx 明明装着且可用。降级时**必须把原始原因留在错误信息里**（`open_fasta` 现在会同时报两种模式的原因）。
13. **不要**在 GFF3 属性列上无条件跑 GTF 正则 `([A-Za-z0-9_]+)\s+"([^"]*)"`：无引号时它在每个起始位置贪婪匹配到串尾再回溯，是 O(n²)。NCBI 属性列约 500 字符 → 单次 57 µs，178 万行 = **101 秒**。先 `'"' in text` 一次扫描就能整段跳过（语义严格等价）。
14. **不要**给每个片段复制整个属性字典：167 万个 `Segment` 各复制一份 6~8 键的 dict 约 600 MB，而全工程只有链方向检查读 `strand`。用 `_segment_attrs()` 显式声明保留哪些键。
14.1 **不要**让 `read_rows` 给每一行都建属性字典。167 万个 exon/CDS 行的字典加属性值字符串实测吃 2.4 GB，而这些行后面只用得到 `Parent` 与 `gene_id`。现在只有 gene/transcript 行（约 10 万行）保留完整 `attrs`，其余行的 `attrs` 恒为空——**拿子特征行的 `attrs` 查东西会拿到空字典，这是设计不是 bug**（见 `FeatureRow` 的注释）。
14.2 `seqid` / `ftype` / `strand` 的取值只有几个，却会为每行各建一个字符串（178 万行 ≈ 400 MB）。用 `sys.intern` 共用同一个对象。
15. `urllib.parse.unquote` 也不是免费的：它在无 `%` 的串上仍走生成器，1150 万次 = 45 秒。先判 `"%" in val`。
16. pyfaidx 的索引属性搬过家：0.6/0.7 是 `fa.index`，0.9+ 改成 `fa.records`（真表在 `fa.faidx.index`）。**只认 `rlen`，不要用 `IndexEntry.length`**——那是文件里的字节跨度，不是碱基数。取长度统一走 `_faidx_entry_lengths()`。
17. `monkeypatch.setattr(模块, "类名", 函数)` 之后调用 `类名(path)` 传的是**普通函数**，第一个参数就是 `path`，不能写成 `(self, path)`；否则报 `missing 1 required positional argument`，而那个 TypeError 会被降级逻辑接住，测试看起来是"别的地方"失败。
18. Qt6 的 `QSplitter.restoreState()` 会把**方向**一并写回，所以每个分隔条必须单独保存/恢复状态；把垂直分隔条的状态恢复到水平分隔条上，会让两列变成上下堆叠（本项目已踩过）。
19. **不要**按"一个碱基一个 run"无脑生成 docx。50 kb 的基因组面板会产生 5 万个 run（约 10 MB），
    Word 打开会明显发卡；同色相邻必须合并（`BlockPainter.runs`）。默认着色方式选"按区域"而不是
    "碱基字母色"，原因就在这张体积表（§9.2）。
20. **不要**指望 `w:pPr/w:rPr` 给序列行当默认字体（§10.3）：Word 认，别的阅读器不一定认，
    而序列行丢字体是整体错版，不是"稍微难看一点"。
21. **不要**在按清单指定角色时只给命中的候选赋值——先清零，否则打分残留的角色会多出来（§7.8）。
22. **不要**用改后缀的方式造 `.doc`：把 docx 改名成 `.doc`，Word 会按纯文本打开，用户看到的是
    "一堆乱码"而不是颜色。程序只产出 `.docx` 与 `.html`，两者都能用 Word 打开。
23. **不要**在 `Bash` 工具里 heredoc 一段含 PowerShell 变量（`$doc` 之类）的脚本文本：
    会被安全策略当成"从 Bash 调 PowerShell"直接拦掉，改用 Edit/Write 落盘再执行。
24. **不要**把 `human_size()`（bp 单位，小写 `bp/kb/Mb`）拿去显示文件大小：`136 bp` 会让人
    以为是 136 个碱基的序列，`863.00 Mb` 会被当成 863 兆碱基的基因组（实际是 863 MiB 的 fasta 文件）。
    字节一律走 `human_bytes()`（`136 B / 20.9 KB / 864.0 MB`）。判据很简单：**数字来自
    `os.path.getsize()` 就用 `human_bytes`，来自 `len(seq)` 就用 `human_size`**。
    §5.6 的打包对话框踩过这个坑——参考蛋白（120 aa / 136 字节）当时显示成 `136 bp`。
25. **不要**写 `dlg.Accepted`（实例访问枚举成员）。PySide6 里这在**实例**上会抛
    `AttributeError: 'QDialog' object has no attribute 'Accepted'`——**类访问 `QDialog.Accepted` 却是好的**，
    所以静态看代码"很像对的"，跑起来必炸。统一写 `QDialog.DialogCode.Accepted`。
    真实事故：`if dlg.exec() != dlg.Accepted:` 让彩色文档导出、数据集打包、压缩包导入、
    批量导出、配色设置**五条路全部静默失效**——对话框正常关闭，槽函数在 `exec()` 之后抛异常中断，
    打包成 exe（无控制台）时异常连 stderr 都看不到，用户只看到"点了没反应、没有文件生成"。
    `tests/test_export_gui.py` 现在有两条打桩 `exec()` 的端到端测试（必须真的落盘）
    与一条源码级防退化断言（正则查 `dlg.Accepted`），注入旧写法会立刻变红。
26. **不要**把"记住路径"当成"恢复状态"。路径只是字符串，数据没加载就搜索不了，
    而界面上文件名齐齐全全，用户会认定是软件坏了（§9.4）。判断"是否真的就绪"必须看
    `provider/model` 这些运行时对象，不能看界面文字；恢复逻辑里也要先确认
    `last_gene in model.genes` 再 `select_gene()`。
27. **不要**在离屏测试里让 `_error()` 弹真模态框：`QMessageBox.exec()` 在无窗口环境里
    会一直等下去，表现为"脚本卡死、没有任何输出"（本机实测：日志停在调用 `select_gene`
    那一行，`faulthandler` 才暴露出卡在 `main_window.py` 的 `_error`）。
    `tests/test_export_gui.py` 的 `no_modal` fixture 就是干这个的；手工脚本里也要打桩，
    或者先用 `model.genes` 确认基因 ID 存在。
28. **不要**用 `grep 1,001` 这种方式判断"文档里有没有坐标列"：正文信息行
    （`chr1:1,001-1,800`）本来就带这个数字，会永远命中。坐标 run 的判据是浅灰
    `w:val="9CA3AF"`（docx）与 `<span class="c">`（html）。
29. **不要**只把"与基因同一侧"的两个端点拿来算窗口边界。真实事故：基因组连续窗口写成
    `region_start = min(tx_start, up_s)` / `region_end = max(tx_end, dn_e)`，
    正链恰好正确（上游在左、下游在右），**负链正好用反**——上游在坐标更大一侧被丢掉、
    下游偏移 `dn_s - region_start` 变成负数，`paint_overlays` 把它裁掉，于是窗口里只剩基因本身。
    `prlrb`（NC_067224.1:14,106,166-14,120,974，负链）就是这样"看不到上下游"的：
    上游/下游面板各 2000 bp 都在，唯独基因组窗口长度仍等于基因长度 14,809。
    正确写法是四个端点全参与：`min(tx_start, up_s, dn_s)` / `max(tx_end, up_e, dn_e)`。
    判据不是"面板有没有序列"，而是**窗口长度 == 上游实际 + 基因长度 + 下游实际**；
    `tests/test_extract.py` 的正/负链各一条用例锁住这个等式。
