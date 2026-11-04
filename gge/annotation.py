"""GFF3 / GTF 注释解析。

支持范围
--------
* GFF3（``##gff-version 3``、``ID=``/``Parent=`` 属性、逗号分隔的多 Parent）
* GFF2 / GTF（``gene_id "x"; transcript_id "y";`` 属性、``frame`` 列）
* 自动识别格式；属性解析同时兼容两种方言，识别失败会记录 warning 而不是崩溃
* ``.gz`` 压缩文件直接读取
* ``##FASTA`` 段自动停止

层级重建
--------
真实注释文件经常不"标准"，解析器显式处理以下情况，并在 ``warnings`` /
``Transcript.warnings`` 中留痕：

1. 只有 gene / transcript 行，没有 exon 行 → 由 CDS 推导外显子，标记
   ``exons_derived=True`` 并提示"外显子范围可能不完整（缺 UTR 部分）"。
2. 有 exon 行、没有 CDS 行 → 视为非编码转录本，CDS 面板显示为空并说明。
3. exon 直接挂在 gene 下（没有 transcript 层）→ 合成一个伪转录本，
   提示"注释缺少转录本层级"。
4. ``Parent`` 指向不存在的 ID → 回退到 ``transcript_id`` / ``gene_id`` 属性，
   仍失败则丢弃该行并计数。
5. 基因行坐标与子特征不一致 → 以子特征合并区间为准并提示。

坐标
----
所有 ``Segment`` 保存 GFF/GTF 原始的 1-based 闭区间；解析阶段不做任何
0-based 转换，避免"转换两次"的经典错误。
"""

from __future__ import annotations

import gzip
import os
import re
import sys
from collections import Counter, OrderedDict
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Tuple
from urllib.parse import unquote

from .errors import AnnotationError, CancelledError, FileFormatError
from .models import Gene, Segment, Transcript

GFF_SUFFIXES = (".gff3", ".gff", ".gtf", ".gff.gz", ".gtf.gz", ".gff3.gz")

# 只保留这些属性，显著降低大注释文件的内存占用
KEEP_ATTRS = {
    "id", "parent", "name", "gene_id", "transcript_id", "gene_name", "gene",
    "gene_biotype", "transcript_biotype", "biotype", "locus_tag", "product",
    "gbkey", "description", "protein_id", "db_xref", "type", "class",
}

GENE_TYPES = {"gene", "pseudogene", "ncrna_gene", "protein_coding_gene"}
TRANSCRIPT_TYPES = {
    "mrna", "transcript", "ncrna", "lncrna", "lincrna", "trna", "rrna", "mirna",
    "snrna", "snorna", "scrna", "pre_mirna", "primary_transcript",
    "pseudogenic_transcript", "processed_transcript", "ncrna_primary_transcript",
    "antisense_rna", "guide_rna", "rnase_p_rna", "srp_rna", "telomerase_rna",
    "vault_rna", "y_rna", "c_gene_segment", "v_gene_segment", "d_gene_segment",
    "j_gene_segment", "unconfirmed_transcript", "Transcript", "mRNA",
    "five_prime_utr",  # 占位，避免被误判（下方单独处理）
}
TRANSCRIPT_TYPES.discard("five_prime_utr")

EXON_TYPES = {"exon"}
CDS_TYPES = {"cds"}
UTR5_TYPES = {"five_prime_utr", "5utr", "5'utr", "utr5"}
UTR3_TYPES = {"three_prime_utr", "3utr", "3'utr", "utr3"}
START_CODON_TYPES = {"start_codon"}
STOP_CODON_TYPES = {"stop_codon"}
# 明确忽略的结构型 feature
IGNORED_TYPES = {
    "region", "chromosome", "supercontig", "contig", "match", "match_part",
    "sequence_feature", "gap", "repeat_region", "direct_repeat", "intron",
    "intergenic_region", "target_site_duplication", "proviral_region",
}


class FeatureRow:
    """一行的紧凑表示（__slots__ 降低内存开销）。

    内存说明（NCBI 草鱼注释 178 万行实测）
    --------------------------------------
    最原始的实现给每一行都建一个属性字典（NCBI 的 exon 行有 6~8 个键），
    ``read_rows`` 一步就吃 **2.4 GB**，14.9 GB / 剩余 3 GB 的机器上直接换页卡死。
    其中 167 万行是子特征（exon/CDS/UTR/密码子），而后续步骤只需要它们的
    ``Parent``（→ ``tx_key``）与 ``gene_id``（→ ``gene_key``），其余属性从来没人读。

    所以现在按类型区别对待：

    * ``gene`` / ``transcript`` 行（约 10 万行）：保留完整 ``attrs``，
      因为这些属性要变成 Gene/Transcript 的 attrs，是展示与导出的一部分；
    * 其他行（约 167 万行）：``attrs`` 一律为空字典，只留 ``tx_key`` / ``gene_key``。

    拿子特征行的 ``attrs`` 去查东西会拿到空字典——这不是 bug，是设计。
    重构前请先读这段。
    """

    __slots__ = ("seqid", "source", "ftype", "start", "end", "strand",
                 "phase", "attrs", "line_no", "invalid", "kind",
                 "tx_key", "gene_key")

    def __init__(self, seqid, source, ftype, start, end, strand, phase, attrs, line_no,
                 kind="", tx_key="", gene_key=""):
        self.seqid = seqid
        self.source = source
        self.ftype = ftype
        self.start = start
        self.end = end
        self.strand = strand
        self.phase = phase
        self.attrs = attrs
        self.line_no = line_no
        self.invalid = start > end
        if self.invalid:
            self.start, self.end = end, start
        self.kind = kind
        self.tx_key = tx_key
        self.gene_key = gene_key


def classify(ftype: str) -> str:
    """把 feature 类型归一化到内部类别。"""
    t = ftype.strip().lower()
    if t in GENE_TYPES or t.endswith("_gene") or t in ("gene",):
        return "gene"
    if t in EXON_TYPES:
        return "exon"
    if t in CDS_TYPES:
        return "cds"
    if t in UTR5_TYPES:
        return "utr5"
    if t in UTR3_TYPES:
        return "utr3"
    if t in START_CODON_TYPES:
        return "start_codon"
    if t in STOP_CODON_TYPES:
        return "stop_codon"
    if t in TRANSCRIPT_TYPES:
        return "transcript"
    if t in IGNORED_TYPES:
        return "ignored"
    if "utr" in t:
        return "utr5" if "five" in t or t.startswith("5") else "utr3"
    if "rna" in t or "transcript" in t:
        return "transcript"
    return "other"


_GTF_ATTR_RE = re.compile(r'([A-Za-z0-9_]+)\s+"([^"]*)"')


def parse_attributes(text: str) -> Dict[str, str]:
    """同时兼容 GFF3 与 GTF 的属性列。

    判定顺序：优先按 GFF3 的 ``key=value`` 解析；若不存在 ``=`` 则按 GTF 的
    ``key "value"`` 解析；两者混合时逐个属性判定。

    性能注意（真实踩坑，NCBI 草鱼 GFF 178 万行实测）
    ---------------------------------------------------
    GTF 正则 ``([A-Za-z0-9_]+)\\s+"([^"]*)"`` 在**不含引号**的 GFF3 属性串上会退化：
    每个起始位置都要把 ``[A-Za-z0-9_]+`` 贪婪匹配到串尾再回溯，对一条 500 字符的
    NCBI 属性列就是 O(n²)，单次 57 µs，178 万次 = **101 秒**。
    而 GFF3 属性列里根本不会有 ``"``，所以先用 ``'"' in text`` 一次字符串扫描把
    正则整个跳过（语义完全等价，正则不可能在没有引号的串上匹配成功）。

    同理 ``unquote`` 只对含 ``%`` 的值调用：NCBI 的值绝大多数不含转义，
    但 unquote 每次都走生成器，1150 万次 = 45 秒。
    """
    text = text.strip()
    if not text or text == ".":
        return {}
    attrs: Dict[str, str] = {}
    if "=" in text and ('"' not in text or not _GTF_ATTR_RE.search(text)):
        for part in text.split(";"):
            if "=" not in part:
                continue
            key, _, val = part.partition("=")
            key = key.strip().lower()
            if key and key in KEEP_ATTRS:
                val = val.strip()
                attrs[key] = unquote(val) if "%" in val else val
        return attrs

    # GTF / GFF2
    for m in _GTF_ATTR_RE.finditer(text):
        key = m.group(1).strip().lower()
        if key in KEEP_ATTRS:
            attrs[key] = m.group(2)
    if attrs:
        return attrs

    # 无引号的 GFF2：key value; key value;
    for part in text.split(";"):
        bits = part.strip().split(None, 1)
        if len(bits) == 2:
            key = bits[0].lower()
            if key in KEEP_ATTRS:
                attrs[key] = bits[1].strip()
    return attrs


def _segment_attrs(attrs: Dict[str, str]) -> Dict[str, str]:
    """外显子/CDS 片段只需要 ``strand`` 一个属性。

    原先这里是 ``dict(row.attrs)``，把整条 GFF3 属性串（NCBI 的 exon 行有 6~8 个键）
    复制进每一个片段。一份全基因组注释有 167 万个片段，仅此一项就吃掉约 600 MB，
    而全工程唯一的消费者是 ``parse_annotation`` 里的链方向一致性检查（只读 strand）。

    说明：``KEEP_ATTRS`` 目前**不含** ``strand``（GFF3 的链方向在第 7 列，不在属性里），
    所以这里通常返回空字典——这与优化前 ``dict(row.attrs)`` 的结果一致，
    那段一致性检查本来就是"有则用、无则跳过"。真要让它生效，得先把 strand 加进
    ``KEEP_ATTRS``，那是另一件事，别在这里偷偷改。
    """
    strand = attrs.get("strand")
    return {"strand": strand} if strand else {}


def detect_format(path: str) -> str:
    """返回 'gff3' / 'gtf' / 'unknown'。"""
    lower = path.lower()
    if lower.endswith((".gtf", ".gtf.gz")):
        return "gtf"
    if lower.endswith((".gff3", ".gff3.gz")):
        return "gff3"
    # 内容嗅探
    opener = gzip.open if lower.endswith(".gz") else open
    gff3_hint = 0
    gtf_hint = 0
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            if i > 400:
                break
            if line.startswith("##gff-version 3"):
                return "gff3"
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            attrs = cols[8]
            if re.search(r"\bID\s*=", attrs) or re.search(r"\bParent\s*=", attrs):
                gff3_hint += 1
            if _GTF_ATTR_RE.search(attrs):
                gtf_hint += 1
    if gff3_hint > gtf_hint:
        return "gff3"
    if gtf_hint > 0:
        return "gtf"
    return "unknown"


def _open_text(path: str):
    if path.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def read_rows(path: str, progress_cb=None, max_rows: int = 0,
              cancel_check=None) -> List[FeatureRow]:
    """逐行读取并解析成 FeatureRow 列表。

    ``cancel_check`` 每 5 万行查一次：500 MB 的 RefSeq 注释要读 70 s，
    没有这个钩子的话"取消"按钮对它完全无效（线程跑完才停）。
    """
    rows: List[FeatureRow] = []
    n_line = 0
    n_bad = 0
    bad_samples: List[str] = []
    opener = _open_text(path)
    # 进度按「已读字符数 / 文件字节数」估算。以前这里固定报 40%，而 NCBI 全基因组
    # GFF 的读取要一两分钟，进度条停着不动会被当成卡死。
    # .gz 的压缩比与解压后的行数无关，没法这样估，退回只报行数。
    try:
        total_chars = 0 if path.lower().endswith(".gz") else os.path.getsize(path)
    except OSError:
        total_chars = 0
    seen_chars = 0
    with opener as fh:
        for line_no, line in enumerate(fh, start=1):
            n_line = line_no
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\r\n").split("\t")
            if len(cols) < 8:
                n_bad += 1
                if len(bad_samples) < 5:
                    bad_samples.append(f"第 {line_no} 行只有 {len(cols)} 列")
                continue
            if len(cols) == 8:
                cols.append(".")  # 缺属性列时补空，宽容处理
            seqid = sys.intern(cols[0].strip())
            source = sys.intern(cols[1].strip())
            ftype = sys.intern(cols[2].strip())
            try:
                start = int(cols[3])
                end = int(cols[4])
            except ValueError:
                n_bad += 1
                if len(bad_samples) < 5:
                    bad_samples.append(f"第 {line_no} 行坐标不是整数：{cols[3]}/{cols[4]}")
                continue
            strand = cols[6].strip()
            if strand not in ("+", "-", "."):
                strand = "."
            else:
                strand = sys.intern(strand)      # 只有 3 种取值，别为每行建一个字符串
            phase_raw = cols[7].strip()
            phase = int(phase_raw) if phase_raw in ("0", "1", "2") else -1
            # 按类型决定保留多少属性（见 FeatureRow 的内存说明）：
            # gene/transcript 行占 6%，保留整字典；其余 94% 只留两个 key。
            kind = classify(cols[2])
            if kind in ("gene", "transcript"):
                attrs = parse_attributes(cols[8])
                tx_key = gene_key = ""
            else:
                parsed = parse_attributes(cols[8])
                attrs = {}
                # 与 parse_annotation 里原来的取法逐字一致：
                # 子特征优先 Parent，没有 Parent 时才用 GTF 的 transcript_id；
                # "other" 类只认 Parent（原实现就是这样，"other" 不看 transcript_id）。
                parent = parsed.get("parent") or ""
                if parent:
                    tx_key = parent.split(",")[0].strip()
                elif kind in ("exon", "cds", "utr5", "utr3", "start_codon", "stop_codon"):
                    tx_key = parsed.get("transcript_id") or ""
                else:
                    tx_key = ""
                gene_key = parsed.get("gene_id") or parsed.get("gene") or ""
            rows.append(FeatureRow(seqid, source, ftype, start, end, strand,
                                   phase, attrs, line_no, kind, tx_key, gene_key))
            seen_chars += len(line)
            if line_no % 50000 == 0:
                if cancel_check and cancel_check():
                    raise CancelledError()
                if progress_cb:
                    if total_chars > 0:
                        pct = 15 + int(55 * min(1.0, seen_chars / float(total_chars)))
                    else:
                        pct = 40
                    progress_cb(min(pct, 70), f"已解析 {line_no:,} 行注释 …")
            if max_rows and len(rows) >= max_rows:
                break
    if n_bad:
        # 行级问题不影响整体解析，但要计入统计
        read_rows.last_bad_count = n_bad  # type: ignore[attr-defined]
        read_rows.last_bad_samples = bad_samples  # type: ignore[attr-defined]
    else:
        read_rows.last_bad_count = 0  # type: ignore[attr-defined]
        read_rows.last_bad_samples = []  # type: ignore[attr-defined]
    return rows


@dataclass
class AnnotationModel:
    """注释文件的内存模型。"""

    path: str = ""
    fmt: str = "unknown"
    genes: "OrderedDict[str, Gene]" = field(default_factory=OrderedDict)
    tx_index: Dict[str, str] = field(default_factory=dict)     # tx_id -> gene_id
    name_index: Dict[str, List[str]] = field(default_factory=dict)
    stats: Dict[str, int] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)
    ignored_types: Counter = field(default_factory=Counter)
    seqids: List[str] = field(default_factory=list)

    # -- 查询 -------------------------------------------------------------- #
    def get_gene(self, gene_id: str) -> Optional[Gene]:
        return self.genes.get(gene_id)

    def gene_by_name(self, name: str) -> List[Gene]:
        out = []
        for gid in self.name_index.get(name.lower(), []):
            g = self.genes.get(gid)
            if g:
                out.append(g)
        return out

    def gene_for_transcript(self, tx_id: str) -> Optional[Gene]:
        gid = self.tx_index.get(tx_id)
        return self.genes.get(gid) if gid else None

    def all_transcripts(self) -> Iterable[Transcript]:
        for gene in self.genes.values():
            for tx in gene.transcripts.values():
                yield tx

    @property
    def n_genes(self) -> int:
        return len(self.genes)

    @property
    def n_transcripts(self) -> int:
        return len(self.tx_index)

    def add_warning(self, msg: str) -> None:
        if msg not in self.warnings:
            self.warnings.append(msg)

    def stat_text(self) -> str:
        s = self.stats
        return (
            f"基因 {s.get('genes', 0):,} | 转录本 {s.get('transcripts', 0):,} | "
            f"外显子 {s.get('exons', 0):,} | CDS 片段 {s.get('cds', 0):,} | "
            f"序列 {len(self.seqids):,} 条"
        )


def _attr(attrs: Dict[str, str], *keys: str, default: str = "") -> str:
    for k in keys:
        v = attrs.get(k.lower())
        if v:
            return v
    return default


def _dedup_segments(segs: List["Segment"]) -> int:
    """合并区间完全相同的重复记录，返回被丢弃的条数。

    真实注释里偶尔会把同一条 exon 重复写两遍（合并注释文件时的常见事故）。
    如果不去重，mRNA 长度会凭空多出一整段，后面所有校验都对不上，
    所以这里按 (start, end) 去重，并把条数报给用户。
    """
    seen = set()
    kept: List["Segment"] = []
    for s in segs:
        key = (int(s.start), int(s.end))
        if key in seen:
            continue
        seen.add(key)
        kept.append(s)
    dropped = len(segs) - len(kept)
    if dropped:
        segs[:] = kept
    return dropped


def parse_annotation(path: str, progress_cb=None, cancel_check=None) -> AnnotationModel:
    """解析 GFF3/GTF 并重建 gene → transcript → exon/CDS 层级。"""
    if not os.path.exists(path):
        raise FileFormatError(f"注释文件不存在：{path}", "确认路径是否正确。")
    if os.path.getsize(path) == 0:
        raise FileFormatError(f"注释文件大小为 0：{path}", "重新下载该注释文件。")

    if progress_cb:
        progress_cb(5, "识别注释格式 …")
    fmt = detect_format(path)
    if progress_cb:
        progress_cb(15, f"格式判定：{fmt.upper()}，开始逐行解析 …")

    rows = read_rows(path, progress_cb=progress_cb, cancel_check=cancel_check)
    if cancel_check and cancel_check():        # 读完到建层级之间还有几十秒
        raise CancelledError()
    if not rows:
        raise AnnotationError(
            f"{os.path.basename(path)} 中没有解析到任何特征行。",
            "确认文件是 GFF3/GTF（Tab 分隔、9 列），且不是空文件。",
        )
    bad = getattr(read_rows, "last_bad_count", 0)
    bad_samples = getattr(read_rows, "last_bad_samples", [])
    if progress_cb:
        progress_cb(72, "建立基因/转录本层级 …")

    model = AnnotationModel(path=path, fmt=fmt)
    if bad:
        model.add_warning(f"{bad} 行因列数不足或坐标非整数被跳过，示例：{'；'.join(bad_samples)}")

    # ------------------------------------------------------------------ #
    # 第一步：分类并登记 gene / transcript 骨架
    # ------------------------------------------------------------------ #
    gene_rows: Dict[str, FeatureRow] = {}
    tx_rows: Dict[str, FeatureRow] = {}
    tx_parent: Dict[str, str] = {}
    child_rows: List[Tuple[str, FeatureRow, str, str]] = []  # (kind, row, tx_key, gene_key)
    counts: Counter = Counter()

    for row in rows:
        kind = row.kind or classify(row.ftype)     # read_rows 已算好，别再算一遍
        counts[kind] += 1
        if kind == "ignored":
            model.ignored_types[row.ftype] += 1
            continue

        attrs = row.attrs
        if kind == "gene":
            gid = _attr(attrs, "id", "gene_id", "name", "locus_tag")
            if not gid:
                gid = f"gene:{row.seqid}:{row.start}-{row.end}"
                model.add_warning(
                    f"第 {row.line_no} 行 gene 缺少 ID/gene_id/Name，已生成临时 ID {gid}"
                )
            gene_rows.setdefault(gid, row)
            continue

        if kind == "transcript":
            tid = _attr(attrs, "id", "transcript_id", "name")
            if not tid:
                tid = f"tx:{row.seqid}:{row.start}-{row.end}"
                model.add_warning(
                    f"第 {row.line_no} 行 transcript 缺少 ID/transcript_id，已生成临时 ID {tid}"
                )
            tx_rows.setdefault(tid, row)
            parent = _attr(attrs, "parent", "gene_id", "gene")
            if parent:
                tx_parent[tid] = parent.split(",")[0].strip()
            continue

        if kind in ("exon", "cds", "utr5", "utr3", "start_codon", "stop_codon"):
            if kind == "cds" and row.phase == -1 and row.ftype.lower() == "cds":
                # phase 为 '.'：不完整模型常见，交给 extract 报告
                pass
            if kind == "cds" and row.invalid:
                model.add_warning(
                    f"第 {row.line_no} 行 CDS 起点大于终点（{row.end}>{row.start}），已按交换后处理"
                )
            child_rows.append((kind, row, row.tx_key, row.gene_key))
            continue

        # other：可能是自定义 feature，尝试当作子特征挂在所属转录本下
        if row.tx_key:
            child_rows.append(("other", row, row.tx_key, row.gene_key))
        else:
            model.ignored_types[row.ftype] += 1

    # ------------------------------------------------------------------ #
    # 第二步：建立 Transcript 对象
    # ------------------------------------------------------------------ #
    transcripts: "OrderedDict[str, Transcript]" = OrderedDict()

    # 先由 transcript 行建立
    for tid, row in tx_rows.items():
        gid = tx_parent.get(tid) or _attr(row.attrs, "gene_id", "gene") or gid_from_tx(tid)
        strand = row.strand if row.strand in ("+", "-") else "."
        tx = Transcript(
            tx_id=tid, gene_id=gid, seqid=row.seqid, strand=strand,
            source=row.source, biotype=row.ftype,
            attrs=dict(row.attrs), line_no=row.line_no,
        )
        transcripts[tid] = tx

    # 子特征指向的 transcript 若不存在（只有 exon 行），补建
    synthesized_tx: List[str] = []
    for kind, row, tx_key, gene_key in child_rows:
        if not tx_key:
            continue
        if tx_key not in transcripts:
            gid = gene_key or gid_from_tx(tx_key)
            strand = row.strand if row.strand in ("+", "-") else "."
            tx = Transcript(
                tx_id=tx_key, gene_id=gid, seqid=row.seqid, strand=strand,
                source=row.source, biotype="transcript",
                attrs={}, line_no=row.line_no,
            )
            tx.add_warning(
                f"转录本 {tx_key} 在注释中没有独立的 transcript/mRNA 行，"
                f"由子特征（第 {row.line_no} 行）反推得到"
            )
            transcripts[tx_key] = tx
            synthesized_tx.append(tx_key)

    # ------------------------------------------------------------------ #
    # 第三步：把子特征挂到转录本
    # ------------------------------------------------------------------ #
    orphans = 0
    orphan_samples: List[str] = []
    if progress_cb:
        progress_cb(85, "把外显子 / CDS 挂到转录本 …")
    for n_done, (kind, row, tx_key, gene_key) in enumerate(child_rows, start=1):
        if n_done % 200000 == 0 and cancel_check and cancel_check():
            raise CancelledError()
        target: Optional[Transcript] = transcripts.get(tx_key) if tx_key else None
        if target is None and gene_key:
            # 有 gene 但无 transcript：挂在基因下唯一/首个转录本
            cands = [t for t in transcripts.values() if t.gene_id == gene_key]
            if len(cands) == 1:
                target = cands[0]
            elif len(cands) > 1:
                # 多转录本时无法判定归属，宁可丢弃也不猜
                orphans += 1
                if len(orphan_samples) < 5:
                    orphan_samples.append(f"第 {row.line_no} 行（gene={gene_key} 有 {len(cands)} 个转录本）")
                continue
        if target is None:
            orphans += 1
            if len(orphan_samples) < 5:
                orphan_samples.append(f"第 {row.line_no} 行（Parent={tx_key or '空'}）")
            continue

        seg = Segment(start=row.start, end=row.end, phase=row.phase,
                      score=None, attrs=_segment_attrs(row.attrs), line_no=row.line_no,
                      invalid=row.invalid)
        if kind == "exon":
            target.exons.append(seg)
        elif kind == "cds":
            target.cds.append(seg)
        elif kind == "utr5":
            target.utr5_features.append(seg)
        elif kind == "utr3":
            target.utr3_features.append(seg)
        elif kind == "start_codon":
            target.start_codon.append(seg)
        elif kind == "stop_codon":
            target.stop_codon.append(seg)

    if orphans:
        model.add_warning(
            f"{orphans} 个特征无法确定所属转录本（缺少 Parent/transcript_id 或父级不唯一），"
            f"已跳过。示例：{'；'.join(orphan_samples)}"
        )
    # 中间行对象到此用完（每条 exon/CDS 行一个，NCBI 草鱼注释 178 万个）。
    # 后面几步只操作 Segment / Transcript，及时释放可以把 GB 级峰值降下来。
    del child_rows
    # 游离的 Parent（典型情况：ID 拼写错误或引用了别的注释文件）
    dangling = [tid for tid in synthesized_tx
                if transcripts[tid].gene_id not in gene_rows]
    if dangling:
        model.add_warning(
            f"{len(dangling)} 个子特征的 Parent 在注释中不存在，已单独建立转录本记录："
            f"{', '.join(dangling[:5])}{' 等' if len(dangling) > 5 else ''}。"
            f"请核对 ID 是否正确（这类问题会导致序列取自错误的转录本）"
        )

    # ------------------------------------------------------------------ #
    # 第四步：补全与自检
    # ------------------------------------------------------------------ #
    if progress_cb:
        progress_cb(92, "补全外显子并自检 …")
    for tx in transcripts.values():
        dup_exon = _dedup_segments(tx.exons)
        dup_cds = _dedup_segments(tx.cds)
        tx.exons.sort(key=lambda s: (s.start, s.end))
        tx.cds.sort(key=lambda s: (s.start, s.end))
        if dup_exon:
            tx.add_warning(
                f"转录本 {tx.tx_id} 有 {dup_exon} 条重复的 exon 记录（区间完全相同），"
                f"已按一条计算；否则 mRNA 长度会凭空变长"
            )
        if dup_cds:
            tx.add_warning(
                f"转录本 {tx.tx_id} 有 {dup_cds} 条重复的 CDS 记录（区间完全相同），已按一条计算"
            )
        if not tx.exons and tx.cds:
            tx.exons = [
                Segment(start=s.start, end=s.end, phase=-1, line_no=s.line_no)
                for s in tx.cds
            ]
            tx.exons_derived = True
            tx.add_warning("注释缺少 exon 行，外显子区间由 CDS 推导，可能缺少 UTR 部分")
        if tx.exons and not tx.cds:
            # 非编码转录本，合法；不报错
            pass
        if not tx.exons and not tx.cds:
            tx.add_warning("该转录本在注释中既没有 exon 也没有 CDS，无法提取序列")
        strand = tx.strand
        if any(s.invalid for s in tx.exons + tx.cds):
            tx.add_warning("存在起点大于终点的片段，已按交换后坐标处理并在报告中标注")
        # 链方向一致性
        child_strands = {s.attrs.get("strand") for s in tx.exons} - {None, ""}
        if strand == "." and child_strands:
            tx.strand = next(iter(child_strands))

    # ------------------------------------------------------------------ #
    # 第五步：建立 Gene
    # ------------------------------------------------------------------ #
    # 由 transcript 反推 gene（gene 行可能缺失或 ID 不同）
    by_gene: Dict[str, List[Transcript]] = {}
    for tx in transcripts.values():
        by_gene.setdefault(tx.gene_id, []).append(tx)

    for gid, row in gene_rows.items():
        strand = row.strand if row.strand in ("+", "-") else "."
        gene = Gene(
            gene_id=gid,
            gene_name=_attr(row.attrs, "name", "gene_name", "gene", "locus_tag", default=gid),
            seqid=row.seqid, strand=strand, start=row.start, end=row.end,
            biotype=_attr(row.attrs, "biotype", "gene_biotype", "type", default=row.ftype),
            source=row.source, attrs=dict(row.attrs), line_no=row.line_no,
        )
        model.genes[gid] = gene

    # 没有 gene 行的转录本，合成基因
    for gid, txs in by_gene.items():
        if gid in model.genes:
            continue
        first = txs[0]
        try:
            gstart, gend = min(t.span()[0] for t in txs), max(t.span()[1] for t in txs)
        except ValueError:
            continue
        gene = Gene(
            gene_id=gid, gene_name=gid, seqid=first.seqid,
            strand=first.strand if first.strand in ("+", "-") else ".",
            start=gstart, end=gend, biotype="", source=first.source,
            attrs={}, line_no=first.line_no,
        )
        gene.add_warning("注释中缺少 gene 行，基因记录由转录本反推得到")
        model.genes[gid] = gene

    # 挂载转录本
    for tx in transcripts.values():
        gene = model.genes.get(tx.gene_id)
        if gene is None:
            gene = Gene(
                gene_id=tx.gene_id, gene_name=tx.gene_id, seqid=tx.seqid,
                strand=tx.strand if tx.strand in ("+", "-") else ".",
                start=0, end=0, source=tx.source, line_no=tx.line_no,
            )
            model.genes[tx.gene_id] = gene
            gene.add_warning("基因记录由转录本反推得到")
        gene.transcripts[tx.tx_id] = tx

    # 丢弃空基因（既无转录本也无坐标）
    for gid in [g for g, gene in model.genes.items() if not gene.transcripts]:
        del model.genes[gid]

    # ------------------------------------------------------------------ #
    # 第六步：一致性检查 + 索引
    # ------------------------------------------------------------------ #
    for gene in model.genes.values():
        if gene.transcripts:
            try:
                tstart = min(t.span()[0] for t in gene.transcripts.values())
                tend = max(t.span()[1] for t in gene.transcripts.values())
            except ValueError:
                continue
            if gene.start <= 0 or gene.end <= 0:
                gene.start, gene.end = tstart, tend
            elif tstart < gene.start or tend > gene.end:
                # 子特征超出基因区间：以子特征为准并提示（真实注释常见）
                if abs(tstart - gene.start) > 0 or abs(tend - gene.end) > 0:
                    gene.add_warning(
                        f"基因坐标 {gene.start}-{gene.end} 未完全包含其转录本"
                        f"（{tstart}-{tend}），显示与提取以转录本坐标为准"
                    )
                gene.start, gene.end = min(gene.start, tstart), max(gene.end, tend)
        for tx in gene.transcripts.values():
            if tx.seqid != gene.seqid:
                tx.add_warning(
                    f"转录本位于 {tx.seqid}，与基因所在 {gene.seqid} 不一致"
                )
            if gene.strand in ("+", "-") and tx.strand in ("+", "-") and tx.strand != gene.strand:
                tx.add_warning(
                    f"转录本链方向 {tx.strand} 与基因链方向 {gene.strand} 不一致"
                )
            if tx.strand == ".":
                tx.strand = gene.strand if gene.strand in ("+", "-") else "+"
                tx.add_warning("转录本链方向缺失，已采用基因链方向")

    for gid, gene in model.genes.items():
        model.name_index.setdefault(gene.gene_name.lower(), []).append(gid)
        model.name_index.setdefault(gid.lower(), []).append(gid)
        for key in ("locus_tag", "name", "gene"):
            v = gene.attrs.get(key)
            if v and v.lower() not in model.name_index:
                model.name_index[v.lower()] = [gid]
        for tid, tx in gene.transcripts.items():
            model.tx_index[tid] = gid

    model.seqids = sorted({g.seqid for g in model.genes.values()})
    model.stats = {
        "genes": len(model.genes),
        "transcripts": len(model.tx_index),
        "exons": sum(len(t.exons) for t in model.all_transcripts()),
        "cds": sum(len(t.cds) for t in model.all_transcripts()),
    }
    if counts.get("other"):
        model.add_warning(
            f"{counts['other']} 个未知类型的特征行被忽略（类型示例："
            f"{', '.join(list(model.ignored_types)[:5]) or '无'}）"
        )
    if progress_cb:
        progress_cb(100, "注释解析完成")
    return model


def gid_from_tx(tx_id: str) -> str:
    """从转录本 ID 猜测基因 ID（仅在注释完全缺失基因层级时使用）。"""
    for sep in (".", "-", "_"):
        if sep in tx_id:
            head = tx_id.rsplit(sep, 1)[0]
            if head:
                return head
    return tx_id
