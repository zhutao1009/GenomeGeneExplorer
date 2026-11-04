"""数据模型。

设计要点
--------
1. ``Segment`` 保存 GFF/GTF 原始坐标（1-based 闭区间），并保留 ``phase`` / 行号，
   使任何结果都能追溯到原始文件。
2. ``Transcript`` 同时保存 ``exons`` 与 ``cds``。两者**不互相推导**，
   缺少 exon 行时由解析器显式补全并写入 warning，避免"看起来正常其实错"。
3. ``Transcript.exons`` 按基因组升序存储；需要转录顺序时调用 ``ordered_exons()``，
   避免把"存储顺序"和"转录顺序"混为一谈。
4. ``Panel`` / ``Block`` / ``Region`` 是纯展示层模型，不含 Qt 依赖，
   因此可以在无 GUI 环境下单元测试。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Tuple

from .coords import order_for_transcript, sort_segments_asc


# --------------------------------------------------------------------------- #
# 注释结构
# --------------------------------------------------------------------------- #
@dataclass
class Segment:
    """一个注释片段（外显子 / CDS / UTR / 密码子）。"""

    start: int                      # 1-based, 含
    end: int                        # 1-based, 含
    phase: int = -1                 # CDS 的 phase/frame：0/1/2，其他为未知
    score: Optional[str] = None
    attrs: Dict[str, str] = field(default_factory=dict)
    line_no: int = -1               # 在注释文件中的行号
    invalid: bool = False           # 原始注释 start > end

    def __post_init__(self) -> None:
        if self.start > self.end:
            # 保留原始值以便报告，但内部按合法区间处理
            self.invalid = True
            self.start, self.end = self.end, self.start

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    @property
    def start0(self) -> int:
        return self.start - 1

    @property
    def end0(self) -> int:
        return self.end

    def overlaps(self, start: int, end: int) -> bool:
        return not (self.end < start or self.start > end)

    def as_interval(self) -> Tuple[int, int]:
        return self.start0, self.end0

    def __str__(self) -> str:  # pragma: no cover - 调试输出
        return f"Segment({self.start}-{self.end}, len={self.length})"


@dataclass
class Transcript:
    """一个转录本（mRNA / ncRNA / 其他 RNA 类型）。"""

    tx_id: str
    gene_id: str
    seqid: str
    strand: str
    source: str = "."
    biotype: str = ""
    exons: List[Segment] = field(default_factory=list)
    cds: List[Segment] = field(default_factory=list)
    start_codon: List[Segment] = field(default_factory=list)
    stop_codon: List[Segment] = field(default_factory=list)
    utr5_features: List[Segment] = field(default_factory=list)
    utr3_features: List[Segment] = field(default_factory=list)
    attrs: Dict[str, str] = field(default_factory=dict)
    line_no: int = -1
    warnings: List[str] = field(default_factory=list)
    exons_derived: bool = False     # 外显子由 CDS 推导而来

    # -- 顺序视图 ---------------------------------------------------------- #
    def ordered_exons(self) -> List[Segment]:
        return order_for_transcript(self.exons, self.strand)

    def ordered_cds(self) -> List[Segment]:
        return order_for_transcript(self.cds, self.strand)

    # -- 区间 -------------------------------------------------------------- #
    def span(self) -> Tuple[int, int]:
        segs = self.exons or self.cds
        if not segs:
            raise ValueError(f"转录本 {self.tx_id} 既无外显子也无 CDS，无法确定区间")
        return min(s.start for s in segs), max(s.end for s in segs)

    def cds_span(self) -> Optional[Tuple[int, int]]:
        if not self.cds:
            return None
        return min(s.start for s in self.cds), max(s.end for s in self.cds)

    @property
    def n_exons(self) -> int:
        return len(self.exons)

    @property
    def n_introns(self) -> int:
        return max(len(self.exons) - 1, 0)

    def mrna_length(self) -> int:
        return sum(s.length for s in self.exons)

    def cds_length(self) -> int:
        return sum(s.length for s in self.cds)

    def display_name(self) -> str:
        """界面显示名：优先 Name / gene_name / product。"""
        for key in ("Name", "gene_name", "gene", "product", "locus_tag"):
            val = self.attrs.get(key)
            if val:
                return f"{self.tx_id} ({val})"
        return self.tx_id

    def add_warning(self, msg: str) -> None:
        if msg not in self.warnings:
            self.warnings.append(msg)


@dataclass
class Gene:
    """一个基因及其所有转录本。"""

    gene_id: str
    gene_name: str
    seqid: str
    strand: str
    start: int
    end: int
    biotype: str = ""
    source: str = "."
    attrs: Dict[str, str] = field(default_factory=dict)
    line_no: int = -1
    transcripts: Dict[str, Transcript] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    @property
    def n_transcripts(self) -> int:
        return len(self.transcripts)

    def span(self) -> Tuple[int, int]:
        """基因区间：优先使用基因行自身坐标，异常时用转录本合并区间。"""
        if self.start <= self.end:
            return self.start, self.end
        bounds = [t.span() for t in self.transcripts.values()]
        if not bounds:
            return self.start, self.end
        return min(b[0] for b in bounds), max(b[1] for b in bounds)

    def tx_list(self) -> List[Transcript]:
        return list(self.transcripts.values())

    def canonical_transcript(self) -> Optional[Transcript]:
        """推荐转录本：CDS 最长者优先，其次外显子总长，最后 ID 排序。

        仅用于"默认选择"，界面必须明确提示已自动选择哪一条，
        并允许用户切换，绝不把多条转录本混合。
        """
        if not self.transcripts:
            return None
        return sorted(
            self.transcripts.values(),
            key=lambda t: (-t.cds_length(), -t.mrna_length(), t.tx_id),
        )[0]

    def add_warning(self, msg: str) -> None:
        if msg not in self.warnings:
            self.warnings.append(msg)


# --------------------------------------------------------------------------- #
# 展示层模型（供 GUI 与导出共用，无 Qt 依赖）
# --------------------------------------------------------------------------- #
REGION_UPSTREAM = "upstream"
REGION_DOWNSTREAM = "downstream"
REGION_EXON = "exon"
REGION_INTRON = "intron"
REGION_CDS = "cds"
REGION_UTR5 = "utr5"
REGION_UTR3 = "utr3"
REGION_FLANK = "flank"
REGION_OTHER = "other"

# 颜色优先级：数字越大越优先，用于重叠注释的着色决策（CDS 覆盖 UTR 覆盖外显子）
REGION_PRIORITY = {
    REGION_OTHER: 0,
    REGION_FLANK: 1,
    REGION_UPSTREAM: 2,
    REGION_DOWNSTREAM: 2,
    REGION_INTRON: 3,
    REGION_EXON: 4,
    REGION_UTR5: 5,
    REGION_UTR3: 5,
    REGION_CDS: 6,
}

REGION_LABELS_CN = {
    REGION_UPSTREAM: "上游",
    REGION_DOWNSTREAM: "下游",
    REGION_EXON: "外显子",
    REGION_INTRON: "内含子",
    REGION_CDS: "CDS",
    REGION_UTR5: "5′ UTR",
    REGION_UTR3: "3′ UTR",
    REGION_FLANK: "侧翼",
    REGION_OTHER: "其他",
}


@dataclass
class Region:
    """块内一段连续的着色区间（0-based 半开），可携带基因组坐标。"""

    start0: int
    end0: int
    kind: str
    label: str = ""
    genomic_start: Optional[int] = None
    genomic_end: Optional[int] = None

    @property
    def length(self) -> int:
        return max(self.end0 - self.start0, 0)


@dataclass
class Block:
    """一个序列块：一段序列 + 其上的着色区间 + 坐标映射信息。"""

    key: str
    title: str
    seq: str
    regions: List[Region] = field(default_factory=list)
    note: str = ""
    meta: Dict[str, object] = field(default_factory=dict)

    @property
    def length(self) -> int:
        return len(self.seq)

    @property
    def strand(self) -> str:
        return str(self.meta.get("strand", "+"))

    @property
    def coord_mode(self) -> str:
        """``genomic``：行号显示基因组坐标；``transcript``：显示转录本内位置。"""
        return str(self.meta.get("coord_mode", "genomic"))

    @property
    def coord_start(self) -> int:
        """坐标轴起点（genomic 模式下为该块的基因组起点或终点，见 coord_reverse）。"""
        return int(self.meta.get("coord_start", 1))

    @property
    def coord_reverse(self) -> bool:
        """负链块：基因组坐标从右向左递减。"""
        return bool(self.meta.get("coord_reverse", False))

    @property
    def genomic_interval(self) -> Optional[Tuple[int, int]]:
        gs = self.meta.get("genomic_start")
        ge = self.meta.get("genomic_end")
        if gs is None or ge is None:
            return None
        return int(gs), int(ge)

    def region_at(self, idx: int) -> Optional[Region]:
        for reg in self.regions:
            if reg.start0 <= idx < reg.end0:
                return reg
        return None

    def first_coord(self, line_idx: int, bases_per_line: int) -> int:
        """第 ``line_idx`` 行第一个碱基对应的坐标（界面左栏与导出文档共用）。

        正链按 ``coord_start`` 递增；负链从 ``coord_start`` 起向小坐标走；
        拼接序列（mRNA / CDS）的坐标是转录本内位置，与基因组无关。
        """
        if self.coord_mode == "transcript" or not self.coord_reverse:
            return self.coord_start + line_idx * bases_per_line
        return self.coord_start + len(self.seq) - 1 - line_idx * bases_per_line


@dataclass
class Panel:
    """一个可独立滚动、复制、导出的序列面板。"""

    key: str
    title: str
    blocks: List[Block] = field(default_factory=list)
    note: str = ""
    seq_type: str = "misc"         # FASTA 头中的序列类型
    empty_text: str = "（无内容）"
    badge: str = ""                # 标题栏右侧的状态标签，如"注释无 UTR（推导）"

    @property
    def sequence(self) -> str:
        return "".join(b.seq for b in self.blocks)

    @property
    def length(self) -> int:
        return sum(b.length for b in self.blocks)

    def is_empty(self) -> bool:
        return self.length == 0

    def kinds_present(self) -> List[str]:
        seen: List[str] = []
        for block in self.blocks:
            for reg in block.regions:
                if reg.kind not in seen:
                    seen.append(reg.kind)
        return seen


@dataclass
class CheckItem:
    """一条自检结果。level: ok / info / warn / error。"""

    level: str
    title: str
    detail: str = ""


@dataclass
class GeneExtraction:
    """一次提取的完整结果（GUI 与导出都以它为唯一数据源）。"""

    # 标识
    gene_id: str
    gene_name: str
    tx_id: str
    seqid: str
    strand: str
    gene_start: int
    gene_end: int

    # 提取区域
    region_start: int = 0
    region_end: int = 0
    upstream_requested: int = 0
    downstream_requested: int = 0
    upstream_actual: int = 0
    downstream_actual: int = 0
    upstream_clipped: bool = False
    downstream_clipped: bool = False

    # 序列
    upstream_seq: str = ""
    downstream_seq: str = ""
    mrna: str = ""
    cds: str = ""
    utr5: str = ""
    utr3: str = ""
    # 注释文件里是否真的给了 UTR 行（RefSeq/NCBI 的 GFF 通常为 False）。
    # 无论有无，本程序的 UTR 都是按「外显子 − CDS」推导的，这个标记只用于如实告知来源。
    utr_in_annotation: bool = False

    # 面板
    panels: Dict[str, Panel] = field(default_factory=dict)

    # 诊断
    warnings: List[str] = field(default_factory=list)
    checks: List[CheckItem] = field(default_factory=list)

    # 版本追溯
    fasta_path: str = ""
    annotation_path: str = ""
    tx_count: int = 1
    n_exons: int = 0
    n_introns: int = 0
    n_cds: int = 0

    def panel(self, key: str) -> Optional[Panel]:
        return self.panels.get(key)

    def warning_text(self) -> str:
        return "\n".join(f"• {w}" for w in self.warnings)

    def summary_line(self) -> str:
        return (
            f"{self.gene_name} ({self.gene_id}) | {self.tx_id} | "
            f"{self.seqid}:{self.gene_start:,}-{self.gene_end:,}({self.strand}) | "
            f"外显子 {self.n_exons} | 内含子 {self.n_introns} | "
            f"mRNA {len(self.mrna):,} bp | CDS {len(self.cds):,} bp"
        )
