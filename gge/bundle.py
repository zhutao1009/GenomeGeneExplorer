"""压缩包导入：从一个 zip / tar 包里自动认出"哪个文件是基因组、哪个是注释、哪个是参考序列"。

设计原则
--------
1. **先看内容，再看文件名**。文件名只是佐证：很多实际发布的包叫
   ``data/seqs.fa``、``annotation_v2.gff3``、``refs/all_transcripts.fasta``，
   光靠扩展名和非 ASCII 关键字判断很容易把 mRNA 当基因组。
2. **不靠"文件顺序"**，也不靠"第一个 FATSA 就是基因组"这种约定；每个候选都独立打分，
   分数最高者当角色，分数接近时标记为"需确认"，界面上必须由人点确认。
3. **不静默吞异常**。认不出的文件进 ``ignored`` 列表并给出原因，缺了必填项会说缺哪个。

判定依据（可解释、可回溯）
--------------------------
* 文本嗅探：首个非注释行是 ``>`` 就是 FASTA，9 列 Tab 就是 GFF/GTF；
* 字母表：核苷酸序列字母几乎全部落在 ACGTUN，蛋白序列会大量出现 E/F/I/L/P/Q；
* 序列名：FASTA 标题 token 命中注释的 seqid → 基因组；命中 transcript_id → 参考序列；
* 长度拟合：把样本记录长度与注释里各转录本的 exon 之和 / CDS 之和对比，
  哪个吻合得更多就是哪一种（这是区分 mRNA 与 CDS 最可靠的一招）；
* 记录数：与注释的转录本数量同量级 → 参考序列；远小于 → 基因组；
* 文件名关键字：仅作加分项（30 分），单靠它不会决定角色。

所有打分都会以中文写进 ``FileCandidate.evidence``，界面直接展示，方便用户判断。

数据集清单（反向能力）
----------------------
"猜"只在别人给的包上才需要。本程序自己导出的包（见 :func:`write_dataset_archive`）
根目录带一份 ``gge_dataset.json``，写明每个文件是哪个角色；扫描时若发现它，
就以清单为准（``_apply_manifest``），不再让打分去猜，也不再为了显示统计量把
几百 MB 的注释整份解析一遍——清单里已经存了基因/转录本/序列条数。
清单版本不匹配或格式不对时静默退回"按内容判断"，不会因为清单坏了打不开包。
"""

from __future__ import annotations

import bz2
import gzip
import io
import lzma
import os
import re
import shutil
import tarfile
import time
import zipfile
from collections import Counter
from dataclasses import dataclass, field
from typing import Callable, Dict, Iterator, List, Optional, Sequence, Tuple

from .annotation import classify, parse_attributes
from .coords import human_bytes, human_size
from .errors import CancelledError, FileFormatError, GGEError

ProgressCB = Optional[Callable[[int, str], None]]
CancelCheck = Optional[Callable[[], bool]]

ROLE_ORDER = ("genome", "annotation", "mrna", "cds", "protein")
REQUIRED_ROLES = ("genome", "annotation")

ROLE_LABELS = {
    "genome": "基因组序列",
    "annotation": "基因注释",
    "mrna": "参考 mRNA",
    "cds": "参考 CDS",
    "protein": "参考蛋白",
    "": "未识别",
}

ARCHIVE_SUFFIXES = (".zip", ".tar", ".tar.gz", ".tgz", ".tar.bz2", ".tbz2",
                    ".tar.xz", ".txz")

#: 这些扩展名明显不是我们要的文本数据，直接跳过，不浪费解压时间
_SKIP_SUFFIXES = (
    ".bam", ".bai", ".cram", ".crai", ".vcf", ".vcf.gz", ".bcf", ".bed", ".bed.gz",
    ".bw", ".bigwig", ".tbi", ".csi", ".fai", ".gzi", ".png", ".jpg", ".jpeg", ".gif",
    ".pdf", ".xlsx", ".xls", ".docx", ".doc", ".pptx", ".db", ".sqlite", ".sqlite3",
    ".py", ".pyc", ".r", ".sh", ".bat", ".md", ".html", ".htm", ".json", ".xml",
    ".idx", ".dict", ".zip", ".tar", ".gz.tmp", ".log", ".csv", ".tsv",
)
_SKIP_NAMES = {".ds_store", "thumbs.db", "desktop.ini", ".gitignore", "readme", "readme.txt"}
_SKIP_DIR_PARTS = {"__macosx", ".git", ".svn", ".idea", "__pycache__"}

_HEAD_BYTES = 256 * 1024          # 每个成员嗅探用多少字节
_FASTA_SAMPLE_BYTES = 32 * 1024 * 1024   # FASTA 统计最多读多少字节
_FASTA_SAMPLE_RECORDS = 5000             # 以及最多读多少条记录
_MAX_TOKENS = 2000
_MAX_SEQIDS = 5000

_NT_LETTERS = set("ACGTUNacgtun")
_PROTEIN_HINT = set("EFILPQXZJOefilpqxzjo")

#: 文件名关键字（只作加分项）
_NAME_HINTS: Dict[str, Tuple[str, ...]] = {
    "genome": ("genome", "chromosome", "assembly", "genomic", "dna", "contig", "scaffold",
               "基因组", "染色体", "组装", "参考序列"),
    "mrna": ("mrna", "transcript", "transcripts", "cdna", "rna", "est", "转录本", "转录组"),
    "cds": ("cds", "coding", "orf", "编码区", "编码序列"),
    "protein": ("protein", "proteins", "pep", "peptide", "prot", "faa", "氨基酸", "蛋白"),
}

_NAME_WEIGHT = 30.0


def _residues(n: int, protein: bool = False) -> str:
    """碱基/氨基酸数的易读写法。蛋白不用 bp，避免"120 bp 的蛋白"这种误导。"""
    if not protein:
        return human_size(n)
    if n >= 1_000_000:
        return f"{n / 1e6:.2f} M aa"
    if n >= 1_000:
        return f"{n / 1e3:.1f} k aa"
    return f"{n:,} aa"


# --------------------------------------------------------------------------- #
# 数据结构
# --------------------------------------------------------------------------- #
@dataclass
class MemberInfo:
    name: str
    size: int
    packed: int = 0

    @property
    def base(self) -> str:
        return os.path.basename(self.name)


@dataclass
class FastaStats:
    n_records: int = 0
    total_bases: int = 0
    tokens: List[str] = field(default_factory=list)
    nt_frac: float = 0.0
    lengths: List[int] = field(default_factory=list)
    complete: bool = False          # 是否读到了文件末尾（否则只是采样）

    @property
    def median_length(self) -> int:
        if not self.lengths:
            return 0
        s = sorted(self.lengths)
        return s[len(s) // 2]

    def divisible_by_3(self) -> float:
        if not self.lengths:
            return 0.0
        return sum(1 for n in self.lengths if n % 3 == 0) / len(self.lengths)

    def text(self, protein: bool = False) -> str:
        if self.n_records <= 0:
            return "没有序列"
        rec = f"{self.n_records:,} 条" + ("" if self.complete else "以上")
        approx = "" if self.complete else "（采样，长度约为）"
        return f"{rec}序列，{approx}{_residues(self.total_bases, protein)}"


@dataclass
class AnnotationSummary:
    """流式统计出的注释概要，只保留判断需要用到的量。"""

    fmt: str = ""
    n_lines: int = 0
    n_genes: int = 0
    n_transcripts: int = 0
    n_exons: int = 0
    n_cds: int = 0
    seqids: List[str] = field(default_factory=list)
    n_seqids_hint: int = 0          # 清单里记的序列条数（没解析注释时用它显示）
    tx_ids: List[str] = field(default_factory=list)
    tx_exon_len: Dict[str, int] = field(default_factory=dict)
    tx_cds_len: Dict[str, int] = field(default_factory=dict)
    complete: bool = False

    def seqid_set(self) -> set:
        return {s for s in self.seqids}

    def tx_id_set(self) -> set:
        return set(self.tx_ids)

    def text(self) -> str:
        return (f"{self.fmt.upper()} · 基因 {self.n_genes:,} · 转录本 {self.n_transcripts:,} · "
                f"外显子 {self.n_exons:,} · CDS {self.n_cds:,} · "
                f"序列 {len(self.seqids) or self.n_seqids_hint:,}")


@dataclass
class FileCandidate:
    member: str
    size: int
    fmt: str = ""                 # FASTA / GFF3 / GTF
    role: str = ""
    score: float = 0.0
    confidence: float = 0.0
    ambiguous: bool = False
    evidence: List[str] = field(default_factory=list)
    stats: Optional[FastaStats] = None
    ann: Optional[AnnotationSummary] = None

    @property
    def base(self) -> str:
        return os.path.basename(self.member)

    def stat_text(self) -> str:
        if self.fmt == "FASTA" and self.stats is not None:
            protein = (self.role == "protein"
                       or (self.stats.nt_frac and self.stats.nt_frac < 0.80))
            return self.stats.text(protein=protein)
        if self.ann is not None:
            return self.ann.text()
        return human_bytes(self.size)

    def evidence_text(self) -> str:
        return "；".join(self.evidence) if self.evidence else "—"


@dataclass
class BundleScan:
    path: str = ""
    kind: str = ""
    archive_size: int = 0
    total_uncompressed: int = 0
    candidates: List[FileCandidate] = field(default_factory=list)
    ignored: List[Tuple[str, str]] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    annotation: Optional[AnnotationSummary] = None
    member_count: int = 0
    manifest: Optional[Dict[str, object]] = None      # 包内数据集清单（若存在）

    def by_role(self) -> Dict[str, Optional[FileCandidate]]:
        out: Dict[str, Optional[FileCandidate]] = {r: None for r in ROLE_ORDER}
        for c in self.candidates:
            if c.role in out and out[c.role] is None:
                out[c.role] = c
        return out

    def missing(self) -> List[str]:
        got = self.by_role()
        return [r for r in REQUIRED_ROLES if got.get(r) is None]

    def recognized_count(self) -> int:
        return len({c.role for c in self.candidates if c.role})

    def name_hint(self) -> str:
        """默认数据集名：清单里的名字优先，其次压缩包文件名。"""
        if self.manifest:
            name = str(self.manifest.get("name") or "").strip()
            if name:
                return name
        base = os.path.basename(self.path or "")
        for ext in (".tar.gz", ".tar.bz2", ".tar.xz"):
            if base.lower().endswith(ext):
                return base[: -len(ext)]
        return os.path.splitext(base)[0]

    def summary_text(self) -> str:
        got = self.by_role()
        parts = [f"{ROLE_LABELS[r]}={'✓' if got.get(r) else '×'}" for r in ROLE_ORDER]
        return " · ".join(parts)

    def text_report(self) -> str:
        lines = [f"压缩包：{os.path.basename(self.path)}（{self.kind}，"
                 f"{human_bytes(self.archive_size)}，{self.member_count} 个文件）"]
        for c in self.candidates:
            lines.append(f"  [{ROLE_LABELS.get(c.role, '?')}] {c.member}  "
                         f"({c.fmt}, {c.stat_text()})  {c.evidence_text()}")
        for name, why in self.ignored:
            lines.append(f"  [跳过] {name}  {why}")
        for w in self.warnings:
            lines.append(f"  ! {w}")
        return "\n".join(lines)


@dataclass
class BundleImport:
    """解包结果：role -> 磁盘路径。"""

    dest: str = ""
    paths: Dict[str, str] = field(default_factory=dict)
    members: Dict[str, str] = field(default_factory=dict)   # role -> 包内原始路径
    warnings: List[str] = field(default_factory=list)

    def missing(self) -> List[str]:
        return [r for r in REQUIRED_ROLES if not self.paths.get(r)]

    def text(self) -> str:
        rows = [f"{ROLE_LABELS[r]}: {os.path.basename(self.paths[r])}"
                for r in ROLE_ORDER if self.paths.get(r)]
        return "；".join(rows)


# --------------------------------------------------------------------------- #
# 压缩包读写
# --------------------------------------------------------------------------- #
def is_archive(path: str) -> bool:
    """按扩展名判断是否是本模块支持的压缩包。"""
    low = (path or "").lower()
    return any(low.endswith(s) for s in ARCHIVE_SUFFIXES)


def archive_kind(path: str) -> str:
    """返回 'zip' / 'tar' / ''（内容优先，扩展名兜底）。"""
    if not os.path.isfile(path):
        return ""
    if zipfile.is_zipfile(path):
        return "zip"
    try:
        if tarfile.is_tarfile(path):
            return "tar"
    except (OSError, tarfile.TarError):
        pass
    return ""


def _open_archive(path: str):
    kind = archive_kind(path)
    if kind == "zip":
        return kind, zipfile.ZipFile(path, "r", allowZip64=True)
    if kind == "tar":
        return kind, tarfile.open(path, "r:*")
    raise FileFormatError(
        f"不是可识别的压缩包：{os.path.basename(path)}",
        "支持 .zip / .tar / .tar.gz / .tgz / .tar.bz2 / .tar.xz；"
        "如果文件是 .7z 或 .rar，请先解压再重新打包成 zip。",
    )


def _list_members(kind: str, archive) -> Tuple[List[MemberInfo], List[Tuple[str, str]], int]:
    members: List[MemberInfo] = []
    ignored: List[Tuple[str, str]] = []
    total = 0
    if kind == "zip":
        infos = archive.infolist()
    else:
        infos = archive.getmembers()
    for info in infos:
        name = info.filename if kind == "zip" else info.name
        name = name.replace("\\", "/")
        is_dir = info.is_dir() if kind == "zip" else info.isdir()
        if is_dir:
            continue
        parts = [p for p in name.split("/") if p and p != "."]
        if not parts:
            continue
        low_parts = [p.lower() for p in parts]
        if any(p in _SKIP_DIR_PARTS for p in low_parts[:-1]):
            ignored.append((name, "系统/版本控制目录"))
            continue
        if low_parts[-1] in _SKIP_NAMES:
            ignored.append((name, "系统文件或说明文件"))
            continue
        if low_parts[-1].startswith("."):
            ignored.append((name, "隐藏文件"))
            continue
        # zip 用 file_size / compress_size，tar 用 size —— 两者都要兼容
        size = int(getattr(info, "file_size", getattr(info, "size", 0)) or 0)
        packed = int(getattr(info, "compress_size", getattr(info, "size", 0)) or 0)
        members.append(MemberInfo(name=name, size=size, packed=packed))
        total += max(size, 0)
    return members, ignored, total


def _inner_suffix(name: str) -> str:
    low = name.lower()
    for ext in (".gz", ".bz2", ".xz"):
        if low.endswith(ext):
            return low[: -len(ext)]
    return low


def _looks_unwanted(name: str) -> Optional[str]:
    inner = _inner_suffix(name)
    for ext in _SKIP_SUFFIXES:
        if inner.endswith(ext):
            return f"{ext} 不是本工具需要的序列/注释格式"
    return None


def _raw_stream(kind: str, archive, member: MemberInfo):
    if kind == "zip":
        return archive.open(member.name, "r")
    fh = archive.extractfile(member.name)
    if fh is None:
        raise FileFormatError(f"无法读取压缩包成员：{member.name}", "压缩包可能已损坏。")
    return fh


def _text_stream(kind: str, archive, member: MemberInfo):
    """打开成员的**解压后文本流**（成员自身是 .gz/.bz2/.xz 时再套一层）。"""
    raw = _raw_stream(kind, archive, member)
    low = member.name.lower()
    if low.endswith(".gz"):
        stream = gzip.GzipFile(fileobj=raw)
    elif low.endswith(".bz2"):
        stream = bz2.BZ2File(raw)
    elif low.endswith(".xz"):
        stream = lzma.LZMAFile(raw)
    else:
        stream = raw
    return stream


class _MemberText:
    """上下文管理器：逐行读取成员文本，退出时按顺序关闭各层句柄。"""

    def __init__(self, kind: str, archive, member: MemberInfo) -> None:
        self._stream = _text_stream(kind, archive, member)
        self._text = io.TextIOWrapper(self._stream, encoding="utf-8",
                                     errors="replace", newline=None)

    def __enter__(self) -> "_MemberText":
        return self

    def __exit__(self, *exc) -> None:
        try:
            self._text.close()
        except Exception:
            pass
        try:
            self._stream.close()
        except Exception:
            pass

    def lines(self) -> Iterator[str]:
        for line in self._text:
            yield line

    def read(self, n: int) -> str:
        return self._text.read(n)


def _read_head(kind: str, archive, member: MemberInfo, limit: int = _HEAD_BYTES) -> bytes:
    raw = _raw_stream(kind, archive, member)
    try:
        low = member.name.lower()
        if low.endswith(".gz"):
            stream = gzip.GzipFile(fileobj=raw)
        elif low.endswith(".bz2"):
            stream = bz2.BZ2File(raw)
        elif low.endswith(".xz"):
            stream = lzma.LZMAFile(raw)
        else:
            stream = raw
        try:
            return stream.read(limit)
        finally:
            if stream is not raw:
                stream.close()
    finally:
        raw.close()


def _is_binary(data: bytes) -> bool:
    return b"\x00" in data


def _is_bgzf(raw: bytes) -> bool:
    """BGZF 是带 ``BC`` 子字段的 gzip，pyfaidx 能直接随机访问它。"""
    if len(raw) < 18 or raw[:3] != b"\x1f\x8b\x08":
        return False
    if not (raw[3] & 0x04):        # FEXTRA
        return False
    xlen = int.from_bytes(raw[10:12], "little")
    extra = raw[12:12 + xlen]
    return b"BC" in extra[:128]


def _member_is_bgzf(kind: str, archive, member: MemberInfo) -> bool:
    if not member.name.lower().endswith(".gz"):
        return False
    raw = _raw_stream(kind, archive, member)
    try:
        return _is_bgzf(raw.read(4096))
    finally:
        raw.close()


# --------------------------------------------------------------------------- #
# 内容嗅探
# --------------------------------------------------------------------------- #
_GFF3_ATTR_RE = re.compile(r"\b(ID|Parent)\s*=")
_GTF_ATTR_RE = re.compile(r'[A-Za-z0-9_]+\s+"[^"]*"')


def _sniff_head(text: str) -> Tuple[str, int, int]:
    """返回 (格式, GFF3 特征数, GTF 特征数)。格式取 'FASTA' / 'GFF3' / 'GTF' / ''。"""
    n_fasta = 0
    gff3 = 0
    gtf = 0
    n_cols = 0
    for i, line in enumerate(text.splitlines()):
        if i > 4000:
            break
        s = line.strip()
        if not s:
            continue
        if s.startswith(">"):
            n_fasta += 1
            if gff3 or gtf or n_cols:
                break
            continue
        if s.startswith("#"):
            continue
        if "\t" in s:
            cols = s.split("\t")
            if len(cols) >= 8:
                n_cols += 1
                attrs = cols[8] if len(cols) >= 9 else ""
                if _GFF3_ATTR_RE.search(attrs) and not _GTF_ATTR_RE.search(attrs):
                    gff3 += 1
                elif _GTF_ATTR_RE.search(attrs):
                    gtf += 1
                continue
        if n_fasta:
            continue
        # 既不是 FASTA 也不是 9 列 Tab：可能是分隔符不对的注释，放弃
        return "", gff3, gtf

    if n_fasta and not (gff3 or gtf or n_cols):
        return "FASTA", gff3, gtf
    if n_cols:
        if gtf > gff3:
            return "GTF", gff3, gtf
        return "GFF3", gff3, gtf
    if n_fasta:
        return "FASTA", gff3, gtf
    return "", gff3, gtf


# --------------------------------------------------------------------------- #
# FASTA 统计
# --------------------------------------------------------------------------- #
def _fasta_stats(kind: str, archive, member: MemberInfo,
                 max_bytes: int = 0, max_records: int = 0) -> FastaStats:
    """采样统计 FASTA：记录数、总长、标题 token、字母构成、各记录长度。

    大文件只读前 ``max_bytes`` 字节 / ``max_records`` 条，此时
    ``FastaStats.complete`` 为 False，界面上会标注"(采样)"。
    """
    max_bytes = max_bytes or _FASTA_SAMPLE_BYTES
    max_records = max_records or _FASTA_SAMPLE_RECORDS
    stats = FastaStats()
    letters: Counter = Counter()
    consumed = 0
    name: Optional[str] = None
    cur_len = 0
    stopped = False
    with _MemberText(kind, archive, member) as mtext:
        for line in mtext.lines():
            consumed += len(line)
            if line.startswith(">"):
                if name is not None:
                    stats.lengths.append(cur_len)
                    stats.total_bases += cur_len
                if len(stats.tokens) < _MAX_TOKENS:
                    token = line[1:].strip().split()[0] if line[1:].strip() else ""
                    stats.tokens.append(token)
                stats.n_records += 1
                name = line[1:].strip()
                cur_len = 0
                if stats.n_records >= max_records or consumed >= max_bytes:
                    stopped = True
                    break
                continue
            seq = line.strip()
            if not seq:
                continue
            cur_len += len(seq)
            letters.update(seq)
    if not stopped and name is not None:      # 正常读完，收尾最后一条
        stats.lengths.append(cur_len)
        stats.total_bases += cur_len
        stats.complete = True
    if stats.n_records == 0:
        return stats
    total_letters = sum(letters.values())
    if total_letters:
        nt = sum(letters[ch] for ch in _NT_LETTERS)
        stats.nt_frac = nt / total_letters
    if not stats.complete and stats.total_bases and consumed:
        # 采样：按"已读字节占比"估算总长，只在界面上标注为约数
        stats.total_bases = int(stats.total_bases / consumed * max(member.size, consumed))
    return stats


# --------------------------------------------------------------------------- #
# 注释概要
# --------------------------------------------------------------------------- #
def _parent_of(attrs: Dict[str, str]) -> str:
    parent = attrs.get("parent") or attrs.get("transcript_id") or ""
    return parent.split(",")[0].strip() if parent else ""


def _summarize_lines(stream, fmt: str, progress_cb: ProgressCB = None,
                     cancel_check: CancelCheck = None,
                     max_lines: int = 0) -> AnnotationSummary:
    """注释流式统计的共用实现（``summarize_annotation*`` 都走这里）。"""
    ann = AnnotationSummary(fmt=(fmt or "GFF3").lower())
    gene_ids: set = set()
    tx_ids: set = set()
    seqids: List[str] = []
    seqid_seen: set = set()
    n_line = 0
    for line in stream:
        n_line += 1
        if max_lines and n_line > max_lines:
            break
        if progress_cb and n_line % 200_000 == 0:
            progress_cb(0, f"扫描注释 {n_line:,} 行 …")
            if cancel_check and cancel_check():
                raise CancelledError()
        if line.startswith("##FASTA"):
            break
        if not line.strip() or line.startswith("#"):
            continue
        cols = line.rstrip("\r\n").split("\t")
        if len(cols) < 8:
            continue
        if len(cols) == 8:
            cols.append(".")
        seqid = cols[0].strip()
        if seqid not in seqid_seen and len(seqids) < _MAX_SEQIDS:
            seqid_seen.add(seqid)
            seqids.append(seqid)
        ckind = classify(cols[2].strip())
        if ckind == "ignored":
            continue
        try:
            start = int(cols[3])
            end = int(cols[4])
        except ValueError:
            continue
        length = abs(end - start) + 1
        attrs = parse_attributes(cols[8])
        # 基因 ID 有三个来源：gene 行、transcript 行的 Parent、GTF 的 gene_id
        gid = attrs.get("gene_id") or ""
        if gid:
            gene_ids.add(gid)
        parent = _parent_of(attrs)
        if ckind == "gene":
            gene_ids.add(attrs.get("id") or attrs.get("gene_id")
                         or attrs.get("name") or f"#{len(gene_ids)}")
            continue
        if ckind == "transcript":
            tid = attrs.get("id") or attrs.get("transcript_id") or attrs.get("name") or ""
            if tid:
                tx_ids.add(tid)
            if attrs.get("parent"):          # GFF3：Parent 就是基因
                gene_ids.add(attrs["parent"].split(",")[0].strip())
            continue
        if not parent:
            continue
        if ckind == "exon":
            ann.n_exons += 1
            if parent not in ann.tx_exon_len:
                tx_ids.add(parent)
            ann.tx_exon_len[parent] = ann.tx_exon_len.get(parent, 0) + length
        elif ckind == "cds":
            ann.n_cds += 1
            if parent not in ann.tx_cds_len:
                tx_ids.add(parent)
            ann.tx_cds_len[parent] = ann.tx_cds_len.get(parent, 0) + length
    ann.n_lines = n_line
    ann.n_genes = len(gene_ids)
    ann.tx_ids = list(tx_ids)
    ann.n_transcripts = len(tx_ids) or max(len(ann.tx_exon_len), len(ann.tx_cds_len))
    ann.seqids = seqids
    ann.complete = not max_lines
    return ann


def summarize_annotation(kind: str, archive, member: MemberInfo,
                         fmt: str = "GFF3", progress_cb: ProgressCB = None,
                         cancel_check: CancelCheck = None,
                         max_lines: int = 0) -> AnnotationSummary:
    """流式统计压缩包内注释文件：基因/转录本/外显子/CDS 数量与每个转录本的长度。

    不构建完整模型，因此几百万行的注释也只占几十 MB 内存。
    """
    with _MemberText(kind, archive, member) as mtext:
        return _summarize_lines(mtext.lines(), fmt, progress_cb, cancel_check, max_lines)


def summarize_annotation_path(path: str, fmt: str = "", progress_cb: ProgressCB = None,
                              cancel_check: CancelCheck = None,
                              max_lines: int = 0) -> AnnotationSummary:
    """对磁盘上的注释文件做同样的统计（供测试、CLI 与导入后复核使用）。"""
    opener = gzip.open if path.lower().endswith(".gz") else open
    guess = fmt or ("gtf" if path.lower().endswith((".gtf", ".gtf.gz")) else "gff3")
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        return _summarize_lines(fh, guess, progress_cb, cancel_check, max_lines)


def _name_hits(member_name: str) -> Dict[str, List[str]]:
    low = member_name.lower()
    hits: Dict[str, List[str]] = {}
    for role, words in _NAME_HINTS.items():
        for w in words:
            if w.isascii():
                if re.search(r"(^|[^a-z])" + re.escape(w) + r"([^a-z]|$)", low):
                    hits.setdefault(role, []).append(w)
            elif w in member_name:
                hits.setdefault(role, []).append(w)
    return hits


def _lookup_token(token: str, table: Dict[str, int]) -> Optional[str]:
    """把参考序列的标题 token 映射回注释里的转录本 ID，容忍常见后缀差异。"""
    if not token:
        return None
    cands = [token, token.split("|")[0], token.split(":")[0]]
    for suffix in ("_cds", ".cds", "-cds", "_CDS", "_mrna", ".mrna", "-mrna", "_transcript",
                   "_protein", "_prot", "_pep", "_aa", ".protein", ".pep"):
        if token.endswith(suffix):
            cands.append(token[: -len(suffix)])
    for c in cands:
        if c in table:
            return c
    low = {k.lower(): k for k in table}
    for c in cands:
        if c.lower() in low:
            return low[c.lower()]
    return None


def _length_fit(pairs: Sequence[Tuple[str, int]], table: Dict[str, int],
                tol: float = 0.02, allow_third: bool = False) -> Tuple[float, int]:
    """返回 (吻合比例, 参与比较的记录数)。"""
    if not table or not pairs:
        return 0.0, 0
    hit = 0
    total = 0
    for token, length in pairs:
        key = _lookup_token(token, table)
        if key is None:
            continue
        expect = table[key]
        tol_abs = max(3.0, tol * max(expect, 1))
        ok = abs(expect - length) <= tol_abs
        if not ok and allow_third:
            # 蛋白长度 ≈ CDS/3（去掉终止密码子后再差 1）
            for cand in (max(1, expect // 3 - 1), max(1, expect // 3)):
                if abs(cand - length) <= max(2.0, 0.03 * cand):
                    ok = True
                    break
        total += 1
        if ok:
            hit += 1
    # 至少 2 条比对得上才算证据：参考序列文件往往只有几条（小基因组/测试集）
    if total < 2:
        return 0.0, total
    return hit / total, total


def _score_candidate(cand: FileCandidate, ann: Optional[AnnotationSummary]) -> None:
    """为候选文件给每个角色打分，结果写回 cand.role/score/evidence。"""
    scores = {r: 0.0 for r in ROLE_ORDER}
    ev = {r: [] for r in ROLE_ORDER}

    stats = cand.stats
    name_hits = _name_hits(cand.member)

    # ---- 1. 字母表：蛋白/核酸 ----
    is_protein = False
    if stats is not None and stats.n_records and stats.nt_frac:
        if stats.nt_frac < 0.80:
            is_protein = True
            scores["protein"] += 70
            ev["protein"].append(f"字母构成像蛋白（ACGTUN 仅占 {stats.nt_frac * 100:.0f}%）")
        else:
            ev["genome"].append(f"核酸字母占 {stats.nt_frac * 100:.1f}%")

    # ---- 2. 文件名（只加分，不决定） ----
    for role, words in name_hits.items():
        if is_protein and role != "protein":
            continue
        scores[role] += _NAME_WEIGHT
        ev[role].append("文件名含 " + "/".join(words[:3]))
    if is_protein:
        scores["genome"] = scores["mrna"] = scores["cds"] = 0.0

    # ---- 3. 与注释的对应关系 ----
    if stats is not None and ann is not None:
        seqids = ann.seqid_set()
        tx_ids = ann.tx_id_set()
        toks = [t for t in stats.tokens if t]
        seq_hit = [t for t in toks if t in seqids]
        tx_hit = [t for t in toks if t in tx_ids]
        if seq_hit:
            scores["genome"] += 55
            ev["genome"].append(f"序列名与注释一致（{seq_hit[0]}"
                                + (f" 等 {len(seq_hit)} 条" if len(seq_hit) > 1 else "") + "）")
        if tx_hit:
            scores["mrna"] += 26
            scores["cds"] += 26
            ev["mrna"].append(f"序列名命中转录本 ID（{tx_hit[0]} 等 {len(tx_hit)} 条）")
            ev["cds"].append(f"序列名命中转录本 ID（{tx_hit[0]} 等 {len(tx_hit)} 条）")

        n_rec = abs(stats.n_records)
        if ann.n_transcripts and n_rec:
            ratio = n_rec / ann.n_transcripts
            if 0.7 <= ratio <= 1.4:
                scores["mrna"] += 20
                scores["cds"] += 20
                ev["mrna"].append(f"记录数 {n_rec:,} ≈ 转录本数 {ann.n_transcripts:,}")
                ev["cds"].append(f"记录数 {n_rec:,} ≈ 转录本数 {ann.n_transcripts:,}")
            elif ratio < 0.15 and n_rec <= 1000:
                scores["genome"] += 22
                ev["genome"].append(f"记录数 {n_rec:,} 远少于转录本数 {ann.n_transcripts:,}")

        pairs = list(zip(stats.tokens, stats.lengths))
        if pairs:
            fit_mrna, n_m = _length_fit(pairs, ann.tx_exon_len)
            fit_cds, n_c = _length_fit(pairs, ann.tx_cds_len)
            if fit_mrna >= 0.5:
                scores["mrna"] += 45 * fit_mrna
                ev["mrna"].append(f"{n_m} 条记录长度与注释 mRNA 长度吻合（{fit_mrna * 100:.0f}%）")
            if fit_cds >= 0.5:
                scores["cds"] += 45 * fit_cds
                ev["cds"].append(f"{n_c} 条记录长度与注释 CDS 长度吻合（{fit_cds * 100:.0f}%）")
            if is_protein:
                fit_p, n_p = _length_fit(pairs, ann.tx_cds_len, allow_third=True)
                if fit_p >= 0.5:
                    scores["protein"] += 45 * fit_p
                    ev["protein"].append(f"{n_p} 条记录长度 ≈ 注释 CDS/3（{fit_p * 100:.0f}%）")

    # ---- 4. 长度能否被 3 整除 ----
    if stats is not None and stats.lengths and not is_protein:
        frac3 = stats.divisible_by_3()
        if frac3 >= 0.9 and len(stats.lengths) >= 3:
            scores["cds"] += 12
            ev["cds"].append(f"{frac3 * 100:.0f}% 的记录长度能被 3 整除")
        elif frac3 <= 0.4 and len(stats.lengths) >= 3:
            scores["mrna"] += 8
            ev["mrna"].append(f"仅 {frac3 * 100:.0f}% 的记录长度能被 3 整除")

    # ---- 5. 没有任何注释时的兜底 ----
    if ann is None and stats is not None and not is_protein:
        n_rec = abs(stats.n_records)
        if cand.size >= 20 * 1024 * 1024 and n_rec <= 1000:
            scores["genome"] += 18
            ev["genome"].append("文件很大而记录数很少，像基因组")

    ranked = sorted(((r, s) for r, s in scores.items() if s > 0), key=lambda x: -x[1])
    if not ranked or ranked[0][1] < 25:
        cand.role = ""
        cand.score = ranked[0][1] if ranked else 0.0
        cand.confidence = 0.0
        cand.evidence = ["证据不足，无法自动判定，请手动指定角色"]
        return
    best_role, best_score = ranked[0]
    second = ranked[1][1] if len(ranked) > 1 else 0.0
    cand.role = best_role
    cand.score = best_score
    cand.confidence = max(0.2, min(0.99, best_score / 100.0))
    cand.ambiguous = (best_score - second) < 12 and best_score < 100
    cand.evidence = ev.get(best_role, []) or ["综合判断"]
    if cand.ambiguous and len(ranked) > 1:
        cand.evidence.append(f"与「{ROLE_LABELS.get(ranked[1][0], ranked[1][0])}」得分接近，建议确认")


def _assign_unique(cands: List[FileCandidate]) -> None:
    """同一角色只能留一个文件，其余标为未识别并说明原因。"""
    best: Dict[str, FileCandidate] = {}
    for c in cands:
        if not c.role:
            continue
        cur = best.get(c.role)
        if cur is None or c.score > cur.score:
            if cur is not None:
                cur.role = ""
                cur.evidence.append(f"与 {os.path.basename(c.member)} 同为"
                                    f"{ROLE_LABELS[c.role]}候选且分数更低，已置为未识别")
            best[c.role] = c
        else:
            role = c.role
            c.role = ""
            c.evidence.append(f"与 {os.path.basename(best[role].member)} 同为"
                              f"{ROLE_LABELS[role]}候选且分数更低，已置为未识别")


def _demoted_warnings(scan: BundleScan) -> None:
    """把"多个文件抢同一个角色"明确说出来，而不是只留一句判定依据。"""
    for c in scan.candidates:
        text = c.evidence_text()
        if not c.role and "同为" in text and "候选" in text:
            scan.warnings.append(
                f"{os.path.basename(c.member)} 也像某种数据但没有被采用（{text}）；"
                "如果判定不对，请在列表中手动改")


def _fallback_single_nucleotide(scan: BundleScan) -> None:
    """包内只有一个核酸 FASTA 又没有基因组时，先当基因组用 —— 但要说清是猜的。

    典型场景是"只给了基因组、没有注释"的包：没有注释做交叉验证，任何判定都是猜。
    与其把唯一候选晾在那里，不如标成低置信度让用户确认。
    """
    if scan.by_role().get("genome") is not None:
        return
    pool = [c for c in scan.candidates
            if c.fmt == "FASTA" and c.stats is not None
            and c.stats.n_records and c.stats.nt_frac >= 0.80 and not c.role]
    if len(pool) != 1:
        return
    c = pool[0]
    c.role = "genome"
    c.score = max(c.score, 20.0)
    c.confidence = 0.3
    c.ambiguous = True
    c.evidence.append("包内唯一的核酸序列文件，暂按基因组处理")
    scan.warnings.append(
        "没有注释可做交叉验证，只能凭「包内唯一的核酸 FASTA」把 "
        f"{os.path.basename(c.member)} 当作基因组，请确认")


# --------------------------------------------------------------------------- #
# 扫描
# --------------------------------------------------------------------------- #
def scan_archive(path: str, progress_cb: ProgressCB = None,
                 cancel_check: CancelCheck = None) -> BundleScan:
    """扫描压缩包并给出每个文件的角色判定。

    只做两件事：读每个成员的头部（判断格式），以及对候选的注释/FASTA 做**有上限**的
    流式统计。不会把整个包解压到磁盘，因此 3 GB 的包也能在几秒内给出结论。
    """
    if not os.path.isfile(path):
        raise FileFormatError(f"文件不存在：{path}", "确认压缩包路径是否正确。")
    if os.path.getsize(path) == 0:
        raise FileFormatError(f"压缩包大小为 0：{os.path.basename(path)}", "重新获取该压缩包。")

    kind, archive = _open_archive(path)
    scan = BundleScan(path=path, kind=kind, archive_size=os.path.getsize(path))
    try:
        members, ignored, total = _list_members(kind, archive)
        scan.ignored = list(ignored)
        scan.total_uncompressed = total
        scan.member_count = len(members)
        if not members:
            raise FileFormatError(
                f"压缩包里没有任何文件：{os.path.basename(path)}",
                "确认压缩包不是空包，且文件都在子目录里也能被读到（本工具会递归查找）。",
            )
        if total > 0 and scan.archive_size > 0 and total / scan.archive_size > 300:
            scan.warnings.append(
                f"压缩比异常（{total / scan.archive_size:.0f}:1），请确认文件来源可信")

        # ---- 第一轮：格式嗅探 ----
        n = len(members)
        for i, m in enumerate(members):
            if cancel_check and cancel_check():
                raise CancelledError()
            if progress_cb:
                progress_cb(int(i / n * 55), f"检查 {m.base}（{i + 1}/{n}）…")
            if m.base.lower() == MANIFEST_NAME:
                continue                     # 数据集清单：不参与嗅探，稍后单独读
            why = _looks_unwanted(m.name)
            if why:
                scan.ignored.append((m.name, why))
                continue
            if m.size == 0:
                scan.ignored.append((m.name, "空文件"))
                continue
            try:
                head = _read_head(kind, archive, m)
            except (OSError, EOFError, zipfile.BadZipFile, tarfile.TarError) as exc:
                scan.ignored.append((m.name, f"读取失败：{exc}"))
                continue
            if not head:
                scan.ignored.append((m.name, "空文件"))
                continue
            if _is_binary(head):
                scan.ignored.append((m.name, "二进制文件"))
                continue
            fmt, n3, ng = _sniff_head(head.decode("utf-8", "replace"))
            if not fmt:
                scan.ignored.append((m.name, "内容既不是 FASTA 也不是 GFF/GTF"))
                continue
            cand = FileCandidate(member=m.name, size=m.size, fmt=fmt)
            scan.candidates.append(cand)

        if not scan.candidates:
            raise FileFormatError(
                f"{os.path.basename(path)} 里没有找到 FASTA 或 GFF/GTF 文件",
                "本工具需要至少一个基因组 FASTA 和一个注释 GFF/GTF。"
                "已跳过：" + ("、".join(x[0] for x in scan.ignored[:6]) or "（无）"),
            )

        # ---- 清单优先：本程序导出的包自带 gge_dataset.json ----
        manifest = _read_manifest_from_archive(kind, archive, members)
        scan.manifest = manifest
        mstats = _manifest_stats(manifest) if manifest else {}

        # ---- 第二轮：注释概要（优先用最大的那个注释）----
        ann_cands = [c for c in scan.candidates if c.fmt in ("GFF3", "GTF")]
        ann_cands.sort(key=lambda c: -c.size)
        if len(ann_cands) > 1:
            scan.warnings.append(
                f"包里有 {len(ann_cands)} 个注释文件（"
                + "、".join(os.path.basename(c.member) for c in ann_cands[:4])
                + "），默认用最大的那个做交叉验证，最终请在列表中确认")
        if mstats.get("n_genes"):
            # 清单里已经记录了统计量：不必再把 500 MB 的注释整份解析一遍
            ann_member = ""
            for entry in manifest.get("files") or []:
                if isinstance(entry, dict) and entry.get("role") == "annotation":
                    ann_member = str(entry.get("member") or "")
                    break
            fmt = "gtf" if ann_member.lower().endswith((".gtf", ".gtf.gz")) else "gff3"
            scan.annotation = AnnotationSummary(
                fmt=fmt,
                n_genes=mstats.get("n_genes", 0),
                n_transcripts=mstats.get("n_transcripts", 0),
                n_seqids_hint=mstats.get("n_seqids", 0),
                complete=True,
            )
            if progress_cb:
                progress_cb(58, "按包内清单读取注释统计（跳过全量解析）")
        elif ann_cands:
            primary = ann_cands[0]
            if progress_cb:
                progress_cb(58, f"统计注释 {os.path.basename(primary.member)} …")
            member = _member_of(members, primary.member)
            try:
                primary.ann = summarize_annotation(
                    kind, archive, member, primary.fmt, progress_cb, cancel_check)
                scan.annotation = primary.ann
            except Exception as exc:      # 注释坏了不应该让扫描整体失败
                scan.warnings.append(f"注释统计失败：{exc}")

        # ---- 第三轮：FASTA 统计 ----
        fasta_cands = [c for c in scan.candidates if c.fmt == "FASTA"]
        for i, c in enumerate(fasta_cands):
            if cancel_check and cancel_check():
                raise CancelledError()
            if progress_cb:
                progress_cb(60 + int(i / max(len(fasta_cands), 1) * 35),
                            f"统计 {os.path.basename(c.member)} …")
            member = _member_of(members, c.member)
            try:
                c.stats = _fasta_stats(kind, archive, member)
            except Exception as exc:
                scan.ignored.append((c.member, f"序列统计失败：{exc}"))
                c.stats = None
                continue
            if c.stats.n_records == 0:
                scan.ignored.append((c.member, "FASTA 中没有序列"))
                c.role = ""
                continue

        # ---- 判分 ----
        for c in scan.candidates:
            if c.fmt == "FASTA" and c.stats is None:
                c.role = ""
                c.evidence = ["无法读取序列内容"]
                continue
            if c.fmt in ("GFF3", "GTF"):
                c.role = "annotation"
                c.score = 200.0
                c.confidence = 0.99
                c.evidence = [f"内容为 {c.fmt} 注释（9 列 Tab + 属性列）"]
                continue
            _score_candidate(c, scan.annotation)
        _assign_unique(scan.candidates)
        _demoted_warnings(scan)
        _fallback_single_nucleotide(scan)
        if manifest:
            _apply_manifest(scan, manifest, members)

        got = scan.by_role()
        if got.get("genome") is None:
            scan.warnings.append("没有识别出基因组 FASTA —— 请在列表中手动指定")
        if got.get("annotation") is None:
            scan.warnings.append("没有识别出 GFF/GTF 注释 —— 请在列表中手动指定")
        for role in ("mrna", "cds", "protein"):
            if got.get(role) is None:
                scan.warnings.append(f"没有找到{ROLE_LABELS[role]}，该校验将不可用（不影响提取）")
        for c in scan.candidates:
            if c.role and c.ambiguous:
                scan.warnings.append(
                    f"{os.path.basename(c.member)} 判为{ROLE_LABELS[c.role]}但证据不够强，请确认")
    finally:
        try:
            archive.close()
        except Exception:
            pass
    if progress_cb:
        progress_cb(100, scan.summary_text())
    return scan


def _member_of(members: Sequence[MemberInfo], name: str) -> MemberInfo:
    for m in members:
        if m.name == name:
            return m
    raise FileFormatError(f"压缩包里找不到成员：{name}", "压缩包可能已损坏。")


# --------------------------------------------------------------------------- #
# 解包
# --------------------------------------------------------------------------- #
_WIN_BAD = re.compile(r'[<>:"/\\|?*\x00-\x1f]')
_RESERVED = {"con", "prn", "aux", "nul", *(f"com{i}" for i in range(1, 10)),
             *(f"lpt{i}" for i in range(1, 10))}


def safe_filename(name: str, fallback: str = "file") -> str:
    """把包内路径压成安全的扁平文件名（同时杜绝 zip-slip）。"""
    base = os.path.basename((name or "").replace("\\", "/"))
    base = _WIN_BAD.sub("_", base).strip().strip(".")
    if not base:
        return fallback
    stem, dot, ext = base.partition(".")
    if stem.lower() in _RESERVED:
        base = f"_{base}"
    return base[:180]


def _unique_name(name: str, used: set) -> str:
    if name not in used:
        used.add(name)
        return name
    stem, dot, ext = name.partition(".")
    i = 2
    while True:
        cand = f"{stem}_{i}{dot}{ext}"
        if cand not in used:
            used.add(cand)
            return cand
        i += 1


def _disk_free(path: str) -> int:
    try:
        return shutil.disk_usage(path).free
    except OSError:
        return 0


def import_bundle(scan: BundleScan, dest_dir: str,
                  roles: Optional[Dict[str, str]] = None,
                  progress_cb: ProgressCB = None,
                  cancel_check: CancelCheck = None) -> BundleImport:
    """按 ``roles``（role -> 包内路径）把文件解到 ``dest_dir``，返回实际路径。

    ``roles`` 为空时使用扫描得到的判定。只解压被选中的文件，其它成员原样留在包里。
    基因组若是普通 gzip（非 bgzf），会顺带解压成 .fa —— pyfaidx 只能对未压缩或
    bgzf 的 FASTA 建索引，普通 .gz 会在加载时报错。
    """
    mapping: Dict[str, str] = {}
    if roles is None:
        for role, cand in scan.by_role().items():
            if cand is not None:
                mapping[role] = cand.member
    else:
        for role, member in roles.items():
            if role in ROLE_ORDER and member:
                mapping[role] = member
    if not mapping:
        raise GGEError("没有选择要导入的文件", "请在列表里至少指定基因组与注释。")

    missing = [r for r in REQUIRED_ROLES if r not in mapping]
    if missing:
        raise GGEError(
            "缺少必填文件：" + "、".join(ROLE_LABELS[r] for r in missing),
            "基因组 FASTA 与 GFF/GTF 注释是必须的；参考序列可选。",
        )

    os.makedirs(dest_dir, exist_ok=True)
    free = _disk_free(dest_dir)
    need = int(scan.total_uncompressed * 0.35) if scan.total_uncompressed else 0
    # 头部空间按需给，不要固定 128 MB：演示包只有 10 KB，却要求 128 MB 空闲，
    # 在"盘快满了但还能写几 MB"的机器上会连一个小包都导不进来。
    headroom = max(16 * 1024 * 1024, min(128 * 1024 * 1024, need))
    if free and need and free < need + headroom:
        raise GGEError(
            f"磁盘空间可能不足：{os.path.dirname(dest_dir)} 所在磁盘可用 {human_bytes(free)}，"
            f"本次预计需要 {human_bytes(need)}（另需 {human_bytes(headroom)} 余量）",
            f"清理该磁盘后重试；或把数据集库目录设到其它分区"
            f"（设置环境变量 GGE_LIBRARY_DIR 指向空间充足的位置，重启程序生效）。",
        )

    kind, archive = _open_archive(scan.path)
    result = BundleImport(dest=dest_dir)
    total_bytes = 0
    plan: List[Tuple[str, MemberInfo, bool]] = []
    try:
        members, _ignored, _total = _list_members(kind, archive)
        for role in ROLE_ORDER:
            member_name = mapping.get(role)
            if not member_name:
                continue
            member = _member_of(members, member_name)
            decompress = False
            if role == "genome" and member.name.lower().endswith(".gz"):
                decompress = not _member_is_bgzf(kind, archive, member)
                if decompress:
                    result.warnings.append(
                        "基因组是普通 gzip（不是 bgzip），已解压为未压缩 FASTA，"
                        "否则无法建立 .fai 索引")
            plan.append((role, member, decompress))
            total_bytes += max(member.size, 1)

        used: set = set()
        done = 0
        for role, member, decompress in plan:
            if cancel_check and cancel_check():
                raise CancelledError()
            base = safe_filename(member.name, fallback=f"{role}.dat")
            if decompress:
                base = re.sub(r"\.gz$", "", base, flags=re.IGNORECASE) or f"{role}.fa"
            target = os.path.join(dest_dir, _unique_name(base, used))
            if progress_cb:
                progress_cb(int(done / max(total_bytes, 1) * 100),
                            f"解包 {base}（{ROLE_LABELS[role]}）…")
            written = _copy_member(kind, archive, member, target, decompress,
                                   total_bytes, done, progress_cb, cancel_check)
            done += max(member.size, 1)
            if written <= 0:
                raise GGEError(f"解包后文件为空：{base}",
                               "压缩包里的这个成员可能是空的。")
            result.paths[role] = target
            result.members[role] = member.name
        if progress_cb:
            progress_cb(100, "解包完成")
    finally:
        try:
            archive.close()
        except Exception:
            pass
    return result


def _copy_member(kind: str, archive, member: MemberInfo, target: str, decompress: bool,
                 total_bytes: int, done: int, progress_cb: ProgressCB,
                 cancel_check: CancelCheck) -> int:
    """把成员写到 ``target``，分块以便报告进度和响应取消。"""
    raw = _raw_stream(kind, archive, member)
    stream = raw
    low = member.name.lower()
    try:
        if low.endswith(".gz"):
            stream = gzip.GzipFile(fileobj=raw)
        elif low.endswith(".bz2"):
            stream = bz2.BZ2File(raw)
        elif low.endswith(".xz"):
            stream = lzma.LZMAFile(raw)
        written = 0
        with open(target, "wb") as out:
            while True:
                if cancel_check and cancel_check():
                    raise CancelledError()
                chunk = stream.read(1024 * 1024)
                if not chunk:
                    break
                out.write(chunk)
                written += len(chunk)
                if progress_cb:
                    pct = int((done + min(written, member.size or written))
                              / max(total_bytes, 1) * 100)
                    progress_cb(min(99, pct), f"解包 {member.base} …")
        return written
    finally:
        if stream is not raw:
            try:
                stream.close()
            except Exception:
                pass
        raw.close()


# --------------------------------------------------------------------------- #
# 便捷入口
# --------------------------------------------------------------------------- #
def import_archive(path: str, dest_dir: str,
                   roles: Optional[Dict[str, str]] = None,
                   progress_cb: ProgressCB = None,
                   cancel_check: CancelCheck = None) -> Tuple[BundleScan, BundleImport]:
    """扫描 + 解包一步到位（CLI 与测试用）。"""
    scan = scan_archive(path, progress_cb=progress_cb, cancel_check=cancel_check)
    result = import_bundle(scan, dest_dir, roles, progress_cb, cancel_check)
    return scan, result


def describe_roles(scan: BundleScan) -> str:
    """一行文字描述扫描结果，供状态栏使用。"""
    got = scan.by_role()
    items = [f"{ROLE_LABELS[r]}={os.path.basename(got[r].member)}" if got.get(r) else ""
             for r in ROLE_ORDER]
    return "，".join([x for x in items if x]) or "未识别出任何文件"


# --------------------------------------------------------------------------- #
# 数据集清单（导出 / 导入认清单）
# --------------------------------------------------------------------------- #
#: 清单文件名。放在包内任意位置都能被认出来，但导出时写在根目录。
MANIFEST_NAME = "gge_dataset.json"
MANIFEST_FORMAT = "gge-dataset"
MANIFEST_VERSION = 1

_MANIFEST_MAX_BYTES = 4 * 1024 * 1024


@dataclass
class BundleExport:
    """一次"打包成数据集"的结果。"""

    path: str = ""
    name: str = ""
    members: Dict[str, str] = field(default_factory=dict)   # role -> 包内路径
    total_bytes: int = 0                                    # 原始总大小
    archive_size: int = 0

    def ratio_text(self) -> str:
        if not self.total_bytes or not self.archive_size:
            return ""
        return f"压缩到 {self.archive_size / self.total_bytes * 100:.0f}%"

    def text(self) -> str:
        rows = [f"{ROLE_LABELS[r]}: {os.path.basename(self.members[r])}"
                for r in ROLE_ORDER if self.members.get(r)]
        return (f"{len(self.members)} 个文件（{human_bytes(self.total_bytes)} → "
                f"{human_bytes(self.archive_size)}）；" + "；".join(rows))


def build_manifest(name: str, paths: Dict[str, str],
                   note: str = "", source: str = "",
                   stats: Optional[Dict[str, int]] = None,
                   members: Optional[Dict[str, str]] = None,
                   created_at: str = "") -> Dict[str, object]:
    """生成 ``gge_dataset.json`` 的内容（纯数据，便于测试与复用）。"""
    from datetime import datetime

    from .version import APP_NAME, __version__

    files = []
    for role in ROLE_ORDER:
        src = paths.get(role, "")
        if not src:
            continue
        size = 0
        try:
            size = os.path.getsize(src)
        except OSError:
            size = 0
        files.append({
            "role": role,
            "member": (members or {}).get(role) or f"data/{safe_filename(os.path.basename(src))}",
            "basename": os.path.basename(src),
            "size": size,
        })
    st = dict(stats or {})
    return {
        "format": MANIFEST_FORMAT,
        "version": MANIFEST_VERSION,
        "name": name or "dataset",
        "note": note or "",
        "source": source or "",
        "created_at": created_at or datetime.now().isoformat(timespec="seconds"),
        "generator": f"{APP_NAME} {__version__}",
        "stats": {k: int(st.get(k, 0) or 0)
                  for k in ("n_genes", "n_transcripts", "n_seqids", "total_bases")},
        "files": files,
    }


def write_dataset_archive(path: str, paths: Dict[str, str], name: str = "",
                          note: str = "", source: str = "",
                          stats: Optional[Dict[str, int]] = None,
                          progress_cb: ProgressCB = None,
                          cancel_check: CancelCheck = None) -> BundleExport:
    """把当前这套文件打成 zip（带清单），供下次一键导入。

    * 只打 ``paths`` 里真实存在的文件，角色由调用方保证；
    * 包内结构：``gge_dataset.json`` + ``data/<原始文件名>``；
    * 分块写入，因此几个 GB 的基因组也能报进度、也能取消；
    * 中途取消/失败会把半成品删掉，不留坏包。
    """
    missing = [r for r in REQUIRED_ROLES if not paths.get(r)]
    if missing:
        raise GGEError("缺少必填文件：" + "、".join(ROLE_LABELS[r] for r in missing),
                       "基因组 FASTA 与 GFF/GTF 注释是必须的。")
    absent = [f"{ROLE_LABELS[r]}：{p}" for r, p in paths.items()
              if p and not os.path.isfile(p)]
    if absent:
        raise GGEError("以下文件不存在或已被移动：\n" + "\n".join(absent),
                       "重新选择文件后再打包。")

    total = 0
    for role in ROLE_ORDER:
        p = paths.get(role)
        if p:
            try:
                total += os.path.getsize(p)
            except OSError:
                pass

    parent = os.path.dirname(os.path.abspath(path)) or "."
    os.makedirs(parent, exist_ok=True)
    free = _disk_free(parent)
    # 压缩包一般只有原文件的 25%–35%，不要按原大小卡磁盘；小文件时至少留 8 MB
    need = int(total * 0.40) + 8 * 1024 * 1024
    if free and free < need:
        raise GGEError(f"磁盘空间可能不足：{parent} 所在磁盘可用 {human_bytes(free)}，"
                       f"预计需要 {human_bytes(need)}",
                       "换到空间充足的分区，或先清理该磁盘。")

    used: set = set()
    members: Dict[str, str] = {}
    for role in ROLE_ORDER:
        src = paths.get(role)
        if not src:
            continue
        members[role] = _unique_name(f"data/{safe_filename(os.path.basename(src))}", used)
    manifest = build_manifest(name, paths, note=note, source=source, stats=stats,
                              members=members)

    done = 0
    import json as _json
    try:
        with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED, allowZip64=True,
                             compresslevel=6) as zf:
            zf.writestr(MANIFEST_NAME,
                        _json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
            for role in ROLE_ORDER:
                src = paths.get(role)
                if not src:
                    continue
                member = members[role]
                if progress_cb:
                    progress_cb(int(done / max(total, 1) * 100),
                                f"打包 {os.path.basename(src)}（{ROLE_LABELS[role]}）…")
                info = zipfile.ZipInfo(member, date_time=time.localtime()[:6])
                info.compress_type = zipfile.ZIP_DEFLATED
                info.external_attr = 0o600 << 16
                with open(src, "rb") as fh, zf.open(info, "w") as dst:
                    while True:
                        if cancel_check and cancel_check():
                            raise CancelledError()
                        chunk = fh.read(4 * 1024 * 1024)
                        if not chunk:
                            break
                        dst.write(chunk)
                        done += len(chunk)
                        if progress_cb:
                            progress_cb(min(99, int(done / max(total, 1) * 100)),
                                        f"打包 {os.path.basename(src)} "
                                        f"（{human_bytes(done)} / {human_bytes(total)}）")
            if progress_cb:
                progress_cb(99, "写入清单与索引 …")
    except BaseException:
        try:
            os.remove(path)          # 不留半截压缩包，否则下次导入会报"包损坏"
        except OSError:
            pass
        raise

    result = BundleExport(path=path, name=str(manifest["name"]), members=members,
                          total_bytes=total)
    try:
        result.archive_size = os.path.getsize(path)
    except OSError:
        result.archive_size = 0
    if progress_cb:
        progress_cb(100, f"打包完成：{os.path.basename(path)}"
                         f"（{human_bytes(result.archive_size)}）")
    return result


def _manifest_member(members: Sequence[MemberInfo]) -> Optional[MemberInfo]:
    for m in members:
        if os.path.basename(m.name).lower() == MANIFEST_NAME:
            return m
    return None


def _read_manifest_from_archive(kind: str, archive,
                                members: Sequence[MemberInfo]) -> Optional[Dict[str, object]]:
    """读包内清单；不是本程序的清单就返回 None（当作普通压缩包继续猜）。"""
    import json as _json

    member = _manifest_member(members)
    if member is None:
        return None
    try:
        raw = _read_head(kind, archive, member, limit=_MANIFEST_MAX_BYTES)
        data = _json.loads(raw.decode("utf-8", "replace"))
    except Exception:
        return None
    if not isinstance(data, dict) or data.get("format") != MANIFEST_FORMAT:
        return None
    return data


def read_manifest(path: str) -> Optional[Dict[str, object]]:
    """从磁盘上的压缩包读清单（CLI / 测试用）。"""
    if not os.path.isfile(path):
        return None
    kind, archive = _open_archive(path)
    try:
        members, _ignored, _total = _list_members(kind, archive)
        return _read_manifest_from_archive(kind, archive, members)
    finally:
        try:
            archive.close()
        except Exception:
            pass


def _apply_manifest(scan: BundleScan, manifest: Dict[str, object],
                    members: Sequence[MemberInfo]) -> None:
    """按清单直接指定角色：可解释、可复现，也不必再让打分去"猜"。"""
    entries = manifest.get("files") or []
    if not isinstance(entries, list):
        return
    by_member = {c.member: c for c in scan.candidates}
    by_base: Dict[str, FileCandidate] = {}
    for c in scan.candidates:
        by_base.setdefault(os.path.basename(c.member).lower(), c)

    # 先清空：清单只列了 3 个角色时，不能让打分残留的第 4 个角色混进来
    for c in scan.candidates:
        c.role = ""

    won: Dict[str, str] = {}
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        role = str(entry.get("role") or "")
        member = str(entry.get("member") or "")
        if role not in ROLE_ORDER or not member or role in won:
            continue
        cand = by_member.get(member) or by_base.get(os.path.basename(member).lower())
        if cand is None:
            scan.ignored.append((member, "清单里有，但压缩包中没有这个文件（可能被改名或删掉）"))
            continue
        if role in won.values():
            continue
        cand.role = role
        cand.score = 300.0
        cand.confidence = 0.99
        cand.ambiguous = False
        note = "按包内 gge_dataset.json 清单指定（本程序导出时记录的角色）"
        name = str(manifest.get("name") or "").strip()
        cand.evidence = [note + (f"；数据集「{name}」" if name else "")]
        won[role] = member

    # 同一个文件被清单写成两个角色时，只保留第一个，避免"一鱼两吃"
    seen: Dict[str, str] = {}
    for cand in scan.candidates:
        if not cand.role:
            continue
        if cand.member in seen:
            cand.evidence.append(f"清单里同一文件已用作{ROLE_LABELS[seen[cand.member]]}，"
                                 f"此处忽略")
            cand.role = ""
            continue
        seen[cand.member] = cand.role


def _manifest_stats(manifest: Dict[str, object]) -> Dict[str, int]:
    stats = manifest.get("stats") or {}
    if not isinstance(stats, dict):
        return {}
    out: Dict[str, int] = {}
    for key in ("n_genes", "n_transcripts", "n_seqids", "total_bases"):
        try:
            out[key] = int(stats.get(key, 0) or 0)
        except (TypeError, ValueError):
            out[key] = 0
    return out
