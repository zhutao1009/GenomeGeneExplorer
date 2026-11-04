"""坐标系统与链方向工具。

约定（全项目唯一约定，任何模块都不得自行假设）
----------------------------------------------
* GFF3 / GTF / 用户界面 / 导出文本：**1-based 闭区间** ``[start, end]``，
  即 ``start`` 与 ``end`` 对应的碱基都包含在内。
* Python 内部切片 / Region 局部偏移：**0-based 半开区间** ``[start0, end0)``。
* 转换：
    - ``to_0based_closed(start, end) -> (start-1, end)``
    - ``to_1based_closed(start0, end0) -> (start0+1, end0)``
* 负链：转录方向与基因组坐标方向相反，因此"按转录顺序"取片段时必须
  **先倒序排列片段，再对每个片段取反向互补**。
"""

from __future__ import annotations

from typing import Iterable, List, Sequence, Tuple

_COMPLEMENT = str.maketrans(
    "ACGTURYSWKMBDHVNacgturyswkmbdhvn",
    "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn",
)

VALID_STRANDS = ("+", "-")


def revcomp(seq: str) -> str:
    """反向互补。非 ACGT 字符（如 N、简并碱基）按 IUPAC 规则保留并提示。"""
    return seq.translate(_COMPLEMENT)[::-1]


def complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)


def to_0based_closed(start: int, end: int) -> Tuple[int, int]:
    """1-based 闭区间 -> 0-based 半开区间。"""
    return start - 1, end


def to_1based_closed(start0: int, end0: int) -> Tuple[int, int]:
    """0-based 半开区间 -> 1-based 闭区间。"""
    return start0 + 1, end0


def clip_interval(start: int, end: int, seq_len: int) -> Tuple[int, int, bool]:
    """把 1-based 闭区间裁剪到 ``[1, seq_len]``。

    返回 ``(start, end, clipped)``；若完全落在区间外则返回 ``(0, 0, True)``，
    调用方应把长度 0 视为"该侧没有序列"而不是错误。
    """
    clipped = False
    if start < 1:
        start = 1
        clipped = True
    if end > seq_len:
        end = seq_len
        clipped = True
    if start > end or start > seq_len or end < 1:
        return 0, 0, True
    return start, end, clipped


def sort_segments_asc(segments: Iterable) -> List:
    """按基因组坐标升序排序（同起点按终点）。"""
    return sorted(segments, key=lambda s: (s.start, s.end))


def order_for_transcript(segments: Iterable, strand: str) -> List:
    """按**转录方向**排序：正链升序，负链降序。"""
    segs = sort_segments_asc(segments)
    if strand == "-":
        segs.reverse()
    return segs


def local_offsets(genomic_start: int, genomic_end: int, block_start: int,
                  block_end: int, strand: str) -> Tuple[int, int]:
    """把块的基因组坐标区间映射为块内 0-based 半开局部区间。

    ``block_start`` / ``block_end`` 是块自身的 1-based 闭区间端点。
    负链块内部是从基因组右端开始书写的，因此需要翻转。
    """
    if strand == "-":
        s0 = block_end - genomic_end
        e0 = block_end - genomic_start + 1
    else:
        s0 = genomic_start - block_start
        e0 = genomic_end - block_start + 1
    if s0 < 0:
        s0 = 0
    if e0 < s0:
        e0 = s0
    return s0, e0


def fmt_interval(start: int, end: int) -> str:
    return f"{start:,}-{end:,}"


def fmt_span(start: int, end: int) -> str:
    return f"{start:,}..{end:,}"


def human_size(n: int) -> str:
    """把 bp 数变成易读字符串（小写 b 表示碱基）。"""
    if n >= 1_000_000_000:
        return f"{n / 1_000_000_000:.2f} Gb"
    if n >= 1_000_000:
        return f"{n / 1_000_000:.2f} Mb"
    if n >= 1_000:
        return f"{n / 1_000:.1f} kb"
    return f"{n} bp"


def human_bytes(n: int) -> str:
    """把文件字节数变成易读字符串（大写 B 表示 Byte）。

    与 :func:`human_size` 分开：文件大小说成 ``136 bp`` / ``863.00 Mb`` 会让人
    当成序列长度，所以字节一律用 B / KB / MB / GB。
    """
    if n >= 1024 ** 3:
        return f"{n / 1024 ** 3:.2f} GB"
    if n >= 1024 ** 2:
        return f"{n / 1024 ** 2:.1f} MB"
    if n >= 1024:
        return f"{n / 1024:.1f} KB"
    return f"{n} B"


def natural_key(text: str) -> Tuple:
    """自然排序键：chr1 < chr2 < chr10。"""
    import re

    return tuple(
        int(part) if part.isdigit() else part.lower()
        for part in re.split(r"(\d+)", str(text))
    )


def iter_windows(seq: str, width: int) -> List[str]:
    """按固定宽度切分（导出 FASTA 用）。"""
    if width <= 0:
        return [seq]
    return [seq[i:i + width] for i in range(0, len(seq), width)] or [""]


def gc_content(seq: str) -> float:
    if not seq:
        return 0.0
    s = seq.upper()
    gc = s.count("G") + s.count("C")
    return gc / len(s) * 100.0


def base_composition(seq: str) -> dict:
    """碱基组成统计；非 ACGT 归入 other，便于发现异常字符。"""
    s = seq.upper()
    comp = {"A": 0, "C": 0, "G": 0, "T": 0, "other": 0}
    for ch in s:
        if ch in "ACGT":
            comp[ch] += 1
        else:
            comp["other"] += 1
    return comp


def unusual_bases(seq: str) -> str:
    """返回序列中出现的非 ACGT 字符集合（排序后）。"""
    return "".join(sorted({c for c in seq.upper() if c not in "ACGT"}))


def wrap_sequence(seq: str, width: int = 70) -> Sequence[str]:
    return [seq[i:i + width] for i in range(0, len(seq), width)]


def merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """合并 0-based 半开区间（用于把相邻的 CDS 片段连通）。"""
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]


def intersect(a: Tuple[int, int], b: Tuple[int, int]) -> Tuple[int, int]:
    """两个 0-based 半开区间求交，空则返回 (0, 0)。"""
    s = max(a[0], b[0])
    e = min(a[1], b[1])
    if e <= s:
        return 0, 0
    return s, e
