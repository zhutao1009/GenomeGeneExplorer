"""序列层面的小工具：翻译、密码子判断、碱基统计。"""

from __future__ import annotations

from typing import List, Optional, Tuple

_BASES = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
STANDARD_TABLE = {
    a + b + c: _AA[i * 16 + j * 4 + k]
    for i, a in enumerate(_BASES)
    for j, b in enumerate(_BASES)
    for k, c in enumerate(_BASES)
}

START_CODONS = {"ATG"}
ALT_START_CODONS = {"ATG", "GTG", "TTG", "CTG", "ATT", "ATC", "ATA"}
STOP_CODONS = {"TAA", "TAG", "TGA"}


def translate(seq: str, stop_at_stop: bool = False) -> str:
    """按标准遗传密码翻译；末尾不足 3 nt 的部分忽略。

    未知密码子（含 N 等简并碱基）翻译为 'X'，不抛异常——真实基因组里
    简并碱基很常见，报错会打断整个流程。
    """
    seq = seq.upper()
    out: List[str] = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        aa = STANDARD_TABLE.get(codon, "X")
        if aa == "*" and stop_at_stop:
            break
        out.append(aa)
    return "".join(out)


def codons(seq: str) -> List[str]:
    s = seq.upper()
    return [s[i:i + 3] for i in range(0, len(s) - 2, 3)]


def first_codon(seq: str) -> str:
    return seq[:3].upper()


def last_codon(seq: str) -> str:
    if len(seq) < 3:
        return seq.upper()
    return seq[-3:].upper()


def is_start_codon(seq: str, allow_alt: bool = False) -> bool:
    return first_codon(seq) in (ALT_START_CODONS if allow_alt else START_CODONS)


def is_stop_codon(seq: str) -> bool:
    return last_codon(seq) in STOP_CODONS


def internal_stops(seq: str) -> List[Tuple[int, str]]:
    """返回内部终止密码子 [(密码子序号从 1 开始, 密码子)]，不含最后一个。"""
    cs = codons(seq)
    out = []
    for i, c in enumerate(cs[:-1] if cs else []):
        if c in STOP_CODONS:
            out.append((i + 1, c))
    return out


def reading_frames(seq: str) -> List[str]:
    """三个正框的翻译结果（用于不完整 CDS 的排查）。"""
    return [translate(seq[off:]) for off in range(3)]


def orf_bounds(seq: str) -> Optional[Tuple[int, int]]:
    """在给定序列中找最长的 ATG..终止 的 ORF，返回 (start0, end0) 半开区间。"""
    seq = seq.upper()
    best: Optional[Tuple[int, int]] = None
    best_len = 0
    for frame in range(3):
        i = frame
        start = -1
        while i + 3 <= len(seq):
            codon = seq[i:i + 3]
            if start < 0 and codon in START_CODONS:
                start = i
            elif start >= 0 and codon in STOP_CODONS:
                length = i + 3 - start
                if length > best_len:
                    best_len = length
                    best = (start, i + 3)
                start = -1
            i += 3
    return best


def gc_percent(seq: str) -> float:
    if not seq:
        return 0.0
    s = seq.upper()
    return (s.count("G") + s.count("C")) / len(s) * 100.0


def molecular_weight_da(seq: str) -> float:
    """双链 DNA 分子量（近似，用于报告展示）。"""
    s = seq.upper()
    return s.count("A") * 313.21 + s.count("T") * 304.2 + \
        s.count("G") * 329.21 + s.count("C") * 289.18
