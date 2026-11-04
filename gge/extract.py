"""序列提取核心算法。

这一层完全不依赖 Qt，可单独单元测试；GUI 只负责把 ``GeneExtraction`` 渲染出来。

坐标与方向的关键约定
--------------------
* ``Segment`` 是 1-based 闭区间（GFF/GTF 原始值）。
* ``fetch`` 返回**正链**序列；负链片段必须 ``revcomp`` 后才等于转录方向序列。
* 负链转录本"按转录顺序"= 基因组坐标降序；因此排序后要整段倒序，
  再逐段反向互补。二者缺一都会得到错误序列（常见 bug）。
* ``Region`` 的 ``start0/end0`` 是**块内**偏移，不是基因组坐标。
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

from .coords import (base_composition, clip_interval, gc_content, local_offsets,
                     order_for_transcript, revcomp, unusual_bases)
from .errors import AnnotationError, CoordinateError, SeqIDMismatchError
from .models import (REGION_CDS, REGION_DOWNSTREAM, REGION_EXON, REGION_FLANK,
                     REGION_INTRON, REGION_OTHER, REGION_PRIORITY, REGION_UPSTREAM,
                     REGION_UTR3, REGION_UTR5, Block, CheckItem, Gene, GeneExtraction,
                     Panel, Region, Segment, Transcript)
from .sequtil import (first_codon, internal_stops, is_start_codon, is_stop_codon,
                      last_codon, translate)

# 面板 key（GUI 与导出共用，顺序即界面顺序）
PANEL_GENOMIC = "genomic"
PANEL_UPSTREAM = "upstream"
PANEL_DOWNSTREAM = "downstream"
PANEL_EXON = "exon"
PANEL_INTRON = "intron"
PANEL_MRNA = "mrna"
PANEL_CDS = "cds"
PANEL_UTR = "utr"

PANEL_ORDER = [PANEL_GENOMIC, PANEL_UPSTREAM, PANEL_DOWNSTREAM, PANEL_EXON,
               PANEL_INTRON, PANEL_MRNA, PANEL_CDS, PANEL_UTR]

PANEL_TITLES = {
    PANEL_GENOMIC: "基因组连续序列",
    PANEL_UPSTREAM: "上游序列",
    PANEL_DOWNSTREAM: "下游序列",
    PANEL_EXON: "外显子序列",
    PANEL_INTRON: "内含子序列",
    PANEL_MRNA: "拼接转录本 mRNA",
    PANEL_CDS: "CDS 编码序列",
    PANEL_UTR: "UTR 序列",
}

SEQ_TYPE_FASTA = {
    PANEL_GENOMIC: "genomic",
    PANEL_UPSTREAM: "upstream",
    PANEL_DOWNSTREAM: "downstream",
    PANEL_EXON: "exon",
    PANEL_INTRON: "intron",
    PANEL_MRNA: "mRNA",
    PANEL_CDS: "CDS",
    PANEL_UTR: "UTR",
}

# 渲染用的类别编码（顺序即优先级，后者覆盖前者）
_KIND_ORDER = [REGION_FLANK, REGION_UPSTREAM, REGION_DOWNSTREAM, REGION_INTRON,
               REGION_EXON, REGION_UTR5, REGION_UTR3, REGION_CDS]
_KIND_CODE = {k: i for i, k in enumerate(_KIND_ORDER)}


# --------------------------------------------------------------------------- #
# 辅助
# --------------------------------------------------------------------------- #
def kinds_to_runs(codes: bytearray, offset: int = 0) -> List[Tuple[int, int, str]]:
    """把逐碱基的类别数组压缩成连续区间。"""
    runs: List[Tuple[int, int, str]] = []
    if not codes:
        return runs
    start = 0
    cur = codes[0]
    for i in range(1, len(codes)):
        if codes[i] != cur:
            runs.append((start + offset, i + offset, _KIND_ORDER[cur]))
            start = i
            cur = codes[i]
    runs.append((start + offset, len(codes) + offset, _KIND_ORDER[cur]))
    return runs


def paint_overlays(length: int, overlays: Sequence[Tuple[int, int, str]],
                   base_kind: str = REGION_FLANK) -> bytearray:
    """按优先级把若干区间"盖"到逐碱基数组上（存在重叠时高优先级胜出）。"""
    codes = bytearray([_KIND_CODE[base_kind]]) * length
    ordered = sorted(overlays, key=lambda o: _KIND_CODE[o[2]])
    for s0, e0, kind in ordered:
        s0 = max(0, min(int(s0), length))
        e0 = max(0, min(int(e0), length))
        if e0 <= s0:
            continue
        codes[s0:e0] = bytes([_KIND_CODE[kind]]) * (e0 - s0)
    return codes


def _regions_from_runs(runs: Sequence[Tuple[int, int, str]],
                       block_genomic_start: int, block_genomic_end: int,
                       strand: str,
                       labels: Optional[Dict[str, str]] = None) -> List[Region]:
    """把块内 run 转成带基因组坐标的 Region 列表。"""
    out: List[Region] = []
    labels = labels or {}
    for s0, e0, kind in runs:
        reg = Region(start0=s0, end0=e0, kind=kind, label=labels.get(kind, ""))
        if strand == "-":
            reg.genomic_start = block_genomic_end - (e0 - 1)
            reg.genomic_end = block_genomic_end - s0
        else:
            reg.genomic_start = block_genomic_start + s0
            reg.genomic_end = block_genomic_start + e0 - 1
        out.append(reg)
    return out


def _fetch_segment(provider, seqid: str, seg: Segment, strand: str) -> str:
    seq = provider.fetch(seqid, seg.start, seg.end)
    return revcomp(seq) if strand == "-" else seq


def _intron_between(a: Segment, b: Segment) -> Optional[Segment]:
    """两个相邻外显子之间的内含子（基因组方向升序）。

    若两个外显子在基因组上重叠（注释异常），则不存在内含子，返回 None，
    由调用方记录 warning，而不是造出一个负长度区间。
    """
    gap_start = min(a.end, b.end) + 1
    gap_end = max(a.start, b.start) - 1
    if gap_end < gap_start:
        return None
    return Segment(start=gap_start, end=gap_end)


def _exon_local_range(exon: Segment, gs: int, ge: int, strand: str) -> Tuple[int, int]:
    """外显子内部基因组区间 -> 外显子块内 0-based 半开区间。"""
    return local_offsets(gs, ge, exon.start, exon.end, strand)


# --------------------------------------------------------------------------- #
# 主函数
# --------------------------------------------------------------------------- #
def extract_gene(provider,
                 gene: Gene,
                 tx: Transcript,
                 upstream: int = 2000,
                 downstream: int = 2000,
                 seqid_override: Optional[str] = None) -> GeneExtraction:
    """提取一个转录本的全部序列与结构信息。"""
    seqid = seqid_override or tx.seqid
    chrom_len = provider.length(seqid)
    if chrom_len is None:
        cands = provider.suggest(seqid)
        raise SeqIDMismatchError(
            f"基因组 FASTA 中找不到序列 “{seqid}”。",
            "注释文件里的序列名与 FASTA 不一致（常见：chr1 vs 1 vs NC_000001.11）。"
            + (f" FASTA 中相近的序列名：{', '.join(cands)}。" if cands else
               " 请确认 FASTA 与注释来自同一组装版本。"),
        )

    exons = tx.ordered_exons()
    if not exons:
        raise AnnotationError(
            f"转录本 {tx.tx_id} 没有任何外显子（也未提供可用 CDS）。",
            "检查注释文件中该转录本的 exon/CDS 行是否被正确解析。",
        )

    tx_start, tx_end = tx.span()
    warnings: List[str] = []
    checks: List[CheckItem] = []

    # ---------------- 坐标合法性 ---------------- #
    if tx_start < 1 or tx_end > chrom_len:
        warnings.append(
            f"转录本坐标 {tx_start}-{tx_end} 超出 {seqid} 长度 {chrom_len:,}，"
            f"超出部分会被截断（原注释可能有误）"
        )
        checks.append(CheckItem("warn", "坐标越界",
                                f"{seqid} 长度 {chrom_len:,}，转录本区间 {tx_start}-{tx_end}"))
    else:
        checks.append(CheckItem("ok", "坐标合法",
                                f"{seqid}:{tx_start:,}-{tx_end:,}（长度 {chrom_len:,}）"))
    tx_start, tx_end, _ = clip_interval(tx_start, tx_end, chrom_len)
    if tx_start == 0:
        raise CoordinateError(
            f"转录本 {tx.tx_id} 的坐标完全落在 {seqid} 之外。",
            "检查注释与 FASTA 是否来自同一组装版本。",
        )

    # ---------------- 上下游区间（按转录方向） ---------------- #
    up_req = max(0, int(upstream))
    dn_req = max(0, int(downstream))
    if tx.strand == "-":
        # 负链：转录起点在基因组右侧 → 上游在 tx_end 的右方，下游在 tx_start 的左方
        up_raw = (tx_end + 1, tx_end + up_req) if up_req else (0, 0)
        dn_raw = (tx_start - dn_req, tx_start - 1) if dn_req else (0, 0)
    else:
        up_raw = (tx_start - up_req, tx_start - 1) if up_req else (0, 0)
        dn_raw = (tx_end + 1, tx_end + dn_req) if dn_req else (0, 0)

    up_s, up_e, up_clip = clip_interval(up_raw[0], up_raw[1], chrom_len) if up_req else (0, 0, False)
    dn_s, dn_e, dn_clip = clip_interval(dn_raw[0], dn_raw[1], chrom_len) if dn_req else (0, 0, False)

    up_seq = ""
    if up_s:
        up_seq = provider.fetch(seqid, up_s, up_e)
        if tx.strand == "-":
            up_seq = revcomp(up_seq)
    dn_seq = ""
    if dn_s:
        dn_seq = provider.fetch(seqid, dn_s, dn_e)
        if tx.strand == "-":
            dn_seq = revcomp(dn_seq)

    if up_req and len(up_seq) < up_req:
        warnings.append(
            f"上游请求 {up_req:,} bp，实际只取到 {len(up_seq):,} bp"
            f"（已到达 {seqid} 边界）"
        )
    if dn_req and len(dn_seq) < dn_req:
        warnings.append(
            f"下游请求 {dn_req:,} bp，实际只取到 {len(dn_seq):,} bp"
            f"（已到达 {seqid} 边界）"
        )

    # ---------------- 区域边界 ---------------- #
    # 四个子区间都要参与取极值：负链基因的上游在基因组右侧、下游在左侧，
    # 只算 (tx_start, up_s) 会把上游丢掉、并让下游偏移变成负数。
    region_start = min([x for x in (tx_start, up_s, dn_s) if x] or [tx_start])
    region_end = max([x for x in (tx_end, up_e, dn_e) if x] or [tx_end])
    region_seq = provider.fetch(seqid, region_start, region_end)

    # ---------------- 逐碱基类别（mRNA 坐标） ---------------- #
    # 1) 外显子在转录本坐标上的偏移
    exon_entries: List[Tuple[Segment, int]] = []
    acc = 0
    for ex in exons:
        exon_entries.append((ex, acc))
        acc += ex.length
    mrna_len = acc
    mrna_codes = bytearray([_KIND_CODE[REGION_EXON]]) * mrna_len

    # 2) CDS 覆盖位置（注意：CDS 与 exon 的交集才属于成熟 mRNA）
    cds_in_exon: List[Tuple[Segment, Segment]] = []  # (cds, exon)
    for cds in tx.ordered_cds():
        for ex, off in exon_entries:
            gs = max(cds.start, ex.start)
            ge = min(cds.end, ex.end)
            if ge < gs:
                continue
            cds_in_exon.append((cds, ex))
            ls, le = _exon_local_range(ex, gs, ge, tx.strand)
            mrna_codes[off + ls:off + le] = bytes([_KIND_CODE[REGION_CDS]]) * (le - ls)

    first_cds_pos = mrna_codes.find(_KIND_CODE[REGION_CDS])
    last_cds_pos = mrna_codes.rfind(_KIND_CODE[REGION_CDS])
    has_cds = first_cds_pos >= 0
    if has_cds:
        if first_cds_pos > 0:
            mrna_codes[0:first_cds_pos] = bytes([_KIND_CODE[REGION_UTR5]]) * first_cds_pos
        if last_cds_pos + 1 < mrna_len:
            n = mrna_len - (last_cds_pos + 1)
            mrna_codes[last_cds_pos + 1:] = bytes([_KIND_CODE[REGION_UTR3]]) * n

    # 3) 组装 mRNA 与各段序列
    exon_blocks: List[Block] = []
    mrna_blocks: List[Block] = []
    mrna_parts: List[str] = []
    for idx, (ex, off) in enumerate(exon_entries, start=1):
        seq = _fetch_segment(provider, seqid, ex, tx.strand)
        mrna_parts.append(seq)
        local_runs = kinds_to_runs(mrna_codes[off:off + ex.length])
        regions = _regions_from_runs(local_runs, ex.start, ex.end, tx.strand,
                                     {REGION_EXON: "外显子", REGION_CDS: "CDS",
                                      REGION_UTR5: "5′UTR", REGION_UTR3: "3′UTR"})
        note = f"外显子 {idx}/{len(exon_entries)} · {seqid}:{ex.start:,}-{ex.end:,} · {ex.length:,} bp"
        exon_blocks.append(Block(
            key=f"exon{idx}", title=f"外显子 {idx}", seq=seq, regions=regions, note=note,
            meta={"strand": tx.strand, "coord_mode": "genomic",
                  "coord_start": ex.start, "coord_reverse": tx.strand == "-",
                  "genomic_start": ex.start, "genomic_end": ex.end,
                  "exon_index": idx, "length": ex.length},
        ))
        mrna_blocks.append(Block(
            key=f"mrna_seg{idx}", title=f"剪接片段 {idx}（外显子 {idx}）", seq=seq,
            regions=[Region(0, len(seq), REGION_EXON, "外显子")],
            note=f"转录本第 {off + 1:,}-{off + ex.length:,} 位 · 来自 {seqid}:{ex.start:,}-{ex.end:,}",
            meta={"strand": tx.strand, "coord_mode": "transcript",
                  "coord_start": off + 1, "coord_reverse": False,
                  "genomic_start": ex.start, "genomic_end": ex.end,
                  "exon_index": idx, "length": ex.length},
        ))
    mrna = "".join(mrna_parts)

    # ---------------- 内含子 ---------------- #
    intron_blocks: List[Block] = []
    for i in range(len(exons) - 1):
        a, b = exons[i], exons[i + 1]
        intron = _intron_between(a, b)
        if intron is None:
            warnings.append(
                f"外显子 {i + 1} 与 {i + 2} 在基因组上重叠，未生成内含子 {i + 1}"
            )
            continue
        seq = _fetch_segment(provider, seqid, intron, tx.strand)
        intron_blocks.append(Block(
            key=f"intron{i + 1}", title=f"内含子 {i + 1}", seq=seq,
            regions=[Region(0, len(seq), REGION_INTRON, "内含子")],
            note=(f"内含子 {i + 1}/{len(exons) - 1} · {seqid}:{intron.start:,}-{intron.end:,} · "
                  f"{intron.length:,} bp" + ("（负链，已反向互补）" if tx.strand == "-" else "")),
            meta={"strand": tx.strand, "coord_mode": "genomic",
                  "coord_start": intron.start, "coord_reverse": tx.strand == "-",
                  "genomic_start": intron.start, "genomic_end": intron.end,
                  "intron_index": i + 1, "length": intron.length},
        ))

    # ---------------- CDS ---------------- #
    cds_blocks: List[Block] = []
    cds_parts: List[str] = []
    cds_ordered = tx.ordered_cds()
    cum = 0
    phase_report: List[str] = []
    for idx, seg in enumerate(cds_ordered, start=1):
        seq = _fetch_segment(provider, seqid, seg, tx.strand)
        cds_parts.append(seq)
        exp_phase = (3 - (cum % 3)) % 3 if cum else 0
        note_bits = [f"CDS {idx}/{len(cds_ordered)}", f"{seqid}:{seg.start:,}-{seg.end:,}",
                     f"{seg.length:,} bp"]
        if tx.strand == "-":
            note_bits.append("负链，已反向互补")
        if seg.phase >= 0:
            note_bits.append(f"phase={seg.phase}")
            if seg.phase != exp_phase and idx > 1:
                phase_report.append(
                    f"CDS 片段 {idx} 注释 phase={seg.phase}，按拼接长度应为 {exp_phase}"
                )
        else:
            note_bits.append("phase 缺失")
        cds_blocks.append(Block(
            key=f"cds{idx}", title=f"CDS 片段 {idx}", seq=seq,
            regions=[Region(0, len(seq), REGION_CDS, "CDS")],
            note=" · ".join(note_bits),
            meta={"strand": tx.strand, "coord_mode": "genomic",
                  "coord_start": seg.start, "coord_reverse": tx.strand == "-",
                  "genomic_start": seg.start, "genomic_end": seg.end,
                  "cds_index": idx, "phase": seg.phase,
                  "frame_offset": cum % 3, "length": seg.length},
        ))
        cum += seg.length
    cds = "".join(cds_parts)
    if phase_report:
        warnings.append("；".join(phase_report[:3]) +
                        ("等" if len(phase_report) > 3 else "") +
                        "。已按注释原样拼接，未做任何碱基裁剪")

    # ---------------- UTR ---------------- #
    # 注释文件（RefSeq/NCBI 的 GFF）通常**没有** UTR 行，这里全部是推导值：
    # 转录本坐标上，首个 CDS 起点之前 = 5′UTR，末个 CDS 终点之后 = 3′UTR。
    utr_blocks: List[Block] = []
    utr5_parts: List[str] = []
    utr3_parts: List[str] = []
    utr_in_annotation = bool(tx.utr5_features or tx.utr3_features)
    if has_cds:
        # 转录方向上的首个 CDS 起点 / 末个 CDS 终点（基因组坐标），用于说明推导依据
        first_cds_bp = cds_ordered[0].start if tx.strand == "+" else cds_ordered[0].end
        last_cds_bp = cds_ordered[-1].end if tx.strand == "+" else cds_ordered[-1].start
        utr_src = ("注释有 UTR 行，本面板仍按 exon−CDS 推导"
                   if utr_in_annotation else "注释中不存在 UTR 行，此为推导值")
        for ex, off in exon_entries:
            slice_codes = mrna_codes[off:off + ex.length]
            # 逐段取 5'UTR / 3'UTR 片段
            for target_kind, store, label in ((REGION_UTR5, utr5_parts, "5′ UTR"),
                                              (REGION_UTR3, utr3_parts, "3′ UTR")):
                code = _KIND_CODE[target_kind]
                i = 0
                n = len(slice_codes)
                while i < n:
                    if slice_codes[i] != code:
                        i += 1
                        continue
                    j = i
                    while j < n and slice_codes[j] == code:
                        j += 1
                    if tx.strand == "-":
                        gs = ex.end - (j - 1)
                        ge = ex.end - i
                    else:
                        gs = ex.start + i
                        ge = ex.start + j - 1
                    seq = provider.fetch(seqid, gs, ge)
                    if tx.strand == "-":
                        seq = revcomp(seq)
                    store.append(seq)
                    basis = (f"末个 CDS 终点 {last_cds_bp:,} + 1 起算"
                             if target_kind == REGION_UTR3
                             else f"首个 CDS 起点 {first_cds_bp:,} − 1 为止")
                    utr_blocks.append(Block(
                        key=f"{target_kind}_{len(utr_blocks)}",
                        title=f"{label}（推导）",
                        seq=seq, regions=[Region(0, len(seq), target_kind, label)],
                        note=f"{label}（推导） · {seqid}:{gs:,}-{ge:,} · {len(seq):,} bp "
                             f"· 转录本第 {off + i + 1:,}-{off + j:,} 位 · "
                             f"{utr_src}；{basis}",
                        meta={"strand": tx.strand, "coord_mode": "genomic",
                              "coord_start": gs, "coord_reverse": tx.strand == "-",
                              "genomic_start": gs, "genomic_end": ge,
                              "length": len(seq), "utr_kind": target_kind},
                    ))
                    i = j
    utr5 = "".join(utr5_parts)
    utr3 = "".join(utr3_parts)

    # ---------------- 基因组连续序列面板 ---------------- #
    overlays: List[Tuple[int, int, str]] = []
    if up_s:
        s0, e0 = up_s - region_start, up_e - region_start + 1
        overlays.append((s0, e0, REGION_UPSTREAM))
    if dn_s:
        s0, e0 = dn_s - region_start, dn_e - region_start + 1
        overlays.append((s0, e0, REGION_DOWNSTREAM))
    for ex in exons:
        s0 = ex.start - region_start
        e0 = ex.end - region_start + 1
        overlays.append((s0, e0, REGION_EXON))
    for i in range(len(exons) - 1):
        intron = _intron_between(exons[i], exons[i + 1])
        if intron:
            overlays.append((intron.start - region_start, intron.end - region_start + 1,
                             REGION_INTRON))
    for cds_seg, ex in cds_in_exon:
        gs, ge = max(cds_seg.start, ex.start), min(cds_seg.end, ex.end)
        overlays.append((gs - region_start, ge - region_start + 1, REGION_CDS))
    if has_cds:
        for utr_kind in (REGION_UTR5, REGION_UTR3):
            for block in utr_blocks:
                if block.meta.get("utr_kind") != utr_kind:
                    continue
                gs = int(block.meta["genomic_start"])
                ge = int(block.meta["genomic_end"])
                overlays.append((gs - region_start, ge - region_start + 1, utr_kind))
    codes = paint_overlays(len(region_seq), overlays, base_kind=REGION_FLANK)
    genomic_regions = _regions_from_runs(kinds_to_runs(codes), region_start,
                                         region_end, "+",
                                         {REGION_FLANK: "基因间", REGION_UPSTREAM: "上游",
                                          REGION_DOWNSTREAM: "下游", REGION_EXON: "外显子",
                                          REGION_INTRON: "内含子", REGION_CDS: "CDS",
                                          REGION_UTR5: "5′UTR", REGION_UTR3: "3′UTR"})
    genomic_block = Block(
        key="genomic", title="基因组窗口", seq=region_seq, regions=genomic_regions,
        note=(f"{seqid}:{region_start:,}-{region_end:,}（{len(region_seq):,} bp，"
              f"以 {seqid} 正链书写；基因在 {tx.strand} 链）"),
        meta={"strand": "+", "coord_mode": "genomic", "coord_start": region_start,
              "coord_reverse": False, "genomic_start": region_start,
              "genomic_end": region_end, "length": len(region_seq)},
    )

    # ---------------- 面板集合 ---------------- #
    panels: Dict[str, Panel] = {}
    panels[PANEL_GENOMIC] = Panel(
        PANEL_GENOMIC, PANEL_TITLES[PANEL_GENOMIC], [genomic_block],
        note=f"含上下游的连续基因组序列（{tx.strand} 链基因，此处按正链参考显示）",
        seq_type="genomic",
    )
    panels[PANEL_UPSTREAM] = Panel(
        PANEL_UPSTREAM, PANEL_TITLES[PANEL_UPSTREAM],
        [Block("up", f"上游 {len(up_seq):,} bp", up_seq,
               [Region(0, len(up_seq), REGION_UPSTREAM, "上游")] if up_seq else [],
               note=(f"{seqid}:{up_raw[0]:,}-{up_raw[1]:,}（{tx.strand} 链转录方向）"
                     if up_seq else "未取到上游序列"),
               meta={"strand": tx.strand, "coord_mode": "genomic",
                     "coord_start": up_s or 1, "coord_reverse": tx.strand == "-",
                     "genomic_start": up_s, "genomic_end": up_e, "length": len(up_seq)})]
        if up_seq else [],
        note="按转录方向定义：负链基因的上游位于基因组坐标更大的一侧，且已反向互补",
        seq_type="upstream", empty_text="上游长度设为 0 或已到达序列边界",
    )
    panels[PANEL_DOWNSTREAM] = Panel(
        PANEL_DOWNSTREAM, PANEL_TITLES[PANEL_DOWNSTREAM],
        [Block("down", f"下游 {len(dn_seq):,} bp", dn_seq,
               [Region(0, len(dn_seq), REGION_DOWNSTREAM, "下游")] if dn_seq else [],
               note=(f"{seqid}:{dn_raw[0]:,}-{dn_raw[1]:,}（{tx.strand} 链转录方向）"
                     if dn_seq else "未取到下游序列"),
               meta={"strand": tx.strand, "coord_mode": "genomic",
                     "coord_start": dn_s or 1, "coord_reverse": tx.strand == "-",
                     "genomic_start": dn_s, "genomic_end": dn_e, "length": len(dn_seq)})]
        if dn_seq else [],
        note="按转录方向定义：负链基因的下游位于基因组坐标更小的一侧，且已反向互补",
        seq_type="downstream", empty_text="下游长度设为 0 或已到达序列边界",
    )
    panels[PANEL_EXON] = Panel(
        PANEL_EXON, PANEL_TITLES[PANEL_EXON], exon_blocks,
        note=f"共 {len(exon_blocks)} 个外显子，按转录方向排序（{tx.strand} 链）",
        seq_type="exon",
    )
    panels[PANEL_INTRON] = Panel(
        PANEL_INTRON, PANEL_TITLES[PANEL_INTRON], intron_blocks,
        note=(f"共 {len(intron_blocks)} 个内含子，按转录方向排序" if intron_blocks
              else "该转录本为单外显子结构，没有内含子"),
        seq_type="intron", empty_text="该转录本为单外显子结构，没有内含子",
    )
    panels[PANEL_MRNA] = Panel(
        PANEL_MRNA, PANEL_TITLES[PANEL_MRNA], mrna_blocks,
        note=f"外显子按 {tx.strand} 链转录方向拼接，共 {len(mrna):,} bp；坐标列为转录本内位置",
        seq_type="mRNA", empty_text="无外显子序列",
    )
    panels[PANEL_CDS] = Panel(
        PANEL_CDS, PANEL_TITLES[PANEL_CDS], cds_blocks,
        note=(f"{len(cds_blocks)} 个 CDS 片段拼接，共 {len(cds):,} bp；"
              f"框长 {len(cds) % 3 if cds else 0}；按注释原样拼接，未按 phase 裁剪碱基"
              if cds else "注释中没有 CDS（该转录本可能为非编码 RNA）"),
        seq_type="CDS", empty_text="该转录本在注释中没有 CDS",
    )
    if not has_cds:
        utr_note = "无 CDS，无法定义 UTR"
    elif not utr_blocks:
        utr_note = ("本转录本的 CDS 覆盖了全部外显子，按「外显子 − CDS」推导没有 UTR 序列"
                    "（注释中也没有 UTR 行）" if not utr_in_annotation else
                    "本转录本的 CDS 覆盖了全部外显子，推导结果没有 UTR 序列")
    elif utr_in_annotation:
        utr_note = (f"注释中另有 UTR 行（5′ {len(tx.utr5_features)} 段 / 3′ "
                    f"{len(tx.utr3_features)} 段），但本面板未使用它们，仍按「外显子 − CDS」推导："
                    f"5′UTR {len(utr5):,} + CDS {len(cds):,} + 3′UTR {len(utr3):,} "
                    f"= mRNA {len(mrna):,} bp")
    else:
        utr_note = (f"⚠ 注释中不存在 UTR 记录（本注释只提供 exon/CDS），以下为推导值："
                    f"5′UTR {len(utr5):,} + CDS {len(cds):,} + 3′UTR {len(utr3):,} "
                    f"= mRNA {len(mrna):,} bp。算法：外显子范围内，首个 CDS 起点之前的碱基记 5′UTR，"
                    f"末个 CDS 终点之后的记 3′UTR")

    panels[PANEL_UTR] = Panel(
        PANEL_UTR, f"{PANEL_TITLES[PANEL_UTR]}（推导）", utr_blocks,
        note=utr_note,
        seq_type="UTR", empty_text="无 CDS 或该转录本没有 UTR",
        badge=("注释有 UTR 行" if utr_in_annotation else "注释无 UTR"),
    )

    # ---------------- 自检项 ---------------- #
    checks.extend(_run_checks(provider, gene, tx, exons, mrna, cds, utr5, utr3,
                             has_cds, uniq_bases=_collect_unusual(
                                 [region_seq, mrna, cds, up_seq, dn_seq]), chrom_len=chrom_len))

    for gmsg in gene.warnings:
        warnings.append(f"基因层面：{gmsg}")
    for tmsg in tx.warnings:
        warnings.append(f"转录本层面：{tmsg}")

    result = GeneExtraction(
        gene_id=gene.gene_id,
        gene_name=gene.gene_name or gene.gene_id,
        tx_id=tx.tx_id,
        seqid=seqid,
        strand=tx.strand,
        gene_start=gene.start or tx_start,
        gene_end=gene.end or tx_end,
        region_start=region_start,
        region_end=region_end,
        upstream_requested=up_req,
        downstream_requested=dn_req,
        upstream_actual=len(up_seq),
        downstream_actual=len(dn_seq),
        upstream_clipped=bool(up_req and len(up_seq) < up_req),
        downstream_clipped=bool(dn_req and len(dn_seq) < dn_req),
        upstream_seq=up_seq,
        downstream_seq=dn_seq,
        mrna=mrna,
        cds=cds,
        utr5=utr5,
        utr3=utr3,
        utr_in_annotation=utr_in_annotation,
        panels=panels,
        warnings=warnings,
        checks=checks,
        fasta_path=getattr(provider, "path", ""),
        annotation_path="",
        tx_count=gene.n_transcripts,
        n_exons=len(exons),
        n_introns=len(intron_blocks),
        n_cds=len(cds_blocks),
    )
    return result


def _collect_unusual(seqs: Sequence[str]) -> str:
    chars = set()
    for s in seqs:
        chars |= set(unusual_bases(s))
    return "".join(sorted(chars))


def _run_checks(provider, gene: Gene, tx: Transcript, exons: List[Segment],
                mrna: str, cds: str, utr5: str, utr3: str, has_cds: bool,
                uniq_bases: str, chrom_len: int) -> List[CheckItem]:
    checks: List[CheckItem] = []

    # 结构
    if len(exons) == 1:
        checks.append(CheckItem("info", "单外显子", "该转录本只有 1 个外显子，没有内含子"))
    else:
        checks.append(CheckItem("ok", "多外显子结构",
                                f"{len(exons)} 个外显子，{len(exons) - 1} 个内含子"))

    expected = sum(e.length for e in exons)
    if len(mrna) == expected:
        checks.append(CheckItem("ok", "mRNA 长度一致",
                                f"拼接长度 {len(mrna):,} = 各外显子之和 {expected:,}"))
    else:
        checks.append(CheckItem("error", "mRNA 长度异常",
                                f"拼接长度 {len(mrna):,} ≠ 外显子之和 {expected:,}"))

    # CDS
    if not has_cds:
        checks.append(CheckItem("info", "无 CDS", "注释中该转录本没有 CDS，可能为非编码转录本"))
    else:
        frame = len(cds) % 3
        if frame == 0:
            checks.append(CheckItem("ok", "CDS 阅读框完整",
                                    f"{len(cds):,} bp = {len(cds) // 3:,} 个密码子"))
        else:
            checks.append(CheckItem(
                "warn", "CDS 长度不是 3 的整数倍",
                f"{len(cds):,} bp，余 {frame} bp。可能是 5′/3′ 端不完整注释或坐标偏差，"
                f"软件未修改原序列",
            ))
        if len(cds) >= 3:
            fc = first_codon(cds)
            if is_start_codon(cds):
                checks.append(CheckItem("ok", "起始密码子", f"ATG（{fc}）"))
            elif is_start_codon(cds, allow_alt=True):
                checks.append(CheckItem("warn", "非典型起始密码子",
                                        f"{fc}，不是标准 ATG（部分物种/注释允许）"))
            else:
                checks.append(CheckItem(
                    "warn", "缺少 ATG 起始密码子",
                    f"实际为 {fc}，提示 CDS 5′ 端不完整、注释坐标偏移或链方向错误",
                ))
            lc = last_codon(cds)
            if is_stop_codon(cds):
                checks.append(CheckItem("ok", "终止密码子", f"{lc}"))
            else:
                checks.append(CheckItem(
                    "warn", "缺少终止密码子",
                    f"实际为 {lc}，提示 CDS 3′ 端不完整或注释坐标偏移",
                ))
            stops = internal_stops(cds)
            if stops:
                sample = "、".join(f"第 {i} 个密码子 {c}" for i, c in stops[:3])
                checks.append(CheckItem(
                    "warn", "存在内部终止密码子",
                    f"共 {len(stops)} 处：{sample}。可能是假基因、注释错误或移码",
                ))
            else:
                checks.append(CheckItem("ok", "无内部终止密码子", "按标准遗传密码翻译连续"))

        # CDS 是否落在 exon 内（逐个片段检查，避免"CDS 跨外显子空缺"被漏掉）
        for seg in tx.cds:
            if not any(seg.overlaps(e.start, e.end) for e in exons):
                checks.append(CheckItem(
                    "error", "CDS 不在外显子内",
                    f"CDS {seg.start:,}-{seg.end:,} 与任何外显子都不重叠，注释可能不一致",
                ))
                break
        if has_cds and utr5 == "" and utr3 == "" and len(cds) == len(mrna):
            checks.append(CheckItem("info", "无 UTR", "CDS 覆盖了全部外显子序列"))

    # 碱基字符
    if uniq_bases:
        checks.append(CheckItem(
            "warn", "存在非 ACGT 字符",
            f"出现：{uniq_bases}。序列按原样保留，请确认是否为简并碱基或软屏蔽区",
        ))
    else:
        checks.append(CheckItem("ok", "碱基字符正常", "仅含 A/C/G/T"))

    # 上下游
    t_start, t_end = tx.span()
    if gene.start and gene.end and (gene.start < t_start or gene.end > t_end):
        checks.append(CheckItem(
            "info", "基因区间大于本转录本",
            f"基因 {gene.start:,}-{gene.end:,}，本转录本 {t_start:,}-{t_end:,}；"
            f"提取与显示均以转录本为准",
        ))
    checks.append(CheckItem(
        "info", "序列来源可追溯",
        f"{getattr(provider, 'path', '')} + 注释 {tx.tx_id}（第 {tx.line_no} 行）",
    ))
    return checks
