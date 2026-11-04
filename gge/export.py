"""导出与复制文本的生成。

所有导出都基于 ``GeneExtraction``，保证"界面上看到的"与"导出的"同源。

FASTA 标题规范
--------------
``>{gene_id}|{tx_id}|{seq_type}|{seqid}:{start}-{end}({strand})|len={n}|{desc}``

* ``seq_type``：genomic / upstream / downstream / exon / intron / mRNA / CDS / UTR
* 拼接序列（mRNA / CDS）没有单一连续区间，使用转录本整体区间并在 ``desc`` 中
  注明 ``spliced`` 与片段数，避免把拼接坐标误当基因组坐标。
* 片段级导出额外带片段序号，例如 ``exon3``、``cds2``，防止不同片段互相覆盖。
"""

from __future__ import annotations

import csv
import io
import os
import re
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .extract import PANEL_ORDER, PANEL_TITLES, SEQ_TYPE_FASTA
from .models import Block, GeneExtraction, Panel

HEADER_SEP = "|"
FASTA_WIDTH = 70


def sanitize(text: str) -> str:
    """FASTA 标题中不允许空白与竖线歧义。"""
    return re.sub(r"\s+", "_", str(text).strip()) or "-"


def _range_text(extraction: GeneExtraction, block: Optional[Block] = None,
                blocks: Optional[Sequence[Block]] = None,
                limit: int = 4) -> str:
    """标题里那段坐标。

    * 单块（或片段级导出）：用它自己的基因组区间，绝不拿基因整体区间顶替——
      上游序列标成基因坐标是错的，会误导别人以为那是基因本体。
    * 多块拼接（完整 mRNA / 完整 CDS）：把所有片段的基因组区间依次列出，
      超过 ``limit`` 个时后面用 ``+N`` 表示数量，保持标题可读。
    * 完全没有基因组区间（理论上不该发生）：退回基因整体区间。
    """
    strand = extraction.strand
    if block is not None:
        gi = block.genomic_interval
        if gi:
            return f"{extraction.seqid}:{gi[0]}-{gi[1]}({strand})"
    if blocks:
        spans = []
        for b in blocks:
            gi = b.genomic_interval
            if gi:
                spans.append(f"{gi[0]}-{gi[1]}")
        if spans:
            if len(spans) > limit:
                shown = ",".join(spans[:limit]) + f",+{len(spans) - limit}"
            else:
                shown = ",".join(spans)
            return f"{extraction.seqid}:{shown}({strand})"
    return f"{extraction.seqid}:{extraction.gene_start}-{extraction.gene_end}({strand})"


def fasta_header(extraction: GeneExtraction, seq_type: str,
                 block: Optional[Block] = None,
                 extra: str = "",
                 blocks: Optional[Sequence[Block]] = None,
                 length: int = 0) -> str:
    parts = [
        sanitize(extraction.gene_id),
        sanitize(extraction.tx_id),
        sanitize(seq_type),
        _range_text(extraction, block, blocks),
    ]
    n = block.length if block is not None else int(length or 0)
    if n:
        parts.append(f"len={n}")
    if extra:
        parts.append(sanitize(extra))
    return ">" + HEADER_SEP.join(parts)


def _wrap(seq: str, width: int = FASTA_WIDTH) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def block_fasta(extraction: GeneExtraction, seq_type: str, block: Block) -> str:
    return f"{fasta_header(extraction, seq_type, block)}\n{_wrap(block.seq)}\n"


def panel_fasta_text(extraction: GeneExtraction, panel_key: str,
                     per_block: bool = True) -> str:
    """把某个面板导出为 FASTA 文本。

    ``per_block=True`` 时每个片段一条记录（外显子/内含子/CDS 片段各自独立）；
    ``False`` 时整个面板拼成一条记录（等价于"完整 mRNA / 完整 CDS"）。
    """
    panel = extraction.panels.get(panel_key)
    if panel is None or panel.is_empty():
        return ""
    seq_type = SEQ_TYPE_FASTA.get(panel_key, panel_key)
    if not per_block or len(panel.blocks) == 1:
        desc = "spliced" if panel_key in ("mrna", "cds") else ""
        if len(panel.blocks) > 1:
            desc = (desc + f" segments={len(panel.blocks)}").strip()
        header = fasta_header(extraction, seq_type, None, extra=desc,
                              blocks=panel.blocks, length=len(panel.sequence))
        return f"{header}\n{_wrap(panel.sequence)}\n"
    out = []
    for i, block in enumerate(panel.blocks, start=1):
        extra = f"segment{i}_of_{len(panel.blocks)}"
        if block.meta.get("exon_index"):
            extra = f"exon{block.meta['exon_index']}"
        elif block.meta.get("intron_index"):
            extra = f"intron{block.meta['intron_index']}"
        elif block.meta.get("cds_index"):
            extra = f"cds{block.meta['cds_index']}"
        header = fasta_header(extraction, seq_type, block, extra=extra)
        out.append(f"{header}\n{_wrap(block.seq)}\n")
    return "".join(out)


def gene_fasta_text(extraction: GeneExtraction,
                    include: Optional[Sequence[str]] = None) -> str:
    """一个基因的全部序列（默认包含全部面板，每个片段单独一条记录）。"""
    keys = include or PANEL_ORDER
    chunks: List[str] = []
    for key in keys:
        text = panel_fasta_text(extraction, key, per_block=True)
        if text:
            chunks.append(text)
    return "".join(chunks)


def gene_basic_info_text(extraction: GeneExtraction, checks: bool = True) -> str:
    """复制"基因基本信息"用。"""
    lines = [
        f"Gene ID: {extraction.gene_id}",
        f"Gene Name: {extraction.gene_name}",
        f"Transcript: {extraction.tx_id}  (该基因共 {extraction.tx_count} 个转录本)",
        f"Location: {extraction.seqid}:{extraction.gene_start:,}-{extraction.gene_end:,} ({extraction.strand})",
        f"Extraction window: {extraction.seqid}:{extraction.region_start:,}-{extraction.region_end:,}",
        f"Upstream: 请求 {extraction.upstream_requested:,} bp，实际 {extraction.upstream_actual:,} bp"
        + ("（已截断到序列边界）" if extraction.upstream_clipped else ""),
        f"Downstream: 请求 {extraction.downstream_requested:,} bp，实际 {extraction.downstream_actual:,} bp"
        + ("（已截断到序列边界）" if extraction.downstream_clipped else ""),
        f"Exons: {extraction.n_exons}   Introns: {extraction.n_introns}   CDS segments: {extraction.n_cds}",
        f"mRNA length: {len(extraction.mrna):,} bp",
        f"CDS length: {len(extraction.cds):,} bp (frame remainder {len(extraction.cds) % 3})",
        f"5'UTR: {len(extraction.utr5):,} bp   3'UTR: {len(extraction.utr3):,} bp"
        + ("   [NOT IN ANNOTATION - computed as exon minus CDS]"
           if not extraction.utr_in_annotation else
           "   [annotation has UTR rows; these values are still computed as exon minus CDS]"),
        f"Genome FASTA: {extraction.fasta_path}",
        f"Annotation: {extraction.annotation_path}",
    ]
    if checks and extraction.checks:
        lines.append("")
        lines.append("Checks:")
        for c in extraction.checks:
            lines.append(f"  [{c.level.upper()}] {c.title}: {c.detail}")
    if extraction.warnings:
        lines.append("")
        lines.append("Warnings:")
        for w in extraction.warnings:
            lines.append(f"  - {w}")
    return "\n".join(lines)


def segments_table(extraction: GeneExtraction) -> List[Dict[str, object]]:
    """把所有片段整理成表格（TSV/CSV/Excel 共用）。"""
    rows: List[Dict[str, object]] = []
    for key in PANEL_ORDER:
        panel = extraction.panels.get(key)
        if not panel:
            continue
        for i, block in enumerate(panel.blocks, start=1):
            gi = block.genomic_interval or (None, None)
            rows.append({
                "panel": key,
                "panel_name": PANEL_TITLES.get(key, key),
                "segment_index": i,
                "segment_title": block.title,
                "seq_type": SEQ_TYPE_FASTA.get(key, key),
                "seqid": extraction.seqid,
                "genomic_start": gi[0] if gi[0] is not None else "",
                "genomic_end": gi[1] if gi[1] is not None else "",
                "strand": extraction.strand,
                "length_bp": block.length,
                "note": block.note,
                "sequence": block.seq,
            })
    return rows


def segments_tsv(extraction: GeneExtraction) -> str:
    rows = segments_table(extraction)
    if not rows:
        return ""
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=list(rows[0].keys()), delimiter="\t",
                            lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buf.getvalue()


def write_text(path: str, text: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text)


def write_csv(path: str, rows: Sequence[Dict[str, object]],
              fields: Optional[Sequence[str]] = None) -> None:
    if not rows:
        write_text(path, "")
        return
    fields = list(fields or rows[0].keys())
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fields})


# --------------------------------------------------------------------------- #
# 校验报告
# --------------------------------------------------------------------------- #
def validation_report_text(extraction: GeneExtraction,
                           results: Dict[str, object],
                           extra_header: str = "") -> str:
    lines: List[str] = []
    lines.append("=" * 78)
    lines.append("序列校验报告 / Sequence validation report")
    lines.append("=" * 78)
    if extra_header:
        lines.append(extra_header)
    lines.append(f"生成时间: {_now()}")
    lines.append(f"基因: {extraction.gene_name} ({extraction.gene_id})")
    lines.append(f"转录本: {extraction.tx_id}")
    lines.append(f"位置: {extraction.seqid}:{extraction.gene_start}-{extraction.gene_end}"
                 f"({extraction.strand})")
    lines.append(f"基因组 FASTA: {extraction.fasta_path}")
    lines.append(f"注释文件: {extraction.annotation_path}")
    lines.append("")

    for kind in ("mRNA", "CDS", "protein"):
        res = results.get(kind)
        if res is None:
            continue
        lines.append("-" * 78)
        lines.append(f"[{kind}] 状态: {res.status_cn}  ({res.status})")
        lines.append(f"  参考 ID      : {res.ref_id or '—'}")
        lines.append(f"  参考长度     : {res.ref_len:,} bp" if res.ref_len else "  参考长度     : —")
        lines.append(f"  提取长度     : {res.obs_len:,}")
        c = res.comparison
        if c:
            lines.append(f"  比对列数     : {c.aligned_columns:,}")
            lines.append(f"  一致 / 不一致: {c.matches:,} / {c.mismatches:,}")
            lines.append(f"  插入 / 缺失  : {c.insertions:,} / {c.deletions:,}  bp")
            lines.append(f"  一致率       : {c.identity:.4f}%")
            if c.events:
                lines.append("  差异事件（最多 50 条）:")
                for ev in c.events[:50]:
                    lines.append(f"    - {ev.describe()}")
                if len(c.events) > 50:
                    lines.append(f"    … 另有 {len(c.events) - 50} 条")
        for level, title, detail in res.coding_checks:
            lines.append(f"  [{level.upper()}] {title}: {detail}")
        for note in res.notes:
            lines.append(f"  * {note}")
        if res.candidates:
            lines.append("  候选参考序列:")
            for cand in res.candidates[:8]:
                lines.append(f"    - {cand.token}  分数 {cand.score:.0f}  {cand.reason_text()}")
        lines.append("")

    if extraction.warnings:
        lines.append("-" * 78)
        lines.append("提取阶段的警告:")
        for w in extraction.warnings:
            lines.append(f"  - {w}")
        lines.append("")
    if extraction.checks:
        lines.append("-" * 78)
        lines.append("提取阶段的自检:")
        for c in extraction.checks:
            lines.append(f"  [{c.level.upper()}] {c.title}: {c.detail}")
        lines.append("")
    lines.append("说明：差异原因仅列出排查方向，不代表已确认的结论；")
    lines.append("      软件不会为了让结果「一致」而修改任何提取到的碱基。")
    lines.append("=" * 78)
    return "\n".join(lines)


def _now() -> str:
    from datetime import datetime

    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def safe_filename(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.\-\u4e00-\u9fff]+", "_", text).strip("_") or "output"
