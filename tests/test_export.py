"""导出与复制文本测试。"""

from __future__ import annotations

import csv
import io
import os

from gge.export import (gene_basic_info_text, gene_fasta_text, panel_fasta_text,
                        segments_table, segments_tsv, validation_report_text)
from gge.extract import PANEL_CDS, PANEL_EXON, PANEL_MRNA, PANEL_UTR

from conftest import extract_one


def _res(provider, ann):
    return extract_one(provider, ann, "g1", "g1.1", up=500, dn=500)


def test_panel_fasta_header_and_wrap(provider, ann_gff3):
    res = _res(provider, ann_gff3)
    text = panel_fasta_text(res, PANEL_EXON, per_block=True)
    lines = text.strip().split("\n")
    header = lines[0]
    assert header.startswith(">g1|g1.1|exon|chr1:")
    assert "(+)" in header and "exon1" in header
    # 每条记录的序列行不超过 70 字符
    for line in lines[1:]:
        assert len(line) <= 70
    assert text.count(">") == 3, "3 个外显子应有 3 条 FASTA 记录"


def test_panel_fasta_joined(provider, ann_gff3):
    res = _res(provider, ann_gff3)
    text = panel_fasta_text(res, PANEL_MRNA, per_block=False)
    seq = "".join(text.strip().split("\n")[1:])
    assert seq == res.mrna
    assert "mRNA" in text.split("\n")[0]
    assert "spliced" in text.split("\n")[0]


def test_gene_fasta_covers_all_panels(provider, ann_gff3):
    res = _res(provider, ann_gff3)
    text = gene_fasta_text(res)
    types = {h.split("|")[2] for h in text.split("\n") if h.startswith(">")}
    assert {"genomic", "upstream", "downstream", "exon", "intron", "mRNA",
            "CDS", "UTR"} <= types
    # 标题里必须能区分不同片段
    assert ">g1|g1.1|intron|" in text
    assert ">g1|g1.1|CDS|" in text


def test_exported_sequence_has_no_decorations(provider, ann_gff3):
    res = _res(provider, ann_gff3)
    text = panel_fasta_text(res, PANEL_CDS, per_block=False)
    body = "".join(line for line in text.strip().split("\n")[1:])
    assert set(body) <= set("ACGTNacgtn")
    assert body == res.cds


def test_basic_info_text(provider, ann_gff3):
    res = _res(provider, ann_gff3)
    res.annotation_path = "demo.gff3"
    info = gene_basic_info_text(res)
    assert "Gene ID: g1" in info
    assert "Transcript: g1.1" in info
    assert "Exons: 3" in info
    assert "CDS length: 360" in info
    assert "demo.gff3" in info


def test_segments_table(provider, ann_gff3):
    res = _res(provider, ann_gff3)
    rows = segments_table(res)
    assert len(rows) == sum(len(p.blocks) for p in res.panels.values())
    first = rows[0]
    for key in ("panel", "seq_type", "genomic_start", "genomic_end", "length_bp", "sequence"):
        assert key in first
    tsv = segments_tsv(res)
    reader = csv.DictReader(io.StringIO(tsv), delimiter="\t")
    parsed = list(reader)
    assert len(parsed) == len(rows)
    assert parsed[0]["length_bp"] == str(len(parsed[0]["sequence"]))


def test_validation_report_text(provider, ann_gff3, mrna_refs, cds_refs):
    from gge.validate import validate_extraction

    gene = ann_gff3.get_gene("g1")
    tx = gene.transcripts["g1.1"]
    res = extract_one(provider, ann_gff3, "g1", "g1.1")
    res.annotation_path = "demo_annotation.gff3"
    v = validate_extraction(res, gene, tx, mrna_refs, cds_refs)
    text = validation_report_text(res, v)
    assert "序列校验报告" in text
    assert "完全一致" in text
    assert "一致率" in text
    assert "g1.1" in text


def test_write_files(provider, ann_gff3, tmp_path):
    from gge.export import write_csv, write_text

    res = _res(provider, ann_gff3)
    p = tmp_path / "sub" / "gene.fa"
    write_text(str(p), gene_fasta_text(res))
    assert p.exists() and p.stat().st_size > 100
    c = tmp_path / "seg.csv"
    write_csv(str(c), segments_table(res))
    assert c.exists()
    with open(c, encoding="utf-8-sig") as fh:
        assert fh.readline().startswith("panel,")


def test_header_uses_own_coordinates_not_gene_span(provider, ann_gff3):
    """上游序列的标题必须是它自己的区间，不能拿基因整体区间顶替。"""
    res = _res(provider, ann_gff3)
    up = panel_fasta_text(res, "upstream", per_block=False).split("\n")[0]
    assert "chr1:501-1000(+)" in up, up
    assert "1001-1800" not in up, "上游标题不应出现基因本体区间"
    down = panel_fasta_text(res, "downstream", per_block=False).split("\n")[0]
    assert "chr1:1801-2300(+)" in down, down


def test_header_lists_all_segments_for_spliced(provider, ann_gff3):
    """拼接序列必须列出各片段区间并带 len，避免把拼接坐标当基因组坐标。"""
    res = _res(provider, ann_gff3)
    head = panel_fasta_text(res, PANEL_MRNA, per_block=False).split("\n")[0]
    for span in ("1001-1200", "1401-1550", "1701-1800"):
        assert span in head, head
    assert "len=450" in head
    assert "spliced" in head


def test_header_len_matches_sequence_length(provider, ann_gff3):
    res = _res(provider, ann_gff3)
    for key in ("genomic", "upstream", "downstream", "exon", "intron", "mRNA",
                "CDS", "UTR"):
        text = panel_fasta_text(res, key, per_block=False)
        if not text:
            continue
        head, *body = text.strip().split("\n")
        n = len("".join(body))
        assert f"len={n}" in head, (key, head, n)
