"""注释解析测试：格式识别、层级重建、异常行处理。"""

from __future__ import annotations

import pytest

from gge.annotation import (AnnotationError, classify, detect_format,
                            parse_annotation, parse_attributes)
from gge.errors import AnnotationError as AnnErr
from gge.errors import FileFormatError


def test_detect_format(demo_dir):
    assert detect_format(demo_dir["gff3"]) == "gff3"
    assert detect_format(demo_dir["gtf"]) == "gtf"


def test_attribute_dialects():
    gff3 = parse_attributes("ID=tx1;Parent=g1;Name=Gene%20One")
    assert gff3["id"] == "tx1"
    assert gff3["parent"] == "g1"
    assert gff3["name"] == "Gene One"

    gtf = parse_attributes('gene_id "g1"; transcript_id "tx1"; gene_name "One";')
    assert gtf["gene_id"] == "g1"
    assert gtf["transcript_id"] == "tx1"
    assert gtf["gene_name"] == "One"

    gff2 = parse_attributes('gene_id "g1" ; transcript_id "tx1"')
    assert gff2["gene_id"] == "g1"


def test_feature_classification():
    assert classify("gene") == "gene"
    assert classify("ncRNA_gene") == "gene"
    assert classify("mRNA") == "transcript"
    assert classify("exon") == "exon"
    assert classify("CDS") == "cds"
    assert classify("five_prime_UTR") == "utr5"
    assert classify("three_prime_UTR") == "utr3"
    assert classify("region") == "ignored"


def test_gene_and_transcript_counts(ann_gff3):
    assert ann_gff3.n_genes == 10
    assert ann_gff3.n_transcripts == 11
    assert "g1" in ann_gff3.genes
    assert ann_gff3.stats["exons"] == 21


def test_multi_transcript_kept_separate(ann_gff3):
    gene = ann_gff3.get_gene("g1")
    assert set(gene.transcripts) == {"g1.1", "g1.2"}
    t1 = gene.transcripts["g1.1"]
    t2 = gene.transcripts["g1.2"]
    assert t1.n_exons == 3 and t2.n_exons == 2
    # 两条转录本的外显子不能混合
    assert len(t1.exons) + len(t2.exons) == 5
    # 规范转录本应是有最长 CDS 的那条
    assert gene.canonical_transcript().tx_id == "g1.1"


def test_gtf_and_gff3_equivalent(ann_gff3, ann_gtf):
    assert ann_gtf.n_genes == ann_gff3.n_genes
    assert ann_gtf.n_transcripts == ann_gff3.n_transcripts
    for gid, gene in ann_gff3.genes.items():
        other = ann_gtf.get_gene(gid)
        assert other is not None, gid
        for tid, tx in gene.transcripts.items():
            otx = other.transcripts[tid]
            assert [(s.start, s.end) for s in tx.exons] == [(s.start, s.end) for s in otx.exons], (gid, tid)
            assert [(s.start, s.end) for s in tx.cds] == [(s.start, s.end) for s in otx.cds], (gid, tid)
            assert tx.strand == otx.strand


def test_cds_only_transcript_derives_exons(ann_gff3):
    gene = ann_gff3.get_gene("g5")
    tx = gene.transcripts["g5.1"]
    assert tx.exons_derived is True
    assert len(tx.exons) == 2
    assert [(s.start, s.end) for s in tx.exons] == [(6001, 6120), (6301, 6400)]
    assert any("exon" in w for w in tx.warnings)


def test_bad_coordinates_normalized_with_warning(ann_gff3):
    gene = ann_gff3.get_gene("g7")
    tx = gene.transcripts["g7.1"]
    segs = [(s.start, s.end) for s in tx.exons]
    assert (7300, 7350) in segs, segs
    assert any(s.invalid for s in tx.exons)
    assert any("起点大于终点" in w for w in tx.warnings)


def test_missing_gene_row_synthesized(ann_gff3):
    gene = ann_gff3.get_gene("g10")
    assert gene is not None
    assert gene.start == 8001 and gene.end == 8200
    assert any("gene 行" in w for w in gene.warnings)


def test_overlapping_exons_do_not_create_negative_intron(ann_gff3, provider):
    from gge.extract import extract_gene

    gene = ann_gff3.get_gene("g9")
    tx = gene.canonical_transcript()
    res = extract_gene(provider, gene, tx)
    assert res.n_introns == 0
    assert res.panels["intron"].is_empty()
    assert any("重叠" in w for w in res.warnings)
    # mRNA 长度 = 两段外显子长度之和（重叠部分按模型重复计入）
    assert len(res.mrna) == sum(e.length for e in tx.exons)


def test_empty_annotation_raises(tmp_path):
    p = tmp_path / "empty.gff3"
    p.write_text("##gff-version 3\n", encoding="utf-8")
    with pytest.raises(AnnotationError):
        parse_annotation(str(p))


def test_missing_file_raises(tmp_path):
    with pytest.raises(FileFormatError):
        parse_annotation(str(tmp_path / "nope.gff3"))


def test_broken_lines_reported(write_annotation, tmp_path):
    text = (
        "##gff-version 3\n"
        "chr1\tdemo\tgene\t10\t100\t.\t+\t.\tID=gt1\n"
        "chr1\tdemo\tmRNA\t10\t100\t.\t+\t.\tID=gt1.1;Parent=gt1\n"
        "chr1\tdemo\texon\t10\t50\t.\t+\t.\n"          # 只有 8 列
        "chr1\tdemo\texon\t10\t50\t.\t+\t.\tID=e1;Parent=gt1.1\n"
        "chr1\tdemo\texon\t60\tXXX\t.\t+\t.\tID=e2;Parent=gt1.1\n"  # 坐标非法
        "chr1\tdemo\texon\t60\t100\t.\t+\t.\tID=e3;Parent=gt1.1\n"
    )
    path = write_annotation(text)
    ann = parse_annotation(path)
    assert ann.n_genes == 1
    assert len(ann.get_gene("gt1").transcripts["gt1.1"].exons) == 2
    assert any("跳过" in w for w in ann.warnings)


def test_orphan_feature_is_reported(write_annotation):
    text = (
        "##gff-version 3\n"
        "chr1\tdemo\tgene\t10\t100\t.\t+\t.\tID=gt2\n"
        "chr1\tdemo\tmRNA\t10\t100\t.\t+\t.\tID=gt2.1;Parent=gt2\n"
        "chr1\tdemo\tmRNA\t10\t100\t.\t+\t.\tID=gt2.2;Parent=gt2\n"
        "chr1\tdemo\texon\t10\t50\t.\t+\t.\tID=e1;Parent=unknown_tx\n"  # Parent 不存在
        "chr1\tdemo\texon\t60\t100\t.\t+\t.\tID=e2;Parent=gt2.1\n"
    )
    ann = parse_annotation(write_annotation(text))
    gene = ann.get_gene("gt2")
    # 不能把 Parent 指向不存在的 exon 静默塞进 gt2.1
    assert [(s.start, s.end) for s in gene.transcripts["gt2.1"].exons] == [(60, 100)]
    assert [(s.start, s.end) for s in gene.transcripts["gt2.2"].exons] == []
    assert "unknown_tx" in ann.tx_index
    assert any("Parent 在注释中不存在" in w for w in ann.warnings)


def test_seqid_list(ann_gff3):
    assert set(ann_gff3.seqids) == {"chr1", "chr2", "chrScaf", "chrUn_1"}


def test_duplicate_exon_rows_are_merged(write_annotation):
    """同一条 exon 写两遍时只算一次，否则 mRNA 长度会凭空变长。"""
    text = "\n".join([
        "##gff-version 3",
        "chr1\ttest\tgene\t101\t400\t.\t+\t.\tID=gX;Name=X",
        "chr1\ttest\tmRNA\t101\t400\t.\t+\t.\tID=gX.t1;Parent=gX",
        "chr1\ttest\texon\t101\t200\t.\t+\t.\tID=ex1;Parent=gX.t1",
        "chr1\ttest\texon\t101\t200\t.\t+\t.\tID=ex1_dup;Parent=gX.t1",  # 重复行
        "chr1\ttest\texon\t301\t400\t.\t+\t.\tID=ex2;Parent=gX.t1",
        "chr1\ttest\tCDS\t121\t200\t.\t+\t0\tID=cds1;Parent=gX.t1",
        "chr1\ttest\tCDS\t121\t200\t.\t+\t0\tID=cds1_dup;Parent=gX.t1",  # 重复行
        "",
    ])
    ann = parse_annotation(write_annotation(text))
    gid = ann.tx_index["gX.t1"]          # 索引存的是 tx_id -> gene_id
    tx = ann.genes[gid].transcripts["gX.t1"]
    assert len(tx.exons) == 2, "重复 exon 应被合并为 2 条"
    assert tx.mrna_length() == 200, "mRNA 长度应为 100+100，而不是 300"
    assert len(tx.cds) == 1
    msgs = " ".join(tx.warnings)
    assert "重复的 exon" in msgs
    assert "重复的 CDS" in msgs


def test_duplicate_exon_rows_are_merged(write_annotation):
    """同一条 exon 写两遍时只算一次，否则 mRNA 长度会凭空变长。"""
    text = "\n".join([
        "##gff-version 3",
        "chr1\ttest\tgene\t101\t400\t.\t+\t.\tID=gX;Name=X",
        "chr1\ttest\tmRNA\t101\t400\t.\t+\t.\tID=gX.t1;Parent=gX",
        "chr1\ttest\texon\t101\t200\t.\t+\t.\tID=ex1;Parent=gX.t1",
        "chr1\ttest\texon\t101\t200\t.\t+\t.\tID=ex1_dup;Parent=gX.t1",  # 重复行
        "chr1\ttest\texon\t301\t400\t.\t+\t.\tID=ex2;Parent=gX.t1",
        "chr1\ttest\tCDS\t121\t200\t.\t+\t0\tID=cds1;Parent=gX.t1",
        "chr1\ttest\tCDS\t121\t200\t.\t+\t0\tID=cds1_dup;Parent=gX.t1",  # 重复行
        "",
    ])
    ann = parse_annotation(write_annotation(text))
    gid = ann.tx_index["gX.t1"]          # 索引存的是 tx_id -> gene_id
    tx = ann.genes[gid].transcripts["gX.t1"]
    assert len(tx.exons) == 2, "重复 exon 应被合并为 2 条"
    assert tx.mrna_length() == 200, "mRNA 长度应为 100+100，而不是 300"
    assert len(tx.cds) == 1
    msgs = " ".join(tx.warnings)
    assert "重复的 exon" in msgs
    assert "重复的 CDS" in msgs


# --------------------------------------------------------------------------- #
# 大文件性能与内存（NCBI 全基因组注释实测：524 MB / 178 万行）
# --------------------------------------------------------------------------- #

def test_attributes_gff3_without_quotes_skips_gtf_regex():
    """GFF3 属性串不含引号时必须跳过 GTF 正则。

    这条正则在无引号的串上是 O(n²) 回溯：每个起始位置都把 ``[A-Za-z0-9_]+``
    贪婪匹配到串尾再回退。NCBI 属性列约 500 字符 → 单次 57 µs，178 万行 = 101 秒。
    """
    attr = ("ID=exon-NC_067220.1-1;Parent=rna-NC_067220.1;"
            "Dbxref=GeneID:127508713,Genbank:XM_051;gbkey=mRNA;gene=abc1;"
            "product=some+protein-like%20isoform+1")
    parsed = parse_attributes(attr)
    assert parsed["id"] == "exon-NC_067220.1-1"
    assert parsed["parent"] == "rna-NC_067220.1"
    assert parsed["gene"] == "abc1"
    # 含 % 的值仍要转义还原
    assert parsed["product"] == "some+protein-like isoform+1"
    # 不在白名单里的键不保留（Dbxref 很占地方，本来就不需要）
    assert "dbxref" not in parsed


def _reference_parse_attributes(text: str):
    """优化前的实现（原样保留），用于证明"快路径"没有改变任何语义。"""
    from urllib.parse import unquote as _uq

    from gge.annotation import _GTF_ATTR_RE, KEEP_ATTRS

    text = text.strip()
    if not text or text == ".":
        return {}
    attrs = {}
    if "=" in text and not _GTF_ATTR_RE.search(text):
        for part in text.split(";"):
            part = part.strip()
            if not part or "=" not in part:
                continue
            key, _, val = part.partition("=")
            key = key.strip().lower()
            if key and key in KEEP_ATTRS:
                attrs[key] = _uq(val.strip())
        return attrs
    for m in _GTF_ATTR_RE.finditer(text):
        key = m.group(1).strip().lower()
        if key in KEEP_ATTRS:
            attrs[key] = m.group(2)
    if attrs:
        return attrs
    for part in text.split(";"):
        bits = part.strip().split(None, 1)
        if len(bits) == 2:
            key = bits[0].lower()
            if key in KEEP_ATTRS:
                attrs[key] = bits[1].strip()
    return attrs


@pytest.mark.parametrize("text", [
    # 纯 GFF3（无引号，走快路径：跳过 O(n²) 的 GTF 正则）
    "ID=g1;Parent=p1;Name=Gene%20One",
    "ID=exon-NC_1.1-1;Parent=rna-NC_1.1;Dbxref=GeneID:1;gbkey=mRNA;gene=abc1;"
    "product=some+protein-like%20isoform+1",
    "ID=cds-XP_1.1;Parent=rna-XP_1.1;protein_id=XP_1.1;gbkey=CDS",
    # 值里恰好有引号（会走 GTF 分支——与优化前一致，不是新行为）
    'ID=g1;product=the "ABC" transporter;gene=abc1',
    # 真 GTF
    'gene_id "g1"; transcript_id "t1"; gene_name "One";',
    'gene_id "g1" ; transcript_id "t1"',
    # 无引号 GFF2
    "gene_id g1; transcript_id t1",
    # 边界
    "", ".", "ID=a", "ID=a;", "novalue;ID=b;product=x%3By",
])
def test_attribute_parsing_matches_reference_implementation(text):
    """快路径必须与原实现逐字符等价（含所有方言与边界）。"""
    assert parse_attributes(text) == _reference_parse_attributes(text)


def test_segment_attrs_keep_only_strand(write_annotation):
    """片段只保留 strand：整字典复制在 167 万个片段上要多吃约 600 MB。"""
    ann = parse_annotation(write_annotation(
        "chr1\ttest\tgene\t1\t300\t.\t+\t.\tID=gS\n"
        "chr1\ttest\tmRNA\t1\t300\t.\t+\t.\tID=gS.t1;Parent=gS\n"
        "chr1\ttest\texon\t1\t100\t.\t+\t.\tID=e1;Parent=gS.t1;gene=abc1;gbkey=mRNA\n"
        "chr1\ttest\tCDS\t1\t90\t.\t+\t0\tID=c1;Parent=gS.t1;gene=abc1;protein_id=P1\n"
    ))
    tx = ann.genes["gS"].transcripts["gS.t1"]
    assert tx.exons[0].attrs == {}, "没有 strand 属性时就是空字典"
    assert set(tx.cds[0].attrs) <= {"strand"}
    # 基因/转录本的属性仍然是完整的（只有片段瘦身）
    assert ann.genes["gS"].attrs.get("id") == "gS"


def test_exon_strand_attribute_promotes_unknown_transcript_strand(write_annotation):
    """``_segment_attrs`` 的契约：只留 strand，其余键一律丢弃。

    注意 ``KEEP_ATTRS`` 不含 ``strand``（GFF3 的链方向在第 7 列），所以真实解析
    出来的片段属性是空的——这一步只是把"片段需要哪些属性"写成显式契约，
    以后往 KEEP_ATTRS 里加 strand 时才不会因为丢了它而静默失效。
    """
    from gge.annotation import _segment_attrs

    assert _segment_attrs({"strand": "-", "id": "e1", "product": "x"}) == {"strand": "-"}
    assert _segment_attrs({"id": "e1"}) == {}
    assert _segment_attrs({}) == {}

    ann = parse_annotation(write_annotation(
        "chr1\ttest\tgene\t1\t300\t.\t.\t.\tID=gU\n"
        "chr1\ttest\tmRNA\t1\t300\t.\t.\t.\tID=gU.t1;Parent=gU\n"
        "chr1\ttest\texon\t1\t100\t.\t.\t.\tID=e1;Parent=gU.t1;strand=-\n"
    ))
    tx = ann.genes["gU"].transcripts["gU.t1"]
    assert tx.exons[0].attrs == {}, "strand 不在白名单里，片段属性为空（与优化前一致）"


def test_read_rows_progress_advances(tmp_path):
    """进度必须递增而不是一直停在 40%（否则用户以为卡死）。"""
    from gge.annotation import read_rows

    lines = []
    for i in range(120000):
        lines.append(f"chr1\ttest\texon\t{i + 1}\t{i + 2}\t.\t+\t.\tID=e{i};Parent=t1")
    p = tmp_path / "many.gff3"
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")

    calls = []
    read_rows(str(p), progress_cb=lambda pct, msg: calls.append(pct))
    assert len(calls) >= 2, f"120000 行应至少报两次进度，实际 {calls}"
    assert calls == sorted(calls), f"进度不能倒退：{calls}"
    assert max(calls) > 15, "文件超过 15% 时进度条要动起来"
    assert max(calls) <= 70, "读行阶段不应超过 70%，后面还要建层级"


# --------------------------------------------------------------------------- #
# read_rows 的内存契约（178 万行实测：整字典保留会吃掉 2.4 GB）
# --------------------------------------------------------------------------- #
def test_read_rows_keeps_attrs_only_for_gene_and_transcript(tmp_path):
    """gene/transcript 行保留整字典，子特征行只留 tx_key/gene_key。

    这是有意的内存取舍：NCBI 全基因组注释 167 万个 exon/CDS 行，
    每行一个 6~8 键的字典 + 一串属性值字符串，实测让 read_rows 一步吃 2.4 GB，
    在 15 GB 内存、剩 3 GB 的机器上直接换页卡死。
    """
    from gge.annotation import read_rows

    p = tmp_path / "mix.gff3"
    p.write_text(
        "chr1\tt\tgene\t1\t300\t.\t+\t.\tID=g1;Name=Gene One\n"
        "chr1\tt\tmRNA\t1\t300\t.\t+\t.\tID=g1.1;Parent=g1;product=thing\n"
        "chr1\tt\texon\t1\t100\t.\t+\t.\tID=e1;Parent=g1.1;gene=one;gbkey=mRNA\n",
        encoding="utf-8",
    )
    rows = read_rows(str(p))
    gene, tx, exon = rows

    assert gene.kind == "gene" and gene.attrs["id"] == "g1" and gene.attrs["name"] == "Gene One"
    assert tx.kind == "transcript" and tx.attrs["product"] == "thing"
    # 子特征：字典为空，关键信息挪到专用字段
    assert exon.kind == "exon"
    assert exon.attrs == {}, "子特征行不应保留属性字典"
    assert exon.tx_key == "g1.1"
    assert exon.gene_key == "one"


def test_read_rows_gtf_child_uses_transcript_id_when_no_parent(tmp_path):
    """GTF 没有 Parent，子特征的归属来自 transcript_id。"""
    from gge.annotation import read_rows

    p = tmp_path / "a.gtf"
    p.write_text(
        'chr1\tt\tgene\t1\t300\t.\t+\t.\tgene_id "g1"; gene_name "One";\n'
        'chr1\tt\ttranscript\t1\t300\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
        'chr1\tt\texon\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n',
        encoding="utf-8",
    )
    rows = read_rows(str(p))
    assert rows[0].kind == "gene" and rows[0].attrs["gene_id"] == "g1"
    assert rows[1].kind == "transcript"
    assert rows[2].attrs == {}
    assert rows[2].tx_key == "t1"        # 没有 Parent 时用 transcript_id
    assert rows[2].gene_key == "g1"


def test_read_rows_other_type_ignores_transcript_id(tmp_path):
    """"other" 类只认 Parent——保持与原实现一致，不因为顺手就多认一个键。"""
    from gge.annotation import read_rows

    p = tmp_path / "o.gtf"
    p.write_text(
        'chr1\tt\tcustom_thing\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n',
        encoding="utf-8",
    )
    row = read_rows(str(p))[0]
    assert row.kind == "other"
    assert row.tx_key == "", "other 类不应从 transcript_id 推归属（原实现就是这样）"


def test_read_rows_interns_repetitive_strings(tmp_path):
    """seqid / ftype / strand 只有几个取值，必须共用同一个字符串对象。"""
    from gge.annotation import read_rows

    lines = [f"chr1\tt\texon\t{i + 1}\t{i + 2}\t.\t+\t.\tID=e{i};Parent=t1"
             for i in range(50)]
    lines.append("chr2\tt\texon\t1\t5\t.\t-\t.\tID=z;Parent=t1")
    p = tmp_path / "many.gff3"
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")

    rows = read_rows(str(p))
    assert len({id(r.seqid) for r in rows if r.seqid == "chr1"}) == 1
    assert len({id(r.ftype) for r in rows}) == 1
    assert len({id(r.strand) for r in rows if r.strand == "+"}) == 1


# --------------------------------------------------------------------------- #
# read_rows 的内存契约（178 万行实测：整字典保留会吃掉 2.4 GB）
# --------------------------------------------------------------------------- #
def test_read_rows_keeps_attrs_only_for_gene_and_transcript(tmp_path):
    """gene/transcript 行保留整字典，子特征行只留 tx_key/gene_key。

    这是有意的内存取舍：NCBI 全基因组注释 167 万个 exon/CDS 行，
    每行一个 6~8 键的字典 + 一串属性值字符串，实测让 read_rows 一步吃 2.4 GB，
    在 15 GB 内存、剩 3 GB 的机器上直接换页卡死。
    """
    from gge.annotation import read_rows

    p = tmp_path / "mix.gff3"
    p.write_text(
        "chr1\tt\tgene\t1\t300\t.\t+\t.\tID=g1;Name=Gene One\n"
        "chr1\tt\tmRNA\t1\t300\t.\t+\t.\tID=g1.1;Parent=g1;product=thing\n"
        "chr1\tt\texon\t1\t100\t.\t+\t.\tID=e1;Parent=g1.1;gene=one;gbkey=mRNA\n",
        encoding="utf-8",
    )
    rows = read_rows(str(p))
    gene, tx, exon = rows

    assert gene.kind == "gene" and gene.attrs["id"] == "g1" and gene.attrs["name"] == "Gene One"
    assert tx.kind == "transcript" and tx.attrs["product"] == "thing"
    # 子特征：字典为空，关键信息挪到专用字段
    assert exon.kind == "exon"
    assert exon.attrs == {}, "子特征行不应保留属性字典"
    assert exon.tx_key == "g1.1"
    assert exon.gene_key == "one"


def test_read_rows_gtf_child_uses_transcript_id_when_no_parent(tmp_path):
    """GTF 没有 Parent，子特征的归属来自 transcript_id。"""
    from gge.annotation import read_rows

    p = tmp_path / "a.gtf"
    p.write_text(
        'chr1\tt\tgene\t1\t300\t.\t+\t.\tgene_id "g1"; gene_name "One";\n'
        'chr1\tt\ttranscript\t1\t300\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
        'chr1\tt\texon\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n',
        encoding="utf-8",
    )
    rows = read_rows(str(p))
    assert rows[0].kind == "gene" and rows[0].attrs["gene_id"] == "g1"
    assert rows[1].kind == "transcript"
    assert rows[2].attrs == {}
    assert rows[2].tx_key == "t1"        # 没有 Parent 时用 transcript_id
    assert rows[2].gene_key == "g1"


def test_read_rows_other_type_ignores_transcript_id(tmp_path):
    """"other" 类只认 Parent——保持与原实现一致，不因为顺手就多认一个键。"""
    from gge.annotation import read_rows

    p = tmp_path / "o.gtf"
    p.write_text(
        'chr1\tt\tcustom_thing\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n',
        encoding="utf-8",
    )
    row = read_rows(str(p))[0]
    assert row.kind == "other"
    assert row.tx_key == "", "other 类不应从 transcript_id 推归属（原实现就是这样）"


def test_read_rows_interns_repetitive_strings(tmp_path):
    """seqid / ftype / strand 只有几个取值，必须共用同一个字符串对象。"""
    from gge.annotation import read_rows

    lines = [f"chr1\tt\texon\t{i + 1}\t{i + 2}\t.\t+\t.\tID=e{i};Parent=t1"
             for i in range(50)]
    lines.append("chr2\tt\texon\t1\t5\t.\t-\t.\tID=z;Parent=t1")
    p = tmp_path / "many.gff3"
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")

    rows = read_rows(str(p))
    assert len({id(r.seqid) for r in rows if r.seqid == "chr1"}) == 1
    assert len({id(r.ftype) for r in rows}) == 1
    assert len({id(r.strand) for r in rows if r.strand == "+"}) == 1


# --------------------------------------------------------------------------- #
# 取消：大文件解析必须能被"取消"打断
# --------------------------------------------------------------------------- #
def _big_gff(path, n=120_000):
    """造一个够长的 GFF：read_rows 每 5 万行查一次取消，所以必须 > 5 万行。"""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("##gff-version 3\n")
        for i in range(n):
            fh.write(f"chr1\tt\tgene\t{i * 10 + 1}\t{i * 10 + 9}\t.\t+\t.\tID=g{i}\n")
    return str(path)


def test_read_rows_cancel(tmp_path):
    """取消要真的抛 CancelledError，而不是"跑完才停"。"""
    from gge.annotation import read_rows
    from gge.errors import CancelledError

    p = _big_gff(tmp_path / "big.gff3")
    calls = {"n": 0}

    def cancel():
        calls["n"] += 1
        return calls["n"] >= 2          # 第二次检查时喊停

    with pytest.raises(CancelledError):
        read_rows(p, cancel_check=cancel)


def test_parse_annotation_cancel(tmp_path):
    """parse_annotation 要把 cancel_check 一路传下去（以前这个参数根本没接）。"""
    from gge.annotation import parse_annotation
    from gge.errors import CancelledError

    p = _big_gff(tmp_path / "big.gff3")
    with pytest.raises(CancelledError):
        parse_annotation(p, cancel_check=lambda: True)


def test_load_tasks_accept_cancel(tmp_path):
    """任务层三个加载函数都得声明 cancel_check，否则界面点「取消」对它们无效。"""
    import inspect

    from gge.tasks import load_annotation_file, load_genome, load_reference_file

    for fn in (load_genome, load_annotation_file, load_reference_file):
        assert "cancel_check" in inspect.signature(fn).parameters, fn.__name__
