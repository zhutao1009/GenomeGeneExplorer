"""序列提取测试：正负链、单/多外显子、边界、UTR、CDS、面板一致性。

断言的是**真实碱基**，不是"函数有没有返回"。
期望值由 ``examples/make_example_data.py`` 用最直白的切片方式独立算出。
"""

from __future__ import annotations

import pytest

from gge.coords import revcomp
from gge.errors import SeqIDMismatchError
from gge.extract import (PANEL_CDS, PANEL_DOWNSTREAM, PANEL_EXON, PANEL_GENOMIC,
                         PANEL_INTRON, PANEL_MRNA, PANEL_UPSTREAM, PANEL_UTR,
                         extract_gene)
from gge.models import (REGION_CDS, REGION_DOWNSTREAM, REGION_EXON,
                        REGION_UPSTREAM, REGION_UTR3, REGION_UTR5)

from conftest import extract_one


# --------------------------------------------------------------------------- #
# 正链
# --------------------------------------------------------------------------- #
def test_plus_strand_exon_and_mrna(provider, ann_gff3, demo_dir):
    res = extract_one(provider, ann_gff3, "g1", "g1.1")
    exp = demo_dir["expected"]["g1"]["g1.1"]
    assert res.strand == "+"
    assert res.n_exons == 3 and res.n_introns == 2
    assert [b.seq for b in res.panels[PANEL_EXON].blocks] == exp["exon_seqs"]
    assert res.mrna == exp["mrna"]
    assert len(res.mrna) == 450
    assert res.cds == exp["cds"]
    assert len(res.cds) == 360


def test_plus_strand_utr(provider, ann_gff3, demo_dir):
    res = extract_one(provider, ann_gff3, "g1", "g1.1")
    genome = demo_dir["genome_sequences"]["chr1"]
    assert res.utr5 == genome[1000:1050]        # 1001..1050
    assert res.utr3 == genome[1760:1800]        # 1761..1800
    assert len(res.utr5) == 50 and len(res.utr3) == 40
    assert res.mrna == res.utr5 + res.cds + res.utr3


@pytest.mark.parametrize("gid,tid", [("g1", "g1.1"), ("g2", "g2.1"),
                                     ("g3", "g3.1"), ("g4", "g4.1")])
def test_utr_blocks_are_exon_minus_cds(provider, ann_gff3, gid, tid):
    """界面上的 UTR 是推导值：注释里**没有** UTR 行，唯一依据是「外显子 − CDS」。

    对每一块断言三件事：
    ① 整块落在某个外显子内（UTR 不可能跑到内含子或基因外）；
    ② 与任何 CDS 片段都不重叠（否则就是 exon 而不是 UTR）；
    ③ 序列能由基因组坐标独立复算（正链直取、负链反向互补）。
    """
    res = extract_one(provider, ann_gff3, gid, tid)
    tx = ann_gff3.genes[gid].transcripts[tid]
    exons, cds = tx.ordered_exons(), tx.ordered_cds()

    for b in res.panels[PANEL_UTR].blocks:
        gs, ge = int(b.meta["genomic_start"]), int(b.meta["genomic_end"])
        assert any(e.start <= gs and ge <= e.end for e in exons), f"{gs}-{ge} 不在任何外显子内"
        assert not any(gs <= c.end and ge >= c.start for c in cds), f"{gs}-{ge} 与 CDS 重叠"
        direct = provider.fetch(res.seqid, gs, ge)
        if tx.strand == "-":
            direct = revcomp(direct)
        assert direct == b.seq


@pytest.mark.parametrize("gid,tid", [("g1", "g1.1"), ("g3", "g3.1")])
def test_utr_length_identity(provider, ann_gff3, gid, tid):
    """注释自洽时：5′UTR + CDS + 3′UTR = mRNA，且 UTR 总量 = 外显子总长 − CDS 长度。"""
    res = extract_one(provider, ann_gff3, gid, tid)
    tx = ann_gff3.genes[gid].transcripts[tid]
    assert len(res.mrna) == len(res.utr5) + len(res.cds) + len(res.utr3)
    assert (len(res.utr5) + len(res.utr3)
            == sum(e.length for e in tx.ordered_exons()) - len(res.cds))


def test_internal_non_cds_stretch_is_not_utr(provider, ann_gff3):
    """g2 的 CDS 在转录本上不连续（150–199、400–496 是非编码空隙）。

    这些空隙位于「第一个 CDS 起点」与「最后一个 CDS 终点」之间，
    属于注释异常，**不能**当作 UTR——UTR 只取 CDS 跨度两端之外的部分。
    """
    res = extract_one(provider, ann_gff3, "g2", "g2.1")
    assert res.panels[PANEL_UTR].is_empty()
    assert res.utr5 == "" and res.utr3 == ""
    assert len(res.mrna) == sum(e.length for e in
                                ann_gff3.genes["g2"].transcripts["g2.1"].ordered_exons())


def test_utr_absence_in_annotation_is_declared(provider, ann_gff3):
    """样例注释没有 UTR 行时，界面与导出必须明说"这是推导值"。

    RefSeq 的 GFF 只给 exon/CDS，用户看到 UTR 面板第一反应就是"注释里没有啊"，
    所以标题、徽标、每段序列说明里都要写清楚来源与算法。
    """
    res = extract_one(provider, ann_gff3, "g1", "g1.1")
    panel = res.panels[PANEL_UTR]
    assert res.utr_in_annotation is False
    assert "推导" in panel.title
    assert "注释无 UTR" in panel.badge
    assert "不存在 UTR" in panel.note
    assert "首个 CDS 起点" in panel.note and "末个 CDS 终点" in panel.note
    for b in panel.blocks:
        assert "（推导）" in b.title
        assert "推导" in b.note
        assert ("CDS 起点" in b.note) or ("CDS 终点" in b.note)


def test_utr_presence_in_annotation_is_declared(provider, demo_dir, tmp_path):
    """注释里**有** UTR 行时也要说明：本程序仍按「外显子 − CDS」推导，不直接采用。"""
    from gge.annotation import parse_annotation

    gff = tmp_path / "with_utr.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\tdemo\tgene\t1001\t1800\t.\t+\t.\tID=g1;Name=GeneOne\n"
        "chr1\tdemo\tmRNA\t1001\t1800\t.\t+\t.\tID=g1.1;Parent=g1\n"
        "chr1\tdemo\texon\t1001\t1200\t.\t+\t.\tID=e1;Parent=g1.1\n"
        "chr1\tdemo\texon\t1401\t1550\t.\t+\t.\tID=e2;Parent=g1.1\n"
        "chr1\tdemo\texon\t1701\t1800\t.\t+\t.\tID=e3;Parent=g1.1\n"
        "chr1\tdemo\tCDS\t1051\t1200\t.\t+\t0\tID=c1;Parent=g1.1\n"
        "chr1\tdemo\tCDS\t1401\t1550\t.\t+\t0\tID=c2;Parent=g1.1\n"
        "chr1\tdemo\tCDS\t1701\t1760\t.\t+\t0\tID=c3;Parent=g1.1\n"
        "chr1\tdemo\tfive_prime_UTR\t1001\t1050\t.\t+\t.\tID=u5;Parent=g1.1\n"
        "chr1\tdemo\tthree_prime_UTR\t1761\t1800\t.\t+\t.\tID=u3;Parent=g1.1\n",
        encoding="utf-8")

    ann = parse_annotation(str(gff))
    res = extract_one(provider, ann, "g1", "g1.1")
    assert res.utr_in_annotation is True
    panel = res.panels[PANEL_UTR]
    assert "注释有 UTR" in panel.badge
    assert "未使用它们" in panel.note
    # 推导结果与注释一致（这份构造数据里两者本就该相同）
    assert res.utr5 == "".join(b.seq for b in panel.blocks
                               if b.meta["utr_kind"] == REGION_UTR5)
    assert len(res.utr5) == 50 and len(res.utr3) == 40



def test_cds_is_not_whole_exon(provider, ann_gff3, demo_dir):
    res = extract_one(provider, ann_gff3, "g1", "g1.1")
    assert len(res.cds) < len(res.mrna)
    assert res.cds != res.mrna
    assert res.cds.startswith("ATG")
    assert res.cds[-3:] == "TAA"


def test_upstream_downstream_plus(provider, ann_gff3, demo_dir):
    res = extract_one(provider, ann_gff3, "g1", "g1.1", up=2000, dn=2000)
    genome = demo_dir["genome_sequences"]["chr1"]
    assert res.upstream_seq == genome[0:1000]          # 只能取到 1000 bp
    assert res.downstream_seq == genome[1800:3800]
    assert res.upstream_requested == 2000
    assert res.upstream_actual == 1000
    assert res.upstream_clipped is True
    assert res.downstream_clipped is False
    assert any("上游请求" in w for w in res.warnings)


def test_plus_strand_genomic_window_covers_flanks(provider, ann_gff3):
    """基因组连续窗口必须横跨「上游 + 基因 + 下游」，正链。"""
    res = extract_one(provider, ann_gff3, "g1", "g1.1", up=2000, dn=2000)
    up = res.panels[PANEL_UPSTREAM].blocks[0].genomic_interval
    dn = res.panels[PANEL_DOWNSTREAM].blocks[0].genomic_interval
    assert (res.region_start, res.region_end) == (min(up[0], dn[0], res.gene_start),
                                                  max(up[1], dn[1], res.gene_end))
    block = res.panels[PANEL_GENOMIC].blocks[0]
    assert len(block.seq) == res.region_end - res.region_start + 1
    assert len(block.seq) == (len(res.upstream_seq) + (res.gene_end - res.gene_start + 1)
                              + len(res.downstream_seq))
    assert block.meta["coord_start"] == res.region_start


# --------------------------------------------------------------------------- #
# 负链
# --------------------------------------------------------------------------- #
def test_minus_strand_sequences(provider, ann_gff3, demo_dir):
    res = extract_one(provider, ann_gff3, "g2", "g2.1")
    exp = demo_dir["expected"]["g2"]["g2.1"]
    assert res.strand == "-"
    assert res.mrna == exp["mrna"]
    assert res.cds == exp["cds"]
    assert res.mrna.startswith("ATG"), "负链转录本应以反向互补后的 ATG 起始"
    assert res.cds.endswith("TAA")


def test_minus_strand_exon_order_descending(provider, ann_gff3, demo_dir):
    res = extract_one(provider, ann_gff3, "g2", "g2.1")
    intervals = [b.genomic_interval for b in res.panels[PANEL_EXON].blocks]
    starts = [iv[0] for iv in intervals]
    assert starts == sorted(starts, reverse=True), starts


def test_minus_strand_introns(provider, ann_gff3, demo_dir):
    res = extract_one(provider, ann_gff3, "g2", "g2.1")
    introns = res.panels[PANEL_INTRON].blocks
    assert len(introns) == 2
    # 转录顺序：先遇到靠右的内含子（4601-4800），再遇到靠左的（4201-4400）
    assert introns[0].genomic_interval == (4601, 4800)
    assert introns[1].genomic_interval == (4201, 4400)
    genome = demo_dir["genome_sequences"]["chr1"]
    assert introns[0].seq == revcomp(genome[4600:4800])
    assert all("反向互补" in b.note for b in introns)


def test_minus_strand_upstream_is_right_side(provider, ann_gff3, demo_dir):
    res = extract_one(provider, ann_gff3, "g2", "g2.1", up=1000, dn=1000)
    genome = demo_dir["genome_sequences"]["chr1"]
    # 负链基因的上游在基因组坐标更大的位置
    assert res.panels[PANEL_UPSTREAM].blocks[0].genomic_interval == (5001, 6000)
    assert res.upstream_seq == revcomp(genome[5000:6000])
    assert res.panels[PANEL_DOWNSTREAM].blocks[0].genomic_interval == (3001, 4000)
    assert res.downstream_seq == revcomp(genome[3000:4000])


def test_minus_strand_genomic_window_covers_flanks(provider, ann_gff3):
    """负链基因的窗口同样要横跨上下游。

    这里曾有一个真实 bug：``region_start`` 只取 ``min(tx_start, up_s)``、
    ``region_end`` 只取 ``max(tx_end, dn_e)``。负链的上游在坐标更大一侧、
    下游在更小一侧，于是上游被丢弃、下游偏移变成负数，
    基因组窗口里只剩下基因本身（prlrb 就是这样"看不到上下游"的）。
    """
    res = extract_one(provider, ann_gff3, "g2", "g2.1", up=1000, dn=1000)
    up = res.panels[PANEL_UPSTREAM].blocks[0].genomic_interval
    dn = res.panels[PANEL_DOWNSTREAM].blocks[0].genomic_interval
    assert up == (5001, 6000) and dn == (3001, 4000)
    # 窗口 = 下游(左) … 基因 … 上游(右)
    assert (res.region_start, res.region_end) == (3001, 6000)
    block = res.panels[PANEL_GENOMIC].blocks[0]
    assert len(block.seq) == 3000
    assert len(block.seq) == (len(res.downstream_seq) + res.gene_end - res.gene_start + 1
                              + len(res.upstream_seq))
    assert block.meta["coord_start"] == 3001

    # 窗口内的着色区间必须落在正坐标上，且能反推回各自的基因组区间
    spans = {r.kind: (r.genomic_start, r.genomic_end) for r in block.regions
             if r.kind in (REGION_UPSTREAM, REGION_DOWNSTREAM)}
    assert spans[REGION_UPSTREAM] == (5001, 6000)
    assert spans[REGION_DOWNSTREAM] == (3001, 4000)
    assert all(r.start0 >= 0 and r.end0 <= len(block.seq) for r in block.regions)



# --------------------------------------------------------------------------- #
# 单外显子 / 边界 / 特殊注释
# --------------------------------------------------------------------------- #
def test_single_exon_no_intron(provider, ann_gff3):
    res = extract_one(provider, ann_gff3, "g3", "g3.1")
    assert res.n_exons == 1
    assert res.n_introns == 0
    assert res.panels[PANEL_INTRON].is_empty()
    assert any(c.title == "单外显子" for c in res.checks)
    assert len(res.mrna) == 400
    assert res.mrna == res.cds + res.utr3


def test_boundary_clipping(provider, ann_gff3, demo_dir):
    """g4 是负链基因，紧贴 chrScaf 起点：
    它的**下游**在左边界外（只剩 50 bp），**上游**朝右反而取满 2000 bp。
    这一条同时验证了"上下游按转录方向定义"而非按坐标方向定义。
    """
    res = extract_one(provider, ann_gff3, "g4", "g4.1", up=2000, dn=2000)
    genome = demo_dir["genome_sequences"]["chrScaf"]
    assert res.upstream_requested == 2000
    assert res.upstream_actual == 2000
    assert res.upstream_clipped is False
    assert res.panels[PANEL_UPSTREAM].blocks[0].genomic_interval == (401, 2400)
    assert res.downstream_requested == 2000
    assert res.downstream_actual == 50
    assert res.downstream_clipped is True
    assert res.downstream_seq == revcomp(genome[0:50])
    assert any("下游请求" in w for w in res.warnings)


def test_cds_only_annotation_extraction(provider, ann_gff3):
    res = extract_one(provider, ann_gff3, "g5", "g5.1")
    assert res.n_exons == 2
    assert res.n_introns == 1
    assert len(res.mrna) == 220
    assert res.mrna == res.cds
    assert any("推导" in w for w in res.warnings)


def test_no_cds_transcript(provider, ann_gff3):
    res = extract_one(provider, ann_gff3, "g6", "g6.1")
    assert res.cds == ""
    assert res.panels[PANEL_CDS].is_empty()
    assert res.panels[PANEL_UTR].is_empty()
    assert any(c.title == "无 CDS" for c in res.checks)


def test_nonstandard_seqid_raises(provider, ann_gff3):
    gene = ann_gff3.get_gene("g8")
    tx = gene.canonical_transcript()
    with pytest.raises(SeqIDMismatchError) as exc:
        extract_gene(provider, gene, tx)
    msg = str(exc.value) + (exc.value.hint or "")
    assert "chrUn_1" in msg
    assert "chr1" in msg or "chrScaf" in msg, "应给出相近序列名建议"


def test_overlapping_exon_no_negative_intron(provider, ann_gff3):
    res = extract_one(provider, ann_gff3, "g9", "g9.1")
    assert res.n_introns == 0
    assert all(b.length > 0 for b in res.panels[PANEL_EXON].blocks)


# --------------------------------------------------------------------------- #
# 面板与坐标
# --------------------------------------------------------------------------- #
def test_genomic_panel_regions_contiguous(provider, ann_gff3):
    res = extract_one(provider, ann_gff3, "g1", "g1.1", up=2000, dn=2000)
    block = res.panels[PANEL_GENOMIC].blocks[0]
    assert sum(r.length for r in block.regions) == block.length
    pos = 0
    for reg in block.regions:
        assert reg.start0 == pos
        pos = reg.end0
    kinds = {r.kind for r in block.regions}
    assert REGION_UPSTREAM in kinds and REGION_DOWNSTREAM in kinds
    # 外显子要么以 exon 出现，要么被子区域（CDS / UTR）完全覆盖
    assert REGION_CDS in kinds
    assert REGION_UTR5 in kinds and REGION_UTR3 in kinds


def test_region_priority_cds_over_exon(provider, ann_gff3):
    res = extract_one(provider, ann_gff3, "g1", "g1.1", up=0, dn=0)
    block = res.panels[PANEL_GENOMIC].blocks[0]
    # 在 CDS 覆盖的基因组位置上，最终类别必须是 CDS 而不是外显子
    cds_reg = [r for r in block.regions if r.kind == REGION_CDS]
    assert cds_reg, "基因组面板应标出 CDS"
    for reg in cds_reg:
        for other in block.regions:
            if other is reg:
                continue
            overlap = min(reg.end0, other.end0) - max(reg.start0, other.start0)
            assert overlap <= 0, "CDS 区间不应与外显子区间在同一位置并列"


def test_block_genomic_coords_consistent(provider, ann_gff3):
    res = extract_one(provider, ann_gff3, "g2", "g2.1")
    for block in res.panels[PANEL_EXON].blocks:
        gs, ge = block.genomic_interval
        assert ge - gs + 1 == block.length
        assert block.coord_reverse is True
    for block in res.panels[PANEL_MRNA].blocks:
        assert block.coord_mode == "transcript"
        assert block.coord_reverse is False


def test_mrna_blocks_concat_equals_mrna(provider, ann_gff3):
    res = extract_one(provider, ann_gff3, "g2", "g2.1")
    assert res.panels[PANEL_MRNA].sequence == res.mrna
    assert res.panels[PANEL_EXON].sequence == res.mrna
    assert res.panels[PANEL_CDS].sequence == res.cds


def test_panel_lengths_reported(provider, ann_gff3):
    res = extract_one(provider, ann_gff3, "g1", "g1.1", up=500, dn=500)
    assert res.panels[PANEL_UPSTREAM].length == 500
    assert res.panels[PANEL_DOWNSTREAM].length == 500
    assert res.panels[PANEL_MRNA].length == 450
    assert res.panels[PANEL_CDS].length == 360
    assert res.panels[PANEL_UTR].length == 90


# --------------------------------------------------------------------------- #
# CDS phase
# --------------------------------------------------------------------------- #
def test_cds_phase_not_used_to_trim_bases(write_fasta, write_annotation):
    """phase 只用于报告，绝不能因此裁掉碱基。"""
    from gge.annotation import parse_annotation
    from gge.fasta import open_fasta

    seq = "ACGT" * 75  # 300 bp
    fa = write_fasta([("chrT", seq)], "t.fa")
    gff = (
        "##gff-version 3\n"
        "chrT\td\tgene\t21\t200\t.\t+\t.\tID=tg\n"
        "chrT\td\tmRNA\t21\t200\t.\t+\t.\tID=tg.1;Parent=tg\n"
        "chrT\td\texon\t21\t80\t.\t+\t.\tID=e1;Parent=tg.1\n"
        "chrT\td\texon\t121\t200\t.\t+\t.\tID=e2;Parent=tg.1\n"
        "chrT\td\tCDS\t21\t80\t.\t+\t0\tID=c1;Parent=tg.1\n"
        "chrT\td\tCDS\t121\t160\t.\t+\t2\tID=c2;Parent=tg.1\n"
    )
    ann = parse_annotation(write_annotation(gff, "t.gff3"))
    p = open_fasta(fa, prefer_memory=True)
    res = extract_gene(p, ann.get_gene("tg"), ann.get_gene("tg").transcripts["tg.1"])
    # 60 + 40，若按 phase 裁剪会得到 58
    assert len(res.cds) == 100
    assert res.cds == seq[20:80] + seq[120:160]
    assert any("phase" in w for w in res.warnings), "phase 不一致应提示但不裁剪"


def test_phase_consistent_no_warning(write_fasta, write_annotation):
    from gge.annotation import parse_annotation
    from gge.fasta import open_fasta

    seq = "ACGT" * 75
    fa = write_fasta([("chrT", seq)], "t2.fa")
    gff = (
        "##gff-version 3\n"
        "chrT\td\tgene\t21\t200\t.\t+\t.\tID=tg2\n"
        "chrT\td\tmRNA\t21\t200\t.\t+\t.\tID=tg2.1;Parent=tg2\n"
        "chrT\td\texon\t21\t80\t.\t+\t.\tID=e1;Parent=tg2.1\n"
        "chrT\td\texon\t121\t200\t.\t+\t.\tID=e2;Parent=tg2.1\n"
        "chrT\td\tCDS\t21\t80\t.\t+\t0\tID=c1;Parent=tg2.1\n"
        "chrT\td\tCDS\t121\t160\t.\t+\t0\tID=c2;Parent=tg2.1\n"
    )
    ann = parse_annotation(write_annotation(gff, "t2.gff3"))
    p = open_fasta(fa, prefer_memory=True)
    res = extract_gene(p, ann.get_gene("tg2"), ann.get_gene("tg2").transcripts["tg2.1"])
    assert not any("phase" in w for w in res.warnings)


def test_unusual_bases_kept_and_reported(write_fasta, write_annotation):
    from gge.annotation import parse_annotation
    from gge.fasta import open_fasta

    seq = "N" * 50 + "ACGT" * 50
    fa = write_fasta([("chrN", seq)], "n.fa")
    gff = (
        "##gff-version 3\n"
        "chrN\td\tgene\t1\t50\t.\t+\t.\tID=ng\n"
        "chrN\td\tmRNA\t1\t50\t.\t+\t.\tID=ng.1;Parent=ng\n"
        "chrN\td\texon\t1\t50\t.\t+\t.\tID=e1;Parent=ng.1\n"
    )
    ann = parse_annotation(write_annotation(gff, "n.gff3"))
    p = open_fasta(fa, prefer_memory=True)
    res = extract_gene(p, ann.get_gene("ng"), ann.get_gene("ng").transcripts["ng.1"])
    assert res.mrna == "N" * 50, "非标准碱基必须原样保留"
    assert any(c.title == "存在非 ACGT 字符" and c.level == "warn" for c in res.checks)
