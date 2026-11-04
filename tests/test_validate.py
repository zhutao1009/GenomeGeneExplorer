"""参考 mRNA / CDS 校验测试。"""

from __future__ import annotations

import pytest

pytest.importorskip("Bio", reason="Biopython 未安装，跳过比对相关测试")

from gge.extract import extract_gene  # noqa: E402
from gge.validate import (STATUS_DIFFERENT, STATUS_IDENTICAL, STATUS_MANUAL,  # noqa: E402
                          STATUS_UNMATCHED, compare_sequences, load_reference,
                          match_reference, validate_extraction)


def _run(provider, ann, gid, refs_m, refs_c, tx_id="", **kw):
    gene = ann.get_gene(gid)
    tx = gene.transcripts[tx_id] if tx_id else gene.canonical_transcript()
    res = extract_gene(provider, gene, tx, kw.pop("up", 2000), kw.pop("dn", 2000))
    return res, validate_extraction(res, gene, tx, refs_m, refs_c, **kw)


# --------------------------------------------------------------------------- #
# 比对算法本身
# --------------------------------------------------------------------------- #
def test_compare_identical():
    c = compare_sequences("ACGTACGT", "ACGTACGT")
    assert c.identical and c.mismatches == 0 and c.identity == 100.0


def test_compare_substitution():
    c = compare_sequences("ACGTACGT", "ACGTTCGT")
    assert not c.identical
    assert c.mismatches == 1 and c.insertions == 0 and c.deletions == 0
    assert c.events[0].kind == "mismatch"
    assert c.events[0].ref_pos == 5


def test_compare_insertion_detected_by_alignment():
    c = compare_sequences("ACGTACGT", "ACGTTACGT")
    assert c.insertions == 1 and c.deletions == 0
    assert c.events[0].kind == "insertion"
    # 剩余 8 个碱基仍然对齐，不能因为一位插入就报 5 处不一致
    assert c.mismatches == 0


def test_compare_deletion():
    c = compare_sequences("ACGTTACGT", "ACGTACGT")
    assert c.deletions == 1 and c.insertions == 0 and c.mismatches == 0


def test_compare_both_indels():
    c = compare_sequences("AAAACCCCGGGG", "AAAACCCGGGGTT")
    assert c.indel_bases == 3  # 1 个缺失 + 2 个插入
    assert c.identity < 100.0


def test_compare_empty():
    c = compare_sequences("", "")
    assert c.identical
    c2 = compare_sequences("ACGT", "")
    assert not c2.identical


# --------------------------------------------------------------------------- #
# mRNA
# --------------------------------------------------------------------------- #
def test_mrna_exact_match(provider, ann_gff3, mrna_refs, cds_refs):
    res, v = _run(provider, ann_gff3, "g1", mrna_refs, cds_refs, "g1.1")
    m = v["mRNA"]
    assert m.status == STATUS_IDENTICAL, m.notes
    assert m.ref_len == 450 and m.obs_len == 450
    assert m.comparison.mismatches == 0


def test_alt_transcript_matched_by_its_own_id(provider, ann_gff3, mrna_refs, cds_refs):
    res, v = _run(provider, ann_gff3, "g1", mrna_refs, cds_refs, "g1.2")
    m = v["mRNA"]
    assert m.status == STATUS_IDENTICAL, (m.ref_id, m.notes)
    assert "g1.2" in m.ref_id


def test_mrna_single_base_deletion(provider, ann_gff3, mrna_refs, cds_refs):
    res, v = _run(provider, ann_gff3, "g2", mrna_refs, cds_refs, "g2.1")
    m = v["mRNA"]
    assert m.status == STATUS_DIFFERENT
    c = m.comparison
    # 参考序列比提取结果少 1 个碱基 → 观测序列相对参考"多出"1 个碱基
    assert c.insertions == 1
    assert c.deletions == 0
    assert c.mismatches == 0, "单碱基缺失不应被报成一串错配"
    assert 99.0 < c.identity < 100.0
    assert any("排查方向" in n for n in m.notes)


def test_mrna_single_substitution(provider, ann_gff3, mrna_refs, cds_refs):
    res, v = _run(provider, ann_gff3, "g3", mrna_refs, cds_refs, "g3.1")
    m = v["mRNA"]
    assert m.status == STATUS_DIFFERENT
    assert m.comparison.mismatches == 1
    assert m.comparison.insertions == 0 and m.comparison.deletions == 0
    assert m.comparison.events[0].ref_pos == 101


def test_mrna_ambiguous_candidates_need_manual(provider, ann_gff3, mrna_refs, cds_refs):
    res, v = _run(provider, ann_gff3, "g6", mrna_refs, cds_refs, "g6.1")
    m = v["mRNA"]
    assert m.status == STATUS_MANUAL
    assert len(m.candidates) >= 2
    assert "得分相同" in " ".join(m.notes) or "无法自动判定" in " ".join(m.notes)


def test_mrna_manual_override(provider, ann_gff3, mrna_refs, cds_refs):
    res, v = _run(provider, ann_gff3, "g6", mrna_refs, cds_refs, "g6.1",
                  manual_mrna_id="g6 copyB")
    m = v["mRNA"]
    assert m.status == STATUS_DIFFERENT
    assert m.ref_token == "g6"
    assert m.comparison.identity < 50.0


def test_mrna_unmatched(provider, ann_gff3, mrna_refs, cds_refs):
    res, v = _run(provider, ann_gff3, "g4", mrna_refs, cds_refs, "g4.1")
    m = v["mRNA"]
    assert m.status == STATUS_UNMATCHED
    assert "没有任何记录能匹配" in " ".join(m.notes)


# --------------------------------------------------------------------------- #
# CDS
# --------------------------------------------------------------------------- #
def test_cds_exact(provider, ann_gff3, mrna_refs, cds_refs):
    res, v = _run(provider, ann_gff3, "g1", mrna_refs, cds_refs, "g1.1")
    c = v["CDS"]
    assert c.status == STATUS_IDENTICAL, c.notes
    assert c.ref_len == 360


def test_cds_with_insertion(provider, ann_gff3, mrna_refs, cds_refs):
    res, v = _run(provider, ann_gff3, "g3", mrna_refs, cds_refs, "g3.1")
    c = v["CDS"]
    assert c.status == STATUS_DIFFERENT
    # 参考 CDS 多出 3 bp → 相对参考，观测序列少了 3 个碱基
    assert c.comparison.deletions == 3
    assert c.comparison.insertions == 0
    assert c.comparison.mismatches == 0


def test_cds_coding_checks_present(provider, ann_gff3, mrna_refs, cds_refs):
    res, v = _run(provider, ann_gff3, "g1", mrna_refs, cds_refs, "g1.1")
    titles = [t for _lvl, t, _d in v["CDS"].coding_checks]
    assert "阅读框" in titles and "起始密码子" in titles and "终止密码子" in titles


def test_no_cds_gives_unmatched(provider, ann_gff3, mrna_refs, cds_refs):
    res, v = _run(provider, ann_gff3, "g6", mrna_refs, cds_refs, "g6.1")
    assert v["CDS"].status == STATUS_UNMATCHED
    assert "没有 CDS" in " ".join(v["CDS"].notes)


# --------------------------------------------------------------------------- #
# 蛋白（扩展校验）
# --------------------------------------------------------------------------- #
def test_protein_validation(provider, ann_gff3, mrna_refs, cds_refs, protein_refs):
    res, v = _run(provider, ann_gff3, "g1", mrna_refs, cds_refs, "g1.1",
                  protein_refs=protein_refs)
    p = v["protein"]
    assert p.status == STATUS_IDENTICAL, p.notes


# --------------------------------------------------------------------------- #
# 匹配逻辑
# --------------------------------------------------------------------------- #
def test_no_order_based_matching(tmp_path, provider, ann_gff3):
    from gge.validate import RefSet, RefRecord

    gene = ann_gff3.get_gene("g1")
    tx = gene.transcripts["g1.1"]
    rs = RefSet(kind="mRNA", records=[
        RefRecord(rid="zzz unrelated", token="zzz", seq="A" * 450, index=0),
        RefRecord(rid="g1.1 correct", token="g1.1", seq="A" * 450, index=1),
    ])
    rs.id_index = {r.rid: r for r in rs.records} | {r.token: r for r in rs.records}
    decision = match_reference(rs, gene, tx, 450)
    assert decision.ref_id == "g1.1 correct"


def test_unmatched_when_ids_do_not_match(provider, ann_gff3):
    from gge.validate import RefSet, RefRecord

    gene = ann_gff3.get_gene("g1")
    tx = gene.transcripts["g1.1"]
    rs = RefSet(kind="mRNA", records=[
        RefRecord(rid="something_else", token="something_else", seq="A" * 100, index=0),
    ])
    rs.id_index = {r.token: r for r in rs.records}
    decision = match_reference(rs, gene, tx, 450)
    assert decision.status == STATUS_UNMATCHED


# --------------------------------------------------------------------------- #
# RefSeq 风格的蛋白匹配（protein.faa 标题只有 XP_ 蛋白 ID）
# --------------------------------------------------------------------------- #
def _refset(kind, pairs):
    """用 (标题, 序列) 快速造一个 RefSet（带 id_index，行为与 load_reference 一致）。"""
    from gge.validate import RefRecord, RefSet

    rs = RefSet(kind=kind)
    for i, (header, seq) in enumerate(pairs):
        token = header.split()[0]
        rs.records.append(RefRecord(rid=header, token=token, seq=seq, index=i))
        rs.id_index.setdefault(header, rs.records[-1])
        rs.id_index.setdefault(token, rs.records[-1])
    return rs


def test_protein_id_taken_from_matched_cds_title():
    """RefSeq 真实标题：CDS 里写着 protein_id=XP_…，蛋白文件里就是 XP_…。"""
    from gge.validate import protein_id_from_cds

    cds = _refset("CDS", [(
        "lcl|NC_067220.1_cds_XP_051755185.1_1 [gene=LOC127515509] "
        "[protein=zinc finger protein 750-like] [protein_id=XP_051755185.1] "
        "[gbkey=CDS]", "ATGAAA")]).records[0]
    prot = _refset("protein", [("XP_051755185.1 zinc finger protein 750-like", "MK"),
                              ("XP_051715495.1 something else", "MKV")])

    assert protein_id_from_cds(cds, prot) == "XP_051755185.1"


def test_protein_id_absent_from_protein_file_returns_empty():
    from gge.validate import protein_id_from_cds

    cds = _refset("CDS", [("lcl|X_cds_YP_9_1 [protein_id=YP_9]", "ATG")]).records[0]
    prot = _refset("protein", [("XP_1 unrelated", "MK")])
    assert protein_id_from_cds(cds, prot) == ""


def test_protein_id_ambiguous_returns_empty():
    """标题里出现两个都存在的蛋白 ID 时不能二选一——交回打分/人工确认。"""
    from gge.validate import protein_id_from_cds

    cds = _refset("CDS", [("[protein_id=XP_1] XP_2", "ATG")]).records[0]
    prot = _refset("protein", [("XP_1 a", "MK"), ("XP_2 b", "MK")])
    assert protein_id_from_cds(cds, prot) == ""


def test_protein_id_without_cds_record_returns_empty():
    from gge.validate import protein_id_from_cds

    prot = _refset("protein", [("XP_1 a", "MK")])
    assert protein_id_from_cds(None, prot) == ""


def test_protein_matches_by_id_from_cds_when_titles_lack_transcript_id(
        provider, ann_gff3, mrna_refs, cds_refs):
    """端到端：蛋白文件只有 XP_ ID 时，仍应自动匹配并判定一致。"""
    from gge.validate import STATUS_IDENTICAL, translate

    gene = ann_gff3.get_gene("g1")
    tx = gene.transcripts["g1.1"]
    res = extract_gene(provider, gene, tx, 2000, 2000)

    # 把 CDS 参考改写成 RefSeq 风格（标题带 protein_id=XP_TEST.1）
    cds2 = _refset("CDS", [
        ("lcl|chr1_cds_XP_TEST.1_1 [gene=g1] [protein_id=XP_TEST.1] [gbkey=CDS]", res.cds),
    ])
    # 蛋白文件标题只有 XP_ ID，完全没有 g1.1 / g1 之类的字样
    prot = _refset("protein", [("XP_TEST.1 PREDICTED: something", translate(res.cds))])

    v = validate_extraction(res, gene, tx, mrna_refs, cds2, protein_refs=prot)
    assert v["CDS"].status == STATUS_IDENTICAL
    assert v["protein"].status == STATUS_IDENTICAL, v["protein"].notes
    assert "protein_id=XP_TEST.1" in " ".join(v["protein"].notes)


def test_protein_length_participates_in_scoring(provider, ann_gff3):
    """没有 protein_id 可依据时，长度至少要让候选浮出来，而不是一条都没有。"""
    from gge.validate import STATUS_MANUAL, translate

    gene = ann_gff3.get_gene("g1")
    tx = gene.transcripts["g1.1"]
    res = extract_gene(provider, gene, tx, 2000, 2000)
    prot_clean = translate(res.cds).replace("*", "")

    prot = _refset("protein", [
        ("XP_1 xxx", prot_clean),           # 长度一致但 ID 无关
        ("XP_2 yyy", "M" * 13),             # 长度差很远
    ])
    # 期望长度用去终止符后的翻译长度（与 validate_extraction 内部一致）
    decision = match_reference(prot, gene, tx, len(prot_clean))
    assert decision.ref_id, decision.notes
    assert decision.ref_id.startswith("XP_1")
    # 唯一命中就自动用；多条同长度才会退回人工确认
    assert decision.status in ("matched", STATUS_MANUAL)
