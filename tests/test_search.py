"""基因检索测试。"""

from __future__ import annotations

from gge.search import resolve_query, search_genes


def test_exact_gene_id(ann_gff3):
    hits = search_genes(ann_gff3, "g1")
    assert hits[0].gene_id == "g1"
    assert hits[0].score == 100
    assert hits[0].matched_by == "Gene ID"


def test_exact_gene_symbol_case_insensitive(ann_gff3):
    hits = search_genes(ann_gff3, "geneone")
    assert hits[0].gene_id == "g1"
    assert hits[0].score == 98
    assert hits[0].n_transcripts == 2


def test_transcript_id_search(ann_gff3):
    hits = search_genes(ann_gff3, "g1.2")
    assert hits[0].gene_id == "g1"
    assert hits[0].tx_id == "g1.2"
    assert hits[0].matched_by == "Transcript ID"


def test_prefix_and_substring(ann_gff3):
    hits = search_genes(ann_gff3, "GeneT")
    assert hits and hits[0].gene_id == "g3"
    assert hits[0].score >= 70

    sub = search_genes(ann_gff3, "eneFiv")
    assert any(h.gene_id == "g5" for h in sub)


def test_fuzzy_search(ann_gff3):
    hits = search_genes(ann_gff3, "GeneOnes")
    assert hits and hits[0].gene_id == "g1"
    assert hits[0].matched_by in ("模糊匹配", "Gene Symbol")


def test_multi_results_returned(ann_gff3):
    hits = search_genes(ann_gff3, "Gene")
    assert len(hits) >= 8
    # 结果按分数降序
    scores = [h.score for h in hits]
    assert scores == sorted(scores, reverse=True)


def test_empty_query(ann_gff3):
    assert search_genes(ann_gff3, "   ") == []


def test_limit(ann_gff3):
    hits = search_genes(ann_gff3, "Gene", limit=3)
    assert len(hits) == 3


def test_resolve_unique(ann_gff3):
    assert resolve_query(ann_gff3, "GeneFour").gene_id == "g4"
    assert resolve_query(ann_gff3, "g9").gene_id == "g9"


def test_resolve_ambiguous_returns_none(ann_gff3):
    # "Gene" 命中多条且分数接近，不能替用户随便挑一个
    assert resolve_query(ann_gff3, "Gene") is None
