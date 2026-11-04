"""基因 / 转录本检索。

排序策略（分数越高越靠前）
--------------------------
* 100  Gene ID 完全一致（忽略大小写）
*  98  Gene Symbol / Name 完全一致
*  96  Transcript ID 完全一致
*  92  别名（locus_tag、gene、Name、product）完全一致
*  80  前缀匹配
*  70  子串匹配
*  60  模糊匹配（difflib 相似度 ≥ 0.72，用于 slc6a6 vs slc6a6a 这类拼写偏差）

同一基因命中多个字段时取最高分，并记录命中的字段名，界面据此说明"为什么这条结果被列出"。
"""

from __future__ import annotations

import difflib
from dataclasses import dataclass
from typing import Dict, List, Optional

from .annotation import AnnotationModel
from .models import Gene


@dataclass
class SearchHit:
    gene_id: str
    gene_name: str
    seqid: str
    start: int
    end: int
    strand: str
    n_transcripts: int
    matched_by: str = ""
    matched_value: str = ""
    score: float = 0.0
    tx_id: str = ""          # 以 Transcript ID 命中时记录，便于直接定位

    def coords_text(self) -> str:
        return f"{self.start:,}-{self.end:,}"

    def label(self) -> str:
        if self.gene_name and self.gene_name != self.gene_id:
            return f"{self.gene_name} ({self.gene_id})"
        return self.gene_id


def _gene_hit(gene: Gene, matched_by: str, matched_value: str, score: float,
              tx_id: str = "") -> SearchHit:
    return SearchHit(
        gene_id=gene.gene_id,
        gene_name=gene.gene_name or gene.gene_id,
        seqid=gene.seqid,
        start=gene.start,
        end=gene.end,
        strand=gene.strand or ".",
        n_transcripts=gene.n_transcripts,
        matched_by=matched_by,
        matched_value=matched_value,
        score=score,
        tx_id=tx_id,
    )


def search_genes(model: AnnotationModel, query: str, limit: int = 300) -> List[SearchHit]:
    """模糊搜索基因，返回按相关度排序的结果。"""
    q = (query or "").strip()
    if not q:
        return []
    ql = q.lower()
    hits: Dict[str, SearchHit] = {}

    def push(gene: Gene, by: str, value: str, score: float, tx_id: str = "") -> None:
        old = hits.get(gene.gene_id)
        if old is None or score > old.score:
            hits[gene.gene_id] = _gene_hit(gene, by, value, score, tx_id)

    for gene in model.genes.values():
        gid = gene.gene_id
        gidl = gid.lower()
        names = {gene.gene_name.lower()} if gene.gene_name else set()
        for key in ("locus_tag", "name", "gene", "product"):
            v = gene.attrs.get(key)
            if v:
                names.add(v.lower())

        if gidl == ql:
            push(gene, "Gene ID", gid, 100)
            continue
        if ql in names:
            push(gene, "Gene Symbol", gene.gene_name or q, 98)
            continue

        # 转录本 ID
        matched_tx = ""
        for tid in gene.transcripts:
            if tid.lower() == ql:
                matched_tx = tid
                break
        if matched_tx:
            push(gene, "Transcript ID", matched_tx, 96, tx_id=matched_tx)
            continue

        if gidl.startswith(ql):
            push(gene, "Gene ID 前缀", gid, 80)
            continue
        if any(n.startswith(ql) for n in names):
            push(gene, "Gene Symbol 前缀", gene.gene_name or q, 80)
            continue
        if ql in gidl:
            push(gene, "Gene ID 子串", gid, 70)
            continue
        if any(ql in n for n in names if n):
            push(gene, "Gene Symbol 子串", gene.gene_name or q, 70)
            continue
        tx_sub = next((t for t in gene.transcripts if ql in t.lower()), "")
        if tx_sub:
            push(gene, "Transcript ID 子串", tx_sub, 68, tx_id=tx_sub)
            continue

        # 模糊匹配（限制长度差，避免噪声）
        best = 0.0
        best_val = ""
        for cand in [gidl] + list(names):
            if not cand or abs(len(cand) - len(ql)) > max(3, len(ql) // 2):
                continue
            ratio = difflib.SequenceMatcher(None, ql, cand).ratio()
            if ratio > best:
                best, best_val = ratio, cand
        if best >= 0.72:
            push(gene, "模糊匹配", best_val, 60 + best * 5)

    out = sorted(hits.values(), key=lambda h: (-h.score, h.gene_name.lower(), h.gene_id))
    return out[:limit]


def resolve_query(model: AnnotationModel, query: str) -> Optional[SearchHit]:
    """把用户输入解析为唯一基因；无法唯一确定时返回 None。"""
    hits = search_genes(model, query, limit=50)
    if not hits:
        return None
    if len(hits) == 1:
        return hits[0]
    top = hits[0]
    if top.score >= 100:
        return top
    if top.score - hits[1].score >= 15:
        return top
    return None
