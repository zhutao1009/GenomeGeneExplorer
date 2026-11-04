"""参考 mRNA / CDS 的匹配与序列校验。

匹配原则
--------
绝不按文件顺序匹配，也绝不在多个候选中默默取第一条。流程是：

1. 用 **Transcript ID 完全一致** → 最高分；
2. **Gene ID 完全一致**；
3. ID 互为子串（``tx1`` 出现在 ``>gene1|tx1|...`` 里）；
4. 标题中出现 Gene Symbol；
5. **长度一致**作为独立加分项（只在前面几条都没命中时起作用）。

分数第一且唯一 → 自动匹配；并列 → 状态 ``manual``，界面必须让用户手动指定。

比对
----
用 Biopython 的 ``PairwiseAligner`` 做全局比对（Needleman–Wunsch），因此
插入/缺失能被正确识别并给出比对坐标，而不是"逐位错位一片"。
超长序列（> 60 kb）会退化为"前后缀锚定"的粗略定位，并在结果中显式标注
``coarse=True``，不会假装做过全序列比对。
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from .coords import revcomp, unusual_bases
from .errors import FileFormatError
from .fasta import read_fasta_records
from .models import Gene, GeneExtraction, Transcript
from .sequtil import (first_codon, internal_stops, is_start_codon, is_stop_codon,
                      last_codon, translate)

STATUS_IDENTICAL = "identical"
STATUS_DIFFERENT = "different"
STATUS_UNMATCHED = "unmatched"
STATUS_MANUAL = "manual"
STATUS_NOT_RUN = "not_run"

STATUS_CN = {
    STATUS_IDENTICAL: "完全一致",
    STATUS_DIFFERENT: "存在差异",
    STATUS_UNMATCHED: "无法匹配",
    STATUS_MANUAL: "需要人工检查",
    STATUS_NOT_RUN: "未校验",
}

STATUS_COLORS = {
    STATUS_IDENTICAL: "#1b7f3b",
    STATUS_DIFFERENT: "#b45309",
    STATUS_UNMATCHED: "#b91c1c",
    STATUS_MANUAL: "#7c3aed",
    STATUS_NOT_RUN: "#6b7280",
}


# --------------------------------------------------------------------------- #
# 参考序列容器
# --------------------------------------------------------------------------- #
@dataclass
class RefRecord:
    rid: str          # 完整标题
    token: str        # 第一个空白前的 token，常用于 ID
    seq: str
    description: str = ""
    index: int = 0


@dataclass
class RefSet:
    path: str = ""
    kind: str = "mRNA"           # 'mRNA' / 'CDS' / 'protein'
    records: List[RefRecord] = field(default_factory=list)
    id_index: Dict[str, RefRecord] = field(default_factory=dict)

    def get(self, rid: str) -> Optional[RefRecord]:
        return self.id_index.get(rid)

    def __len__(self) -> int:
        return len(self.records)

    def ids(self) -> List[str]:
        return [r.rid for r in self.records]

    def summary(self) -> str:
        total = sum(len(r.seq) for r in self.records)
        return f"{len(self.records):,} 条记录，共 {total:,} bp"


def load_reference(path: str, kind: str, progress_cb=None,
                  cancel_check=None) -> RefSet:
    """读取参考 FASTA 并建立 ID 索引。"""
    records = read_fasta_records(path, progress_cb=progress_cb,
                                cancel_check=cancel_check)
    rs = RefSet(path=path, kind=kind)
    for i, (header, seq) in enumerate(records):
        token = header.split()[0] if header.split() else header
        desc = header[len(token):].strip()
        rec = RefRecord(rid=header, token=token, seq=seq, description=desc, index=i)
        rs.records.append(rec)
        rs.id_index.setdefault(header, rec)
        rs.id_index.setdefault(token, rec)
    return rs


# --------------------------------------------------------------------------- #
# ID 匹配
# --------------------------------------------------------------------------- #
@dataclass
class Candidate:
    rid: str
    token: str
    score: float
    reasons: List[str] = field(default_factory=list)
    length: int = 0

    def reason_text(self) -> str:
        return "；".join(self.reasons)


@dataclass
class MatchDecision:
    status: str
    ref_id: str = ""
    ref_token: str = ""
    candidates: List[Candidate] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)

    def candidate_text(self, limit: int = 8) -> str:
        if not self.candidates:
            return "（无候选）"
        rows = []
        for c in self.candidates[:limit]:
            rows.append(f"{c.token}  分数 {c.score:.0f}  {c.reason_text()}")
        if len(self.candidates) > limit:
            rows.append(f"… 另有 {len(self.candidates) - limit} 条候选")
        return "\n".join(rows)


def _split_title(text: str) -> List[str]:
    """把标题拆成可用于匹配的片段（支持 | 、空格、逗号分隔）。"""
    out = set()
    for sep in ("|", " ", ",", ";", "\t"):
        for part in text.replace(sep, "\x00").split("\x00"):
            part = part.strip()
            if part:
                out.add(part)
    return sorted(out)


_PROTEIN_ID_RE = re.compile(r"protein_id=([^\s;\],]+)")


def protein_id_from_cds(cds_rec: Optional[RefRecord],
                        protein_refs: RefSet) -> str:
    """从已匹配的参考 CDS 标题里取出蛋白 ID，并在参考蛋白文件里确认存在。

    为什么需要这一步
    ----------------
    RefSeq 的 ``protein.faa`` 标题只有蛋白 ID（``>XP_051755185.1 …``），**不含转录本 ID**，
    而原来的蛋白匹配只靠「转录本 ID / 基因 ID / 基因名」，且 ``expected_len`` 传 0
    （长度项不加分）——两条路都断，结果是蛋白参考文件永远"无法匹配"，
    用户只能从 5 万条候选里手动挑一条。

    但参考 CDS 的标题里写着答案：
    ``lcl|NC_067220.1_cds_XP_051755185.1_1 [gene=LOC…] [protein_id=XP_051755185.1] …``
    CDS 已经匹配成功（那是逐碱基验证过的），它指向的蛋白就是唯一正确的对照。

    只在**唯一命中**时才返回：多个候选说明这条 CDS 标题有歧义，
    宁可返回空字符串（退回打分 + 人工确认），也不猜。
    """
    if cds_rec is None or not protein_refs or not protein_refs.records:
        return ""
    found: List[str] = []
    m = _PROTEIN_ID_RE.search(cds_rec.rid)
    if m:
        found.append(m.group(1).strip("[]"))
    # 退一步：标题里出现的词如果正好是参考蛋白的 ID 也算（Ensembl 之类没有 protein_id= 的写法）
    known = {r.token for r in protein_refs.records}
    for tok in _split_title(cds_rec.rid):
        tok = tok.strip("[]")
        if tok in known and tok not in found:
            found.append(tok)
    hits = [c for c in found if protein_refs.get(c) is not None]
    return hits[0] if len(hits) == 1 else ""


def score_candidates(refset: RefSet, gene: Gene, tx: Transcript,
                     expected_len: int) -> List[Candidate]:
    """对所有参考记录打分，返回降序候选列表。"""
    tx_id = tx.tx_id
    gene_id = gene.gene_id
    gene_name = gene.gene_name or ""
    tx_l = tx_id.lower()
    gene_l = gene_id.lower()
    name_l = gene_name.lower()

    scored: List[Candidate] = []
    for rec in refset.records:
        score = 0.0
        reasons: List[str] = []
        tokens = [t.lower() for t in _split_title(rec.rid)]
        token_l = rec.token.lower()

        if token_l == tx_l or any(t == tx_l for t in tokens):
            score += 100
            reasons.append("Transcript ID 完全一致")
        elif token_l == gene_l or any(t == gene_l for t in tokens):
            score += 85
            reasons.append("Gene ID 完全一致")
        else:
            if tx_l and tx_l in rec.rid.lower():
                score += 60
                reasons.append("标题包含 Transcript ID")
            elif gene_l and gene_l in rec.rid.lower():
                score += 50
                reasons.append("标题包含 Gene ID")
            elif name_l and len(name_l) >= 3 and name_l in rec.rid.lower():
                score += 45
                reasons.append("标题包含 Gene Symbol")

        if expected_len and len(rec.seq) == expected_len:
            score += 40
            reasons.append(f"长度一致（{expected_len:,} bp）")
        elif expected_len and abs(len(rec.seq) - expected_len) <= max(3, expected_len // 100):
            score += 12
            reasons.append("长度接近（差异 ≤ 1%）")

        if score > 0:
            scored.append(Candidate(rid=rec.rid, token=rec.token, score=score,
                                    reasons=reasons, length=len(rec.seq)))
    scored.sort(key=lambda c: (-c.score, c.token))
    return scored


def match_reference(refset: RefSet, gene: Gene, tx: Transcript,
                    expected_len: int, manual_id: Optional[str] = None,
                    reason: str = "") -> MatchDecision:
    """决定使用哪条参考序列。"""
    if not refset or not refset.records:
        return MatchDecision(STATUS_NOT_RUN, notes=["未加载参考序列文件"])

    if manual_id:
        rec = refset.get(manual_id)
        if rec is None:
            for r in refset.records:
                if r.token == manual_id:
                    rec = r
                    break
        if rec is None:
            return MatchDecision(STATUS_UNMATCHED,
                                 notes=[f"手动指定的参考序列 {manual_id} 不在文件中"])
        return MatchDecision("matched", ref_id=rec.rid, ref_token=rec.token,
                             notes=[reason or "按用户手动指定的对应关系进行校验"])

    cands = score_candidates(refset, gene, tx, expected_len)
    if not cands:
        notes = [f"参考文件中没有任何记录能匹配 {tx.tx_id} / {gene.gene_id}"]
        if refset.records:
            notes.append("前几条参考标题示例：" +
                         "；".join(r.token for r in refset.records[:3]))
        return MatchDecision(STATUS_UNMATCHED, notes=notes)

    best = cands[0]
    top_ties = [c for c in cands if abs(c.score - best.score) < 1e-6]
    if len(top_ties) > 1:
        return MatchDecision(
            STATUS_MANUAL, candidates=cands,
            notes=[f"{len(top_ties)} 条参考序列得分相同（{best.score:.0f}），"
                   f"无法自动判定，请手动确认对应关系"],
        )
    if best.score < 40:
        return MatchDecision(
            STATUS_MANUAL, candidates=cands,
            notes=["候选得分过低，无法可靠自动匹配，请手动确认"],
        )
    return MatchDecision("matched", ref_id=best.rid, ref_token=best.token,
                         candidates=cands)


# --------------------------------------------------------------------------- #
# 序列比对
# --------------------------------------------------------------------------- #
@dataclass
class DiffEvent:
    kind: str            # mismatch / insertion / deletion
    ref_pos: int         # 1-based；插入事件为插入点前一个碱基的位置
    qry_pos: int
    ref_base: str
    qry_base: str

    def describe(self) -> str:
        if self.kind == "mismatch":
            return f"ref {self.ref_pos:,}: {self.ref_base} → obs {self.qry_base}"
        if self.kind == "insertion":
            return f"obs 在 {self.qry_pos:,} 后多出 {self.qry_base!r}（ref 无对应碱基）"
        return f"ref 在 {self.ref_pos:,} 处缺失 {self.ref_base!r}"


@dataclass
class Comparison:
    identical: bool = False
    ref_len: int = 0
    qry_len: int = 0
    aligned_columns: int = 0
    matches: int = 0
    mismatches: int = 0
    insertions: int = 0      # 观测序列相对参考多出的碱基数
    deletions: int = 0       # 参考相对观测多出的碱基数
    identity: float = 0.0
    events: List[DiffEvent] = field(default_factory=list)
    aligned_ref: str = ""
    aligned_qry: str = ""
    coarse: bool = False
    notes: List[str] = field(default_factory=list)
    first_diff_ref: int = -1
    first_diff_qry: int = -1

    @property
    def indel_bases(self) -> int:
        return self.insertions + self.deletions


def _make_aligner():
    from Bio.Align import PairwiseAligner

    al = PairwiseAligner()
    al.mode = "global"
    al.match_score = 2.0
    al.mismatch_score = -3.0
    al.open_gap_score = -6.0
    al.extend_gap_score = -1.0
    al.end_gap_score = -0.5
    return al


def _gapped_strings(alignment) -> Tuple[str, str]:
    """从 Alignment 对象取出带 gap 的两条字符串。"""
    a = alignment[0]
    b = alignment[1]
    if not isinstance(a, str):
        a = str(a)
    if not isinstance(b, str):
        b = str(b)
    return a, b


def compare_sequences(ref: str, qry: str, max_bases: int = 60_000) -> Comparison:
    """全局比对两条序列并统计差异。"""
    ref = (ref or "").strip().upper()
    qry = (qry or "").strip().upper()
    res = Comparison(ref_len=len(ref), qry_len=len(qry))

    if ref == qry:
        res.identical = True
        res.matches = len(ref)
        res.aligned_columns = len(ref)
        res.identity = 100.0
        res.aligned_ref = ref
        res.aligned_qry = qry
        return res

    if not ref or not qry:
        res.notes.append("其中一条序列为空，无法比对")
        res.identity = 0.0
        return res

    if max(len(ref), len(qry)) > max_bases:
        return _coarse_compare(ref, qry, res)

    aligner = _make_aligner()
    try:
        alignments = aligner.align(ref, qry)
        alignment = alignments[0]
    except Exception as exc:  # pragma: no cover - 极端情况
        res.notes.append(f"比对失败：{exc}")
        return res

    ga, gb = _gapped_strings(alignment)
    res.aligned_ref = ga
    res.aligned_qry = gb
    res.aligned_columns = len(ga)

    r_pos = 0
    q_pos = 0
    matches = 0
    mism = 0
    ins = 0
    dele = 0
    events: List[DiffEvent] = []
    for i in range(len(ga)):
        cb = gb[i] if i < len(gb) else "-"
        ca = ga[i]
        if ca == "-":
            ins += 1
            if res.first_diff_qry < 0:
                res.first_diff_qry = q_pos
                res.first_diff_ref = r_pos
            if len(events) < 5000:
                events.append(DiffEvent("insertion", r_pos, q_pos + 1, "-", cb))
            q_pos += 1
        elif cb == "-":
            dele += 1
            if res.first_diff_ref < 0:
                res.first_diff_ref = r_pos
            if len(events) < 5000:
                events.append(DiffEvent("deletion", r_pos + 1, q_pos, ca, "-"))
            r_pos += 1
        else:
            if ca == cb:
                matches += 1
            else:
                mism += 1
                if res.first_diff_ref < 0:
                    res.first_diff_ref = r_pos
                    res.first_diff_qry = q_pos
                if len(events) < 5000:
                    events.append(DiffEvent("mismatch", r_pos + 1, q_pos + 1, ca, cb))
            r_pos += 1
            q_pos += 1

    res.matches = matches
    res.mismatches = mism
    res.insertions = ins
    res.deletions = dele
    res.events = events
    res.identity = matches / res.aligned_columns * 100.0 if res.aligned_columns else 0.0
    res.identical = (mism == 0 and ins == 0 and dele == 0)
    if len(events) >= 5000:
        res.notes.append("差异事件超过 5000 条，仅展示前 5000 条")
    return res


def _coarse_compare(ref: str, qry: str, res: Comparison) -> Comparison:
    """超长序列的粗略比较：公共前后缀 + 中段长度差。"""
    res.coarse = True
    n = min(len(ref), len(qry))
    pre = 0
    while pre < n and ref[pre] == qry[pre]:
        pre += 1
    suf = 0
    while suf < n - pre and ref[len(ref) - 1 - suf] == qry[len(qry) - 1 - suf]:
        suf += 1
    mid_ref = ref[pre:len(ref) - suf]
    mid_qry = qry[pre:len(qry) - suf]
    res.matches = pre + suf
    res.aligned_columns = max(len(ref), len(qry))
    res.identity = res.matches / res.aligned_columns * 100.0
    res.first_diff_ref = pre
    res.first_diff_qry = pre
    if len(mid_qry) > len(mid_ref):
        res.insertions = len(mid_qry) - len(mid_ref)
        res.events.append(DiffEvent("insertion", pre, pre + 1, "-",
                                    mid_qry[len(mid_ref):len(mid_qry)][:20]))
    elif len(mid_ref) > len(mid_qry):
        res.deletions = len(mid_ref) - len(mid_qry)
        res.events.append(DiffEvent("deletion", pre + 1, pre, mid_ref[len(mid_qry):][:20], "-"))
    else:
        res.mismatches = sum(1 for a, b in zip(mid_ref, mid_qry) if a != b)
        for i, (a, b) in enumerate(zip(mid_ref, mid_qry)):
            if a != b:
                res.events.append(DiffEvent("mismatch", pre + i + 1, pre + i + 1, a, b))
                break
    res.identical = not (res.mismatches or res.insertions or res.deletions)
    res.notes.append(
        f"序列长度超过 {60_000:,} bp，未做全序列动态规划比对；"
        f"已定位差异区间 ref {pre:,}-{len(ref) - suf:,} / obs {pre:,}-{len(qry) - suf:,}"
    )
    return res


def has_non_acgt(seq: str) -> str:
    return unusual_bases(seq)


# --------------------------------------------------------------------------- #
# 校验结果
# --------------------------------------------------------------------------- #
@dataclass
class CheckResult:
    kind: str                        # 'mRNA' / 'CDS' / 'protein'
    status: str = STATUS_NOT_RUN
    ref_id: str = ""
    ref_token: str = ""
    ref_len: int = 0
    obs_len: int = 0
    comparison: Optional[Comparison] = None
    notes: List[str] = field(default_factory=list)
    candidates: List[Candidate] = field(default_factory=list)
    coding_checks: List[Tuple[str, str, str]] = field(default_factory=list)
    unused_ids: List[str] = field(default_factory=list)
    manual_id: str = ""

    @property
    def status_cn(self) -> str:
        return STATUS_CN.get(self.status, self.status)

    @property
    def color(self) -> str:
        return STATUS_COLORS.get(self.status, "#6b7280")

    def metrics(self) -> List[Tuple[str, str]]:
        c = self.comparison
        rows = [
            ("参考序列", self.ref_id or "—"),
            ("参考长度", f"{self.ref_len:,} bp" if self.ref_len else "—"),
            ("提取长度", f"{self.obs_len:,} bp" if self.obs_len else "—"),
        ]
        if c:
            rows += [
                ("比对列数", f"{c.aligned_columns:,}"),
                ("一致碱基", f"{c.matches:,}"),
                ("不一致碱基", f"{c.mismatches:,}"),
                ("插入（obs 多出）", f"{c.insertions:,}"),
                ("缺失（ref 多出）", f"{c.deletions:,}"),
                ("一致率", f"{c.identity:.4f}%"),
            ]
        return rows

    def headline(self) -> str:
        c = self.comparison
        if self.status == STATUS_IDENTICAL:
            return f"完全一致（{self.ref_len:,} bp）"
        if self.status == STATUS_DIFFERENT and c:
            return (f"存在差异：{c.mismatches} 处替换，{c.insertions} 个插入碱基，"
                    f"{c.deletions} 个缺失碱基，一致率 {c.identity:.3f}%")
        if self.status == STATUS_UNMATCHED:
            return "无法匹配：参考文件中找不到对应记录"
        if self.status == STATUS_MANUAL:
            return "需要人工检查：候选不唯一或置信度不足"
        return "未校验"


def observe_cds_checks(cds: str) -> List[Tuple[str, str, str]]:
    """CDS 层面的阅读框/密码子检查（与参考序列是否存在无关）。"""
    out: List[Tuple[str, str, str]] = []
    if not cds:
        out.append(("info", "无 CDS", "提取结果中没有 CDS"))
        return out
    frame = len(cds) % 3
    out.append(("ok" if frame == 0 else "warn", "阅读框",
                f"{len(cds):,} bp，余 {frame} bp"))
    out.append(("ok" if is_start_codon(cds) else "warn", "起始密码子", first_codon(cds)))
    out.append(("ok" if is_stop_codon(cds) else "warn", "终止密码子", last_codon(cds)))
    stops = internal_stops(cds)
    out.append(("ok" if not stops else "warn", "内部终止密码子", f"{len(stops)} 处"))
    # 翻译产物（截断展示）
    prot = translate(cds)
    out.append(("info", "翻译产物",
                f"{len(prot)} aa，前 40 aa：{prot[:40]}{'…' if len(prot) > 40 else ''}"))
    return out


def validate_extraction(extraction: GeneExtraction,
                        gene: Gene,
                        tx: Transcript,
                        mrna_refs: Optional[RefSet] = None,
                        cds_refs: Optional[RefSet] = None,
                        protein_refs: Optional[RefSet] = None,
                        manual_mrna_id: str = "",
                        manual_cds_id: str = "",
                        manual_protein_id: str = "") -> Dict[str, CheckResult]:
    """对一次提取结果做 mRNA / CDS / 蛋白三重校验。"""
    results: Dict[str, CheckResult] = {}

    # ---- mRNA ----
    mr = CheckResult(kind="mRNA", obs_len=len(extraction.mrna))
    if mrna_refs is None:
        mr.notes.append("未加载参考 mRNA 文件，跳过校验")
        mr.coding_checks = []
    else:
        decision = match_reference(mrna_refs, gene, tx, len(extraction.mrna),
                                  manual_id=manual_mrna_id or None)
        mr.candidates = decision.candidates
        mr.notes.extend(decision.notes)
        if decision.ref_id:
            rec = mrna_refs.get(decision.ref_id)
            if rec is None:
                mr.status = STATUS_UNMATCHED
                mr.notes.append("内部错误：匹配到的参考记录已失效")
            else:
                mr.ref_id = rec.rid
                mr.ref_token = rec.token
                mr.ref_len = len(rec.seq)
                mr.manual_id = manual_mrna_id
                cmp = compare_sequences(rec.seq, extraction.mrna)
                mr.comparison = cmp
                mr.status = STATUS_IDENTICAL if cmp.identical else STATUS_DIFFERENT
                n_acgt = has_non_acgt(rec.seq)
                if n_acgt:
                    mr.notes.append(f"参考序列含非 ACGT 字符：{n_acgt}（按原样比对）")
                mr.unused_ids = [r.token for r in mrna_refs.records
                                 if r.rid != rec.rid][:20]
                if not cmp.identical:
                    mr.notes.append(
                        "排查方向（仅为可能性，不是结论）：注释坐标与参考 mRNA 是否来自同一"
                        "组装版本；外显子边界是否有 1-2 bp 偏移；链方向是否正确；"
                        "Transcript ID 是否指向了不同转录本；是否存在可变剪接导致的外显子差异"
                    )
        else:
            mr.status = decision.status
    results["mRNA"] = mr

    # ---- CDS ----
    cr = CheckResult(kind="CDS", obs_len=len(extraction.cds))
    cr.coding_checks = observe_cds_checks(extraction.cds)
    if cds_refs is None:
        cr.notes.append("未加载参考 CDS 文件，跳过校验")
    elif not extraction.cds:
        cr.status = STATUS_UNMATCHED
        cr.notes.append("该转录本在注释中没有 CDS，无法与参考 CDS 比较")
    else:
        decision = match_reference(cds_refs, gene, tx, len(extraction.cds),
                                  manual_id=manual_cds_id or None)
        cr.candidates = decision.candidates
        cr.notes.extend(decision.notes)
        if decision.ref_id:
            rec = cds_refs.get(decision.ref_id)
            if rec is None:
                cr.status = STATUS_UNMATCHED
            else:
                cr.ref_id = rec.rid
                cr.ref_token = rec.token
                cr.ref_len = len(rec.seq)
                cr.manual_id = manual_cds_id
                cmp = compare_sequences(rec.seq, extraction.cds)
                cr.comparison = cmp
                cr.status = STATUS_IDENTICAL if cmp.identical else STATUS_DIFFERENT
                if cr.status == STATUS_IDENTICAL:
                    cr.notes.append("提取 CDS 与参考 CDS 逐碱基完全一致")
                else:
                    cr.notes.append(
                        "提取 CDS 与参考存在差异；若差异集中在 5′ 端，优先怀疑起始坐标/"
                        "phase 注释问题；若存在插入缺失，需核对是否有额外外显子"
                    )
        else:
            cr.status = decision.status
    results["CDS"] = cr

    # ---- 蛋白（可选） ----
    if protein_refs is not None:
        pr = CheckResult(kind="protein")
        prot = translate(extraction.cds) if extraction.cds else ""
        pr.obs_len = len(prot)
        if not prot:
            pr.status = STATUS_UNMATCHED
            pr.notes.append("没有 CDS，无法翻译出蛋白")
        else:
            # 期望长度用翻译结果（去掉终止符）——参考蛋白 FASTA 通常不写终止符。
            # 以前这里传 0（"蛋白没有 transcript 长度可比较"），等于把长度这条唯一的
            # 非 ID 证据也关掉了：RefSeq 的 protein.faa 标题只有 XP_ 蛋白 ID，不含
            # 转录本 ID，于是永远匹配不上。现在两条路一起走：
            #   ① 从已匹配的参考 CDS 标题里取 protein_id（确定性证据）；
            #   ② 取不到就靠长度打分，把候选摆出来让人确认（绝不硬猜一条）。
            prot_clean = prot.replace("*", "")
            forced = manual_protein_id or ""
            hint_note = ""
            if not forced:
                cds_rec = cds_refs.get(cr.ref_id) if (cds_refs and cr.ref_id) else None
                forced = protein_id_from_cds(cds_rec, protein_refs)
                if forced:
                    hint_note = (f"参考 CDS 标题里写着 protein_id={forced}，"
                                 f"该 ID 在参考蛋白文件中存在，据此确定对应关系")
            decision = match_reference(protein_refs, gene, tx, len(prot_clean),
                                      manual_id=forced or None, reason=hint_note)
            pr.candidates = decision.candidates
            if hint_note:
                pr.notes.insert(0, hint_note)
            pr.notes.extend(n for n in decision.notes if n != hint_note)
            if decision.ref_id:
                rec = protein_refs.get(decision.ref_id)
                if rec:
                    pr.ref_id = rec.rid
                    pr.ref_token = rec.token
                    pr.ref_len = len(rec.seq)
                    cmp = compare_sequences(rec.seq.replace("*", ""),
                                            prot.replace("*", ""))
                    pr.comparison = cmp
                    pr.status = STATUS_IDENTICAL if cmp.identical else STATUS_DIFFERENT
                    pr.notes.append("蛋白序列由提取 CDS 直接翻译，未做任何手工修改")
                else:
                    pr.status = STATUS_UNMATCHED
            else:
                pr.status = decision.status
        results["protein"] = pr

    return results
