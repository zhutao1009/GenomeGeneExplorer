"""生成示例数据集（同时作为测试数据源）。

运行：``python examples/make_example_data.py``
会在本文件所在目录生成：

    demo_genome.fa        合成基因组（chr1 12 kb / chr2 6 kb / chrScaf 3 kb）
    demo_annotation.gff3  覆盖 10 种典型情形的注释
    demo_annotation.gtf   与上面等价的 GTF 版本
    demo_mrna.fa          参考 mRNA（含完全匹配、单碱基缺失、候选二义）
    demo_cds.fa           参考 CDS（含完全匹配、3 bp 插入）
    demo_protein.fa       参考蛋白
    demo_README.txt       数据说明

设计目的：让**每一个**自动化测试都能验证真实碱基内容，而不是只验证
"程序没崩"。因此每条序列都由固定随机种子生成，坐标与参考序列都可复现。
"""

from __future__ import annotations

import os
import random
import sys

SEED = 20261028
CHROMS = [("chr1", 12000), ("chr2", 6000), ("chrScaf", 3000)]

# --------------------------------------------------------------------------- #
# 基因定义
# --------------------------------------------------------------------------- #
# 每个基因：id, name, seqid, strand, gene_span, exons, cds
GENES = [
    dict(gid="g1", name="GeneOne", seqid="chr1", strand="+", span=(1001, 1800),
         exons=[(1001, 1200), (1401, 1550), (1701, 1800)],
         cds=[(1051, 1200), (1401, 1550), (1701, 1760)],
         biotype="protein_coding",
         desc="正链、3 外显子、含 5′/3′ UTR；另有可变剪接转录本 g1.2"),
    dict(gid="g1", name="GeneOne", seqid="chr1", strand="+", span=(1001, 1800),
         exons=[(1001, 1200), (1701, 1800)],
         cds=[(1051, 1200), (1701, 1760)],
         tx="g1.2", biotype="protein_coding",
         desc="g1 的第二转录本（跳过外显子 2）"),
    dict(gid="g2", name="GeneTwo", seqid="chr1", strand="-", span=(4001, 5000),
         exons=[(4001, 4200), (4401, 4600), (4801, 5000)],
         cds=[(4001, 4103), (4401, 4600), (4851, 5000)],
         biotype="protein_coding",
         desc="负链、3 外显子，用于验证反向互补与片段倒序"),
    dict(gid="g3", name="GeneThree", seqid="chr2", strand="+", span=(501, 900),
         exons=[(501, 900)], cds=[(501, 860)],
         biotype="protein_coding",
         desc="单外显子（应报告无内含子）"),
    dict(gid="g4", name="GeneFour", seqid="chrScaf", strand="-", span=(51, 400),
         exons=[(51, 150), (251, 400)], cds=[(51, 150), (251, 400)],
         biotype="protein_coding",
         desc="贴近 scaffold 起点，上游只能取到 50 bp（边界截断）"),
    dict(gid="g5", name="GeneFive", seqid="chr1", strand="+", span=(6001, 6400),
         exons=[], cds=[(6001, 6120), (6301, 6400)],
         biotype="protein_coding",
         desc="注释只有 CDS、没有 exon（外显子需由 CDS 推导）"),
    dict(gid="g6", name="GeneSix", seqid="chr2", strand="+", span=(2001, 2200),
         exons=[(2001, 2200)], cds=[],
         biotype="lncRNA",
         desc="非编码转录本（无 CDS，不应报错）"),
    dict(gid="g7", name="GeneSeven", seqid="chr1", strand="+", span=(7001, 7400),
         exons=[(7001, 7100), (7300, 7350)], cds=[(7001, 7100), (7300, 7350)],
         biotype="protein_coding", bad_exon=True,
         desc="含 start>end 的异常 exon 行（解析器应交换并警告）"),
    dict(gid="g8", name="GeneEight", seqid="chrUn_1", strand="+", span=(1, 500),
         exons=[(1, 500)], cds=[(1, 300)],
         biotype="protein_coding",
         desc="位于 FASTA 中不存在的序列 chrUn_1（用于验证 ID 不匹配报错）"),
    dict(gid="g9", name="GeneNine", seqid="chr2", strand="+", span=(2501, 2700),
         exons=[(2501, 2600), (2550, 2700)], cds=[(2501, 2600)],
         biotype="protein_coding",
         desc="两个外显子在基因组上重叠（注释异常，不应造出负长度内含子）"),
    dict(gid="g10", name="GeneTen", seqid="chr1", strand="+", span=(0, 0),
         exons=[(8001, 8100), (8150, 8200)], cds=[(8001, 8100)],
         biotype="protein_coding", no_gene_row=True,
         desc="注释缺少 gene 行（由转录本反推基因区间）"),
]


def _random_genome() -> dict:
    rng = random.Random(SEED)
    return {name: "".join(rng.choice("ACGT") for _ in range(n)) for name, n in CHROMS}


def _patch(seq: str, start1: int, new: str) -> str:
    """把 1-based 位置 start1 起始的碱基替换为 new。"""
    i = start1 - 1
    return seq[:i] + new + seq[i + len(new):]


def _rc(seq: str) -> str:
    table = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(table)[::-1]


def _apply_codon_patches(genome: dict) -> None:
    """让每个 CDS 以 ATG 起始、以 TAA 结束（仅改动必要的 3 nt）。"""
    for g in GENES:
        if not g["cds"]:
            continue
        seqid, strand = g["seqid"], g["strand"]
        seq = genome.get(seqid)
        if seq is None:
            continue
        exons = g["exons"] or g["cds"]
        # 起始密码子
        first = g["cds"][-1] if strand == "-" else g["cds"][0]
        if strand == "+":
            start1 = first[0]
            seq = _patch(seq, start1, "ATG")
        else:
            start1 = first[1] - 2
            seq = _patch(seq, start1, _rc("ATG"))
        genome[seqid] = seq

        # 终止密码子
        last = g["cds"][0] if strand == "-" else g["cds"][-1]
        cds_len = sum(e - s + 1 for s, e in g["cds"])
        if cds_len % 3 != 0:
            continue
        if strand == "+":
            stop1 = last[1] - 2
            seq = _patch(genome[seqid], stop1, "TAA")
        else:
            stop1 = last[0]
            seq = _patch(genome[seqid], stop1, _rc("TAA"))
        genome[seqid] = seq


def expected_for(g: dict, genome: dict) -> dict:
    """独立于 gge 包，用最直白的方式算出期望序列。"""
    seqid, strand = g["seqid"], g["strand"]
    seq = genome.get(seqid, "")
    exons = g["exons"] or g["cds"]
    ordered = sorted(exons, reverse=(strand == "-"))
    exon_seqs = []
    for s, e in ordered:
        piece = seq[s - 1:e]
        exon_seqs.append(_rc(piece) if strand == "-" else piece)
    mrna = "".join(exon_seqs)
    cds_ordered = sorted(g["cds"], reverse=(strand == "-"))
    cds_seq = ""
    for s, e in cds_ordered:
        piece = seq[s - 1:e] if seq else ""
        cds_seq += _rc(piece) if strand == "-" else piece
    return {
        "exon_order": ordered,
        "exon_seqs": exon_seqs,
        "mrna": mrna,
        "cds": cds_seq,
        "seqid": seqid,
        "strand": strand,
    }


def _gff3(g: dict, tx_id: str) -> list:
    """返回该转录本的 GFF3 行。"""
    out = []
    seqid, strand = g["seqid"], g["strand"]
    if not g.get("no_gene_row"):
        out.append(
            f"{seqid}\tdemo\tgene\t{g['span'][0]}\t{g['span'][1]}\t.\t{strand}\t.\t"
            f"ID={g['gid']};Name={g['name']};biotype={g['biotype']}"
        )
    out.append(
        f"{seqid}\tdemo\tmRNA\t{g['span'][0] if g['span'][0] else min(s for s, e in (g['exons'] or g['cds']))}"
        f"\t{g['span'][1] if g['span'][1] else max(e for s, e in (g['exons'] or g['cds']))}"
        f"\t.\t{strand}\t.\tID={tx_id};Parent={g['gid']};biotype={g['biotype']}"
    )
    if g.get("bad_exon") and len(g["exons"]) > 1:
        # 第一个外显子正常，第二个故意写成 start > end
        s0, e0 = g["exons"][0]
        out.append(f"{seqid}\tdemo\texon\t{s0}\t{e0}\t.\t{strand}\t.\tID=exon:{tx_id}:1;Parent={tx_id}")
        s, e = g["exons"][1]
        out.append(f"{seqid}\tdemo\texon\t{e}\t{s}\t.\t{strand}\t.\tID=exon:{tx_id}:2;Parent={tx_id}")
    else:
        for i, (s, e) in enumerate(sorted(g["exons"]), start=1):
            out.append(f"{seqid}\tdemo\texon\t{s}\t{e}\t.\t{strand}\t.\tID=exon:{tx_id}:{i};Parent={tx_id}")
    for i, (s, e) in enumerate(sorted(g["cds"]), start=1):
        out.append(f"{seqid}\tdemo\tCDS\t{s}\t{e}\t.\t{strand}\t0\tID=cds:{tx_id}:{i};Parent={tx_id}")
    return out


def _gtf_line(seqid, ftype, start, end, strand, attrs: str) -> str:
    return f"{seqid}\tdemo\t{ftype}\t{start}\t{end}\t.\t{strand}\t.\t{attrs}"


def build_example_data(outdir: str, verbose: bool = False) -> dict:
    genome = _random_genome()
    _apply_codon_patches(genome)

    os.makedirs(outdir, exist_ok=True)
    fa_path = os.path.join(outdir, "demo_genome.fa")
    with open(fa_path, "w", encoding="utf-8", newline="\n") as fh:
        for name, _n in CHROMS:
            fh.write(f">{name} synthetic demo chromosome\n")
            s = genome[name]
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")

    # ---------------- GFF3 ---------------- #
    gff_lines = ["##gff-version 3", "##sequence-region chr1 1 12000"]
    tx_ids = {}
    for g in GENES:
        if "tx" in g:
            tx = g["tx"]
        else:
            tx = f"{g['gid']}.1"
        if g["gid"] in tx_ids:
            tx_ids[g["gid"]].append(tx)
        else:
            tx_ids[g["gid"]] = [tx]
        gff_lines.extend(_gff3(g, tx))
    gff_path = os.path.join(outdir, "demo_annotation.gff3")
    with open(gff_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(gff_lines) + "\n")

    # ---------------- GTF ---------------- #
    gtf_lines = ["#!genome-build demo"]
    for g in GENES:
        tx = g.get("tx", f"{g['gid']}.1")
        if not g.get("no_gene_row"):
            gtf_lines.append(_gtf_line(
                g["seqid"], "gene", g["span"][0], g["span"][1], g["strand"],
                f'gene_id "{g["gid"]}"; gene_name "{g["name"]}"; gene_biotype "{g["biotype"]}";'))
        span = g["span"]
        t_start = span[0] or min(s for s, e in (g["exons"] or g["cds"]))
        t_end = span[1] or max(e for s, e in (g["exons"] or g["cds"]))
        gtf_lines.append(_gtf_line(
            g["seqid"], "transcript", t_start, t_end, g["strand"],
            f'gene_id "{g["gid"]}"; transcript_id "{tx}"; gene_name "{g["name"]}";'))
        if g.get("bad_exon") and len(g["exons"]) > 1:
            s0, e0 = g["exons"][0]
            gtf_lines.append(_gtf_line(g["seqid"], "exon", s0, e0, g["strand"],
                                       f'gene_id "{g["gid"]}"; transcript_id "{tx}"; exon_number "1";'))
            s, e = g["exons"][1]
            gtf_lines.append(_gtf_line(g["seqid"], "exon", e, s, g["strand"],
                                       f'gene_id "{g["gid"]}"; transcript_id "{tx}"; exon_number "2";'))
        else:
            for i, (s, e) in enumerate(sorted(g["exons"]), start=1):
                gtf_lines.append(_gtf_line(g["seqid"], "exon", s, e, g["strand"],
                                           f'gene_id "{g["gid"]}"; transcript_id "{tx}"; exon_number "{i}";'))
        for i, (s, e) in enumerate(sorted(g["cds"]), start=1):
            phase = 0 if i == 1 else "."
            gtf_lines.append(_gtf_line(g["seqid"], "CDS", s, e, g["strand"],
                                       f'gene_id "{g["gid"]}"; transcript_id "{tx}"; frame {phase};'))
    gtf_path = os.path.join(outdir, "demo_annotation.gtf")
    with open(gtf_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(gtf_lines) + "\n")

    # ---------------- 期望值 ---------------- #
    exp = {}
    for g in GENES:
        tx = g.get("tx", f"{g['gid']}.1")
        e = expected_for(g, genome)
        e["tx"] = tx
        e["desc"] = g["desc"]
        exp.setdefault(g["gid"], {})[tx] = e

    # ---------------- 参考 mRNA ---------------- #
    mrna_records = []
    primary = exp["g1"]["g1.1"]
    mrna_records.append(("g1.1 g1 primary mRNA (exact)", primary["mrna"]))
    mrna_records.append(("g1.2 g1 alternative mRNA (exact)", exp["g1"]["g1.2"]["mrna"]))

    # g2：单碱基缺失
    g2 = exp["g2"]["g2.1"]["mrna"]
    pos = len(g2) // 2
    mrna_records.append(("g2 GeneTwo mRNA (1 bp deleted at %d)" % (pos + 1),
                         g2[:pos] + g2[pos + 1:]))

    # g3：单碱基替换
    g3 = exp["g3"]["g3.1"]["mrna"]
    pos3 = 100
    alt = {"A": "C", "C": "G", "G": "T", "T": "A"}[g3[pos3]]
    mrna_records.append(("g3 GeneThree mRNA (1 substitution at %d)" % (pos3 + 1),
                         g3[:pos3] + alt + g3[pos3 + 1:]))

    # g6：两条等长候选，制造"需要人工检查"
    g6 = exp["g6"]["g6.1"]["mrna"]
    mrna_records.append(("g6 copyA", g6))
    mrna_records.append(("g6 copyB", "A" * len(g6)))
    mrna_records.append(("unrelated_sequence no_gene_match", "ACGT" * 50))

    mrna_path = os.path.join(outdir, "demo_mrna.fa")
    with open(mrna_path, "w", encoding="utf-8", newline="\n") as fh:
        for header, seq in mrna_records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")

    # ---------------- 参考 CDS ---------------- #
    cds_records = [
        ("g1.1_cds exact", exp["g1"]["g1.1"]["cds"]),
        ("g2.1_cds exact", exp["g2"]["g2.1"]["cds"]),
        ("g3.1_cds with 3bp insertion",
         exp["g3"]["g3.1"]["cds"][:60] + "GGG" + exp["g3"]["g3.1"]["cds"][60:]),
    ]
    cds_path = os.path.join(outdir, "demo_cds.fa")
    with open(cds_path, "w", encoding="utf-8", newline="\n") as fh:
        for header, seq in cds_records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")

    # ---------------- 参考蛋白 ---------------- #
    codon_table = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
        "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
        "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
        "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
        "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R",
        "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    }
    cds1 = exp["g1"]["g1.1"]["cds"]
    prot = "".join(codon_table.get(cds1[i:i + 3], "X") for i in range(0, len(cds1) - 2, 3))
    prot_path = os.path.join(outdir, "demo_protein.fa")
    with open(prot_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(">g1.1_protein\n")
        for i in range(0, len(prot), 60):
            fh.write(prot[i:i + 60] + "\n")

    # ---------------- 说明 ---------------- #
    readme = os.path.join(outdir, "demo_README.txt")
    lines = [
        "示例数据说明 / Demo dataset",
        "=" * 60,
        f"随机种子 {SEED}，基因组完全合成，不含任何真实物种序列。",
        f"染色体长度：{', '.join(f'{n}={l}' for n, l in CHROMS)}",
        "",
        "每个基因覆盖的情形：",
    ]
    for g in GENES:
        tx = g.get("tx", f"{g['gid']}.1")
        lines.append(f"  {g['gid']:4s} {tx:8s} {g['seqid']:8s}{g['strand']}  {g['desc']}")
    lines += [
        "",
        "参考 mRNA 文件里 g2 故意少 1 个碱基、g3 故意替换 1 个碱基；",
        "g6 有两条等长候选，用来触发『需要人工检查』。",
        "参考 CDS 文件里 g3 故意插入 3 bp，用于验证 indel 定位。",
        "",
        "原始坐标：GFF3/GTF 为 1-based 闭区间。",
    ]
    with open(readme, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")

    if verbose:
        print(f"示例数据已生成于 {outdir}")
        print(f"  基因 {len(GENES)} 个 / 转录本 {sum(len(v) for v in tx_ids.values())} 个")
    return {
        "genome": fa_path,
        "gff3": gff_path,
        "gtf": gtf_path,
        "mrna": mrna_path,
        "cds": cds_path,
        "protein": prot_path,
        "genome_sequences": genome,
        "expected": exp,
        "chroms": dict(CHROMS),
    }


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    build_example_data(here, verbose=True)
