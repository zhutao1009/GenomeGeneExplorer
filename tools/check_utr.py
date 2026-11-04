# -*- coding: utf-8 -*-
"""核查 UTR 推导是否正确（真实数据实测）。

用法：python tools/check_utr.py gdf9
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gge.tasks import load_genome, load_annotation_file
from gge.extract import extract_gene, PANEL_UTR
from gge.models import REGION_UTR5, REGION_UTR3

GENOME = r"D:\grasscarp_gene_extract\GCF_019924925.1_HZGC01_genomic.fna"
GFF = r"D:\grasscarp_gene_extract\GCF_019924925.1_HZGC01_genomic.gff"


def independent_utr(tx):
    """完全不依赖程序逻辑：对每个 exon，减去与 CDS 的交集，剩下的就是 UTR。

    返回 [(kind, start, end)]，start/end 为 1-based 闭区间基因组坐标。
    """
    out = []
    cds = tx.ordered_cds()
    cds_first = min(c.start for c in cds)
    cds_last = max(c.end for c in cds)
    for ex in tx.ordered_exons():
        # 该 exon 内被 CDS 覆盖的区间
        covered = [(max(c.start, ex.start), min(c.end, ex.end)) for c in cds
                   if c.start <= ex.end and c.end >= ex.start]
        # exon 减去这些覆盖区间
        pieces = [(ex.start, ex.end)]
        for cs, ce in covered:
            nxt = []
            for s, e in pieces:
                if ce < s or cs > e:
                    nxt.append((s, e)); continue
                if s < cs:
                    nxt.append((s, cs - 1))
                if e > ce:
                    nxt.append((ce + 1, e))
            pieces = nxt
        for s, e in pieces:
            if s > e:
                continue
            if cds_first <= s <= cds_last:
                # 落在 CDS 跨度之内却没被 CDS 覆盖 —— 内含子样区域，不该出现
                kind = "???"
            elif e < cds_first:
                kind = REGION_UTR5
            else:
                kind = REGION_UTR3
            out.append((kind, s, e))
    return out


def main():
    name = sys.argv[1] if len(sys.argv) > 1 else "gdf9"
    print("加载基因组与注释（约 70 s）…", flush=True)
    info = load_genome(GENOME)
    model = load_annotation_file(GFF)
    provider = info.provider

    genes = [g for gid, g in model.genes.items()
             if gid == f"gene-{name}" or g.gene_name == name]
    if not genes:
        print("未找到该基因"); return
    gene = genes[0]
    genome = provider.fetch(gene.seqid, 1, provider.length(gene.seqid))

    print(f"\n基因 {gene.gene_id} ({gene.gene_name})  "
          f"{gene.seqid}:{gene.start:,}-{gene.end:,} ({gene.strand})  "
          f"转录本 {gene.n_transcripts} 个")
    for t in gene.tx_list():
        cds = t.ordered_cds()
        print(f"  {t.tx_id}  外显子 {len(t.exons)}  CDS片段 {len(cds)}  "
              f"CDS总长 {sum(c.length for c in cds)}")

    tx = gene.canonical_transcript()
    print(f"\n默认（canonical）转录本：{tx.tx_id}")

    # 独立计算
    ind = independent_utr(tx)
    print("\n--- 独立计算结果（exon 减去 CDS 的交集）---")
    for kind, s, e in ind:
        print(f"  {kind:9s} {gene.seqid}:{s:,}-{e:,}  {e - s + 1:,} bp")

    # 程序输出
    ex = extract_gene(provider, gene, tx, 2000, 2000)
    print("\n--- 程序 UTR 面板 ---")
    panel = ex.panels[PANEL_UTR]
    for b in panel.blocks:
        print(f"  {b.title:8s} {b.note}")
    print(f"  汇总：5'UTR {len(ex.utr5):,} bp / 3'UTR {len(ex.utr3):,} bp")

    # 逐条比对坐标与序列
    print("\n--- 一致性核对 ---")
    got = []
    for b in panel.blocks:
        gs, ge = int(b.meta["genomic_start"]), int(b.meta["genomic_end"])
        got.append((b.meta["utr_kind"], gs, ge, b.seq))
    ok = True
    for kind, s, e in ind:
        hit = [g for g in got if g[0] == kind and g[1] == s and g[2] == e]
        if not hit:
            print(f"  [不一致] {kind} {s}-{e} 在程序输出里找不到"); ok = False
            continue
        seq_prog = hit[0][3]
        seq_direct = genome[s - 1:e]          # 正链基因直接用切片
        if tx.strand == "-":
            from gge.coords import revcomp
            seq_direct = revcomp(seq_direct)
        same = seq_prog == seq_direct
        print(f"  {'OK ' if same else 'ERR'} {kind} {s:,}-{e:,}  "
              f"程序 {len(seq_prog)} bp / 直接切片 {len(seq_direct)} bp  序列相同={same}")
        ok &= same
    print(f"\n总体：{'全部一致' if ok else '存在不一致'}")

    # 终止密码子与 CDS 末端
    cds = tx.ordered_cds()
    last = max(cds, key=lambda c: c.end) if tx.strand == "+" else min(cds, key=lambda c: c.start)
    cds_last_bp = last.end if tx.strand == "+" else last.start
    codon = genome[cds_last_bp - 3:cds_last_bp]
    if tx.strand == "-":
        from gge.coords import revcomp
        codon = revcomp(codon)
    print(f"\nCDS 末端 3 bp（{gene.seqid}:{cds_last_bp - 2:,}-{cds_last_bp:,}）：{codon}  "
          f"→ {'终止密码子' if codon.upper() in ('TAA', 'TAG', 'TGA') else '不是终止密码子'}")
    print(f"CDS 之后紧接着的就是 3′UTR 起点 {cds_last_bp + 1:,}")

    # 拼接检验：exon2 序列 == CDS2 + 3'UTR
    if len(tx.ordered_exons()) >= 2:
        ex2 = tx.ordered_exons()[-1]
        print(f"\n检验：最后一个外显子 {ex2.start:,}-{ex2.end:,} 长 {ex2.length:,} bp；"
              f"CDS 部分 {ex2.length - len(ex.utr3):,} bp + 3′UTR {len(ex.utr3):,} bp "
              f"= {ex2.length - len(ex.utr3) + len(ex.utr3):,} bp "
              f"{'（吻合）' if ex2.length == (ex2.length - len(ex.utr3)) + len(ex.utr3) else ''}")


if __name__ == "__main__":
    main()
