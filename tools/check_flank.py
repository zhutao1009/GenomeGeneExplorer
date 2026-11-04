# -*- coding: utf-8 -*-
"""核查上下游侧翼序列与基因组连续窗口（真实数据实测）。

用法：
    python tools/check_flank.py prlrb 2000
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gge.tasks import load_genome, load_annotation_file
from gge.extract import extract_gene, PANEL_UPSTREAM, PANEL_DOWNSTREAM, PANEL_GENOMIC

GENOME = r"D:\grasscarp_gene_extract\GCF_019924925.1_HZGC01_genomic.fna"
GFF = r"D:\grasscarp_gene_extract\GCF_019924925.1_HZGC01_genomic.gff"


def main():
    name = sys.argv[1] if len(sys.argv) > 1 else "prlrb"
    req = int(sys.argv[2]) if len(sys.argv) > 2 else 2000

    t0 = time.time()
    print(f"加载基因组 … {os.path.basename(GENOME)}", flush=True)
    info = load_genome(GENOME)
    provider = info.provider
    print(f"  完成 {time.time()-t0:.1f}s，序列数 {len(info.seqids)}", flush=True)

    t1 = time.time()
    print(f"加载注释 … {os.path.basename(GFF)}", flush=True)
    model = load_annotation_file(GFF)
    print(f"  完成 {time.time()-t1:.1f}s，基因数 {len(model.genes):,}", flush=True)

    # 按 gene_id 与 gene name 两种方式找
    cands = [g for gid, g in model.genes.items()
             if gid == f"gene-{name}" or g.gene_name == name
             or gid.endswith(f"-{name}")]
    if not cands:
        cands = model.gene_by_name(name)
    print(f"\n命中基因 {len(cands)} 个：", flush=True)
    for g in cands[:5]:
        print(f"  {g.gene_id}  name={g.gene_name}  {g.seqid}:{g.start:,}-{g.end:,} "
              f"({g.strand})  转录本 {g.n_transcripts}", flush=True)
    if not cands:
        print("未找到该基因"); return

    gene = cands[0]
    tx = gene.canonical_transcript() or gene.tx_list()[0]
    print(f"\n使用转录本：{tx.tx_id}  {tx.seqid}:{tx.span()[0]:,}-{tx.span()[1]:,} "
          f"({tx.strand})  外显子 {len(tx.exons)}", flush=True)

    print(f"\n提取：上游 {req:,} bp / 下游 {req:,} bp", flush=True)
    t2 = time.time()
    ex = extract_gene(provider, gene, tx, req, req)
    print(f"  完成 {time.time()-t2:.1f}s", flush=True)

    print("\n--- 结果 ---")
    print(f"seqid          : {ex.seqid}  (长度 {provider.length(ex.seqid):,})")
    print(f"strand         : {ex.strand}")
    print(f"gene span      : {ex.gene_start:,}-{ex.gene_end:,}")
    print(f"region         : {ex.region_start:,}-{ex.region_end:,}  "
          f"长度 {ex.region_end - ex.region_start + 1:,}")
    print(f"上游 requested : {ex.upstream_requested:,}")
    print(f"上游 actual    : {ex.upstream_actual:,}   clipped={ex.upstream_clipped}")
    print(f"下游 requested : {ex.downstream_requested:,}")
    print(f"下游 actual    : {ex.downstream_actual:,}   clipped={ex.downstream_clipped}")

    up = ex.panels.get(PANEL_UPSTREAM)
    dn = ex.panels.get(PANEL_DOWNSTREAM)
    gm = ex.panels.get(PANEL_GENOMIC)
    print("\n--- 面板 ---")
    print(f"upstream  面板: {len(up.blocks) if up else 0} 块  "
          f"长度 {sum(len(b.seq) for b in up.blocks) if up else 0:,}  "
          f"empty_text={getattr(up, 'empty_text', '')!r}")
    print(f"downstream面板: {len(dn.blocks) if dn else 0} 块  "
          f"长度 {sum(len(b.seq) for b in dn.blocks) if dn else 0:,}  "
          f"empty_text={getattr(dn, 'empty_text', '')!r}")
    print(f"genomic   面板: {len(gm.blocks) if gm else 0} 块  "
          f"长度 {sum(len(b.seq) for b in gm.blocks) if gm else 0:,}")

    if up and up.blocks:
        b = up.blocks[0]
        print(f"\n上游块 note: {b.note}")
        print(f"上游前 60 bp: {b.seq[:60]}")
        print(f"上游后 60 bp: {b.seq[-60:]}")
    if dn and dn.blocks:
        b = dn.blocks[0]
        print(f"\n下游块 note: {b.note}")
        print(f"下游前 60 bp: {b.seq[:60]}")
        print(f"下游后 60 bp: {b.seq[-60:]}")

    if gm and gm.blocks:
        b = gm.blocks[0]
        print(f"\n基因组窗口 note: {b.note}")
        print(f"窗口长度 {len(b.seq):,}  = 上游 {ex.upstream_actual:,} + "
              f"基因 {ex.gene_end - ex.gene_start + 1:,} + 下游 {ex.downstream_actual:,}")
        # 校验窗口里基因两端能对上
        gs_off = ex.gene_start - ex.region_start
        print(f"窗口中基因起点偏移 {gs_off:,}；该处序列 {b.seq[max(0,gs_off-10):gs_off+10]}")

    # 独立复核：直接用 pyfaidx 风格取序列比对
    print("\n--- 独立复核（直接从 FASTA 取）---")
    ts, te = tx.span()
    if tx.strand == "+":
        up_raw = (ts - req, ts - 1)
        dn_raw = (te + 1, te + req)
    else:
        up_raw = (te + 1, te + req)
        dn_raw = (ts - req, ts - 1)
    print(f"上游应为 {ex.seqid}:{up_raw[0]:,}-{up_raw[1]:,}（长 {up_raw[1]-up_raw[0]+1:,}）")
    print(f"下游应为 {ex.seqid}:{dn_raw[0]:,}-{dn_raw[1]:,}（长 {dn_raw[1]-dn_raw[0]+1:,}）")
    raw_up = provider.fetch(ex.seqid, up_raw[0], up_raw[1])
    raw_dn = provider.fetch(ex.seqid, dn_raw[0], dn_raw[1])
    if tx.strand == "-":          # 负链面板存的是转录方向的序列，需要反向互补后再比
        from gge.coords import revcomp
        raw_up, raw_dn = revcomp(raw_up), revcomp(raw_dn)
    print(f"直接取到 上游 {len(raw_up):,} bp / 下游 {len(raw_dn):,} bp")
    print(f"与结果一致：上游 {raw_up == (ex.upstream_seq or '')}  "
          f"下游 {raw_dn == (ex.downstream_seq or '')}")


if __name__ == "__main__":
    main()
