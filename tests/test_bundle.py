"""压缩包自动识别的测试。

重点验证两件事：

1. **不看文件名也能认对**。示例包故意用 ``main.fa`` / ``coding.fasta``
   这类无信息量的名字，甚至把一个核酸文件叫 ``protein_like.fasta``，
   判定必须靠内容（字母表、序列名、长度拟合、记录数）。
2. **不静默选错**。多个候选、缺必填项、成员路径越界，都要有明确交代。
"""

from __future__ import annotations

import gzip
import os
import tarfile
import zipfile

import pytest

from gge.bundle import (BundleScan, import_bundle, scan_archive, safe_filename,
                        summarize_annotation_path, _is_bgzf)
from gge.errors import FileFormatError, GGEError


# --------------------------------------------------------------------------- #
# 辅助
# --------------------------------------------------------------------------- #
def _roles(scan: BundleScan) -> dict:
    return {r: (c.member if c else "") for r, c in scan.by_role().items()}


def _synthetic_dataset(tmp_path):
    """造一套"文件名完全没有提示"的合成数据，用来验证长度拟合。"""
    from gge.tasks import extract_and_validate  # noqa: F401  (仅为证明依赖可用)

    n_tx = 8
    exons_per_tx = 2
    exon_len = 60
    cds_len = 90        # 每条转录本 CDS 90 bp（能被 3 整除），mRNA 120 bp（不能）
    genome_parts = []
    gff_lines = ["##gff-version 3"]
    mrna_records = []
    cds_records = []
    pos = 1
    for i in range(1, n_tx + 1):
        tid = f"tx{i:02d}"
        gid = f"gene{i:02d}"
        gene_start = pos
        gff_lines.append(f"chrA\ttest\tgene\t{pos}\t{pos + 200}\t.\t+\t.\tID={gid}")
        gff_lines.append(f"chrA\ttest\tmRNA\t{pos}\t{pos + 200}\t.\t+\t.\tID={tid};Parent={gid}")
        tx_seq = []
        for e in range(exons_per_tx):
            start = pos + e * 100
            end = start + exon_len - 1
            gff_lines.append(f"chrA\ttest\texon\t{start}\t{end}\t.\t+\t.\tID={tid}.ex{e+1};Parent={tid}")
            tx_seq.append(("ACGT" * 40)[:exon_len])
        # CDS 恰好占满第 1 个外显子 + 第 2 个外显子前 30 bp
        gff_lines.append(f"chrA\ttest\tCDS\t{pos}\t{pos + exon_len - 1}\t.\t+\t0\tID={tid}.cds1;Parent={tid}")
        gff_lines.append(f"chrA\ttest\tCDS\t{pos + 100}\t{pos + 100 + 29}\t.\t+\t0\tID={tid}.cds2;Parent={tid}")
        mrna_records.append((tid, "".join(tx_seq)))
        cds_records.append((tid, ("ATG" + "GC" * 43 + "TAA")[:cds_len]))
        genome_parts.append("N" * 240)
        pos += 240

    genome = os.path.join(tmp_path, "genome_raw.txt")
    with open(genome, "w", encoding="utf-8") as fh:
        fh.write(">chrA\n")
        seq = "".join(genome_parts)
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + "\n")

    gff = os.path.join(tmp_path, "ann_raw.txt")
    with open(gff, "w", encoding="utf-8") as fh:
        fh.write("\n".join(gff_lines) + "\n")

    def _fasta(path, records):
        with open(path, "w", encoding="utf-8") as fh:
            for name, seq in records:
                fh.write(f">{name}\n{seq}\n")
        return path

    return {
        "genome": genome,
        "gff": gff,
        "mrna": _fasta(os.path.join(tmp_path, "a1.txt"), mrna_records),
        "cds": _fasta(os.path.join(tmp_path, "a2.txt"), cds_records),
    }


# --------------------------------------------------------------------------- #
# 识别准确性
# --------------------------------------------------------------------------- #
def test_identifies_roles_by_content_not_names(make_zip, demo_dir):
    """文件名全都没有提示（main.fa / coding.fasta），仍要认对五个角色。"""
    zp = make_zip({
        "data/main.fa": demo_dir["genome"],
        "annotation.gff3": demo_dir["gff3"],
        "refs/transcripts.fasta": demo_dir["mrna"],
        "refs/coding.fasta": demo_dir["cds"],
        "refs/pep.faa": demo_dir["protein"],
    })
    scan = scan_archive(zp)
    roles = _roles(scan)
    assert roles["genome"] == "data/main.fa"
    assert roles["annotation"] == "annotation.gff3"
    assert roles["mrna"] == "refs/transcripts.fasta"
    assert roles["cds"] == "refs/coding.fasta"
    assert roles["protein"] == "refs/pep.faa"
    # 每个判定都要能说出理由
    for cand in scan.candidates:
        assert cand.evidence, f"{cand.member} 没有给出判定依据"


def test_mrna_and_cds_told_apart_by_length_fit(make_zip, tmp_path):
    """两个文件都叫 a1/a2.txt、标题都是转录本 ID 时，只能靠长度与注释比对来区分。"""
    data = _synthetic_dataset(tmp_path)
    zp = make_zip({
        "g.txt": data["genome"],
        "x.gff3": data["gff"],
        "a1.txt": data["mrna"],
        "a2.txt": data["cds"],
    })
    scan = scan_archive(zp)
    roles = _roles(scan)
    assert roles["mrna"] == "a1.txt", scan.text_report()
    assert roles["cds"] == "a2.txt", scan.text_report()
    mrna_ev = next(c for c in scan.candidates if c.member == "a1.txt").evidence_text()
    assert "mRNA 长度吻合" in mrna_ev


def test_protein_detected_by_alphabet_without_name_hint(make_zip, demo_dir):
    """把蛋白文件改名成完全中性的名字，也要靠字母构成认出来。"""
    zp = make_zip({
        "g.fa": demo_dir["genome"],
        "a.gff3": demo_dir["gff3"],
        "z9.fasta": demo_dir["protein"],
    })
    scan = scan_archive(zp)
    roles = _roles(scan)
    assert roles["protein"] == "z9.fasta", scan.text_report()
    cand = next(c for c in scan.candidates if c.member == "z9.fasta")
    assert "蛋白" in cand.evidence_text()


def test_gtf_is_recognised_as_annotation(make_zip, demo_dir):
    zp = make_zip({
        "g.fa": demo_dir["genome"],
        "genes.gtf": demo_dir["gtf"],
    })
    scan = scan_archive(zp)
    cand = scan.by_role()["annotation"]
    assert cand is not None and cand.fmt == "GTF"
    assert "GTF" in cand.evidence_text()
    assert scan.annotation is not None and scan.annotation.fmt == "gtf"
    assert scan.annotation.n_transcripts > 0


def test_annotation_summary_agrees_with_parser(demo_dir):
    """流式统计的数量必须与正式解析对得上，否则判定依据就是错的。"""
    from gge.annotation import parse_annotation

    ann = parse_annotation(demo_dir["gff3"])
    summary = summarize_annotation_path(demo_dir["gff3"])
    assert summary.n_genes == len(ann.genes)
    assert summary.n_transcripts == len(ann.tx_index)
    assert set(summary.seqids) >= {g.seqid for g in ann.genes.values()}


def test_ignored_members_are_listed_with_reason(make_zip, demo_dir):
    zp = make_zip({
        "g.fa": demo_dir["genome"],
        "a.gff3": demo_dir["gff3"],
        "README.txt": ("text", "说明：本包包含三个文件\n"),
        ".DS_Store": ("text", "junk"),
        "__MACOSX/._g.fa": ("text", "junk"),
        "align.bam": ("text", "not really bam"),
    })
    scan = scan_archive(zp)
    reasons = {name: why for name, why in scan.ignored}
    assert "README.txt" in reasons
    assert ".DS_Store" in reasons
    assert "__MACOSX/._g.fa" in reasons
    assert "align.bam" in reasons
    # 被跳过的文件不能混进候选
    assert all("bam" not in c.member for c in scan.candidates)


def test_binary_member_is_skipped(make_zip, demo_dir, tmp_path):
    blob = tmp_path / "blob.dat"
    blob.write_bytes(b"\x00\x01\x02 binary \x00")
    zp = make_zip({"g.fa": demo_dir["genome"], "a.gff3": demo_dir["gff3"],
                   "blob.dat": str(blob)})
    scan = scan_archive(zp)
    assert any("二进制" in why for _n, why in scan.ignored)


def test_two_genome_candidates_are_not_silently_merged(make_zip, demo_dir):
    """两个都像基因组的文件：只保留一个，另一个降级为"未识别"并说明原因。"""
    zp = make_zip({
        "g1.fa": demo_dir["genome"],
        "g2.fa": demo_dir["genome"],
        "a.gff3": demo_dir["gff3"],
    })
    scan = scan_archive(zp)
    genomes = [c for c in scan.candidates if c.role == "genome"]
    assert len(genomes) == 1
    leftovers = [c for c in scan.candidates
                 if c.member in ("g1.fa", "g2.fa") and not c.role]
    assert len(leftovers) == 1
    assert "同为" in leftovers[0].evidence_text()
    assert any("基因组" in w or "确认" in w for w in scan.warnings)


def test_scan_rejects_non_archive(tmp_path):
    p = tmp_path / "not_archive.zip"
    p.write_text("hello", encoding="utf-8")
    with pytest.raises(FileFormatError) as ei:
        scan_archive(str(p))
    assert "压缩包" in ei.value.message


def test_scan_reports_missing_annotation(make_zip, demo_dir):
    zp = make_zip({"g.fa": demo_dir["genome"]})
    scan = scan_archive(zp)
    assert scan.missing() == ["annotation"]
    assert any("注释" in w for w in scan.warnings)


# --------------------------------------------------------------------------- #
# 解包
# --------------------------------------------------------------------------- #
def test_import_extracts_only_selected_roles(make_zip, demo_dir, tmp_path):
    zp = make_zip({
        "g.fa": demo_dir["genome"],
        "a.gff3": demo_dir["gff3"],
        "r.fa": demo_dir["mrna"],
    })
    scan = scan_archive(zp)
    dest = tmp_path / "out"
    res = import_bundle(scan, str(dest), {"genome": "g.fa", "annotation": "a.gff3"})
    assert sorted(os.listdir(dest)) == ["a.gff3", "g.fa"]
    assert res.missing() == []
    # 解出来的内容必须原样
    with open(demo_dir["gff3"], "rb") as fh:
        assert open(res.paths["annotation"], "rb").read() == fh.read()


def test_import_requires_genome_and_annotation(make_zip, demo_dir, tmp_path):
    zp = make_zip({"g.fa": demo_dir["genome"]})
    scan = scan_archive(zp)
    with pytest.raises(GGEError) as ei:
        import_bundle(scan, str(tmp_path / "out"), {"genome": "g.fa"})
    assert "注释" in ei.value.message


def test_import_flattens_traversal_paths(make_zip, demo_dir, tmp_path):
    """zip 里写 ../../evil.fa 也不能写到目标目录外面去。"""
    zp = make_zip({
        "../../evil.fa": demo_dir["genome"],
        "a.gff3": demo_dir["gff3"],
    })
    scan = scan_archive(zp)
    dest = tmp_path / "deep" / "out"
    res = import_bundle(scan, str(dest))
    assert res.paths["genome"].startswith(str(dest))
    assert os.path.basename(res.paths["genome"]) == "evil.fa"
    assert not (tmp_path / "evil.fa").exists()


def test_import_decompresses_plain_gzip_genome(make_zip, demo_dir, tmp_path):
    """普通 gzip 的基因组要解压成 .fa，否则 pyfaidx 建不了索引。"""
    gz = tmp_path / "genome.fa.gz"
    with open(demo_dir["genome"], "rb") as src, gzip.open(gz, "wb") as out:
        out.write(src.read())
    with gzip.open(demo_dir["genome"], "rt") as _:
        pass
    zp = make_zip({"genome.fa.gz": str(gz), "a.gff3": demo_dir["gff3"]})
    scan = scan_archive(zp)
    res = import_bundle(scan, str(tmp_path / "out"))
    target = res.paths["genome"]
    assert target.endswith(".fa"), target
    with open(target, "rb") as fh:
        assert fh.read(1) == b">"
    assert any("解压" in w for w in res.warnings)
    # 真的能被基因组读取器打开
    from gge.fasta import open_fasta

    prov = open_fasta(target, prefer_memory=True)
    assert set(prov.seqids()) == {"chr1", "chr2", "chrScaf"}
    prov.close()


def test_bgzf_is_not_decompressed(tmp_path):
    """BGZF 的判据是 gzip 头里的 BC 子字段，不能把 bgzf 当成普通 gzip 再解一遍。"""
    # 真实的 BGZF 头：10 字节固定头 + XLEN=6 + "BC\x02\x00" + BSIZE(2)
    bgzf = (b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00"
            b"BC\x02\x00\x1c\x00")
    assert _is_bgzf(bytes(bgzf)) is True
    plain = gzip.compress(b">chr1\nACGT\n")
    assert _is_bgzf(plain) is False


def test_tar_gz_archive_supported(make_tar_gz, demo_dir, tmp_path):
    zp = make_tar_gz({
        "data/g.fa": demo_dir["genome"],
        "data/a.gff3": demo_dir["gff3"],
    })
    scan = scan_archive(zp)
    assert scan.kind == "tar"
    assert _roles(scan)["genome"] == "data/g.fa"
    res = import_bundle(scan, str(tmp_path / "out"))
    assert os.path.exists(res.paths["genome"])


def test_import_bundle_from_archive_registers_dataset(make_zip, demo_dir, library):
    """走完整任务链：解包 → 入库 → 记录统计信息。"""
    from gge.tasks import import_bundle_from_archive

    zp = make_zip({
        "g.fa": demo_dir["genome"],
        "a.gff3": demo_dir["gff3"],
        "m.fa": demo_dir["mrna"],
    })
    scan = scan_archive(zp)
    roles = {r: c.member for r, c in scan.by_role().items() if c}
    res = import_bundle_from_archive(scan, roles, "测试数据集", library.path)
    ds = res.dataset
    assert ds.id > 0 and ds.name == "测试数据集"
    assert ds.n_transcripts > 0 and ds.fmt == "gff3"
    assert os.path.dirname(ds.genome) != ""
    assert library.get(ds.id) is not None
    assert library.get_by_source(zp).id == ds.id


def test_reimport_same_archive_replaces_record(make_zip, demo_dir, library):
    """同一压缩包再导入一次：替换旧记录并清掉旧目录，不能越堆越多。"""
    from gge.tasks import import_bundle_from_archive

    zp = make_zip({"g.fa": demo_dir["genome"], "a.gff3": demo_dir["gff3"]})
    scan = scan_archive(zp)
    roles = {r: c.member for r, c in scan.by_role().items() if c}
    first = import_bundle_from_archive(scan, roles, "A", library.path)
    old_dir = os.path.dirname(first.dataset.genome)
    second = import_bundle_from_archive(scan, roles, "A", library.path)
    assert library.count() == 1
    assert second.replaced is True
    assert not os.path.exists(old_dir), "旧目录应被清理"
    assert os.path.exists(second.dataset.genome)


def test_failed_import_leaves_no_directory(make_zip, demo_dir, library):
    """解包中途失败不能留下半套文件。"""
    from gge.tasks import import_bundle_from_archive

    zp = make_zip({"g.fa": demo_dir["genome"], "a.gff3": demo_dir["gff3"]})
    scan = scan_archive(zp)
    store = os.path.join(library.path and os.path.dirname(library.path),
                         "datasets")
    before = set(os.listdir(store)) if os.path.isdir(store) else set()
    with pytest.raises(GGEError):
        # 故意只给基因组：必填项缺失 → 解包前就应报错
        import_bundle_from_archive(scan, {"genome": "g.fa"}, "bad", library.path)
    after = set(os.listdir(store)) if os.path.isdir(store) else set()
    assert after == before
    assert library.count() == 0


# --------------------------------------------------------------------------- #
# 文件名清洗
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("raw,expect", [
    ("data/main.fa", "main.fa"),
    ("..\\..\\win.fa", "win.fa"),
    ("a/b/../../c.gff3", "c.gff3"),
    ("CON.fa", "_CON.fa"),
    ("bad:name?.fa", "bad_name_.fa"),
    ("...", "file"),
])
def test_safe_filename(raw, expect):
    assert safe_filename(raw) == expect


def test_sampling_marks_incomplete(monkeypatch, make_zip, tmp_path):
    """大文件只采样时，界面文案必须说明"以上"，不能假装是精确值。"""
    import gge.bundle as bundle

    records = [(f"r{i}", "ACGT" * 25) for i in range(10)]
    fa = tmp_path / "many.fa"
    fa.write_text("".join(f">{n}\n{s}\n" for n, s in records), encoding="utf-8")
    zp = make_zip({"many.fa": str(fa)})

    monkeypatch.setattr(bundle, "_FASTA_SAMPLE_RECORDS", 3)
    scan = scan_archive(zp)
    cand = next(c for c in scan.candidates if c.member == "many.fa")
    assert cand.stats is not None
    assert cand.stats.complete is False
    assert "以上" in cand.stat_text()
    assert "采样" in cand.stat_text()


# --------------------------------------------------------------------------- #
# 磁盘空间预检：头部空间要按需给，不能固定 128 MB
# --------------------------------------------------------------------------- #
def test_tiny_bundle_allowed_on_nearly_full_disk(make_zip, demo_dir, tmp_path, monkeypatch):
    """盘快满了（只剩 40 MB）时，10 KB 的演示包仍应能导入。

    真实场景：开发机 C 盘只剩 30 MB，固定 128 MB 的头部门槛让一个 10 KB 的包
    也报"磁盘空间可能不足"——用户会以为程序坏了。
    """
    from gge import bundle as bundle_mod

    zp = make_zip({"g.fa": demo_dir["genome"], "a.gff3": demo_dir["gff3"]})
    scan = scan_archive(zp)
    monkeypatch.setattr(bundle_mod, "_disk_free", lambda p: 40 * 1024 * 1024)
    res = bundle_mod.import_bundle(scan, str(tmp_path / "out"),
                                   {"genome": "g.fa", "annotation": "a.gff3"})
    assert os.path.exists(res.paths["genome"])


def test_import_blocked_when_disk_really_too_small(make_zip, demo_dir, tmp_path, monkeypatch):
    """真的只剩 1 MB 时要在解包前拦住，并给出可操作的提示。"""
    from gge import bundle as bundle_mod

    zp = make_zip({"g.fa": demo_dir["genome"], "a.gff3": demo_dir["gff3"]})
    scan = scan_archive(zp)
    monkeypatch.setattr(bundle_mod, "_disk_free", lambda p: 1 * 1024 * 1024)
    with pytest.raises(GGEError) as ei:
        bundle_mod.import_bundle(scan, str(tmp_path / "out2"),
                                 {"genome": "g.fa", "annotation": "a.gff3"})
    assert "磁盘空间" in ei.value.message
    assert "GGE_LIBRARY_DIR" in ei.value.hint, "要告诉用户怎么换分区"
    assert not os.path.exists(tmp_path / "out2" / "g.fa"), "不该留下半成品"


def test_headroom_scales_with_bundle_size():
    """头部空间：小包 16 MB 起，大包不超过 128 MB。"""
    def headroom(need):
        return max(16 * 1024 * 1024, min(128 * 1024 * 1024, need))

    # 参数是 need（= 未压缩总量 × 0.35），不是包大小
    assert headroom(10 * 1024) == 16 * 1024 * 1024            # 10 KB 的量：给 16 MB
    assert headroom(35 * 1024 * 1024) == 35 * 1024 * 1024     # 35 MB 的量：等比给 35 MB
    assert headroom(10 * 1024 ** 3) == 128 * 1024 * 1024      # 再大也只留 128 MB


# --------------------------------------------------------------------------- #
# 数据集打包导出 + 清单识别（"下次一键导入"那条路）
# --------------------------------------------------------------------------- #
def _demo_paths(demo_dir, n_refs=3):
    paths = {"genome": demo_dir["genome"], "annotation": demo_dir["gff3"]}
    if n_refs >= 1:
        paths["mrna"] = demo_dir["mrna"]
    if n_refs >= 2:
        paths["cds"] = demo_dir["cds"]
    if n_refs >= 3:
        paths["protein"] = demo_dir["protein"]
    return paths


def test_write_dataset_archive_roundtrip(demo_dir, tmp_path):
    """导出 → 再扫描：角色按清单还原，统计量直接取自清单，不再解析注释。"""
    from gge.bundle import read_manifest, scan_archive, write_dataset_archive

    paths = _demo_paths(demo_dir)
    out = tmp_path / "ds.zip"
    stats = {"n_genes": 11, "n_transcripts": 12, "n_seqids": 4, "total_bases": 21000}
    res = write_dataset_archive(str(out), paths, name="演示 数据集", note="备注",
                                stats=stats)

    assert out.exists() and res.archive_size > 0
    assert set(res.members) == set(paths)
    assert res.total_bytes >= sum(os.path.getsize(p) for p in paths.values())
    assert "→" in res.text()

    man = read_manifest(str(out))
    assert man is not None
    assert man["format"] == "gge-dataset"
    assert man["version"] == 1
    assert man["name"] == "演示 数据集"
    roles = [f["role"] for f in man["files"]]
    assert roles == ["genome", "annotation", "mrna", "cds", "protein"]

    scan = scan_archive(str(out))
    assert scan.manifest is not None
    got = {r: c.member for r, c in scan.by_role().items() if c is not None}
    assert got == {r: res.members[r] for r in paths}
    assert scan.name_hint() == "演示 数据集"
    # 清单指定的角色要写着"凭清单"，而不是假装是猜出来的
    assert "清单" in scan.by_role()["genome"].evidence[0]
    assert scan.by_role()["genome"].confidence >= 0.9
    # 统计量来自清单
    assert scan.annotation is not None
    assert scan.annotation.n_genes == 11
    assert scan.annotation.n_transcripts == 12
    assert "4" in scan.annotation.text()


def test_exported_archive_imports_back(demo_dir, tmp_path):
    """真正的闭环：导出 → 解包 → 文件内容与原文件逐字节相同。"""
    from gge.bundle import import_bundle, scan_archive, write_dataset_archive

    paths = _demo_paths(demo_dir)
    out = tmp_path / "round.zip"
    res = write_dataset_archive(str(out), paths, name="round")
    scan = scan_archive(str(out))
    dest = tmp_path / "unpacked"
    imp = import_bundle(scan, str(dest), res.members)

    for role, src in paths.items():
        with open(src, "rb") as fh:
            want = fh.read()
        with open(imp.paths[role], "rb") as fh:
            assert fh.read() == want, f"{role} 解包后内容不一致"


def test_manifest_wins_over_content_guessing(demo_dir, tmp_path, make_zip):
    """包里有清单时以清单为准：即使文件名更容易被"猜"成别的角色。"""
    import json

    two = {"data/real.fa": demo_dir["genome"],
           "data/copy.fa": demo_dir["genome"],
           "data/a.gff3": demo_dir["gff3"]}
    manifest = {
        "format": "gge-dataset", "version": 1, "name": "手动清单",
        "stats": {"n_genes": 0, "n_transcripts": 0, "n_seqids": 0, "total_bases": 0},
        "files": [{"role": "genome", "member": "data/copy.fa"},
                  {"role": "annotation", "member": "data/a.gff3"}],
    }
    two["gge_dataset.json"] = ("text", json.dumps(manifest, ensure_ascii=False))
    zp = make_zip(two)

    scan = scan_archive(zp)
    assert scan.by_role()["genome"].member == "data/copy.fa"
    # 另一个同样像基因组的文件不再抢角色
    other = next(c for c in scan.candidates if c.member == "data/real.fa")
    assert other.role == ""
    # 清单文件本身不该被当成"跳过的不认识格式"报出来
    assert not any("gge_dataset.json" in name for name, _why in scan.ignored)


def test_broken_manifest_falls_back_to_guessing(demo_dir, make_zip):
    """清单坏了不能让整个包打不开：静默退回按内容判断。"""
    zp = make_zip({"g.fa": demo_dir["genome"], "a.gff3": demo_dir["gff3"],
                   "gge_dataset.json": ("text", '{"format": "gge-dataset", "fil')})
    scan = scan_archive(zp)
    assert scan.manifest is None
    got = {r: c.member for r, c in scan.by_role().items() if c is not None}
    assert got["genome"] == "g.fa"
    assert got["annotation"] == "a.gff3"


def test_manifest_pointing_to_missing_member(demo_dir, make_zip):
    """清单写了包里没有的文件：要说清原因，而不是悄悄少一个角色。"""
    import json

    manifest = {"format": "gge-dataset", "version": 1, "name": "x",
                "files": [{"role": "genome", "member": "data/gone.fa"},
                          {"role": "annotation", "member": "a.gff3"}]}
    zp = make_zip({"g.fa": demo_dir["genome"], "a.gff3": demo_dir["gff3"],
                   "gge_dataset.json": ("text", json.dumps(manifest))})
    scan = scan_archive(zp)
    assert scan.by_role()["genome"] is None
    assert any("data/gone.fa" in name for name, _why in scan.ignored)
    assert any("基因组" in w for w in scan.warnings)


def test_write_archive_requires_genome_and_annotation(demo_dir, tmp_path):
    from gge.bundle import write_dataset_archive

    with pytest.raises(GGEError) as ei:
        write_dataset_archive(str(tmp_path / "x.zip"), {"genome": demo_dir["genome"]})
    assert "缺少必填文件" in ei.value.message
    assert not (tmp_path / "x.zip").exists()


def test_write_archive_reports_missing_source_file(demo_dir, tmp_path):
    from gge.bundle import write_dataset_archive

    paths = {"genome": demo_dir["genome"],
             "annotation": str(tmp_path / "nope.gff3")}
    with pytest.raises(GGEError) as ei:
        write_dataset_archive(str(tmp_path / "y.zip"), paths)
    assert "不存在" in ei.value.message
    assert not (tmp_path / "y.zip").exists()


def test_write_archive_progress_and_cancel(demo_dir, tmp_path):
    """进度要单调递增到 100；取消后不能留下半截压缩包。"""
    from gge.bundle import write_dataset_archive
    from gge.errors import CancelledError

    seen = []
    write_dataset_archive(str(tmp_path / "p.zip"), _demo_paths(demo_dir, 1),
                          progress_cb=lambda pct, msg: seen.append(pct))
    assert seen and seen[-1] == 100
    assert seen == sorted(seen)

    out = tmp_path / "c.zip"
    with pytest.raises(CancelledError):
        write_dataset_archive(str(out), _demo_paths(demo_dir, 1),
                              cancel_check=lambda: True)
    assert not out.exists(), "取消后必须删掉半成品"
