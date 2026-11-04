"""无 GUI 依赖的"任务函数"。

GUI 只负责在后台线程里调用它们并显示进度，因此这些函数也可以直接在
脚本或测试中使用（``tests`` 里就有用例）。所有长任务都支持
``progress_cb(percent, message)`` 与 ``cancel_check() -> bool``。
"""

from __future__ import annotations

import os
import traceback
from dataclasses import dataclass, field
from typing import Callable, Dict, Iterable, List, Optional, Sequence

from .annotation import AnnotationModel, parse_annotation
from .coords import human_bytes, human_size
from .errors import CancelledError, GGEError
from .export import (gene_basic_info_text, gene_fasta_text, safe_filename,
                     segments_table, validation_report_text, write_csv, write_text)
from .extract import PANEL_ORDER, extract_gene
from .fasta import open_fasta
from .models import GeneExtraction
from .validate import RefSet, load_reference, validate_extraction

ProgressCB = Optional[Callable[[int, str], None]]
CancelCheck = Optional[Callable[[], bool]]


def _check(cancel_check: CancelCheck) -> None:
    if cancel_check and cancel_check():
        raise CancelledError()


@dataclass
class GenomeInfo:
    path: str
    provider: object
    seqids: List[str]
    lengths: Dict[str, int]
    mode: str = ""

    @property
    def total_bases(self) -> int:
        return sum(self.lengths.values())

    def summary(self) -> str:
        return (f"{len(self.seqids):,} 条序列，共 {human_size(self.total_bases)}"
                f"（{self.mode}）")


def load_genome(path: str, progress_cb: ProgressCB = None,
                cancel_check: CancelCheck = None) -> GenomeInfo:
    """打开 FASTA 并建立索引。"""
    if progress_cb:
        progress_cb(10, "打开 FASTA …")
    provider = open_fasta(path, progress_cb=progress_cb)
    if progress_cb:
        progress_cb(90, "统计序列信息 …")
    lengths = provider.lengths()
    if not lengths:
        raise GGEError("FASTA 中没有任何序列记录。", "确认文件没有损坏。")
    if progress_cb:
        progress_cb(100, "索引就绪")
    return GenomeInfo(path=path, provider=provider, seqids=provider.seqids(),
                      lengths=lengths, mode=getattr(provider, "mode", ""))


def load_annotation_file(path: str, progress_cb: ProgressCB = None,
                         cancel_check: CancelCheck = None) -> AnnotationModel:
    """解析注释。

    ``cancel_check`` 不是摆设：500 MB 的 RefSeq GFF 要读 70 s，没有它的话
    界面上的「取消」要等线程自己跑完才生效（旧版就是这样）。
    """
    return parse_annotation(path, progress_cb=progress_cb, cancel_check=cancel_check)


def load_reference_file(path: str, kind: str, progress_cb: ProgressCB = None,
                        cancel_check: CancelCheck = None) -> RefSet:
    if progress_cb:
        progress_cb(30, f"读取参考 {kind} …")
    rs = load_reference(path, kind, progress_cb=progress_cb, cancel_check=cancel_check)
    if progress_cb:
        progress_cb(100, f"参考 {kind}：{rs.summary()}")
    return rs


@dataclass
class ExtractionBundle:
    """一次"提取 + 校验"的完整产物。"""

    extraction: GeneExtraction
    validation: Dict[str, object] = field(default_factory=dict)
    gene: object = None
    tx: object = None


def extract_and_validate(provider,
                         model: AnnotationModel,
                         gene_id: str,
                         tx_id: str = "",
                         upstream: int = 2000,
                         downstream: int = 2000,
                         mrna_refs: Optional[RefSet] = None,
                         cds_refs: Optional[RefSet] = None,
                         protein_refs: Optional[RefSet] = None,
                         manual_mrna_id: str = "",
                         manual_cds_id: str = "",
                         manual_protein_id: str = "",
                         progress_cb: ProgressCB = None,
                         cancel_check: CancelCheck = None) -> ExtractionBundle:
    gene = model.get_gene(gene_id)
    if gene is None:
        raise GGEError(f"注释中找不到基因 {gene_id}", "重新搜索并选择一个基因。")
    tx = gene.transcripts.get(tx_id) if tx_id else gene.canonical_transcript()
    if tx is None:
        raise GGEError(
            f"基因 {gene_id} 中找不到转录本 {tx_id}",
            "该转录本可能只有注释行、没有 exon/CDS；请换一个转录本。",
        )
    if progress_cb:
        progress_cb(10, f"提取 {gene_id} / {tx.tx_id} …")
    _check(cancel_check)
    extraction = extract_gene(provider, gene, tx, upstream, downstream)
    extraction.annotation_path = model.path
    if progress_cb:
        progress_cb(55, "提取完成，准备校验 …")
    _check(cancel_check)

    validation = {}
    if mrna_refs or cds_refs or protein_refs:
        if progress_cb:
            progress_cb(65, "比对参考序列 …")
        validation = validate_extraction(
            extraction, gene, tx, mrna_refs, cds_refs, protein_refs,
            manual_mrna_id=manual_mrna_id,
            manual_cds_id=manual_cds_id,
            manual_protein_id=manual_protein_id,
        )
    if progress_cb:
        progress_cb(100, "完成")
    return ExtractionBundle(extraction=extraction, validation=validation, gene=gene, tx=tx)


def export_gene(extraction: GeneExtraction, outdir: str,
                validation: Optional[Dict[str, object]] = None,
                write_table: bool = True,
                write_info: bool = True) -> List[str]:
    """导出一个基因的全部文件，返回生成的文件路径列表。"""
    os.makedirs(outdir, exist_ok=True)
    base = safe_filename(f"{extraction.gene_id}_{extraction.tx_id}")
    written: List[str] = []
    fa = os.path.join(outdir, f"{base}.fasta")
    write_text(fa, gene_fasta_text(extraction))
    written.append(fa)
    if write_table:
        rows = segments_table(extraction)
        if rows:
            csv_path = os.path.join(outdir, f"{base}.segments.csv")
            write_csv(csv_path, rows)
            written.append(csv_path)
    if write_info:
        info_path = os.path.join(outdir, f"{base}.info.txt")
        write_text(info_path, gene_basic_info_text(extraction))
        written.append(info_path)
    if validation:
        rep = os.path.join(outdir, f"{base}.validation.txt")
        write_text(rep, validation_report_text(extraction, validation))
        written.append(rep)
    return written


@dataclass
class BatchResult:
    total: int = 0
    ok: int = 0
    failed: int = 0
    written: int = 0
    problems: List[tuple] = field(default_factory=list)   # (gene_id, message)
    rows: List[dict] = field(default_factory=list)

    def summary(self) -> str:
        return (f"成功 {self.ok} / {self.total} 个基因，生成 {self.written} 个文件，"
                f"失败 {self.failed} 个")


def batch_export(provider,
                 model: AnnotationModel,
                 gene_ids: Sequence[str],
                 outdir: str,
                 upstream: int = 2000,
                 downstream: int = 2000,
                 all_transcripts: bool = False,
                 write_table: bool = True,
                 mrna_refs: Optional[RefSet] = None,
                 cds_refs: Optional[RefSet] = None,
                 progress_cb: ProgressCB = None,
                 cancel_check: CancelCheck = None) -> BatchResult:
    """批量导出多个基因。单个基因失败不会中断整体，但会被记录下来。"""
    result = BatchResult(total=len(gene_ids))
    os.makedirs(outdir, exist_ok=True)
    for i, gid in enumerate(gene_ids):
        _check(cancel_check)
        if progress_cb:
            progress_cb(int(i / max(len(gene_ids), 1) * 100), f"导出 {gid}（{i + 1}/{len(gene_ids)}）")
        gene = model.get_gene(gid)
        if gene is None:
            result.failed += 1
            result.problems.append((gid, "注释中找不到该基因"))
            continue
        txs = gene.tx_list() if all_transcripts else [gene.canonical_transcript()]
        gene_ok = False
        for tx in txs:
            if tx is None:
                continue
            try:
                extraction = extract_gene(provider, gene, tx, upstream, downstream)
                extraction.annotation_path = model.path
                validation = {}
                if mrna_refs or cds_refs:
                    validation = validate_extraction(extraction, gene, tx,
                                                     mrna_refs, cds_refs)
                files = export_gene(extraction, outdir, validation, write_table=write_table)
                result.written += len(files)
                gene_ok = True
                result.rows.append({
                    "gene_id": extraction.gene_id,
                    "gene_name": extraction.gene_name,
                    "transcript": extraction.tx_id,
                    "seqid": extraction.seqid,
                    "start": extraction.gene_start,
                    "end": extraction.gene_end,
                    "strand": extraction.strand,
                    "exons": extraction.n_exons,
                    "cds_len": len(extraction.cds),
                    "mrna_len": len(extraction.mrna),
                    "mRNA_status": validation.get("mRNA").status_cn if validation.get("mRNA") else "",
                    "CDS_status": validation.get("CDS").status_cn if validation.get("CDS") else "",
                })
            except (GGEError, ValueError) as exc:
                result.failed += 1
                msg = exc.full_text() if isinstance(exc, GGEError) else str(exc)
                result.problems.append((f"{gid}/{tx.tx_id}", msg))
            except Exception as exc:  # pragma: no cover - 兜底
                result.failed += 1
                result.problems.append((f"{gid}/{tx.tx_id}", f"未预期错误：{exc}"))
        if gene_ok:
            result.ok += 1
    if result.rows:
        write_csv(os.path.join(outdir, "_batch_summary.csv"), result.rows)
    if result.problems:
        lines = ["gene\ttranscript\tproblem"]
        for gid, msg in result.problems:
            lines.append(f"{gid}\t-\t{msg.replace(chr(9), ' ').replace(chr(10), ' ')}")
        write_text(os.path.join(outdir, "_batch_problems.tsv"), "\n".join(lines) + "\n")
    if progress_cb:
        progress_cb(100, result.summary())
    return result


# --------------------------------------------------------------------------- #
# 压缩包导入（数据集库）
# --------------------------------------------------------------------------- #
@dataclass
class BundleImportResult:
    """一次性把"解包 + 登记进库"的结果带回来。"""

    dataset: object = None            # gge.library.Dataset
    scan: object = None               # gge.bundle.BundleScan
    import_result: object = None      # gge.bundle.BundleImport
    replaced: bool = False            # 是否覆盖了同一来源的旧记录

    def summary(self) -> str:
        ds = self.dataset
        name = getattr(ds, "name", "")
        return f"已导入数据集「{name}」：{self.import_result.text() if self.import_result else ''}"


def import_bundle_from_archive(scan,
                               roles: Dict[str, str],
                               name: str,
                               library_path: Optional[str] = None,
                               progress_cb: ProgressCB = None,
                               cancel_check: CancelCheck = None) -> BundleImportResult:
    """把扫描结果解包到库目录并登记。

    失败时会把刚建好的目录清掉，避免留下半套文件；同一个压缩包重复导入时，
    会替换旧记录（旧目录一并删除），不会越堆越多。
    """
    from .bundle import import_bundle as _import_bundle
    from .library import Library

    lib = Library(library_path)
    _check(cancel_check)
    dest = lib.new_dataset_dir(name or os.path.basename(scan.path))
    try:
        if progress_cb:
            progress_cb(5, "准备解包 …")
        result = _import_bundle(scan, dest, roles, progress_cb=progress_cb,
                                cancel_check=cancel_check)
        _check(cancel_check)

        ann = scan.annotation
        genome_stats = None
        for cand in scan.candidates:
            if cand.member == result.members.get("genome"):
                genome_stats = cand.stats
        stats = {
            "n_genes": getattr(ann, "n_genes", 0) if ann else 0,
            "n_transcripts": getattr(ann, "n_transcripts", 0) if ann else 0,
            "n_seqids": (len(getattr(ann, "seqids", []) or [])
                         or getattr(ann, "n_seqids_hint", 0)) if ann else 0,
            "total_bases": getattr(genome_stats, "total_bases", 0) if genome_stats else 0,
        }
        # 清单里记的总碱基数比"采样估算"更准，有就用它
        mstats = {}
        if getattr(scan, "manifest", None):
            from .bundle import _manifest_stats
            mstats = _manifest_stats(scan.manifest)
        if mstats.get("total_bases") and not stats["total_bases"]:
            stats["total_bases"] = mstats["total_bases"]
        old = lib.get_by_source(scan.path)
        if progress_cb:
            progress_cb(95, "写入数据集库 …")
        if old is not None:
            lib.remove(old.id, delete_files=True)
        dataset = lib.add(name=name or os.path.basename(scan.path), paths=result.paths,
                          source=scan.path, fmt=(ann.fmt if ann else ""),
                          stats=stats, members=result.members,
                          note="；".join(result.warnings))
        if progress_cb:
            progress_cb(100, "导入完成")
        return BundleImportResult(dataset=dataset, scan=scan, import_result=result,
                                  replaced=old is not None)
    except BaseException:
        # 中途失败：清掉半成品目录（它一定位于库目录内，不会碰到用户的原始文件）
        import shutil
        shutil.rmtree(dest, ignore_errors=True)
        raise


# --------------------------------------------------------------------------- #
# 数据集打包（把"分别导入的几个文件"固化成一个可再次导入的包）
# --------------------------------------------------------------------------- #
@dataclass
class DatasetExportResult:
    path: str = ""
    name: str = ""
    members: Dict[str, str] = field(default_factory=dict)
    total_bytes: int = 0
    archive_size: int = 0

    def summary(self) -> str:
        return (f"{len(self.members)} 个文件 → {os.path.basename(self.path)}"
                f"（{human_bytes(self.total_bytes)} → {human_bytes(self.archive_size)}）")


def export_dataset_archive(path: str, paths: Dict[str, str], name: str = "",
                           note: str = "", source: str = "",
                           stats: Optional[Dict[str, int]] = None,
                           progress_cb: ProgressCB = None,
                           cancel_check: CancelCheck = None) -> DatasetExportResult:
    """把当前这套文件打成一个带清单的 zip（后台任务入口）。"""
    from .bundle import write_dataset_archive

    if progress_cb:
        progress_cb(0, "准备打包 …")
    res = write_dataset_archive(path, paths, name=name, note=note, source=source,
                               stats=stats, progress_cb=progress_cb,
                               cancel_check=cancel_check)
    return DatasetExportResult(path=res.path, name=res.name, members=res.members,
                               total_bytes=res.total_bytes,
                               archive_size=res.archive_size)


# --------------------------------------------------------------------------- #
# 彩色序列文档（Word / 网页）
# --------------------------------------------------------------------------- #
@dataclass
class DocumentExportResult:
    files: List[str] = field(default_factory=list)
    bases: int = 0
    panels: int = 0

    def summary(self) -> str:
        names = "、".join(os.path.basename(f) for f in self.files)
        return f"已导出 {names}（{self.panels} 个面板 / {self.bases:,} bp）"


def export_sequence_document(path: str, extraction: GeneExtraction,
                             panel_keys: Optional[Sequence[str]] = None,
                             mode: str = "region", bpl: int = 60,
                             font_pt: float = 8.5,
                             palette: Optional[Dict[str, Dict[str, str]]] = None,
                             include_header: bool = True,
                             include_coord: bool = True,
                             progress_cb: ProgressCB = None,
                             cancel_check: CancelCheck = None) -> DocumentExportResult:
    """导出带颜色的序列文档；失败不修改任何原始文件。"""
    from .seqdoc import write_document

    keys = list(panel_keys or [k for k, p in extraction.panels.items() if not p.is_empty()])
    bases = sum(extraction.panels[k].length for k in keys if k in extraction.panels)
    panels = len([k for k in keys if k in extraction.panels])
    if progress_cb:
        progress_cb(5, f"排版 {panels} 个面板 …")
    _check(cancel_check)
    if progress_cb:
        progress_cb(40, "写入文档 …")
    files = write_document(path, extraction, keys, mode=mode, bpl=bpl,
                          font_pt=font_pt, palette=palette,
                          include_header=include_header,
                          include_coord=include_coord)
    _check(cancel_check)
    if progress_cb:
        progress_cb(100, f"已导出 {os.path.basename(files[0])}")
    return DocumentExportResult(files=files, bases=bases, panels=panels)
