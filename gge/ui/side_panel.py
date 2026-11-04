"""左侧栏：文件导入、基因搜索、基因信息与自检。"""

from __future__ import annotations

from typing import Dict, List, Optional

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor, QFont
from PySide6.QtWidgets import (QAbstractItemView, QComboBox, QFileDialog, QFrame,
                               QGridLayout, QGroupBox, QHBoxLayout, QHeaderView,
                               QLabel, QLineEdit, QListWidget, QListWidgetItem,
                               QPushButton, QSizePolicy, QTableWidget,
                               QTableWidgetItem, QToolButton, QVBoxLayout, QWidget)

from ..models import CheckItem, GeneExtraction
from ..search import SearchHit
from .theme import ui_font

FILE_SPECS = [
    ("genome", "基因组序列", "FASTA (*.fa *.fasta *.fna *.fas *.fa.gz *.fna.gz);;所有文件 (*)", True),
    ("annotation", "基因注释", "GFF/GTF (*.gff3 *.gff *.gtf *.gff3.gz *.gtf.gz);;所有文件 (*)", True),
    ("mrna", "参考 mRNA", "FASTA (*.fa *.fasta *.fna *.gz);;所有文件 (*)", False),
    ("cds", "参考 CDS", "FASTA (*.fa *.fasta *.fna *.gz);;所有文件 (*)", False),
    ("protein", "参考蛋白", "FASTA (*.fa *.fasta *.faa *.gz);;所有文件 (*)", False),
]

ROLE_LABELS = {
    "genome": "基因组", "annotation": "注释", "mrna": "mRNA",
    "cds": "CDS", "protein": "蛋白",
}


class FileRow(QWidget):
    browseRequested = Signal(str)
    clearRequested = Signal(str)

    def __init__(self, role: str, label: str, filters: str, required: bool = False,
                 parent=None) -> None:
        super().__init__(parent)
        self.role = role
        self.filters = filters
        lay = QGridLayout(self)
        lay.setContentsMargins(0, 2, 0, 2)
        lay.setHorizontalSpacing(4)
        lay.setVerticalSpacing(1)

        title = QLabel(label + (" *" if required else ""))
        title.setStyleSheet("color:#374151; font-weight:600; font-size:11px;")
        title.setToolTip("必填" if required else "可选")
        self.path_edit = QLineEdit()
        self.path_edit.setReadOnly(True)
        self.path_edit.setPlaceholderText("尚未选择文件")
        self.path_edit.setStyleSheet("font-size:11px; padding:2px 4px;")

        self.btn = QToolButton()
        self.btn.setText("选择")
        self.btn.setToolTip("选择文件")
        self.btn_clear = QToolButton()
        self.btn_clear.setText("×")
        self.btn_clear.setToolTip("清除该文件")
        self.btn_clear.setStyleSheet("QToolButton{color:#9ca3af;font-size:13px;}")

        self.status = QLabel("未加载")
        self.status.setStyleSheet("color:#9ca3af; font-size:11px;")
        self.status.setWordWrap(True)

        lay.addWidget(title, 0, 0, 1, 3)
        lay.addWidget(self.path_edit, 1, 0, 1, 2)
        lay.addWidget(self.btn, 1, 2)
        lay.addWidget(self.status, 2, 0, 1, 2)
        lay.addWidget(self.btn_clear, 2, 2)
        lay.setColumnStretch(0, 1)

        self.btn.clicked.connect(lambda: self.browseRequested.emit(self.role))
        self.btn_clear.clicked.connect(lambda: self.clearRequested.emit(self.role))

    def set_state(self, path: str, status: str = "", level: str = "ok") -> None:
        self.path_edit.setText(path)
        self.path_edit.setToolTip(path)
        colors = {"ok": "#15803d", "warn": "#b45309", "error": "#b91c1c", "idle": "#9ca3af"}
        prefix = {"ok": "✓ ", "warn": "⚠ ", "error": "✕ ", "idle": ""}[level]
        self.status.setText(prefix + status)
        self.status.setStyleSheet(f"color:{colors.get(level, '#6b7280')}; font-size:11px;")

    def clear(self) -> None:
        self.path_edit.clear()
        self.set_state("", "未加载", "idle")


class FilePanel(QGroupBox):
    """五个独立文件选择行 + 压缩包/数据集库入口。"""

    browseRequested = Signal(str)
    clearRequested = Signal(str)
    importRequested = Signal()
    libraryRequested = Signal()
    exportRequested = Signal()

    def __init__(self, parent=None) -> None:
        super().__init__("1 · 数据文件", parent)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(8, 6, 8, 8)
        lay.setSpacing(2)

        # 一行两个入口：整包导入（推荐） / 浏览已导入的数据集
        top = QHBoxLayout()
        top.setSpacing(4)
        self.btn_import = QToolButton()
        self.btn_import.setText("导入压缩包…")
        self.btn_import.setToolTip(
            "选择一个 zip/tar 包，自动识别里面的基因组、注释与参考序列\n"
            "支持 .zip / .tar.gz / .tgz 等（Ctrl+Shift+I）")
        self.btn_import.setStyleSheet(
            "QToolButton{border:1px solid #c7d2fe;background:#eef2ff;border-radius:4px;"
            "padding:3px 8px;color:#3730a3;font-size:11px;}"
            "QToolButton:hover{background:#e0e7ff;}")
        self.btn_library = QToolButton()
        self.btn_library.setText("数据集库")
        self.btn_library.setToolTip("查看已导入的数据集，一键切换（Ctrl+L）")
        self.btn_library.setStyleSheet(
            "QToolButton{border:1px solid #e5e7eb;border-radius:4px;padding:3px 8px;"
            "color:#374151;font-size:11px;}"
            "QToolButton:hover{background:#f3f4f6;}")
        self.btn_export_ds = QToolButton()
        self.btn_export_ds.setText("打包导出")
        self.btn_export_ds.setToolTip(
            "把当前这几个文件打包成一个 zip（内含 gge_dataset.json 清单）\n"
            "下次「导入压缩包」时按清单直接认回角色，不用再逐个选（Ctrl+Shift+E）")
        self.btn_export_ds.setStyleSheet(
            "QToolButton{border:1px solid #e5e7eb;border-radius:4px;padding:3px 8px;"
            "color:#374151;font-size:11px;}"
            "QToolButton:hover{background:#f3f4f6;}")
        top.addWidget(self.btn_import)
        top.addWidget(self.btn_library)
        top.addWidget(self.btn_export_ds)
        top.addStretch(1)
        lay.addLayout(top)

        self.dataset_label = QLabel("")
        self.dataset_label.setStyleSheet("color:#4b5563; font-size:10px;")
        self.dataset_label.setWordWrap(True)
        self.dataset_label.hide()
        lay.addWidget(self.dataset_label)

        self.rows: Dict[str, FileRow] = {}
        for role, label, filters, required in FILE_SPECS:
            row = FileRow(role, label, filters, required)
            row.browseRequested.connect(self.browseRequested)
            row.clearRequested.connect(self.clearRequested)
            lay.addWidget(row)
            self.rows[role] = row
        hint = QLabel("基因组与注释为必填；参考序列用于校验，可稍后再加载。")
        hint.setStyleSheet("color:#9ca3af; font-size:10px;")
        hint.setWordWrap(True)
        lay.addWidget(hint)

        self.btn_import.clicked.connect(self.importRequested)
        self.btn_library.clicked.connect(self.libraryRequested)
        self.btn_export_ds.clicked.connect(self.exportRequested)

    def set_dataset(self, name: str, detail: str = "") -> None:
        """显示"当前来自哪个数据集"，名称为空时隐藏。"""
        if not name:
            self.dataset_label.hide()
            self.dataset_label.setText("")
            return
        text = f"当前数据集：{name}"
        if detail:
            text += f"（{detail}）"
        self.dataset_label.setText(text)
        self.dataset_label.setToolTip(text)
        self.dataset_label.show()


class SearchPanel(QGroupBox):
    """基因搜索 + 结果列表。"""

    searchRequested = Signal(str)
    geneChosen = Signal(str)     # gene_id

    HEADERS = ["基因", "位置", "链", "转录本", "命中方式"]

    def __init__(self, parent=None) -> None:
        super().__init__("2 · 基因搜索", parent)
        self._hits: List[SearchHit] = []
        lay = QVBoxLayout(self)
        lay.setContentsMargins(8, 6, 8, 8)
        lay.setSpacing(4)

        top = QHBoxLayout()
        top.setSpacing(4)
        self.edit = QLineEdit()
        self.edit.setPlaceholderText("Gene Symbol / Gene ID / Transcript ID")
        self.edit.setClearButtonEnabled(True)
        self.btn = QPushButton("搜索")
        self.btn.setObjectName("primary")
        self.btn.setFixedWidth(58)
        top.addWidget(self.edit, 1)
        top.addWidget(self.btn)
        lay.addLayout(top)

        self.table = QTableWidget(0, len(self.HEADERS))
        self.table.setHorizontalHeaderLabels(self.HEADERS)
        self.table.verticalHeader().setVisible(False)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setAlternatingRowColors(True)
        self.table.setMinimumHeight(120)
        hh = self.table.horizontalHeader()
        hh.setSectionResizeMode(0, QHeaderView.Stretch)
        for i in range(1, len(self.HEADERS)):
            hh.setSectionResizeMode(i, QHeaderView.ResizeToContents)
        lay.addWidget(self.table, 1)

        self.info = QLabel("输入基因名或 ID 后回车")
        self.info.setStyleSheet("color:#9ca3af; font-size:10px;")
        lay.addWidget(self.info)

        self.btn.clicked.connect(self._emit_search)
        self.edit.returnPressed.connect(self._emit_search)
        self.table.cellClicked.connect(self._on_cell)
        self.table.itemActivated.connect(lambda _i: self._on_cell(self.table.currentRow(), 0))

    def _emit_search(self) -> None:
        self.searchRequested.emit(self.edit.text().strip())

    def set_hits(self, hits: List[SearchHit]) -> None:
        self._hits = hits
        self.table.setRowCount(len(hits))
        for r, h in enumerate(hits):
            name = h.gene_name if h.gene_name and h.gene_name != h.gene_id else h.gene_id
            cells = [name, f"{h.seqid}:{h.coords_text()}", h.strand or ".",
                     str(h.n_transcripts), h.matched_by]
            for c, text in enumerate(cells):
                item = QTableWidgetItem(text)
                item.setToolTip(f"{h.gene_id}  {h.seqid}:{h.start}-{h.end}({h.strand})\n"
                                f"命中：{h.matched_by} = {h.matched_value}")
                if c == 0:
                    f = item.font()
                    f.setBold(True)
                    item.setFont(f)
                if c == 1:
                    item.setFont(_mono())
                self.table.setItem(r, c, item)
        if hits:
            self.info.setText(f"{len(hits)} 条结果 · 双击选择")
            self.table.selectRow(0)
        else:
            self.info.setText("没有匹配结果——可尝试更短的片段或基因 ID")

    def current_gene(self) -> Optional[str]:
        r = self.table.currentRow()
        if 0 <= r < len(self._hits):
            return self._hits[r].gene_id
        return None

    def current_hit(self) -> Optional[SearchHit]:
        r = self.table.currentRow()
        if 0 <= r < len(self._hits):
            return self._hits[r]
        return None

    def _on_cell(self, row: int, _col: int) -> None:
        if 0 <= row < len(self._hits):
            self.table.selectRow(row)
            self.geneChosen.emit(self._hits[row].gene_id)


def _mono():
    f = QFont("Consolas")
    f.setStyleHint(QFont.Monospace)
    f.setPointSize(8)
    return f


class InfoPanel(QGroupBox):
    """基因基本信息、转录本选择、自检结果。"""

    transcriptChanged = Signal(str)

    def __init__(self, parent=None) -> None:
        super().__init__("3 · 基因与转录本", parent)
        self._updating = False
        lay = QVBoxLayout(self)
        lay.setContentsMargins(8, 6, 8, 8)
        lay.setSpacing(4)

        self.title = QLabel("未选择基因")
        f = ui_font(10)
        f.setBold(True)
        self.title.setFont(f)
        self.title.setWordWrap(True)
        lay.addWidget(self.title)

        self.meta = QLabel("")
        self.meta.setStyleSheet("color:#4b5563; font-size:11px;")
        self.meta.setWordWrap(True)
        lay.addWidget(self.meta)

        row = QHBoxLayout()
        row.setSpacing(4)
        lab = QLabel("转录本")
        lab.setStyleSheet("color:#374151; font-size:11px;")
        self.combo = QComboBox()
        self.combo.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        self.combo.setStyleSheet("font-size:11px; padding:1px 4px;")
        row.addWidget(lab)
        row.addWidget(self.combo, 1)
        lay.addLayout(row)

        self.tx_note = QLabel("")
        self.tx_note.setStyleSheet("color:#b45309; font-size:10px;")
        self.tx_note.setWordWrap(True)
        lay.addWidget(self.tx_note)

        self.checks = QListWidget()
        self.checks.setMinimumHeight(110)
        self.checks.setStyleSheet("font-size:11px; border:1px solid #eef2f7;")
        self.checks.setWordWrap(True)
        lay.addWidget(self.checks, 1)

        self.warn_label = QLabel("")
        self.warn_label.setStyleSheet("color:#b45309; font-size:10px;")
        self.warn_label.setWordWrap(True)
        lay.addWidget(self.warn_label)

        self.combo.currentIndexChanged.connect(self._on_combo)

    def set_gene(self, gene, canonical_tx_id: str = "") -> None:
        self._updating = True
        try:
            self.title.setText(f"{gene.gene_name}  ({gene.gene_id})")
            txs = gene.tx_list()
            self.meta.setText(
                f"{gene.seqid}:{gene.start:,}-{gene.end:,} ({gene.strand}) · "
                f"{gene.length:,} bp · {gene.biotype or '类型未知'}\n"
                f"转录本 {len(txs)} 条"
            )
            self.combo.clear()
            ordered = sorted(txs, key=lambda t: (-t.cds_length(), -t.mrna_length(), t.tx_id))
            for tx in ordered:
                label = f"{tx.tx_id}  ·  {tx.n_exons} 外显子 · CDS {tx.cds_length():,} bp"
                self.combo.addItem(label, tx.tx_id)
            if canonical_tx_id:
                idx = self.combo.findData(canonical_tx_id)
                if idx >= 0:
                    self.combo.setCurrentIndex(idx)
                self.tx_note.setText(
                    f"已自动选择 CDS 最长的转录本 {canonical_tx_id}；"
                    f"如需查看其它转录本请在上方切换（不会混合不同转录本的结构）"
                )
        finally:
            self._updating = False

    def _on_combo(self, _idx: int) -> None:
        if self._updating:
            return
        tx_id = self.combo.currentData()
        if tx_id:
            self.transcriptChanged.emit(tx_id)

    def current_tx(self) -> str:
        return self.combo.currentData() or ""

    def set_checks(self, checks: List[CheckItem], warnings: List[str]) -> None:
        self.checks.clear()
        order = {"error": 0, "warn": 1, "info": 2, "ok": 3}
        colors = {"error": "#b91c1c", "warn": "#b45309", "info": "#2563eb", "ok": "#15803d"}
        icons = {"error": "✕", "warn": "⚠", "info": "i", "ok": "✓"}
        for c in sorted(checks, key=lambda x: order.get(x.level, 9)):
            item = QListWidgetItem(f"{icons.get(c.level, '·')} {c.title}")
            item.setForeground(QColor(colors.get(c.level, "#374151")))
            item.setToolTip(c.detail)
            if c.detail:
                item.setText(f"{icons.get(c.level, '·')} {c.title} — {c.detail}")
            self.checks.addItem(item)
        if warnings:
            self.warn_label.setText(f"⚠ {len(warnings)} 条警告：" + warnings[0]
                                    + ("…" if len(warnings) > 1 else ""))
            self.warn_label.setToolTip("\n".join(warnings))
        else:
            self.warn_label.setText("")
            self.warn_label.setToolTip("")

    def clear(self) -> None:
        self.title.setText("未选择基因")
        self.meta.setText("")
        self.combo.clear()
        self.tx_note.setText("")
        self.checks.clear()
        self.warn_label.setText("")
