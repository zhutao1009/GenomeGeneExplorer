"""对话框：差异详情、配色、批量导出、关于。"""

from __future__ import annotations

import html
import os
from typing import Dict, List, Optional, Sequence

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor, QFont
from PySide6.QtWidgets import (QAbstractItemView, QCheckBox, QColorDialog, QComboBox,
                               QDialog, QDialogButtonBox, QDoubleSpinBox, QFileDialog,
                               QFormLayout, QGridLayout, QGroupBox, QHBoxLayout,
                               QHeaderView, QLabel, QLineEdit, QMessageBox,
                               QPlainTextEdit, QPushButton, QSpinBox, QTableWidget,
                               QTableWidgetItem, QTabWidget, QTextEdit, QVBoxLayout,
                               QWidget)

from ..colors import DEFAULT_PALETTE, LEGEND_ORDER, ensure_palette, region_label
from ..coords import human_bytes
from ..extract import PANEL_ORDER, PANEL_TITLES
from ..models import GeneExtraction
from ..seqdoc import (DEFAULT_BPL, MODE_LABELS, MODE_REGION, WARN_BASES, estimate)
from ..validate import CheckResult, STATUS_CN
from .theme import mono_font, ui_font


def _badge(text: str, color: str) -> QLabel:
    lab = QLabel(text)
    lab.setStyleSheet(
        f"background:{color}; color:white; border-radius:8px;"
        f"padding:2px 10px; font-size:11px; font-weight:600;"
    )
    lab.setAlignment(Qt.AlignCenter)
    return lab


class DiffDialog(QDialog):
    """mRNA / CDS 校验详情：指标、差异列表、比对视图。"""

    def __init__(self, result: CheckResult, extraction: Optional[GeneExtraction] = None,
                 parent=None) -> None:
        super().__init__(parent)
        self.result = result
        self.extraction = extraction
        self.setWindowTitle(f"{result.kind} 校验详情")
        self.resize(880, 620)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(12, 12, 12, 10)
        lay.setSpacing(8)

        head = QHBoxLayout()
        head.addWidget(_badge(result.status_cn, result.color))
        title = QLabel(result.headline())
        f = ui_font(10)
        f.setBold(True)
        title.setFont(f)
        title.setWordWrap(True)
        head.addWidget(title, 1)
        lay.addLayout(head)

        tabs = QTabWidget()
        tabs.addTab(self._metrics_tab(), "指标")
        tabs.addTab(self._events_tab(), "差异位点")
        tabs.addTab(self._alignment_tab(), "比对视图")
        if result.candidates or result.notes:
            tabs.addTab(self._notes_tab(), "候选与说明")
        lay.addWidget(tabs, 1)

        buttons = QDialogButtonBox(QDialogButtonBox.Close)
        self.btn_export = QPushButton("导出校验报告…")
        buttons.addButton(self.btn_export, QDialogButtonBox.ActionRole)
        buttons.rejected.connect(self.reject)
        buttons.accepted.connect(self.accept)
        self.btn_export.clicked.connect(self._export)
        lay.addWidget(buttons)

    # ------------------------------------------------------------------ #
    def _metrics_tab(self) -> QWidget:
        w = QWidget()
        form = QFormLayout(w)
        form.setContentsMargins(10, 10, 10, 10)
        form.setVerticalSpacing(5)
        for key, val in self.result.metrics():
            k = QLabel(key)
            k.setStyleSheet("color:#4b5563;")
            v = QLabel(val)
            v.setTextInteractionFlags(Qt.TextSelectableByMouse)
            if key in ("参考序列",):
                v.setFont(mono_font(9))
            form.addRow(k, v)
        if self.result.comparison and self.result.comparison.coarse:
            warn = QLabel("⚠ 序列过长，本次仅做前后缀锚定的粗略定位，未做全序列动态规划比对")
            warn.setStyleSheet("color:#b45309; font-size:11px;")
            warn.setWordWrap(True)
            form.addRow(warn)
        return w

    def _events_tab(self) -> QWidget:
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setContentsMargins(8, 8, 8, 8)
        c = self.result.comparison
        events = list(c.events) if c else []
        if not events:
            lab = QLabel("没有差异位点" if self.result.status == "identical"
                         else "没有可定位的差异（可能是长度差异或需要人工检查）")
            lab.setStyleSheet("color:#6b7280;")
            lay.addWidget(lab)
            return w
        table = QTableWidget(len(events[:2000]), 4)
        table.setHorizontalHeaderLabels(["类型", "参考位置", "参考碱基", "观测碱基"])
        table.verticalHeader().setVisible(False)
        table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        kind_cn = {"mismatch": "替换", "insertion": "插入", "deletion": "缺失"}
        for r, ev in enumerate(events[:2000]):
            vals = [kind_cn.get(ev.kind, ev.kind), f"{ev.ref_pos:,}", ev.ref_base, ev.qry_base]
            for cidx, text in enumerate(vals):
                item = QTableWidgetItem(text)
                if cidx in (1, 2, 3):
                    item.setFont(mono_font(8))
                if ev.kind != "mismatch":
                    item.setBackground(QColor("#fef3c7"))
                table.setItem(r, cidx, item)
        lay.addWidget(table)
        if len(events) > 2000:
            lab = QLabel(f"仅显示前 2000 条，共 {len(events):,} 条")
            lab.setStyleSheet("color:#9ca3af; font-size:10px;")
            lay.addWidget(lab)
        return w

    def _alignment_tab(self) -> QWidget:
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setContentsMargins(8, 8, 8, 8)
        c = self.result.comparison
        if not c or not c.aligned_ref:
            lab = QLabel("没有可显示的比对结果")
            lab.setStyleSheet("color:#6b7280;")
            lay.addWidget(lab)
            return w
        view = QTextEdit()
        view.setReadOnly(True)
        view.setLineWrapMode(QTextEdit.NoWrap)
        view.setFont(mono_font(9))
        view.setHtml(self._alignment_html(c))
        lay.addWidget(view)
        return w

    def _alignment_html(self, c) -> str:
        ref = c.aligned_ref
        obs = c.aligned_qry
        n = len(ref)
        window = 160
        center = max(0, c.first_diff_ref)
        start = max(0, center - window // 2)
        end = min(n, start + window)
        start = max(0, end - window)
        piece_r = ref[start:end]
        piece_o = obs[start:end]

        def colorize(a: str, b: str, for_ref: bool) -> str:
            out = []
            for i in range(len(a)):
                ca = a[i]
                cb = b[i] if i < len(b) else "-"
                if ca == "-" or cb == "-":
                    bg = "#fde68a"
                elif ca != cb:
                    bg = "#fca5a5"
                else:
                    bg = "#ffffff"
                ch = ca if for_ref else cb
                out.append(f'<span style="background:{bg};">{html.escape(ch)}</span>')
            return "".join(out)

        head = (
            f"<div style='font-family:Consolas,monospace;font-size:11px;'>"
            f"列 {start + 1:,}–{end:,} / 共 {n:,} 列 &nbsp; "
            f"一致 {c.matches:,} · 替换 {c.mismatches:,} · "
            f"插入 {c.insertions:,} · 缺失 {c.deletions:,} · 一致率 {c.identity:.3f}%"
            f"</div><br>"
        )
        body = (
            f"<div style='font-family:Consolas,monospace;font-size:12px;line-height:1.35;'>"
            f"<b>ref</b> {colorize(piece_r, piece_o, True)}<br>"
            f"<b>obs</b> {colorize(piece_r, piece_o, False)}"
            f"</div>"
        )
        note = ""
        if start > 0 or end < n:
            note = ("<br><span style='color:#6b7280;font-size:11px;'>"
                    "为使差异可见，仅显示差异附近的一小段；完整比对可通过导出报告查看。</span>")
        color_legend = ("<br><span style='font-size:11px;color:#4b5563;'>"
                        "<span style='background:#fca5a5;'>&nbsp;&nbsp;</span> 替换 &nbsp; "
                        "<span style='background:#fde68a;'>&nbsp;&nbsp;</span> 插入/缺失</span>")
        return head + body + color_legend + note

    def _notes_tab(self) -> QWidget:
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setContentsMargins(8, 8, 8, 8)
        if self.result.candidates:
            lab = QLabel("候选参考序列（按得分降序，不自动取第一条）")
            lab.setStyleSheet("color:#374151; font-weight:600;")
            lay.addWidget(lab)
            table = QTableWidget(len(self.result.candidates), 3)
            table.setHorizontalHeaderLabels(["参考 ID", "得分", "命中依据"])
            table.verticalHeader().setVisible(False)
            table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
            table.horizontalHeader().setSectionResizeMode(2, QHeaderView.Stretch)
            table.setEditTriggers(QAbstractItemView.NoEditTriggers)
            for r, cand in enumerate(self.result.candidates[:50]):
                table.setItem(r, 0, QTableWidgetItem(cand.rid))
                table.setItem(r, 1, QTableWidgetItem(f"{cand.score:.0f}"))
                table.setItem(r, 2, QTableWidgetItem(cand.reason_text()))
            table.setMaximumHeight(200)
            lay.addWidget(table)
        txt = QPlainTextEdit()
        txt.setReadOnly(True)
        txt.setPlainText("\n".join(f"• {n}" for n in self.result.notes))
        lay.addWidget(txt, 1)
        return w

    def _export(self) -> None:
        from ..export import validation_report_text, write_text

        if self.extraction is None:
            QMessageBox.information(self, "无法导出", "缺少提取结果上下文。")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "保存校验报告", f"{self.extraction.gene_id}_validation.txt",
            "文本 (*.txt);;Markdown (*.md)")
        if not path:
            return
        write_text(path, validation_report_text(self.extraction, {self.result.kind: self.result}))
        QMessageBox.information(self, "已导出", f"报告已保存到：\n{path}")


class ColorSettingsDialog(QDialog):
    """配色自定义。"""

    def __init__(self, palette: Dict[str, Dict[str, str]], parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("配色设置")
        self.resize(460, 420)
        self._palette = ensure_palette(palette)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(12, 12, 12, 10)

        hint = QLabel("点击色块修改背景色；文字色会自动取同色系深色以保证可读。")
        hint.setStyleSheet("color:#6b7280; font-size:11px;")
        hint.setWordWrap(True)
        lay.addWidget(hint)

        self.table = QTableWidget(len(LEGEND_ORDER), 3)
        self.table.setHorizontalHeaderLabels(["区域", "背景色", "预览"])
        self.table.verticalHeader().setVisible(False)
        self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        self.table.horizontalHeader().setSectionResizeMode(2, QHeaderView.Stretch)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        for r, kind in enumerate(LEGEND_ORDER):
            conf = self._palette.get(kind, DEFAULT_PALETTE.get(kind, {}))
            self.table.setItem(r, 0, QTableWidgetItem(region_label(kind, self._palette)))
            btn = QPushButton(conf.get("bg", "#ffffff"))
            btn.setStyleSheet(f"background:{conf.get('bg', '#ffffff')}; border:1px solid #d1d5db;")
            btn.clicked.connect(lambda _c=False, k=kind, b=btn: self._pick(k, b))
            self.table.setCellWidget(r, 1, btn)
            prev = QLabel("ACGTACGT")
            prev.setAlignment(Qt.AlignCenter)
            prev.setFont(mono_font(9))
            prev.setStyleSheet(f"background:{conf.get('bg', '#ffffff')}; color:{conf.get('fg', '#111827')};")
            self.table.setCellWidget(r, 2, prev)
        lay.addWidget(self.table, 1)

        row = QHBoxLayout()
        btn_reset = QPushButton("恢复默认")
        btn_reset.clicked.connect(self._reset)
        row.addWidget(btn_reset)
        row.addStretch(1)
        box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        box.button(QDialogButtonBox.Ok).setText("确定")
        box.button(QDialogButtonBox.Cancel).setText("取消")
        box.accepted.connect(self.accept)
        box.rejected.connect(self.reject)
        row.addWidget(box)
        lay.addLayout(row)

    def _pick(self, kind: str, button: QPushButton) -> None:
        cur = QColor(self._palette[kind]["bg"])
        col = QColorDialog.getColor(cur, self, f"选择 {region_label(kind, self._palette)} 背景色")
        if not col.isValid():
            return
        self._palette[kind]["bg"] = col.name()
        self._palette[kind]["fg"] = _contrast_fg(col.name())
        button.setText(col.name())
        button.setStyleSheet(f"background:{col.name()}; border:1px solid #d1d5db;")
        r = LEGEND_ORDER.index(kind)
        prev = self.table.cellWidget(r, 2)
        if prev:
            prev.setStyleSheet(f"background:{col.name()}; color:{self._palette[kind]['fg']};")

    def _reset(self) -> None:
        self._palette = {k: dict(v) for k, v in DEFAULT_PALETTE.items()}
        for r, kind in enumerate(LEGEND_ORDER):
            conf = self._palette[kind]
            btn = self.table.cellWidget(r, 1)
            btn.setText(conf["bg"])
            btn.setStyleSheet(f"background:{conf['bg']}; border:1px solid #d1d5db;")
            prev = self.table.cellWidget(r, 2)
            prev.setStyleSheet(f"background:{conf['bg']}; color:{conf['fg']};")

    def palette(self) -> Dict[str, Dict[str, str]]:
        return self._palette


def _contrast_fg(bg_hex: str) -> str:
    c = QColor(bg_hex)
    lum = 0.299 * c.red() + 0.587 * c.green() + 0.114 * c.blue()
    return "#111827" if lum > 140 else "#f9fafb"


class BatchExportDialog(QDialog):
    """多基因批量导出。"""

    def __init__(self, gene_ids: Sequence[str], outdir: str,
                 upstream: int, downstream: int, parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("批量导出")
        self.resize(520, 460)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(12, 12, 12, 10)
        lay.setSpacing(8)

        lab = QLabel("每行一个 Gene ID（可先用搜索结果填入）：")
        lab.setStyleSheet("color:#374151; font-size:11px;")
        lay.addWidget(lab)
        self.edit = QPlainTextEdit()
        self.edit.setPlainText("\n".join(gene_ids))
        self.edit.setFont(mono_font(9))
        lay.addWidget(self.edit, 1)

        box = QGroupBox("参数")
        form = QFormLayout(box)
        self.spin_up = QSpinBox()
        self.spin_up.setRange(0, 1_000_000)
        self.spin_up.setSingleStep(500)
        self.spin_up.setValue(upstream)
        self.spin_dn = QSpinBox()
        self.spin_dn.setRange(0, 1_000_000)
        self.spin_dn.setSingleStep(500)
        self.spin_dn.setValue(downstream)
        form.addRow("上游 (bp)", self.spin_up)
        form.addRow("下游 (bp)", self.spin_dn)
        self.chk_all_tx = QCheckBox("导出该基因的全部转录本")
        self.chk_table = QCheckBox("同时导出片段明细表 (CSV)")
        self.chk_table.setChecked(True)
        form.addRow(self.chk_all_tx)
        form.addRow(self.chk_table)
        lay.addWidget(box)

        row = QHBoxLayout()
        self.dir_edit = QLineEdit(outdir)
        btn_dir = QPushButton("输出目录…")
        btn_dir.clicked.connect(self._pick_dir)
        row.addWidget(self.dir_edit, 1)
        row.addWidget(btn_dir)
        lay.addLayout(row)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.button(QDialogButtonBox.Ok).setText("开始导出")
        buttons.button(QDialogButtonBox.Cancel).setText("取消")
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        lay.addWidget(buttons)

    def _pick_dir(self) -> None:
        d = QFileDialog.getExistingDirectory(self, "选择输出目录", self.dir_edit.text())
        if d:
            self.dir_edit.setText(d)

    def values(self) -> dict:
        ids = [line.strip() for line in self.edit.toPlainText().splitlines() if line.strip()]
        return {
            "gene_ids": ids,
            "outdir": self.dir_edit.text().strip(),
            "upstream": self.spin_up.value(),
            "downstream": self.spin_dn.value(),
            "all_transcripts": self.chk_all_tx.isChecked(),
            "write_table": self.chk_table.isChecked(),
        }


class SequenceDocDialog(QDialog):
    """导出带颜色的序列文档（Word ``.docx`` / 网页 ``.html``）。"""

    FORMATS = [("Word 文档 (*.docx)", ".docx"), ("网页 (*.html)", ".html")]

    def __init__(self, extraction: GeneExtraction, preset_panels: Sequence[str],
                 default_path: str, bpl: int = 60, font_pt: float = 8.5,
                 palette: Optional[Dict[str, Dict[str, str]]] = None,
                 include_coord: bool = False, parent=None) -> None:
        super().__init__(parent)
        self.extraction = extraction
        self.palette_map = ensure_palette(palette)
        self.setWindowTitle("导出彩色序列文档")
        self.resize(620, 660)

        lay = QVBoxLayout(self)
        lay.setContentsMargins(12, 12, 12, 10)
        lay.setSpacing(8)

        tip = QLabel("导出的序列带底色（与界面上看到的一致），可直接在 Word 里编辑、"
                     "加批注或截图使用。")
        tip.setWordWrap(True)
        tip.setStyleSheet("color:#4b5563; font-size:11px;")
        lay.addWidget(tip)

        # ---------------- 面板 ---------------- #
        box_panels = QGroupBox("包含哪些面板")
        grid = QGridLayout(box_panels)
        grid.setContentsMargins(10, 8, 10, 8)
        grid.setHorizontalSpacing(12)
        grid.setVerticalSpacing(4)
        self.panel_checks: Dict[str, QCheckBox] = {}
        keys = [k for k in PANEL_ORDER if k in extraction.panels]
        for i, key in enumerate(keys):
            panel = extraction.panels[key]
            chk = QCheckBox(f"{PANEL_TITLES.get(key, key)}（{panel.length:,} bp）")
            chk.setChecked(key in preset_panels and not panel.is_empty())
            if panel.is_empty():
                chk.setEnabled(False)
                chk.setToolTip("这个面板没有内容")
            chk.toggled.connect(self._update_hint)
            grid.addWidget(chk, i // 2, i % 2)
            self.panel_checks[key] = chk
        row = QHBoxLayout()
        btn_all = QPushButton("全选")
        btn_none = QPushButton("全不选")
        btn_all.clicked.connect(lambda: self._set_all(True))
        btn_none.clicked.connect(lambda: self._set_all(False))
        row.addWidget(btn_all)
        row.addWidget(btn_none)
        row.addStretch(1)
        grid.addLayout(row, (len(keys) + 1) // 2, 0, 1, 2)
        lay.addWidget(box_panels)

        # ---------------- 参数 ---------------- #
        box = QGroupBox("排版与着色")
        form = QFormLayout(box)
        form.setContentsMargins(10, 8, 10, 8)
        form.setVerticalSpacing(5)
        self.cmb_mode = QComboBox()
        for key, label in MODE_LABELS.items():
            self.cmb_mode.addItem(label, key)
        self.cmb_mode.setCurrentIndex(max(0, self.cmb_mode.findData(MODE_REGION)))
        self.cmb_mode.currentIndexChanged.connect(self._update_hint)
        form.addRow("着色方式", self.cmb_mode)

        self.spin_bpl = QSpinBox()
        self.spin_bpl.setRange(20, 200)
        self.spin_bpl.setSingleStep(10)
        self.spin_bpl.setValue(max(20, min(200, int(bpl or DEFAULT_BPL))))
        self.spin_bpl.valueChanged.connect(self._update_hint)
        form.addRow("每行碱基数", self.spin_bpl)

        self.spin_font = QDoubleSpinBox()
        self.spin_font.setRange(6.0, 14.0)
        self.spin_font.setSingleStep(0.5)
        self.spin_font.setDecimals(1)
        self.spin_font.setValue(float(font_pt or 8.5))
        form.addRow("字号 (pt)", self.spin_font)

        self.chk_head = QCheckBox("包含标题、基因信息与颜色图例")
        self.chk_head.setChecked(True)
        self.chk_head.toggled.connect(self._update_hint)
        form.addRow(self.chk_head)

        self.chk_coord = QCheckBox("每行左侧显示坐标列")
        self.chk_coord.setChecked(bool(include_coord))
        self.chk_coord.setToolTip(
            "关闭（默认）时文档里只有序列本身，左侧不会出现位置数字，"
            "贴到论文或报告里更干净；需要按位置核对时再打开。")
        form.addRow(self.chk_coord)
        lay.addWidget(box)

        # ---------------- 输出 ---------------- #
        box_out = QGroupBox("输出")
        out = QFormLayout(box_out)
        out.setContentsMargins(10, 8, 10, 8)
        out.setVerticalSpacing(5)
        self.cmb_format = QComboBox()
        for label, _suffix in self.FORMATS:
            self.cmb_format.addItem(label)
        self.cmb_format.currentIndexChanged.connect(self._on_format_changed)
        out.addRow("格式", self.cmb_format)
        prow = QHBoxLayout()
        self.path_edit = QLineEdit(default_path)
        btn_path = QPushButton("浏览…")
        btn_path.clicked.connect(self._pick_path)
        prow.addWidget(self.path_edit, 1)
        prow.addWidget(btn_path)
        out.addRow("保存到", prow)
        lay.addWidget(box_out)

        self.hint = QLabel("")
        self.hint.setWordWrap(True)
        self.hint.setStyleSheet("color:#6b7280; font-size:11px;")
        lay.addWidget(self.hint)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.button(QDialogButtonBox.Ok).setText("导出")
        buttons.button(QDialogButtonBox.Cancel).setText("取消")
        buttons.accepted.connect(self._on_ok)
        buttons.rejected.connect(self.reject)
        lay.addWidget(buttons)
        self._update_hint()

    # ------------------------------------------------------------------ #
    def _set_all(self, on: bool) -> None:
        for chk in self.panel_checks.values():
            if chk.isEnabled():
                chk.setChecked(on)

    def selected_panels(self) -> List[str]:
        return [k for k, chk in self.panel_checks.items() if chk.isChecked()]

    def mode(self) -> str:
        return str(self.cmb_mode.currentData() or MODE_REGION)

    def _on_format_changed(self) -> None:
        suffix = self.FORMATS[max(0, self.cmb_format.currentIndex())][1]
        stem = os.path.splitext(self.path_edit.text().strip())[0]
        if stem:
            self.path_edit.setText(stem + suffix)
        self._update_hint()

    def _pick_path(self) -> None:
        suffix = self.FORMATS[max(0, self.cmb_format.currentIndex())][1]
        flt = "Word 文档 (*.docx)" if suffix == ".docx" else "网页 (*.html)"
        current = self.path_edit.text().strip()
        if current and not current.lower().endswith(suffix):
            current = os.path.splitext(current)[0] + suffix
        path, _ = QFileDialog.getSaveFileName(self, "保存彩色序列文档",
                                              current or "sequence.docx", flt)
        if path:
            if not path.lower().endswith((".docx", ".html", ".htm")):
                path += suffix
            self.path_edit.setText(path)

    def _update_hint(self) -> None:
        keys = self.selected_panels()
        if not keys:
            self.hint.setText("请至少选择一个面板。")
            self.hint.setStyleSheet("color:#b45309; font-size:11px;")
            return
        est = estimate(self.extraction, keys, self.spin_bpl.value(), self.mode())
        text = (f"共 {est['panels']} 个面板 / {est['bases']:,} bp，"
                f"约 {est['lines']:,} 行（{est['pages']:,} 页），"
                f"Word 文件约 {human_bytes(est['docx_bytes'])}")
        color = "#6b7280"
        if est["bases"] > WARN_BASES:
            text += ("\n⚠ 序列很长，Word 打开与导出都会比较慢；"
                     "建议改用「按区域着色」或选更少的面板。")
            color = "#b45309"
        self.hint.setText(text)
        self.hint.setStyleSheet(f"color:{color}; font-size:11px;")

    def _on_ok(self) -> None:
        if not self.selected_panels():
            QMessageBox.warning(self, "没有选择面板", "请至少勾选一个面板。")
            return
        if not self.path_edit.text().strip():
            QMessageBox.warning(self, "没有指定保存位置", "请填写要保存的文件路径。")
            return
        self.accept()

    def values(self) -> dict:
        path = self.path_edit.text().strip()
        suffix = self.FORMATS[max(0, self.cmb_format.currentIndex())][1]
        if not path.lower().endswith((".docx", ".html", ".htm")):
            path += suffix
        return {
            "path": path,
            "panels": self.selected_panels(),
            "mode": self.mode(),
            "bpl": self.spin_bpl.value(),
            "font_pt": self.spin_font.value(),
            "include_header": self.chk_head.isChecked(),
            "include_coord": self.chk_coord.isChecked(),
        }


class DatasetExportDialog(QDialog):
    """把当前分别导入的几个文件打包成一个"下次能直接导入"的数据集压缩包。"""

    def __init__(self, paths: Dict[str, str], name: str = "", note: str = "",
                 default_path: str = "", parent=None) -> None:
        super().__init__(parent)
        from ..bundle import REQUIRED_ROLES, ROLE_LABELS, ROLE_ORDER
        self._required = set(REQUIRED_ROLES)
        self.setWindowTitle("导出数据集压缩包")
        self.resize(660, 480)

        lay = QVBoxLayout(self)
        lay.setContentsMargins(12, 12, 12, 10)
        lay.setSpacing(8)

        tip = QLabel("把当前这套文件打成一个 zip：里面带一份 gge_dataset.json 清单，"
                     "下次「文件 → 导入压缩包」时会按清单直接认回角色，不用再猜。")
        tip.setWordWrap(True)
        tip.setStyleSheet("color:#4b5563; font-size:11px;")
        lay.addWidget(tip)

        box = QGroupBox("数据集")
        form = QFormLayout(box)
        form.setContentsMargins(10, 8, 10, 8)
        form.setVerticalSpacing(5)
        self.name_edit = QLineEdit(name)
        self.note_edit = QLineEdit(note)
        self.note_edit.setPlaceholderText("可留空，例如：GCF_019924925.1 草鱼 HZGC01 组装")
        form.addRow("名称", self.name_edit)
        form.addRow("备注", self.note_edit)
        lay.addWidget(box)

        self.table = QTableWidget(len(ROLE_ORDER), 4)
        self.table.setHorizontalHeaderLabels(["包含", "角色", "文件", "大小"])
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setSelectionMode(QAbstractItemView.NoSelection)
        head = self.table.horizontalHeader()
        head.setSectionResizeMode(0, QHeaderView.ResizeToContents)
        head.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        head.setSectionResizeMode(2, QHeaderView.Stretch)
        head.setSectionResizeMode(3, QHeaderView.ResizeToContents)
        self._rows: Dict[str, int] = {}
        for r, role in enumerate(ROLE_ORDER):
            path = paths.get(role, "")
            chk = QTableWidgetItem()
            chk.setFlags(Qt.ItemIsEnabled | Qt.ItemIsUserCheckable)
            chk.setCheckState(Qt.Checked if path else Qt.Unchecked)
            self.table.setItem(r, 0, chk)
            self.table.setItem(r, 1, QTableWidgetItem(ROLE_LABELS[role]))
            item = QTableWidgetItem(os.path.basename(path) if path else "（未加载）")
            item.setToolTip(path)
            self.table.setItem(r, 2, item)
            size = 0
            if path and os.path.isfile(path):
                try:
                    size = os.path.getsize(path)
                except OSError:
                    size = 0
            self.table.setItem(r, 3, QTableWidgetItem(human_bytes(size) if size else "—"))
            self._rows[role] = r
            if role in self._required and path:
                chk.setToolTip("基因组与注释是必填项，不能取消")
        lay.addWidget(self.table, 1)

        prow = QHBoxLayout()
        prow.addWidget(QLabel("保存到"))
        self.path_edit = QLineEdit(default_path)
        btn_path = QPushButton("浏览…")
        btn_path.clicked.connect(self._pick_path)
        prow.addWidget(self.path_edit, 1)
        prow.addWidget(btn_path)
        lay.addLayout(prow)

        self.hint = QLabel("")
        self.hint.setWordWrap(True)
        self.hint.setStyleSheet("color:#6b7280; font-size:11px;")
        lay.addWidget(self.hint)
        for r in range(len(ROLE_ORDER)):
            self.table.item(r, 0).setCheckState(self.table.item(r, 0).checkState())
        self.table.itemChanged.connect(lambda _i: self._update_hint())
        self._update_hint()

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.button(QDialogButtonBox.Ok).setText("开始打包")
        buttons.button(QDialogButtonBox.Cancel).setText("取消")
        buttons.accepted.connect(self._on_ok)
        buttons.rejected.connect(self.reject)
        lay.addWidget(buttons)

    # ------------------------------------------------------------------ #
    def _checked_roles(self) -> List[str]:
        from ..bundle import ROLE_ORDER

        out = []
        for role in ROLE_ORDER:
            item = self.table.item(self._rows[role], 0)
            if item.checkState() == Qt.Checked:
                out.append(role)
        return out

    def _update_hint(self) -> None:
        from ..bundle import ROLE_LABELS, ROLE_ORDER

        roles = self._checked_roles()
        total = 0
        for role in roles:
            item = self.table.item(self._rows[role], 2)
            path = item.toolTip()
            if path and os.path.isfile(path):
                total += os.path.getsize(path)
        missing = [ROLE_LABELS[r] for r in self._required if r not in roles]
        if missing:
            self.hint.setText("缺少必填文件：" + "、".join(missing) +
                              "。请先加载，或改用在左侧逐个选择文件。")
            self.hint.setStyleSheet("color:#b91c1c; font-size:11px;")
            return
        self.hint.setText(
            f"共 {len(roles)} 个文件，原始 {human_bytes(total)}；"
            f"压缩后通常只剩 25%–35%（FASTA 压得很好）。"
            f"压缩过程按字节推进，可以在进度条旁随时取消。")
        self.hint.setStyleSheet("color:#6b7280; font-size:11px;")

    def _pick_path(self) -> None:
        path, _ = QFileDialog.getSaveFileName(
            self, "保存数据集压缩包", self.path_edit.text().strip() or "dataset.zip",
            "ZIP 压缩包 (*.zip)")
        if path:
            if not path.lower().endswith(".zip"):
                path += ".zip"
            self.path_edit.setText(path)

    def _on_ok(self) -> None:
        from ..bundle import ROLE_LABELS

        roles = self._checked_roles()
        missing = [ROLE_LABELS[r] for r in self._required if r not in roles]
        if missing:
            QMessageBox.warning(self, "缺少必填文件",
                                "基因组与注释是必须的：缺少 " + "、".join(missing))
            return
        if not self.path_edit.text().strip():
            QMessageBox.warning(self, "没有指定保存位置", "请填写压缩包的保存路径。")
            return
        self.accept()

    def values(self) -> dict:
        paths: Dict[str, str] = {}
        for role, row in self._rows.items():
            if self.table.item(row, 0).checkState() == Qt.Checked:
                paths[role] = self.table.item(row, 2).toolTip()
        path = self.path_edit.text().strip()
        if not path.lower().endswith(".zip"):
            path += ".zip"
        return {
            "path": path,
            "name": self.name_edit.text().strip() or "dataset",
            "note": self.note_edit.text().strip(),
            "paths": paths,
        }


class AboutDialog(QDialog):
    def __init__(self, text: str, parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("关于 / 使用说明")
        self.resize(640, 560)
        lay = QVBoxLayout(self)
        view = QTextEdit()
        view.setReadOnly(True)
        view.setHtml(text)
        lay.addWidget(view)
        box = QDialogButtonBox(QDialogButtonBox.Close)
        box.rejected.connect(self.reject)
        box.accepted.connect(self.accept)
        lay.addWidget(box)
