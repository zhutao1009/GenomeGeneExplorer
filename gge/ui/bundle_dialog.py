"""压缩包导入对话框 + 数据集库对话框。

导入对话框把"自动识别"的结果**摊开给用户看**：每个文件出现在一行，写明判成什么、
依据是什么，并且可以逐行改。程序从不擅自决定"哪个文件是基因组"就一路跑下去——
必填项缺失时"开始导入"按钮是灰的。
"""

from __future__ import annotations

import os
from typing import Dict, List, Optional

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (QAbstractItemView, QComboBox, QDialog,
                               QDialogButtonBox, QFrame, QHBoxLayout, QHeaderView,
                               QInputDialog, QLabel, QLineEdit, QMessageBox,
                               QPushButton, QTableWidget, QTableWidgetItem,
                               QTextEdit, QVBoxLayout, QWidget)

from ..bundle import (REQUIRED_ROLES, ROLE_LABELS, ROLE_ORDER, BundleScan,
                      describe_roles)
from ..coords import human_bytes
from ..library import Dataset, Library
from .theme import mono_font, ui_font

_ROLE_ITEMS = [("", "未识别 / 不导入")] + [(r, ROLE_LABELS[r]) for r in ROLE_ORDER]


def _archive_stem(path: str) -> str:
    name = os.path.basename(path or "")
    for ext in (".tar.gz", ".tar.bz2", ".tar.xz", ".tgz", ".tbz2", ".txz", ".zip", ".tar"):
        if name.lower().endswith(ext):
            return name[: -len(ext)]
    return name


class BundleImportDialog(QDialog):
    """确认压缩包里的文件识别结果。"""

    def __init__(self, scan: BundleScan, parent=None,
                 default_name: str = "") -> None:
        super().__init__(parent)
        self.scan = scan
        self.setWindowTitle("导入压缩包")
        self.resize(1180, 520)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(12, 10, 12, 10)
        lay.setSpacing(8)

        # ---------------- 顶部摘要 ---------------- #
        head = QFrame()
        head.setStyleSheet("QFrame{background:#f8fafc;border:1px solid #e5e7eb;"
                           "border-radius:6px;}")
        hv = QVBoxLayout(head)
        hv.setContentsMargins(10, 8, 10, 8)
        hv.setSpacing(4)

        title = QLabel(f"{os.path.basename(scan.path)}")
        f = ui_font(10)
        f.setBold(True)
        title.setFont(f)
        hv.addWidget(title)

        meta = QLabel(f"{scan.kind.upper()} 压缩包 · {human_bytes(scan.archive_size)} · "
                      f"解压后约 {human_bytes(scan.total_uncompressed)} · "
                      f"{scan.member_count} 个文件（已跳过 {len(scan.ignored)} 个无关文件）")
        meta.setStyleSheet("color:#6b7280; font-size:11px;")
        meta.setWordWrap(True)
        hv.addWidget(meta)

        row = QHBoxLayout()
        row.setSpacing(6)
        lab = QLabel("数据集名称")
        lab.setStyleSheet("color:#374151; font-size:11px;")
        # 本程序导出的包带清单，清单里的名字比压缩包文件名更准确
        self.name_edit = QLineEdit(default_name or scan.name_hint()
                                   or _archive_stem(scan.path))
        self.name_edit.setPlaceholderText("给这套数据起个名字，例如 斑马鱼 GRCz11")
        row.addWidget(lab)
        row.addWidget(self.name_edit, 1)
        hv.addLayout(row)
        lay.addWidget(head)

        # ---------------- 识别结果表 ---------------- #
        cols = ["文件", "格式", "内容概况", "识别为", "判断依据"]
        self.table = QTableWidget(0, len(cols))
        self.table.setHorizontalHeaderLabels(cols)
        self.table.verticalHeader().setVisible(False)
        self.table.setAlternatingRowColors(True)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        hh = self.table.horizontalHeader()
        # 判断依据最长、最需要空间，所以让最后一列吃掉剩余宽度；
        # 文件名与内容概况给固定初值（用户仍可拖动），否则两列 Stretch 会把依据挤成省略号。
        hh.setSectionResizeMode(0, QHeaderView.Interactive)
        self.table.setColumnWidth(0, 280)
        hh.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        hh.setSectionResizeMode(2, QHeaderView.Interactive)
        self.table.setColumnWidth(2, 225)
        hh.setSectionResizeMode(3, QHeaderView.Fixed)
        self.table.setColumnWidth(3, 120)
        hh.setSectionResizeMode(4, QHeaderView.Stretch)
        lay.addWidget(self.table, 1)

        self._combos: List[QComboBox] = []
        self._fill_table()

        # ---------------- 底部 ---------------- #
        self.summary = QLabel("")
        self.summary.setStyleSheet("font-size:11px;")
        self.summary.setWordWrap(True)
        lay.addWidget(self.summary)

        self.warn_label = QLabel("")
        self.warn_label.setStyleSheet("color:#b45309; font-size:10px;")
        self.warn_label.setWordWrap(True)
        lay.addWidget(self.warn_label)

        bottom = QHBoxLayout()
        btn_report = QPushButton("查看扫描详情")
        btn_report.setToolTip("列出所有文件（含被跳过的）及判定依据")
        btn_report.clicked.connect(self._show_report)
        bottom.addWidget(btn_report)
        bottom.addStretch(1)
        self.buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        self.buttons.button(QDialogButtonBox.Ok).setText("开始导入")
        self.buttons.button(QDialogButtonBox.Cancel).setText("取消")
        self.buttons.accepted.connect(self._on_ok)
        self.buttons.rejected.connect(self.reject)
        bottom.addWidget(self.buttons)
        lay.addLayout(bottom)

        self._refresh_summary()
        self._refresh_warnings()

    # ------------------------------------------------------------------ #
    def _fill_table(self) -> None:
        cands = self.scan.candidates
        self.table.setRowCount(len(cands))
        for r, cand in enumerate(cands):
            name_item = QTableWidgetItem(cand.member)
            name_item.setToolTip(cand.member)
            self.table.setItem(r, 0, name_item)

            fmt_item = QTableWidgetItem(cand.fmt or "—")
            fmt_item.setFont(mono_font(8))
            self.table.setItem(r, 1, fmt_item)

            stat_item = QTableWidgetItem(cand.stat_text())
            stat_item.setToolTip(cand.stat_text())
            self.table.setItem(r, 2, stat_item)

            combo = QComboBox()
            for value, label in _ROLE_ITEMS:
                combo.addItem(label, value)
            idx = combo.findData(cand.role)
            combo.setCurrentIndex(max(0, idx))
            combo.setProperty("member", cand.member)
            combo.currentIndexChanged.connect(
                lambda _i, c=combo: self._on_role_changed(c))
            self.table.setCellWidget(r, 3, combo)
            self._combos.append(combo)

            ev = QTableWidgetItem(cand.evidence_text())
            ev.setToolTip(cand.evidence_text())
            if cand.ambiguous:
                ev.setForeground(QColor("#b45309"))
            self.table.setItem(r, 4, ev)
            self.table.setRowHeight(r, 26)

    def _role_of(self, combo: QComboBox) -> str:
        return combo.currentData() or ""

    def _on_role_changed(self, combo: QComboBox) -> None:
        role = self._role_of(combo)
        if role:
            # 一个角色只能对应一个文件：把占用同一角色的其它行清掉
            for other in self._combos:
                if other is not combo and self._role_of(other) == role:
                    other.blockSignals(True)
                    other.setCurrentIndex(0)
                    other.blockSignals(False)
        self._refresh_summary()

    def _refresh_summary(self) -> None:
        got = self.roles()
        parts = []
        for role in ROLE_ORDER:
            ok = bool(got.get(role))
            mark = "✓" if ok else ("×" if role in REQUIRED_ROLES else "—")
            color = "#15803d" if ok else ("#b91c1c" if role in REQUIRED_ROLES else "#9ca3af")
            parts.append(f"<span style='color:{color}'>{ROLE_LABELS[role]} {mark}</span>")
        missing = [r for r in REQUIRED_ROLES if not got.get(r)]
        extra = ""
        if missing:
            extra = ("<br><span style='color:#b91c1c'>还缺："
                     + "、".join(ROLE_LABELS[r] for r in missing)
                     + " —— 请在下拉框里手动指定</span>")
        self.summary.setText("识别结果：" + " · ".join(parts) + extra)
        ok = self.buttons.button(QDialogButtonBox.Ok)
        ok.setEnabled(not missing)
        ok.setToolTip("必填项未指定" if missing else "解包并加入数据集库")

    def _refresh_warnings(self) -> None:
        warns = list(self.scan.warnings)
        if not warns:
            self.warn_label.setText("")
            return
        self.warn_label.setText(f"⚠ {len(warns)} 条提示：" + warns[0]
                                + ("…（详见「查看扫描详情」）" if len(warns) > 1 else ""))
        self.warn_label.setToolTip("\n".join(warns))

    def _show_report(self) -> None:
        dlg = QDialog(self)
        dlg.setWindowTitle("扫描详情")
        dlg.resize(820, 520)
        v = QVBoxLayout(dlg)
        view = QTextEdit()
        view.setReadOnly(True)
        view.setFont(mono_font(9))
        view.setPlainText(self.scan.text_report())
        v.addWidget(view, 1)
        box = QDialogButtonBox(QDialogButtonBox.Close)
        box.button(QDialogButtonBox.Close).setText("关闭")
        box.rejected.connect(dlg.reject)
        box.accepted.connect(dlg.accept)
        v.addWidget(box)
        dlg.exec()

    # ------------------------------------------------------------------ #
    def roles(self) -> Dict[str, str]:
        """role -> 包内路径（只含被指定的）。"""
        out: Dict[str, str] = {}
        for combo in self._combos:
            role = self._role_of(combo)
            member = combo.property("member")
            if role and member and role not in out:
                out[role] = member
        return out

    def dataset_name(self) -> str:
        return self.name_edit.text().strip() or _archive_stem(self.scan.path)

    def _on_ok(self) -> None:
        missing = [r for r in REQUIRED_ROLES if r not in self.roles()]
        if missing:
            QMessageBox.warning(self, "还缺必填文件",
                                "请为以下角色指定文件：" +
                                "、".join(ROLE_LABELS[r] for r in missing))
            return
        self.accept()


class LibraryDialog(QDialog):
    """数据集库浏览：载入、重命名、移除。"""

    def __init__(self, library: Library, parent=None,
                 current_id: int = 0) -> None:
        super().__init__(parent)
        self.lib = library
        self.wants_import = False
        self._datasets: List[Dataset] = []
        self.setWindowTitle("数据集库")
        self.resize(940, 520)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(12, 10, 12, 10)
        lay.setSpacing(8)

        top = QHBoxLayout()
        self.search = QLineEdit()
        self.search.setPlaceholderText("按名称或来源筛选…")
        self.search.setClearButtonEnabled(True)
        self.search.textChanged.connect(self._reload)
        btn_refresh = QPushButton("刷新")
        btn_refresh.clicked.connect(self._reload)
        btn_import = QPushButton("导入压缩包…")
        btn_import.setObjectName("primary")
        btn_import.clicked.connect(self._on_import)
        top.addWidget(self.search, 1)
        top.addWidget(btn_refresh)
        top.addWidget(btn_import)
        lay.addLayout(top)

        cols = ["名称", "概况", "注释", "来源", "最近使用", "状态"]
        self.table = QTableWidget(0, len(cols))
        self.table.setHorizontalHeaderLabels(cols)
        self.table.verticalHeader().setVisible(False)
        self.table.setAlternatingRowColors(True)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        hh = self.table.horizontalHeader()
        # 「概况」最需要弹性（内容长），名称给固定初值，其余按内容收紧
        hh.setSectionResizeMode(0, QHeaderView.Interactive)
        self.table.setColumnWidth(0, 190)
        hh.setSectionResizeMode(1, QHeaderView.Stretch)
        for i in (2, 3, 4, 5):
            hh.setSectionResizeMode(i, QHeaderView.ResizeToContents)
        self.table.doubleClicked.connect(lambda _i: self._load_selected())
        self.table.itemSelectionChanged.connect(self._update_buttons)
        lay.addWidget(self.table, 1)

        self.hint = QLabel("")
        self.hint.setStyleSheet("color:#9ca3af; font-size:10px;")
        self.hint.setWordWrap(True)
        lay.addWidget(self.hint)

        bottom = QHBoxLayout()
        self.btn_load = QPushButton("载入")
        self.btn_load.setObjectName("primary")
        self.btn_load.clicked.connect(self._load_selected)
        self.btn_rename = QPushButton("重命名")
        self.btn_rename.clicked.connect(self._rename_selected)
        self.btn_remove = QPushButton("移除")
        self.btn_remove.clicked.connect(self._remove_selected)
        self.btn_dir = QPushButton("打开数据目录")
        self.btn_dir.clicked.connect(self._open_dir)
        bottom.addWidget(self.btn_load)
        bottom.addWidget(self.btn_rename)
        bottom.addWidget(self.btn_remove)
        bottom.addWidget(self.btn_dir)
        bottom.addStretch(1)
        box = QDialogButtonBox(QDialogButtonBox.Close)
        box.button(QDialogButtonBox.Close).setText("关闭")
        box.rejected.connect(self.reject)
        bottom.addWidget(box)
        lay.addLayout(bottom)

        self._reload()
        if current_id:
            self._select_id(current_id)
        elif self._datasets:
            # 默认选中第一行：否则「载入/重命名/移除」全是灰的，看着像坏了
            self.table.selectRow(0)

    # ------------------------------------------------------------------ #
    def _reload(self) -> None:
        self._datasets = self.lib.list(self.search.text().strip())
        self.table.setRowCount(len(self._datasets))
        for r, ds in enumerate(self._datasets):
            for c, text in enumerate(ds.row()):
                item = QTableWidgetItem(text)
                item.setToolTip(f"{ds.name}\n{ds.summary()}\n"
                                f"基因组：{ds.genome}\n注释：{ds.annotation}\n"
                                f"来源：{ds.source}")
                if c == 0:
                    f = item.font()
                    f.setBold(True)
                    item.setFont(f)
                self.table.setItem(r, c, item)
            missing = ds.missing_required()
            broken = ds.missing()
            status = "就绪" if not broken else f"缺 {len(broken)} 个文件"
            item = QTableWidgetItem(status)
            item.setForeground(QColor("#15803d" if not broken else "#b91c1c"))
            if broken:
                item.setToolTip("找不到：" + "\n".join(ds.path(r) for r in broken))
            elif missing:
                item.setToolTip("缺少必需角色：" + "、".join(missing))
            self.table.setItem(r, 5, item)
            self.table.setRowHeight(r, 26)
        self.hint.setText(f"共 {len(self._datasets)} 套数据。"
                          "双击任一行即可载入；数据目录："
                          f"{self.lib.path}")
        self._update_buttons()

    def _select_id(self, dataset_id: int) -> None:
        for r, ds in enumerate(self._datasets):
            if ds.id == dataset_id:
                self.table.selectRow(r)
                return

    def _current(self) -> Optional[Dataset]:
        r = self.table.currentRow()
        if 0 <= r < len(self._datasets):
            return self._datasets[r]
        return None

    def _update_buttons(self) -> None:
        ds = self._current()
        for b in (self.btn_load, self.btn_rename, self.btn_remove, self.btn_dir):
            b.setEnabled(ds is not None)

    # ------------------------------------------------------------------ #
    def _on_import(self) -> None:
        self.wants_import = True
        self.accept()

    def _load_selected(self) -> None:
        ds = self._current()
        if ds is None:
            return
        if ds.missing_required():
            QMessageBox.warning(
                self, "文件缺失",
                "这套数据的文件已经不在原来的位置：\n" +
                "\n".join(ds.path(r) or ROLE_LABELS[r] for r in ds.missing_required()) +
                "\n\n可能被移动或删除了，请重新导入压缩包。")
            return
        if ds.missing():
            QMessageBox.information(
                self, "部分参考文件缺失",
                "以下参考序列找不到，将无法用于校验（不影响提取）：\n" +
                "\n".join(ds.path(r) for r in ds.missing()))
        self.accept()

    def _rename_selected(self) -> None:
        ds = self._current()
        if ds is None:
            return
        text, ok = QInputDialog.getText(self, "重命名数据集", "新名称：", text=ds.name)
        if ok and text.strip():
            self.lib.rename(ds.id, text.strip())
            self._reload()

    def _remove_selected(self) -> None:
        ds = self._current()
        if ds is None:
            return
        owned = self.lib.owns_files(ds)
        detail = ("会同时删除库目录内的解包文件（原始压缩包不受影响）。"
                  if owned else "只删除库里的记录，磁盘上的文件保持不动。")
        ret = QMessageBox.question(
            self, "移除数据集",
            f"确定移除「{ds.name}」吗？\n{detail}",
            QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if ret != QMessageBox.Yes:
            return
        removed = self.lib.remove(ds.id, delete_files=owned)
        self._reload()
        if removed:
            self.hint.setText(f"已移除「{ds.name}」，删除目录：{removed[0]}")

    def _open_dir(self) -> None:
        ds = self._current()
        if ds is None:
            return
        path = os.path.dirname(ds.genome or ds.annotation)
        if path and os.path.isdir(path):
            try:
                os.startfile(path)      # type: ignore[attr-defined]
            except Exception:
                QMessageBox.information(self, "数据目录", path)
        else:
            QMessageBox.information(self, "数据目录", "文件已不在原位置。")

    def selected_dataset(self) -> Optional[Dataset]:
        return self._current()
