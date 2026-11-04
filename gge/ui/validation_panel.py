"""校验结果面板（mRNA / CDS / 蛋白）。"""

from __future__ import annotations

from typing import Dict, Optional

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (QComboBox, QFrame, QGridLayout, QHBoxLayout,
                               QLabel, QPushButton, QScrollArea, QSizePolicy,
                               QVBoxLayout, QWidget)

from ..validate import STATUS_CN, CheckResult
from .theme import ui_font


class ResultCard(QFrame):
    """单个校验结果卡片。"""

    detailRequested = Signal(str)
    manualChanged = Signal(str, str)

    def __init__(self, kind: str, parent=None) -> None:
        super().__init__(parent)
        self.kind = kind
        self._result: Optional[CheckResult] = None
        self.setObjectName("resultCard")
        self.setStyleSheet(
            "QFrame#resultCard { background:#ffffff; border:1px solid #e5e7eb;"
            " border-radius:6px; }"
        )
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(10, 8, 10, 8)
        lay.setSpacing(4)

        head = QHBoxLayout()
        self.title = QLabel(kind)
        f = ui_font(9)
        f.setBold(True)
        self.title.setFont(f)
        self.badge = QLabel("未校验")
        self.badge.setStyleSheet(
            "background:#6b7280; color:white; border-radius:7px; padding:1px 8px;"
            " font-size:10px; font-weight:600;"
        )
        head.addWidget(self.title)
        head.addWidget(self.badge)
        head.addStretch(1)
        self.btn_detail = QPushButton("详情")
        self.btn_detail.setFixedHeight(22)
        self.btn_detail.setStyleSheet("font-size:11px; padding:1px 8px;")
        self.btn_detail.setEnabled(False)
        head.addWidget(self.btn_detail)
        lay.addLayout(head)

        self.headline = QLabel("—")
        self.headline.setStyleSheet("color:#374151; font-size:11px;")
        self.headline.setWordWrap(True)
        lay.addWidget(self.headline)

        self.metrics = QGridLayout()
        self.metrics.setHorizontalSpacing(14)
        self.metrics.setVerticalSpacing(2)
        lay.addLayout(self.metrics)

        row = QHBoxLayout()
        lab = QLabel("参考")
        lab.setStyleSheet("color:#6b7280; font-size:10px;")
        self.combo = QComboBox()
        self.combo.setStyleSheet("font-size:10px; padding:0 4px;")
        self.combo.setMinimumWidth(120)
        self.combo.addItem("（自动匹配）", "")
        self.combo.currentIndexChanged.connect(self._on_combo)
        row.addWidget(lab)
        row.addWidget(self.combo, 1)
        lay.addLayout(row)

        self.btn_detail.clicked.connect(lambda: self.detailRequested.emit(self.kind))

    def _on_combo(self, _idx: int) -> None:
        self.manualChanged.emit(self.kind, self.combo.currentData() or "")

    def set_result(self, result: Optional[CheckResult], ref_ids=None) -> None:
        self._result = result
        self._blocking = True
        try:
            if result is None:
                self.badge.setText("未校验")
                self.badge.setStyleSheet(
                    "background:#6b7280;color:white;border-radius:7px;padding:1px 8px;"
                    "font-size:10px;font-weight:600;")
                self.headline.setText("尚无可显示的校验结果")
                self._clear_metrics()
                self.btn_detail.setEnabled(False)
                return
            self.badge.setText(result.status_cn)
            self.badge.setStyleSheet(
                f"background:{result.color};color:white;border-radius:7px;"
                f"padding:1px 8px;font-size:10px;font-weight:600;")
            self.headline.setText(result.headline())
            self._clear_metrics()
            keep = {"一致率", "不一致碱基", "插入（obs 多出）", "缺失（ref 多出）"}
            col = 0
            row = 0
            for key, val in result.metrics():
                if key not in keep:
                    continue
                k = QLabel(key)
                k.setStyleSheet("color:#6b7280; font-size:10px;")
                v = QLabel(val)
                v.setStyleSheet("color:#111827; font-size:11px; font-weight:600;")
                self.metrics.addWidget(k, row, col * 2)
                self.metrics.addWidget(v, row, col * 2 + 1)
                col += 1
                if col == 2:
                    col = 0
                    row += 1
            self.btn_detail.setEnabled(True)

            self.combo.blockSignals(True)
            self.combo.clear()
            self.combo.addItem("（自动匹配）", "")
            for rid in (ref_ids or []):
                self.combo.addItem(rid[:60], rid)
            if result.manual_id:
                idx = self.combo.findData(result.manual_id)
                if idx >= 0:
                    self.combo.setCurrentIndex(idx)
            self.combo.blockSignals(False)
        finally:
            self._blocking = False

    def _clear_metrics(self) -> None:
        while self.metrics.count():
            item = self.metrics.takeAt(0)
            w = item.widget()
            if w:
                w.deleteLater()


class ValidationPanel(QFrame):
    """mRNA / CDS / 蛋白三张卡片。"""

    detailRequested = Signal(str)
    manualChanged = Signal(str, str)
    revalidateRequested = Signal()

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.setObjectName("validationPanel")
        self.setStyleSheet(
            "QFrame#validationPanel { background:#ffffff; border:1px solid #e5e7eb;"
            " border-radius:6px; }"
        )
        lay = QVBoxLayout(self)
        lay.setContentsMargins(8, 6, 8, 8)
        lay.setSpacing(6)

        head = QHBoxLayout()
        title = QLabel("序列校验")
        f = ui_font(9)
        f.setBold(True)
        title.setFont(f)
        head.addWidget(title)
        self.hint = QLabel("加载参考 mRNA/CDS 后自动校验")
        self.hint.setStyleSheet("color:#9ca3af; font-size:10px;")
        head.addWidget(self.hint)
        head.addStretch(1)
        self.btn_revalidate = QPushButton("重新校验")
        self.btn_revalidate.setFixedHeight(22)
        self.btn_revalidate.setStyleSheet("font-size:11px; padding:1px 8px;")
        self.btn_revalidate.clicked.connect(self.revalidateRequested)
        head.addWidget(self.btn_revalidate)
        lay.addLayout(head)

        # 卡片放进滚动区：窗口变小时滚动，而不是把内容压没或顶高整个窗口
        self.cards: Dict[str, ResultCard] = {}
        host = QWidget()
        host.setStyleSheet("background:transparent;")
        host_lay = QVBoxLayout(host)
        host_lay.setContentsMargins(0, 0, 0, 0)
        host_lay.setSpacing(6)
        for kind in ("mRNA", "CDS", "protein"):
            card = ResultCard(kind)
            card.detailRequested.connect(self.detailRequested)
            card.manualChanged.connect(self.manualChanged)
            self.cards[kind] = card
            host_lay.addWidget(card)
        host_lay.addStretch(1)

        self._scroll = QScrollArea()
        self._scroll.setWidgetResizable(True)
        self._scroll.setFrameShape(QFrame.NoFrame)
        self._scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self._scroll.setStyleSheet("QScrollArea{background:transparent;}")
        self._scroll.setWidget(host)
        lay.addWidget(self._scroll, 1)
        self.setMinimumHeight(112)

    def set_results(self, results: Dict[str, CheckResult],
                    ref_ids: Dict[str, list], has_refs: bool) -> None:
        for kind, card in self.cards.items():
            res = results.get(kind)
            if kind == "protein" and res is None:
                card.setVisible(False)
                continue
            card.setVisible(True)
            card.set_result(res, ref_ids.get(kind, []))
        if has_refs:
            self.hint.setText("已按 ID 打分匹配；如需指定请用卡片下方下拉框")
        else:
            self.hint.setText("未加载参考序列文件，校验未运行")

    def set_blocking(self, blocking: bool) -> None:
        """校验结果为空时给出占位说明。"""
        for card in self.cards.values():
            card.btn_detail.setEnabled(not blocking and card._result is not None)
