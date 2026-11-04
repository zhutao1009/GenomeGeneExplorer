"""彩色序列渲染控件与序列面板。

为什么自己画而不用 QTextEdit
---------------------------
一个 2 Mb 的基因区间如果塞进富文本控件，Qt 会为每个字符建一个 text fragment，
直接卡死。这里用 ``QAbstractScrollArea`` 自绘：

* 只绘制可视行（虚拟化），与序列长度基本无关；
* 每个碱基的颜色来自 ``Block.regions`` 压缩成的 bytearray，绘制时按"同色连续段"
  批量填矩形；
* 坐标、选区、搜索高亮、横向/纵向滚动全部自己控制。

面板（``SequencePanel``）在其上加了标题、独立复制/导出按钮、颜色图例和折叠，
以满足"多个序列同时可见、各自独立操作"的要求。
"""

from __future__ import annotations

import math
import re
from typing import Dict, List, Optional, Sequence, Tuple

from PySide6.QtCore import QPointF, QRectF, QSize, Qt, Signal
from PySide6.QtGui import (QColor, QFontMetrics, QGuiApplication, QKeySequence,
                           QPainter, QPen, QShortcut)
from PySide6.QtWidgets import (QAbstractScrollArea, QFrame, QHBoxLayout, QLabel,
                               QMenu, QSizePolicy, QToolButton, QVBoxLayout,
                               QWidget)

from ..colors import LEGEND_ORDER, ensure_palette, region_label
from ..models import Block, Panel
from .theme import mono_font, ui_font


def _qc(hex_color: str, alpha: int = 255) -> QColor:
    c = QColor(hex_color)
    if alpha < 255:
        c.setAlpha(alpha)
    return c


class SequenceCanvas(QAbstractScrollArea):
    """虚拟化彩色序列视图。"""

    cursorInfo = Signal(str)
    selectionChanged = Signal()
    zoomRequested = Signal(int)
    richCopyRequested = Signal()      # 复制为彩色（HTML）序列，可直接粘到 Word
    docExportRequested = Signal()     # 导出彩色序列文档

    def __init__(self, palette: Dict[str, Dict[str, str]] | None = None,
                 ui_colors: Dict[str, str] | None = None, parent=None) -> None:
        super().__init__(parent)
        self._palette = ensure_palette(palette)
        self._ui = dict(ui_colors or {})
        self._blocks: List[Block] = []
        self._codes: List[bytearray] = []
        self._kind_names: List[str] = []
        self._kind_ids: Dict[str, int] = {}
        self._bpl = 60
        self._pt = 10
        self._font = mono_font(self._pt)
        self._fm = QFontMetrics(self._font)
        self._cw = max(5, self._fm.horizontalAdvance("A"))
        self._gap = max(3, int(round(self._cw * 0.7)))
        self._line_h = self._fm.height() + 2
        self._gutter = 82
        self._layout: List[Tuple[int, int, int, int]] = []
        self._total_h = 0
        self._content_w = 0
        self._sel: Optional[Tuple[int, int, int]] = None
        self._anchor: Optional[Tuple[int, int]] = None
        self._dragging = False
        self._hits: List[Tuple[int, int, int]] = []
        self._hit_pos = -1
        self._query = ""
        self._auto_fit = False
        self._pending_goto: Optional[Tuple[int, int]] = None

        self.setFrameShape(QFrame.NoFrame)
        self.setFocusPolicy(Qt.StrongFocus)
        self.setMouseTracking(True)
        self.viewport().setCursor(Qt.IBeamCursor)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self._recalc()

    # ------------------------------------------------------------------ #
    # 基本设置
    # ------------------------------------------------------------------ #
    def set_palette_map(self, palette: Dict[str, Dict[str, str]],
                        ui_colors: Dict[str, str] | None = None) -> None:
        self._palette = ensure_palette(palette)
        if ui_colors:
            self._ui.update(ui_colors)
        self.viewport().update()

    def set_font_size(self, pt: int) -> None:
        self._pt = max(6, min(28, int(pt)))
        self._font = mono_font(self._pt)
        self._recalc()
        self._relayout()
        self.viewport().update()

    def font_size(self) -> int:
        return self._pt

    def set_bases_per_line(self, n: int) -> None:
        self._bpl = max(10, min(240, int(n)))
        self._recalc()
        self._relayout()
        self.viewport().update()

    def bases_per_line(self) -> int:
        return self._bpl

    def set_auto_fit(self, enabled: bool) -> None:
        """自适应：每行碱基数跟着面板宽度走，避免序列被右侧裁掉。"""
        self._auto_fit = bool(enabled)
        if self._auto_fit:
            self._fit_bpl()
            self.viewport().update()

    def _fit_bpl(self) -> bool:
        if not self._auto_fit:
            return False
        vw = self.viewport().width()
        if vw < 80:
            return False
        usable = vw - self._gutter - 14
        per_base = self._cw + self._gap / 10.0
        n = int(usable / per_base)
        n = max(20, (n // 10) * 10)
        n = min(240, n)
        if n != self._bpl:
            self._bpl = n
            self._recalc()
            self._relayout()
            return True
        return False

    def _recalc(self) -> None:
        self._fm = QFontMetrics(self._font)
        self._cw = max(5, self._fm.horizontalAdvance("A"))
        self._gap = max(3, int(round(self._cw * 0.7)))
        self._line_h = self._fm.height() + 2
        groups = math.ceil(self._bpl / 10)
        self._content_w = self._gutter + 8 + self._bpl * self._cw + max(0, groups - 1) * self._gap + 10

    # ------------------------------------------------------------------ #
    # 数据
    # ------------------------------------------------------------------ #
    def _kind_id(self, kind: str) -> int:
        if kind not in self._kind_ids:
            self._kind_ids[kind] = len(self._kind_names)
            self._kind_names.append(kind)
        return self._kind_ids[kind]

    def set_blocks(self, blocks: Sequence[Block]) -> None:
        self._blocks = list(blocks)
        self._codes = []
        for b in self._blocks:
            idx = bytearray(len(b.seq))
            for reg in b.regions:
                s = max(0, int(reg.start0))
                e = min(len(b.seq), int(reg.end0))
                if e <= s:
                    continue
                idx[s:e] = bytes([self._kind_id(reg.kind)]) * (e - s)
            self._codes.append(idx)
        self._sel = None
        self._hits = []
        self._hit_pos = -1
        self._query = ""
        self._relayout()
        self._scroll_to_interesting()
        self.viewport().update()

    def _scroll_to_interesting(self) -> None:
        """单块面板（基因组窗口）默认滚到第一个外显子/ CDS，而不是停在上游。"""
        if len(self._blocks) != 1:
            self.verticalScrollBar().setValue(0)
            return
        block = self._blocks[0]
        target = 0
        for reg in block.regions:
            if reg.kind in ("exon", "cds"):
                target = max(0, reg.start0 - self._bpl)
                break
        if target <= 0:
            self.verticalScrollBar().setValue(0)
            return
        self._pending_goto = (0, target)
        self.goto_block(0, target, align="top")

    def clear(self) -> None:
        self.set_blocks([])

    def block_count(self) -> int:
        return len(self._blocks)

    def all_text(self) -> str:
        return "".join(b.seq for b in self._blocks)

    def _relayout(self) -> None:
        self._layout = []
        y = 0
        for bi, b in enumerate(self._blocks):
            header_h = self._line_h + 8 if (b.title or b.note) else 0
            n_lines = max(1, math.ceil(len(b.seq) / self._bpl)) if b.seq else 1
            h = header_h + n_lines * self._line_h + 12
            self._layout.append((y, header_h, n_lines, h))
            y += h
        self._total_h = y
        self._update_scrollbars()

    def _update_scrollbars(self) -> None:
        vh = max(1, self.viewport().height())
        vw = max(1, self.viewport().width())
        vmax = max(0, self._total_h - vh)
        hmax = max(0, self._content_w - vw)
        self.verticalScrollBar().setRange(0, vmax)
        self.verticalScrollBar().setPageStep(vh)
        self.verticalScrollBar().setSingleStep(self._line_h)
        self.horizontalScrollBar().setRange(0, hmax)
        self.horizontalScrollBar().setPageStep(vw)
        self.horizontalScrollBar().setSingleStep(self._cw * 4)

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        self._fit_bpl()
        if self._pending_goto and self.viewport().height() > 40:
            bi, idx = self._pending_goto
            self._pending_goto = None
            self.goto_block(bi, idx, align="top")
        self._update_scrollbars()

    # ------------------------------------------------------------------ #
    # 绘制
    # ------------------------------------------------------------------ #
    def _region_color(self, bi: int, idx: int) -> Tuple[QColor, QColor]:
        kind = self._kind_names[self._codes[bi][idx]] if idx < len(self._codes[bi]) else "other"
        conf = self._palette.get(kind) or self._palette["other"]
        return _qc(conf["bg"]), _qc(conf["fg"])

    def _line_first_coord(self, b: Block, line: int) -> int:
        if b.coord_mode == "transcript":
            return b.coord_start + line * self._bpl
        if b.coord_reverse:
            return b.coord_start + len(b.seq) - 1 - line * self._bpl
        return b.coord_start + line * self._bpl

    def paintEvent(self, event) -> None:  # noqa: N802
        p = QPainter(self.viewport())
        p.setFont(self._font)
        p.fillRect(self.viewport().rect(), _qc(self._ui.get("canvas_bg", "#ffffff")))
        if not self._blocks:
            p.setPen(_qc(self._ui.get("gutter_fg", "#6b7280")))
            p.drawText(self.viewport().rect().adjusted(12, 12, -12, -12),
                       Qt.AlignTop | Qt.AlignLeft | Qt.TextWordWrap,
                       "（当前面板没有内容）")
            p.end()
            return

        voff = self.verticalScrollBar().value()
        hoff = self.horizontalScrollBar().value()
        vh = self.viewport().height()
        x0 = self._gutter + 8 - hoff

        # 左侧坐标槽背景
        p.fillRect(QRectF(0, 0, self._gutter, vh), _qc(self._ui.get("gutter_bg", "#f9fafb")))
        p.setPen(QPen(_qc("#e5e7eb"), 1))
        p.drawLine(QPointF(self._gutter - 0.5, 0), QPointF(self._gutter - 0.5, vh))

        base_pen_color = _qc(self._ui.get("text", "#111827"))
        muted = _qc(self._ui.get("gutter_fg", "#6b7280"))

        for bi, b in enumerate(self._blocks):
            y, header_h, n_lines, h = self._layout[bi]
            top = y - voff
            if top + h < 0 or top > vh:
                continue

            # 标题行
            if header_h:
                hy = top + 2
                p.setFont(self._font)
                p.setPen(_qc("#1f2937"))
                title_rect = QRectF(6, hy, self.viewport().width() - 12, header_h - 2)
                p.drawText(title_rect, Qt.AlignLeft | Qt.AlignVCenter,
                           f"▍{b.title}" if b.title else "")
                if b.note:
                    p.setFont(ui_font(8))
                    p.setPen(muted)
                    meta_w = 380
                    p.drawText(QRectF(max(120, self.viewport().width() - meta_w - 8), hy,
                                      meta_w, header_h - 2),
                               Qt.AlignRight | Qt.AlignVCenter, b.note)
                    p.setFont(self._font)

            first_line = max(0, (voff - y - header_h) // self._line_h)
            last_line = min(n_lines - 1, (voff + vh - y - header_h) // self._line_h)

            for line in range(first_line, last_line + 1):
                ly = top + header_h + line * self._line_h
                line_start = line * self._bpl
                line_end = min(len(b.seq), line_start + self._bpl)

                # 坐标
                coord = self._line_first_coord(b, line)
                p.setFont(self._font)
                p.setPen(muted)
                p.drawText(QRectF(4, ly, self._gutter - 10, self._line_h),
                           Qt.AlignRight | Qt.AlignVCenter, f"{coord:,}")

                if line_start >= line_end:
                    continue

                # 背景色块（按同色连续段批量填充）
                j = line_start
                codes = self._codes[bi]
                while j < line_end:
                    k = codes[j]
                    k2 = j + 1
                    while k2 < line_end and codes[k2] == k:
                        k2 += 1
                    bx = x0 + (j - line_start) * self._cw + ((j - line_start) // 10) * self._gap
                    ex = (x0 + (k2 - 1 - line_start) * self._cw
                          + ((k2 - 1 - line_start) // 10) * self._gap + self._cw)
                    conf = self._palette.get(self._kind_names[k]) or self._palette["other"]
                    p.fillRect(QRectF(bx, ly, ex - bx, self._line_h - 1),
                               _qc(conf["bg"]))
                    j = k2

                # 碱基文本（按 run 分段，段内再按 10 bp 分组绘制）
                j = line_start
                while j < line_end:
                    k = codes[j]
                    k2 = j + 1
                    while k2 < line_end and codes[k2] == k:
                        k2 += 1
                    conf = self._palette.get(self._kind_names[k]) or self._palette["other"]
                    p.setPen(_qc(conf["fg"]))
                    g = j
                    while g < k2:
                        group_end = min(k2, ((g - line_start) // 10 + 1) * 10 + line_start)
                        gx = x0 + (g - line_start) * self._cw + ((g - line_start) // 10) * self._gap
                        p.drawText(QPointF(gx, ly + self._line_h - 3), b.seq[g:group_end])
                        g = group_end
                    j = k2

                # 搜索高亮
                if self._hits:
                    for hi, (hb, hs, he) in enumerate(self._hits):
                        if hb != bi:
                            continue
                        if he <= line_start or hs >= line_end:
                            continue
                        s = max(hs, line_start)
                        e = min(he, line_end)
                        rx = x0 + (s - line_start) * self._cw + ((s - line_start) // 10) * self._gap
                        rw = (e - s) * self._cw + ((e - 1 - line_start) // 10 - (s - line_start) // 10) * self._gap
                        color = QColor("#f59e0b")
                        color.setAlpha(220 if hi == self._hit_pos else 120)
                        p.setPen(QPen(color, 2))
                        p.setBrush(Qt.NoBrush)
                        p.drawRect(QRectF(rx, ly, rw, self._line_h - 1))

                # 选区
                if self._sel and self._sel[0] == bi:
                    s, e = self._sel[1], self._sel[2]
                    if not (e <= line_start or s >= line_end):
                        s2 = max(s, line_start)
                        e2 = min(e, line_end)
                        sx = x0 + (s2 - line_start) * self._cw + ((s2 - line_start) // 10) * self._gap
                        sw = ((e2 - s2) * self._cw
                              + ((e2 - 1 - line_start) // 10 - (s2 - line_start) // 10) * self._gap)
                        p.fillRect(QRectF(sx, ly, sw, self._line_h - 1),
                                   QColor(37, 99, 235, 60))
            p.setPen(QPen(_qc("#f1f5f9"), 1))
            p.drawLine(QPointF(0, top + h - 1), QPointF(self.viewport().width(), top + h - 1))

        p.end()

    # ------------------------------------------------------------------ #
    # 命中测试与交互
    # ------------------------------------------------------------------ #
    def _pos_to_index(self, pos) -> Tuple[Optional[int], Optional[int]]:
        voff = self.verticalScrollBar().value()
        hoff = self.horizontalScrollBar().value()
        x = pos.x() + hoff
        y = pos.y() + voff
        if x < self._gutter:
            return None, None
        for bi, b in enumerate(self._blocks):
            top, header_h, n_lines, h = self._layout[bi]
            if not (top + header_h <= y <= top + h):
                continue
            line = int((y - top - header_h) // self._line_h)
            if line < 0 or line >= n_lines:
                return None, None
            rel_x = x - self._gutter - 8
            # 逐个碱基定位（考虑每 10 bp 的间隙）
            col = 0
            acc = 0.0
            bpl = self._bpl
            while col < bpl:
                w = self._cw
                if acc + w > rel_x:
                    break
                acc += w
                col += 1
                if col % 10 == 0 and col < bpl:
                    acc += self._gap
            idx = line * bpl + min(col, bpl - 1)
            if idx >= len(b.seq):
                idx = len(b.seq) - 1
            return bi, max(0, idx)
        return None, None

    def mousePressEvent(self, event) -> None:  # noqa: N802
        if event.button() != Qt.LeftButton:
            return
        bi, idx = self._pos_to_index(event.position())
        if bi is None:
            return
        self.setFocus(Qt.MouseFocusReason)
        self._anchor = (bi, idx)
        self._sel = (bi, idx, idx + 1)
        self._dragging = True
        self.selectionChanged.emit()
        self.viewport().update()

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
        bi, idx = self._pos_to_index(event.position())
        if bi is not None:
            b = self._blocks[bi]
            coord = b.coord_start + idx if not b.coord_reverse else b.coord_start + len(b.seq) - 1 - idx
            kind = ""
            if bi < len(self._codes) and idx < len(self._codes[bi]) and self._codes[bi]:
                kind = self._kind_names[self._codes[bi][idx]]
            label = region_label(kind, self._palette)
            self.cursorInfo.emit(f"{b.title} · 坐标 {coord:,} · 碱基 {b.seq[idx]} · {label}")
            if self._dragging and self._anchor and self._anchor[0] == bi:
                a = self._anchor[1]
                self._sel = (bi, min(a, idx), max(a, idx) + 1)
                self.selectionChanged.emit()
                self.viewport().update()

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802
        self._dragging = False

    def mouseDoubleClickEvent(self, event) -> None:  # noqa: N802
        bi, idx = self._pos_to_index(event.position())
        if bi is None:
            return
        self._sel = (bi, 0, len(self._blocks[bi].seq))
        self.selectionChanged.emit()
        self.viewport().update()

    def wheelEvent(self, event) -> None:  # noqa: N802
        if event.modifiers() & Qt.ControlModifier:
            delta = 1 if event.angleDelta().y() > 0 else -1
            self.zoomRequested.emit(delta)
            event.accept()
            return
        if event.modifiers() & Qt.ShiftModifier:
            bar = self.horizontalScrollBar()
            bar.setValue(bar.value() - event.angleDelta().y())
            event.accept()
            return
        super().wheelEvent(event)

    def keyPressEvent(self, event) -> None:  # noqa: N802
        if event.matches(QKeySequence.Copy):
            self.copy_selection()
            event.accept()
            return
        if event.key() == Qt.Key_A and event.modifiers() & Qt.ControlModifier:
            if self._blocks:
                self._sel = (0, 0, len(self._blocks[0].seq))
                self.selectionChanged.emit()
                self.viewport().update()
            event.accept()
            return
        super().keyPressEvent(event)

    def contextMenuEvent(self, event) -> None:  # noqa: N802
        menu = QMenu(self)
        act_copy = menu.addAction("复制选中序列")
        act_copy.setEnabled(bool(self._sel))
        act_all = menu.addAction("复制本面板全部序列")
        menu.addSeparator()
        act_rich = menu.addAction("复制为彩色序列（可直接粘贴到 Word）")
        act_doc = menu.addAction("导出彩色序列文档（Word / 网页）…")
        menu.addSeparator()
        act_find = menu.addAction("在此面板中查找…")
        chosen = menu.exec(event.globalPos())
        if chosen is act_copy:
            self.copy_selection()
        elif chosen is act_all:
            self.copy_to_clipboard(self.all_text())
        elif chosen is act_rich:
            self.richCopyRequested.emit()
        elif chosen is act_doc:
            self.docExportRequested.emit()
        elif chosen is act_find:
            self.zoomRequested.emit(0)  # 复用信号通知外层聚焦搜索框

    # ------------------------------------------------------------------ #
    # 复制 / 选区
    # ------------------------------------------------------------------ #
    def selected_text(self) -> str:
        if not self._sel:
            return ""
        bi, s, e = self._sel
        return self._blocks[bi].seq[s:e]

    def copy_selection(self) -> str:
        text = self.selected_text()
        if text:
            self.copy_to_clipboard(text)
        return text

    def copy_to_clipboard(self, text: str) -> None:
        if text:
            QGuiApplication.clipboard().setText(text)

    def selection_label(self) -> str:
        if not self._sel:
            return ""
        bi, s, e = self._sel
        return f"{self._blocks[bi].title}  选中 {e - s:,} bp"

    # ------------------------------------------------------------------ #
    # 查找与定位
    # ------------------------------------------------------------------ #
    def find_all(self, query: str, case_sensitive: bool = False) -> int:
        self._query = query or ""
        self._hits = []
        self._hit_pos = -1
        if not query:
            self.viewport().update()
            return 0
        q = query if case_sensitive else query.upper()
        total = 0
        for bi, b in enumerate(self._blocks):
            seq = b.seq if case_sensitive else b.seq.upper()
            start = 0
            while True:
                pos = seq.find(q, start)
                if pos < 0:
                    break
                self._hits.append((bi, pos, pos + len(q)))
                total += 1
                start = pos + max(1, len(q))
                if total > 5000:
                    break
            if total > 5000:
                break
        if self._hits:
            self._hit_pos = 0
            self._scroll_to_hit(0)
        self.viewport().update()
        return total

    def next_hit(self) -> int:
        if not self._hits:
            return -1
        self._hit_pos = (self._hit_pos + 1) % len(self._hits)
        self._scroll_to_hit(self._hit_pos)
        self.viewport().update()
        return self._hit_pos

    def prev_hit(self) -> int:
        if not self._hits:
            return -1
        self._hit_pos = (self._hit_pos - 1) % len(self._hits)
        self._scroll_to_hit(self._hit_pos)
        self.viewport().update()
        return self._hit_pos

    def clear_search(self) -> None:
        self._hits = []
        self._hit_pos = -1
        self._query = ""
        self.viewport().update()

    def _scroll_to_hit(self, i: int) -> None:
        bi, s, _e = self._hits[i]
        line = s // self._bpl
        y, header_h, _n, _h = self._layout[bi]
        target = y + header_h + line * self._line_h
        bar = self.verticalScrollBar()
        if target < bar.value() or target > bar.value() + self.viewport().height() - self._line_h * 2:
            bar.setValue(max(0, target - self.viewport().height() // 2))

    def goto_block(self, bi: int, idx: int = 0, align: str = "center") -> None:
        """滚动到第 ``bi`` 块的 ``idx`` 处。

        ``align="top"`` 让目标行贴着视口顶端（行首对齐，不会出现半行被切）；
        ``align="center"`` 尽量把目标行放在视口上三分之一处，同样按行高对齐。
        """
        if not (0 <= bi < len(self._blocks)):
            return
        line = idx // self._bpl
        y, header_h, _n, _h = self._layout[bi]
        base = y + header_h
        target = base + line * self._line_h
        if align == "top":
            value = target
        else:
            value = target - self.viewport().height() // 3
            if self._line_h > 0 and value > base:
                value = base + ((value - base) // self._line_h) * self._line_h
        self.verticalScrollBar().setValue(max(0, value))

    def goto_coordinate(self, coord: int) -> bool:
        """按基因组坐标定位（仅对 genomic 模式块有效）。"""
        for bi, b in enumerate(self._blocks):
            gi = b.genomic_interval
            if not gi or b.coord_mode != "genomic":
                continue
            lo, hi = min(gi), max(gi)
            if lo <= coord <= hi:
                offset = abs(coord - lo) if not b.coord_reverse else abs(hi - coord)
                self.goto_block(bi, offset)
                self._sel = (bi, offset, min(offset + 1, len(b.seq)))
                self.selectionChanged.emit()
                self.viewport().update()
                return True
        return False


class LegendBar(QWidget):
    """全局颜色图例（一条，放在序列区上方，避免每个面板重复占高度）。"""

    def __init__(self, palette: Dict[str, Dict[str, str]] | None = None, parent=None) -> None:
        super().__init__(parent)
        self._palette = ensure_palette(palette)
        self.setFixedHeight(24)
        self._lay = QHBoxLayout(self)
        self._lay.setContentsMargins(4, 0, 4, 0)
        self._lay.setSpacing(10)
        self._rebuild()

    def _rebuild(self) -> None:
        while self._lay.count():
            item = self._lay.takeAt(0)
            w = item.widget()
            if w:
                w.deleteLater()
        for kind in LEGEND_ORDER:
            conf = self._palette.get(kind, self._palette["other"])
            dot = QLabel()
            dot.setFixedSize(11, 11)
            dot.setStyleSheet(
                f"background:{conf['bg']}; border:1px solid #cbd5e1; border-radius:2px;")
            text = QLabel(region_label(kind, self._palette))
            text.setStyleSheet("color:#4b5563; font-size:11px;")
            box = QWidget()
            bl = QHBoxLayout(box)
            bl.setContentsMargins(0, 0, 0, 0)
            bl.setSpacing(4)
            bl.addWidget(dot)
            bl.addWidget(text)
            self._lay.addWidget(box)
        self._lay.addStretch(1)
        tip = QLabel("颜色优先级：CDS > UTR > 外显子 > 内含子 > 上下游")
        tip.setStyleSheet("color:#9ca3af; font-size:10px;")
        self._lay.addWidget(tip)

    def set_palette_map(self, palette) -> None:
        self._palette = ensure_palette(palette)
        self._rebuild()


class SequencePanel(QFrame):
    """带标题、独立操作按钮与图例的序列面板。"""

    requestExport = Signal(str)
    requestCopyFasta = Signal(str)
    requestRichCopy = Signal(str)
    requestDocExport = Signal(str)
    requestMaximize = Signal(str)
    cursorInfo = Signal(str)

    def __init__(self, key: str, title: str,
                 palette: Dict[str, Dict[str, str]] | None = None,
                 ui_colors: Dict[str, str] | None = None,
                 parent=None) -> None:
        super().__init__(parent)
        self.key = key
        self._panel: Optional[Panel] = None
        self._collapsed = False
        self._auto_fit = False
        self._palette = ensure_palette(palette)
        self._ui = dict(ui_colors or {})
        self.setObjectName("seqPanel")
        self.setFrameShape(QFrame.NoFrame)
        self.setStyleSheet(
            "QFrame#seqPanel { background:#ffffff; border:1px solid #e5e7eb; border-radius:6px; }"
        )
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        outer = QVBoxLayout(self)
        outer.setContentsMargins(1, 1, 1, 1)
        outer.setSpacing(0)

        header = QWidget()
        header.setStyleSheet("background:transparent;")
        hlay = QHBoxLayout(header)
        hlay.setContentsMargins(7, 2, 5, 2)
        hlay.setSpacing(3)
        self.title_label = QLabel(title)
        f = ui_font(9)
        f.setBold(True)
        self.title_label.setFont(f)
        self.title_label.setStyleSheet("color:#111827;")
        self.info_label = QLabel("")
        self.info_label.setStyleSheet("color:#6b7280; font-size:11px;")
        # 面板变窄时允许信息标签被挤压（而不是把整行顶宽）
        self.info_label.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Preferred)
        self.info_label.setMinimumWidth(0)
        hlay.addWidget(self.title_label)
        hlay.addWidget(self.info_label, 1)

        self.btn_copy = self._tool("复制", "复制选中序列（无选中时复制整面板）")
        self.btn_fasta = self._tool("FASTA", "复制为 FASTA 格式（含标题）")
        self.btn_export = self._tool("导出", "导出该面板为 FASTA 文件")
        self.btn_doc = self._tool("彩色文档导出", "导出带底色的序列文档（Word / 网页），"
                                                "颜色与界面一致")
        self.btn_max = self._tool("放大", "把该面板放大到整个序列区（再点一次还原）")
        self.btn_collapse = self._tool("折叠", "折叠/展开该面板")
        for b in (self.btn_copy, self.btn_fasta, self.btn_export, self.btn_doc,
                  self.btn_max, self.btn_collapse):
            hlay.addWidget(b)
        outer.addWidget(header)

        self.canvas = SequenceCanvas(self._palette, self._ui)
        outer.addWidget(self.canvas, 1)

        self.placeholder = QLabel("（无内容）")
        self.placeholder.setAlignment(Qt.AlignCenter)
        self.placeholder.setStyleSheet("color:#9ca3af; font-size:12px; padding:12px;")
        self.placeholder.hide()
        self.placeholder.hide()
        outer.addWidget(self.placeholder, 1)

        # 最小尺寸：保证 4 个面板 + 校验区能在 1080p 里同屏；再小就靠面板折叠/拖分隔条
        self.setMinimumHeight(74)
        self.setMinimumWidth(300)

        self.btn_copy.clicked.connect(self.copy_plain)
        self.btn_fasta.clicked.connect(lambda: self.requestCopyFasta.emit(self.key))
        self.btn_export.clicked.connect(lambda: self.requestExport.emit(self.key))
        self.btn_doc.clicked.connect(lambda: self.requestDocExport.emit(self.key))
        self.canvas.richCopyRequested.connect(lambda: self.requestRichCopy.emit(self.key))
        self.canvas.docExportRequested.connect(lambda: self.requestDocExport.emit(self.key))
        self.btn_max.clicked.connect(lambda: self.requestMaximize.emit(self.key))
        self.btn_collapse.clicked.connect(self.toggle_collapse)
        self.canvas.cursorInfo.connect(self.cursorInfo)
        self.canvas.zoomRequested.connect(self._on_zoom)

    def _tool(self, text: str, tip: str) -> QToolButton:
        b = QToolButton()
        b.setText(text)
        b.setToolTip(tip)
        b.setAutoRaise(True)
        b.setStyleSheet(
            "QToolButton { border:1px solid transparent; border-radius:4px; padding:2px 6px;"
            " color:#374151; font-size:11px; }"
            "QToolButton:hover { background:#eef2ff; border-color:#c7d2fe; }"
        )
        return b

    # ------------------------------------------------------------------ #
    def _on_zoom(self, delta: int) -> None:
        if delta == 0:
            return
        self.canvas.set_font_size(self.canvas.font_size() + delta)

    def set_panel(self, panel: Optional[Panel]) -> None:
        self._panel = panel
        if panel is None:
            self.clear()
            return
        self.title_label.setText(panel.title)
        # badge 放最左：标签栏变窄时被裁掉的是右侧的长度数字，而不是"注释无 UTR"这个关键提示
        info = f"{panel.badge} · " if panel.badge else ""
        info += f"{panel.length:,} bp · {len(panel.blocks)} 段" if panel.blocks else "空"
        self.info_label.setText(info)
        self.title_label.setToolTip(panel.note)
        self.info_label.setToolTip(panel.note)
        empty = panel.is_empty()
        self.canvas.setVisible(not empty)
        self.placeholder.setVisible(empty)
        if empty:
            self.placeholder.setText(panel.empty_text or "（无内容）")
        self.canvas.set_blocks(panel.blocks)
        self.canvas.set_auto_fit(self._auto_fit)

    def clear(self) -> None:
        self._panel = None
        self.canvas.clear()
        self.canvas.setVisible(True)
        self.placeholder.setText("（尚未提取）")
        self.placeholder.setVisible(True)
        self.info_label.setText("")
        self.title_label.setToolTip("")

    def toggle_collapse(self) -> None:
        self._collapsed = not self._collapsed
        for w in (self.canvas, self.placeholder):
            w.setVisible(not self._collapsed)
        self.btn_collapse.setText("展开" if self._collapsed else "折叠")
        if self._collapsed:
            self.setMaximumHeight(self.sizeHint().height())
        else:
            self.setMaximumHeight(16777215)

    def is_collapsed(self) -> bool:
        return self._collapsed

    # ------------------------------------------------------------------ #
    def copy_plain(self) -> str:
        text = self.canvas.selected_text() or self.canvas.all_text()
        self.canvas.copy_to_clipboard(text)
        return text

    def copy_all(self) -> str:
        text = self.canvas.all_text()
        self.canvas.copy_to_clipboard(text)
        return text

    def set_palette_map(self, palette, ui_colors=None) -> None:
        self._palette = ensure_palette(palette)
        self.canvas.set_palette_map(self._palette, ui_colors)

    def set_auto_fit(self, enabled: bool) -> None:
        self._auto_fit = bool(enabled)
        self.canvas.set_auto_fit(self._auto_fit)

    def set_font_size(self, pt: int) -> None:
        self.canvas.set_font_size(pt)

    def set_bases_per_line(self, n: int) -> None:
        self.canvas.set_bases_per_line(n)

    def find(self, text: str) -> int:
        return self.canvas.find_all(text)

    def next_hit(self) -> int:
        return self.canvas.next_hit()

    def clear_search(self) -> None:
        self.canvas.clear_search()
