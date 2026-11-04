"""基因结构示意图。

用 ``QPainter`` 直接画，导出 PNG/SVG 与屏幕显示走同一段绘制代码，
因此"屏幕上看到的"和"导出的图"不会不一致。

图形约定
--------
* 横轴 = 基因组坐标（正链从左到右递增；负链基因在左上角标注转录方向箭头）；
* 灰色粗线 = 内含子，线上的小箭头指向**转录方向**；
* 绿色块 = 外显子；橙色粗块 = CDS；紫色/青色块 = 5′/3′ UTR；
* 左右两侧淡色带 = 上下游侧翼，宽度与用户设置的 bp 数一致。
"""

from __future__ import annotations

import os
from typing import Dict, Optional, Tuple

from PySide6.QtCore import QPointF, QRectF, QSize, Qt
from PySide6.QtGui import (QColor, QFontMetrics, QImage, QPainter, QPen,
                           QPixmap, QPolygonF)
from PySide6.QtWidgets import QSizePolicy, QWidget

from ..colors import ensure_palette, UI_COLORS
from ..models import (REGION_CDS, REGION_DOWNSTREAM, REGION_EXON, REGION_INTRON,
                      REGION_UPSTREAM, REGION_UTR3, REGION_UTR5, GeneExtraction)
from .theme import mono_font, ui_font


def _q(hex_color: str, alpha: int = 255) -> QColor:
    c = QColor(hex_color)
    c.setAlpha(alpha)
    return c


class GeneStructureView(QWidget):
    """可缩放/平移的基因结构图。"""

    def __init__(self, palette: Dict[str, Dict[str, str]] | None = None,
                 parent=None) -> None:
        super().__init__(parent)
        self._palette = ensure_palette(palette)
        self._ui = dict(UI_COLORS)
        self._ex: Optional[GeneExtraction] = None
        self._scale = 0.02          # px / bp
        self._offset = 1            # 左边界基因组坐标
        self._drag_x = None
        self._drag_offset = 1

        self.setMinimumHeight(132)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
        self.setMouseTracking(True)
        self.setToolTip("滚轮缩放，拖动平移；双击恢复整段显示")

    # ------------------------------------------------------------------ #
    def set_extraction(self, ex: Optional[GeneExtraction]) -> None:
        self._ex = ex
        self.fit()
        self.update()

    def set_palette_map(self, palette) -> None:
        self._palette = ensure_palette(palette)
        self.update()

    def extraction(self) -> Optional[GeneExtraction]:
        return self._ex

    def span(self) -> Tuple[int, int]:
        if not self._ex:
            return 1, 1
        return self._ex.region_start, self._ex.region_end

    def fit(self) -> None:
        if not self._ex:
            return
        s, e = self.span()
        w = max(80, self.width() - 2 * self._margin())
        self._offset = s
        self._scale = w / max(1, (e - s + 1))
        self.update()

    def _margin(self) -> int:
        return 24

    # ------------------------------------------------------------------ #
    # 坐标转换
    # ------------------------------------------------------------------ #
    def _x(self, pos: float, width: Optional[float] = None) -> float:
        return self._margin() + (pos - self._offset) * self._scale

    def _bp_at(self, x: float) -> float:
        return self._offset + (x - self._margin()) / self._scale if self._scale else self._offset

    # ------------------------------------------------------------------ #
    # 交互
    # ------------------------------------------------------------------ #
    def wheelEvent(self, event) -> None:  # noqa: N802
        if not self._ex:
            return
        factor = 1.25 if event.angleDelta().y() > 0 else 1 / 1.25
        anchor_bp = self._bp_at(event.position().x())
        new_scale = max(1e-5, min(8.0, self._scale * factor))
        # 保持鼠标下的碱基不动
        self._offset = anchor_bp - (event.position().x() - self._margin()) / new_scale
        self._scale = new_scale
        self._clamp()
        self.update()

    def mousePressEvent(self, event) -> None:  # noqa: N802
        if event.button() == Qt.LeftButton:
            self._drag_x = event.position().x()
            self._drag_offset = self._offset
            self.setCursor(Qt.ClosedHandCursor)

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
        if self._drag_x is not None and self._scale:
            dx = event.position().x() - self._drag_x
            self._offset = self._drag_offset - dx / self._scale
            self._clamp()
            self.update()

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802
        self._drag_x = None
        self.setCursor(Qt.ArrowCursor)

    def mouseDoubleClickEvent(self, event) -> None:  # noqa: N802
        self.fit()

    def _clamp(self) -> None:
        if not self._ex:
            return
        s, e = self.span()
        visible = self.width() / max(self._scale, 1e-9)
        min_off = s - visible * 0.5
        max_off = max(min_off, e - visible * 0.5)
        self._offset = max(min_off, min(max_off, self._offset))

    # ------------------------------------------------------------------ #
    # 绘制
    # ------------------------------------------------------------------ #
    def paintEvent(self, event) -> None:  # noqa: N802
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing, True)
        p.fillRect(self.rect(), _q(self._ui.get("structure_bg", "#ffffff")))
        self._draw(p, self.width(), self.height())
        p.end()

    def _draw(self, p: QPainter, w: int, h: int) -> None:
        margin = self._margin()
        if not self._ex:
            p.setPen(_q("#9ca3af"))
            p.setFont(ui_font(9))
            p.drawText(QRectF(0, 0, w, h), Qt.AlignCenter,
                       "选择基因后在此显示结构图")
            return

        ex = self._ex
        pal = self._palette

        # ---- 标题 ----
        p.setPen(_q("#111827"))
        f = ui_font(10)
        f.setBold(True)
        p.setFont(f)
        title = f"{ex.gene_name}  ({ex.gene_id})   {ex.tx_id}"
        p.drawText(QPointF(margin, 18), title)
        p.setFont(mono_font(8))
        p.setPen(_q("#6b7280"))
        meta = (f"{ex.seqid}:{ex.gene_start:,}-{ex.gene_end:,} ({ex.strand})   "
                f"外显子 {ex.n_exons} · 内含子 {ex.n_introns} · "
                f"mRNA {len(ex.mrna):,} bp · CDS {len(ex.cds):,} bp   "
                f"窗口 {ex.region_start:,}-{ex.region_end:,}")
        p.drawText(QPointF(margin, 33), meta)
        p.setPen(_q("#9ca3af"))
        p.drawText(QRectF(w - margin - 150, 21, 150, 12), Qt.AlignRight,
                   f"{self._scale:.4g} px/bp")

        top = 46
        bottom = h - 42
        mid = (top + bottom) / 2
        if bottom - top < 40:
            mid = top + 20

        # ---- 上下游色带 ----
        if ex.upstream_actual:
            gi = self._panel_interval("upstream")
            if gi:
                self._band(p, gi, top, bottom, pal[REGION_UPSTREAM]["bg"], "上游")
        if ex.downstream_actual:
            gi = self._panel_interval("downstream")
            if gi:
                self._band(p, gi, top, bottom, pal[REGION_DOWNSTREAM]["bg"], "下游")

        # ---- 内含子线 ----
        line_pen = QPen(_q(self._ui.get("structure_line", "#374151")), max(1.5, min(6.0, self._scale * 30)))
        p.setPen(line_pen)
        gin = self._gene_interval()
        if gin:
            x1, x2 = self._x(gin[0]), self._x(gin[1])
            p.drawLine(QPointF(x1, mid), QPointF(x2, mid))
            self._draw_arrows(p, x1, x2, mid, ex.strand)

        # ---- 外显子 / CDS / UTR ----
        blocks = ex.panels["exon"].blocks if "exon" in ex.panels else []
        for block in blocks:
            for reg in block.regions:
                if reg.genomic_start is None or reg.genomic_end is None:
                    continue
                kind = reg.kind
                conf = pal.get(kind, pal[REGION_EXON])
                if kind == REGION_CDS:
                    hh = 13
                elif kind in (REGION_UTR5, REGION_UTR3):
                    hh = 6
                else:
                    hh = 9
                x1 = self._x(reg.genomic_start)
                x2 = self._x(reg.genomic_end + 1)
                if x2 - x1 < 1.0:
                    x2 = x1 + 1.0
                rect = QRectF(x1, mid - hh, max(1.0, x2 - x1), hh * 2)
                p.setBrush(_q(conf["bg"]))
                p.setPen(QPen(_q(conf["fg"]), 0.8))
                p.drawRect(rect)
            # 外显子编号
            gi = block.genomic_interval
            if gi:
                idx = block.meta.get("exon_index", "")
                xc = (self._x(gi[0]) + self._x(gi[1] + 1)) / 2
                if self._x(gi[1] + 1) - self._x(gi[0]) > 14:
                    p.setFont(mono_font(7))
                    p.setPen(_q("#4b5563"))
                    p.drawText(QRectF(xc - 20, mid - 30, 40, 12), Qt.AlignCenter, f"{idx}")

        # ---- 坐标轴 ----
        axis_y = max(mid + 34, bottom + 6)
        p.setPen(QPen(_q("#d1d5db"), 1))
        p.drawLine(QPointF(margin, axis_y), QPointF(w - margin, axis_y))
        s, e = self.span()
        ticks = self._ticks(s, e)
        p.setFont(mono_font(7))
        for val in ticks:
            x = self._x(val)
            if x < margin - 1 or x > w - margin + 1:
                continue
            p.setPen(QPen(_q("#9ca3af"), 1))
            p.drawLine(QPointF(x, axis_y - 3), QPointF(x, axis_y + 3))
            p.setPen(_q("#6b7280"))
            tx = min(max(x - 45, margin), w - margin - 90)
            p.drawText(QRectF(tx, axis_y + 4, 90, 12), Qt.AlignCenter, f"{val:,}")

        # ---- 图例 ----
        legend = [(REGION_UPSTREAM, "上游"), (REGION_EXON, "外显子"),
                  (REGION_CDS, "CDS"), (REGION_UTR5, "UTR"), (REGION_DOWNSTREAM, "下游")]
        lx = margin
        ly = axis_y + 18
        p.setFont(ui_font(8))
        for kind, label in legend:
            conf = pal.get(kind, pal[REGION_EXON])
            p.setBrush(_q(conf["bg"]))
            p.setPen(QPen(_q(conf["fg"]), 0.8))
            p.drawRect(QRectF(lx, ly - 8, 10, 8))
            p.setPen(_q("#4b5563"))
            p.drawText(QRectF(lx + 13, ly - 10, 46, 12), Qt.AlignLeft | Qt.AlignVCenter, label)
            lx += 62


    def _draw_arrows(self, p: QPainter, x1: float, x2: float, mid: float, strand: str) -> None:
        """在内含子线上画转录方向箭头。"""
        if x2 - x1 < 30:
            return
        step = 46.0
        size = 3.2
        p.setBrush(_q(self._ui.get("structure_arrow", "#1f2937")))
        p.setPen(Qt.NoPen)
        x = x1 + step / 2
        while x < x2 - 4:
            if strand == "-":
                poly = QPolygonF([QPointF(x, mid), QPointF(x + size * 2, mid - size * 1.6),
                                  QPointF(x + size * 2, mid + size * 1.6)])
            else:
                poly = QPolygonF([QPointF(x + size * 2, mid), QPointF(x, mid - size * 1.6),
                                  QPointF(x, mid + size * 1.6)])
            p.drawPolygon(poly)
            x += step

    def _band(self, p: QPainter, gi: Tuple[int, int], top: float, bottom: float,
              color: str, label: str) -> None:
        x1, x2 = self._x(gi[0]), self._x(gi[1] + 1)
        if x2 <= x1:
            return
        p.fillRect(QRectF(x1, top, x2 - x1, bottom - top), _q(color))
        if x2 - x1 > 40:
            p.setPen(_q("#64748b"))
            p.setFont(ui_font(8))
            p.drawText(QRectF(x1, top + 2, min(x2 - x1, 90), 12),
                       Qt.AlignLeft | Qt.AlignVCenter, label)

    def _panel_interval(self, key: str) -> Optional[Tuple[int, int]]:
        if not self._ex:
            return None
        panel = self._ex.panels.get(key)
        if not panel or not panel.blocks:
            return None
        return panel.blocks[0].genomic_interval

    def _gene_interval(self) -> Optional[Tuple[int, int]]:
        if not self._ex:
            return None
        return self._ex.gene_start, self._ex.gene_end

    def _ticks(self, start: int, end: int, target: int = 6) -> list:
        span = max(1, end - start)
        raw = span / target
        mag = 10 ** max(0, int(f"{raw:e}".split("e")[1]))
        for mult in (1, 2, 5, 10):
            step = mag * mult
            if step >= raw:
                break
        else:
            step = mag * 10
        step = max(1, int(step))
        first = (int(start) // step + 1) * step
        out = [int(start)] + list(range(first, int(end), step)) + [int(end)]
        return sorted(set(out))

    # ------------------------------------------------------------------ #
    # 导出
    # ------------------------------------------------------------------ #
    def render_image(self, width: int = 1600, height: int = 420,
                     scale: float = 2.0) -> QImage:
        img = QImage(int(width * scale), int(height * scale), QImage.Format_ARGB32_Premultiplied)
        img.setDevicePixelRatio(scale)
        img.fill(_q("#ffffff"))
        p = QPainter(img)
        p.setRenderHint(QPainter.Antialiasing, True)
        saved = (self._scale, self._offset)
        # 导出时先把可视范围适配到整段，再绘制
        self._offset, self._scale = self.span()[0], (width - 2 * self._margin()) / max(
            1, self.span()[1] - self.span()[0] + 1)
        self._draw(p, width, height)
        self._scale, self._offset = saved
        p.end()
        return img

    def export_png(self, path: str, width: int = 1600, height: int = 420) -> bool:
        img = self.render_image(width, height)
        return img.save(path, "PNG")

    def export_svg(self, path: str, width: int = 1600, height: int = 420) -> bool:
        try:
            from PySide6.QtSvg import QSvgGenerator
        except ImportError:  # pragma: no cover
            return False
        gen = QSvgGenerator()
        gen.setFileName(path)
        gen.setSize(QSize(width, height))
        gen.setViewBox(QRectF(0, 0, width, height))
        gen.setTitle(f"{self._ex.tx_id if self._ex else 'gene'} structure")
        gen.setDescription("Generated by Genome Gene Explorer")
        p = QPainter(gen)
        p.setRenderHint(QPainter.Antialiasing, True)
        saved = (self._scale, self._offset)
        self._offset, self._scale = self.span()[0], (width - 2 * self._margin()) / max(
            1, self.span()[1] - self.span()[0] + 1)
        self._draw(p, width, height)
        self._scale, self._offset = saved
        p.end()
        return os.path.exists(path)

    def export(self, path: str) -> bool:
        if path.lower().endswith(".svg"):
            return self.export_svg(path)
        return self.export_png(path)
