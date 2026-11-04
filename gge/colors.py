"""配色方案。

默认颜色贴近需求文档的建议（上游浅蓝 / 外显子绿 / 内含子浅灰 / CDS 橙 /
5′UTR 紫 / 3′UTR 蓝绿 / 下游浅红），并把前景色一起定义，保证浅色背景上
黑字可读。所有颜色都可在界面中修改并持久化到配置文件。
"""

from __future__ import annotations

from typing import Dict

from .models import (REGION_CDS, REGION_DOWNSTREAM, REGION_EXON, REGION_FLANK,
                     REGION_INTRON, REGION_LABELS_CN, REGION_OTHER, REGION_UPSTREAM,
                     REGION_UTR3, REGION_UTR5)

DEFAULT_PALETTE: Dict[str, Dict[str, str]] = {
    REGION_UPSTREAM:   {"bg": "#dbeafe", "fg": "#1e3a8a", "label": "上游"},
    REGION_EXON:       {"bg": "#bbf7d0", "fg": "#14532d", "label": "外显子"},
    REGION_INTRON:     {"bg": "#e5e7eb", "fg": "#374151", "label": "内含子"},
    REGION_CDS:        {"bg": "#fed7aa", "fg": "#7c2d12", "label": "CDS"},
    REGION_UTR5:       {"bg": "#e9d5ff", "fg": "#4c1d95", "label": "5′ UTR"},
    REGION_UTR3:       {"bg": "#99f6e4", "fg": "#134e4a", "label": "3′ UTR"},
    REGION_DOWNSTREAM: {"bg": "#fecaca", "fg": "#7f1d1d", "label": "下游"},
    REGION_FLANK:      {"bg": "#f8fafc", "fg": "#475569", "label": "侧翼/基因间"},
    REGION_OTHER:      {"bg": "#ffffff", "fg": "#111827", "label": "其他"},
}

# 界面其它颜色
UI_COLORS = {
    "canvas_bg": "#ffffff",
    "gutter_bg": "#f9fafb",
    "gutter_fg": "#6b7280",
    "text": "#111827",
    "selection": "#3b82f6",
    "search_hit": "#f59e0b",
    "diff_same": "#ffffff",
    "diff_mismatch": "#fca5a5",
    "diff_gap": "#fde68a",
    "structure_bg": "#ffffff",
    "structure_line": "#374151",
    "structure_arrow": "#1f2937",
    "structure_tick": "#6b7280",
    "structure_cds": "#c2410c",
    "structure_exon": "#15803d",
    "structure_utr": "#7e22ce",
    "structure_up": "#bfdbfe",
    "structure_down": "#fecaca",
}

DIFF_STATUS_COLORS = {
    "identical": "#15803d",
    "different": "#b45309",
    "unmatched": "#b91c1c",
    "manual": "#7c3aed",
    "not_run": "#6b7280",
}

LEGEND_ORDER = [REGION_UPSTREAM, REGION_EXON, REGION_INTRON, REGION_CDS,
                REGION_UTR5, REGION_UTR3, REGION_DOWNSTREAM, REGION_FLANK]


def region_label(kind: str, palette: Dict[str, Dict[str, str]] | None = None) -> str:
    if palette and kind in palette and palette[kind].get("label"):
        return palette[kind]["label"]
    return REGION_LABELS_CN.get(kind, kind)


def ensure_palette(palette: Dict[str, Dict[str, str]] | None) -> Dict[str, Dict[str, str]]:
    """合并默认配色与用户配色，缺失的键用默认值补齐（兼容旧配置文件）。"""
    out = {k: dict(v) for k, v in DEFAULT_PALETTE.items()}
    if palette:
        for kind, conf in palette.items():
            if kind in out:
                out[kind].update({k2: v2 for k2, v2 in conf.items() if k2 in ("bg", "fg", "label")})
            else:
                out[kind] = dict(conf)
    return out
