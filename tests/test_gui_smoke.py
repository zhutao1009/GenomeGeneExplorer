"""GUI 冒烟测试（离屏运行，不需要显示器）。

只验证"数据能正确到达控件并被画出来"，不依赖鼠标交互；
任何会弹出模态对话框的路径都在测试里绕开（离屏下模态会阻塞）。
"""

from __future__ import annotations

import os
import sys

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PySide6.QtWidgets", reason="PySide6 未安装，跳过 GUI 测试")

from PySide6.QtCore import QSize  # noqa: E402
from PySide6.QtGui import QImage  # noqa: E402

from gge.config import Settings  # noqa: E402
from gge.extract import PANEL_ORDER  # noqa: E402
from gge.tasks import extract_and_validate  # noqa: E402


@pytest.fixture()
def settings(tmp_path):
    return Settings(path=str(tmp_path / "settings.json"))


@pytest.fixture()
def window(qapp, settings):
    from gge.ui.main_window import MainWindow

    w = MainWindow(settings=settings)
    w.resize(1400, 900)
    w.show()
    yield w
    w.shutdown()      # 先等后台线程退出，否则 Qt 析构时会崩掉整个进程
    w.close()


def _load(window, provider, ann, gid="g1", tx_id="g1.1", up=500, dn=500,
          mrna_refs=None, cds_refs=None):
    bundle = extract_and_validate(provider, ann, gid, tx_id, up, dn,
                                  mrna_refs, cds_refs)
    window.provider = provider
    window.model = ann
    window.current_gene_id = gid
    window.current_tx_id = tx_id
    window.spin_up.setValue(up)
    window.spin_dn.setValue(dn)
    # 与真实主窗口一致：结果先入缓存再渲染
    window._store_and_show(window._cache_key(), bundle)
    return bundle


def test_window_builds(window):
    assert len(window.panels) == 8
    for key in PANEL_ORDER:
        assert key in window.panels
    assert window.validation_panel is not None
    assert window.structure is not None
    # 空状态下所有面板都在，只是没有内容
    assert all(p.canvas.block_count() == 0 for p in window.panels.values())


def test_panels_filled_after_extraction(window, provider, ann_gff3, mrna_refs, cds_refs):
    bundle = _load(window, provider, ann_gff3, mrna_refs=mrna_refs, cds_refs=cds_refs)
    ex = bundle.extraction
    assert window.extraction is ex
    assert window.panels["exon"].canvas.block_count() == 3
    assert window.panels["intron"].canvas.block_count() == 2
    assert window.panels["mrna"].canvas.all_text() == ex.mrna
    assert window.panels["cds"].canvas.all_text() == ex.cds
    assert window.panels["genomic"].canvas.block_count() == 1
    # UTR 面板应有 5′ 与 3′ 两段
    assert window.panels["utr"].canvas.block_count() == 2


def test_canvas_renders_pixels(window, provider, ann_gff3):
    _load(window, provider, ann_gff3)
    canvas = window.panels["genomic"].canvas
    canvas.resize(QSize(900, 400))
    img = QImage(canvas.viewport().size(), QImage.Format_ARGB32)
    img.fill(0xFFFFFFFF)
    canvas.viewport().render(img)
    assert not img.isNull()
    # 至少要画出一些非白像素（文字/色块）
    colored = 0
    for y in range(0, img.height(), 3):
        for x in range(0, img.width(), 3):
            if img.pixel(x, y) != 0xFFFFFFFF:
                colored += 1
    assert colored > 50, f"序列视图几乎是空白（非白像素 {colored}）"


def test_structure_renders_and_exports(window, provider, ann_gff3, tmp_path):
    _load(window, provider, ann_gff3)
    win = window.structure
    win.resize(QSize(1000, 260))
    img = QImage(win.size(), QImage.Format_ARGB32)
    img.fill(0xFFFFFFFF)
    win.render(img)
    assert not img.isNull()

    png = tmp_path / "structure.png"
    svg = tmp_path / "structure.svg"
    assert win.export_png(str(png))
    assert png.stat().st_size > 3000
    assert win.export_svg(str(svg))
    assert svg.stat().st_size > 1000
    text = svg.read_text(encoding="utf-8", errors="ignore")
    assert "<svg" in text


def test_panel_copy_and_fasta(window, provider, ann_gff3):
    _load(window, provider, ann_gff3)
    panel = window.panels["cds"]
    copied = panel.copy_plain()
    assert set(copied) <= set("ACGTN")
    assert copied == window.extraction.cds
    # 未选中时复制整面板
    assert panel.canvas.selected_text() == ""


def test_search_results_table(window, ann_gff3):
    from gge.search import search_genes

    window.model = ann_gff3
    hits = search_genes(ann_gff3, "Gene")
    window.search_panel.set_hits(hits)
    assert window.search_panel.table.rowCount() == len(hits)
    assert window.search_panel.current_gene() == hits[0].gene_id


def test_bpl_and_font_updates(window, provider, ann_gff3):
    _load(window, provider, ann_gff3)
    window._apply_bpl(30)
    assert window.panels["genomic"].canvas.bases_per_line() == 30
    window._apply_font_size(13)
    assert window.panels["genomic"].canvas.font_size() == 13


def test_find_across_panels(window, provider, ann_gff3):
    _load(window, provider, ann_gff3)
    window.find_edit.setText("ATG")
    window._do_find()
    assert "未找到" not in window.find_info.text()
    assert "面板" in window.find_info.text()


def test_panel_visibility_toggle(window, provider, ann_gff3):
    _load(window, provider, ann_gff3)
    window._set_panel_visible("intron", False)
    assert not window.panels["intron"].isVisible()
    window._set_panel_visible("intron", True)
    assert window.panels["intron"].isVisible()


def test_minus_strand_window(window, provider, ann_gff3):
    _load(window, provider, ann_gff3, "g2", "g2.1")
    assert window.extraction.strand == "-"
    first_exon = window.extraction.panels["exon"].blocks[0]
    assert first_exon.genomic_interval[0] > window.extraction.panels["exon"].blocks[1].genomic_interval[0]
    assert window.panels["mrna"].canvas.all_text().startswith("ATG")


def test_validation_panel_cards(window, provider, ann_gff3, mrna_refs, cds_refs):
    _load(window, provider, ann_gff3, mrna_refs=mrna_refs, cds_refs=cds_refs)
    window._refresh_validation()
    card = window.validation_panel.cards["mRNA"]
    assert card.badge.text() == "完全一致"
    assert card._result is not None
    assert card.btn_detail.isEnabled()


def test_extraction_cache_hit(window, provider, ann_gff3):
    """参数不变时第二次提取应命中缓存，不再启动后台任务。"""
    _load(window, provider, ann_gff3)
    assert window.extraction is not None
    assert len(window._cache) == 1

    key = window._cache_key()
    cached = window._cache.get(key)
    assert cached is not None and cached.extraction is window.extraction

    window.extract_current()          # 应走缓存分支：不排队、不起线程
    assert window._job is None
    assert len(window._queue) == 0
    assert "缓存" in window.status_msg.text()


def test_cache_key_changes_with_parameters(window, provider, ann_gff3):
    """改上下游长度 / 换转录本 / 手动指定参考，都必须让缓存键失效。"""
    _load(window, provider, ann_gff3)
    k0 = window._cache_key()

    window.spin_up.setValue(300)
    k1 = window._cache_key()
    assert k1 != k0
    assert window._cache.get(k1) is None      # 新参数没有旧结果可用

    window.spin_up.setValue(500)
    assert window._cache_key() == k0          # 改回来又能命中

    window._manual["mRNA"] = "g1.1_ref"
    assert window._cache_key() != k0


def test_cache_cleared_on_new_annotation(window, provider, ann_gff3, tmp_path):
    """换了注释文件就清空缓存，避免把旧模型的结果当新结果。"""
    _load(window, provider, ann_gff3)
    assert window._cache
    window._annotation_loaded(ann_gff3)
    assert not window._cache


def test_demo_autosearch_waits_for_loading(window):
    """示例数据是排队异步加载的：自动搜索必须等加载完再发。

    早期版本写死 ``QTimer.singleShot(1200, do_search)``，机器慢一点就会在注释
    还没解析完时搜索，弹出一个「还没有加载注释文件」的模态框并卡在那里。
    """
    win = window
    called = []
    win.do_search = lambda text="": called.append(text)   # 只记录，不真的搜

    win.model = None
    win._job = None
    win._demo_autosearch(tries=0)
    assert called == [], "没有注释时不应发起搜索"

    win._job = object()                 # 假装还有后台任务在跑
    try:
        win._demo_autosearch(tries=0)   # tries=0：不再排重试，直接返回
    finally:
        win._job = None
    assert called == [], "后台任务没结束就不该搜索"

    win.model = object()                # 已经加载完
    win._demo_autosearch(tries=0)
    assert called == ["GeneOne"]

    win.model = None        # 用完还原，否则 closeEvent 会去读 model.path
