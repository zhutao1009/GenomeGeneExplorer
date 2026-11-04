"""彩色序列导出与数据集打包导出的界面测试（离屏运行）。

只构造控件、点按钮的槽函数，不调用 ``exec()``（离屏下模态会卡住）；
真正"导出"的链路走 ``tasks`` 里的函数，与用户点确定后执行的路径一致。
"""

from __future__ import annotations

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PySide6.QtWidgets", reason="PySide6 未安装，跳过 GUI 测试")

from PySide6.QtCore import QEventLoop, QTimer  # noqa: E402
from PySide6.QtWidgets import QDialog, QDialogButtonBox, QMessageBox  # noqa: E402

from gge.config import Settings  # noqa: E402
from gge.extract import PANEL_CDS, PANEL_GENOMIC, PANEL_MRNA  # noqa: E402

from conftest import extract_one


@pytest.fixture(autouse=True)
def no_modal(monkeypatch):
    shown = []

    def _fake_exec(self):
        shown.append((self.windowTitle(), self.text()))
        return QMessageBox.Ok

    def _fake_warning(parent, title, text, *a, **kw):
        shown.append((title, text))
        return QMessageBox.Ok

    monkeypatch.setattr(QMessageBox, "exec", _fake_exec)
    monkeypatch.setattr(QMessageBox, "warning", _fake_warning)
    return shown


@pytest.fixture()
def extraction(provider, ann_gff3):
    return extract_one(provider, ann_gff3, "g1", "g1.1", up=400, dn=400)


def _wait(ms: int = 50) -> None:
    loop = QEventLoop()
    QTimer.singleShot(ms, loop.quit)
    loop.exec()


def _wait_idle(win, timeout: int = 30000) -> None:
    waited = 0
    while (win._job is not None or win._queue) and waited < timeout:
        _wait(50)
        waited += 50
    _wait(150)


# --------------------------------------------------------------------------- #
# 彩色序列文档对话框
# --------------------------------------------------------------------------- #
def test_doc_dialog_defaults(qapp, extraction, tmp_path):
    from gge.ui.dialogs import SequenceDocDialog

    dlg = SequenceDocDialog(extraction, [PANEL_MRNA], str(tmp_path / "a.docx"), 60, 8.5)
    vals = dlg.values()
    assert vals["path"].endswith(".docx")
    assert vals["panels"] == [PANEL_MRNA]
    assert vals["mode"] == "region"          # 默认与界面一致、文件最小
    assert vals["bpl"] == 60 and vals["include_header"] is True
    dlg.close()


def test_doc_dialog_format_switches_suffix(qapp, extraction, tmp_path):
    from gge.ui.dialogs import SequenceDocDialog

    dlg = SequenceDocDialog(extraction, [PANEL_MRNA], str(tmp_path / "a.docx"), 60, 8.5)
    dlg.cmb_format.setCurrentIndex(1)        # 网页
    assert dlg.path_edit.text().endswith(".html")
    dlg.cmb_format.setCurrentIndex(0)
    assert dlg.path_edit.text().endswith(".docx")
    dlg.close()


def test_doc_dialog_select_all_and_hint(qapp, extraction, tmp_path):
    from gge.ui.dialogs import SequenceDocDialog

    dlg = SequenceDocDialog(extraction, [PANEL_MRNA], str(tmp_path / "a.docx"), 60, 8.5)
    assert "bp" in dlg.hint.text()
    dlg._set_all(True)
    assert PANEL_GENOMIC in dlg.selected_panels()
    assert PANEL_CDS in dlg.selected_panels() or extraction.panels[PANEL_CDS].is_empty()
    dlg._set_all(False)
    assert dlg.selected_panels() == []
    assert "至少选择一个面板" in dlg.hint.text()
    dlg._on_ok()                              # 空选择必须拦住
    assert dlg.result() != QDialog.DialogCode.Accepted
    dlg.close()


def test_doc_dialog_warns_on_huge_selection(qapp, extraction, tmp_path):
    """超过 20 万 bp 要给出"会慢"的提示，而不是默默让人等。"""
    from gge.ui.dialogs import SequenceDocDialog

    dlg = SequenceDocDialog(extraction, [PANEL_MRNA], str(tmp_path / "a.docx"), 60, 8.5)
    import gge.ui.dialogs as dlg_mod

    old = dlg_mod.WARN_BASES
    dlg_mod.WARN_BASES = 10                   # 人为压低阈值
    try:
        dlg.cmb_mode.setCurrentIndex(dlg.cmb_mode.findData("both"))
        dlg._update_hint()
        assert "很长" in dlg.hint.text()
    finally:
        dlg_mod.WARN_BASES = old
    dlg.close()


# --------------------------------------------------------------------------- #
# 数据集打包对话框
# --------------------------------------------------------------------------- #
def test_dataset_dialog_lists_roles(qapp, demo_dir, tmp_path):
    from gge.ui.dialogs import DatasetExportDialog

    paths = {"genome": demo_dir["genome"], "annotation": demo_dir["gff3"],
             "mrna": demo_dir["mrna"]}
    dlg = DatasetExportDialog(paths, "演示", "", str(tmp_path / "d.zip"))
    vals = dlg.values()
    assert vals["name"] == "演示"
    assert vals["path"].endswith(".zip")
    assert set(vals["paths"]) == {"genome", "annotation", "mrna"}
    assert "原始" in dlg.hint.text()
    dlg.close()


def test_dataset_dialog_uses_byte_units(qapp, demo_dir, tmp_path):
    """文件大小必须用 B/KB，不能说成 bp——蛋白文件 136 字节不是 136 bp。"""
    from gge.ui.dialogs import DatasetExportDialog

    paths = {"genome": demo_dir["genome"], "annotation": demo_dir["gff3"],
             "protein": demo_dir["protein"]}
    dlg = DatasetExportDialog(paths, "演示", "", str(tmp_path / "d.zip"))
    texts = [dlg.table.item(dlg._rows[role], 3).text() for role in paths]
    assert not any("bp" in t for t in texts), texts
    assert all(t.endswith(("B", "KB", "MB", "GB")) for t in texts), texts
    dlg.close()


def test_dialogs_localize_buttons(qapp, extraction, demo_dir, tmp_path):
    """中文界面里不该冒出英文 Cancel。"""
    from gge.ui.dialogs import DatasetExportDialog, SequenceDocDialog

    doc = SequenceDocDialog(extraction, ["genomic"], str(tmp_path / "a.docx"))
    txt = [b.text() for b in doc.findChildren(QDialogButtonBox)[0].buttons()]
    assert "取消" in txt and "Cancel" not in txt, txt
    doc.close()

    paths = {"genome": demo_dir["genome"], "annotation": demo_dir["gff3"]}
    ds = DatasetExportDialog(paths, "x", "", str(tmp_path / "b.zip"))
    txt = [b.text() for b in ds.findChildren(QDialogButtonBox)[0].buttons()]
    assert "取消" in txt and "Cancel" not in txt, txt
    ds.close()


def test_dataset_dialog_blocks_without_annotation(qapp, demo_dir, tmp_path):
    from PySide6.QtCore import Qt

    from gge.ui.dialogs import DatasetExportDialog

    paths = {"genome": demo_dir["genome"], "annotation": demo_dir["gff3"]}
    dlg = DatasetExportDialog(paths, "x", "", str(tmp_path / "d.zip"))
    row = dlg._rows["annotation"]
    dlg.table.item(row, 0).setCheckState(Qt.Unchecked)
    assert "缺少必填文件" in dlg.hint.text()
    dlg._on_ok()
    assert dlg.result() != QDialog.DialogCode.Accepted
    dlg.close()


# --------------------------------------------------------------------------- #
# 主窗口接线
# --------------------------------------------------------------------------- #
@pytest.fixture()
def window(qapp, tmp_path):
    from gge.ui.main_window import MainWindow

    w = MainWindow(settings=Settings(path=str(tmp_path / "settings.json")))
    w.resize(1400, 900)
    w.show()
    yield w
    w.shutdown()
    w.close()


def test_window_exposes_role_paths_and_rich_copy(window):
    """分别导入 5 个文件后，能拿到路径打包；面板能复制成彩色 HTML。"""
    window.load_demo()
    _wait_idle(window)
    assert window.model is not None and window.provider is not None

    paths = window.current_role_paths()
    assert {"genome", "annotation"} <= set(paths)
    stats = window.current_stats()
    assert stats.get("n_genes", 0) > 0
    assert stats.get("total_bases", 0) > 0

    window.do_search("GeneOne")
    _wait_idle(window)
    assert window.search_panel.table.rowCount() > 0
    gene_id = window.search_panel.current_gene()
    assert gene_id
    window.select_gene(gene_id)
    _wait_idle(window)
    assert window.extraction is not None
    assert window.act_export_doc.isEnabled()
    assert window.act_export_dataset.isEnabled()

    window.copy_panel_rich(PANEL_MRNA)
    from PySide6.QtGui import QGuiApplication

    mime = QGuiApplication.clipboard().mimeData()
    assert mime.hasHtml() and "background-color" in mime.html()
    assert mime.hasText() and len(mime.text()) > 10


# --------------------------------------------------------------------------- #
# 回归：点「导出」必须真的落盘
# --------------------------------------------------------------------------- #
def test_doc_export_really_writes_file(window, tmp_path, monkeypatch):
    """回归 dlg.Accepted 事故：对话框 accepted 后任务根本没启动，界面却毫无提示。

    PySide6 里 ``dlg.Accepted`` 在实例上访问会抛 AttributeError（只有
    ``QDialog.Accepted`` / ``QDialog.DialogCode.Accepted`` 可用），异常发生在槽函数里，
    无控制台的 exe 上表现为"点了没反应"。这里把 exec 打桩成"用户点了导出"，
    断言文件真的生成、且内容是 docx。
    """
    from PySide6.QtWidgets import QDialog

    from gge.ui.dialogs import SequenceDocDialog

    out = tmp_path / "GeneOne_colored.docx"

    def fake_exec(self):
        self.path_edit.setText(str(out))
        return QDialog.DialogCode.Accepted

    monkeypatch.setattr(SequenceDocDialog, "exec", fake_exec)
    window.load_demo()
    _wait_idle(window)
    window.do_search("GeneOne")
    _wait_idle(window)
    assert window.extraction is not None

    window.export_sequence_document_dialog()
    _wait_idle(window)

    assert out.exists(), "对话框点了导出，但没有任何文件生成"
    assert out.stat().st_size > 1000
    assert out.read_bytes()[:2] == b"PK"


def test_dataset_export_really_writes_file(window, tmp_path, monkeypatch):
    """同上，数据集打包那条路。"""
    from PySide6.QtWidgets import QDialog

    from gge.ui.dialogs import DatasetExportDialog

    out = tmp_path / "demo_dataset.zip"

    def fake_exec(self):
        self.path_edit.setText(str(out))
        return QDialog.DialogCode.Accepted

    monkeypatch.setattr(DatasetExportDialog, "exec", fake_exec)
    window.load_demo()
    _wait_idle(window)

    window.export_dataset_dialog()
    _wait_idle(window)

    assert out.exists(), "点了打包导出，但没有生成 zip"
    assert out.read_bytes()[:2] == b"PK"


def test_no_instance_level_dialog_code():
    """防退化：``某对话框实例.Accepted`` 一定炸，必须用 ``QDialog.DialogCode.Accepted``。"""
    import re
    from pathlib import Path

    src = Path(__file__).resolve().parent.parent / "gge" / "ui" / "main_window.py"
    text = src.read_text(encoding="utf-8")
    bad = re.findall(r"\bdlg\.Accepted\b", text)
    assert not bad, f"main_window.py 里还有 {len(bad)} 处 dlg.Accepted（实例访问会抛 AttributeError）"


def test_default_export_dir_is_usable(window, tmp_path):
    """默认导出目录必须是绝对、存在、可写的目录——空串会变成相对路径，
    exe 双击启动时相对的是安装目录，用户就会"导出成功但找不到文件"。"""
    d = window.default_export_dir()
    assert os.path.isabs(d), d
    assert os.path.isdir(d), d
    assert os.access(d, os.W_OK), d

    window.settings["export_dir"] = str(tmp_path)
    assert window.default_export_dir() == str(tmp_path)

    # 记录里的目录已被删掉时不能原样返回
    gone = tmp_path / "gone"
    window.settings["export_dir"] = str(gone)
    window.settings["last_dir"] = str(gone)
    assert window.default_export_dir() != str(gone)


def test_doc_dialog_coord_off_by_default(qapp, extraction, tmp_path):
    """坐标列默认关闭（贴到论文里更干净），但用户仍可以手动打开。"""
    from gge.ui.dialogs import SequenceDocDialog

    dlg = SequenceDocDialog(extraction, [PANEL_MRNA], str(tmp_path / "a.docx"), 60, 8.5)
    assert dlg.values()["include_coord"] is False
    dlg.chk_coord.setChecked(True)
    assert dlg.values()["include_coord"] is True
    dlg.close()

    dlg2 = SequenceDocDialog(extraction, [PANEL_MRNA], str(tmp_path / "b.docx"), 60, 8.5,
                             include_coord=True)
    assert dlg2.values()["include_coord"] is True, "上次的选择要能带回来"
    dlg2.close()


# --------------------------------------------------------------------------- #
# 重启后自动恢复（旧版只填路径不加载，用户以为坏了）
# --------------------------------------------------------------------------- #
def _wait_until(cond, timeout: int = 30000) -> bool:
    waited = 0
    while waited < timeout:
        if cond():
            return True
        _wait(50)
        waited += 50
    return bool(cond())


def test_reopen_loads_previous_files(tmp_path):
    """再打开时不仅路径要在，数据必须真的加载好、能直接搜索。"""
    from gge.ui.main_window import MainWindow

    cfg = str(tmp_path / "settings.json")
    w1 = MainWindow(settings=Settings(path=cfg))
    w1.resize(1200, 800)
    w1.show()
    w1.load_demo()
    _wait_idle(w1)
    assert w1.model is not None and w1.provider is not None
    w1.do_search("GeneOne")
    _wait_idle(w1)
    w1._save_settings()
    w1.shutdown()
    w1.close()

    w2 = MainWindow(settings=Settings(path=cfg))
    w2.resize(1200, 800)
    w2.show()
    try:
        ok = _wait_until(lambda: w2.provider is not None and w2.model is not None)
        assert ok, "重新打开后没有自动加载上次的文件"
        assert _wait_until(lambda: w2.mrna_refs is not None), "参考 mRNA 也该恢复"
        row = w2.file_panel.rows["genome"]
        assert row.status.text().startswith("✓"), row.status.text()
        # 数据可用：搜索能直接出结果，不需要用户按 Ctrl+R
        _wait_idle(w2)
        w2.do_search("GeneOne")
        _wait_idle(w2)
        assert w2.search_panel.table.rowCount() > 0
        assert w2.extraction is not None
    finally:
        w2.shutdown()
        w2.close()


def test_auto_restore_can_be_turned_off(tmp_path):
    """关掉开关后只恢复路径、不加载（给超大基因组/注释的用户留出口）。"""
    from gge.ui.main_window import MainWindow

    cfg = str(tmp_path / "settings.json")
    w1 = MainWindow(settings=Settings(path=cfg))
    w1.show()
    w1.load_demo()
    _wait_idle(w1)
    w1._save_settings()
    w1.shutdown()
    w1.close()

    s = Settings(path=cfg)
    s["auto_restore"] = False
    s.save()

    w2 = MainWindow(settings=Settings(path=cfg))
    w2.show()
    try:
        _wait(600)                      # 留足时间让定时器有机会触发
        assert w2.provider is None and w2.model is None
        assert w2.file_panel.rows["genome"].path_edit.text(), "路径仍然要显示出来"
        assert not w2.act_auto_restore.isChecked()
    finally:
        w2.shutdown()
        w2.close()


def test_auto_restore_toggle_persists(tmp_path):
    from gge.ui.main_window import MainWindow

    cfg = str(tmp_path / "settings.json")
    w = MainWindow(settings=Settings(path=cfg))
    w.show()
    try:
        assert w.act_auto_restore.isChecked()
        w.act_auto_restore.setChecked(False)
        assert Settings(path=cfg).get("auto_restore") is False
        w.act_auto_restore.setChecked(True)
        assert Settings(path=cfg).get("auto_restore") is True
    finally:
        w.shutdown()
        w.close()
