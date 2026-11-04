"""压缩包导入相关界面的测试（离屏运行）。

只构造控件并检查状态，不调用 ``exec()``（离屏下模态会阻塞）。
真正"解包 + 加载"的链路走 ``import_bundle_from_archive`` + ``MainWindow._bundle_imported``，
与用户点按钮时执行的路径完全一致。
"""

from __future__ import annotations

import os
import sys

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PySide6.QtWidgets", reason="PySide6 未安装，跳过 GUI 测试")

from PySide6.QtCore import QEventLoop, QTimer  # noqa: E402
from PySide6.QtWidgets import QDialogButtonBox  # noqa: E402

from gge.bundle import scan_archive  # noqa: E402
from gge.config import Settings  # noqa: E402
from gge.tasks import import_bundle_from_archive  # noqa: E402


@pytest.fixture(autouse=True)
def no_modal(monkeypatch):
    """把 QMessageBox.exec 换成立即返回，避免离屏下真的弹窗把测试挂住。

    返回的列表里会记录所有弹过的框，需要时可以直接断言。
    """
    from PySide6.QtWidgets import QMessageBox

    shown = []

    def _fake_exec(self):
        shown.append((self.windowTitle(), self.text()))
        return QMessageBox.Ok

    monkeypatch.setattr(QMessageBox, "exec", _fake_exec)
    return shown


@pytest.fixture()
def window(qapp, tmp_path):
    from gge.ui.main_window import MainWindow

    w = MainWindow(settings=Settings(path=str(tmp_path / "settings.json")))
    w.resize(1400, 900)
    w.show()
    yield w
    # 先等后台线程退出再销毁窗口：QThread 在运行中被回收会直接崩溃
    w.shutdown()
    w.close()


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


def _ok_button(dlg):
    return dlg.buttons.button(QDialogButtonBox.Ok)


# --------------------------------------------------------------------------- #
# 导入对话框
# --------------------------------------------------------------------------- #
def test_dialog_shows_every_member_and_roles(qapp, make_zip, demo_dir):
    from gge.ui.bundle_dialog import BundleImportDialog

    zp = make_zip({
        "data/main.fa": demo_dir["genome"],
        "annotation.gff3": demo_dir["gff3"],
        "refs/transcripts.fasta": demo_dir["mrna"],
    })
    dlg = BundleImportDialog(scan_archive(zp))
    assert dlg.table.rowCount() == 3
    assert dlg.dataset_name() == "bundle"        # 去掉扩展名的包名
    roles = dlg.roles()
    assert roles["genome"] == "data/main.fa"
    assert roles["annotation"] == "annotation.gff3"
    assert roles["mrna"] == "refs/transcripts.fasta"
    assert _ok_button(dlg).isEnabled()
    # 每行都要有"依据"文字，不能是空白
    for r in range(dlg.table.rowCount()):
        assert dlg.table.item(r, 4).text().strip()
    dlg.close()


def test_dialog_blocks_import_when_required_role_missing(qapp, make_zip, demo_dir):
    from gge.ui.bundle_dialog import BundleImportDialog

    # 完整的两件套：默认可导入
    zp = make_zip({"g.fa": demo_dir["genome"], "a.gff3": demo_dir["gff3"]})
    dlg = BundleImportDialog(scan_archive(zp))
    assert _ok_button(dlg).isEnabled()

    combos = {c.property("member"): c for c in dlg._combos}
    combos["g.fa"].setCurrentIndex(combos["g.fa"].findData(""))   # 取消基因组
    assert not _ok_button(dlg).isEnabled()
    assert "还缺" in dlg.summary.text()

    combos["g.fa"].setCurrentIndex(combos["g.fa"].findData("genome"))
    assert _ok_button(dlg).isEnabled()
    assert "还缺" not in dlg.summary.text()
    dlg.close()


def test_dialog_cannot_pick_same_file_twice(qapp, make_zip, demo_dir):
    """只给一个文件时，无论怎么选都不可能同时满足基因组与注释。"""
    from gge.ui.bundle_dialog import BundleImportDialog

    zp = make_zip({"g.fa": demo_dir["genome"]})
    dlg = BundleImportDialog(scan_archive(zp))
    combo = dlg._combos[0]
    for role in ("annotation", "genome"):
        combo.setCurrentIndex(combo.findData(role))
        assert not _ok_button(dlg).isEnabled()
    assert "还缺" in dlg.summary.text()
    dlg.close()


def test_dialog_role_combo_is_exclusive(qapp, make_zip, demo_dir):
    """两个文件不能同时被指定成基因组，后选的会顶掉先选的。"""
    from gge.ui.bundle_dialog import BundleImportDialog

    zp = make_zip({
        "g.fa": demo_dir["genome"],
        "a.gff3": demo_dir["gff3"],
        "m.fa": demo_dir["mrna"],
    })
    dlg = BundleImportDialog(scan_archive(zp))
    rows = {c.property("member"): c for c in dlg._combos}
    rows["m.fa"].setCurrentIndex(rows["m.fa"].findData("genome"))
    roles = dlg.roles()
    assert roles["genome"] == "m.fa"
    assert "g.fa" not in roles.values()
    dlg.close()


def test_dialog_name_defaults_to_archive_stem(qapp, make_zip, demo_dir):
    from gge.ui.bundle_dialog import BundleImportDialog

    zp = make_zip({"g.fa": demo_dir["genome"], "a.gff3": demo_dir["gff3"]},
                  name="Danio_rerio_GRCz11.zip")
    dlg = BundleImportDialog(scan_archive(zp))
    assert dlg.dataset_name() == "Danio_rerio_GRCz11"
    dlg.name_edit.setText("我的数据")
    assert dlg.dataset_name() == "我的数据"
    dlg.close()


# --------------------------------------------------------------------------- #
# 数据集库对话框
# --------------------------------------------------------------------------- #
def test_library_dialog_lists_and_marks_broken(qapp, library, tmp_path):
    from gge.ui.bundle_dialog import LibraryDialog

    good = tmp_path / "g.fa"
    ann = tmp_path / "a.gff3"
    good.write_text(">chr1\nACGT\n", encoding="utf-8")
    ann.write_text("##gff-version 3\n", encoding="utf-8")
    ds = library.add("好数据", {"genome": str(good), "annotation": str(ann)})
    library.add("坏数据", {"genome": str(tmp_path / "gone.fa"),
                          "annotation": str(ann)})

    dlg = LibraryDialog(library)
    assert dlg.table.rowCount() == 2
    status = {dlg.table.item(r, 0).text(): dlg.table.item(r, 5).text()
              for r in range(2)}
    assert status["好数据"] == "就绪"
    assert "缺" in status["坏数据"]

    dlg._select_id(ds.id)
    assert dlg.selected_dataset().id == ds.id
    assert dlg.btn_load.isEnabled()
    dlg.close()


# --------------------------------------------------------------------------- #
# 主窗口：导入 → 加载 → 切换
# --------------------------------------------------------------------------- #
def test_main_window_imports_and_loads_dataset(window, make_zip, demo_dir, library):
    zp = make_zip({
        "main.fa": demo_dir["genome"],
        "ann.gff3": demo_dir["gff3"],
        "t.fa": demo_dir["mrna"],
        "c.fa": demo_dir["cds"],
        "p.faa": demo_dir["protein"],
    })
    scan = scan_archive(zp)
    roles = {r: c.member for r, c in scan.by_role().items() if c}
    result = import_bundle_from_archive(scan, roles, "演示", library.path)
    window.library = library
    window._bundle_imported(result)
    _wait_idle(window)

    assert window.provider is not None and window.model is not None
    assert window.mrna_refs is not None and window.cds_refs is not None
    assert window.protein_refs is not None
    assert window.current_dataset is not None and window._dataset_id == result.dataset.id
    assert "演示" in window.file_panel.dataset_label.text()
    for role in ("genome", "annotation", "mrna"):
        text = window.file_panel.rows[role].path_edit.text()
        assert text and os.path.exists(text)

    window.search_panel.edit.setText("GeneOne")
    window.do_search("GeneOne")
    _wait_idle(window)
    assert window.extraction is not None
    assert window.extraction.gene_id == "g1"
    assert window.validation.get("mRNA") is not None
    assert window.validation["mRNA"].status == "identical"


def test_switching_dataset_clears_stale_roles(window, library, demo_dir, no_modal):
    """从"有 mRNA 参考"的数据集切到"没有"的数据集，旧的参考必须被清掉。"""
    base = {"genome": demo_dir["genome"], "annotation": demo_dir["gff3"]}
    ds1 = library.add("带参考", dict(base, mrna=demo_dir["mrna"]))
    ds2 = library.add("无参考", dict(base))

    window.library = library
    window.load_dataset(ds1)
    _wait_idle(window)
    assert window.mrna_refs is not None
    assert not no_modal, f"不应弹窗：{no_modal}"

    window.load_dataset(ds2)
    _wait_idle(window)
    assert window.mrna_refs is None, "切数据集后不应该还留着上一套的参考序列"
    assert window.provider is not None


def test_reset_view_clears_previous_gene(window, provider, ann_gff3, demo_dir,
                                        library, no_modal):
    """换数据前必须把上一条基因的结果清干净，否则两套数据的结果会混在一起。"""
    from gge.tasks import extract_and_validate

    window.provider = provider
    window.model = ann_gff3
    bundle = extract_and_validate(provider, ann_gff3, "g1", "g1.1", 500, 500)
    window._store_and_show(window._cache_key(), bundle)
    assert window.extraction is not None
    assert window.panels["exon"].canvas.block_count() > 0

    ds = library.add("另一套", {"genome": demo_dir["genome"],
                              "annotation": demo_dir["gff3"]})

    window._apply_dataset(ds)
    assert window.extraction is None
    assert window.current_gene_id == ""
    assert all(p.canvas.block_count() == 0 for p in window.panels.values())
    assert window.search_panel.table.rowCount() == 0
    assert not no_modal, f"换数据不应弹窗：{no_modal}"
    window._cancel_job()


def test_drag_drop_archive_routes_to_import(window, tmp_path, monkeypatch):
    """拖一个 zip 进窗口，应该走导入流程而不是"无法识别的文件类型"。"""
    zp = tmp_path / "x.zip"
    zp.write_bytes(b"PK\x03\x04")
    called = {}
    monkeypatch.setattr(window, "import_archive_dialog",
                        lambda path="": called.setdefault("path", path))

    class _Url:
        def toLocalFile(self):
            return str(zp)

    class _Mime:
        def urls(self):
            return [_Url()]

    class _Event:
        def mimeData(self):
            return _Mime()

    window.dropEvent(_Event())
    assert called["path"] == str(zp)


def test_if_gff_member_is_gtf_dialog_says_gtf(qapp, make_zip, demo_dir):
    from gge.ui.bundle_dialog import BundleImportDialog

    zp = make_zip({"g.fa": demo_dir["genome"], "genes.gtf": demo_dir["gtf"]})
    dlg = BundleImportDialog(scan_archive(zp))
    formats = [dlg.table.item(r, 1).text() for r in range(dlg.table.rowCount())]
    assert "GTF" in formats
    dlg.close()
