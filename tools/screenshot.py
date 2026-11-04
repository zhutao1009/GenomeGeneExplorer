"""开发辅助：在真实窗口下截屏，用于人工核对界面。

``python tools/screenshot.py out_dir``

之所以不用 ``QT_QPA_PLATFORM=offscreen``：离屏平台没有字体数据库
（``QFontDatabase.families()`` 返回空），文字会渲染成方块，
只能验证布局不能验证显示。这里用真实窗口 + ``QWidget.grab()``。
"""

from __future__ import annotations

import os
import sys
import tempfile

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

# 截图用的数据集库与设置文件都放临时目录：既不往用户真实的 %APPDATA% 库里塞演示数据，
# 也避免"上一次截图留下的 dataset_id 让这次启动去载入一个已经不存在的数据集"。
_SHOT_TMP = tempfile.mkdtemp(prefix="gge_shot_")
os.environ["GGE_LIBRARY_DIR"] = os.path.join(_SHOT_TMP, "library")


def main() -> int:
    from PySide6.QtCore import QEventLoop, QTimer
    from PySide6.QtWidgets import QApplication, QComboBox

    from gge.ui.app import create_app
    from gge.ui.main_window import MainWindow

    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if "--out" in sys.argv:
        i = sys.argv.index("--out")
        args = [sys.argv[i + 1]] if i + 1 < len(sys.argv) else []
    outdir = args[0] if args else os.path.join(HERE, "screenshots")
    os.makedirs(outdir, exist_ok=True)

    from gge.config import Settings

    app = create_app([sys.argv[0]])
    win = MainWindow(Settings(os.path.join(_SHOT_TMP, "settings.json")))
    win.resize(1680, 1020)
    win.show()

    def wait(ms: int) -> None:
        loop = QEventLoop()
        QTimer.singleShot(ms, loop.quit)
        loop.exec()

    def wait_idle(timeout: int = 30000) -> None:
        waited = 0
        while (win._job is not None or win._queue) and waited < timeout:
            wait(50)
            waited += 50
        wait(200)

    win.load_demo()
    wait_idle()
    win.do_search("GeneOne")
    wait_idle()
    wait(300)
    print("[shot] 主窗口就绪", flush=True)

    win.grab().save(os.path.join(outdir, "01_main_window.png"))

    # 只看结构图
    win.structure.grab().save(os.path.join(outdir, "02_structure.png"))

    # 单个面板放大截屏
    for key in ("genomic", "exon", "mrna", "cds", "utr"):
        panel = win.panels.get(key)
        if panel and panel.canvas.isVisible():
            panel.grab().save(os.path.join(outdir, f"03_panel_{key}.png"))

    # 校验卡片区域
    win.validation_panel.grab().save(os.path.join(outdir, "04_validation.png"))

    # 差异详情对话框（g3 有 1 处替换）
    win.search_panel.edit.setText("GeneThree")
    win.do_search("GeneThree")
    wait_idle()
    if win.validation.get("mRNA"):
        from gge.ui.dialogs import DiffDialog

        dlg = DiffDialog(win.validation["mRNA"], win.extraction, win)
        dlg.resize(880, 620)
        dlg.show()
        wait(300)
        dlg.grab().save(os.path.join(outdir, "05_diff_dialog.png"))
        dlg.close()

    # 配色对话框
    from gge.ui.dialogs import ColorSettingsDialog

    dlg2 = ColorSettingsDialog(win.palette_map, win)
    dlg2.show()
    wait(250)
    dlg2.grab().save(os.path.join(outdir, "06_colors.png"))
    dlg2.close()

    # 彩色序列文档导出对话框（用上一次搜索的基因，此时面板里确实有序列）
    if win.extraction is not None:
        from gge.ui.dialogs import SequenceDocDialog

        keys = [k for k in ("genomic", "exon", "intron", "mrna", "cds", "utr")
                if win.extraction.panels.get(k) is not None]
        dlg5 = SequenceDocDialog(
            win.extraction, keys,
            os.path.join(outdir, "GeneOne_colored.docx"),
            60, 8.5, win.palette_map, include_coord=False, parent=win)
        dlg5.show()
        wait(350)
        dlg5.grab().save(os.path.join(outdir, "09_doc_export.png"))
        dlg5.close()

    # 数据集打包导出对话框（用示例目录里那几个文件，模拟"逐个手动导入"后的状态）
    demo_paths = {}
    for role, fn in (("genome", "demo_genome.fa"), ("annotation", "demo_annotation.gff3"),
                     ("mrna", "demo_mrna.fa"), ("cds", "demo_cds.fa"),
                     ("protein", "demo_protein.fa")):
        p = os.path.join(HERE, "examples", fn)
        if os.path.exists(p):
            demo_paths[role] = p
    if demo_paths:
        from gge.ui.dialogs import DatasetExportDialog

        dlg6 = DatasetExportDialog(demo_paths, "演示数据集", "截图用",
                                   os.path.join(outdir, "demo_dataset.zip"), win)
        dlg6.show()
        wait(350)
        dlg6.grab().save(os.path.join(outdir, "10_dataset_export.png"))
        dlg6.close()

    # 压缩包导入对话框 + 数据集库（用自带的演示包，文件名故意无信息量）
    demo_zip = os.path.join(HERE, "examples", "demo_bundle.zip")
    if os.path.exists(demo_zip):
        from gge import bundle as bundle_mod
        from gge.library import Library
        from gge.tasks import import_bundle_from_archive
        from gge.ui.bundle_dialog import BundleImportDialog, LibraryDialog

        print("[shot] 生成导入对话框…", flush=True)
        scan = bundle_mod.scan_archive(demo_zip)
        dlg3 = BundleImportDialog(scan, win)
        dlg3.resize(1240, 500)
        dlg3.show()
        wait(400)
        dlg3.grab().save(os.path.join(outdir, "07_bundle_import.png"))
        dlg3.close()

        lib = Library(os.path.join(os.environ["GGE_LIBRARY_DIR"], "library.sqlite"))
        roles = {r: c.member for r, c in scan.by_role().items() if c}
        import_bundle_from_archive(scan, roles, "演示数据集", lib.path)
        dlg4 = LibraryDialog(lib, win)
        dlg4.resize(960, 420)
        dlg4.show()
        wait(400)
        dlg4.grab().save(os.path.join(outdir, "08_library.png"))
        dlg4.close()

    print("截图已保存到", outdir, flush=True)
    win._save_settings()
    win.shutdown()
    win.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
