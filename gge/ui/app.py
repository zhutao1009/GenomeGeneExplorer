"""应用入口。

除了正常启动 GUI，还提供 ``--self-test``：在无显示器环境（offscreen）下
加载示例数据、跑一遍完整流程并把各面板渲染成图片，用于自动化验证。
"""

from __future__ import annotations

import argparse
import os
import sys
import traceback

from ..version import APP_NAME, APP_NAME_CN, ORG_NAME, APP_ID, __version__


def create_app(argv=None):
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import QApplication

    os.environ.setdefault("QT_ENABLE_HIGHDPI_SCALING", "1")
    app = QApplication(argv if argv is not None else sys.argv)
    app.setApplicationName(APP_NAME)
    app.setApplicationDisplayName(APP_NAME_CN)
    app.setOrganizationName(ORG_NAME)
    app.setApplicationVersion(__version__)
    from .theme import apply_theme

    apply_theme(app)
    return app


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=f"{APP_NAME} {__version__}")
    parser.add_argument("--genome", help="启动后立即加载的基因组 FASTA")
    parser.add_argument("--annotation", help="启动后立即加载的注释 GFF/GTF")
    parser.add_argument("--gene", help="启动后立即搜索/提取的基因")
    parser.add_argument("--import-bundle", metavar="ZIP",
                        help="启动后导入指定压缩包（自动识别文件）并加载")
    parser.add_argument("--list-datasets", action="store_true",
                        help="在命令行列出数据集库内容后退出（不需要图形界面）")
    parser.add_argument("--self-test", action="store_true",
                        help="离屏自检：加载示例数据、提取、渲染、导出，然后退出")
    parser.add_argument("--self-test-out", default="", help="自检输出目录")
    args, qt_args = parser.parse_known_args(argv)

    if args.list_datasets:
        return _list_datasets()

    if args.self_test:
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
        # 自检必须与用户环境隔离：它会把示例数据导进数据集库、并写 settings.json。
        # 早期版本直接用了真实的 MainWindow()，于是每次自检都在用户库里多出一套
        # 「演示数据集」，还把用户上次的真实文件路径覆盖成示例数据——用户再打开
        # 就看到一堆演示文件、自己的数据反而不见了。
        outdir = args.self_test_out or os.path.join(os.getcwd(), "selftest_out")
        os.environ["GGE_LIBRARY_DIR"] = os.path.join(outdir, "_library")
        os.makedirs(outdir, exist_ok=True)

    app = create_app([sys.argv[0]] + qt_args)
    from ..config import Settings
    from .main_window import MainWindow

    # 自检/命令行已经指定了要加载什么，就别再自动恢复上次的文件去抢队列
    explicit = bool(args.self_test or args.import_bundle or args.genome or args.annotation)
    if args.self_test:
        win = MainWindow(settings=Settings(
            path=os.path.join(args.self_test_out or
                              os.path.join(os.getcwd(), "selftest_out"),
                              "selftest_settings.json")),
            auto_restore=False)
    else:
        win = MainWindow(auto_restore=not explicit)

    if args.self_test:
        return _run_self_test(win, args, app)

    win.show()
    if args.import_bundle:
        win.import_archive_dialog(args.import_bundle)
    if args.genome:
        win.load_role("genome", args.genome)
    if args.annotation:
        win.load_role("annotation", args.annotation)
    if args.gene:
        win.search_panel.edit.setText(args.gene)
        win.do_search(args.gene)
    return app.exec()


def _file_size(n: int) -> str:
    """文件大小（不是碱基数，所以单位用 KB/MB）。"""
    if n >= 1024 ** 3:
        return f"{n / 1024 ** 3:.2f} GB"
    if n >= 1024 ** 2:
        return f"{n / 1024 ** 2:.1f} MB"
    if n >= 1024:
        return f"{n / 1024:.1f} KB"
    return f"{n} B"


def _pad(text: str, width: int) -> str:
    """按显示宽度补空格：中文占两列，直接 ljust 会错位。"""
    shown = sum(2 if ord(ch) > 0x2E80 else 1 for ch in text)
    return text + " " * max(0, width - shown)


def _list_datasets() -> int:
    """命令行列出数据集库，便于在服务器/脚本里核对导入结果。"""
    from ..library import Library, ROLE_LABELS

    lib = Library()
    rows = lib.list()
    print(f"数据集库：{lib.path}")
    print(f"共 {len(rows)} 套数据")
    for ds in rows:
        state = "就绪" if not ds.missing() else f"缺 {len(ds.missing())} 个文件"
        print(f"\n[{ds.id}] {ds.name}  （{state}）")
        print(f"     {' · '.join(x for x in (ds.summary(), ds.fmt.upper()) if x)}")
        if ds.source:
            print(f"     来源：{ds.source}")
        for role in ("genome", "annotation", "mrna", "cds", "protein"):
            p = ds.path(role)
            if p:
                ok = "✓" if os.path.exists(p) else "✕"
                size = _file_size(os.path.getsize(p)) if os.path.exists(p) else "-"
                print(f"     {ok} {_pad(ROLE_LABELS[role], 8)}{size:>9}  {p}")
    return 0


def _self_test_bundle(outdir: str) -> str:
    """优先用随程序发布的演示压缩包；没有就从示例文件现造一个。"""
    import zipfile

    from .main_window import _demo_dir

    base = _demo_dir()
    if base:
        packed = os.path.join(base, "demo_bundle.zip")
        if os.path.isfile(packed):
            return packed
    if not base:
        return ""
    target = os.path.join(outdir, "self_test_bundle.zip")
    pairs = [("data/assembly.fna", "demo_genome.fa"),
             ("data/annotation.gff3", "demo_annotation.gff3"),
             ("data/seq_a.fna", "demo_mrna.fa"),
             ("data/seq_b.fna", "demo_cds.fa"),
             ("data/pep.faa", "demo_protein.fa")]
    with zipfile.ZipFile(target, "w", zipfile.ZIP_DEFLATED) as zf:
        for member, src in pairs:
            path = os.path.join(base, src)
            if os.path.isfile(path):
                zf.write(path, member)
    return target if os.path.isfile(target) else ""


def _run_self_test(win, args, app) -> int:
    """离屏跑完整链路，返回 0 表示通过。"""
    from PySide6.QtCore import QEventLoop, QTimer
    from PySide6.QtGui import QImage, QPainter

    outdir = args.self_test_out or os.path.join(os.getcwd(), "selftest_out")
    os.makedirs(outdir, exist_ok=True)
    win.resize(1600, 1000)
    win.show()

    def wait(ms: int) -> None:
        loop = QEventLoop()
        QTimer.singleShot(ms, loop.quit)
        loop.exec()

    def wait_idle(timeout_ms: int = 20000) -> bool:
        """等到"确实没有任务在跑"，并且这个状态维持一个窗口期。

        不能只等一次：示例数据是 5 个任务排队加载，加载完还会触发
        ``_demo_autosearch`` → ``do_search`` → 提取，于是"空了一下的队列"马上又忙起来。
        旧写法把这种"紧接着又起了一个任务"误判成加载超时（自检会偶发失败，
        而总耗时只有 9 s，看着完全不像超时）。现在只要有新任务就把计时清零，
        只有连续 150 ms 都空闲才算完成。
        """
        waited = 0
        while True:
            busy = win._job is not None or win._queue
            if not busy:
                wait(150)
                if win._job is None and not win._queue:
                    return True
                waited = 0
                continue
            if waited >= timeout_ms:
                return False
            wait(50)
            waited += 50

    failures = []
    try:
        win.load_demo()
        if not wait_idle(60000):
            failures.append("示例数据加载超时")
        win.do_search("GeneOne")
        wait_idle()
        if win.extraction is None:
            failures.append("没有提取到 GeneOne")
        win.extraction and win.extraction.panels
        wait(100)

        # 各面板渲染成图片（真正调用 paintEvent）
        for key, panel in win.panels.items():
            canvas = panel.canvas
            if canvas.width() <= 0:
                continue
            img = QImage(canvas.viewport().width(), canvas.viewport().height(),
                         QImage.Format_ARGB32)
            img.fill(0xFFFFFFFF)
            canvas.viewport().render(img)
            img.save(os.path.join(outdir, f"panel_{key}.png"))

        struct = QImage(win.structure.width(), win.structure.height(), QImage.Format_ARGB32)
        struct.fill(0xFFFFFFFF)
        win.structure.render(struct)
        struct.save(os.path.join(outdir, "structure.png"))

        # 导出结构图与序列
        win.structure.export_png(os.path.join(outdir, "structure_export.png"))
        win.structure.export_svg(os.path.join(outdir, "structure_export.svg"))
        from ..export import gene_fasta_text, write_text

        write_text(os.path.join(outdir, "GeneOne.fasta"), gene_fasta_text(win.extraction))
        win._copy_gene_info()
        win.find_edit.setText("ATG")
        win._do_find()
        if "未找到" in win.find_info.text():
            failures.append("在所有面板中都找不到 ATG")

        if not win.validation:
            failures.append("校验结果为空（参考文件应已加载）")
        else:
            for kind, res in win.validation.items():
                if res.status not in ("identical", "different"):
                    failures.append(f"{kind} 校验状态异常：{res.status}")
                print(f"[self-test] {kind}: {res.status} {res.headline()}")
        if len(win.extraction.panels) != 8:
            failures.append(f"面板数量异常：{len(win.extraction.panels)}")

        # ---------------- 彩色序列文档（Word / 网页） ----------------
        from ..seqdoc import write_document

        doc_docx = os.path.join(outdir, "GeneOne_colored.docx")
        doc_html = os.path.join(outdir, "GeneOne_colored.html")
        write_document(doc_docx, win.extraction, ["genomic", "mrna", "cds"], mode="region",
                       include_coord=False)
        write_document(doc_html, win.extraction, ["genomic", "mrna", "cds"], mode="both")
        for path in (doc_docx, doc_html):
            if not os.path.isfile(path) or os.path.getsize(path) < 1000:
                failures.append(f"彩色文档导出失败：{path}")
        try:
            import zipfile

            with zipfile.ZipFile(doc_docx) as zf:
                xml = zf.read("word/document.xml").decode("utf-8")
            if "w:shd" not in xml:
                failures.append("彩色文档里没有字符底纹（w:shd），Word 会显示成黑白")
            if 'w:val="9CA3AF"' in xml:
                failures.append("显式关掉坐标列后，文档里仍出现坐标 run")
        except Exception:
            failures.append("彩色文档不是合法的 docx（zip 包不完整）")
        print("[self-test] color doc: .docx / .html 已导出（无坐标列）")

        # ---------------- 数据集压缩包（导出 → 再识别） ----------------
        from ..bundle import scan_archive as _scan_archive_once
        from ..bundle import write_dataset_archive

        role_paths = win.current_role_paths()
        ds_zip = os.path.join(outdir, "self_test_dataset.zip")
        try:
            exp = write_dataset_archive(ds_zip, role_paths, name="自检数据集",
                                        stats=win.current_stats())
            again = _scan_archive_once(ds_zip)
            back = {r: c.member for r, c in again.by_role().items() if c}
            if set(back) != set(role_paths):
                failures.append(
                    f"数据集压缩包再识别不一致：{sorted(back)} vs {sorted(role_paths)}")
            if again.manifest is None:
                failures.append("数据集压缩包里没有读回 gge_dataset.json 清单")
            print(f"[self-test] dataset: {len(exp.members)} 个文件 → "
                  f"{exp.archive_size / 1024:.0f} KB，清单已读回")
        except Exception as exc:
            failures.append(f"数据集压缩包导出失败：{exc}")

        # ---------------- 压缩包导入（自动识别）全链路 ----------------
        from ..bundle import scan_archive
        from ..library import Library
        from ..tasks import import_bundle_from_archive

        bundle_path = _self_test_bundle(outdir)
        if not bundle_path:
            failures.append("找不到/造不出演示压缩包")
        else:
            scan = scan_archive(bundle_path)
            roles = {r: c.member for r, c in scan.by_role().items() if c}
            print(f"[self-test] bundle: {scan.summary_text()}")
            for role in ("genome", "annotation", "mrna", "cds", "protein"):
                if role not in roles:
                    failures.append(f"压缩包识别失败：{role} 未被识别")
            lib = Library(os.path.join(outdir, "_library", "self_test.sqlite"))
            bres = import_bundle_from_archive(scan, roles, "自检数据集", lib.path)
            win.library = lib
            win._bundle_imported(bres)
            if not wait_idle(90000):
                failures.append("导入数据集后加载超时")
            if win.provider is None or win.model is None:
                failures.append("导入数据集后没有加载基因组/注释")
            else:
                win.search_panel.edit.setText("GeneOne")
                win.do_search("GeneOne")
                wait_idle()
                if win.extraction is None:
                    failures.append("导入数据集后无法提取 GeneOne")
                else:
                    res = win.validation.get("mRNA")
                    if res is None or res.status != "identical":
                        failures.append(
                            f"导入数据集的 mRNA 校验异常：{res.status if res else '无结果'}")
            if len(lib.list()) != 1:
                failures.append(f"数据集库记录数异常：{len(lib.list())}")
    except Exception:
        failures.append(traceback.format_exc())
    finally:
        win._save_settings()
        win.close()

    if failures:
        print("[self-test] 失败:")
        for f in failures:
            print("  -", f)
        return 1
    print(f"[self-test] 通过，输出位于 {outdir}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
