#!/usr/bin/env python
"""启动脚本（Windows 下双击即用，也可 `python run_app.py`）。

如果启动失败，脚本会尝试在窗口里显示错误；控制台不可见时也不会"闪退无提示"。
"""

from __future__ import annotations

import os
import sys
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)


def _fatal(message: str, detail: str = "") -> None:
    """在无控制台的环境下也能看到错误。"""
    print(message)
    if detail:
        print(detail)
    try:
        from PySide6.QtWidgets import QApplication, QMessageBox

        app = QApplication.instance() or QApplication(sys.argv)
        box = QMessageBox()
        box.setIcon(QMessageBox.Critical)
        box.setWindowTitle("Genome Gene Explorer 启动失败")
        box.setText(message)
        if detail:
            box.setDetailedText(detail)
        box.exec()
    except Exception:
        pass


def main() -> int:
    try:
        from gge.ui.app import main as app_main
    except ImportError as exc:
        _fatal(
            "缺少依赖，无法启动。",
            "\n".join([
                f"当前解释器：{sys.executable}",
                f"Python 版本：{sys.version.split()[0]}",
                f"错误信息：{exc}",
                "",
                "推荐做法：在项目目录双击 start.bat（它会自动挑一个装好依赖的 Python）。",
                "",
                f"也可以手动安装到当前解释器：\n    \"{sys.executable}\" -m pip install -r requirements.txt",
                f"若提示没有 pip：\n    \"{sys.executable}\" -m ensurepip --upgrade",
                "",
                "查看环境详情：  python tools\\doctor.py",
                "或者直接用打包好的程序（不需要 Python）：  dist\\GenomeGeneExplorer\\GenomeGeneExplorer.exe",
            ]),
        )
        return 2
    try:
        return app_main(sys.argv[1:])
    except Exception:
        _fatal("程序运行时发生未处理的错误。", traceback.format_exc())
        return 1


if __name__ == "__main__":
    sys.exit(main())
