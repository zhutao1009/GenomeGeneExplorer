"""界面主题。

只做一件事：把 Qt 默认外观调成"科研软件"该有的克制风格——浅色、细边框、
等宽序列字体、不要把注意力抢走。深浅两套配色由 ``is_dark`` 决定。
"""

from __future__ import annotations

from PySide6.QtGui import QColor, QFont, QFontDatabase, QPalette
from PySide6.QtWidgets import QApplication

UI_FONT_CANDIDATES = ["Microsoft YaHei UI", "Microsoft YaHei", "Segoe UI", "Noto Sans CJK SC"]
MONO_FONT_CANDIDATES = ["Consolas", "Cascadia Mono", "DejaVu Sans Mono", "Courier New"]

ACCENT = "#2563eb"
CARD_BG = "#ffffff"
BORDER = "#e5e7eb"
TEXT = "#111827"
MUTED = "#6b7280"


def pick_font(candidates) -> str:
    families = set(QFontDatabase.families())
    for name in candidates:
        if name in families:
            return name
    return candidates[-1]


def ui_font(size: int = 9) -> QFont:
    f = QFont(pick_font(UI_FONT_CANDIDATES))
    f.setPointSize(size)
    return f


def mono_font(size: int = 10) -> QFont:
    f = QFont(pick_font(MONO_FONT_CANDIDATES))
    f.setPointSize(size)
    f.setStyleHint(QFont.Monospace)
    f.setFixedPitch(True)
    f.setKerning(False)
    return f


STYLESHEET = """
QWidget { color: #111827; }
QMainWindow, QDialog { background: #f8fafc; }
QGroupBox {
    border: 1px solid #e5e7eb; border-radius: 6px; margin-top: 14px;
    background: #ffffff; padding-top: 6px;
}
QGroupBox::title {
    subcontrol-origin: margin; left: 10px; padding: 0 4px;
    color: #374151; font-weight: 600;
}
QPushButton {
    background: #ffffff; border: 1px solid #d1d5db; border-radius: 5px;
    padding: 5px 10px; color: #111827;
}
QPushButton:hover { background: #f3f4f6; border-color: #9ca3af; }
QPushButton:pressed { background: #e5e7eb; }
QPushButton:disabled { color: #9ca3af; background: #f9fafb; border-color: #e5e7eb; }
QPushButton#primary {
    background: #2563eb; border: 1px solid #1d4ed8; color: white; font-weight: 600;
}
QPushButton#primary:hover { background: #1d4ed8; }
QPushButton#primary:disabled { background: #bfdbfe; border-color: #bfdbfe; color: #eff6ff; }
QToolButton {
    background: transparent; border: 1px solid transparent; border-radius: 4px;
    padding: 3px 6px; color: #374151;
}
QToolButton:hover { background: #eef2ff; border-color: #c7d2fe; }
QToolButton:disabled { color: #c0c4cc; }
QLineEdit, QComboBox, QSpinBox, QPlainTextEdit, QTextEdit {
    background: #ffffff; border: 1px solid #d1d5db; border-radius: 5px;
    padding: 4px 6px; selection-background-color: #bfdbfe; selection-color: #111827;
}
QLineEdit:focus, QComboBox:focus, QSpinBox:focus { border-color: #2563eb; }
QComboBox::drop-down { border: none; width: 18px; }
QHeaderView::section {
    background: #f3f4f6; border: none; border-right: 1px solid #e5e7eb;
    border-bottom: 1px solid #e5e7eb; padding: 4px 6px; font-weight: 600; color: #374151;
}
QTableWidget, QTreeWidget, QListWidget {
    background: #ffffff; border: 1px solid #e5e7eb; border-radius: 6px;
    gridline-color: #f1f5f9; alternate-background-color: #fafafa;
}
QTableWidget::item:selected, QListWidget::item:selected {
    background: #dbeafe; color: #111827;
}
QTabWidget::pane { border: 1px solid #e5e7eb; border-radius: 6px; background: #ffffff; }
QTabBar::tab {
    background: #f3f4f6; border: 1px solid #e5e7eb; border-bottom: none;
    padding: 5px 12px; margin-right: 2px; border-top-left-radius: 5px;
    border-top-right-radius: 5px; color: #4b5563;
}
QTabBar::tab:selected { background: #ffffff; color: #111827; font-weight: 600; }
QScrollBar:vertical { background: #f8fafc; width: 12px; margin: 0; }
QScrollBar::handle:vertical { background: #cbd5e1; border-radius: 6px; min-height: 24px; }
QScrollBar::handle:vertical:hover { background: #94a3b8; }
QScrollBar:horizontal { background: #f8fafc; height: 12px; margin: 0; }
QScrollBar::handle:horizontal { background: #cbd5e1; border-radius: 6px; min-width: 24px; }
QScrollBar::handle:horizontal:hover { background: #94a3b8; }
QScrollBar::add-line, QScrollBar::sub-line { width: 0; height: 0; }
QSplitter::handle { background: #e5e7eb; }
QSplitter::handle:hover { background: #c7d2fe; }
QStatusBar { background: #ffffff; border-top: 1px solid #e5e7eb; }
QStatusBar QLabel { color: #4b5563; }
QMenuBar { background: #ffffff; border-bottom: 1px solid #e5e7eb; }
QMenuBar::item { padding: 5px 10px; background: transparent; }
QMenuBar::item:selected { background: #eef2ff; }
QMenu { background: #ffffff; border: 1px solid #e5e7eb; padding: 4px; }
QMenu::item { padding: 5px 22px 5px 20px; border-radius: 4px; }
QMenu::item:selected { background: #eef2ff; }
QProgressBar {
    border: 1px solid #e5e7eb; border-radius: 4px; background: #f3f4f6;
    text-align: center; height: 14px; color: #374151;
}
QProgressBar::chunk { background: #60a5fa; border-radius: 3px; }
QToolTip {
    background: #111827; color: #f9fafb; border: none; padding: 4px 6px;
    border-radius: 4px;
}
"""


def apply_theme(app: QApplication) -> None:
    app.setStyle("Fusion")
    app.setFont(ui_font(9))
    app.setStyleSheet(STYLESHEET)
    pal = QPalette()
    pal.setColor(QPalette.Window, QColor("#f8fafc"))
    pal.setColor(QPalette.Base, QColor("#ffffff"))
    pal.setColor(QPalette.Text, QColor(TEXT))
    pal.setColor(QPalette.WindowText, QColor(TEXT))
    pal.setColor(QPalette.Highlight, QColor("#bfdbfe"))
    pal.setColor(QPalette.HighlightedText, QColor(TEXT))
    app.setPalette(pal)
