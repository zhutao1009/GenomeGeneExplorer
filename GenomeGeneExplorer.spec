# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller 打包配置。

用法：
    python -m PyInstaller --noconfirm --clean GenomeGeneExplorer.spec
产物：
    dist/GenomeGeneExplorer/GenomeGeneExplorer.exe   （目录模式，启动快）
如需单文件 exe，把文件末尾的 COLLECT 段替换为注释里的 ONEFILE 段即可。
"""

import os
import sys

ROOT = os.path.abspath(os.getcwd())

# 示例数据一起打包，方便用户"开箱即用"验证功能
datas = []
ex_dir = os.path.join(ROOT, "examples")
if os.path.isdir(ex_dir):
    datas.append((ex_dir, "examples"))

# Qt 里用不到的模块排除掉，可显著减小体积
excludes = [
    "tkinter", "matplotlib", "scipy", "pandas", "notebook", "IPython",
    "PySide6.QtWebEngineCore", "PySide6.QtWebEngineWidgets", "PySide6.QtWebChannel",
    "PySide6.QtQuick", "PySide6.QtQuick3D", "PySide6.QtQml", "PySide6.QtQuickWidgets",
    "PySide6.Qt3DCore", "PySide6.Qt3DRender", "PySide6.Qt3DInput", "PySide6.Qt3DLogic",
    "PySide6.Qt3DAnimation", "PySide6.Qt3DExtras", "PySide6.QtCharts",
    "PySide6.QtDataVisualization", "PySide6.QtMultimedia", "PySide6.QtMultimediaWidgets",
    "PySide6.QtNetworkAuth", "PySide6.QtPdf", "PySide6.QtPdfWidgets",
    "PySide6.QtTest", "PySide6.QtDesigner", "PySide6.QtHelp", "PySide6.QtSql",
    "PySide6.QtBluetooth", "PySide6.QtNfc", "PySide6.QtPositioning",
    "PySide6.QtSerialPort", "PySide6.QtSensors", "PySide6.QtSerialBus",
    "PySide6.QtWebSockets", "PySide6.QtRemoteObjects", "PySide6.QtScxml",
    "PySide6.QtStateMachine", "PySide6.QtTextToSpeech", "PySide6.QtUiTools",
]

a = Analysis(
    ["run_app.py"],
    pathex=[ROOT],
    binaries=[],
    datas=datas,
    hiddenimports=["Bio", "Bio.Align", "pyfaidx", "PySide6.QtSvg"],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=excludes,
    noarchive=False,
    optimize=0,
)

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="GenomeGeneExplorer",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,          # 无控制台窗口；启动失败时会弹窗提示
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=None,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    upx_exclude=[],
    name="GenomeGeneExplorer",
)

# ---------------------------------------------------------------- #
# 单文件模式（把上面 COLLECT 换成下面这段）：
#
# exe = EXE(
#     pyz, a.scripts, a.binaries, a.datas, [],
#     name="GenomeGeneExplorer", debug=False, strip=False, upx=False,
#     console=False, onefile=True,
# )
# ---------------------------------------------------------------- #
