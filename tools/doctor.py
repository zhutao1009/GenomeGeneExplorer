#!/usr/bin/env python
"""环境自检：告诉你当前用的是哪个 Python、依赖齐不齐、路径对不对。

用法：

    python tools/doctor.py            # 完整报告
    python tools/doctor.py --quiet    # 只有一行，供 start.bat 调用

退出码：0 = 依赖齐全可以启动；1 = 缺东西（并打印修复命令）。
"""

from __future__ import annotations

import argparse
import os
import platform
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

REQUIRED = [("PySide6", "PySide6"), ("Bio", "biopython"), ("pyfaidx", "pyfaidx")]
OPTIONAL = [("pytest", "pytest"), ("PyInstaller", "pyinstaller")]


def _version(mod_name: str) -> str:
    try:
        mod = __import__(mod_name)
    except Exception as exc:  # noqa: BLE001 - 缺依赖就是缺，不细分原因
        return f"缺失（{type(exc).__name__}）"
    for attr in ("__version__", "version", "VERSION"):
        val = getattr(mod, attr, None)
        if isinstance(val, str) and val:
            return val
    return "已安装"


def _settings_path() -> str:
    base = os.environ.get("APPDATA") or os.path.expanduser("~")
    return os.path.join(base, "GenomeGeneExplorer", "settings.json")


def collect() -> dict:
    info = {
        "exe": sys.executable,
        "version": platform.python_version(),
        "bits": 64 if sys.maxsize > 2**32 else 32,
        "impl": platform.python_implementation(),
        "frozen": bool(getattr(sys, "frozen", False)),
        "venv": sys.prefix != getattr(sys, "base_prefix", sys.prefix),
        "deps": [(label, _version(mod)) for mod, label in REQUIRED],
        "extra": [(label, _version(mod)) for mod, label in OPTIONAL],
    }
    info["missing"] = [label for label, ver in info["deps"] if "缺失" in ver]
    info["root"] = ROOT
    info["examples"] = os.path.isdir(os.path.join(ROOT, "examples"))
    info["exe_dist"] = os.path.join(
        ROOT, "dist", "GenomeGeneExplorer", "GenomeGeneExplorer.exe")
    info["settings"] = _settings_path()
    return info


def report(info: dict) -> str:
    lines = []
    lines.append("Genome Gene Explorer 环境自检")
    lines.append("=" * 46)
    lines.append(f"解释器    : {info['exe']}")
    lines.append(f"版本      : Python {info['version']} ({info['bits']} 位, {info['impl']})")
    lines.append(f"虚拟环境  : {'是' if info['venv'] else '否（直接用系统解释器）'}")
    lines.append(f"项目目录  : {info['root']}")
    lines.append("")
    lines.append("必需依赖：")
    for label, ver in info["deps"]:
        mark = "OK  " if "缺失" not in ver else "缺失"
        lines.append(f"  [{mark}] {label:<10} {ver}")
    lines.append("")
    lines.append("开发依赖（打包/测试才需要）：")
    for label, ver in info["extra"]:
        mark = "OK  " if "缺失" not in ver else "--  "
        lines.append(f"  [{mark}] {label:<12} {ver}")
    lines.append("")
    lines.append(f"示例数据  : examples/ {'存在' if info['examples'] else '缺失'}")
    lines.append(f"已打包程序: {'存在' if os.path.exists(info['exe_dist']) else '未构建'}  {info['exe_dist']}")
    lines.append(f"设置文件  : {info['settings']}")
    if info["missing"]:
        lines.append("")
        lines.append("修复方法（在项目目录执行）：")
        lines.append(f'  "{info["exe"]}" -m pip install -r requirements.txt')
        lines.append("若提示没有 pip，先执行：")
        lines.append(f'  "{info["exe"]}" -m ensurepip --upgrade')
        lines.append("国内网络较慢时可加镜像：")
        lines.append("  -i https://pypi.tuna.tsinghua.edu.cn/simple")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="检查运行环境")
    ap.add_argument("--quiet", action="store_true", help="只打印一行摘要")
    args = ap.parse_args(argv)

    info = collect()
    if args.quiet:
        if info["missing"]:
            print(f"[env] Python {info['version']} @ {info['exe']}  -> 缺少 {', '.join(info['missing'])}")
        else:
            parts = [f"Python {info['version']}"]
            parts += [f"{label} {ver}" for label, ver in info["deps"]]
            print("[env] " + " | ".join(parts))
        return 1 if info["missing"] else 0

    text = report(info)
    try:
        print(text)
    except UnicodeEncodeError:      # 老式控制台兜底
        sys.stdout.buffer.write(text.encode("utf-8", "replace"))
        sys.stdout.write("\n")
    return 1 if info["missing"] else 0


if __name__ == "__main__":
    sys.exit(main())
