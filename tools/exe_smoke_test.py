#!/usr/bin/env python
"""打包产物冒烟测试：用 dist 里的 exe 跑一遍 --self-test，确认打包没漏东西。

用法：  python tools/exe_smoke_test.py [--keep]
退出码：0 = 通过；1 = 失败。
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
EXE = os.path.join(ROOT, "dist", "GenomeGeneExplorer", "GenomeGeneExplorer.exe")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="测试打包后的 exe")
    ap.add_argument("--keep", action="store_true", help="保留自检输出目录")
    ap.add_argument("--timeout", type=int, default=600)
    args = ap.parse_args(argv)

    if not os.path.exists(EXE):
        print(f"[FAIL] 没有找到打包产物：{EXE}")
        print("       请先执行 build_exe.bat")
        return 1

    out = tempfile.mkdtemp(prefix="gge_exetest_")
    env = dict(os.environ)
    env["QT_QPA_PLATFORM"] = "offscreen"
    env.pop("PYTHONPATH", None)          # 模拟干净机器：不借用源码目录
    env.pop("PYTHONHOME", None)

    size = os.path.getsize(EXE) / 2**20
    print(f"exe        : {EXE}  ({size:.1f} MB)")
    print(f"输出目录   : {out}")
    try:
        rc = subprocess.call([EXE, "--self-test", "--self-test-out", out],
                             env=env, timeout=args.timeout)
        files = sorted(os.listdir(out)) if os.path.isdir(out) else []
        pngs = [f for f in files if f.endswith(".png")]
        print(f"退出码     : {rc}")
        print(f"产出文件   : {len(files)} 个（其中 PNG {len(pngs)} 张）")
        for f in files[:12]:
            print(f"  - {f}")
        ok = rc == 0 and pngs
        print("[PASS] 打包产物可用" if ok else "[FAIL] 打包产物的自检没有跑通")
        return 0 if ok else 1
    except subprocess.TimeoutExpired:
        print(f"[FAIL] 超过 {args.timeout} 秒仍未结束")
        return 1
    finally:
        if args.keep:
            print(f"已保留：{out}")
        else:
            shutil.rmtree(out, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
