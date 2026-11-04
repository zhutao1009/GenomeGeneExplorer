"""打包与启动脚本的回归测试。

这里守的是一个真实踩过的坑：`start.bat` 被保存成 LF 换行后，
cmd.exe 解析 `if (...)` 块会错位，用户看到的是
`'cho' is not recognized as an internal or external command`
和 `The binding handle is invalid`——完全看不出真正原因。
所以「批处理必须是 CRLF + 纯 ASCII」值得用测试钉死。
"""

from __future__ import annotations

import pathlib
import re

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
BATS = sorted(ROOT.glob("*.bat"))


def test_batch_files_exist():
    names = {p.name for p in BATS}
    assert {"start.bat", "build_exe.bat", "_find_python.bat"} <= names, names


@pytest.mark.parametrize("path", BATS, ids=lambda p: p.name)
def test_batch_uses_crlf(path: pathlib.Path):
    data = path.read_bytes()
    lf = data.count(b"\n")
    crlf = data.count(b"\r\n")
    assert lf == crlf, f"{path.name} 有 {lf - crlf} 处裸 LF 换行，cmd.exe 会解析错位"


@pytest.mark.parametrize("path", BATS, ids=lambda p: p.name)
def test_batch_is_ascii_only(path: pathlib.Path):
    data = path.read_bytes()
    bad = [(i, b) for i, b in enumerate(data) if b > 127]
    assert not bad, f"{path.name} 含非 ASCII 字节（首处偏移 {bad[0][0] if bad else -1}）"


@pytest.mark.parametrize("path", BATS, ids=lambda p: p.name)
def test_batch_has_no_stray_chcp(path: pathlib.Path):
    """不再改控制台代码页：纯 ASCII 批处理不需要它，改了反而会让别的工具输出乱码。"""
    text = path.read_text(encoding="ascii")
    assert "chcp" not in text.lower()


def test_gitattributes_pins_batch_to_crlf():
    attr = ROOT / ".gitattributes"
    assert attr.exists(), "缺 .gitattributes，clone 后 *.bat 可能又变回 LF"
    text = attr.read_text(encoding="utf-8")
    assert re.search(r"^\*\.bat\s+text\s+eol=crlf", text, re.M), text


def test_start_bat_forwards_arguments():
    """start.bat 必须能把 --self-test 之类的参数透传给 run_app.py。"""
    text = (ROOT / "start.bat").read_text(encoding="ascii")
    assert "run_app.py" in text
    assert "%*" in text
    assert "_find_python.bat" in text


def test_build_script_runs_packaged_smoke_test():
    text = (ROOT / "build_exe.bat").read_text(encoding="ascii")
    assert "exe_smoke_test.py" in text, "打包脚本没有验证产物，可能产出跑不起来的目录"


def test_doctor_reports_dependencies():
    """doctor.py 是给用户看的诊断入口，依赖齐的时候必须返回 0。"""
    import importlib.util

    spec = importlib.util.spec_from_file_location("gge_doctor", ROOT / "tools" / "doctor.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    info = mod.collect()
    assert info["missing"] == [], info["deps"]
    assert mod.main(["--quiet"]) == 0
    assert "PySide6" in mod.report(info)


def test_exe_smoke_test_reports_missing_build(monkeypatch, tmp_path):
    """没打包时不应抛异常，而是明确告诉你先去 build_exe.bat。"""
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "gge_smoke", ROOT / "tools" / "exe_smoke_test.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    monkeypatch.setattr(mod, "EXE", str(tmp_path / "nope.exe"))
    assert mod.main([]) == 1
