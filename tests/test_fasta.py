"""基因组 FASTA 随机访问的测试（`gge/fasta.py`）。

为什么单独一个文件
------------------
其他测试的 ``provider`` 夹具一律用 ``open_fasta(..., prefer_memory=True)``，
于是 **pyfaidx 那条路径在测试里从未被执行过**。真实后果：pyfaidx 0.9 把
``Fasta.index`` 改名成 ``Fasta.records``，包装层读 ``.index`` 直接
``AttributeError``，而 ``open_fasta`` 用裸 ``except Exception: pass`` 把它吞了，
最后报出"文件过大，请安装 pyfaidx"——用户装的 pyfaidx 明明好好的。
一个 900 MB 的草鱼基因组打不开，就是这么来的。

这里专门锁住两件事：
1. 磁盘模式（pyfaidx）能真的读长度和序列，且不依赖某个版本的私有属性名；
2. 磁盘模式失败时，**真实原因**必须出现在报错里，不能被内存模式的"文件过大"顶掉。
"""

from __future__ import annotations

import os

import pytest


def _write(path, text):
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text)
    return str(path)


@pytest.fixture()
def small_fasta(tmp_path):
    """两条序列，长度故意不相等且都不是 3 的倍数。"""
    rec1 = "ACGTACGTAC" * 7 + "ACG"          # 73 bp
    rec2 = "TTTTGGGGCC" * 5 + "TT"           # 52 bp
    return _write(tmp_path / "mini.fa",
                  f">chr1 first sequence\n{rec1}\n>chr2 second sequence\n{rec2}\n"), rec1, rec2


def test_pyfaidx_provider_reads_lengths(small_fasta):
    """磁盘模式必须能读出每条序列的长度（回归：pyfaidx 0.9 删掉了 Fasta.index）。"""
    path, rec1, rec2 = small_fasta
    from gge.fasta import PyfaidxProvider

    prov = PyfaidxProvider(path)
    try:
        assert prov.seqids() == ["chr1", "chr2"]
        assert prov.lengths() == {"chr1": len(rec1), "chr2": len(rec2)}
        assert prov.length("chr1") == len(rec1)
        assert prov.length("nope") is None
    finally:
        prov.close()


def test_pyfaidx_provider_fetch_returns_real_bases(small_fasta):
    """取出来的必须是文件里真实的碱基，不是占位符长度。"""
    path, rec1, rec2 = small_fasta
    from gge.fasta import PyfaidxProvider

    prov = PyfaidxProvider(path)
    try:
        assert prov.fetch("chr1", 1, 10) == rec1[:10]
        assert prov.fetch("chr1", 5, 5) == rec1[4]
        # 超界裁剪
        assert prov.fetch("chr2", 1, 9999) == rec2
        # 负链取反向互补
        assert prov.fetch_oriented("chr1", 1, 10, "-") == _revcomp(rec1[:10])
    finally:
        prov.close()


def test_open_fasta_uses_disk_mode_for_normal_file(small_fasta):
    """默认路径就该走 pyfaidx，而不是悄悄退回内存模式。"""
    path, _, _ = small_fasta
    from gge.fasta import open_fasta

    prov = open_fasta(path)
    try:
        assert "pyfaidx" in prov.mode, prov.mode
    finally:
        prov.close()


def test_faidx_lengths_survive_index_rename(monkeypatch):
    """索引属性改名时也不能崩：``.index`` / ``.records`` / ``.faidx.index`` 都试。"""
    from gge import fasta as fasta_mod

    class Entry:
        def __init__(self, rlen):
            self.rlen = rlen

    class FakeFasta:
        """只有 records（模拟 pyfaidx 0.9+），没有 index。"""

        def __init__(self):
            self.records = {"a": Entry(11), "b": Entry(22)}

        def keys(self):
            raise AssertionError("有 records 就不该退回到逐条 len()")

    assert fasta_mod._faidx_entry_lengths(FakeFasta()) == {"a": 11, "b": 22}

    class OldFasta:
        """模拟 pyfaidx 0.6/0.7：只有 index。"""

        def __init__(self):
            self.index = {"x": Entry(7)}

    assert fasta_mod._faidx_entry_lengths(OldFasta()) == {"x": 7}

    class Bare:
        """两个属性都没有：退回 len(fa[name])，一个都不能少。"""

        def keys(self):
            return ["k"]

        def __getitem__(self, name):
            return "ABC"

    assert fasta_mod._faidx_entry_lengths(Bare()) == {"k": 3}


def test_open_fasta_reports_disk_failure_not_size_excuse(tmp_path, monkeypatch):
    """磁盘模式失败时要把真实原因报出来。

    真实事故：pyfaidx 0.9 不兼容 → AttributeError 被吞 → 只剩"文件过大"，
    用户被引导去"安装 pyfaidx"，而它装着且能用。这条测试就是防这个。
    """
    from gge import fasta as fasta_mod
    from gge.errors import FileFormatError

    path = _write(tmp_path / "big.fa", ">chr1\nACGT\n")

    def boom(path, progress_cb=None):
        # 注意：monkeypatch 到模块属性上的是**普通函数**，不会像类方法那样自动绑定，
        # 所以第一个参数就是 path，不能写成 (self, path)。
        raise AttributeError("'Fasta' object has no attribute 'index'")

    def too_big(path, max_bytes=64 * 1024 * 1024):
        raise FileFormatError("文件过大（904.4 MB），无法使用内存模式。")

    monkeypatch.setattr(fasta_mod, "PyfaidxProvider", boom)
    monkeypatch.setattr(fasta_mod, "PlainFastaProvider", too_big)

    with pytest.raises(FileFormatError) as ei:
        fasta_mod.open_fasta(path)
    msg = str(ei.value)
    assert "磁盘模式" in msg and "内存模式" in msg
    assert "AttributeError" in msg, "必须带上真实异常类型，否则又变成误导"
    assert "'Fasta' object has no attribute 'index'" in msg


def test_open_fasta_missing_file_message(tmp_path):
    from gge.errors import FileFormatError
    from gge.fasta import open_fasta

    with pytest.raises(FileFormatError) as ei:
        open_fasta(str(tmp_path / "nope.fa"))
    assert "不存在" in str(ei.value)


def test_open_fasta_empty_file_message(tmp_path):
    from gge.errors import FileFormatError
    from gge.fasta import open_fasta

    p = tmp_path / "empty.fa"
    p.write_text("", encoding="utf-8")
    with pytest.raises(FileFormatError) as ei:
        open_fasta(str(p))
    assert "大小为 0" in str(ei.value)


def _revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]
