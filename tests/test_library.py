"""数据集库（SQLite）的测试。"""

from __future__ import annotations

import os

import pytest

from gge.library import (Dataset, Library, ROLE_ORDER, default_db_path, library_dir,
                         store_root, slugify)


def _fake_files(tmp_path, roles=("genome", "annotation")):
    """造几个假文件当数据集内容。"""
    paths = {}
    for role in roles:
        p = tmp_path / f"{role}.dat"
        p.write_text(f"{role}\n", encoding="utf-8")
        paths[role] = str(p)
    return paths


def test_add_and_get(library, tmp_path):
    paths = _fake_files(tmp_path)
    ds = library.add("测试集", paths, source=str(tmp_path / "a.zip"),
                     fmt="gff3", stats={"n_genes": 12, "n_transcripts": 30,
                                        "n_seqids": 3, "total_bases": 1000})
    assert ds.id > 0
    assert ds.name == "测试集" and ds.fmt == "gff3"
    assert ds.n_genes == 12 and ds.n_transcripts == 30
    assert ds.genome == os.path.abspath(paths["genome"])
    assert ds.missing() == []
    again = library.get(ds.id)
    assert again is not None and again.name == ds.name
    assert len(again.files) == 2


def test_persists_across_instances(library, tmp_path):
    paths = _fake_files(tmp_path)
    ds = library.add("持久化", paths)
    other = Library(library.path)
    assert other.count() == 1
    assert other.get(ds.id).name == "持久化"


def test_list_orders_by_last_used(library, tmp_path):
    a = library.add("A", _fake_files(tmp_path, ("genome",)))
    b = library.add("B", _fake_files(tmp_path, ("annotation",)))
    names = [d.name for d in library.list()]
    assert names[0] == "B"           # 新建的排前面
    library.touch(a.id)
    assert [d.name for d in library.list()][0] == "A"


def test_search_filters_by_name_and_source(library, tmp_path):
    library.add("斑马鱼 GRCz11", _fake_files(tmp_path, ("genome",)),
                source=str(tmp_path / "zebra.zip"))
    library.add("小鼠 GRCm39", _fake_files(tmp_path, ("genome",)),
                source=str(tmp_path / "mouse.zip"))
    assert [d.name for d in library.list("斑马")] == ["斑马鱼 GRCz11"]
    assert [d.name for d in library.list("mouse")] == ["小鼠 GRCm39"]
    assert len(library.list()) == 2


def test_get_by_source(library, tmp_path):
    src = str(tmp_path / "pkg.zip")
    ds = library.add("X", _fake_files(tmp_path), source=src)
    assert library.get_by_source(src).id == ds.id
    assert library.get_by_source(str(tmp_path / "none.zip")) is None
    assert library.get_by_source("") is None


def test_rename_and_touch(library, tmp_path):
    ds = library.add("旧名", _fake_files(tmp_path))
    library.rename(ds.id, "新名")
    assert library.get(ds.id).name == "新名"
    library.rename(ds.id, "   ")          # 空名字应被忽略
    assert library.get(ds.id).name == "新名"
    library.touch(ds.id)
    assert library.get(ds.id).last_used


def test_remove_deletes_only_owned_files(library, tmp_path):
    """库目录内的文件随记录删除；用户自己指定到别处的文件必须留着。"""
    owned = library.new_dataset_dir("owned")
    inside = os.path.join(owned, "genome.fa")
    with open(inside, "w", encoding="utf-8") as fh:
        fh.write(">chr1\nACGT\n")
    outside = tmp_path / "my_own_genome.fa"
    outside.write_text(">chr1\nACGT\n", encoding="utf-8")

    ds1 = library.add("库内", {"genome": inside})
    ds2 = library.add("库外", {"genome": str(outside)})

    removed = library.remove(ds1.id, delete_files=True)
    assert removed and not os.path.exists(inside)
    assert library.get(ds1.id) is None

    removed2 = library.remove(ds2.id, delete_files=True)
    assert removed2 == []                 # 不在库目录里 → 不删
    assert outside.exists(), "绝不能删用户自己的文件"
    assert library.count() == 0


def test_remove_keeps_files_when_asked(library):
    owned = library.new_dataset_dir("keep")
    inside = os.path.join(owned, "a.fa")
    open(inside, "w").close()
    ds = library.add("保留文件", {"genome": inside})
    library.remove(ds.id, delete_files=False)
    assert os.path.exists(inside)


def test_missing_files_detected(library, tmp_path):
    paths = _fake_files(tmp_path)
    ds = library.add("会丢文件", paths)
    os.remove(paths["annotation"])
    reloaded = library.get(ds.id)
    assert "annotation" in reloaded.missing()
    assert reloaded.missing_required() == ["annotation"]


def test_files_rows_record_member_and_size(library, tmp_path):
    paths = _fake_files(tmp_path)
    ds = library.add("带成员名", paths, members={"genome": "data/main.fa"})
    f = ds.files["genome"]
    assert f.member == "data/main.fa"
    assert f.size == os.path.getsize(paths["genome"])


def test_update_paths_replaces_old_rows(library, tmp_path):
    ds = library.add("更新", _fake_files(tmp_path, ("genome",)))
    new_paths = _fake_files(tmp_path, ("genome", "cds"))
    library.update(ds.id, paths=new_paths, stats={"n_genes": 5})
    reloaded = library.get(ds.id)
    assert set(reloaded.paths()) == {"genome", "cds"}
    assert reloaded.n_genes == 5


def test_new_dataset_dir_is_unique(library):
    a = library.new_dataset_dir("same")
    b = library.new_dataset_dir("same")
    assert a != b and os.path.isdir(a) and os.path.isdir(b)
    assert os.path.dirname(a) == store_root()


def test_dataset_row_and_summary(library, tmp_path):
    ds = library.add("行显示", _fake_files(tmp_path),
                     source=str(tmp_path / "pkg.zip"),
                     stats={"n_genes": 3, "n_transcripts": 4})
    row = ds.row()
    assert row[0] == "行显示"
    assert "基因 3" in row[1]
    assert row[3] == "pkg.zip"
    assert ds.source_name() == "pkg.zip"


def test_slugify_keeps_ascii_and_falls_back():
    assert slugify("Zebrafish GRCz11") == "Zebrafish_GRCz11"
    assert slugify("斑马鱼") == "dataset"
    assert slugify("", fallback="x") == "x"
    assert len(slugify("a" * 100)) <= 40


def test_library_paths_follow_env(monkeypatch, tmp_path):
    monkeypatch.setenv("GGE_LIBRARY_DIR", str(tmp_path / "libroot"))
    monkeypatch.delenv("GGE_LIBRARY_DB", raising=False)
    assert library_dir() == str(tmp_path / "libroot")
    assert default_db_path() == os.path.join(str(tmp_path / "libroot"), "library.sqlite")
    assert store_root() == os.path.join(str(tmp_path / "libroot"), "datasets")


def test_role_order_is_stable():
    assert ROLE_ORDER == ("genome", "annotation", "mrna", "cds", "protein")


def test_empty_library_count(library):
    assert library.count() == 0
    assert library.list() == []
    assert library.get(999) is None
    assert library.remove(999, delete_files=True) == []
