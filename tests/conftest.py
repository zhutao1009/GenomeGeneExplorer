"""pytest 公共夹具。

测试数据用 ``examples/make_example_data.py`` 现场生成到临时目录，
因此测试不依赖仓库里已存在的文件，也不会污染示例数据。
"""

from __future__ import annotations

import os
import sys
import zipfile

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXAMPLES = os.path.join(ROOT, "examples")
for p in (ROOT, EXAMPLES):
    if p not in sys.path:
        sys.path.insert(0, p)

import make_example_data  # noqa: E402


@pytest.fixture(scope="session", autouse=True)
def isolated_library(tmp_path_factory):
    """把数据集库重定向到临时目录。

    否则跑一次测试就会往用户 ``%APPDATA%`` 里塞一堆数据集，
    而且测试之间会互相干扰。
    """
    d = tmp_path_factory.mktemp("gge_library")
    old_dir = os.environ.get("GGE_LIBRARY_DIR")
    old_db = os.environ.get("GGE_LIBRARY_DB")
    os.environ["GGE_LIBRARY_DIR"] = str(d)
    os.environ.pop("GGE_LIBRARY_DB", None)
    yield str(d)
    if old_dir is None:
        os.environ.pop("GGE_LIBRARY_DIR", None)
    else:
        os.environ["GGE_LIBRARY_DIR"] = old_dir
    if old_db is not None:
        os.environ["GGE_LIBRARY_DB"] = old_db


@pytest.fixture(scope="session")
def demo_dir(tmp_path_factory):
    d = tmp_path_factory.mktemp("demo_data")
    return make_example_data.build_example_data(str(d))


@pytest.fixture(scope="session")
def demo_meta(demo_dir):
    return demo_dir


@pytest.fixture(scope="session")
def provider(demo_dir):
    from gge.fasta import open_fasta

    p = open_fasta(demo_dir["genome"], prefer_memory=True)
    yield p
    p.close()


@pytest.fixture(scope="session")
def ann_gff3(demo_dir):
    from gge.annotation import parse_annotation

    return parse_annotation(demo_dir["gff3"])


@pytest.fixture(scope="session")
def ann_gtf(demo_dir):
    from gge.annotation import parse_annotation

    return parse_annotation(demo_dir["gtf"])


@pytest.fixture(scope="session")
def mrna_refs(demo_dir):
    from gge.validate import load_reference

    return load_reference(demo_dir["mrna"], "mRNA")


@pytest.fixture(scope="session")
def cds_refs(demo_dir):
    from gge.validate import load_reference

    return load_reference(demo_dir["cds"], "CDS")


@pytest.fixture(scope="session")
def protein_refs(demo_dir):
    from gge.validate import load_reference

    return load_reference(demo_dir["protein"], "protein")


@pytest.fixture()
def write_annotation(tmp_path):
    """把注释文本写进临时文件，返回路径。"""

    def _write(text: str, name: str = "custom.gff3") -> str:
        p = tmp_path / name
        p.write_text(text, encoding="utf-8")
        return str(p)

    return _write


@pytest.fixture()
def write_fasta(tmp_path):
    def _write(records, name: str = "custom.fa") -> str:
        p = tmp_path / name
        with open(p, "w", encoding="utf-8", newline="\n") as fh:
            for header, seq in records:
                fh.write(f">{header}\n")
                for i in range(0, len(seq), 60):
                    fh.write(seq[i:i + 60] + "\n")
        return str(p)

    return _write


def extract_one(provider, ann, gene_id: str, tx_id: str = "", up=2000, dn=2000):
    from gge.extract import extract_gene

    gene = ann.get_gene(gene_id)
    assert gene is not None, f"注释中找不到基因 {gene_id}"
    tx = gene.transcripts.get(tx_id) if tx_id else gene.canonical_transcript()
    assert tx is not None
    return extract_gene(provider, gene, tx, up, dn)


@pytest.fixture()
def library(isolated_library):
    """每个测试一个全新的库（同一目录、不同文件名），互不影响。"""
    from gge.library import Library

    name = f"lib_{os.urandom(4).hex()}.sqlite"
    return Library(os.path.join(isolated_library, name))


@pytest.fixture()
def make_zip(tmp_path):
    """把 {包内路径: 源文件路径 或 ("text", 内容)} 打成一个 zip，返回路径。"""

    def _make(entries, name: str = "bundle.zip", compression=zipfile.ZIP_DEFLATED) -> str:
        path = tmp_path / name
        with zipfile.ZipFile(path, "w", compression) as zf:
            for member, src in entries.items():
                if isinstance(src, tuple) and src[0] == "text":
                    zf.writestr(member, src[1])
                else:
                    zf.write(src, member)
        return str(path)

    return _make


@pytest.fixture()
def make_tar_gz(tmp_path):
    import tarfile

    def _make(entries, name: str = "bundle.tar.gz") -> str:
        path = tmp_path / name
        with tarfile.open(path, "w:gz") as tf:
            for member, src in entries.items():
                tf.add(src, arcname=member)
        return str(path)

    return _make


@pytest.fixture(scope="session")
def qapp():
    """整个测试会话只允许存在一个 QApplication。

    每个测试模块各建一个会触发 "Please destroy the QApplication singleton"，
    因此统一定义在这里，两个 GUI 测试模块共用。
    """
    pytest.importorskip("PySide6.QtWidgets", reason="PySide6 未安装，跳过 GUI 测试")
    from gge.ui.app import create_app

    app = create_app([sys.argv[0]])
    yield app
