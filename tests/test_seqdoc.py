"""彩色序列文档（Word / 网页）导出的测试。

这一层最怕的不是"导出失败"，而是**导出的东西和界面对不上**：

* 颜色对了，碱基少了/多了；
* 行对得上，坐标错了（负链最容易出错）；
* Word 打开发现每个碱基都成了单独的 run（几十万个 run，打开要转半天）。

因此这里不只检查"文件生成成功"，而是把 document.xml 拆开，逐条比对
段落里的文字与原始序列是否逐字符相等。
"""

from __future__ import annotations

import os
import re
import xml.etree.ElementTree as ET
import zipfile

import pytest

from gge import seqdoc
from gge.extract import PANEL_CDS, PANEL_EXON, PANEL_GENOMIC, PANEL_MRNA

from conftest import extract_one

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"


def _res(provider, ann, up=400, dn=400):
    return extract_one(provider, ann, "g1", "g1.1", up=up, dn=dn)


def _docx_parts(data: bytes) -> dict:
    import io

    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        return {name: zf.read(name) for name in zf.namelist()}


def _paragraphs(xml: str):
    """返回 [(整段文字, 是否含底纹)]。"""
    root = ET.fromstring(xml)
    out = []
    for p in root.iter(f"{{{W}}}p"):
        runs = list(p.iter(f"{{{W}}}r"))
        text = "".join(t.text or "" for r in runs for t in r.iter(f"{{{W}}}t"))
        shaded = any(
            (r.find(f"{{{W}}}rPr/{{{W}}}shd") is not None) for r in runs)
        out.append((text, shaded))
    return out


def _seq_paragraphs(xml: str):
    """只挑"坐标 + 序列"这种段落。"""
    return [t for t, _s in _paragraphs(xml) if re.match(r"^\s*[\d,]+\s\s\S", t)]


def _bases(text: str) -> str:
    """去掉坐标前缀、分组空格与首尾空白，只留碱基。"""
    body = text.strip()
    body = re.sub(r"^[\d,]+\s+", "", body)
    return body.replace(" ", "")


# --------------------------------------------------------------------------- #
# 行模型
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("mode", list(seqdoc.MODE_LABELS))
@pytest.mark.parametrize("panel", [PANEL_GENOMIC, PANEL_MRNA, PANEL_CDS, PANEL_EXON])
def test_runs_reconstruct_exact_sequence(provider, ann_gff3, panel, mode):
    """把每一行的 run 拼回去，必须与 Block 里的序列逐字符相同。

    这条是整套导出的地基：颜色可以调，但碱基一个都不能少、不能串行。
    """
    res = _res(provider, ann_gff3)
    blk = res.panels[panel].blocks[0]
    bpl = 60
    painter = seqdoc.BlockPainter(blk, mode=mode)
    got = []
    for line in range((len(blk.seq) + bpl - 1) // bpl):
        start = line * bpl
        end = min(len(blk.seq), start + bpl)
        painter.line_start = start
        text = "".join(t for t, _f, _c in painter.runs(start, end))
        got.append(text.replace(" ", ""))
    assert "".join(got) == blk.seq


def test_group_separator_every_ten_bases(provider, ann_gff3):
    res = _res(provider, ann_gff3)
    blk = res.panels[PANEL_GENOMIC].blocks[0]
    painter = seqdoc.BlockPainter(blk, mode=seqdoc.MODE_REGION)
    painter.line_start = 0
    runs = list(painter.runs(0, 60))
    spaces = [i for i, (t, _f, _c) in enumerate(runs) if t == " "]
    # 60 bp 一行、每 10 bp 一个空格 → 5 个分隔空格（不含行尾）
    assert len(spaces) == 5
    assert all(fill == "" for t, fill, _c in runs if t == " "), "分组空格不能带底纹"


def test_region_mode_merges_runs(provider, ann_gff3):
    """区域着色时，一段连续的外显子只能是 1 个 run（不能被拆成每个碱基一个）。"""
    res = _res(provider, ann_gff3)
    blk = res.panels[PANEL_GENOMIC].blocks[0]
    painter = seqdoc.BlockPainter(blk, mode=seqdoc.MODE_REGION)
    painter.line_start = 0
    runs = [r for r in painter.runs(0, min(60, len(blk.seq)))]
    # 60 个碱基整块同色时：2 个 run + 5 个空格
    assert len(runs) <= 12, runs


def test_base_mode_colors_four_bases(provider, ann_gff3):
    res = _res(provider, ann_gff3)
    blk = res.panels[PANEL_GENOMIC].blocks[0]
    painter = seqdoc.BlockPainter(blk, mode=seqdoc.MODE_BASE)
    painter.line_start = 0
    seen = {}
    for t, fill, color in painter.runs(0, 200):
        if t.strip():
            assert fill == "", "只按碱基着色时不应有底纹"
        for ch in set(t.upper()):
            if ch in "ACGT":
                seen.setdefault(ch, color)
    assert len(set(seen.values())) >= 3, f"四种碱基应至少有 3 种颜色，实际 {seen}"


def test_coord_column_follows_strand(provider, ann_gff3):
    """坐标列的起点必须来自 Block，而不是"从 1 数"。"""
    res = _res(provider, ann_gff3)
    lines = seqdoc.document_lines(res, [PANEL_GENOMIC], bpl=60)
    seq_lines = [ln for ln in lines if ln.kind == "line"]
    blk = res.panels[PANEL_GENOMIC].blocks[0]
    assert seq_lines[0].coord == f"{blk.coord_start:,}".rjust(len(seq_lines[0].coord))
    assert seq_lines[1].coord.strip() == f"{blk.first_coord(1, 60):,}"


def test_minus_strand_coord_counts_down(provider, ann_gff3):
    """负链块的坐标应逐行变小 —— 用人工 Block 钉死，不依赖示例基因的方向。"""
    from gge.models import Block, Region, REGION_CDS

    blk = Block("b", "t", "ACGT" * 30,
                [Region(0, 120, REGION_CDS, "CDS")],
                meta={"coord_start": 5000, "coord_reverse": True,
                      "coord_mode": "genomic"})
    assert blk.first_coord(0, 60) == 5000 + 120 - 1
    assert blk.first_coord(1, 60) == 5000 + 120 - 1 - 60
    assert blk.first_coord(2, 60) < blk.first_coord(1, 60)


# --------------------------------------------------------------------------- #
# docx
# --------------------------------------------------------------------------- #
def test_docx_is_minimal_valid_package(provider, ann_gff3):
    res = _res(provider, ann_gff3)
    lines = seqdoc.document_lines(res, [PANEL_GENOMIC], mode=seqdoc.MODE_REGION)
    parts = _docx_parts(seqdoc.render_docx(lines))
    assert set(parts) == {"[Content_Types].xml", "_rels/.rels", "word/document.xml"}
    for name, blob in parts.items():
        ET.fromstring(blob)          # 每一份都必须是良构 XML
        assert name
    doc = parts["word/document.xml"].decode("utf-8")
    assert doc.count("<w:body>") == 1 and doc.count("</w:body>") == 1
    assert "<w:sectPr>" in doc          # 没有节属性，Word 会用奇怪的页面尺寸


def test_docx_text_matches_sequence(provider, ann_gff3):
    res = _res(provider, ann_gff3)
    lines = seqdoc.document_lines(res, [PANEL_MRNA], mode=seqdoc.MODE_REGION, bpl=60)
    doc = _docx_parts(seqdoc.render_docx(lines))["word/document.xml"].decode("utf-8")
    got = "".join(_bases(t) for t in _seq_paragraphs(doc))
    assert got == res.panels[PANEL_MRNA].sequence


def test_docx_has_character_shading(provider, ann_gff3):
    res = _res(provider, ann_gff3)
    lines = seqdoc.document_lines(res, [PANEL_GENOMIC], mode=seqdoc.MODE_REGION)
    doc = _docx_parts(seqdoc.render_docx(lines))["word/document.xml"].decode("utf-8")
    assert "w:shd" in doc
    fills = set(re.findall(r'w:fill="([0-9A-F]{6})"', doc))
    # 基因组面板含上游/外显子/内含子等，至少 3 种底色
    assert len(fills) >= 3, fills


def test_docx_runs_carry_font(provider, ann_gff3):
    """每个序列 run 都要自带等宽字体，不能指望段落继承（WPS/预览器不一致）。"""
    res = _res(provider, ann_gff3)
    lines = seqdoc.document_lines(res, [PANEL_MRNA], mode=seqdoc.MODE_REGION)
    doc = _docx_parts(seqdoc.render_docx(lines))["word/document.xml"].decode("utf-8")
    assert doc.count("Consolas") >= doc.count("<w:r><w:rPr>")


def test_base_mode_makes_bigger_file(provider, ann_gff3):
    res = _res(provider, ann_gff3)
    keys = [PANEL_GENOMIC, PANEL_MRNA]
    region = seqdoc.render_docx(seqdoc.document_lines(res, keys, mode="region"))
    both = seqdoc.render_docx(seqdoc.document_lines(res, keys, mode="both"))
    assert len(both) > len(region)


def test_title_is_xml_escaped(provider, ann_gff3):
    """基因名里出现 & < > 时不能把 document.xml 弄坏。"""
    res = _res(provider, ann_gff3)
    res.gene_name = "A&B <test>"
    lines = seqdoc.document_lines(res, [PANEL_CDS])
    doc = _docx_parts(seqdoc.render_docx(lines))["word/document.xml"].decode("utf-8")
    assert "A&amp;B &lt;test&gt;" in doc
    ET.fromstring(doc)


def test_legend_lists_regions(provider, ann_gff3):
    res = _res(provider, ann_gff3)
    lines = seqdoc.document_lines(res, [PANEL_GENOMIC], mode=seqdoc.MODE_REGION)
    legend = [ln for ln in lines if ln.kind == "legend"]
    assert legend, "区域着色时必须有图例"
    text = "".join(t for t, _f, _c in legend[0].runs)
    for word in ("上游", "外显子", "内含子"):
        assert word in text


# --------------------------------------------------------------------------- #
# html / 落盘
# --------------------------------------------------------------------------- #
def test_html_fragment_has_backgrounds(provider, ann_gff3):
    res = _res(provider, ann_gff3)
    lines = seqdoc.document_lines(res, [PANEL_GENOMIC], mode=seqdoc.MODE_REGION)
    html = seqdoc.render_html_fragment(lines)
    assert html.startswith("<!DOCTYPE html>")
    assert "background-color:#" in html
    assert "<pre class=\"seq\"" in html
    # 碱基必须原样出现：先摘掉坐标列，再去标签、去分组空格，逐字符比对
    body = re.sub(r'<span class="c">[^<]*</span>', "", html)
    plain = re.sub(r"<[^>]+>", "", body).replace(" ", "").replace("\n", "")
    seq = res.panels[PANEL_GENOMIC].blocks[0].seq
    assert seq in plain


def test_clipboard_html_is_a_complete_fragment(provider, ann_gff3):
    res = _res(provider, ann_gff3)
    html = seqdoc.clipboard_html(res, [PANEL_MRNA], mode=seqdoc.MODE_BOTH)
    assert html.startswith("<html>") and html.rstrip().endswith("</body></html>")
    assert "<!DOCTYPE" not in html, "剪贴板里不需要 DOCTYPE"
    assert "<pre" in html and "<span" in html


def test_write_document_by_suffix(provider, ann_gff3, tmp_path):
    res = _res(provider, ann_gff3)
    docx = tmp_path / "a.docx"
    files = seqdoc.write_document(str(docx), res, [PANEL_MRNA])
    assert files == [str(docx)] and docx.exists() and docx.stat().st_size > 500
    html = tmp_path / "b.html"
    files = seqdoc.write_document(str(html), res, [PANEL_MRNA])
    assert files == [str(html)]
    assert "background-color" in html.read_text(encoding="utf-8")


def test_estimate_matches_line_count(provider, ann_gff3):
    res = _res(provider, ann_gff3)
    est = seqdoc.estimate(res, [PANEL_EXON], bpl=50)
    n = sum((len(b.seq) + 49) // 50 for b in res.panels[PANEL_EXON].blocks)
    assert est["lines"] == n
    assert est["bases"] == res.panels[PANEL_EXON].length
    big = seqdoc.estimate(res, [PANEL_EXON], bpl=50, mode=seqdoc.MODE_BOTH)
    small = seqdoc.estimate(res, [PANEL_EXON], bpl=50, mode=seqdoc.MODE_REGION)
    assert big["docx_bytes"] > small["docx_bytes"] * 5


def test_empty_panel_is_skipped(provider, ann_gff3):
    res = _res(provider, ann_gff3)
    lines = seqdoc.document_lines(res, ["not_a_panel", PANEL_MRNA])
    headings = [ln.text for ln in lines if ln.kind == "heading"]
    assert len(headings) == 1 and "mRNA" in headings[0]


def test_human_bytes_vs_human_size():
    """文件大小必须说 B/KB/MB，碱基数才说 bp/kb/Mb。"""
    from gge.coords import human_bytes, human_size

    assert human_bytes(136) == "136 B"
    assert human_bytes(21_948) == "21.4 KB"
    assert human_bytes(905_969_664) == "864.0 MB"
    assert human_size(136) == "136 bp"
    assert human_size(21_948) == "21.9 kb"
    # 同一个数字，两种单位不能混用
    assert human_bytes(136) != human_size(136)


# --------------------------------------------------------------------------- #
# 坐标列可关闭（论文/报告里通常不需要）
# --------------------------------------------------------------------------- #
def test_document_lines_without_coord(provider, ann_gff3):
    """关掉坐标后，序列行的 ``coord`` 必须为空，但碱基片段不能有任何变化。"""
    res = _res(provider, ann_gff3)
    keys = [PANEL_GENOMIC, PANEL_MRNA]
    with_c = [ln for ln in seqdoc.document_lines(res, keys, include_coord=True)
              if ln.kind == "line"]
    without = [ln for ln in seqdoc.document_lines(res, keys, include_coord=False)
               if ln.kind == "line"]
    assert with_c and len(with_c) == len(without)
    assert all(ln.coord.strip() for ln in with_c)
    assert all(ln.coord == "" for ln in without)
    assert [ln.runs for ln in with_c] == [ln.runs for ln in without], \
        "关掉坐标不能顺手改掉颜色或碱基分组"


def test_docx_without_coord(provider, ann_gff3):
    """docx 里不能再有灰色的坐标 run；每行文字 = 原行去掉坐标前缀。"""
    res = _res(provider, ann_gff3)
    keys = [PANEL_GENOMIC, PANEL_MRNA]
    on = _docx_parts(seqdoc.render_docx(
        seqdoc.document_lines(res, keys, include_coord=True)))["word/document.xml"].decode()
    off = _docx_parts(seqdoc.render_docx(
        seqdoc.document_lines(res, keys, include_coord=False)))["word/document.xml"].decode()

    assert 'w:val="9CA3AF"' in on, "带坐标时应该有浅灰的坐标 run"
    assert 'w:val="9CA3AF"' not in off, "关掉坐标后不该再出现坐标 run"

    p_on = [t for t, _s in _paragraphs(on)]
    p_off = [t for t, _s in _paragraphs(off)]
    assert len(p_on) == len(p_off), "段落数应与带坐标时完全一致"

    pairs = 0
    for t_on, t_off in zip(p_on, p_off):
        if re.match(r"^\s*[\d,]+\s\s", t_on):
            pairs += 1
            assert re.sub(r"^\s*[\d,]+\s+", "", t_on) == t_off
    assert pairs > 0, "没找到任何序列行，测试无效"


def test_html_without_coord(provider, ann_gff3):
    res = _res(provider, ann_gff3)
    html_on = seqdoc.render_html_fragment(
        seqdoc.document_lines(res, [PANEL_GENOMIC], include_coord=True), full_page=False)
    html_off = seqdoc.render_html_fragment(
        seqdoc.document_lines(res, [PANEL_GENOMIC], include_coord=False), full_page=False)
    assert '<span class="c">' in html_on
    assert '<span class="c">' not in html_off
    # 底色片段不能跟着丢
    assert html_on.count("background-color") == html_off.count("background-color")
