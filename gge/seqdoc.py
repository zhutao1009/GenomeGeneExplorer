"""带颜色的序列文档导出（Word ``.docx`` / 网页 ``.html``）。

为什么自己拼 OOXML
------------------
需求是"把界面上那张彩色序列表贴进 Word 里"。FASTA/TSV 是纯文本，天生带不了
底色；``python-docx`` 又只是把 OOXML 包了一层，为了几个色块再拖一个依赖进
PyInstaller 不划算。Word 真正需要的只有三个部件：

* ``[Content_Types].xml``（声明部件类型）
* ``_rels/.rels``（指向正文）
* ``word/document.xml``（正文）

文字级底色用 ``<w:shd w:val="clear" w:fill="RRGGBB"/>`` 表示，Word 2007
之后的版本都认（这也是"字符底纹"按钮存盘时写的东西）。

两种格式共用一份"行模型"
------------------------
:func:`document_lines` 把一次提取结果摊平成 ``DocLine`` 列表：标题、信息、
图例、面板标题、序列行（坐标 + 若干 ``(文本, 底色, 字色)`` 片段）。docx 与
html 两个渲染器都只消费这个列表，因此**颜色、坐标、10 bp 分组位置完全一致**；
以后要加 PDF 或剪贴板富文本，照抄一个渲染器即可。

三个着色模式
------------
* ``region``（默认）：只看区域（上游/外显子/CDS/UTR…），背景与字色都取界面配色。
  颜色变化次数等于"结构变化次数"，所以 run 极少、文件最小，和界面看到的一模一样；
* ``base``：只看碱基（A 绿 / T 红 / G 蓝 / C 橙），不带背景；
* ``both``：背景按区域、字母按碱基 —— 更好读，但每个碱基都可能自成一个 run，
  文件大小约是 ``region`` 的十几倍，几十 kb 以上要有心理准备。

一个容易被忽略的细节
--------------------
大序列不是"一个碱基一个 run"。同样颜色的连续碱基必须合并成一段，否则一个
50 kb 的面板会产生 5 万个 ``<w:r>``，Word 打开要转半天。合并后的 run 数量
与"颜色变化次数"成正比，而不是与碱基数成正比。
"""

from __future__ import annotations

import os
import zipfile
from dataclasses import dataclass, field
from datetime import datetime
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
from xml.sax.saxutils import escape as _xml_escape

from .colors import DEFAULT_PALETTE, LEGEND_ORDER, ensure_palette, region_label
from .coords import human_bytes
from .models import REGION_OTHER, Block, GeneExtraction, Panel
from .version import APP_NAME, __version__

# --------------------------------------------------------------------------- #
# 常量
# --------------------------------------------------------------------------- #
MODE_REGION = "region"
MODE_BASE = "base"
MODE_BOTH = "both"

MODE_LABELS = {
    MODE_REGION: "按区域着色（与界面一致，文件最小）",
    MODE_BOTH: "区域背景 + 碱基字母色（更好读，文件更大）",
    MODE_BASE: "只按碱基字母着色（A/T/G/C）",
}

#: 碱基配色（字母色）。U 用于 RNA，N/其他用灰，保证"有问题的碱基"不抢眼。
BASE_COLORS = {
    "A": "#15803d",
    "T": "#b91c1c",
    "U": "#0e7490",
    "G": "#1d4ed8",
    "C": "#b45309",
    "N": "#6b7280",
}
BASE_FALLBACK = "#111827"

GROUP_SIZE = 10          # 每 10 bp 一个空格，方便数位置
DEFAULT_BPL = 60
MIN_BPL, MAX_BPL = 20, 200
DEFAULT_FONT_PT = 8.5

#: 超过这个碱基数就提醒：不是不能导，而是 Word 会明显变慢
WARN_BASES = 200_000

DOCX_CONTENT_TYPES = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
<Default Extension="xml" ContentType="application/xml"/>
<Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>
</Types>"""

DOCX_RELS = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>
</Relationships>"""

_W_NS = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
_R_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
MONO_FONT = "Consolas"
CJK_FONT = "Microsoft YaHei"


def base_color(base: str) -> str:
    return BASE_COLORS.get((base or "").upper(), BASE_FALLBACK)


# --------------------------------------------------------------------------- #
# 行模型
# --------------------------------------------------------------------------- #
@dataclass
class DocLine:
    """文档里的一行。

    ``runs`` 里每项是 ``(文本, 底色, 字色)``；底色为空表示不加底纹。
    """

    kind: str = "line"          # title/info/legend/heading/subhead/note/line/spacer
    text: str = ""
    coord: str = ""
    runs: List[Tuple[str, str, str]] = field(default_factory=list)
    note: str = ""


def _palette_conf(palette: Dict[str, Dict[str, str]], kind: str) -> Dict[str, str]:
    return palette.get(kind) or palette.get(REGION_OTHER) or DEFAULT_PALETTE[REGION_OTHER]


class BlockPainter:
    """把一个 :class:`Block` 的着色结果切成 ``(文本, 底色, 字色)`` 片段。

    区间查找用游标而不是每行从头扫：调用方总是从左到右逐行取，游标只前进，
    因此整体是 O(区间数 + 字符数)，与"行数 × 区间数"无关。
    """

    def __init__(self, block: Block, mode: str = MODE_REGION,
                 palette: Optional[Dict[str, Dict[str, str]]] = None,
                 line_start: int = 0) -> None:
        self.block = block
        self.mode = mode if mode in MODE_LABELS else MODE_BOTH
        self.palette = ensure_palette(palette)
        self.regions = sorted(block.regions, key=lambda r: (r.start0, r.end0))
        self.line_start = line_start
        self._cursor = 0

    # -- 区间 -------------------------------------------------------------- #
    def spans(self, start: int, end: int) -> Iterable[Tuple[int, int, str]]:
        """把 ``[start, end)`` 切成 ``(起, 止, 区域类型)``，缺口补 other。"""
        regs = self.regions
        n = len(regs)
        if not n:
            yield start, end, REGION_OTHER
            return
        i = min(self._cursor, n - 1)
        while i > 0 and regs[i].start0 > start:
            i -= 1
        pos = start
        while pos < end:
            if i < n and regs[i].start0 <= pos:
                stop = min(end, regs[i].end0)
                kind = regs[i].kind
                if stop <= pos:                 # 空区间，防御性跳过
                    i += 1
                    continue
                if stop >= regs[i].end0:
                    i += 1
                yield pos, stop, kind
                pos = stop
            else:
                stop = min(end, regs[i].start0) if i < n else end
                if stop <= pos:
                    stop = end
                yield pos, stop, REGION_OTHER
                pos = stop
        self._cursor = max(0, i - 1)

    # -- 颜色 -------------------------------------------------------------- #
    def colors(self, kind: str) -> Tuple[str, str]:
        conf = _palette_conf(self.palette, kind)
        bg = conf.get("bg", "#ffffff")
        fg = conf.get("fg", "#111827")
        if self.mode == MODE_BASE:
            return "", fg
        return bg, fg

    def runs(self, start: int, end: int) -> Iterable[Tuple[str, str, str]]:
        """按"同色连续"合并后的片段序列（含每 10 bp 的分组空格）。"""
        seq = self.block.seq
        for s, e, kind in self.spans(start, end):
            bg, fg = self.colors(kind)
            x = s
            while x < e:
                off = x - self.line_start
                if off and off % GROUP_SIZE == 0:
                    yield " ", "", ""            # 分组空格，不带底纹
                group_end = min(e, self.line_start + (off // GROUP_SIZE + 1) * GROUP_SIZE)
                if self.mode == MODE_REGION:
                    yield seq[x:group_end], bg, fg
                else:
                    a = x
                    while a < group_end:
                        col = base_color(seq[a])
                        b = a + 1
                        while b < group_end and base_color(seq[b]) == col:
                            b += 1
                        yield seq[a:b], bg, col
                        a = b
                x = group_end


# --------------------------------------------------------------------------- #
# 行模型生成
# --------------------------------------------------------------------------- #
def info_rows(extraction: GeneExtraction) -> List[Tuple[str, str]]:
    """文档开头的元信息（与"复制基因基本信息"同源，但更短）。"""
    rows = [
        ("基因", f"{extraction.gene_name}（{extraction.gene_id}）"),
        ("转录本", f"{extraction.tx_id}（该基因共 {extraction.tx_count} 个转录本）"),
        ("位置", f"{extraction.seqid}:{extraction.gene_start:,}-{extraction.gene_end:,}"
                 f"（{extraction.strand} 链）"),
        ("结构", f"外显子 {extraction.n_exons} · 内含子 {extraction.n_introns} · "
                 f"CDS 片段 {extraction.n_cds}"),
        ("长度", f"mRNA {len(extraction.mrna):,} bp · CDS {len(extraction.cds):,} bp · "
                 f"5′UTR {len(extraction.utr5):,} bp · 3′UTR {len(extraction.utr3):,} bp"
                 + ("（UTR 非注释内容：注释中无 UTR 行，按「外显子 − CDS」推导）"
                    if not extraction.utr_in_annotation else
                    "（UTR 按「外显子 − CDS」推导，未直接采用注释中的 UTR 行）")),
        ("上下游", f"上游 {extraction.upstream_actual:,} bp · "
                   f"下游 {extraction.downstream_actual:,} bp"),
    ]
    if extraction.fasta_path:
        rows.append(("基因组", extraction.fasta_path))
    if extraction.annotation_path:
        rows.append(("注释", extraction.annotation_path))
    return rows


def _coord_width(extraction: GeneExtraction, panels: Sequence[Panel]) -> int:
    """坐标列宽：按最长坐标定，保证每一行的数字右对齐。"""
    top = max(extraction.region_end or 0, extraction.gene_end or 0, 1)
    for panel in panels:
        for block in panel.blocks:
            gi = block.genomic_interval
            if gi:
                top = max(top, max(gi))
            top = max(top, block.coord_start + len(block.seq))
    return max(7, len(f"{top:,}"))


def legend_runs(mode: str = MODE_BOTH,
                palette: Optional[Dict[str, Dict[str, str]]] = None
                ) -> List[Tuple[str, str, str]]:
    pal = ensure_palette(palette)
    runs: List[Tuple[str, str, str]] = []
    if mode in (MODE_REGION, MODE_BOTH):
        for kind in LEGEND_ORDER:
            conf = _palette_conf(pal, kind)
            runs.append(("███", conf["bg"], conf["bg"]))
            runs.append((f" {region_label(kind, pal)}   ", "", "#4b5563"))
    if mode in (MODE_BASE, MODE_BOTH):
        if runs:
            runs.append(("    ", "", ""))
        for base in ("A", "T", "G", "C", "N"):
            runs.append((base, "", BASE_COLORS[base]))
            runs.append((" ", "", ""))
        runs.append(("碱基配色", "", "#4b5563"))
    return runs


def document_lines(extraction: GeneExtraction,
                   panel_keys: Optional[Sequence[str]] = None,
                   mode: str = MODE_REGION,
                   bpl: int = DEFAULT_BPL,
                   palette: Optional[Dict[str, Dict[str, str]]] = None,
                   include_header: bool = True,
                   include_legend: bool = True,
                   generated_at: str = "",
                   include_coord: bool = True) -> List[DocLine]:
    """把一次提取摊平成行模型（docx / html / 剪贴板共用）。

    ``include_coord=False`` 时不输出每行左侧的坐标列（``DocLine.coord`` 留空），
    渲染器据此跳过坐标与它后面的两个空格。复制到论文/报告里时坐标往往用不上，
    却会占掉最左边一列、还让整行偏移 —— 所以界面上默认关闭。
    """
    bpl = max(MIN_BPL, min(MAX_BPL, int(bpl or DEFAULT_BPL)))
    pal = ensure_palette(palette)
    keys = list(panel_keys or [k for k in extraction.panels])
    panels = [extraction.panels[k] for k in keys
              if k in extraction.panels and not extraction.panels[k].is_empty()]
    out: List[DocLine] = []

    out.append(DocLine("title",
                       text=f"{extraction.gene_name}（{extraction.gene_id}）序列彩色图谱"))
    if include_header:
        for label, value in info_rows(extraction):
            out.append(DocLine("info", text=f"{label}：{value}",
                               note=label))
        out.append(DocLine("info", text=f"生成时间：{generated_at or _now()}"
                                       f"　·　{APP_NAME} v{__version__}"))
    if include_legend:
        out.append(DocLine("legend", runs=legend_runs(mode, pal)))

    total = sum(p.length for p in panels)
    out.append(DocLine("info", text=f"共 {len(panels)} 个面板 / {total:,} bp"
                                   f"　·　每行 {bpl} bp"))

    width = _coord_width(extraction, panels) if include_coord else 0
    for panel in panels:
        out.append(DocLine("heading",
                           text=f"{panel.title}（{panel.length:,} bp）"))
        if panel.note:
            out.append(DocLine("note", text=panel.note))
        multi = len(panel.blocks) > 1
        for block in panel.blocks:
            if not block.seq:
                continue
            if multi:
                head = f"▍{block.title}"
                if block.note:
                    head += f"　{block.note}"
                out.append(DocLine("subhead", text=head))
            elif block.note and not panel.note:
                out.append(DocLine("note", text=block.note))
            n_lines = (len(block.seq) + bpl - 1) // bpl
            painter = BlockPainter(block, mode, pal)
            for line in range(n_lines):
                start = line * bpl
                end = min(len(block.seq), start + bpl)
                painter.line_start = start
                coord = block.first_coord(line, bpl)
                out.append(DocLine(
                    "line",
                    coord=(f"{coord:>{width},}" if include_coord else ""),
                    runs=list(painter.runs(start, end)),
                    note=f"{len(block.seq):,}",
                ))
        out.append(DocLine("spacer"))
    return out


def _now() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M")


# --------------------------------------------------------------------------- #
# Word（.docx）
# --------------------------------------------------------------------------- #
def _hex(color: str) -> str:
    return (color or "").lstrip("#").upper()[:6]


def _default_rpr(size_pt: float, font: str = MONO_FONT) -> str:
    """段落默认运行属性（保留给需要"整段同格式"的场合）。

    注意：序列行**没有**用它。``w:pPr/w:rPr`` 按规范是"段落标记"的字符格式，
    Word 一般会继承，但不同版本/不同阅读器（WPS、LibreOffice、各种预览器）
    对它的继承并不一致；序列行即使丢字体也会整体错版，风险不值得省这点字节。
    因此每个碱基 run 都自带 ``rFonts`` + ``sz``，宁可文件大一点也要保证打开就对上。
    """
    half = max(2, int(round(size_pt * 2)))
    return (f'<w:rPr><w:rFonts w:ascii="{font}" w:hAnsi="{font}" '
            f'w:eastAsia="{CJK_FONT}" w:cs="{font}"/>'
            f'<w:sz w:val="{half}"/><w:szCs w:val="{half}"/></w:rPr>')


def _run(text: str, *, fill: str = "", color: str = "", size_pt: float = DEFAULT_FONT_PT,
         bold: bool = False, italic: bool = False, font: str = MONO_FONT) -> str:
    """完整的 run（字体/字号全写），用于标题、说明等零散文本。"""
    if not text:
        return ""
    props = [f'<w:rFonts w:ascii="{font}" w:hAnsi="{font}" '
             f'w:eastAsia="{CJK_FONT}" w:cs="{font}"/>']
    if bold:
        props.append("<w:b/><w:bCs/>")
    if italic:
        props.append("<w:i/><w:iCs/>")
    if color:
        props.append(f'<w:color w:val="{_hex(color)}"/>')
    half = max(2, int(round(size_pt * 2)))
    props.append(f'<w:sz w:val="{half}"/><w:szCs w:val="{half}"/>')
    if fill:
        props.append(f'<w:shd w:val="clear" w:color="auto" w:fill="{_hex(fill)}"/>')
    return (f"<w:r><w:rPr>{''.join(props)}</w:rPr>"
            f'<w:t xml:space="preserve">{_xml_escape(text)}</w:t></w:r>')


def _srun(text: str, fill: str = "", color: str = "",
          size_pt: float = DEFAULT_FONT_PT) -> str:
    """序列行专用 run：自带字体与字号（不依赖段落继承），其余按需写。"""
    if not text:
        return ""
    half = max(2, int(round(size_pt * 2)))
    props = [f'<w:rFonts w:ascii="{MONO_FONT}" w:hAnsi="{MONO_FONT}" '
             f'w:eastAsia="{CJK_FONT}"/>',
             f'<w:sz w:val="{half}"/>']
    if color:
        props.append(f'<w:color w:val="{_hex(color)}"/>')
    if fill:
        props.append(f'<w:shd w:val="clear" w:color="auto" w:fill="{_hex(fill)}"/>')
    return (f"<w:r><w:rPr>{''.join(props)}</w:rPr>"
            f'<w:t xml:space="preserve">{_xml_escape(text)}</w:t></w:r>')


def _par(runs: str, *, before: float = 0, after: float = 0, line_pt: float = 0,
         align: str = "", indent: float = 0, keep_next: bool = False,
         default_rpr: str = "") -> str:
    bits = []
    if keep_next:
        bits.append("<w:keepNext/>")
    if before or after or line_pt:
        line = f' w:line="{int(round(line_pt * 20))}" w:lineRule="exact"' if line_pt else ""
        bits.append(f'<w:spacing w:before="{int(round(before * 20))}"'
                    f' w:after="{int(round(after * 20))}"{line}/>')
    if indent:
        bits.append(f'<w:ind w:left="{int(round(indent * 20))}"/>')
    if align:
        bits.append(f'<w:jc w:val="{align}"/>')
    if default_rpr:                      # schema 要求 w:rPr 是 pPr 的最后一个子元素
        bits.append(default_rpr)
    ppr = f"<w:pPr>{''.join(bits)}</w:pPr>" if bits else ""
    return f"<w:p>{ppr}{runs}</w:p>"


def render_docx(lines: Sequence[DocLine], *, font_pt: float = DEFAULT_FONT_PT,
                title: str = "") -> bytes:
    """把行模型渲染成 ``.docx``（返回 zip 字节流）。"""
    body: List[str] = []
    for line in lines:
        if line.kind == "title":
            body.append(_par(_run(line.text, size_pt=15, bold=True, font=CJK_FONT,
                                  color="#111827"), align="center", after=4))
        elif line.kind == "info":
            body.append(_par(_run(line.text, size_pt=9, font=CJK_FONT, color="#374151"),
                             after=1))
        elif line.kind == "legend":
            runs = "".join(_run(t, fill=f, color=c, size_pt=9,
                                font=CJK_FONT if not t.strip("█") else MONO_FONT)
                           for t, f, c in line.runs)
            body.append(_par(runs, before=2, after=6))
        elif line.kind == "heading":
            body.append(_par(_run(line.text, size_pt=11.5, bold=True, font=CJK_FONT,
                                  color="#111827"), before=8, after=2, keep_next=True))
        elif line.kind == "subhead":
            body.append(_par(_run(line.text, size_pt=9, bold=True, font=CJK_FONT,
                                  color="#374151"), before=3, after=1, keep_next=True))
        elif line.kind == "note":
            body.append(_par(_run(line.text, size_pt=8, italic=True, font=CJK_FONT,
                                  color="#6b7280"), after=2))
        elif line.kind == "spacer":
            body.append(_par(_run(" ", size_pt=6)))
        else:                                    # 序列行
            runs = ""
            if line.coord:                       # 取消了坐标列时这里整段不写
                runs += _srun(line.coord, color="#9ca3af", size_pt=font_pt)
                runs += _srun("  ", size_pt=font_pt)
            for text, fill, color in line.runs:
                runs += _srun(text, fill=fill, color=color or "#111827", size_pt=font_pt)
            body.append(_par(runs, line_pt=font_pt * 1.25))

    sect = ('<w:sectPr><w:pgSz w:w="11906" w:h="16838"/>'
            '<w:pgMar w:top="1020" w:right="1020" w:bottom="1020" w:left="1020"'
            ' w:header="567" w:footer="567" w:gutter="0"/></w:sectPr>')
    document = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                f'<w:document xmlns:w="{_W_NS}" xmlns:r="{_R_NS}">'
                f'<w:body>{"".join(body)}{sect}</w:body></w:document>')
    return _zip_bytes({"[Content_Types].xml": DOCX_CONTENT_TYPES,
                       "_rels/.rels": DOCX_RELS,
                       "word/document.xml": document})


def _zip_bytes(parts: Dict[str, str]) -> bytes:
    import io

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        for name, text in parts.items():
            zf.writestr(name, text.encode("utf-8"))
    return buf.getvalue()


# --------------------------------------------------------------------------- #
# 网页（.html）
# --------------------------------------------------------------------------- #
HTML_STYLE = """
body { margin: 24px; background:#ffffff; color:#111827;
       font-family: "Microsoft YaHei", "PingFang SC", sans-serif; }
h1 { font-size: 17pt; margin: 0 0 6px; }
h2 { font-size: 12pt; margin: 22px 0 4px; padding-bottom: 3px;
     border-bottom: 1px solid #e5e7eb; }
h3 { font-size: 10pt; margin: 10px 0 3px; color:#374151; }
.meta { font-size: 9pt; color:#374151; line-height: 1.55; }
.meta b { display:inline-block; min-width: 4.6em; color:#6b7280; font-weight:600; }
.panel-note { font-size: 8.5pt; color:#6b7280; margin: 2px 0 4px; }
.legend { font-size: 9pt; margin: 8px 0 2px; }
pre.seq { font-family: Consolas, "Courier New", monospace; font-size: 8.5pt;
          line-height: 1.35; margin: 0; white-space: pre; }
pre.seq .c { color:#9ca3af; }
.foot { margin-top: 18px; font-size: 8pt; color:#9ca3af; }
""".strip()


def render_html_fragment(lines: Sequence[DocLine], *, font_pt: float = DEFAULT_FONT_PT,
                        full_page: bool = True, title: str = "") -> str:
    """把行模型渲染成 HTML。

    ``full_page=False`` 时只返回 ``<body>`` 内部片段 —— 这正是剪贴板富文本
    （``text/html``）需要的形态，Word / WPS 粘贴时按内联样式还原颜色。
    """
    out: List[str] = []
    in_seq = False
    for line in lines:
        if line.kind != "line" and in_seq:
            out.append("</pre>")
            in_seq = False
        if line.kind == "title":
            out.append(f"<h1>{_html_escape(line.text)}</h1>")
        elif line.kind == "info":
            label = line.note or ""
            text = line.text
            if label and text.startswith(label):
                out.append(f'<div class="meta"><b>{_html_escape(label)}</b>'
                           f'{_html_escape(text[len(label):].lstrip("：: "))}</div>')
            else:
                out.append(f'<div class="meta">{_html_escape(text)}</div>')
        elif line.kind == "legend":
            spans = "".join(_html_run(t, f, c) for t, f, c in line.runs)
            out.append(f'<div class="legend">{spans}</div>')
        elif line.kind == "heading":
            out.append(f"<h2>{_html_escape(line.text)}</h2>")
        elif line.kind == "subhead":
            out.append(f"<h3>{_html_escape(line.text)}</h3>")
        elif line.kind == "note":
            out.append(f'<div class="panel-note">{_html_escape(line.text)}</div>')
        elif line.kind == "spacer":
            out.append('<div style="height:10px"></div>')
        else:
            if not in_seq:
                out.append(f'<pre class="seq" style="font-size:{font_pt:g}pt">')
                in_seq = True
            runs = (f'<span class="c">{_html_escape(line.coord)}</span>  '
                    if line.coord else "")
            runs += "".join(_html_run(t, f, c) for t, f, c in line.runs)
            out.append(runs + "\n")
    if in_seq:
        out.append("</pre>")
    body = "".join(out)
    if not full_page:
        return body
    head = (f'<meta charset="utf-8"><title>{_html_escape(title or "序列彩色图谱")}</title>'
            f"<style>{HTML_STYLE}</style>")
    foot = ('<div class="foot">由 ' + APP_NAME + " v" + __version__ +
            " 生成；序列与注释文件原样读取，未做任何修改。</div>")
    return (f'<!DOCTYPE html><html lang="zh-CN"><head>{head}</head>'
            f"<body>{body}{foot}</body></html>")


def _html_escape(text: str) -> str:
    return (text or "").replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def _html_run(text: str, fill: str, color: str) -> str:
    if not text:
        return ""
    style = []
    # 写 background-color 而不是简写 background：Word 粘贴 HTML 时只认前者
    if fill:
        style.append(f"background-color:{fill}")
    if color:
        style.append(f"color:{color}")
    if not style:
        return _html_escape(text)
    return f'<span style="{";".join(style)}">{_html_escape(text)}</span>'


# --------------------------------------------------------------------------- #
# 落地
# --------------------------------------------------------------------------- #
DOC_SUFFIXES = (".docx", ".html", ".htm")


def suffix_of(path: str) -> str:
    low = os.path.splitext(path)[1].lower()
    return low if low in DOC_SUFFIXES else ".docx"


def estimate(extraction: GeneExtraction, panel_keys: Optional[Sequence[str]],
             bpl: int = DEFAULT_BPL, mode: str = MODE_REGION) -> Dict[str, int]:
    """粗略估算工作量，供界面提示"这次会导出多大"。"""
    keys = list(panel_keys or extraction.panels)
    bases = 0
    lines = 0
    panels = 0
    for k in keys:
        panel = extraction.panels.get(k)
        if panel is None or panel.is_empty():
            continue
        panels += 1
        bases += panel.length
        for block in panel.blocks:
            lines += (len(block.seq) + bpl - 1) // bpl
    # 每碱基的 XML 开销（实测）：区域着色只受"结构变化次数"影响，约 4 B/bp；
    # 按碱基着色要一个碱基一个 run，约 200 B/bp。docx 是 zip，再按 4 倍压缩折算。
    factor = {MODE_REGION: 4.0, MODE_BOTH: 210.0, MODE_BASE: 200.0}.get(mode, 4.0)
    raw = int(bases * factor + lines * 420)
    return {"panels": panels, "bases": bases, "lines": lines,
            "pages": max(1, round(lines / 62)),
            "docx_bytes": max(1024, int(raw * 0.25))}


def write_document(path: str, extraction: GeneExtraction,
                   panel_keys: Optional[Sequence[str]] = None,
                   mode: str = MODE_REGION, bpl: int = DEFAULT_BPL,
                   font_pt: float = DEFAULT_FONT_PT,
                   palette: Optional[Dict[str, Dict[str, str]]] = None,
                   include_header: bool = True,
                   include_legend: bool = True,
                   generated_at: str = "",
                   include_coord: bool = True) -> List[str]:
    """按扩展名写 ``.docx`` 或 ``.html``；返回实际写出的文件（可能两个）。"""
    suffix = suffix_of(path)
    stem = os.path.splitext(path)[0]
    lines = document_lines(extraction, panel_keys, mode=mode, bpl=bpl, palette=palette,
                           include_header=include_header, include_legend=include_legend,
                           generated_at=generated_at, include_coord=include_coord)
    title = f"{extraction.gene_id}_{extraction.tx_id} 序列彩色图谱"
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    written: List[str] = []
    if suffix == ".docx":
        with open(path, "wb") as fh:
            fh.write(render_docx(lines, font_pt=font_pt, title=title))
        written.append(path)
    else:
        html_path = path if suffix in (".html", ".htm") else stem + ".html"
        with open(html_path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(render_html_fragment(lines, font_pt=font_pt, full_page=True,
                                          title=title))
        written.append(html_path)
    return written


def clipboard_html(extraction: GeneExtraction,
                   panel_keys: Optional[Sequence[str]] = None,
                   mode: str = MODE_REGION, bpl: int = DEFAULT_BPL,
                   font_pt: float = DEFAULT_FONT_PT,
                   palette: Optional[Dict[str, Dict[str, str]]] = None,
                   include_coord: bool = True) -> str:
    """剪贴板用 HTML 片段：粘到 Word 里就是带底色的表格。"""
    lines = document_lines(extraction, panel_keys, mode=mode, bpl=bpl, palette=palette,
                           include_coord=include_coord)
    body = render_html_fragment(lines, font_pt=font_pt, full_page=False)
    return ("<html><head><meta charset=\"utf-8\"><style>" + HTML_STYLE +
            "</style></head><body>" + body + "</body></html>")


def describe(extraction: GeneExtraction, panel_keys: Optional[Sequence[str]],
             mode: str = MODE_REGION, bpl: int = DEFAULT_BPL) -> str:
    """一行摘要，用于成功提示。"""
    est = estimate(extraction, panel_keys, bpl, mode)
    return (f"{est['panels']} 个面板 / {est['bases']:,} bp / 约 {est['pages']:,} 页"
            f"，Word 文件约 {human_bytes(est['docx_bytes'])}")
