"""主窗口：把文件、搜索、提取、可视化、校验、导出串起来。"""

from __future__ import annotations

import os
import sys
from collections import OrderedDict
from typing import Dict, List, Optional

from PySide6.QtCore import (QByteArray, QMimeData, QSettings, QSize, QStandardPaths,
                            Qt, QTimer)
from PySide6.QtGui import QAction, QGuiApplication, QIcon, QKeySequence, QShortcut
from PySide6.QtWidgets import (QApplication, QCheckBox, QComboBox, QDialog,
                               QFileDialog, QFrame, QHBoxLayout, QLabel, QLineEdit,
                               QMainWindow, QMessageBox, QProgressBar, QPushButton,
                               QScrollArea, QSizePolicy, QSpinBox, QSplitter,
                               QToolButton, QVBoxLayout, QWidget)

from ..bundle import ROLE_LABELS, ROLE_ORDER, describe_roles, is_archive, scan_archive
from ..colors import LEGEND_ORDER, UI_COLORS, ensure_palette
from ..config import Settings
from ..errors import CancelledError, GGEError, NotReadyError
from ..library import Dataset, Library
from ..export import (gene_basic_info_text, gene_fasta_text, panel_fasta_text,
                      safe_filename, validation_report_text, write_text)
from ..extract import PANEL_TITLES
from ..models import GeneExtraction
from ..search import SearchHit, search_genes
from ..seqdoc import MODE_REGION, clipboard_html
from ..tasks import (BatchResult, BundleImportResult, DatasetExportResult,
                     DocumentExportResult, ExtractionBundle, GenomeInfo,
                     batch_export, export_dataset_archive, export_gene,
                     export_sequence_document, extract_and_validate,
                     import_bundle_from_archive, load_annotation_file, load_genome,
                     load_reference_file)
from ..version import APP_NAME, APP_NAME_CN, __version__
from ..validate import CheckResult
from .bundle_dialog import BundleImportDialog, LibraryDialog
from .dialogs import (AboutDialog, BatchExportDialog, ColorSettingsDialog,
                      DatasetExportDialog, DiffDialog, SequenceDocDialog)
from .sequence_view import LegendBar, SequencePanel
from .side_panel import FilePanel, InfoPanel, SearchPanel
from .structure_view import GeneStructureView
from .theme import mono_font, ui_font
from .validation_panel import ValidationPanel
from .workers import Job

PANEL_LAYOUT = [
    ["genomic", "exon", "intron", "cds"],
    ["upstream", "downstream", "mrna", "utr"],
]

ABOUT_HTML = """
<h2>Genome Gene Explorer</h2>
<p style="color:#6b7280">基因组序列提取与注释校验工具 v{ver}</p>
<p>纯本地运行，不联网、不依赖在线数据库。</p>
<h3>技术栈</h3>
<ul>
<li>界面：PySide6（Qt 6）</li>
<li>大基因组随机访问：pyfaidx（.fai 索引）</li>
<li>序列比对：Biopython PairwiseAligner（全局比对）</li>
<li>数据集库：标准库 sqlite3（不额外装东西）</li>
<li>自动化测试：pytest（{tests} 项）</li>
</ul>
<h3>目录</h3>
<p>源码与文档见程序目录下的 README.md 与 docs/ALGORITHMS.md。</p>
<p style="color:#9ca3af">设置文件：%APPDATA%\\GenomeGeneExplorer\\settings.json<br>
数据集库：%APPDATA%\\GenomeGeneExplorer\\library\\</p>
"""

HELP_HTML = """
<h2>基因组序列提取与注释校验工具</h2>
<p style="color:#6b7280">Genome Gene Explorer v{ver}</p>
<h3>基本流程</h3>
<ol>
<li><b>选择基因组 FASTA</b>（支持 .fa/.fasta/.fna，自动建立 .fai 索引，几十 GB 也只占极小内存）。</li>
<li><b>选择注释文件</b>（GFF3 / GTF，支持 .gz）。</li>
<li>在搜索框输入 Gene Symbol / Gene ID / Transcript ID，双击结果选中基因。</li>
<li>在上方结构图与下方 8 个序列面板中查看结果；每个面板可独立复制、导出。</li>
<li>如需校验，加载参考 mRNA / CDS / 蛋白 FASTA，程序会按 ID 打分匹配并做全局比对。</li>
<li>要把带底色的序列发给别人：面板右上角「彩色文档导出」→ 选 Word 或网页；
   也可以右键「复制为彩色序列」，直接粘进 Word。文档默认<b>不带左侧坐标列</b>，
   需要按位置核对时在导出对话框里勾选「每行左侧显示坐标列」。</li>
</ol>
<p><b>下次打开</b>：程序会自动把上次用的文件加载好，并回到上次查看的基因，
不用再按 Ctrl+R。想改成"只填路径、不自动加载"（比如超大基因组），
在「文件」菜单里取消勾选「启动时自动加载上次的文件」。</p>
<p><b>懒人路线</b>：把上面这几个文件打成一个 zip / tar.gz，菜单「文件 → 导入压缩包」，
程序会按内容（不是文件名）自动认出谁是基因组、谁是注释、谁是 mRNA/CDS/蛋白，解包并登记到数据集库；
下次从「文件 → 数据集库」双击即可载入。</p>
<p><b>反过来也行</b>：手动逐个选好文件后，用「文件 → 导出数据集压缩包」把它们打成一个 zip
（内含 gge_dataset.json 清单）。下次「导入压缩包」选中它，会按清单直接认回角色，
不必再猜、也不必重新统计注释，打开更快。</p>
<h3>快捷键</h3>
<ul>
<li><b>Ctrl+O</b> 打开基因组；<b>Ctrl+Shift+O</b> 打开注释</li>
<li><b>Ctrl+Shift+I</b> 导入数据压缩包；<b>Ctrl+L</b> 打开数据集库</li>
<li><b>Ctrl+F</b> 在所有序列面板中查找；<b>F3</b> 下一个匹配</li>
<li><b>Ctrl+Shift+C</b> 复制当前面板全部序列；<b>Ctrl+Shift+D</b> 导出彩色序列文档</li>
<li><b>Ctrl+Shift+E</b> 把当前这套文件导出成数据集压缩包</li>
<li><b>F5</b> 重新提取 / 重新校验</li>
<li>序列视图内：<b>Ctrl+滚轮</b> 缩放字号，<b>Shift+滚轮</b> 横向滚动，<b>拖动</b> 选择碱基，<b>Ctrl+C</b> 复制选区</li>
</ul>
<h3>颜色优先级</h3>
<p>CDS &gt; UTR &gt; 外显子 &gt; 内含子 &gt; 上下游 &gt; 基因间。因此落在外显子内部的 CDS 会显示为
CDS 颜色，而不是外显子颜色。</p>
<h3>数据正确性</h3>
<ul>
<li>注释坐标是 1-based 闭区间，内部切片转成 0-based 半开区间，只在读取处转换一次。</li>
<li>负链基因的上下游按<b>转录方向</b>定义并反向互补，不是按坐标方向。</li>
<li>CDS 的 phase 只用于报告与一致性检查，<b>不会据此裁剪碱基</b>。</li>
<li>参考序列不匹配时只报告差异，<b>绝不修改</b>提取到的序列。</li>
<li>FASTA 与注释的序列名对不上时直接报错并给出相近名称，不会静默返回空序列。</li>
</ul>
<h3>示例数据</h3>
<p>菜单「帮助 → 加载示例数据」会用程序自带的小基因组（约 21 kb）演示全部功能。</p>
<p>想试压缩包导入，用程序目录下的 <code>examples/demo_bundle.zip</code>——里面 5 个文件的文件名是故意起得没有信息量的。</p>
"""


class MainWindow(QMainWindow):
    def __init__(self, settings: Optional[Settings] = None,
                 auto_restore: bool = True) -> None:
        super().__init__()
        self.settings = settings or Settings()
        self.palette_map = ensure_palette(self.settings.palette())
        #: 启动时是否把上次的文件自动加载回来（``--self-test`` / 命令行指定文件时关掉）
        self._auto_restore_on_start = bool(auto_restore)
        #: 自动恢复期间不弹"注释有提醒"的对话框，免得开机就挡用户一脸
        self._restoring = False

        self.provider = None
        self.genome_info: Optional[GenomeInfo] = None
        self.model = None
        self.mrna_refs = None
        self.cds_refs = None
        self.protein_refs = None

        self.extraction: Optional[GeneExtraction] = None
        self.validation: Dict[str, CheckResult] = {}
        self.current_gene_id = ""
        self.current_tx_id = ""
        self._manual: Dict[str, str] = {"mRNA": "", "CDS": "", "protein": ""}

        self.library = Library()
        self.current_dataset: Optional[Dataset] = None
        self._dataset_id = 0

        self._cache: "OrderedDict[tuple, ExtractionBundle]" = OrderedDict()
        self._cache_limit = 8
        self._job: Optional[Job] = None
        self._job_done = None
        self._queue: List[dict] = []
        self._maximized: Optional[str] = None
        self._col_sizes: List[List[int]] = []
        self._search_order: List[str] = []
        self._search_idx = 0

        self.setWindowTitle(f"{APP_NAME} · {APP_NAME_CN}")
        self.resize(1560, 940)
        self.setAcceptDrops(True)
        self._build_ui()
        self._build_actions()
        self._restore_settings()
        self._update_actions()
        self._schedule_auto_restore()

    # ================================================================== #
    # 构建界面
    # ================================================================== #
    def _build_ui(self) -> None:
        splitter = QSplitter(Qt.Horizontal)
        splitter.setChildrenCollapsible(False)
        self.setCentralWidget(splitter)

        # ---------------- 左栏 ---------------- #
        left = QWidget()
        left.setMinimumWidth(330)
        left.setMaximumWidth(520)
        lv = QVBoxLayout(left)
        lv.setContentsMargins(8, 8, 4, 8)
        lv.setSpacing(8)

        self.file_panel = FilePanel()
        self.file_panel.browseRequested.connect(self._browse_role)
        self.file_panel.clearRequested.connect(self._clear_role)
        self.file_panel.importRequested.connect(lambda: self.import_archive_dialog())
        self.file_panel.libraryRequested.connect(self.open_library)
        self.file_panel.exportRequested.connect(self.export_dataset_dialog)
        lv.addWidget(self.file_panel)

        self.search_panel = SearchPanel()
        self.search_panel.searchRequested.connect(self.do_search)
        self.search_panel.geneChosen.connect(self.select_gene)
        lv.addWidget(self.search_panel, 1)

        self.info_panel = InfoPanel()
        self.info_panel.transcriptChanged.connect(self._on_transcript_changed)
        lv.addWidget(self.info_panel, 1)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(left)
        scroll.setFrameShape(QFrame.NoFrame)
        scroll.setMinimumWidth(330)
        splitter.addWidget(scroll)

        # ---------------- 右栏 ---------------- #
        right = QWidget()
        rv = QVBoxLayout(right)
        rv.setContentsMargins(4, 8, 8, 8)
        rv.setSpacing(6)

        rv.addWidget(self._build_control_bar())
        self.legend_bar = LegendBar(self.palette_map)
        rv.addWidget(self.legend_bar)

        self.vsplit = QSplitter(Qt.Vertical)
        self.vsplit.setChildrenCollapsible(False)

        # 结构图
        struct_frame = QFrame()
        struct_frame.setObjectName("structFrame")
        struct_frame.setStyleSheet(
            "QFrame#structFrame{background:#ffffff;border:1px solid #e5e7eb;border-radius:6px;}")
        sv = QVBoxLayout(struct_frame)
        sv.setContentsMargins(4, 4, 4, 2)
        sv.setSpacing(0)
        sbar = QHBoxLayout()
        sbar.setContentsMargins(6, 2, 4, 0)
        t = QLabel("基因结构")
        f = ui_font(9)
        f.setBold(True)
        t.setFont(f)
        sbar.addWidget(t)
        self.struct_hint = QLabel("滚轮缩放 · 拖动平移 · 双击复位")
        self.struct_hint.setStyleSheet("color:#9ca3af; font-size:10px;")
        sbar.addWidget(self.struct_hint)
        sbar.addStretch(1)
        btn_fit = QToolButton()
        btn_fit.setText("适配")
        btn_fit.setToolTip("完整显示整个窗口区间")
        btn_fit.clicked.connect(lambda: self.structure.fit())
        btn_png = QToolButton()
        btn_png.setText("导出图")
        btn_png.setToolTip("导出为 PNG / SVG")
        btn_png.clicked.connect(self.export_structure)
        for b in (btn_fit, btn_png):
            b.setStyleSheet(
                "QToolButton{border:1px solid transparent;border-radius:4px;padding:2px 6px;"
                "color:#374151;font-size:11px;}"
                "QToolButton:hover{background:#eef2ff;border-color:#c7d2fe;}")
            sbar.addWidget(b)
        sv.addLayout(sbar)
        self.structure = GeneStructureView(self.palette_map)
        sv.addWidget(self.structure, 1)
        self.vsplit.addWidget(struct_frame)

        # 序列面板
        self.panels_splitter = QSplitter(Qt.Orientation.Horizontal)
        self.panels_splitter.setOrientation(Qt.Orientation.Horizontal)
        self.panels_splitter.setChildrenCollapsible(False)
        self.panels: Dict[str, SequencePanel] = {}
        self._col_splitters: List[QSplitter] = []
        for keys in PANEL_LAYOUT:
            col = QSplitter(Qt.Orientation.Vertical)
            col.setOrientation(Qt.Orientation.Vertical)
            col.setChildrenCollapsible(False)
            for key in keys:
                panel = SequencePanel(key, PANEL_TITLES[key], self.palette_map, UI_COLORS)
                panel.requestExport.connect(self.export_panel)
                panel.requestCopyFasta.connect(self.copy_panel_fasta)
                panel.requestRichCopy.connect(self.copy_panel_rich)
                panel.requestDocExport.connect(self.export_sequence_document_dialog)
                panel.requestMaximize.connect(self._toggle_maximize)
                panel.cursorInfo.connect(self._on_cursor_info)
                col.addWidget(panel)
                self.panels[key] = panel
            self._col_splitters.append(col)
            self.panels_splitter.addWidget(col)

        self.validation_panel = ValidationPanel()
        self.validation_panel.detailRequested.connect(self._show_diff_dialog)
        self.validation_panel.manualChanged.connect(self._on_manual_ref)
        self.validation_panel.revalidateRequested.connect(self.extract_current)
        self._col_splitters[1].addWidget(self.validation_panel)

        for sp in [self.vsplit, self.panels_splitter] + self._col_splitters:
            sp.setHandleWidth(6)
        self.vsplit.addWidget(self.panels_splitter)
        self._apply_default_sizes()
        rv.addWidget(self.vsplit, 1)
        splitter.addWidget(right)
        splitter.setSizes([380, 1180])

        # ---------------- 状态栏 ---------------- #
        sb = self.statusBar()
        self.status_msg = QLabel("就绪。请先选择基因组 FASTA 与注释文件。")
        self.status_msg.setStyleSheet("color:#4b5563; font-size:11px;")
        sb.addWidget(self.status_msg, 1)
        self.progress = QProgressBar()
        self.progress.setFixedWidth(180)
        self.progress.setTextVisible(True)
        self.progress.hide()
        sb.addPermanentWidget(self.progress)
        self.btn_cancel = QPushButton("取消")
        self.btn_cancel.setFixedWidth(56)
        self.btn_cancel.hide()
        self.btn_cancel.clicked.connect(self._cancel_job)
        sb.addPermanentWidget(self.btn_cancel)

        # 默认空状态
        for p in self.panels.values():
            p.clear()
        self._apply_font_size(self.settings.get("font_size", 10))
        self._apply_bpl(self.settings.get("bases_per_line", 60))
        self._apply_auto_fit(bool(self.settings.get("auto_fit", True)))

    def _build_control_bar(self) -> QWidget:
        bar = QFrame()
        bar.setObjectName("ctrlBar")
        bar.setStyleSheet(
            "QFrame#ctrlBar{background:#ffffff;border:1px solid #e5e7eb;border-radius:6px;}"
        )
        lay = QHBoxLayout(bar)
        lay.setContentsMargins(8, 4, 8, 4)
        lay.setSpacing(6)

        def lab(text: str) -> QLabel:
            l = QLabel(text)
            l.setStyleSheet("color:#6b7280; font-size:11px;")
            return l

        self.spin_up = QSpinBox()
        self.spin_up.setRange(0, 1_000_000)
        self.spin_up.setSingleStep(500)
        self.spin_up.setValue(int(self.settings.get("upstream", 2000)))
        self.spin_up.setFixedWidth(76)
        self.spin_up.setToolTip("基因上游长度（按转录方向定义）")
        self.spin_dn = QSpinBox()
        self.spin_dn.setRange(0, 1_000_000)
        self.spin_dn.setSingleStep(500)
        self.spin_dn.setValue(int(self.settings.get("downstream", 2000)))
        self.spin_dn.setFixedWidth(76)
        self.spin_dn.setToolTip("基因下游长度（按转录方向定义）")

        self.btn_extract = QPushButton("提取")
        self.btn_extract.setObjectName("primary")
        self.btn_extract.setFixedWidth(72)
        self.btn_extract.setToolTip("按当前上下游长度重新提取（F5）")
        self.btn_extract.clicked.connect(self.extract_current)

        lay.addWidget(lab("上游"))
        lay.addWidget(self.spin_up)
        lay.addWidget(lab("下游"))
        lay.addWidget(self.spin_dn)
        lay.addWidget(self.btn_extract)

        line = QFrame()
        line.setFrameShape(QFrame.VLine)
        line.setStyleSheet("color:#e5e7eb;")
        lay.addWidget(line)

        self.find_edit = QLineEdit()
        self.find_edit.setPlaceholderText("在所有面板中查找序列 (Ctrl+F)")
        self.find_edit.setFixedWidth(160)
        self.find_edit.returnPressed.connect(self._do_find)
        self.btn_find_prev = QToolButton()
        self.btn_find_prev.setText("↑")
        self.btn_find_prev.setToolTip("上一个匹配 (Shift+F3)")
        self.btn_find_prev.clicked.connect(lambda: self._step_find(-1))
        self.btn_find_next = QToolButton()
        self.btn_find_next.setText("↓")
        self.btn_find_next.setToolTip("下一个匹配 (F3)")
        self.btn_find_next.clicked.connect(lambda: self._step_find(1))
        self.find_info = QLabel("")
        self.find_info.setStyleSheet("color:#6b7280; font-size:10px;")
        lay.addWidget(self.find_edit)
        lay.addWidget(self.btn_find_prev)
        lay.addWidget(self.btn_find_next)
        lay.addWidget(self.find_info)

        lay.addStretch(1)

        self.chk_auto = QCheckBox("自适应")
        self.chk_auto.setToolTip("每行碱基数按面板宽度自动计算，避免长行被右侧裁掉")
        self.chk_auto.setChecked(bool(self.settings.get("auto_fit", True)))
        self.chk_auto.toggled.connect(self._apply_auto_fit)
        lay.addWidget(self.chk_auto)
        self.spin_bpl = QSpinBox()
        self.spin_bpl.setRange(20, 200)
        self.spin_bpl.setSingleStep(10)
        self.spin_bpl.setValue(int(self.settings.get("bases_per_line", 60)))
        self.spin_bpl.setFixedWidth(56)
        self.spin_bpl.setToolTip("每行显示的碱基数")
        self.spin_bpl.valueChanged.connect(self._apply_bpl)
        self.spin_bpl.setEnabled(not self.chk_auto.isChecked())
        self.spin_font = QSpinBox()
        self.spin_font.setRange(7, 22)
        self.spin_font.setValue(int(self.settings.get("font_size", 10)))
        self.spin_font.setFixedWidth(48)
        self.spin_font.setToolTip("序列字号")
        self.spin_font.valueChanged.connect(self._apply_font_size)
        lay.addWidget(lab("每行"))
        lay.addWidget(self.spin_bpl)
        lay.addWidget(lab("字号"))
        lay.addWidget(self.spin_font)
        return bar

    # ================================================================== #
    # 动作
    # ================================================================== #
    def _build_actions(self) -> None:
        mbar = self.menuBar()

        m_file = mbar.addMenu("文件(&F)")
        self.act_genome = QAction("打开基因组 FASTA…", self)
        self.act_genome.setShortcut("Ctrl+O")
        self.act_genome.triggered.connect(lambda: self._browse_role("genome"))
        self.act_ann = QAction("打开注释 GFF/GTF…", self)
        self.act_ann.setShortcut("Ctrl+Shift+O")
        self.act_ann.triggered.connect(lambda: self._browse_role("annotation"))
        m_file.addAction(self.act_genome)
        m_file.addAction(self.act_ann)
        m_file.addSeparator()
        self.act_import_archive = QAction("导入压缩包（自动识别文件）…", self)
        self.act_import_archive.setShortcut("Ctrl+Shift+I")
        self.act_import_archive.triggered.connect(lambda: self.import_archive_dialog())
        m_file.addAction(self.act_import_archive)
        self.act_library = QAction("数据集库…", self)
        self.act_library.setShortcut("Ctrl+L")
        self.act_library.triggered.connect(self.open_library)
        m_file.addAction(self.act_library)
        self.act_export_dataset = QAction("导出数据集压缩包（方便下次导入）…", self)
        self.act_export_dataset.setShortcut("Ctrl+Shift+E")
        self.act_export_dataset.triggered.connect(self.export_dataset_dialog)
        m_file.addAction(self.act_export_dataset)
        m_file.addSeparator()
        for role, label in (("mrna", "参考 mRNA FASTA…"), ("cds", "参考 CDS FASTA…"),
                            ("protein", "参考蛋白 FASTA…")):
            act = QAction(label, self)
            act.triggered.connect(lambda _c=False, r=role: self._browse_role(r))
            m_file.addAction(act)
        m_file.addSeparator()
        self.act_reload = QAction("重新加载上次的文件", self)
        self.act_reload.setShortcut("Ctrl+R")
        self.act_reload.triggered.connect(self._reload_saved)
        m_file.addAction(self.act_reload)
        self.act_auto_restore = QAction("启动时自动加载上次的文件", self)
        self.act_auto_restore.setCheckable(True)
        self.act_auto_restore.setChecked(self._auto_restore_enabled())
        self.act_auto_restore.setToolTip(
            "打开时自动把上次使用的 5 个文件加载好；关闭后只填回路径，需要时按 Ctrl+R。")
        self.act_auto_restore.toggled.connect(self._toggle_auto_restore)
        m_file.addAction(self.act_auto_restore)
        m_file.addSeparator()
        self.act_export_gene = QAction("导出当前基因全部序列…", self)
        self.act_export_gene.triggered.connect(self.export_current_gene)
        self.act_export_doc = QAction("导出彩色序列文档（Word / 网页）…", self)
        self.act_export_doc.setShortcut("Ctrl+Shift+D")
        self.act_export_doc.triggered.connect(lambda: self.export_sequence_document_dialog())
        self.act_batch = QAction("批量导出多个基因…", self)
        self.act_batch.triggered.connect(self.batch_export_dialog)
        self.act_export_report = QAction("导出校验报告…", self)
        self.act_export_report.triggered.connect(self.export_report)
        m_file.addAction(self.act_export_gene)
        m_file.addAction(self.act_export_doc)
        m_file.addAction(self.act_batch)
        m_file.addAction(self.act_export_report)
        m_file.addSeparator()
        act_quit = QAction("退出", self)
        act_quit.setShortcut("Ctrl+Q")
        act_quit.triggered.connect(self.close)
        m_file.addAction(act_quit)

        m_edit = mbar.addMenu("编辑(&E)")
        self.act_copy_panel = QAction("复制当前面板全部序列", self)
        self.act_copy_panel.setShortcut("Ctrl+Shift+C")
        self.act_copy_panel.triggered.connect(self._copy_focused_panel)
        m_edit.addAction(self.act_copy_panel)
        act_find = QAction("查找序列…", self)
        act_find.setShortcut("Ctrl+F")
        act_find.triggered.connect(self._focus_find)
        m_edit.addAction(act_find)
        act_next = QAction("下一个匹配", self)
        act_next.setShortcut("F3")
        act_next.triggered.connect(lambda: self._step_find(1))
        m_edit.addAction(act_next)
        act_prev = QAction("上一个匹配", self)
        act_prev.setShortcut("Shift+F3")
        act_prev.triggered.connect(lambda: self._step_find(-1))
        m_edit.addAction(act_prev)
        self.act_copy_info = QAction("复制基因基本信息", self)
        self.act_copy_info.triggered.connect(self._copy_gene_info)
        m_edit.addAction(self.act_copy_info)

        m_view = mbar.addMenu("视图(&V)")
        m_panels = m_view.addMenu("显示面板")
        self._panel_actions = {}
        for key in list(self.panels) + ["validation"]:
            act = QAction(PANEL_TITLES.get(key, "序列校验"), self)
            act.setCheckable(True)
            act.setChecked(True)
            act.toggled.connect(lambda on, k=key: self._set_panel_visible(k, on))
            m_panels.addAction(act)
            self._panel_actions[key] = act
        act_colors = QAction("配色设置…", self)
        act_colors.triggered.connect(self._edit_colors)
        m_view.addAction(act_colors)
        act_reset = QAction("恢复默认布局", self)
        act_reset.triggered.connect(self._reset_layout)
        m_view.addAction(act_reset)

        m_tools = mbar.addMenu("工具(&T)")
        act_re = QAction("重新提取当前转录本", self)
        act_re.setShortcut("F5")
        act_re.triggered.connect(self.extract_current)
        m_tools.addAction(act_re)
        act_val = QAction("重新校验参考序列", self)
        act_val.triggered.connect(self.revalidate)
        m_tools.addAction(act_val)
        act_struct = QAction("导出基因结构图…", self)
        act_struct.triggered.connect(self.export_structure)
        m_tools.addAction(act_struct)

        m_help = mbar.addMenu("帮助(&H)")
        act_demo = QAction("加载示例数据", self)
        act_demo.triggered.connect(self.load_demo)
        m_help.addAction(act_demo)
        act_help = QAction("使用说明", self)
        act_help.triggered.connect(self._show_help)
        m_help.addAction(act_help)
        act_about = QAction("关于", self)
        act_about.triggered.connect(self._show_about)
        m_help.addAction(act_about)

        QShortcut(QKeySequence("F3"), self, activated=lambda: self._step_find(1))

    # ================================================================== #
    # 设置
    # ================================================================== #
    def _restore_settings(self) -> None:
        for role in ("genome", "annotation", "mrna", "cds", "protein"):
            path = self.settings.get(f"{role}_path", "")
            if path and os.path.exists(path):
                self.file_panel.rows[role].set_state(path, "上次使用（未加载）", "idle")
        geo = self.settings.get("window_geometry", "")
        if geo:
            try:
                self.restoreGeometry(QByteArray.fromBase64(geo.encode()))
            except Exception:
                pass
        state = self.settings.get("splitter_states", {})
        if isinstance(state, dict):
            self._restore_splitter("vsplit", self.vsplit, state)
            self._restore_splitter("panels", self.panels_splitter, state)
            for i, col in enumerate(self._col_splitters):
                self._restore_splitter(f"col{i}", col, state)
        # 方向是布局的硬约束：任何恢复动作之后都重申一次，防止状态串用导致串栏
        self.panels_splitter.setOrientation(Qt.Orientation.Horizontal)
        for col in self._col_splitters:
            col.setOrientation(Qt.Orientation.Vertical)
        self._apply_font_size(int(self.settings.get("font_size", 10)))
        self._apply_bpl(int(self.settings.get("bases_per_line", 60)))
        self._apply_auto_fit(bool(self.settings.get("auto_fit", True)))
        ds_id = int(self.settings.get("dataset_id", 0) or 0)
        if ds_id:
            ds = self.library.get(ds_id)
            if ds is not None:
                self._dataset_id = ds_id
                self.current_dataset = ds
                self.file_panel.set_dataset(ds.name, "上次使用，尚未加载")
        if self.settings.existing_path("genome_path"):
            if self._auto_restore_on_start and bool(self.settings.get("auto_restore", True)):
                self.status_msg.setText("正在恢复上次使用的文件…")
            else:
                self.status_msg.setText(
                    "已恢复上次的文件路径。点击「文件 → 重新加载上次的文件」"
                    "(Ctrl+R) 即可加载。")

    # ---------------------------- 启动自动恢复 ---------------------------- #
    def _schedule_auto_restore(self) -> None:
        """窗口建好后再排自动加载，避免界面还没画出来就开始弹进度条。"""
        if not self._auto_restore_on_start:
            return
        if not bool(self.settings.get("auto_restore", True)):
            return
        if not self.settings.existing_path("genome_path"):
            return
        QTimer.singleShot(120, self._auto_restore_files)

    def _auto_restore_files(self) -> None:
        """把上次用过的文件真正加载回来。

        旧版重启后只把路径填回左侧列表（状态写"上次使用（未加载）"），
        用户看到文件名都在、以为已经就绪，一搜索却毫无反应 —— 得自己找到
        「文件 → 重新加载上次的文件」才能用，很容易被当成软件坏了。
        现在窗口一出现就把存在的文件排队加载，重新打开即可直接用。
        """
        if (self.provider is not None or self.model is not None
                or self._job is not None or self._queue):
            return          # 已经在加载（命令行指定文件、示例数据…），不要抢
        roles = [r for r in ("genome", "annotation", "mrna", "cds", "protein")
                 if self.settings.existing_path(f"{r}_path")]
        if not roles:
            return
        self._restoring = True
        for role in roles:
            self.load_role(role, self.settings.get(f"{role}_path", ""))
        self.status_msg.setText(f"正在恢复上次使用的 {len(roles)} 个文件…")

    def _auto_restore_enabled(self) -> bool:
        return bool(self.settings.get("auto_restore", True))

    def _toggle_auto_restore(self, on: bool) -> None:
        self.settings["auto_restore"] = bool(on)
        self.settings.save()
        self.status_msg.setText("已开启：启动时自动加载上次的文件"
                                if on else
                                "已关闭：启动时只恢复路径，需要时按 Ctrl+R 加载")

    def _apply_default_sizes(self) -> None:
        """默认分配：结构图约 1/4 高，4 个序列面板均分，校验区略高一点。"""
        self.vsplit.setSizes([232, 700])
        self.panels_splitter.setSizes([1, 1])
        for col in self._col_splitters:
            n = col.count()
            if n >= 5:      # 右列：4 个序列面板 + 校验区
                col.setSizes([108] * (n - 1) + [230])
            else:
                col.setSizes([150] * max(n, 1))

    def _restore_splitter(self, key: str, sp: QSplitter, store: Dict[str, str]) -> None:
        """按名字恢复单个分隔条的尺寸。

        早期版本只存了 ``vsplit`` 一份状态，却同时恢复到水平分栏上，
        而 Qt6 的 ``restoreState`` 会把方向一并写回，结果两列变成了上下叠放。
        现在每个分隔条独立存一份，并且恢复后再确认一次方向。
        """
        blob = store.get(key, "")
        if not blob:
            return
        try:
            sp.restoreState(QByteArray.fromBase64(blob.encode()))
        except Exception:
            return
        sp.setOrientation(Qt.Orientation.Horizontal if key == "panels"
                          else Qt.Orientation.Vertical)

    def _save_settings(self) -> None:
        s = self.settings
        if self.provider is not None:
            s["genome_path"] = getattr(self.provider, "path", "")
        if self.model is not None:
            s["annotation_path"] = self.model.path
        s["upstream"] = self.spin_up.value()
        s["downstream"] = self.spin_dn.value()
        s["bases_per_line"] = self.spin_bpl.value()
        s["auto_fit"] = self.chk_auto.isChecked()
        s["font_size"] = self.spin_font.value()
        s["palette"] = self.palette_map
        s["dataset_id"] = int(self._dataset_id or 0)
        s["last_gene"] = self.current_gene_id or ""
        try:
            s["window_geometry"] = bytes(self.saveGeometry().toBase64()).decode()
            s["splitter_states"] = {
                "vsplit": bytes(self.vsplit.saveState().toBase64()).decode(),
                "panels": bytes(self.panels_splitter.saveState().toBase64()).decode(),
                **{f"col{i}": bytes(col.saveState().toBase64()).decode()
                   for i, col in enumerate(self._col_splitters)},
            }
        except Exception:
            pass
        s.save()

    def closeEvent(self, event) -> None:  # noqa: N802
        running = self._job is not None and self._job.isRunning()
        if running:
            ret = QMessageBox.question(self, "任务进行中", "有任务正在执行，确定退出吗？")
            if ret != QMessageBox.Yes:
                event.ignore()
                return
        self.shutdown()
        self._save_settings()
        super().closeEvent(event)

    def shutdown(self) -> None:
        """停掉后台任务、关闭文件句柄。

        ``QThread`` 在运行中被销毁会让整个进程崩溃（"Destroyed while thread is
        still running"），所以退出前必须等线程真正结束；只有等不到时才强制
        terminate（仅可能发生在第三方库卡在系统调用里的极端情况）。
        """
        self._queue.clear()
        job, self._job = self._job, None
        if job is not None:
            job.cancel()
            if not job.wait(5000):
                job.terminate()
                job.wait(2000)
        if self.provider is not None:
            try:
                self.provider.close()
            except Exception:
                pass
            self.provider = None

    # ================================================================== #
    # 文件加载
    # ================================================================== #
    ROLE_TITLES = {
        "genome": "基因组 FASTA",
        "annotation": "注释 GFF/GTF",
        "mrna": "参考 mRNA FASTA",
        "cds": "参考 CDS FASTA",
        "protein": "参考蛋白 FASTA",
    }

    def _browse_role(self, role: str) -> None:
        row = self.file_panel.rows[role]
        start_dir = (self.settings.existing_path(f"{role}_path")
                     or self.settings.get("last_dir", ""))
        title = f"选择{self.ROLE_TITLES.get(role, '文件')}"
        if start_dir and os.path.isfile(start_dir):
            start_dir = os.path.dirname(start_dir)
        path, _ = QFileDialog.getOpenFileName(self, title, start_dir, row.filters)
        if not path:
            return
        self.settings["last_dir"] = os.path.dirname(path)
        self.load_role(role, path)

    def load_role(self, role: str, path: str) -> None:
        if role == "genome":
            self._start_job(load_genome, path, busy="正在打开基因组并建立索引…",
                            title="加载基因组", on_done=self._genome_loaded)
        elif role == "annotation":
            self._start_job(load_annotation_file, path, busy="正在解析注释文件…",
                            title="加载注释", on_done=self._annotation_loaded)
        else:
            kind = {"mrna": "mRNA", "cds": "CDS", "protein": "protein"}[role]
            self._start_job(load_reference_file, path, kind, busy=f"正在读取参考 {kind} …",
                            title=f"加载参考 {kind}",
                            on_done=lambda rs, r=role, k=kind: self._reference_loaded(r, k, rs))

    def _genome_loaded(self, info: GenomeInfo) -> None:
        self._cache.clear()
        self.provider = info.provider
        self.genome_info = info
        self.file_panel.rows["genome"].set_state(
            info.path, f"{info.summary()}", "ok")
        self.settings["genome_path"] = info.path
        self.status_msg.setText(f"基因组已加载：{info.summary()}")
        self._update_actions()
        self._auto_extract_if_ready()

    def _annotation_loaded(self, model) -> None:
        self._cache.clear()
        self.model = model
        row = self.file_panel.rows["annotation"]
        level = "ok" if not model.warnings else "warn"
        row.set_state(model.path, model.stat_text(), level)
        self.settings["annotation_path"] = model.path
        self.status_msg.setText(f"注释已加载：{model.stat_text()}")
        if model.warnings:
            if self._restoring:
                # 启动自动恢复时不打断用户：提醒留在状态栏与校验面板里
                self.status_msg.setText(
                    f"注释已加载：{model.stat_text()}"
                    f"（{len(model.warnings)} 条提醒，见「数据校验」）")
            else:
                self._warn("注释解析完成，但有需要留意的地方", "\n".join(model.warnings[:20]))
        self._update_actions()
        # 预填搜索框，降低"不知道输什么"的门槛
        if model.genes:
            first = next(iter(model.genes))
            if not self.search_panel.edit.text():
                self.search_panel.edit.setPlaceholderText(
                    f"例如：{first}  或  {model.genes[first].gene_name}")
        self._auto_extract_if_ready()

    def _reference_loaded(self, role: str, kind: str, rs) -> None:
        self._cache.clear()
        if role == "mrna":
            self.mrna_refs = rs
        elif role == "cds":
            self.cds_refs = rs
        else:
            self.protein_refs = rs
        self.file_panel.rows[role].set_state(rs.path, f"{len(rs)} 条 | {rs.summary()}", "ok")
        self.settings[f"{role}_path"] = rs.path
        self.status_msg.setText(f"参考 {kind} 已加载：{rs.summary()}")
        self._update_actions()
        if self.extraction is not None:
            self.revalidate()

    def _clear_role(self, role: str) -> None:
        self._cache.clear()
        self.file_panel.rows[role].clear()
        if role == "genome":
            self.provider = None
            self.genome_info = None
        elif role == "annotation":
            self.model = None
        elif role == "mrna":
            self.mrna_refs = None
        elif role == "cds":
            self.cds_refs = None
        elif role == "protein":
            self.protein_refs = None
        self.settings[f"{role}_path"] = ""
        self._update_actions()

    def _reload_saved(self) -> None:
        loaded = []
        for role in ("genome", "annotation", "mrna", "cds", "protein"):
            path = self.settings.get(f"{role}_path", "")
            if path and os.path.exists(path) and not self._role_loaded(role):
                self.load_role(role, path)
                loaded.append(role)
        if not loaded:
            self._info("没有需要重新加载的文件",
                       "可用「文件」菜单逐个选择文件，或先把文件路径填进左侧列表。")

    def _role_loaded(self, role: str) -> bool:
        return {
            "genome": self.provider is not None,
            "annotation": self.model is not None,
            "mrna": self.mrna_refs is not None,
            "cds": self.cds_refs is not None,
            "protein": self.protein_refs is not None,
        }.get(role, False)

    def load_demo(self) -> None:
        base = _demo_dir()
        if not base:
            self._error("找不到示例数据",
                        "示例数据应位于程序目录的 examples/ 下；"
                        "也可运行 python examples/make_example_data.py 重新生成。")
            return
        self.load_role("genome", os.path.join(base, "demo_genome.fa"))
        self.load_role("annotation", os.path.join(base, "demo_annotation.gff3"))
        self.load_role("mrna", os.path.join(base, "demo_mrna.fa"))
        self.load_role("cds", os.path.join(base, "demo_cds.fa"))
        self.load_role("protein", os.path.join(base, "demo_protein.fa"))
        self.search_panel.edit.setText("GeneOne")
        # 不能写死一个延时就去搜索：5 个文件是排队异步加载的，
        # 加载没完时搜到一半会弹「还没有加载注释文件」。这里改成"加载完再搜"。
        QTimer.singleShot(200, self._demo_autosearch)

    def _demo_autosearch(self, tries: int = 100) -> None:
        """等示例数据全部加载完再自动搜索（每 200 ms 检查一次，最多 20 s）。"""
        if self._job is not None or self._queue:
            if tries > 0:
                QTimer.singleShot(200, lambda: self._demo_autosearch(tries - 1))
            return
        if self.model is not None:
            self.do_search("GeneOne")
        # 加载失败时不再追一次搜索：真正的错误已经由加载流程报过了

    # ================================================================== #
    # 压缩包导入 / 数据集库
    # ================================================================== #
    ARCHIVE_FILTER = ("压缩包 (*.zip *.tar *.tar.gz *.tgz *.tar.bz2 *.tbz2 "
                      "*.tar.xz *.txz);;所有文件 (*)")

    def import_archive_dialog(self, path: str = "") -> None:
        """选一个压缩包 → 自动识别 → 人工确认 → 解包入库 → 加载。"""
        if not path:
            start = (self.settings.get("last_archive_dir", "")
                     or self.settings.get("last_dir", ""))
            path, _ = QFileDialog.getOpenFileName(
                self, "选择数据压缩包", start, self.ARCHIVE_FILTER)
        if not path:
            return
        self.settings["last_archive_dir"] = os.path.dirname(path)
        scan = self._scan_archive_gui(path)
        if scan is None:
            return
        dlg = BundleImportDialog(scan, self)
        if dlg.exec() != QDialog.DialogCode.Accepted:
            self.status_msg.setText("已取消导入（没有解包任何文件）")
            return
        roles = dlg.roles()
        name = dlg.dataset_name()
        self.status_msg.setText(f"正在解包「{name}」…")
        self._start_job(
            import_bundle_from_archive, scan, roles, name, self.library.path,
            busy=f"正在解包「{name}」…", title="导入压缩包",
            on_done=self._bundle_imported, show_percent=True,
        )

    def _scan_archive_gui(self, path: str):
        """带进度条与取消按钮地扫描压缩包（只读头部，不解包）。"""
        from PySide6.QtWidgets import QProgressDialog

        box = QProgressDialog(f"正在扫描 {os.path.basename(path)} …", "取消",
                              0, 100, self)
        box.setWindowTitle("扫描压缩包")
        box.setWindowModality(Qt.WindowModal)
        box.setMinimumDuration(0)
        box.setAutoClose(False)
        box.setAutoReset(False)
        box.setValue(0)
        state = {"cancel": False}

        def cb(pct: int, msg: str) -> None:
            box.setValue(max(0, min(100, int(pct))))
            if msg:
                box.setLabelText(msg)
            QApplication.processEvents()
            if box.wasCanceled():
                state["cancel"] = True

        try:
            return scan_archive(path, progress_cb=cb,
                                cancel_check=lambda: state["cancel"])
        except CancelledError:
            self.status_msg.setText("已取消扫描")
            return None
        except GGEError as exc:
            self._error(exc.message, exc.hint or "")
            return None
        except Exception as exc:               # pragma: no cover - 兜底
            self._error(f"扫描压缩包时出错：{exc}",
                        "确认压缩包没有损坏；也可以先手动解压，再逐个选择文件。")
            return None
        finally:
            box.close()

    def _bundle_imported(self, result: BundleImportResult) -> None:
        ds = result.dataset
        if ds is None:
            self._error("导入结果为空", "请重试；若反复出现请把压缩包提供给开发者。")
            return
        self._apply_dataset(ds)
        msg = f"已导入「{ds.name}」：{ds.summary()}"
        if result.replaced:
            msg += "（替换了同一压缩包的旧记录）"
        self.status_msg.setText(msg)
        notes = list(getattr(result.import_result, "warnings", []) or [])
        if notes:
            self._warn("导入完成，但有以下提示", "\n".join(notes))
        elif self._dataset_id:
            self.file_panel.set_dataset(ds.name, ds.summary())

    def open_library(self) -> None:
        dlg = LibraryDialog(self.library, self, current_id=self._dataset_id)
        dlg.exec()
        if dlg.wants_import:
            self.import_archive_dialog()
            return
        ds = dlg.selected_dataset()
        if ds is not None and dlg.result() == QDialog.DialogCode.Accepted:
            self.load_dataset(ds)

    # ------------------------------------------------------------------ #
    # 反向：把"分别导入的几个文件"打包成一个可再次导入的数据集
    # ------------------------------------------------------------------ #
    def current_role_paths(self) -> Dict[str, str]:
        """当前真正加载着的文件（没加载的角色不猜、不塞）。"""
        out: Dict[str, str] = {}
        if self.genome_info is not None:
            out["genome"] = self.genome_info.path
        if self.model is not None:
            out["annotation"] = self.model.path
        for role, refs in (("mrna", self.mrna_refs), ("cds", self.cds_refs),
                           ("protein", self.protein_refs)):
            path = getattr(refs, "path", "") if refs is not None else ""
            if path:
                out[role] = path
        return {k: v for k, v in out.items() if v and os.path.isfile(v)}

    def current_stats(self) -> Dict[str, int]:
        stats: Dict[str, int] = {}
        if self.model is not None:
            stats["n_genes"] = getattr(self.model, "n_genes", 0)
            stats["n_transcripts"] = getattr(self.model, "n_transcripts", 0)
            stats["n_seqids"] = len(getattr(self.model, "seqids", []) or [])
        if self.genome_info is not None:
            stats["total_bases"] = int(self.genome_info.total_bases)
        return stats

    def export_dataset_dialog(self) -> None:
        paths = self.current_role_paths()
        if not paths.get("genome") or not paths.get("annotation"):
            self._error("还不能打包",
                        "至少要先加载基因组 FASTA 与注释 GFF/GTF；"
                        "参考 mRNA/CDS/蛋白 可选。")
            return
        name = (self.current_dataset.name if self.current_dataset else "") or \
            os.path.splitext(os.path.basename(paths["genome"]))[0]
        base = self.default_export_dir()
        default = os.path.join(base, safe_filename(name) + ".zip")
        dlg = DatasetExportDialog(paths, name, "", default, self)
        if dlg.exec() != QDialog.DialogCode.Accepted:
            return
        vals = dlg.values()
        self._start_job(
            export_dataset_archive, vals["path"], vals["paths"], vals["name"],
            vals["note"], "", self.current_stats(),
            busy=f"正在打包「{vals['name']}」…", title="导出数据集压缩包",
            on_done=self._dataset_exported, show_percent=True,
        )

    def _dataset_exported(self, result: DatasetExportResult) -> None:
        self.settings["export_dir"] = os.path.dirname(result.path)
        self.settings["last_archive_dir"] = os.path.dirname(result.path)
        self.status_msg.setText(f"数据集压缩包已导出：{result.path}")
        self._info(
            "数据集压缩包已导出",
            f"{result.summary()}\n\n{result.path}\n\n"
            "包里带一份 gge_dataset.json 清单，记录每个文件是哪个角色。"
            "下次用「文件 → 导入压缩包」选中它，会按清单直接认回来"
            "（不再靠内容猜，也跳过全量注释统计，打开更快）。")

    def load_dataset(self, ds: Dataset) -> None:
        """从库里载入一套数据：先把不属于它的旧文件清掉，再逐个加载。"""
        miss = ds.missing_required()
        if miss:
            self._error("这套数据的必需文件不见了",
                        "\n".join(ds.path(r) or ROLE_LABELS[r] for r in miss)
                        + "\n\n文件可能被移动或删除，请重新导入压缩包。")
            return
        self.library.touch(ds.id)
        self._apply_dataset(ds)
        self.status_msg.setText(f"正在加载数据集「{ds.name}」（{ds.summary()}）…")

    def _apply_dataset(self, ds: Dataset) -> None:
        self.current_dataset = ds
        self._dataset_id = ds.id
        self.settings["dataset_id"] = ds.id
        self.file_panel.set_dataset(ds.name, ds.summary())
        paths = ds.paths()
        self._reset_view_for_new_data()
        self._clear_missing_roles(set(paths))
        for role in ROLE_ORDER:
            p = paths.get(role, "")
            if p and os.path.exists(p):
                self.load_role(role, p)
            elif p:
                self.status_msg.setText(f"{ROLE_LABELS[role]}文件不存在：{p}")

    def _reset_view_for_new_data(self) -> None:
        """换数据前把上一条基因的结果清干净，避免把两个数据源的结果混在一起。"""
        self.extraction = None
        self.validation = {}
        self.current_gene_id = ""
        self.current_tx_id = ""
        self._manual = {"mRNA": "", "CDS": "", "protein": ""}
        self._cache.clear()
        for p in self.panels.values():
            p.clear()
        self.structure.set_extraction(None)
        self.info_panel.clear()
        self.search_panel.set_hits([])
        self.search_panel.edit.clear()
        self.validation_panel.set_results({}, {"mRNA": [], "CDS": [], "protein": []}, False)
        self.setWindowTitle(f"{APP_NAME} · {APP_NAME_CN}")
        self._update_actions()

    def _clear_missing_roles(self, keep: set) -> None:
        for role in ROLE_ORDER:
            if role not in keep and self._role_loaded(role):
                self._clear_role(role)

    # ================================================================== #
    # 搜索与选择
    # ================================================================== #
    def do_search(self, text: str = "") -> None:
        text = text or self.search_panel.edit.text().strip()
        if not text:
            return
        if self.model is None:
            self._error("还没有加载注释文件", "请先选择 GFF3/GTF 注释文件。")
            return
        hits: List[SearchHit] = search_genes(self.model, text)
        self.search_panel.set_hits(hits)
        if not hits:
            self.status_msg.setText(f"没有找到与 “{text}” 匹配的基因")
            return
        self.status_msg.setText(f"找到 {len(hits)} 个基因，选择一条查看。")
        if len(hits) >= 1:
            # 高置信度（完全一致）时自动选中，其余让用户自己点
            if hits[0].score >= 96:
                self.select_gene(hits[0].gene_id)
            else:
                self.search_panel.table.selectRow(0)

    def select_gene(self, gene_id: str) -> None:
        if self.model is None:
            return
        gene = self.model.get_gene(gene_id)
        if gene is None:
            self._error(f"注释中找不到基因 {gene_id}")
            return
        self.current_gene_id = gene_id
        canon = gene.canonical_transcript()
        self.info_panel.set_gene(gene, canon.tx_id if canon else "")
        self.settings.add_recent_gene(gene_id)
        self.current_tx_id = self.info_panel.current_tx()
        if not gene.transcripts:
            self._error(f"基因 {gene_id} 没有任何转录本",
                        "该注释记录只有 gene 行，缺少 exon/CDS，无法提取序列。")
            return
        self.extract_current()

    def _on_transcript_changed(self, tx_id: str) -> None:
        self.current_tx_id = tx_id
        self.extract_current()

    def _auto_extract_if_ready(self) -> None:
        if self.provider is None or self.model is None:
            return
        if self.current_gene_id:
            return
        if self._restoring:
            # 自动恢复时顺手把上次看的基因也打开，重新启动就回到原处
            gene = str(self.settings.get("last_gene", "") or "")
            if gene and gene in self.model.genes:
                self.search_panel.edit.setText(gene)
                self.select_gene(gene)
                return
        self.status_msg.setText("文件就绪，请输入基因名或 ID 开始搜索。")

    # ================================================================== #
    # 提取
    # ================================================================== #
    def _cache_key(self) -> tuple:
        """缓存键：基因/转录本/上下游 + 全部输入文件 + 手动指定的参考 ID。

        任一环节变化（换了文件、改了长度、手动指定了别的参考序列）都视为不同结果，
        避免把旧结果当成新结果展示。
        """
        def stamp(path: str) -> tuple:
            try:
                st = os.stat(path)
                return (path, int(st.st_mtime), st.st_size)
            except OSError:
                return (path, 0, 0)

        return (
            self.current_gene_id, self.current_tx_id,
            int(self.spin_up.value()), int(self.spin_dn.value()),
            stamp(getattr(self.provider, "path", "")),
            stamp(getattr(self.model, "path", "")),
            stamp(getattr(self.mrna_refs, "path", "")),
            stamp(getattr(self.cds_refs, "path", "")),
            stamp(getattr(self.protein_refs, "path", "")),
            self._manual["mRNA"], self._manual["CDS"], self._manual["protein"],
            bool(self.mrna_refs is not None), bool(self.cds_refs is not None),
            bool(self.protein_refs is not None),
        )

    def extract_current(self) -> None:
        if self.provider is None:
            self._error("还没有加载基因组 FASTA", "请先选择基因组序列文件。")
            return
        if self.model is None:
            self._error("还没有加载注释文件", "请先选择 GFF3/GTF 文件。")
            return
        if not self.current_gene_id:
            gid = self.search_panel.current_gene()
            if not gid:
                self._error("请先搜索并选择一个基因")
                return
            self.current_gene_id = gid

        key = self._cache_key()
        cached = self._cache.get(key)
        if cached is not None:
            # 参数没变就直接复用，不重跑 IO 与比对
            self._cache.move_to_end(key)
            self._extraction_done(cached)
            self.status_msg.setText(
                self.status_msg.text() + "　（命中缓存，未重新计算）")
            return

        self._start_job(
            extract_and_validate,
            self.provider, self.model, self.current_gene_id,
            self.current_tx_id,
            self.spin_up.value(), self.spin_dn.value(),
            self.mrna_refs, self.cds_refs, self.protein_refs,
            self._manual["mRNA"], self._manual["CDS"], self._manual["protein"],
            busy=f"正在提取 {self.current_gene_id} …",
            title="提取序列",
            on_done=lambda bundle, k=key: self._store_and_show(k, bundle),
        )

    def _store_and_show(self, key: tuple, bundle: ExtractionBundle) -> None:
        """缓存最近若干次提取结果（顺序移动 + FIFO 淘汰），再渲染。"""
        self._cache[key] = bundle
        self._cache.move_to_end(key)
        while len(self._cache) > self._cache_limit:
            self._cache.popitem(last=False)
        self._extraction_done(bundle)

    def revalidate(self) -> None:
        if self.extraction is None:
            self._info("还没有可校验的结果", "请先选择一个基因并提取。")
            return
        self.extract_current()

    def _extraction_done(self, bundle: ExtractionBundle) -> None:
        ex = bundle.extraction
        self.extraction = ex
        self.validation = bundle.validation or {}
        self.current_gene_id = ex.gene_id
        self.current_tx_id = ex.tx_id

        self.structure.set_extraction(ex)
        for key, panel in self.panels.items():
            panel.set_panel(ex.panels.get(key))
        self.info_panel.set_checks(ex.checks, ex.warnings)
        self._refresh_validation()
        self.status_msg.setText(ex.summary_line())
        self.setWindowTitle(
            f"{APP_NAME} · {ex.gene_name} ({ex.gene_id}) / {ex.tx_id}")
        if ex.warnings:
            self.status_msg.setToolTip("\n".join(ex.warnings))
        self._update_actions()

    def _refresh_validation(self) -> None:
        ref_ids = {
            "mRNA": [r.rid for r in self.mrna_refs.records] if self.mrna_refs else [],
            "CDS": [r.rid for r in self.cds_refs.records] if self.cds_refs else [],
            "protein": [r.rid for r in self.protein_refs.records] if self.protein_refs else [],
        }
        has_refs = bool(self.mrna_refs or self.cds_refs or self.protein_refs)
        self.validation_panel.set_results(self.validation, ref_ids, has_refs)

    def _on_manual_ref(self, kind: str, ref_id: str) -> None:
        self._manual[kind] = ref_id
        if self.extraction is not None:
            self.extract_current()

    def _show_diff_dialog(self, kind: str) -> None:
        res = self.validation.get(kind)
        if res is None:
            return
        dlg = DiffDialog(res, self.extraction, self)
        dlg.exec()

    # ================================================================== #
    # 面板操作
    # ================================================================== #
    def _apply_auto_fit(self, enabled: bool) -> None:
        for p in self.panels.values():
            p.set_auto_fit(bool(enabled))
        self.spin_bpl.setEnabled(not enabled)
        if not enabled:
            for p in self.panels.values():
                p.set_bases_per_line(self.spin_bpl.value())
        self.settings["auto_fit"] = bool(enabled)

    def _apply_bpl(self, n: int) -> None:
        for p in self.panels.values():
            p.set_bases_per_line(int(n))
        if hasattr(self, "spin_bpl") and self.spin_bpl.value() != int(n):
            self.spin_bpl.setValue(int(n))

    def _apply_font_size(self, pt: int) -> None:
        for p in self.panels.values():
            p.set_font_size(int(pt))
        if hasattr(self, "spin_font") and self.spin_font.value() != int(pt):
            self.spin_font.setValue(int(pt))

    def _on_cursor_info(self, text: str) -> None:
        self.status_msg.setText(text)

    def _focused_panel(self) -> Optional[SequencePanel]:
        w = QApplication.focusWidget()
        while w is not None:
            if isinstance(w, SequencePanel):
                return w
            w = w.parentWidget()
        # 没有焦点时用基因组面板
        return self.panels.get("genomic")

    def _copy_focused_panel(self) -> None:
        panel = self._focused_panel()
        if panel is None:
            return
        text = panel.copy_plain()
        self.status_msg.setText(f"已复制 {len(text):,} bp 到剪贴板")

    def _copy_gene_info(self) -> None:
        if self.extraction is None:
            return
        QGuiApplication.clipboard().setText(gene_basic_info_text(self.extraction))
        self.status_msg.setText("已复制基因基本信息")

    def copy_panel_fasta(self, key: str) -> None:
        if self.extraction is None:
            return
        text = panel_fasta_text(self.extraction, key, per_block=True)
        if not text:
            self._info("该面板没有可导出的序列")
            return
        QGuiApplication.clipboard().setText(text)
        self.status_msg.setText(f"已复制 {PANEL_TITLES.get(key, key)} 的 FASTA 文本")

    def export_panel(self, key: str) -> None:
        if self.extraction is None:
            return
        panel = self.extraction.panels.get(key)
        if not panel or panel.is_empty():
            self._info("该面板没有可导出的序列")
            return
        default = safe_filename(
            f"{self.extraction.gene_id}_{self.extraction.tx_id}_{key}.fasta")
        base = self.default_export_dir()
        path, _ = QFileDialog.getSaveFileName(self, "导出面板序列",
                                              os.path.join(base, default),
                                              "FASTA (*.fa *.fasta);;所有文件 (*)")
        if not path:
            return
        write_text(path, panel_fasta_text(self.extraction, key, per_block=True))
        self.settings["export_dir"] = os.path.dirname(path)
        self.status_msg.setText(f"已导出：{path}")

    def default_export_dir(self) -> str:
        """导出的默认落点：上次用过的目录 → 文档目录 → 用户主目录。

        **不能返回空串**：``os.path.join("", name)`` 得到的是相对路径，exe 双击启动时
        相对的是程序目录（安装位置往往只读），用户就会"点了导出、却不知道文件去哪了"。
        """
        for key in ("export_dir", "last_dir"):
            d = str(self.settings.get(key) or "").strip()
            if d and os.path.isdir(d) and os.access(d, os.W_OK):
                return d
        docs = QStandardPaths.writableLocation(QStandardPaths.DocumentsLocation)
        for d in (docs, os.path.expanduser("~")):
            if d and os.path.isdir(d):
                return d
        return os.getcwd()

    def export_current_gene(self) -> None:
        if self.extraction is None:
            self._info("还没有提取结果", "请先选择一个基因。")
            return
        base = self.default_export_dir()
        outdir = QFileDialog.getExistingDirectory(self, "选择导出目录", base)
        if not outdir:
            return
        try:
            files = export_gene(self.extraction, outdir, self.validation)
        except Exception as exc:  # pragma: no cover
            self._error(f"导出失败：{exc}")
            return
        self.settings["export_dir"] = outdir
        self._info("导出完成", "\n".join(files))

    # ------------------------------------------------------------------ #
    # 彩色序列文档（Word / 网页）
    # ------------------------------------------------------------------ #
    def export_sequence_document_dialog(self, focused_key: str = "") -> None:
        """把面板里的彩色序列导出成 Word / 网页（颜色与界面同源）。"""
        if self.extraction is None:
            self._info("还没有提取结果", "请先选择一个基因，再导出彩色序列文档。")
            return
        preset: List[str] = []
        if focused_key and focused_key in self.extraction.panels:
            preset = [focused_key]
        else:
            panel = self._focused_panel()
            if panel is not None:
                preset = [panel.key]
        if not preset:
            preset = ["genomic"]
        for key in ("mrna", "cds"):
            p = self.extraction.panels.get(key)
            if p is not None and not p.is_empty() and key not in preset:
                preset.append(key)
        base = self.default_export_dir()
        default = os.path.join(base, safe_filename(
            f"{self.extraction.gene_id}_{self.extraction.tx_id}_sequence") + ".docx")
        bpl = int(self.spin_bpl.value()) if hasattr(self, "spin_bpl") else 60
        font = float(self.settings.get("font_size", 9)) * 0.75
        dlg = SequenceDocDialog(self.extraction, preset, default, bpl,
                                max(6.0, min(14.0, font)), self.palette_map,
                                bool(self.settings.get("doc_include_coord", False)), self)
        if dlg.exec() != QDialog.DialogCode.Accepted:
            return
        vals = dlg.values()
        self.settings["doc_include_coord"] = bool(vals["include_coord"])
        self._start_job(
            export_sequence_document, vals["path"], self.extraction, vals["panels"],
            vals["mode"], vals["bpl"], vals["font_pt"], self.palette_map,
            vals["include_header"], vals["include_coord"],
            busy="正在生成彩色序列文档…", title="导出彩色序列文档",
            on_done=self._document_exported, show_percent=True,
        )

    def _document_exported(self, result: DocumentExportResult) -> None:
        if result.files:
            self.settings["export_dir"] = os.path.dirname(result.files[0])
        self.status_msg.setText(result.summary())
        self._info("彩色序列文档已导出", "\n".join(result.files) +
                   "\n\nWord 文档里的底色是文字级底纹，可以直接编辑；"
                   "想贴到已有文档里，也可以在序列面板上右键「复制为彩色序列」。")

    def copy_panel_rich(self, key: str) -> None:
        """把某个面板按彩色 HTML 放进剪贴板：粘到 Word / WPS 就是带底色的表。"""
        if self.extraction is None:
            return
        panel = self.extraction.panels.get(key)
        if panel is None or panel.is_empty():
            self._info("该面板没有可复制的序列")
            return
        bpl = int(self.spin_bpl.value()) if hasattr(self, "spin_bpl") else 60
        html_text = clipboard_html(self.extraction, [key], mode=MODE_REGION, bpl=bpl,
                                   palette=self.palette_map,
                                   include_coord=bool(
                                       self.settings.get("doc_include_coord", False)))
        mime = QMimeData()
        mime.setHtml(html_text)
        mime.setText(panel.sequence)
        QGuiApplication.clipboard().setMimeData(mime)
        self.status_msg.setText(
            f"已复制 {PANEL_TITLES.get(key, key)} 的彩色序列（{panel.length:,} bp），"
            f"可直接粘贴到 Word")

    def export_report(self) -> None:
        if self.extraction is None or not self.validation:
            self._info("还没有校验结果", "请先加载参考 mRNA/CDS 并完成校验。")
            return
        default = safe_filename(f"{self.extraction.gene_id}_{self.extraction.tx_id}_validation.txt")
        base = self.default_export_dir()
        path, _ = QFileDialog.getSaveFileName(self, "导出校验报告",
                                              os.path.join(base, default),
                                              "文本 (*.txt);;Markdown (*.md)")
        if not path:
            return
        write_text(path, validation_report_text(self.extraction, self.validation))
        self.status_msg.setText(f"校验报告已导出：{path}")

    def export_structure(self) -> None:
        if self.extraction is None:
            self._info("还没有结构图", "请先选择一个基因。")
            return
        default = safe_filename(f"{self.extraction.gene_id}_{self.extraction.tx_id}_structure.png")
        base = self.default_export_dir()
        path, flt = QFileDialog.getSaveFileName(self, "导出基因结构图",
                                               os.path.join(base, default),
                                               "PNG 图片 (*.png);;SVG 矢量图 (*.svg)")
        if not path:
            return
        if flt.startswith("SVG") and not path.lower().endswith(".svg"):
            path += ".svg"
        ok = self.structure.export(path)
        if ok:
            self.settings["export_dir"] = os.path.dirname(path)
            self.status_msg.setText(f"结构图已导出：{path}")
        else:
            self._error("导出结构图失败", "请检查目标目录是否可写，或换一个文件名。")

    def batch_export_dialog(self) -> None:
        if self.provider is None or self.model is None:
            self._error("请先加载基因组与注释文件")
            return
        preset: List[str] = []
        if self.extraction is not None:
            preset = [self.extraction.gene_id]
        hits = search_genes(self.model, self.search_panel.edit.text().strip()) \
            if self.search_panel.edit.text().strip() else []
        if hits:
            preset = [h.gene_id for h in hits]
        dlg = BatchExportDialog(preset, self.settings.get("export_dir", "") or
                                self.settings.get("last_dir", ""),
                                self.spin_up.value(), self.spin_dn.value(), self)
        if dlg.exec() != QDialog.DialogCode.Accepted:
            return
        vals = dlg.values()
        if not vals["gene_ids"]:
            self._error("没有指定基因", "在列表里每行填一个 Gene ID。")
            return
        if not vals["outdir"]:
            self._error("没有指定输出目录")
            return
        self._start_job(
            batch_export, self.provider, self.model, vals["gene_ids"], vals["outdir"],
            vals["upstream"], vals["downstream"], vals["all_transcripts"],
            vals["write_table"], self.mrna_refs, self.cds_refs,
            busy="正在批量导出…", title="批量导出", on_done=self._batch_done,
            show_percent=True,
        )

    def _batch_done(self, result: BatchResult) -> None:
        self.settings["export_dir"] = self.settings.get("export_dir", "")
        msg = result.summary()
        if result.problems:
            detail = "\n".join(f"{gid}: {m}" for gid, m in result.problems[:15])
            self._warn(f"批量导出完成（{msg}）", detail)
        else:
            self._info("批量导出完成", msg)
        self.status_msg.setText(f"批量导出：{msg}")

    # ================================================================== #
    # 查找
    # ================================================================== #
    def _focus_find(self) -> None:
        self.find_edit.setFocus(Qt.ShortcutFocusReason)
        self.find_edit.selectAll()

    def _do_find(self) -> None:
        text = self.find_edit.text().strip()
        if not text:
            for p in self.panels.values():
                p.clear_search()
            self.find_info.setText("")
            return
        order: List[str] = []
        total = 0
        for key, panel in self.panels.items():
            n = panel.find(text)
            if n:
                order.append(key)
                total += n
        self._search_order = order
        self._search_idx = 0
        if not order:
            self.find_info.setText("未找到")
            self.status_msg.setText(f"所有面板中都没有找到 “{text}”")
            return
        self.find_info.setText(f"{len(order)} 个面板 / {total} 处")
        self.status_msg.setText(
            f"在 {len(order)} 个面板中找到 {total} 处匹配；F3 跳下一个")
        self.panels[order[0]].canvas.goto_block(0)

    def _step_find(self, direction: int) -> None:
        if not self._search_order:
            if self.find_edit.text().strip():
                self._do_find()
            return
        key = self._search_order[self._search_idx % len(self._search_order)]
        panel = self.panels[key]
        if self._search_idx > 0 or direction < 0:
            panel.next_hit() if direction > 0 else panel.prev_hit()
        self._search_idx = (self._search_idx + (1 if direction > 0 else -1)) % len(self._search_order)
        self.status_msg.setText(f"匹配位于：{PANEL_TITLES.get(key, key)}")

    # ================================================================== #
    # 视图
    # ================================================================== #
    def _set_panel_visible(self, key: str, visible: bool) -> None:
        if key == "validation":
            self.validation_panel.setVisible(visible)
            return
        panel = self.panels.get(key)
        if panel is not None:
            panel.setVisible(visible)

    def _toggle_maximize(self, key: str) -> None:
        panel = self.panels.get(key)
        if panel is None:
            return
        if self._maximized == key:
            for k, p in self.panels.items():
                p.setVisible(self._panel_actions[k].isChecked())
            for i, sp in enumerate(self._col_splitters):
                sp.setVisible(True)
            self.validation_panel.setVisible(self._panel_actions["validation"].isChecked())
            self._maximized = None
            self._apply_default_sizes()
            return
        col_idx = 0 if key in PANEL_LAYOUT[0] else 1
        for i, sp in enumerate(self._col_splitters):
            sp.setVisible(i == col_idx)
        for k, p in self.panels.items():
            p.setVisible(k == key)
        self.validation_panel.setVisible(False)
        self._maximized = key

    def _edit_colors(self) -> None:
        dlg = ColorSettingsDialog(self.palette_map, self)
        if dlg.exec() != QDialog.DialogCode.Accepted:
            return
        self.palette_map = dlg.palette()
        self.settings.set_palette(self.palette_map)
        for p in self.panels.values():
            p.set_palette_map(self.palette_map, UI_COLORS)
        self.structure.set_palette_map(self.palette_map)
        self.status_msg.setText("配色已更新")

    def _reset_layout(self) -> None:
        for k, act in self._panel_actions.items():
            act.setChecked(True)
        if self._maximized:
            self._toggle_maximize(self._maximized)
        self._apply_default_sizes()

    # ================================================================== #
    # 任务与提示
    # ================================================================== #
    def _start_job(self, fn, *args, busy: str = "", title: str = "",
                   on_done=None, show_percent: bool = False, **kwargs) -> None:
        params = {
            "fn": fn, "args": args, "busy": busy, "title": title,
            "on_done": on_done, "show_percent": show_percent, "kwargs": kwargs,
        }
        if self._job is not None:
            # 允许连续请求（例如"加载示例数据"要连开 5 个文件），排队执行
            self._queue.append(params)
            self.status_msg.setText(f"已排队：{title or '任务'}（前面还有 {len(self._queue)} 个）")
            return
        self._launch_job(params)

    def _launch_job(self, params: dict) -> None:
        job = Job(params["fn"], *params["args"], **params["kwargs"])
        self._job = job
        self._job_done = params["on_done"]
        self._job_title = params["title"]
        self._show_percent = params["show_percent"]
        job.progressed.connect(self._on_progress)
        job.succeeded.connect(self._on_job_success)
        job.failed.connect(self._on_job_failed)
        job.wasCancelled.connect(self._on_job_cancelled)
        job.finished.connect(self._on_job_finished)
        self.progress.setRange(0, 0 if not params["show_percent"] else 100)
        self.progress.setValue(0)
        self.progress.show()
        self.btn_cancel.show()
        self.status_msg.setText(params["busy"] or "处理中…")
        self._set_busy(True)
        job.start()

    def _on_progress(self, percent: int, message: str) -> None:
        if self._show_percent:
            self.progress.setRange(0, 100)
            self.progress.setValue(max(0, min(100, percent)))
        else:
            if self.progress.maximum() == 0:
                self.progress.setRange(0, 0)
        if message:
            self.status_msg.setText(message)

    def _on_job_success(self, result) -> None:
        handler = self._job_done
        self._reset_job_ui()
        if handler is not None:
            try:
                handler(result)
            except GGEError as exc:
                self._error(exc.message, exc.hint or "")
            except Exception as exc:  # pragma: no cover
                self._error(f"处理结果时出错：{exc}")

    def _on_job_failed(self, message: str, hint: str) -> None:
        self._reset_job_ui()
        self._error(message, hint)

    def _on_job_cancelled(self) -> None:
        self._reset_job_ui()
        self.status_msg.setText("已取消")

    def _on_job_finished(self) -> None:
        self._job = None
        if self._queue:
            nxt = self._queue.pop(0)
            QTimer.singleShot(0, lambda p=nxt: self._launch_job(p))
        else:
            self._restoring = False
            self._update_actions()

    def _reset_job_ui(self) -> None:
        self.progress.hide()
        self.progress.setRange(0, 100)
        self.btn_cancel.hide()
        self._set_busy(False)
        self._job_done = None

    def _cancel_job(self) -> None:
        if self._queue:
            n = len(self._queue)
            self._queue.clear()
            self.status_msg.setText(f"已取消 {n} 个排队任务")
        if self._job is not None:
            self._job.cancel()
            self.status_msg.setText("正在取消…")

    def _set_busy(self, busy: bool) -> None:
        self.btn_extract.setEnabled(not busy)
        self.search_panel.btn.setEnabled(not busy)
        self.act_reload.setEnabled(not busy and self.provider is None)
        if not busy:
            self._update_actions()

    def _update_actions(self) -> None:
        has_ex = self.extraction is not None
        ready = self.provider is not None and self.model is not None
        self.act_export_gene.setEnabled(has_ex)
        self.act_export_doc.setEnabled(has_ex)
        self.act_export_report.setEnabled(bool(self.validation))
        self.act_batch.setEnabled(ready)
        self.act_export_dataset.setEnabled(self.provider is not None
                                           and self.model is not None
                                           and self._job is None)
        self.btn_extract.setEnabled(ready and bool(self.current_gene_id)
                                    and self._job is None and not self._queue)
        self.act_copy_info.setEnabled(has_ex)
        self.act_reload.setEnabled(self._job is None)

    def _error(self, message: str, hint: str = "") -> None:
        box = QMessageBox(self)
        box.setIcon(QMessageBox.Critical)
        box.setWindowTitle("出错了")
        box.setText(message)
        if hint:
            box.setInformativeText(hint if len(hint) < 600 else hint[:600] + "…")
            box.setDetailedText(hint)
        box.exec()

    def _warn(self, message: str, detail: str = "") -> None:
        box = QMessageBox(self)
        box.setIcon(QMessageBox.Warning)
        box.setWindowTitle("提示")
        box.setText(message)
        if detail:
            box.setDetailedText(detail)
        box.exec()

    def _info(self, message: str, detail: str = "") -> None:
        box = QMessageBox(self)
        box.setIcon(QMessageBox.Information)
        box.setWindowTitle("提示")
        box.setText(message)
        if detail:
            box.setDetailedText(detail)
        box.exec()

    def _show_help(self) -> None:
        AboutDialog(HELP_HTML.format(ver=__version__), self).exec()

    def _show_about(self) -> None:
        AboutDialog(ABOUT_HTML.format(ver=__version__, tests=243), self).exec()

    # ================================================================== #
    # 拖放
    # ================================================================== #
    def dragEnterEvent(self, event) -> None:  # noqa: N802
        if event.mimeData().hasUrls():
            event.acceptProposedAction()

    def dropEvent(self, event) -> None:  # noqa: N802
        for url in event.mimeData().urls():
            path = url.toLocalFile()
            if not os.path.isfile(path):
                continue
            if is_archive(path) or path.lower().endswith((".7z", ".rar")):
                self.import_archive_dialog(path)
                return
            low = path.lower()
            if low.endswith((".fa", ".fasta", ".fna", ".fas", ".fa.gz", ".fna.gz")):
                role = "genome" if self.provider is None else (
                    "mrna" if self.mrna_refs is None else (
                        "cds" if self.cds_refs is None else "protein"))
                self.load_role(role, path)
            elif low.endswith((".gff3", ".gff", ".gtf", ".gff3.gz", ".gtf.gz")):
                self.load_role("annotation", path)
            else:
                self.status_msg.setText(f"无法识别的文件类型：{os.path.basename(path)}")


def _demo_dir() -> Optional[str]:
    """定位示例数据目录（兼容 PyInstaller 打包后的路径）。"""
    candidates = []
    if getattr(sys, "frozen", False):
        base = os.path.dirname(sys.executable)
        candidates.append(os.path.join(base, "examples"))
        candidates.append(os.path.join(getattr(sys, "_MEIPASS", base), "examples"))
    here = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    candidates.append(os.path.join(here, "examples"))
    candidates.append(os.path.join(os.getcwd(), "examples"))
    for c in candidates:
        if os.path.isfile(os.path.join(c, "demo_genome.fa")):
            return c
    return None
