"""Genome Gene Explorer (GGE)

基因组序列提取、可视化与注释校验桌面软件的核心库。

模块划分
--------
gge.models       数据模型（Gene / Transcript / Segment / Panel / Block / Region）
gge.errors       统一的异常类型
gge.fasta        基因组 FASTA 随机访问（pyfaidx 优先，小文件回退到内存索引）
gge.annotation   GFF3 / GTF 解析与层级重建
gge.search       基因 / 转录本的模糊检索
gge.coords       坐标与链方向工具（1-based 闭区间 <-> 0-based 半开区间）
gge.extract      序列提取核心算法
gge.validate     参考 mRNA / CDS 匹配与比对校验
gge.export       FASTA / TSV / CSV / 报告导出
gge.config       参数持久化
gge.colors       配色方案
gge.ui.*         PySide6 图形界面
"""

from .version import __version__, APP_NAME

__all__ = ["__version__", "APP_NAME"]
