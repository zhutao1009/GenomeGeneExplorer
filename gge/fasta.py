"""基因组 FASTA 随机访问。

策略
----
1. 优先使用 ``pyfaidx.Fasta``：它用 ``.fai`` 索引做 ``seek``，几十 GB 的基因组
   只占用极小内存，这是我们处理大基因组的首选。
2. 若 pyfaidx 不可用，或文件过小（< 64 MB）且索引建立失败，退回到
   ``PlainFastaProvider``：一次性读入内存并建立 (offset, length) 索引。
   此时会在状态里标注 "内存模式"，让用户知道内存占用。
3. 索引比 FASTA 旧时（重新下载序列的常见情形）自动重建一次并提示。

线程
----
pyfaidx 的文件句柄不是线程安全的，这里用 ``threading.RLock`` 串行化所有读取。
提取任务本身在后台线程执行，因此这个锁是必需的。
"""

from __future__ import annotations

import os
import re
import threading
from typing import Dict, Iterable, List, Optional, Tuple

from .coords import revcomp
from .errors import FileFormatError

FASTA_SUFFIXES = (".fa", ".fasta", ".fna", ".fas", ".fa.gz", ".fasta.gz", ".fna.gz")


def looks_like_fasta(path: str) -> bool:
    return path.lower().endswith(FASTA_SUFFIXES)


def _sanitize_seqid(name: str) -> str:
    """取 FASTA 标题的第一个 token 作为序列 ID。"""
    return name.strip().split()[0] if name.strip() else name.strip()


class BaseFastaProvider:
    """序列访问接口。所有实现都必须满足：fetch 返回**正链**基因组序列。"""

    path: str = ""
    mode: str = "unknown"

    def seqids(self) -> List[str]:
        raise NotImplementedError

    def length(self, seqid: str) -> Optional[int]:
        raise NotImplementedError

    def lengths(self) -> Dict[str, int]:
        raise NotImplementedError

    def fetch(self, seqid: str, start: int, end: int) -> str:
        """取 1-based 闭区间 ``[start, end]`` 的正链序列，超界自动裁剪。"""
        raise NotImplementedError

    def fetch_oriented(self, seqid: str, start: int, end: int, strand: str) -> str:
        seq = self.fetch(seqid, start, end)
        return revcomp(seq) if strand == "-" else seq

    def close(self) -> None:
        pass

    # 便于测试的辅助
    def has(self, seqid: str) -> bool:
        return self.length(seqid) is not None

    def suggest(self, seqid: str, limit: int = 5) -> List[str]:
        """当 seqid 找不到时，给出相似 ID 建议（chr1 vs 1 vs NC_xxx）。"""
        target = _normalize_seqid(seqid)
        scored: List[Tuple[int, str]] = []
        for cand in self.seqids():
            score = _seqid_similarity(target, _normalize_seqid(cand))
            if score > 0:
                scored.append((score, cand))
        scored.sort(key=lambda x: (-x[0], x[1]))
        return [c for _, c in scored[:limit]]


def _normalize_seqid(s: str) -> str:
    s = s.strip().lower()
    s = re.sub(r"^chr(?=[0-9mxy])", "", s)
    s = s.split(".")[0]
    s = s.replace("_", "").replace("-", "")
    return s


def _seqid_similarity(a: str, b: str) -> int:
    if not a or not b:
        return 0
    if a == b:
        return 100
    if a in b or b in a:
        return 60 - abs(len(a) - len(b))
    # 公共前缀
    n = 0
    for x, y in zip(a, b):
        if x != y:
            break
        n += 1
    return n * 3 if n >= 2 else 0


def _faidx_entry_lengths(fa) -> Dict[str, int]:
    """取出每条序列的长度，兼容 pyfaidx 0.6 ~ 0.9+。

    pyfaidx 的索引对象搬过两次家：
      * 0.6/0.7：``fa.index``（dict，值带 ``rlen`` 字段）
      * 0.9+ ：``fa.index`` 被删除，改名 ``fa.records``；真正的索引表在 ``fa.faidx.index``

    注意 ``IndexEntry.length`` 是**文件里的字节跨度**，不是碱基数，所以只认 ``rlen``。
    拿不到属性时退回 ``len(fa[name])``——pyfaidx 的 ``FastaRecord.__len__`` 直接读
    ``rlen``，不会把整条序列读进内存（31 条染色体实测 < 1 ms）。
    """
    lengths: Dict[str, int] = {}

    # 路径一：直接从索引条目取字段（最快，且不构造任何记录对象）
    for holder in (getattr(fa, "faidx", None), fa):
        if holder is None:
            continue
        for attr in ("index", "records"):
            table = getattr(holder, attr, None)
            if not isinstance(table, dict):
                continue
            for name, entry in table.items():
                rlen = getattr(entry, "rlen", None)
                if isinstance(rlen, int):
                    lengths[name] = rlen
            if lengths:
                return lengths

    # 路径二：任何版本都能用的兜底
    for name in fa.keys():
        try:
            lengths[name] = len(fa[name])
        except Exception:
            continue
    return lengths


class PyfaidxProvider(BaseFastaProvider):
    """基于 pyfaidx 的磁盘随机访问实现。"""

    mode = "pyfaidx（磁盘随机访问）"

    def __init__(self, path: str, progress_cb=None) -> None:
        import pyfaidx  # 延迟导入，缺依赖时给出明确提示

        self.path = path
        self._lock = threading.RLock()
        if progress_cb:
            progress_cb(5, "检查 FASTA 索引 …")
        try:
            self._fa = pyfaidx.Fasta(path, as_raw=True, sequence_always_upper=False,
                                     duplicate_action="first")
        except Exception as exc:  # 索引缺失/过期
            if progress_cb:
                progress_cb(15, "索引不可用，正在重建 .fai …")
            try:
                self._fa = pyfaidx.Fasta(path, as_raw=True, sequence_always_upper=False,
                                         duplicate_action="first", rebuild=True)
            except Exception as exc2:
                raise FileFormatError(
                    f"无法为 FASTA 建立索引：{exc2}",
                    "确认文件未被其他程序占用且磁盘可写；也可先手动运行 "
                    "`samtools faidx <file.fa>`。",
                ) from exc2
        self._lengths = _faidx_entry_lengths(self._fa)
        if not self._lengths:
            raise FileFormatError(
                f"FASTA 索引为空：{os.path.basename(path)}",
                "文件可能只有标题行没有序列；确认后用 `samtools faidx` 手动建索引再试。",
            )
        # pyfaidx 的 keys() 可能带 '>name' 形式，这里统一清洗
        cleaned: Dict[str, int] = {}
        for k, v in self._lengths.items():
            cleaned[_sanitize_seqid(k)] = v
        self._lengths = cleaned
        self._seqids = list(cleaned.keys())

    def seqids(self) -> List[str]:
        return list(self._seqids)

    def length(self, seqid: str) -> Optional[int]:
        if seqid in self._lengths:
            return self._lengths[seqid]
        return None

    def lengths(self) -> Dict[str, int]:
        return dict(self._lengths)

    def fetch(self, seqid: str, start: int, end: int) -> str:
        n = self.length(seqid)
        if n is None:
            raise KeyError(seqid)
        if start < 1:
            start = 1
        if end > n:
            end = n
        if start > end:
            return ""
        with self._lock:
            rec = self._fa[seqid]
            seq = rec[start - 1:end]
            if hasattr(seq, "seq"):
                seq = str(seq.seq)
            else:
                seq = str(seq)
        return seq

    def close(self) -> None:
        try:
            self._fa.close()
        except Exception:
            pass


class PlainFastaProvider(BaseFastaProvider):
    """内存索引实现：适合测试与小型基因组（< 64 MB）。"""

    mode = "内存索引"

    def __init__(self, path: str, max_bytes: int = 64 * 1024 * 1024) -> None:
        self.path = path
        size = os.path.getsize(path)
        if size > max_bytes:
            raise FileFormatError(
                f"文件过大（{size / 1e6:.1f} MB），无法使用内存模式。",
                "请安装 pyfaidx：`pip install pyfaidx`，以便对大基因组做磁盘随机访问。",
            )
        self._seqs: Dict[str, str] = {}
        self._lengths: Dict[str, int] = {}
        name: Optional[str] = None
        chunks: List[str] = []
        opener = open
        if path.lower().endswith(".gz"):
            import gzip

            opener = gzip.open
        try:
            with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    if line.startswith(">"):
                        if name is not None:
                            self._seqs[name] = "".join(chunks)
                        header = _sanitize_seqid(line[1:])
                        name = header
                        if name in self._seqs:
                            # 重复 ID：保留第一条，与 pyfaidx duplicate_action='first' 一致
                            name = name + "#dup"
                        chunks = []
                    else:
                        if name is None:
                            raise FileFormatError(
                                "FASTA 文件在首个 '>' 标题行之前就出现了序列内容。",
                                "确认选择的是 FASTA 文件而不是其他格式。",
                            )
                        chunks.append(line.strip())
                if name is not None:
                    self._seqs[name] = "".join(chunks)
        except UnicodeDecodeError as exc:
            raise FileFormatError(f"FASTA 文件不是文本格式：{exc}") from exc
        if not self._seqs:
            raise FileFormatError(
                "FASTA 文件中没有找到任何序列。",
                "检查文件是否为空，或标题行是否以 '>' 开头。",
            )
        self._lengths = {k: len(v) for k, v in self._seqs.items()}

    def seqids(self) -> List[str]:
        return list(self._seqs.keys())

    def length(self, seqid: str) -> Optional[int]:
        return self._lengths.get(seqid)

    def lengths(self) -> Dict[str, int]:
        return dict(self._lengths)

    def fetch(self, seqid: str, start: int, end: int) -> str:
        seq = self._seqs.get(seqid)
        if seq is None:
            raise KeyError(seqid)
        n = len(seq)
        start = max(1, start)
        end = min(n, end)
        if start > end:
            return ""
        return seq[start - 1:end]


def open_fasta(path: str, progress_cb=None, prefer_memory: bool = False) -> BaseFastaProvider:
    """打开基因组 FASTA 并建立索引，自动在两种实现之间选择。

    磁盘模式失败的**真实原因必须保留**：以前这里写的是 ``except Exception: pass``，
    于是 pyfaidx 版本的任何不兼容（例如 0.9 删掉了 ``Fasta.index``）都会被吞掉，
    最后只报出内存模式的"文件过大，请安装 pyfaidx"——明明 pyfaidx 装着且能用，
    用户却被指向一个完全错误的方向。现在两种模式都失败时，两个原因一起报出来。
    """
    if not os.path.exists(path):
        raise FileFormatError(
            f"文件不存在：{path}",
            "确认路径是否正确，或文件是否被移动/删除。",
        )
    if os.path.getsize(path) == 0:
        raise FileFormatError(f"文件大小为 0：{path}", "重新下载该 FASTA 文件。")

    disk_error = None          # (描述, 原始异常)
    if not prefer_memory:
        try:
            return PyfaidxProvider(path, progress_cb=progress_cb)
        except ImportError as exc:
            disk_error = ("pyfaidx 未安装", exc)
        except FileFormatError:
            raise                   # 已经是带"怎么办"的明确错误，别降级掩盖
        except Exception as exc:
            disk_error = (f"磁盘索引失败（{type(exc).__name__}: {exc}）", exc)

    try:
        return PlainFastaProvider(path)
    except FileFormatError as exc:
        if disk_error is None:
            raise
        raise FileFormatError(
            f"无法打开 FASTA：{os.path.basename(path)}\n"
            f"　· 磁盘模式：{disk_error[0]}\n"
            f"　· 内存模式：{exc}",
            "多数情况是上面「磁盘模式」的原因导致的（内存模式只适合 < 64 MB 的文件）。"
            "可以先把这条错误发出来；也可以手动运行 `samtools faidx <文件.fa>` 生成 .fai 后再打开。",
        ) from disk_error[1]


def read_fasta_records(path: str, progress_cb=None,
                       cancel_check=None) -> "list[tuple[str, str]]":
    """读取普通 FASTA（参考 mRNA / CDS / 蛋白），返回 [(id, seq)]，保持文件顺序。

    这里不用 pyfaidx：参考序列文件通常很小，而且我们需要**全部**记录的内容
    用于 ID 匹配与"未使用的参考序列"统计。
    """
    if not os.path.exists(path):
        raise FileFormatError(f"文件不存在：{path}", "确认路径是否正确。")
    records: List[Tuple[str, str]] = []
    name: Optional[str] = None
    chunks: List[str] = []
    opener = open
    if path.lower().endswith(".gz"):
        import gzip

        opener = gzip.open
    line_no = 0
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line_no += 1
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name = line[1:].strip()
                chunks = []
            elif line.startswith(";"):
                continue  # 旧式注释行
            else:
                if name is None:
                    raise FileFormatError(
                        f"{os.path.basename(path)} 第 {line_no} 行：序列出现在标题行之前。",
                        "确认文件是 FASTA 格式（以 '>' 开头）。",
                    )
                chunks.append(line.strip())
            if line_no % 20000 == 0:
                if cancel_check and cancel_check():     # 250 MB 的参考也要能中途取消
                    from .errors import CancelledError

                    raise CancelledError()
                if progress_cb:
                    progress_cb(50, f"已读取 {line_no:,} 行 …")
        if name is not None:
            records.append((name, "".join(chunks)))
    if not records:
        raise FileFormatError(
            f"{os.path.basename(path)} 中没有解析到任何序列记录。",
            "检查文件是否为空或格式是否被破坏。",
        )
    # 去除空记录（只有标题没有序列）
    cleaned = [(n, s) for n, s in records if s.strip()]
    if not cleaned:
        raise FileFormatError(
            f"{os.path.basename(path)} 中所有记录都没有序列内容。",
            "重新导出参考序列文件。",
        )
    return cleaned
