"""本地数据集库（SQLite）。

把"一套数据"（基因组 + 注释 + 可选的 mRNA/CDS/蛋白参考）登记成一条记录，
下次启动时不用再翻文件夹。使用标准库 ``sqlite3``，不引入额外依赖，也不需要
管理员权限：

* Windows: ``%APPDATA%\\GenomeGeneExplorer\\library.sqlite``
* 其它平台: ``~/.config/GenomeGeneExplorer/library.sqlite``

数据库里只存**元数据与文件路径**，序列本身仍在解包目录里。这样：

* 打开软件不必读几十 GB 的基因组；
* 数据可以放在任意磁盘（甚至移动硬盘），只需路径仍然有效；
* 删掉数据库不会丢数据，删掉数据目录则会在库里显示"文件缺失"。

表结构
------
``datasets``  一条记录 = 一套数据；``files``  每个角色一行，记录包内原始文件名和大小。
``meta``      记录 schema 版本，便于以后升级。

删除数据集时只会递归删除**本程序自己的存储目录**（``library/datasets/``）下的文件，
用户自己指定到别处的路径一律不删，只删记录。
"""

from __future__ import annotations

import os
import re
import shutil
import sqlite3
from contextlib import contextmanager
from dataclasses import dataclass, field
from datetime import datetime
from typing import Dict, Iterable, Iterator, List, Optional

from .config import config_dir

SCHEMA_VERSION = 1

ROLE_ORDER = ("genome", "annotation", "mrna", "cds", "protein")
ROLE_LABELS = {
    "genome": "基因组",
    "annotation": "注释",
    "mrna": "mRNA",
    "cds": "CDS",
    "protein": "蛋白",
}
REQUIRED_ROLES = ("genome", "annotation")


def library_dir() -> str:
    """库根目录（数据库 + 解包后的数据）。"""
    override = os.environ.get("GGE_LIBRARY_DIR", "")
    if override:
        path = override
    else:
        path = os.path.join(config_dir(), "library")
    os.makedirs(path, exist_ok=True)
    return path


def default_db_path() -> str:
    override = os.environ.get("GGE_LIBRARY_DB", "")
    if override:
        os.makedirs(os.path.dirname(override) or ".", exist_ok=True)
        return override
    return os.path.join(library_dir(), "library.sqlite")


def store_root() -> str:
    """导入的数据解包到这里，删除数据集时也只会动这里。"""
    path = os.path.join(library_dir(), "datasets")
    os.makedirs(path, exist_ok=True)
    return path


def _now() -> str:
    # 精确到微秒：同一秒内连续导入/使用时排序才有意义（显示时只取到分钟）
    return datetime.now().isoformat(timespec="microseconds")


def slugify(text: str, fallback: str = "dataset") -> str:
    """把任意名字变成可用作目录名的 ASCII 片段。

    中文名不丢：Windows 下中文目录是可用的，但为了跨工具兼容，这里保留
    ASCII 字母数字，其余替换成 ``_``；若全部被替换掉则退回 ``fallback``。
    """
    keep = re.sub(r"[^0-9A-Za-z._-]+", "_", text or "").strip("._-")
    return (keep or fallback)[:40]


@dataclass
class DatasetFile:
    role: str
    path: str
    member: str = ""
    size: int = 0


@dataclass
class Dataset:
    id: int = 0
    name: str = ""
    source: str = ""
    fmt: str = ""
    n_genes: int = 0
    n_transcripts: int = 0
    n_seqids: int = 0
    total_bases: int = 0
    note: str = ""
    created_at: str = ""
    updated_at: str = ""
    last_used: str = ""
    files: Dict[str, DatasetFile] = field(default_factory=dict)

    # -- 便捷访问 ---------------------------------------------------------- #
    def path(self, role: str) -> str:
        f = self.files.get(role)
        return f.path if f else ""

    @property
    def genome(self) -> str:
        return self.path("genome")

    @property
    def annotation(self) -> str:
        return self.path("annotation")

    def paths(self) -> Dict[str, str]:
        return {r: f.path for r, f in self.files.items() if f.path}

    def missing(self) -> List[str]:
        return [r for r, f in self.files.items() if f.path and not os.path.exists(f.path)]

    def missing_required(self) -> List[str]:
        return [r for r in REQUIRED_ROLES
                if not self.path(r) or not os.path.exists(self.path(r))]

    def total_size(self) -> int:
        return sum(f.size for f in self.files.values() if f.size)

    def summary(self) -> str:
        bits = [f"基因 {self.n_genes:,}" if self.n_genes else "",
                f"转录本 {self.n_transcripts:,}" if self.n_transcripts else "",
                f"参考 {len([r for r in ('mrna', 'cds', 'protein') if self.path(r)])} 个"]
        return " · ".join([b for b in bits if b])

    def row(self) -> List[str]:
        return [
            self.name,
            self.summary() or "—",
            self.fmt.upper() or "—",
            self.source_name() or "（手动指定）",
            (self.last_used or self.updated_at or "")[:16].replace("T", " "),
        ]

    def source_name(self) -> str:
        return os.path.basename(self.source) if self.source else ""


class Library:
    """数据集库。所有方法都可以在 GUI 线程直接调用（内部按调用开关连接）。"""

    def __init__(self, db_path: Optional[str] = None) -> None:
        self.path = os.path.abspath(db_path or default_db_path())
        os.makedirs(os.path.dirname(self.path) or ".", exist_ok=True)
        self._init_schema()

    # -- 连接 -------------------------------------------------------------- #
    @contextmanager
    def _conn(self) -> Iterator[sqlite3.Connection]:
        conn = sqlite3.connect(self.path, timeout=10)
        conn.row_factory = sqlite3.Row
        try:
            conn.execute("PRAGMA foreign_keys = ON")
            yield conn
            conn.commit()
        finally:
            conn.close()

    def _init_schema(self) -> None:
        with self._conn() as conn:
            conn.executescript(
                """
                CREATE TABLE IF NOT EXISTS meta (
                    key   TEXT PRIMARY KEY,
                    value TEXT NOT NULL
                );
                CREATE TABLE IF NOT EXISTS datasets (
                    id            INTEGER PRIMARY KEY AUTOINCREMENT,
                    name          TEXT NOT NULL,
                    source        TEXT NOT NULL DEFAULT '',
                    fmt           TEXT NOT NULL DEFAULT '',
                    n_genes       INTEGER NOT NULL DEFAULT 0,
                    n_transcripts INTEGER NOT NULL DEFAULT 0,
                    n_seqids      INTEGER NOT NULL DEFAULT 0,
                    total_bases   INTEGER NOT NULL DEFAULT 0,
                    note          TEXT NOT NULL DEFAULT '',
                    created_at    TEXT NOT NULL,
                    updated_at    TEXT NOT NULL,
                    last_used     TEXT NOT NULL DEFAULT ''
                );
                CREATE TABLE IF NOT EXISTS files (
                    id         INTEGER PRIMARY KEY AUTOINCREMENT,
                    dataset_id INTEGER NOT NULL REFERENCES datasets(id) ON DELETE CASCADE,
                    role       TEXT NOT NULL,
                    member     TEXT NOT NULL DEFAULT '',
                    path       TEXT NOT NULL,
                    size       INTEGER NOT NULL DEFAULT 0,
                    UNIQUE(dataset_id, role)
                );
                CREATE INDEX IF NOT EXISTS idx_files_dataset ON files(dataset_id);
                CREATE INDEX IF NOT EXISTS idx_datasets_source ON datasets(source);
                """
            )
            conn.execute("INSERT OR REPLACE INTO meta(key, value) VALUES('schema_version', ?)",
                         (str(SCHEMA_VERSION),))

    # -- 读 ----------------------------------------------------------------- #
    def _row_to_dataset(self, row: sqlite3.Row, files: Iterable[sqlite3.Row]) -> Dataset:
        ds = Dataset(
            id=row["id"], name=row["name"], source=row["source"], fmt=row["fmt"],
            n_genes=row["n_genes"], n_transcripts=row["n_transcripts"],
            n_seqids=row["n_seqids"], total_bases=row["total_bases"], note=row["note"],
            created_at=row["created_at"], updated_at=row["updated_at"],
            last_used=row["last_used"],
        )
        for f in files:
            ds.files[f["role"]] = DatasetFile(role=f["role"], path=f["path"],
                                             member=f["member"], size=f["size"])
        return ds

    def list(self, search: str = "", limit: int = 500) -> List[Dataset]:
        """按"最近使用"排序；从未使用过的排在后面（按创建时间倒序）。"""
        with self._conn() as conn:
            if search:
                like = f"%{search}%"
                rows = conn.execute(
                    "SELECT * FROM datasets WHERE name LIKE ? OR source LIKE ? "
                    "ORDER BY COALESCE(NULLIF(last_used,''), created_at) DESC, id DESC "
                    "LIMIT ?",
                    (like, like, limit)).fetchall()
            else:
                rows = conn.execute(
                    "SELECT * FROM datasets "
                    "ORDER BY COALESCE(NULLIF(last_used,''), created_at) DESC, id DESC "
                    "LIMIT ?",
                    (limit,)).fetchall()
            out = []
            for row in rows:
                files = conn.execute("SELECT * FROM files WHERE dataset_id=?",
                                     (row["id"],)).fetchall()
                out.append(self._row_to_dataset(row, files))
            return out

    def get(self, dataset_id: int) -> Optional[Dataset]:
        with self._conn() as conn:
            row = conn.execute("SELECT * FROM datasets WHERE id=?", (dataset_id,)).fetchone()
            if row is None:
                return None
            files = conn.execute("SELECT * FROM files WHERE dataset_id=?",
                                 (dataset_id,)).fetchall()
            return self._row_to_dataset(row, files)

    def count(self) -> int:
        with self._conn() as conn:
            return int(conn.execute("SELECT COUNT(*) FROM datasets").fetchone()[0])

    def get_by_source(self, source: str) -> Optional[Dataset]:
        """同一个压缩包重复导入时，用它判断是"更新"而不是"再存一份"。"""
        if not source:
            return None
        with self._conn() as conn:
            row = conn.execute(
                "SELECT * FROM datasets WHERE source=? ORDER BY id DESC LIMIT 1",
                (os.path.abspath(source),)).fetchone()
            if row is None:
                return None
            files = conn.execute("SELECT * FROM files WHERE dataset_id=?",
                                 (row["id"],)).fetchall()
            return self._row_to_dataset(row, files)

    # -- 写 ----------------------------------------------------------------- #
    def add(self, name: str, paths: Dict[str, str], source: str = "", fmt: str = "",
            stats: Optional[Dict[str, int]] = None, note: str = "",
            members: Optional[Dict[str, str]] = None) -> Dataset:
        stats = stats or {}
        members = members or {}
        now = _now()
        with self._conn() as conn:
            cur = conn.execute(
                "INSERT INTO datasets(name, source, fmt, n_genes, n_transcripts, n_seqids,"
                " total_bases, note, created_at, updated_at, last_used)"
                " VALUES(?,?,?,?,?,?,?,?,?,?,?)",
                (name, os.path.abspath(source) if source else "", fmt,
                 int(stats.get("n_genes", 0)), int(stats.get("n_transcripts", 0)),
                 int(stats.get("n_seqids", 0)), int(stats.get("total_bases", 0)),
                 note, now, now, now))
            dataset_id = int(cur.lastrowid)
            self._insert_files(conn, dataset_id, paths, members)
        ds = self.get(dataset_id)
        assert ds is not None
        return ds

    def update(self, dataset_id: int, paths: Optional[Dict[str, str]] = None,
               name: Optional[str] = None, source: Optional[str] = None,
               fmt: Optional[str] = None, stats: Optional[Dict[str, int]] = None,
               note: Optional[str] = None,
               members: Optional[Dict[str, str]] = None) -> Optional[Dataset]:
        fields = ["updated_at=?"]
        values: List[object] = [_now()]
        if name is not None:
            fields.append("name=?")
            values.append(name)
        if source is not None:
            fields.append("source=?")
            values.append(os.path.abspath(source) if source else "")
        if fmt is not None:
            fields.append("fmt=?")
            values.append(fmt)
        if note is not None:
            fields.append("note=?")
            values.append(note)
        for key in ("n_genes", "n_transcripts", "n_seqids", "total_bases"):
            if stats and key in stats:
                fields.append(f"{key}=?")
                values.append(int(stats[key]))
        values.append(dataset_id)
        with self._conn() as conn:
            conn.execute(f"UPDATE datasets SET {', '.join(fields)} WHERE id=?", values)
            if paths is not None:
                conn.execute("DELETE FROM files WHERE dataset_id=?", (dataset_id,))
                self._insert_files(conn, dataset_id, paths, members or {})
        return self.get(dataset_id)

    def _insert_files(self, conn: sqlite3.Connection, dataset_id: int,
                      paths: Dict[str, str], members: Dict[str, str]) -> None:
        for role in ROLE_ORDER:
            path = paths.get(role, "")
            if not path:
                continue
            size = 0
            try:
                size = os.path.getsize(path)
            except OSError:
                pass
            conn.execute(
                "INSERT OR REPLACE INTO files(dataset_id, role, member, path, size)"
                " VALUES(?,?,?,?,?)",
                (dataset_id, role, members.get(role, os.path.basename(path)),
                 os.path.abspath(path), size))

    def rename(self, dataset_id: int, name: str) -> None:
        name = (name or "").strip()
        if not name:
            return
        with self._conn() as conn:
            conn.execute("UPDATE datasets SET name=?, updated_at=? WHERE id=?",
                         (name, _now(), dataset_id))

    def touch(self, dataset_id: int) -> None:
        with self._conn() as conn:
            conn.execute("UPDATE datasets SET last_used=? WHERE id=?", (_now(), dataset_id))

    def remove(self, dataset_id: int, delete_files: bool = True) -> List[str]:
        """删记录；``delete_files`` 为真时删除**库目录内**的文件。返回被删的路径。"""
        ds = self.get(dataset_id)
        if ds is None:
            return []
        removed: List[str] = []
        if delete_files:
            for f in ds.files.values():
                target = self._owned_dir(f.path)
                if target and os.path.isdir(target) and target not in removed:
                    shutil.rmtree(target, ignore_errors=True)
                    removed.append(target)
        with self._conn() as conn:
            conn.execute("DELETE FROM files WHERE dataset_id=?", (dataset_id,))
            conn.execute("DELETE FROM datasets WHERE id=?", (dataset_id,))
        return removed

    def owns_files(self, dataset: Dataset) -> bool:
        """这套数据的文件是否就在本库目录里（决定移除时要不要删文件）。"""
        return any(self._owned_dir(f.path) for f in dataset.files.values())

    def _owned_dir(self, file_path: str) -> str:
        """只有位于库存储目录内的文件才允许删除，避免误删用户的原始数据。"""
        root = os.path.abspath(store_root())
        target = os.path.abspath(os.path.dirname(file_path or ""))
        try:
            if os.path.commonpath([root, target]).lower() == root.lower() and target != root:
                return target
        except ValueError:      # 不同盘符
            pass
        return ""

    # -- 目录 -------------------------------------------------------------- #
    def new_dataset_dir(self, name: str) -> str:
        """为一次导入分配一个新目录，名字带时间戳，避免撞名。"""
        stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        base = f"{slugify(name)}_{stamp}"
        path = os.path.join(store_root(), base)
        i = 2
        while os.path.exists(path):
            path = os.path.join(store_root(), f"{base}_{i}")
            i += 1
        os.makedirs(path, exist_ok=False)
        return path
