"""参数与路径持久化。

用 JSON 存在用户目录下（不需要管理员权限，也不写注册表）：

* Windows: ``%APPDATA%\\GenomeGeneExplorer\\settings.json``
* 其它平台: ``~/.config/GenomeGeneExplorer/settings.json``

只保存"下次启动该恢复什么"，不保存任何序列数据；写入采用临时文件 + 替换，
避免断电时留下半个文件导致下次启动崩溃。
"""

from __future__ import annotations

import json
import os
import sys
import tempfile
from typing import Any, Dict, List, Optional

from .colors import DEFAULT_PALETTE, ensure_palette
from .version import APP_ID

SETTINGS_VERSION = 5


def config_dir() -> str:
    if sys.platform.startswith("win"):
        base = os.environ.get("APPDATA") or os.path.expanduser("~")
    elif sys.platform == "darwin":
        base = os.path.expanduser("~/Library/Application Support")
    else:
        base = os.environ.get("XDG_CONFIG_HOME") or os.path.expanduser("~/.config")
    path = os.path.join(base, APP_ID)
    os.makedirs(path, exist_ok=True)
    return path


def settings_path() -> str:
    return os.path.join(config_dir(), "settings.json")


def default_settings() -> Dict[str, Any]:
    return {
        "version": SETTINGS_VERSION,
        "genome_path": "",
        "annotation_path": "",
        "mrna_path": "",
        "cds_path": "",
        "protein_path": "",
        "last_dir": "",
        "upstream": 2000,
        "downstream": 2000,
        "bases_per_line": 60,
        "font_size": 11,
        "panel_columns": 2,
        "palette": {k: dict(v) for k, v in DEFAULT_PALETTE.items()},
        "window_geometry": "",
        "splitter_states": {},
        "auto_fit": True,
        "dataset_id": 0,
        "last_archive_dir": "",
        "recent_genes": [],
        "auto_validate": True,
        "wrap_sequence": False,
        "export_dir": "",
        # 启动时把上次的文件真正加载回来（只填路径会让用户以为已就绪）
        "auto_restore": True,
        # 导出的彩色文档默认不带左侧坐标列：贴进论文/报告更干净
        "doc_include_coord": False,
        # 上次查看的基因：重启自动恢复时一并打开
        "last_gene": "",
    }


class Settings:
    """轻量设置容器：``settings["upstream"]`` 读写，自动补默认值。"""

    def __init__(self, path: Optional[str] = None) -> None:
        self.path = path or settings_path()
        self._data: Dict[str, Any] = default_settings()
        self.load()

    # -- 读写 -------------------------------------------------------------- #
    def load(self) -> None:
        if not os.path.exists(self.path):
            return
        try:
            with open(self.path, "r", encoding="utf-8") as fh:
                data = json.load(fh)
        except (OSError, json.JSONDecodeError):
            # 配置文件损坏不应阻止软件启动
            return
        if not isinstance(data, dict):
            return
        defaults = default_settings()
        for key, val in data.items():
            if key in defaults and isinstance(defaults[key], type(val)):
                self._data[key] = val
            elif key in defaults and defaults[key] is None:
                self._data[key] = val
            elif key in defaults and isinstance(defaults[key], (int, float)) and isinstance(val, (int, float)):
                self._data[key] = val
        # 兼容：palette 只做了部分覆盖时补齐
        self._data["palette"] = ensure_palette(self._data.get("palette"))
        if not isinstance(self._data.get("recent_genes"), list):
            self._data["recent_genes"] = []
        if not isinstance(self._data.get("splitter_states"), dict):
            self._data["splitter_states"] = {}
        # 旧的单份分隔条状态（vsplit）已被 splitter_states 取代，直接丢弃
        self._data.pop("splitter_state", None)

    def save(self) -> bool:
        data = dict(self._data)
        data["version"] = SETTINGS_VERSION
        payload = json.dumps(data, ensure_ascii=False, indent=2)
        try:
            folder = os.path.dirname(self.path) or "."
            os.makedirs(folder, exist_ok=True)
            fd, tmp = tempfile.mkstemp(prefix=".settings_", dir=folder, text=True)
            with os.fdopen(fd, "w", encoding="utf-8") as fh:
                fh.write(payload)
            os.replace(tmp, self.path)
            return True
        except OSError:
            return False

    # -- 字典接口 ---------------------------------------------------------- #
    def __getitem__(self, key: str) -> Any:
        return self._data.get(key, default_settings().get(key))

    def __setitem__(self, key: str, value: Any) -> None:
        self._data[key] = value

    def get(self, key: str, default: Any = None) -> Any:
        return self._data.get(key, default)

    def update(self, **kwargs: Any) -> None:
        self._data.update(kwargs)

    def as_dict(self) -> Dict[str, Any]:
        return dict(self._data)

    # -- 便捷 -------------------------------------------------------------- #
    def add_recent_gene(self, gene_id: str, limit: int = 15) -> None:
        items: List[str] = [g for g in self._data.get("recent_genes", []) if g != gene_id]
        items.insert(0, gene_id)
        self._data["recent_genes"] = items[:limit]

    def palette(self) -> Dict[str, Dict[str, str]]:
        return ensure_palette(self._data.get("palette"))

    def set_palette(self, palette: Dict[str, Dict[str, str]]) -> None:
        self._data["palette"] = ensure_palette(palette)

    def existing_path(self, key: str) -> str:
        """返回仍然存在的路径，否则返回空串（避免启动时报"文件不存在"）。"""
        p = self._data.get(key) or ""
        return p if p and os.path.exists(p) else ""
