"""统一异常类型。

每个异常都带 `hint` 字段，GUI 在弹窗中把 message 与 hint 一起显示，
避免出现"静默失败"或只有一行无信息量的报错。
"""

from __future__ import annotations

from typing import Optional


class GGEError(Exception):
    """所有业务异常的基类。"""

    def __init__(self, message: str, hint: Optional[str] = None):
        super().__init__(message)
        self.message = message
        self.hint = hint

    def full_text(self) -> str:
        if self.hint:
            return f"{self.message}\n\n建议：{self.hint}"
        return self.message


class FileFormatError(GGEError):
    """文件无法读取或格式不符合预期。"""


class AnnotationError(GGEError):
    """注释文件内容无法解析。"""


class SeqIDMismatchError(GGEError):
    """FASTA 与注释文件的序列 ID 对不上。"""


class CoordinateError(GGEError):
    """坐标非法（越界、起止颠倒、长度为 0 等）。"""


class NotReadyError(GGEError):
    """尚未加载必需文件就发起了操作。"""


class CancelledError(Exception):
    """用户主动取消。

    不继承 ``GGEError``：取消是正常流程，不该弹出"出错了"对话框。
    ``gge.tasks`` 里为兼容旧代码保留了同名别名。
    """
