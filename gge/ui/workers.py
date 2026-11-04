"""后台任务线程。

所有耗时操作（建索引、解析注释、提取、比对、批量导出）都走这里，主线程只
接收信号，因此大基因组加载时界面不会"未响应"。

取消：``Job.cancel()`` 只是置标志，任务函数在安全检查点抛出 ``CancelledError``。
这样不会强杀线程（强杀会让 pyfaidx 的文件句柄处于未知状态）。
"""

from __future__ import annotations

import inspect
import traceback
from typing import Any, Callable

from PySide6.QtCore import QThread, Signal

from ..errors import GGEError
from ..tasks import CancelledError


class Job(QThread):
    """把任意函数放到工作线程执行。"""

    progressed = Signal(int, str)
    succeeded = Signal(object)
    failed = Signal(str, str)
    wasCancelled = Signal()

    def __init__(self, fn: Callable, *args, parent=None, **kwargs) -> None:
        super().__init__(parent)
        self._fn = fn
        self._args = args
        self._kwargs = kwargs
        self._cancelled = False
        self._finished_normally = False
        sig = None
        try:
            sig = inspect.signature(fn)
        except (TypeError, ValueError):
            sig = None
        self._accepts_progress = sig is None or "progress_cb" in sig.parameters
        self._accepts_cancel = sig is None or "cancel_check" in sig.parameters

    # -- 控制 -------------------------------------------------------------- #
    def cancel(self) -> None:
        self._cancelled = True

    @property
    def is_cancelled(self) -> bool:
        return self._cancelled

    # -- 内部 -------------------------------------------------------------- #
    def _emit(self, percent: int, message: str) -> None:
        if not self._cancelled:
            self.progressed.emit(int(percent), str(message))

    def run(self) -> None:  # noqa: D401 - QThread 接口
        try:
            kwargs = dict(self._kwargs)
            if self._accepts_progress:
                kwargs["progress_cb"] = self._emit
            if self._accepts_cancel:
                kwargs["cancel_check"] = lambda: self._cancelled
            result = self._fn(*self._args, **kwargs)
        except CancelledError:
            self.wasCancelled.emit()
        except GGEError as exc:
            self.failed.emit(exc.message, exc.hint or "")
        except Exception as exc:  # pragma: no cover - 兜底
            self.failed.emit(f"{type(exc).__name__}: {exc}",
                             traceback.format_exc(limit=12))
        else:
            self._finished_normally = True
            self.succeeded.emit(result)


def wait_for(job: Job, timeout_ms: int = 120000) -> bool:
    """测试用：等待任务结束。"""
    return job.wait(timeout_ms)
