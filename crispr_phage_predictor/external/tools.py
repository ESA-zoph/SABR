from __future__ import annotations

import shutil


def tool_available(command: str) -> bool:
    return shutil.which(command) is not None


def missing_tool_message(command: str) -> str:
    return (
        f"External tool '{command}' was not found on PATH. Install it or select the "
        "internal Python backend."
    )

