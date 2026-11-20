"""Small per-user memory of the interface: last mode and folders used.
Stored in <home>/.msparser/state.json (or in $MSPARSER_HOME)."""
from __future__ import annotations

import json
import os
from pathlib import Path


def _file() -> Path:
    return Path(os.environ.get("MSPARSER_HOME") or Path.home() / ".msparser") / "state.json"


def load() -> dict:
    try:
        return json.loads(_file().read_text())
    except (OSError, ValueError):
        return {}


def save(state: dict) -> None:
    try:
        _file().parent.mkdir(parents=True, exist_ok=True)
        _file().write_text(json.dumps(state, indent=2))
    except OSError:
        pass  # the memory is a convenience only


def remember_folder(state: dict, folder: str, key: str = "recent_folders", keep: int = 8) -> None:
    recent = [f for f in state.get(key, []) if f != folder]
    state[key] = [folder] + recent[: keep - 1]
