"""Single source of truth: VERSION.txt at the repo root.

`pyproject.toml` reads VERSION.txt at build time, so once the package is
installed (editable or wheel) `importlib.metadata.version("irfinder-mdl")`
returns the same string.  The `VERSION.txt`-on-disk fallback below covers the
"run straight out of a checkout without pip install" case.
"""
from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version as _pkg_version
from pathlib import Path


def _read_version_file() -> str:
    # VERSION.txt sits two levels up from this file: <repo>/VERSION.txt
    candidate = Path(__file__).resolve().parent.parent / "VERSION.txt"
    if candidate.is_file():
        return candidate.read_text(encoding="utf-8").strip()
    return "0+unknown"


try:
    __version__ = _pkg_version("irfinder-mdl")
except PackageNotFoundError:
    __version__ = _read_version_file()
