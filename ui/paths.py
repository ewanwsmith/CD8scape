"""
paths.py — UI-only helpers for the folder picker.

Intentionally tiny. The UI is a *wrapper* around CD8scape; it must not
reimplement pipeline logic. In particular we deliberately do NOT:

    * detect variant file types (.vcf / .aa / .out) — that's CD8scape's job
    * decide when to pass --aa — that's the user's call
    * interpret reading-frame files — that's CD8scape's job

All we answer here are plain UI questions: is the path the user typed a
real folder on their machine, and (optionally) what files are in it so we
can show a read-only listing.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple


@dataclass
class FolderCheck:
    """Result of a basic filesystem sanity check on a user-typed path."""

    path: Path
    ok: bool
    message: str = ""


def resolve_folder(raw_path: str) -> FolderCheck:
    """Expand, resolve, and sanity-check a path the user typed.

    The only questions we answer are:
        - was the path non-empty?
        - does it exist?
        - is it a directory?

    Anything about the folder's *contents* is left to CD8scape.
    """
    stripped = (raw_path or "").strip().strip('"').strip("'")
    if not stripped:
        return FolderCheck(path=Path("."), ok=False, message="Please enter a folder path.")

    path = Path(stripped).expanduser()
    try:
        path = path.resolve()
    except OSError:
        # Extremely unusual — e.g. filesystem loops — but we still want to
        # fail clearly rather than explode.
        return FolderCheck(
            path=path,
            ok=False,
            message=f"Could not resolve '{stripped}'. Check the path and try again.",
        )

    if not path.exists():
        return FolderCheck(
            path=path,
            ok=False,
            message=f"The folder '{path}' doesn't exist. Check the path and try again.",
        )
    if not path.is_dir():
        return FolderCheck(
            path=path,
            ok=False,
            message=(
                f"'{path}' is a file, not a folder. Pick the folder that contains "
                "your CD8scape inputs."
            ),
        )

    return FolderCheck(path=path, ok=True, message="")


def list_folder(folder: Path, limit: Optional[int] = 50) -> Tuple[List[str], int]:
    """Return a plain list of entries inside *folder* for read-only display.

    We sort alphabetically (folders first, then files) so the listing is
    stable across runs. We do NOT interpret filenames — that's on CD8scape.

    Returns (entries, total) where `entries` is capped at `limit` and
    `total` is the full count (so the UI can say "+N more").
    """
    try:
        raw = list(folder.iterdir())
    except OSError:
        return [], 0

    def key(p: Path) -> Tuple[int, str]:
        # 0 for directories, 1 for files — puts folders first.
        return (0 if p.is_dir() else 1, p.name.lower())

    raw.sort(key=key)
    names = [f"{p.name}/" if p.is_dir() else p.name for p in raw]
    total = len(names)
    if limit is not None and total > limit:
        names = names[:limit]
    return names, total


# -- CLI smoke test -----------------------------------------------------------

if __name__ == "__main__":
    import sys as _sys

    paths = _sys.argv[1:] or ["."]
    bad = 0
    for raw in paths:
        check = resolve_folder(raw)
        print(f"=== {raw} ===")
        print(f"  resolved: {check.path}")
        print(f"  ok:       {check.ok}")
        if check.message:
            print(f"  msg:      {check.message}")
        if check.ok:
            entries, total = list_folder(check.path)
            print(f"  entries ({total}):")
            for name in entries:
                print(f"    {name}")
        if not check.ok:
            bad += 1
    _sys.exit(1 if bad else 0)
