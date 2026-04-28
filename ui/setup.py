"""
setup.py — Environment configuration helpers for the CD8scape UI.

Handles everything the user would otherwise need to do manually in a text
editor or terminal before running CD8scape:

    1. Reading the current netMHCpan path from src/settings.txt
    2. Validating that path (exists on disk, is an executable file)
    3. Writing an updated path back to settings.txt (preserving comments)
    4. Checking that Perl is available on PATH

No Streamlit dependency — this module is pure Python / stdlib so it can be
tested without a running app.

settings.txt format (either form is accepted by env.jl):
    NETMHCPAN=/full/path/to/netMHCpan      ← key=value (preferred)
    /full/path/to/netMHCpan                 ← bare path (legacy fallback)
    # comment lines and blank lines are ignored
"""

from __future__ import annotations

import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

# Absolute path to the settings file — one level up from ui/, into src/.
SETTINGS_PATH: Path = Path(__file__).resolve().parent.parent / "src" / "settings.txt"

# Strings that identify a placeholder / unconfigured path.
_PLACEHOLDER_MARKERS = ("full/path", "<", ">")


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------


def read_netmhcpan_path() -> Optional[str]:
    """Return the NETMHCPAN value from settings.txt.

    Returns None when the file does not exist, the key is absent, or the
    value still looks like the example placeholder.
    """
    if not SETTINGS_PATH.exists():
        return None

    for raw in SETTINGS_PATH.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue

        if "=" in line:
            key, _, value = line.partition("=")
            if key.strip().upper() == "NETMHCPAN":
                value = value.strip()
                if value and not any(m in value for m in _PLACEHOLDER_MARKERS):
                    return value
        else:
            # Bare-path form — treat the whole non-comment line as the path.
            if not any(m in line for m in _PLACEHOLDER_MARKERS):
                return line

    return None


# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------


def write_netmhcpan_path(new_path: str) -> None:
    """Update the NETMHCPAN= line in settings.txt.

    Behaviour:
        * If settings.txt already contains a NETMHCPAN= line (or a bare
          non-comment path line), replace it in-place so that comments and
          other keys are preserved.
        * If no such line exists, append one.
        * If the file does not exist, create it with just the key=value line.

    The written value is always in ``NETMHCPAN=<path>`` form.
    """
    new_line = f"NETMHCPAN={new_path}"

    if not SETTINGS_PATH.exists():
        SETTINGS_PATH.write_text(new_line + "\n", encoding="utf-8")
        return

    lines = SETTINGS_PATH.read_text(encoding="utf-8").splitlines()
    result = []
    replaced = False

    for raw in lines:
        stripped = raw.strip()
        if not stripped or stripped.startswith("#"):
            result.append(raw)
            continue

        if "=" in stripped:
            key, _, _ = stripped.partition("=")
            if key.strip().upper() == "NETMHCPAN":
                result.append(new_line)
                replaced = True
                continue
        else:
            # Bare-path line — replace it too.
            result.append(new_line)
            replaced = True
            continue

        result.append(raw)

    if not replaced:
        result.append(new_line)

    SETTINGS_PATH.write_text("\n".join(result) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


@dataclass
class PathCheck:
    ok: bool
    message: str
    resolved: Optional[Path] = None


def validate_netmhcpan(raw_path: str) -> PathCheck:
    """Check whether *raw_path* points to an executable netMHCpan file."""
    if not raw_path or not raw_path.strip():
        return PathCheck(ok=False, message="No path provided.")

    p = Path(raw_path.strip()).expanduser()
    if not p.is_absolute():
        p = (SETTINGS_PATH.parent / p).resolve()
    else:
        p = p.resolve()

    if not p.exists():
        return PathCheck(
            ok=False,
            message=f"Not found: {p}\nDouble-check the path to your netMHCpan installation.",
            resolved=p,
        )
    if not p.is_file():
        return PathCheck(
            ok=False,
            message=f"{p} is a directory, not a file. Select the `netMHCpan` executable itself.",
            resolved=p,
        )
    if not os.access(p, os.X_OK):
        return PathCheck(
            ok=False,
            message=f"{p} is not executable.\nRun: chmod +x {p}",
            resolved=p,
        )

    return PathCheck(ok=True, message=f"Executable found: {p}", resolved=p)


@dataclass
class EnvStatus:
    """Snapshot of the CD8scape environment readiness."""

    netmhcpan_path: Optional[str]
    netmhcpan_ok: bool
    netmhcpan_msg: str
    perl_ok: bool
    perl_path: Optional[str]

    @property
    def ready(self) -> bool:
        return self.netmhcpan_ok and self.perl_ok


def check_env() -> EnvStatus:
    """Return a complete snapshot of the environment status."""
    raw_path = read_netmhcpan_path()
    if raw_path:
        pcheck = validate_netmhcpan(raw_path)
        nmhc_ok = pcheck.ok
        nmhc_msg = pcheck.message
    else:
        nmhc_ok = False
        nmhc_msg = "netMHCpan path not configured in settings.txt."

    perl_exe = shutil.which("perl")

    return EnvStatus(
        netmhcpan_path=raw_path,
        netmhcpan_ok=nmhc_ok,
        netmhcpan_msg=nmhc_msg,
        perl_ok=perl_exe is not None,
        perl_path=perl_exe,
    )
