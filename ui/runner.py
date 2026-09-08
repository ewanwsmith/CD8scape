"""
runner.py — cross-platform subprocess launcher for CD8scape.

Design goals (Phase 1):
    * No shell invocation (shell=False) — avoids shell-injection and
      shell-specific quoting/escaping differences between bash/zsh/cmd/PowerShell.
    * Locate the Julia interpreter portably via shutil.which, so the shebang
      line of CD8scape.jl is irrelevant. Julia's location is looked up the
      same way on macOS, Linux, and Windows.
    * Paths are built with pathlib.Path so separators are never hard-coded.
    * stdout and stderr are merged and streamed line-by-line to a caller-supplied
      callback, so the UI can render them live.
    * The function is a generator: the UI stays responsive and we never buffer
      the whole log in memory before showing it.

This module has NO Streamlit dependency — it can be unit-tested and reused
from any front-end (CLI, Streamlit, Flask, Tk, etc.).
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, List, Optional  # Optional kept for RunResult


# ---------------------------------------------------------------------------
# Repo root resolution
# ---------------------------------------------------------------------------

REPO_ROOT: Path = Path(__file__).resolve().parent.parent
CD8SCAPE_SCRIPT: Path = REPO_ROOT / "CD8scape.jl"

# ---------------------------------------------------------------------------
# Pinned Julia version
# ---------------------------------------------------------------------------
# CD8scape's Manifest.toml is resolved for the Julia 1.11 LTS series. Newer
# Julia releases (1.12+, and the 1.13 pre-releases) removed internal symbols
# (e.g. Base.memhash_seed) that the pinned dependency versions still call, so
# running the pipeline on them crashes inside CSV.jl/InlineStrings. Until the
# dependencies are upgraded, we pin every invocation to 1.11.
#
# When Julia is managed by juliaup (the usual case), the "+<channel>" selector
# forces that channel per-invocation without touching the user's global
# default. Non-juliaup installs don't understand "+<channel>", so we only add
# it when juliaup is actually present; julia_launch_prefix() raises a clear
# error if a bare install is on the wrong version.
#
# CD8scape.jl propagates this pin to its worker subprocesses via
# Sys.BINDIR + Base.julia_exename(), so the whole pipeline stays on 1.11.
JULIA_CHANNEL: str = "1.11"


def get_repo_root() -> Path:
    """Return the CD8scape repository root (parent of ui/)."""
    return REPO_ROOT


@dataclass
class RunResult:
    """Structured summary of a completed CD8scape invocation."""

    command: List[str]
    returncode: int
    success: bool


class CD8scapeNotFoundError(RuntimeError):
    """Raised when CD8scape.jl cannot be located."""


class JuliaNotFoundError(RuntimeError):
    """Raised when the Julia interpreter cannot be located on PATH."""


def find_julia() -> str:
    """Locate the Julia executable in a cross-platform way.

    Returns the absolute path as a string (suitable for subprocess).
    Raises JuliaNotFoundError if Julia is not on PATH.

    Notes:
        * shutil.which respects PATH and PATHEXT, so it finds ``julia.exe``
          on Windows and ``julia`` on macOS/Linux without any OS-specific code.
        * We do NOT rely on the shebang line in CD8scape.jl, because Windows
          does not honour shebangs.
        * The launcher (.command / .bat) inherits the user's shell PATH, so
          juliaup, Homebrew, and conda installations are all visible here.
    """
    julia = shutil.which("julia")
    if julia is None:
        raise JuliaNotFoundError(
            "The 'julia' executable was not found on your PATH.\n\n"
            "Install Julia 1.11 from https://julialang.org/downloads/ "
            "(or run `juliaup add 1.11`) and make sure the 'julia' command "
            "works in a new terminal."
        )
    return julia


def _juliaup_available() -> bool:
    """Return True if the juliaup version multiplexer is on PATH.

    Only juliaup understands the ``julia +<channel>`` selector we use to pin
    the version. A direct (non-juliaup) install would treat ``+1.11`` as a
    script path and fail, so we must not pass it in that case.
    """
    return shutil.which("juliaup") is not None


def julia_launch_prefix() -> List[str]:
    """Return the argv prefix used to launch the pinned Julia version.

    * With juliaup: ``[julia, "+1.11"]`` — forces the 1.11 channel for this
      invocation only, regardless of the user's global default channel.
    * Without juliaup: ``[julia]`` — plus a best-effort version check that
      raises a clear error if a bare install is on an unsupported version.
    """
    julia = find_julia()
    if _juliaup_available():
        return [julia, f"+{JULIA_CHANNEL}"]

    # Bare install: we can't switch versions, so verify the one we have.
    try:
        out = subprocess.run(
            [julia, "--version"],
            capture_output=True,
            text=True,
            timeout=30,
        ).stdout
    except Exception:
        # If the check itself fails, don't block the run — let Julia speak.
        return [julia]

    # Output looks like "julia version 1.11.9".
    version = out.strip().split()[-1] if out.strip() else ""
    if not version.startswith(JULIA_CHANNEL + "."):
        raise JuliaNotFoundError(
            f"CD8scape is pinned to the Julia {JULIA_CHANNEL} LTS series, but "
            f"the 'julia' on your PATH is version {version or 'unknown'}.\n\n"
            f"Newer Julia versions crash the pipeline with the current pinned "
            f"dependencies. Install Julia {JULIA_CHANNEL} from "
            f"https://julialang.org/downloads/, or install juliaup and run "
            f"`juliaup add {JULIA_CHANNEL}`."
        )
    return [julia]


def check_cd8scape_script() -> Path:
    """Confirm that CD8scape.jl exists in the configured repo root."""
    root = get_repo_root()
    script = root / "CD8scape.jl"
    if not script.is_file():
        raise CD8scapeNotFoundError(
            f"Could not find CD8scape.jl at {script}.\n\n"
            "Please check the CD8scape Repository path in Setup."
        )
    return script


def build_command(cd8scape_args: List[str]) -> List[str]:
    """Build the argv list for a CD8scape invocation.

    Using a list (never a string) means subprocess hands the arguments to the
    OS verbatim — no shell parsing, no cross-platform quoting surprises.
    """
    prefix = julia_launch_prefix()
    script = check_cd8scape_script()
    return [*prefix, str(script), *cd8scape_args]


def stream_cd8scape(cd8scape_args: List[str]) -> Iterator[str]:
    """Run CD8scape and yield merged stdout/stderr lines as they are produced.

    The final yielded item is a sentinel dict-like string of the form
    ``"__CD8SCAPE_EXIT__:<returncode>"`` so callers can detect completion and
    the exit status from a single stream.

    Cross-platform notes:
        * ``shell=False`` (default when argv is a list) avoids every known
          shell-quoting pitfall on Windows, macOS, and Linux.
        * ``text=True`` + ``encoding='utf-8'`` + ``errors='replace'`` keeps
          output decoding deterministic regardless of the user's locale.
        * ``bufsize=1`` (line-buffered) plus merging stderr into stdout
          ensures we see logs in the order they were emitted.
        * ``cwd=REPO_ROOT`` means CD8scape's relative include()s work
          exactly as they do from the command line.
    """
    argv = build_command(cd8scape_args)

    # Merge stderr into stdout so ordering is preserved and we only
    # need to poll one pipe.
    proc = subprocess.Popen(
        argv,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        cwd=str(get_repo_root()),
        text=True,
        encoding="utf-8",
        errors="replace",
        bufsize=1,
        # env=os.environ.copy() lets CD8SCAPE_* env vars flow through if set.
        env=os.environ.copy(),
    )

    # Stream lines as they arrive. readline() returns '' only at EOF.
    if proc.stdout is None:
        raise RuntimeError(
            "subprocess stdout is None — was Popen called without stdout=PIPE?"
        )
    try:
        for line in proc.stdout:
            # Preserve line content; strip only the trailing newline so the UI
            # can join with its own separator.
            yield line.rstrip("\r\n")
    finally:
        proc.stdout.close()
        returncode = proc.wait()

    yield f"__CD8SCAPE_EXIT__:{returncode}"


def run_cd8scape(cd8scape_args: List[str]) -> RunResult:
    """Blocking helper that collects the stream and returns a RunResult.

    Mostly useful for tests and the CLI smoke-test (``python runner.py``).
    The Streamlit UI uses :func:`stream_cd8scape` directly so it can render
    lines live.
    """
    lines: List[str] = []
    returncode: Optional[int] = None
    for item in stream_cd8scape(cd8scape_args):
        if item.startswith("__CD8SCAPE_EXIT__:"):
            returncode = int(item.split(":", 1)[1])
        else:
            lines.append(item)
    if returncode is None:
        raise RuntimeError(
            "stream_cd8scape() ended without emitting an __CD8SCAPE_EXIT__ sentinel"
        )
    print("\n".join(lines))
    return RunResult(
        command=build_command(cd8scape_args),
        returncode=returncode,
        success=returncode == 0,
    )


if __name__ == "__main__":
    # Tiny CLI smoke test:
    #     python ui/runner.py            -> runs `CD8scape.jl --help`
    #     python ui/runner.py read data/Example_data --aa
    args = sys.argv[1:] or ["--help"]
    try:
        result = run_cd8scape(args)
    except (JuliaNotFoundError, CD8scapeNotFoundError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
    sys.exit(0 if result.success else 1)
