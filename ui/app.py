"""
app.py — Streamlit UI for CD8scape.

This is a *wrapper*. It renders CLI arguments as UI controls, spawns
`julia CD8scape.jl …` via a cross-platform subprocess, and streams the
output. It must not know anything about CD8scape's input files, output
files, or pipeline logic — if CD8scape's behaviour changes, this UI
keeps working unchanged.

Current scope:
    * Pick a data folder (basic existence check only).
    * Choose the CD8scape subcommand (Phase 2 ships with `read`).
    * Optional --aa toggle — a direct mirror of the CLI flag.
    * Run; stream merged stdout/stderr live; report the exit code.
"""

from __future__ import annotations

import platform
import sys
from pathlib import Path
from typing import List

import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent))

from runner import (  # noqa: E402
    CD8scapeNotFoundError,
    JuliaNotFoundError,
    REPO_ROOT,
    build_command,
    stream_cd8scape,
)
from paths import FolderCheck, list_folder, resolve_folder  # noqa: E402


# -- Page setup ---------------------------------------------------------------

st.set_page_config(page_title="CD8scape", layout="centered")
st.title("CD8scape")
st.caption(
    "Launcher for the CD8scape pipeline. "
    "This UI is a thin wrapper — CD8scape does all the work."
)


DEFAULT_DATA_FOLDER = str(REPO_ROOT / "data" / "Example_data")


# -- Session state ------------------------------------------------------------

def _init_state() -> None:
    st.session_state.setdefault("data_folder", DEFAULT_DATA_FOLDER)
    st.session_state.setdefault("folder_check", None)  # FolderCheck | None
    st.session_state.setdefault("use_aa", False)
    st.session_state.setdefault("log_lines", [])
    st.session_state.setdefault("last_status", None)   # "success"|"failure"|"error"|None
    st.session_state.setdefault("last_message", "")


_init_state()


# -- Helpers ------------------------------------------------------------------


def build_args(folder: Path, use_aa: bool) -> List[str]:
    """Build the CD8scape argv tail. Direct, 1:1 mapping from UI to CLI."""
    args: List[str] = ["read", str(folder)]
    if use_aa:
        args.append("--aa")
    return args


def render_folder_feedback(check: FolderCheck | None) -> None:
    if check is None:
        st.info("Enter a folder above and click **Check folder**.")
        return
    if check.ok:
        st.success(f"Folder found: `{check.path}`")
        entries, total = list_folder(check.path)
        with st.expander(f"What's in this folder ({total} items)?"):
            if not entries:
                st.caption("This folder is empty.")
            else:
                st.code("\n".join(entries), language="text")
                if total > len(entries):
                    st.caption(f"(showing {len(entries)} of {total})")
        st.caption(
            "The UI does not inspect what these files are — CD8scape decides "
            "how to handle them when you click Run."
        )
    else:
        st.error(check.message)


def render_logs(slot) -> None:
    log_text = "\n".join(st.session_state.log_lines) or "(no output yet)"
    slot.code(log_text, language="text")


def render_status(slot) -> None:
    status = st.session_state.last_status
    msg = st.session_state.last_message
    if status == "success":
        slot.success(msg)
    elif status == "failure":
        slot.error(msg)
    elif status == "error":
        slot.warning(msg)
    else:
        slot.empty()


# -- Sidebar ------------------------------------------------------------------


def render_sidebar() -> None:
    with st.sidebar:
        st.subheader("Environment")
        st.text(f"OS:     {platform.system()} {platform.release()}")
        st.text(f"Python: {platform.python_version()}")
        st.text(f"Repo:   {REPO_ROOT}")

        st.subheader("Next command")
        check: FolderCheck | None = st.session_state.folder_check
        if check is None or not check.ok:
            st.caption("Pick a folder to preview the full command.")
            return
        try:
            argv = build_command(build_args(check.path, st.session_state.use_aa))
            st.code("\n".join(argv), language="text")
            st.caption(
                "This is the exact command the Run button will execute. "
                "No shell is involved."
            )
        except (JuliaNotFoundError, CD8scapeNotFoundError) as exc:
            st.warning(str(exc))


render_sidebar()


# -- Main: folder picker ------------------------------------------------------

st.subheader("1. Choose your data folder")
st.caption(
    "Enter the path to a folder containing your CD8scape inputs. "
    "The UI only checks that the path exists and is a folder; "
    "CD8scape decides which files inside it to use."
)

folder_text = st.text_input(
    "Data folder",
    value=st.session_state.data_folder,
    key="folder_text_input",
    help="Tip: drag a folder from Finder / File Explorer into this field to paste its path.",
)
st.session_state.data_folder = folder_text

check_col, reset_col = st.columns([1, 1])
check_clicked = check_col.button("Check folder", use_container_width=True)
reset_clicked = reset_col.button("Reset to example", use_container_width=True)

if reset_clicked:
    st.session_state.data_folder = DEFAULT_DATA_FOLDER
    st.session_state.folder_check = None
    st.session_state.log_lines = []
    st.session_state.last_status = None
    st.session_state.last_message = ""
    st.rerun()

if check_clicked:
    st.session_state.folder_check = resolve_folder(st.session_state.data_folder)

render_folder_feedback(st.session_state.folder_check)


# -- Main: CLI options --------------------------------------------------------

st.subheader("2. Options")
st.caption(
    "Each control here is a direct mirror of a CD8scape CLI flag. "
    "The UI does not choose flags for you."
)

st.session_state.use_aa = st.checkbox(
    "Pass `--aa` (variants are amino-acid format, .aa file)",
    value=st.session_state.use_aa,
    help=(
        "Tick this only if your variant file is a .aa amino-acid file. "
        "CD8scape documents this flag — see ./CD8scape.jl --help."
    ),
)


# -- Main: run ---------------------------------------------------------------

st.subheader("3. Run CD8scape")

status_slot = st.empty()
log_slot = st.empty()
render_status(status_slot)
render_logs(log_slot)

check: FolderCheck | None = st.session_state.folder_check
can_run = check is not None and check.ok

run_col, clear_col = st.columns([1, 1])
run_clicked = run_col.button(
    "Run CD8scape",
    type="primary",
    use_container_width=True,
    disabled=not can_run,
    help=(
        None
        if can_run
        else "Click 'Check folder' first. The Run button unlocks once the folder path is valid."
    ),
)
clear_clicked = clear_col.button("Clear logs", use_container_width=True)

if clear_clicked:
    st.session_state.log_lines = []
    st.session_state.last_status = None
    st.session_state.last_message = ""
    render_status(status_slot)
    render_logs(log_slot)

if run_clicked and can_run:
    # Re-check at launch time in case the path was edited after the last check.
    live_check = resolve_folder(st.session_state.data_folder)
    st.session_state.folder_check = live_check
    if not live_check.ok:
        render_folder_feedback(live_check)
        st.session_state.last_status = "error"
        st.session_state.last_message = (
            "The folder path is no longer valid. Re-check it before running."
        )
        render_status(status_slot)
        st.stop()

    args = build_args(live_check.path, st.session_state.use_aa)
    st.session_state.log_lines = []
    st.session_state.last_status = None
    st.session_state.last_message = ""
    render_status(status_slot)
    render_logs(log_slot)

    try:
        stream = stream_cd8scape(args)
    except (JuliaNotFoundError, CD8scapeNotFoundError) as exc:
        st.session_state.last_status = "error"
        st.session_state.last_message = str(exc)
        render_status(status_slot)
        st.stop()

    exit_code: int | None = None
    with st.spinner(f"Running `CD8scape.jl {' '.join(args)}` …"):
        for item in stream:
            if item.startswith("__CD8SCAPE_EXIT__:"):
                exit_code = int(item.split(":", 1)[1])
                break
            st.session_state.log_lines.append(item)
            render_logs(log_slot)

    if exit_code == 0:
        st.session_state.last_status = "success"
        st.session_state.last_message = (
            "CD8scape finished successfully. Any output files it produced are "
            f"inside `{live_check.path}` — see CD8scape's own docs for what to expect."
        )
    else:
        st.session_state.last_status = "failure"
        st.session_state.last_message = (
            f"CD8scape exited with a non-zero status (exit code {exit_code}). "
            "The log above is CD8scape's own output — use it to diagnose the run."
        )
    render_status(status_slot)
    render_logs(log_slot)
