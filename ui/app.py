"""
app.py — Streamlit UI for CD8scape.

The UI is a *wrapper*. It renders CLI arguments as UI controls, spawns
`julia CD8scape.jl …` via a cross-platform subprocess, and streams the
output. It deliberately knows nothing about CD8scape's input files,
output files, or pipeline logic — if CD8scape changes its internals,
this UI keeps working unchanged.

Current scope (Phases 1–3):
    * Select a CD8scape subcommand from a dropdown.
    * Pick a data folder (basic filesystem check only). The picker is
      hidden for `prep`, which takes no folder.
    * Configure CLI flags through controls that are a 1:1 mirror of
      CD8scape's own --help output.
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
from commands import (  # noqa: E402
    COMMAND_HELP,
    COMMANDS,
    needs_folder,
    render_options,
)


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
    st.session_state.setdefault("command", "read")
    st.session_state.setdefault("data_folder", DEFAULT_DATA_FOLDER)
    st.session_state.setdefault("folder_check", None)     # FolderCheck | None
    st.session_state.setdefault("log_lines", [])
    st.session_state.setdefault("last_status", None)      # "success"|"failure"|"error"|None
    st.session_state.setdefault("last_message", "")


_init_state()


# -- Helpers ------------------------------------------------------------------


def build_cd8scape_args(
    command: str,
    folder: Path | None,
    option_args: List[str],
) -> List[str]:
    """Build the full CD8scape argv tail: subcommand + folder (optional) + options."""
    argv: List[str] = [command]
    if folder is not None:
        argv.append(str(folder))
    argv.extend(option_args)
    return argv


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


def render_sidebar(preview_args: List[str]) -> None:
    with st.sidebar:
        st.subheader("Environment")
        st.text(f"OS:     {platform.system()} {platform.release()}")
        st.text(f"Python: {platform.python_version()}")
        st.text(f"Repo:   {REPO_ROOT}")

        st.subheader("Next command")
        if not preview_args:
            st.caption("Configure the options to preview the full command.")
            return
        try:
            argv = build_command(preview_args)
            st.code("\n".join(argv), language="text")
            st.caption(
                "This is the exact command the Run button will execute. "
                "No shell is involved."
            )
        except (JuliaNotFoundError, CD8scapeNotFoundError) as exc:
            st.warning(str(exc))


# -- Main: command picker -----------------------------------------------------

st.subheader("1. Choose a CD8scape command")

command = st.selectbox(
    "Subcommand",
    options=COMMANDS,
    index=COMMANDS.index(st.session_state.command),
    key="command",
    help=(
        "Direct mirror of the subcommands listed by `./CD8scape.jl --help`. "
        "Each one maps to a command you'd normally type in a terminal."
    ),
)
st.caption(COMMAND_HELP[command])


# -- Main: folder picker (skipped for prep) -----------------------------------

folder_check: FolderCheck | None = None
if needs_folder(command):
    st.subheader("2. Choose your data folder")
    st.caption(
        "Enter the path to a folder containing your CD8scape inputs. "
        "The UI only checks that the path exists and is a folder; "
        "CD8scape decides which files inside it to use."
    )

    folder_text = st.text_input(
        "Data folder",
        value=st.session_state.data_folder,
        key="folder_text_input",
        help=(
            "Tip: drag a folder from Finder / File Explorer into this field "
            "to paste its path."
        ),
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

    folder_check = st.session_state.folder_check
    render_folder_feedback(folder_check)
else:
    # `prep` has no folder; clear any stale check so it doesn't influence preview.
    st.info("`prep` doesn't take a folder — skip straight to Run.")


# -- Main: options ------------------------------------------------------------

step_no = 3 if needs_folder(command) else 2
st.subheader(f"{step_no}. Options")
st.caption(
    "Each control here is a direct mirror of a CD8scape CLI flag. "
    "The UI does not choose flags for you."
)
option_args = render_options(command)


# -- Main: run ---------------------------------------------------------------

run_step = step_no + 1
st.subheader(f"{run_step}. Run CD8scape")

status_slot = st.empty()
log_slot = st.empty()
render_status(status_slot)
render_logs(log_slot)

# Build the preview argv (after options have been rendered) so the sidebar
# shows the live command. For non-prep commands we only include the folder
# if it checked out OK; otherwise we leave it absent and show "configure..."
preview_folder = folder_check.path if (folder_check is not None and folder_check.ok) else None
can_run_folder = (not needs_folder(command)) or (folder_check is not None and folder_check.ok)

if can_run_folder:
    preview_args = build_cd8scape_args(command, preview_folder, option_args)
else:
    preview_args = []

render_sidebar(preview_args)

run_col, clear_col = st.columns([1, 1])
run_clicked = run_col.button(
    "Run CD8scape",
    type="primary",
    use_container_width=True,
    disabled=not can_run_folder,
    help=(
        None
        if can_run_folder
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

if run_clicked and can_run_folder:
    # For folder-requiring commands, re-check the path at launch time so a
    # stale green tick from before editing the field doesn't slip through.
    live_folder: Path | None = None
    if needs_folder(command):
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
        live_folder = live_check.path

    args = build_cd8scape_args(command, live_folder, option_args)

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
        if live_folder is not None:
            trailer = (
                f"Any output files it produced are inside `{live_folder}` — "
                "see CD8scape's own docs for what to expect."
            )
        else:
            trailer = "See CD8scape's own docs for what this command does."
        st.session_state.last_status = "success"
        st.session_state.last_message = f"CD8scape finished successfully. {trailer}"
    else:
        st.session_state.last_status = "failure"
        st.session_state.last_message = (
            f"CD8scape exited with a non-zero status (exit code {exit_code}). "
            "The log above is CD8scape's own output — use it to diagnose the run."
        )
    render_status(status_slot)
    render_logs(log_slot)
