"""
app.py — Phase 1 Streamlit launcher for CD8scape.

Scope (deliberately tiny):
    * One button: "Run CD8scape"
    * One hardcoded invocation: ``CD8scape.jl --help``
    * Live stdout/stderr stream
    * Success / failure banner at the end

Phase 1 proves the plumbing — nothing else. File pickers, parameters,
results viewers, and safety rails arrive in later phases.

Run:
    streamlit run ui/app.py

Cross-platform notes are in the docstring of runner.stream_cd8scape.
"""

from __future__ import annotations

import platform
import sys
from pathlib import Path

import streamlit as st

# Ensure this folder is importable when Streamlit is launched from the
# repo root (`streamlit run ui/app.py`).
sys.path.insert(0, str(Path(__file__).resolve().parent))

from runner import (  # noqa: E402  (import after sys.path mutation)
    CD8scapeNotFoundError,
    JuliaNotFoundError,
    REPO_ROOT,
    build_command,
    stream_cd8scape,
)


# -- Page setup ---------------------------------------------------------------

st.set_page_config(
    page_title="CD8scape",
    page_icon=None,
    layout="centered",
)

st.title("CD8scape")
st.caption(
    "A friendly launcher for the CD8scape pipeline. "
    "Phase 1: smoke-test the pipeline with one click."
)


# -- Hardcoded Phase 1 test input --------------------------------------------
#
# We invoke `CD8scape.jl --help` because it exercises the full call path
# (find julia → spawn subprocess → stream output → read exit code) without
# requiring netMHCpan, Julia packages, or any user data. Later phases will
# let the user choose the command and pass their own inputs.
PHASE1_ARGS = ["--help"]


# -- Session state ------------------------------------------------------------

if "log_lines" not in st.session_state:
    st.session_state.log_lines = []
if "last_status" not in st.session_state:
    # None | "success" | "failure" | "error"
    st.session_state.last_status = None
if "last_message" not in st.session_state:
    st.session_state.last_message = ""


# -- Sidebar: environment snapshot (useful for diagnosing cross-platform issues)

with st.sidebar:
    st.subheader("Environment")
    st.text(f"OS:      {platform.system()} {platform.release()}")
    st.text(f"Python:  {platform.python_version()}")
    st.text(f"Repo:    {REPO_ROOT}")
    try:
        preview_cmd = build_command(PHASE1_ARGS)
        st.text("Command:")
        # Show each argv entry on its own line — makes it obvious that we
        # are NOT building a shell string.
        st.code("\n".join(preview_cmd), language="text")
    except (JuliaNotFoundError, CD8scapeNotFoundError) as exc:
        st.warning(str(exc))


# -- Main controls ------------------------------------------------------------

col_run, col_clear = st.columns([1, 1])
run_clicked = col_run.button("Run CD8scape", type="primary", use_container_width=True)
clear_clicked = col_clear.button("Clear logs", use_container_width=True)

if clear_clicked:
    st.session_state.log_lines = []
    st.session_state.last_status = None
    st.session_state.last_message = ""

status_slot = st.empty()
log_slot = st.empty()


def render_logs() -> None:
    """Render the full log buffer in a fixed-height code block."""
    log_text = "\n".join(st.session_state.log_lines) or "(no output yet)"
    log_slot.code(log_text, language="text")


def render_status() -> None:
    """Render the most recent run status, if any."""
    status = st.session_state.last_status
    msg = st.session_state.last_message
    if status == "success":
        status_slot.success(msg)
    elif status == "failure":
        status_slot.error(msg)
    elif status == "error":
        status_slot.warning(msg)
    else:
        status_slot.empty()


# Render existing state first so "Clear" and page reruns look right.
render_status()
render_logs()


# -- Run action ---------------------------------------------------------------

if run_clicked:
    st.session_state.log_lines = []
    st.session_state.last_status = None
    st.session_state.last_message = ""
    render_status()
    render_logs()

    try:
        stream = stream_cd8scape(PHASE1_ARGS)
    except (JuliaNotFoundError, CD8scapeNotFoundError) as exc:
        st.session_state.last_status = "error"
        st.session_state.last_message = str(exc)
        render_status()
        st.stop()

    exit_code: int | None = None
    # Iterate the generator. stream_cd8scape yields one line per read, and a
    # final sentinel of the form "__CD8SCAPE_EXIT__:<code>".
    with st.spinner("Running CD8scape…"):
        for item in stream:
            if item.startswith("__CD8SCAPE_EXIT__:"):
                exit_code = int(item.split(":", 1)[1])
                break
            st.session_state.log_lines.append(item)
            render_logs()

    if exit_code == 0:
        st.session_state.last_status = "success"
        st.session_state.last_message = "CD8scape finished successfully (exit code 0)."
    else:
        st.session_state.last_status = "failure"
        st.session_state.last_message = (
            f"CD8scape exited with a non-zero status (exit code {exit_code}). "
            "Scroll the log above for details."
        )

    render_status()
    render_logs()
