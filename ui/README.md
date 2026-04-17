# CD8scape UI — Phase 1

A minimal, cross-platform graphical launcher for CD8scape. Phase 1 focuses
on **one thing only**: proving we can safely run CD8scape from a UI on
macOS, Linux, and Windows, and stream its logs live.

Later phases will add file pickers, parameter controls, result viewers,
and polish. See the project plan for the phase-by-phase roadmap.

---

## What Phase 1 does

- A single-page Streamlit app with:
  - one **Run CD8scape** button
  - a live log panel (stdout + stderr, merged and streamed)
  - a **success / failure** banner when the run ends
- Hardcoded invocation: `CD8scape.jl --help`
  - Exercises the full launch path without requiring netMHCpan or
    any Julia packages to be installed.

## What Phase 1 deliberately does NOT do

- No file pickers
- No parameter controls
- No results viewer
- No configuration saving

---

## Requirements

| Tool | Version | Notes |
|------|---------|-------|
| Python | 3.9+ | For Streamlit and the runner module. |
| Julia | 1.11+ | Required by CD8scape itself. Install from <https://julialang.org/downloads/>. |
| Streamlit | 1.32+ | Installed via `requirements.txt`. |

You do **not** need `netMHCpan` or Perl installed to run the Phase 1 smoke
test — `--help` does not touch them. You will need them for later phases.

---

## Install

From the repo root (the folder that contains `CD8scape.jl`):

### macOS / Linux (bash/zsh)

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r ui/requirements.txt
```

### Windows (PowerShell)

```powershell
py -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r ui\requirements.txt
```

### Windows (cmd.exe)

```cmd
py -m venv .venv
.venv\Scripts\activate.bat
pip install -r ui\requirements.txt
```

---

## Run the UI

From the repo root, with the venv active:

```bash
streamlit run ui/app.py
```

Streamlit will open a browser tab at <http://localhost:8501>. Click
**Run CD8scape**; you should see the CD8scape help text stream in and
a green success banner.

---

## Run the headless smoke test

Useful for CI and for debugging before touching Streamlit:

```bash
python ui/runner.py            # runs `CD8scape.jl --help`
python ui/runner.py read data/Example_data --aa
```

Exit codes:

- `0` — CD8scape ran and exited successfully
- `1` — CD8scape ran but exited non-zero
- `2` — Environment error (Julia or CD8scape.jl could not be located)

---

## Cross-platform design — why this works on every OS

This is the single most important property of Phase 1, so it is worth
spelling out.

1. **No shell invocation.** `subprocess.Popen` is called with an argv
   list and `shell=False` (the default for list args). That means there
   is no bash, zsh, cmd, or PowerShell parser between us and the OS, and
   no quoting/escaping differences to worry about.
2. **Julia located by PATH, not by shebang.** We use `shutil.which("julia")`.
   On Windows this finds `julia.exe` via `PATHEXT`; on macOS/Linux it
   finds `julia` on `PATH`. Windows does not honour `#!/usr/bin/env julia`,
   so shelling out to the script directly would fail there — we always
   invoke as `julia CD8scape.jl …`.
3. **Paths via `pathlib.Path`.** Never a hardcoded `/` or `\`. The repo
   root is computed from `__file__`, so the UI works from any check-out
   location.
4. **UTF-8 with replacement.** `text=True, encoding="utf-8",
   errors="replace"` means decoding never crashes regardless of the
   user's locale (cp1252 on Windows, UTF-8 on macOS/Linux, etc.).
5. **Line-buffered merged pipe.** `stderr=subprocess.STDOUT, bufsize=1`
   preserves ordering and prevents us from waiting for a full block of
   output before rendering.
6. **Working directory is the repo root.** `cwd=REPO_ROOT` means
   CD8scape's relative `include(...)` calls resolve identically to how
   they do when a user runs it manually.

## Known limitations (by design, for Phase 1)

- Only the hardcoded `--help` command can be run from the UI.
- The sidebar only shows the resolved command — it doesn't offer to
  edit it (that's Phase 3).
- Log lines are rendered into a single code block; very long runs
  (thousands of lines) will be slow to redraw. Phase 4/5 will add a
  virtualised log viewer if needed.
- Closing the browser tab does not cancel a running subprocess yet;
  Phase 5 will add explicit cancellation.
