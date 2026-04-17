# CD8scape UI — Phases 1 · 2 · 3

A minimal, cross-platform graphical launcher for CD8scape.

- **Phase 1 (shipped):** cross-platform subprocess launch + live log stream.
- **Phase 2 (shipped):** pick your data folder; the UI checks only that
  the path is a folder and forwards everything to CD8scape.
- **Phase 3 (shipped):** full CLI surface (all six subcommands, every flag)
  exposed as UI controls that are 1:1 mirrors of CD8scape's `--help`.
- **Phase 4+ (planned):** results viewer, polish.

---

## Design principle — the UI is a wrapper, nothing more

The UI renders CLI arguments as UI controls, spawns `julia CD8scape.jl …`,
and streams the output. It deliberately does **not** know anything about
CD8scape's input files, output files, or pipeline logic. If CD8scape's
internals change, this UI keeps working unchanged.

Consequences, in practice:

- The folder picker checks only that the path exists and is a directory.
  Discovery of variant files, reading-frame files, etc. is entirely up to
  CD8scape.
- Every toggle, textbox and dropdown in the UI is a 1:1 mirror of a
  CD8scape CLI flag.
- Errors shown to the user come from CD8scape itself (the log panel),
  not from a parallel validator in the UI.

## What the UI currently does

- A four-section Streamlit page:
  1. **Choose a CD8scape command.** Dropdown listing every subcommand in
     `./CD8scape.jl --help`: `prep`, `read`, `simulate`, `run`,
     `run_supertype`, `percentile`. A short description is shown below it.
  2. **Choose your data folder.** Text input (drag-and-drop from Finder /
     File Explorer works). A **Check folder** button does a basic
     filesystem sanity check and — on success — shows a read-only listing
     of the folder's contents. No interpretation. This section is hidden
     for `prep`, which takes no folder.
  3. **Options.** Controls that map directly onto the CLI flags of the
     selected subcommand. The set of controls changes with the command.
  4. **Run CD8scape.** Invokes `julia CD8scape.jl <command> [<folder>]
     [flags…]` through a cross-platform subprocess. Stdout and stderr
     stream into a live log panel, followed by a green / red banner and
     the exit code.
- The sidebar shows OS, Python version, the repo root, and the exact
  `argv` list that will be spawned.

## CLI ↔ UI mapping

| Command          | UI controls (all are direct mirrors of CLI flags) |
|------------------|---------------------------------------------------|
| `prep`           | — (no options) |
| `read`           | `--aa` checkbox · `--suffix` text · `--latest`/`--no-latest` radio |
| `simulate`       | sampling radio (all / `--n` / `--p`) · `--seed` number · `--suffix` (default *simulated*) · `--latest`/`--no-latest` radio |
| `run`            | `--t` number or `--t max` checkbox · `--per-allele` · `--verbose` · `--suffix` · `--latest`/`--no-latest` |
| `run_supertype`  | same as `run` |
| `percentile`     | `--per-allele` · `--s` text (simulated file) · `--o` text (observed file) |

The UI does not interpret what these flags do; it simply assembles them
into an argv list and hands the whole thing to CD8scape. The "Next
command" block in the sidebar always shows the full argv that will be
executed, so you can see the translation for yourself before clicking
Run.

## What the UI deliberately does NOT do

- No duplication of CD8scape's file-discovery logic.
- No auto-setting of CLI flags based on folder contents — that's the
  user's call.
- No results viewer — Phase 4 will render the CSVs inline.
- No run cancellation. Phase 5 will add a Cancel button and progress bar.
- No config save/load.

---

## Requirements

| Tool | Version | Notes |
|------|---------|-------|
| Python | 3.9+ | For Streamlit and the runner module. |
| Julia | 1.11+ | Required by CD8scape itself. Install from <https://julialang.org/downloads/>. |
| Streamlit | 1.32+ | Installed via `requirements.txt`. |

For Phase 2's `read` stage you need Julia and the CD8scape Julia packages
installed (run `./CD8scape.jl prep` once from a terminal before using the UI).
You do **not** need `netMHCpan` or Perl until Phase 3's run stages.

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
python ui/runner.py                          # runs `CD8scape.jl --help`
python ui/runner.py read data/Example_data --aa
python ui/paths.py data/Example_data         # UI-only folder sanity check
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

## Folder check — what it actually does

The UI's **Check folder** button only answers filesystem questions:

- Did you enter something?
- Does it exist?
- Is it a directory?

If so, it shows a plain listing of what's in the folder as a convenience.
It does **not** classify files or decide what CD8scape will do with them.
That's CD8scape's job, and if CD8scape can't use the folder it will say
so in the run log.

## Known limitations

- The Run button currently invokes the `read` stage only. The other
  CD8scape commands (`simulate`, `run`, `run_supertype`, `percentile`)
  arrive in Phase 3.
- No CLI parameters (thread count, suffix, `--per-allele`, etc.) yet.
- Log lines are rendered into a single code block; very long runs
  (thousands of lines) can be slow to redraw.
- Closing the browser tab does not cancel a running subprocess yet;
  Phase 5 will add explicit cancellation.
