# CD8scape UI

Python source for the CD8scape desktop app.

## Files

| File | Purpose |
|------|---------|
| `app_qt.py` | Main PyQt6 application — the five-step wizard |
| `runner.py` | Locates Julia, builds the subprocess command, streams stdout/stderr |
| `workflow.py` | Defines the available prep and run workflows and their CLI flags |
| `setup.py` | Reads and writes `src/settings.txt`; validates the netMHCpan path |
| `paths.py` | Folder resolution and listing helpers for the data folder picker |
| `progress_estimator.py` | Estimates run time from log output for the progress bar |
| `requirements.txt` | Python dependencies (`PyQt6`) |

## Design

The UI is a thin wrapper around the CD8scape CLI. It assembles `julia CD8scape.jl <command> [flags…]` argument lists from widget state and streams the subprocess output live. It does not duplicate any pipeline logic — file discovery, variant parsing, and all error messages come from CD8scape itself.

## Running directly

From the repo root:

```bash
python launch.py
```

Or invoke the smoke tests:

```bash
python ui/runner.py                        # runs CD8scape.jl --help
python ui/runner.py read data/Example_data
python ui/paths.py data/Example_data       # folder sanity check only
```
