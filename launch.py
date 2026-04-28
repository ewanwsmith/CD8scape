"""
launch.py — Start the CD8scape desktop app.

Usage:
    python launch.py

What it does:
    1. Checks that Python 3.8+ is being used.
    2. Installs PyQt6 if it isn't already present.
    3. Sets QT_QPA_PLATFORM_PLUGIN_PATH so the native platform plugin (cocoa
       on macOS, xcb on Linux) is always found — needed when PyQt6 is
       installed via pip rather than the system package manager.
    4. Launches the CD8scape native desktop application.

No arguments needed. Run this from anywhere inside the CD8scape folder.
"""

import os
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
UI_SCRIPT = REPO_ROOT / "ui" / "app_qt.py"

# ---------------------------------------------------------------------------
# Python version check
# ---------------------------------------------------------------------------

if sys.version_info < (3, 8):
    sys.exit(
        f"CD8scape requires Python 3.8 or newer.\n"
        f"You are running Python {sys.version}.\n"
        f"Download a newer version from https://www.python.org/downloads/"
    )

# ---------------------------------------------------------------------------
# PyQt6: install if missing
# ---------------------------------------------------------------------------

try:
    import PyQt6  # noqa: F401
except ImportError:
    print("PyQt6 not found — installing now (this only happens once)...")
    subprocess.check_call(
        [sys.executable, "-m", "pip", "install", "PyQt6>=6.4"],
    )
    print("PyQt6 installed.\n")

try:
    import setproctitle  # noqa: F401
except ImportError:
    try:
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install", "setproctitle"],
        )
    except Exception:
        pass  # Not critical — the app will still run without it

# ---------------------------------------------------------------------------
# Fix platform-plugin path (macOS pip install quirk)
#
# When PyQt6 is installed via pip the Qt6 plugins directory is inside the
# Python package rather than in a system-wide Qt installation.  Qt itself
# looks for platform plugins in directories listed in the environment variable
# QT_QPA_PLATFORM_PLUGIN_PATH; if the variable is empty Qt prints the "could
# not find platform plugin 'cocoa'" error and aborts.
#
# We locate the plugins/platforms directory that pip put inside the PyQt6
# package and inject it into the environment before we spawn the child process.
# ---------------------------------------------------------------------------

def _find_qt_plugins() -> str | None:
    """Return the path to PyQt6's Qt6/plugins directory, or None."""
    try:
        import importlib.util
        spec = importlib.util.find_spec("PyQt6")
        if spec is None or not spec.submodule_search_locations:
            return None
        pyqt6_dir = Path(list(spec.submodule_search_locations)[0])
        plugins = pyqt6_dir / "Qt6" / "plugins"
        if plugins.is_dir():
            return str(plugins)
        # Older PyQt6 wheel layout
        plugins_alt = pyqt6_dir / "Qt" / "plugins"
        if plugins_alt.is_dir():
            return str(plugins_alt)
    except Exception:
        pass
    return None


if "QT_QPA_PLATFORM_PLUGIN_PATH" not in os.environ:
    qt_plugins = _find_qt_plugins()
    if qt_plugins:
        # QT_QPA_PLATFORM_PLUGIN_PATH should point to the *platforms* sub-dir.
        platforms_dir = os.path.join(qt_plugins, "platforms")
        if os.path.isdir(platforms_dir):
            os.environ["QT_QPA_PLATFORM_PLUGIN_PATH"] = platforms_dir
        # QT_PLUGIN_PATH is a broader search path that Qt also respects.
        os.environ.setdefault("QT_PLUGIN_PATH", qt_plugins)

# ---------------------------------------------------------------------------
# Launch
# ---------------------------------------------------------------------------

print("Starting CD8scape…")

env = os.environ.copy()

try:
    subprocess.run(
        [sys.executable, str(UI_SCRIPT)],
        cwd=str(REPO_ROOT),
        env=env,
        check=False,
    )
except KeyboardInterrupt:
    print("\nStopped.")
