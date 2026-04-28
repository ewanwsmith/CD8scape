#!/bin/bash
# Launch CD8scape.command — macOS double-click launcher.
#
# Double-click this file in Finder to start CD8scape.
# A Terminal window will open while the app is running; close it to quit.

cd "$(dirname "$0")"

# ── Find Python 3.8+ ──────────────────────────────────────────────────────────
# Check conda prefix first (works when conda is active in the parent shell or
# when the user's default shell activates conda automatically).  Then fall
# through to common fixed locations so it works even from Finder.

find_python() {
    local candidate ok

    # 1. Currently active conda environment
    if [ -n "${CONDA_PREFIX:-}" ] && [ -x "$CONDA_PREFIX/bin/python3" ]; then
        echo "$CONDA_PREFIX/bin/python3"; return
    fi

    # 2. Common conda / Python install locations
    for candidate in \
        /opt/anaconda3/bin/python3 \
        /opt/miniconda3/bin/python3 \
        "$HOME/anaconda3/bin/python3" \
        "$HOME/miniconda3/bin/python3" \
        /opt/homebrew/bin/python3 \
        /usr/local/bin/python3 \
        python3
    do
        local resolved
        resolved=$(command -v "$candidate" 2>/dev/null || echo "")
        [ -z "$resolved" ] && continue
        ok=$("$resolved" -c "import sys; print(sys.version_info >= (3,8))" 2>/dev/null)
        [ "$ok" = "True" ] && { echo "$resolved"; return; }
    done
}

PYTHON=$(find_python)

if [ -z "$PYTHON" ]; then
    osascript -e 'display alert "Python not found" message "CD8scape requires Python 3.8 or newer.\n\nInstall it from https://www.python.org/downloads/ and try again."' 2>/dev/null \
        || echo "ERROR: Python 3.8+ not found. Install from https://www.python.org/downloads/"
    exit 1
fi

echo "Python: $PYTHON ($("$PYTHON" --version 2>&1))"
echo ""
echo "Starting CD8scape…"
echo ""

# ── Run ───────────────────────────────────────────────────────────────────────
"$PYTHON" launch.py
EXIT_CODE=$?

if [ $EXIT_CODE -ne 0 ]; then
    echo ""
    echo "────────────────────────────────────────────"
    echo "  CD8scape exited with an error (code $EXIT_CODE)."
    echo "  Scroll up to read the error message."
    echo "────────────────────────────────────────────"
    echo "Press Enter to close this window."
    read -r
fi
