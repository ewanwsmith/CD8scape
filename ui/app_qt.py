#!/usr/bin/env python3
"""
app_qt.py — PyQt6 native desktop application for CD8scape.

Five-step wizard:
  0  Setup    — configure netMHCpan path, check Perl, install Julia deps
  1  Data     — choose a working folder (example datasets or your own)
  2  Prepare  — choose input type and advanced options
  3  Run      — choose analysis type and advanced options
  4  Execute  — live log stream, then view / download / delete output files

Run with:
    python ui/app_qt.py          (from the CD8scape repo root)

Or via the launcher:
    python launch.py
"""
from __future__ import annotations

import os
import platform
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Fix Qt platform-plugin path BEFORE importing PyQt6.
#
# When PyQt6 is installed via pip (rather than via a system Qt), the platform
# plugins (libqcocoa.dylib on macOS, libqxcb.so on Linux) live inside the
# Python package directory.  Qt won't find them unless QT_QPA_PLATFORM_PLUGIN_PATH
# is set.  We do this here so `python ui/app_qt.py` also works without the
# launcher.
# ---------------------------------------------------------------------------

def _ensure_qt_plugin_path() -> None:
    if "QT_QPA_PLATFORM_PLUGIN_PATH" in os.environ:
        return
    try:
        import importlib.util
        spec = importlib.util.find_spec("PyQt6")
        if spec is None or not spec.submodule_search_locations:
            return
        pyqt6_dir = Path(list(spec.submodule_search_locations)[0])
        for candidate in ("Qt6", "Qt"):
            plugins = pyqt6_dir / candidate / "plugins"
            if plugins.is_dir():
                platforms = plugins / "platforms"
                if platforms.is_dir():
                    os.environ["QT_QPA_PLATFORM_PLUGIN_PATH"] = str(platforms)
                os.environ.setdefault("QT_PLUGIN_PATH", str(plugins))
                return
    except Exception:
        pass

_ensure_qt_plugin_path()

# ---------------------------------------------------------------------------
# Rename the process in the OS process list (Activity Monitor, ps, htop).
# When launched via CD8scape.app, the menu bar / Dock name comes from the
# bundle's Info.plist — no extra work needed there.  setproctitle is
# optional; it only affects the BSD process title column.
# ---------------------------------------------------------------------------

try:
    import setproctitle as _spt  # type: ignore
    _spt.setproctitle("CD8scape")
except ImportError:
    try:
        sys.argv[0] = "CD8scape"
    except Exception:
        pass

# ---------------------------------------------------------------------------
# Make sure the ui/ directory is importable regardless of cwd
# ---------------------------------------------------------------------------
_UI_DIR = Path(__file__).resolve().parent
_REPO_ROOT = _UI_DIR.parent
if str(_UI_DIR) not in sys.path:
    sys.path.insert(0, str(_UI_DIR))

from PyQt6.QtCore import Qt, QThread, QTimer, QUrl, pyqtSignal
from PyQt6.QtGui import QColor, QDesktopServices, QFont, QTextCursor
from PyQt6.QtWidgets import (
    QApplication,
    QButtonGroup,
    QCheckBox,
    QDoubleSpinBox,
    QFileDialog,
    QFrame,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMainWindow,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QRadioButton,
    QScrollArea,
    QSizePolicy,
    QSpinBox,
    QStackedWidget,
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from paths import list_folder, resolve_folder
from progress_estimator import ProgressEstimator, infer_threads_from_args
from runner import (
    CD8scapeNotFoundError,
    JuliaNotFoundError,
    REPO_ROOT,
    stream_cd8scape,
)
from setup import check_env, read_netmhcpan_path, validate_netmhcpan, write_netmhcpan_path
from workflow import PREP_CHOICES, RUN_CHOICES, WorkflowChoice, choice_by_key

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

STEP_NAMES = ["Setup", "Data", "Prepare", "Run", "Execute", "Output"]

EXAMPLE_DATASETS: Dict[str, Path] = {
    "Synthetic Example": _REPO_ROOT / "data" / "Example_data",
    "Stanevich et al.": _REPO_ROOT / "data" / "Stanevich_et_al",
}

FILE_REQUIREMENTS = {
    "nucleotide": (
        "Required files for nucleotide input:\n\n"
        "  Variant file (one of):\n"
        "    • A .vcf or .vcf.gz file\n"
        "    • A single_locus_trajectories…out file (Samfire)\n\n"
        "  Reading-frame files (one of these pairs):\n"
        "    • sequences.fasta + Consensus.fa  (NCBI-style data)\n"
        "    • Reading_Frames.dat              (Samfire-style data)\n\n"
        "  Allele files (one or both):\n"
        "    • alleles.txt         (one HLA allele per line, e.g. HLA-A03:01)\n"
        "    • supertype_panel.csv (columns: Allele, Frequency)"
    ),
    "amino_acid": (
        "Required files for amino-acid input:\n\n"
        "  Variant file:\n"
        "    • A .aa file (pre-called amino-acid variants)\n\n"
        "  Reading-frame files (one of these pairs):\n"
        "    • sequences.fasta + Consensus.fa  (NCBI-style data)\n"
        "    • Reading_Frames.dat              (Samfire-style data)\n\n"
        "  Allele files (one or both):\n"
        "    • alleles.txt         (one HLA allele per line, e.g. HLA-A03:01)\n"
        "    • supertype_panel.csv (columns: Allele, Frequency)"
    ),
    "simulation": (
        "Required files for simulation:\n\n"
        "  Reading-frame files only (no variant file needed):\n"
        "    • sequences.fasta + Consensus.fa  (NCBI-style data)\n"
        "    • Reading_Frames.dat              (Samfire-style data)\n\n"
        "  Allele files (one or both):\n"
        "    • alleles.txt         (one HLA allele per line, e.g. HLA-A03:01)\n"
        "    • supertype_panel.csv (columns: Allele, Frequency)"
    ),
}

# ---------------------------------------------------------------------------
# Platform font — pick a font that actually exists on each OS so Qt doesn't
# spend 60 ms populating alias tables for missing families.
# ---------------------------------------------------------------------------

_PLATFORM = platform.system()
if _PLATFORM == "Darwin":
    _UI_FONT = '"Helvetica Neue", "Arial"'
elif _PLATFORM == "Windows":
    _UI_FONT = '"Segoe UI", "Arial"'
else:  # Linux / FreeBSD / etc.
    _UI_FONT = '"Ubuntu", "DejaVu Sans", "Arial"'

# ---------------------------------------------------------------------------
# Stylesheet
# ---------------------------------------------------------------------------

STYLE = """
/* ── overall background ── */
QMainWindow { background-color: #f5f5f7; }
QWidget { font-family: _UI_FONT_PLACEHOLDER_; }

/* ── step bar ── */
QWidget#step_bar { background-color: #ffffff; }

/* ── page content ── */
QWidget#page_root { background-color: #f5f5f7; }
QScrollArea { background-color: #f5f5f7; border: none; }
QScrollArea > QWidget > QWidget { background-color: #f5f5f7; }

/* ── nav bar ── */
QWidget#nav_bar { background-color: #ffffff; }

/* ── cards / group boxes ── */
QGroupBox {
    border: 1px solid #d1d1d6;
    border-radius: 10px;
    margin-top: 10px;
    padding: 12px 12px 8px 12px;
    background-color: #ffffff;
    font-size: 13px;
    font-weight: 600;
    color: #1d1d1f;
}
QGroupBox::title {
    subcontrol-origin: margin;
    subcontrol-position: top left;
    padding: 0 6px;
    left: 10px;
}

/* ── inputs ── */
QLineEdit {
    border: 1px solid #d1d1d6;
    border-radius: 7px;
    padding: 6px 10px;
    background-color: #ffffff;
    font-size: 13px;
    color: #1d1d1f;
}
QLineEdit:focus { border-color: #0071e3; }

QTextEdit {
    border: 1px solid #d1d1d6;
    border-radius: 7px;
    background-color: #1e1e1e;
    color: #d4d4d4;
    font-family: "Menlo", "Monaco", "Consolas", monospace;
    font-size: 12px;
}

QListWidget {
    border: 1px solid #d1d1d6;
    border-radius: 7px;
    background-color: #ffffff;
    font-size: 13px;
    outline: none;
}
QListWidget::item { padding: 5px 10px; }
QListWidget::item:selected {
    background-color: #e3f0fd;
    color: #0071e3;
}

QSpinBox, QDoubleSpinBox, QComboBox {
    border: 1px solid #d1d1d6;
    border-radius: 7px;
    padding: 5px 8px;
    background-color: #ffffff;
    font-size: 13px;
}

QRadioButton { font-size: 13px; color: #1d1d1f; spacing: 8px; }
QCheckBox    { font-size: 13px; color: #1d1d1f; spacing: 8px; }

/* ── separator ── */
QFrame#sep { background-color: #d1d1d6; max-height: 1px; min-height: 1px; }

/* ── buttons ── */
QPushButton {
    border-radius: 8px;
    padding: 7px 18px;
    font-size: 13px;
}
QPushButton#btn_primary {
    background-color: #0071e3;
    color: #ffffff;
    border: none;
    font-weight: 600;
}
QPushButton#btn_primary:hover   { background-color: #0077ed; }
QPushButton#btn_primary:pressed { background-color: #005bb5; }
QPushButton#btn_primary:disabled { background-color: #a8a8ae; color: #ffffff; }

QPushButton#btn_secondary {
    background-color: #e5e5ea;
    color: #1d1d1f;
    border: none;
}
QPushButton#btn_secondary:hover   { background-color: #d8d8dc; }
QPushButton#btn_secondary:disabled { color: #a8a8ae; }

QPushButton#btn_outline {
    background-color: #ffffff;
    color: #0071e3;
    border: 1.5px solid #0071e3;
}
QPushButton#btn_outline:hover   { background-color: #f0f7ff; }
QPushButton#btn_outline:disabled { color: #a8a8ae; border-color: #a8a8ae; }

QPushButton#btn_danger {
    background-color: #ffffff;
    color: #ff3b30;
    border: 1.5px solid #ff3b30;
}
QPushButton#btn_danger:hover   { background-color: #fff5f5; }

QPushButton#btn_link {
    background: transparent;
    border: none;
    color: #0071e3;
    font-size: 12px;
    padding: 2px 4px;
    text-decoration: underline;
}
QPushButton#btn_link:hover { color: #0055b3; }

/* ── status labels ── */
QLabel#lbl_ok   { color: #34c759; font-size: 12px; }
QLabel#lbl_err  { color: #ff3b30; font-size: 12px; }
QLabel#lbl_warn { color: #ff9f0a; font-size: 12px; }
QLabel#lbl_info { color: #636366; font-size: 12px; }

/* ── info/warn callout boxes ── */
QLabel#callout_info {
    background-color: #e3f0fd;
    color: #004080;
    border: 1px solid #b0d0f0;
    border-radius: 7px;
    padding: 8px 12px;
    font-size: 12px;
}
QLabel#callout_warn {
    background-color: #fff8e1;
    color: #664d00;
    border: 1px solid #ffd060;
    border-radius: 7px;
    padding: 8px 12px;
    font-size: 12px;
}

/* ── tabs ── */
QTabWidget::pane {
    border: none;
    background: transparent;
}
QTabBar {
    background: transparent;
}
QTabBar::tab {
    background: #e5e5ea;
    color: #636366;
    border: none;
    padding: 7px 20px;
    font-size: 13px;
    border-top-left-radius: 6px;
    border-top-right-radius: 6px;
    margin-right: 2px;
}
QTabBar::tab:selected {
    background: #ffffff;
    color: #1d1d1f;
    font-weight: 600;
}
QTabBar::tab:hover:!selected { background: #d8d8dc; }

/* ── progress bar ── */
QProgressBar {
    border: none;
    border-radius: 5px;
    background-color: #e5e5ea;
    max-height: 10px;
    min-height: 10px;
}
QProgressBar::chunk {
    background-color: #0071e3;
    border-radius: 5px;
}

/* ── page / section titles ── */
QLabel#title_page    { font-size: 22px; font-weight: 700; color: #1d1d1f; }
QLabel#title_section { font-size: 15px; font-weight: 600; color: #1d1d1f; }
""".replace("_UI_FONT_PLACEHOLDER_", _UI_FONT)


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

def _sep() -> QFrame:
    f = QFrame()
    f.setObjectName("sep")
    f.setFrameShape(QFrame.Shape.HLine)
    return f


def _btn(text: str, style: str = "btn_primary") -> QPushButton:
    b = QPushButton(text)
    b.setObjectName(style)
    return b


def _label(text: str, obj: str = "") -> QLabel:
    lbl = QLabel(text)
    if obj:
        lbl.setObjectName(obj)
    return lbl


def _scrolled(widget: QWidget) -> QScrollArea:
    """Wrap a page widget in a borderless scroll area.

    * setWidgetResizable(True) — lets the inner widget expand to fill the
      viewport when there is room, and allows scrolling when there isn't.
    * Horizontal scrollbar is hidden; the app has a 940 px minimum width so
      horizontal overflow should never happen in practice.
    * No frame so the scroll area is invisible — the page looks the same as
      before, just gains a vertical scrollbar when the window is too short.
    """
    sa = QScrollArea()
    sa.setWidget(widget)
    sa.setWidgetResizable(True)
    sa.setFrameShape(QFrame.Shape.NoFrame)
    sa.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
    sa.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
    return sa


# ---------------------------------------------------------------------------
# Background threads
# ---------------------------------------------------------------------------

class _StreamThread(QThread):
    """Run one or two CD8scape argv tails sequentially, emitting log lines."""

    line = pyqtSignal(str)            # one log line
    phase_done = pyqtSignal(int, int) # (phase_index, returncode)
    finished_ok = pyqtSignal(bool)    # True = all phases succeeded

    def __init__(self, steps: List[List[str]], parent: Optional[QWidget] = None):
        super().__init__(parent)
        self._steps = steps

    def run(self) -> None:
        for i, args in enumerate(self._steps):
            label = " ".join(args[:3]) + ("…" if len(args) > 3 else "")
            self.line.emit(f"\n{'─'*60}")
            self.line.emit(f"Step {i+1}: julia CD8scape.jl {label}")
            self.line.emit(f"{'─'*60}\n")
            rc = 0
            try:
                for token in stream_cd8scape(args):
                    if token.startswith("__CD8SCAPE_EXIT__:"):
                        rc = int(token.split(":", 1)[1])
                    else:
                        self.line.emit(token)
            except (JuliaNotFoundError, CD8scapeNotFoundError) as exc:
                self.line.emit(f"\nERROR: {exc}\n")
                self.phase_done.emit(i, 1)
                self.finished_ok.emit(False)
                return
            self.phase_done.emit(i, rc)
            if rc != 0:
                self.line.emit(f"\nStep {i+1} failed (exit code {rc}).")
                self.finished_ok.emit(False)
                return
        self.finished_ok.emit(True)


# ---------------------------------------------------------------------------
# Step bar
# ---------------------------------------------------------------------------

class StepBar(QWidget):
    """Step pills at the top of the window. Each pill is clickable."""

    step_clicked = pyqtSignal(int)

    def __init__(self, names: List[str], parent: Optional[QWidget] = None):
        super().__init__(parent)
        self.setObjectName("step_bar")
        self._pills: List[QPushButton] = []
        self._locked: bool = False
        self._active: int = 0
        layout = QHBoxLayout(self)
        layout.setContentsMargins(24, 10, 24, 10)
        layout.setSpacing(4)

        for i, name in enumerate(names):
            if i > 0:
                arrow = QLabel("›")
                arrow.setStyleSheet("color: #a8a8ae; font-size: 16px;")
                layout.addWidget(arrow)
            pill = QPushButton(f"  {name}  ")
            pill.setFixedHeight(28)
            pill.setCursor(Qt.CursorShape.PointingHandCursor)
            pill.setFlat(True)
            pill.clicked.connect(lambda checked, idx=i: self._on_pill_clicked(idx))
            self._pills.append(pill)
            layout.addWidget(pill)
        layout.addStretch()
        self.set_step(0)

    def _on_pill_clicked(self, idx: int) -> None:
        if self._locked:
            QMessageBox.information(
                self,
                "Analysis running",
                "Navigation is locked while an analysis is running.\n\n"
                "Wait for the current run to finish before switching steps.",
            )
            return
        self.step_clicked.emit(idx)

    def set_locked(self, locked: bool) -> None:
        """Lock or unlock the step bar. Locked pills show a forbidden cursor
        and display a message when clicked instead of navigating."""
        self._locked = locked
        cursor = Qt.CursorShape.ForbiddenCursor if locked else Qt.CursorShape.PointingHandCursor
        for pill in self._pills:
            pill.setCursor(cursor)
        # Re-apply styling so locked pills get their visual treatment
        self.set_step(self._active)

    def set_step(self, active: int) -> None:
        self._active = active
        for i, pill in enumerate(self._pills):
            if self._locked:
                # All pills get a uniform dimmed/striped look while locked
                is_active = (i == active)
                pill.setStyleSheet(
                    f"QPushButton {{ background:{'#4a90d9' if is_active else '#c8c8ce'};"
                    f" color:{'#cce0f5' if is_active else '#999999'}; border-radius:12px;"
                    " padding:2px 10px; font-size:12px; border:none;"
                    f"{'font-weight:bold;' if is_active else ''} }}"
                )
            elif i < active:
                pill.setStyleSheet(
                    "QPushButton { background:#d9eaf9; color:#0071e3; border-radius:12px;"
                    " padding:2px 10px; font-size:12px; border:none; }"
                    "QPushButton:hover { background:#c5dff5; }"
                )
            elif i == active:
                pill.setStyleSheet(
                    "QPushButton { background:#0071e3; color:#ffffff; border-radius:12px;"
                    " padding:2px 10px; font-size:12px; font-weight:bold; border:none; }"
                )
            else:
                pill.setStyleSheet(
                    "QPushButton { background:#e5e5ea; color:#888888; border-radius:12px;"
                    " padding:2px 10px; font-size:12px; border:none; }"
                    "QPushButton:hover { background:#d8d8dc; color:#555555; }"
                )


# ---------------------------------------------------------------------------
# Shared option sub-widgets
# ---------------------------------------------------------------------------

class _SuffixLatestWidget(QWidget):
    """--suffix and --latest / --no-latest controls, shared across commands."""

    def __init__(self, default_suffix: str = "", parent: Optional[QWidget] = None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)

        # Suffix row
        row = QHBoxLayout()
        row.addWidget(QLabel("Output suffix (--suffix):"))
        self._suffix = QLineEdit()
        self._suffix.setPlaceholderText("leave blank for default")
        self._suffix.setText(default_suffix)
        self._suffix.setMaximumWidth(200)
        row.addWidget(self._suffix)
        row.addStretch()
        layout.addLayout(row)

        # Latest/no-latest
        latest_lbl = QLabel("When multiple input candidates exist:")
        layout.addWidget(latest_lbl)
        self._latest = QRadioButton("Pick the most recent (--latest, CD8scape default)")
        self._nlatest = QRadioButton("Error on ambiguity (--no-latest)")
        self._latest.setChecked(True)
        layout.addWidget(self._latest)
        layout.addWidget(self._nlatest)

    def suffix(self) -> str:
        return self._suffix.text().strip()

    def args(self) -> List[str]:
        result: List[str] = []
        s = self.suffix()
        if s:
            result += ["--suffix", s]
        if self._nlatest.isChecked():
            result.append("--no-latest")
        return result


class ReadOptionsWidget(QWidget):
    def __init__(self, parent: Optional[QWidget] = None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)
        self._sl = _SuffixLatestWidget(parent=self)
        layout.addWidget(self._sl)

    def args(self) -> List[str]:
        return self._sl.args()


class SimulateOptionsWidget(QWidget):
    def __init__(self, parent: Optional[QWidget] = None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)

        # Sampling mode
        layout.addWidget(QLabel("Sampling mode:"))
        self._grp = QButtonGroup(self)
        self._rb_all = QRadioButton("Write all variants (no limit)")
        self._rb_n   = QRadioButton("Sample a fixed number (--n)")
        self._rb_p   = QRadioButton("Sample a proportion (--p)")
        self._rb_n.setChecked(True)
        for rb in (self._rb_all, self._rb_n, self._rb_p):
            self._grp.addButton(rb)
            layout.addWidget(rb)

        # n / p inputs in a stacked-ish way — just show/hide
        self._n_spin = QSpinBox()
        self._n_spin.setRange(1, 10_000_000)
        self._n_spin.setValue(1000)
        self._n_spin.setMaximumWidth(160)
        n_row = QHBoxLayout()
        n_row.addWidget(QLabel("Number of variants (--n):"))
        n_row.addWidget(self._n_spin)
        n_row.addStretch()
        self._n_widget = QWidget()
        self._n_widget.setLayout(n_row)
        layout.addWidget(self._n_widget)

        self._p_spin = QDoubleSpinBox()
        self._p_spin.setRange(0.001, 0.999)
        self._p_spin.setSingleStep(0.01)
        self._p_spin.setValue(0.10)
        self._p_spin.setDecimals(3)
        self._p_spin.setMaximumWidth(160)
        p_row = QHBoxLayout()
        p_row.addWidget(QLabel("Proportion to sample (--p):"))
        p_row.addWidget(self._p_spin)
        p_row.addStretch()
        self._p_widget = QWidget()
        self._p_widget.setLayout(p_row)
        layout.addWidget(self._p_widget)
        self._p_widget.hide()

        self._rb_all.toggled.connect(self._update_visibility)
        self._rb_n.toggled.connect(self._update_visibility)
        self._rb_p.toggled.connect(self._update_visibility)

        # Seed
        seed_row = QHBoxLayout()
        seed_row.addWidget(QLabel("Random seed (--seed):"))
        self._seed = QSpinBox()
        self._seed.setRange(0, 99_999_999)
        self._seed.setValue(1320)
        self._seed.setMaximumWidth(160)
        seed_row.addWidget(self._seed)
        seed_row.addStretch()
        layout.addLayout(seed_row)

        layout.addWidget(_sep())
        self._sl = _SuffixLatestWidget(default_suffix="simulated", parent=self)
        layout.addWidget(self._sl)

    def _update_visibility(self) -> None:
        self._n_widget.setVisible(self._rb_n.isChecked())
        self._p_widget.setVisible(self._rb_p.isChecked())

    def args(self) -> List[str]:
        result: List[str] = []
        if self._rb_n.isChecked():
            result += ["--n", str(self._n_spin.value())]
        elif self._rb_p.isChecked():
            result += ["--p", f"{self._p_spin.value():g}"]
        result += ["--seed", str(self._seed.value())]
        result += self._sl.args()
        return result

    def sim_relative_size(self) -> float:
        """Estimated simulation size relative to a --n 1000 baseline.

        Used by ProgressEstimator to weight the run-simulated step fairly.
        The run-simulated step scales roughly linearly with variant count, so a
        --n 100 run should occupy ~10% of the bar that a --n 1000 run would.

        For --p (proportion), we assume a typical dataset produces ~5 000 total
        variant positions, giving a nominal count of p × 5000.  This is a
        deliberate heuristic — the actual count is data-dependent — but it puts
        the estimate in the right order of magnitude.
        """
        if self._rb_n.isChecked():
            return self._n_spin.value() / 1000.0
        elif self._rb_p.isChecked():
            # p × 5000 variants / 1000 baseline = p × 5
            return self._p_spin.value() * 5.0
        else:
            # "all variants" — no limit; treat like --p 1.0
            return 5.0


class RunOptionsWidget(QWidget):
    def __init__(self, show_supertype_note: bool = False, parent: Optional[QWidget] = None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)

        if show_supertype_note:
            note = QLabel(
                "Note: supertype-panel alleles are population-frequency surrogates, "
                "not an individual's genotype."
            )
            note.setWordWrap(True)
            note.setObjectName("lbl_info")
            layout.addWidget(note)

        # Threads
        self._tmax = QCheckBox("Use maximum parallel threads (--t max)")
        self._tmax.setChecked(False)
        layout.addWidget(self._tmax)

        t_row = QHBoxLayout()
        t_row.addWidget(QLabel("Parallel chunks (--t):"))
        self._t = QSpinBox()
        self._t.setRange(1, 256)
        self._t.setValue(1)
        self._t.setMaximumWidth(100)
        t_row.addWidget(self._t)
        t_row.addStretch()
        self._t_widget = QWidget()
        self._t_widget.setLayout(t_row)
        layout.addWidget(self._t_widget)
        self._tmax.toggled.connect(lambda v: self._t_widget.setVisible(not v))

        # per-allele / verbose
        self._per_allele = QCheckBox(
            "Per-allele output (--per-allele): writes per_allele_best_ranks.csv"
        )
        self._verbose = QCheckBox(
            "Verbose mode (--verbose): preserve temp files and per-allele logs"
        )
        layout.addWidget(self._per_allele)
        layout.addWidget(self._verbose)

        layout.addWidget(_sep())
        self._sl = _SuffixLatestWidget(parent=self)
        layout.addWidget(self._sl)

    def suffix(self) -> str:
        return self._sl.suffix()

    def per_allele(self) -> bool:
        return self._per_allele.isChecked()

    def args(self) -> List[str]:
        result: List[str] = []
        if self._tmax.isChecked():
            result += ["--t", "max"]
        else:
            result += ["--t", str(self._t.value())]
        if self._per_allele.isChecked():
            result.append("--per-allele")
        if self._verbose.isChecked():
            result.append("--verbose")
        result += self._sl.args()
        return result

    def effective_threads(self) -> int:
        """Return the numeric thread count the user has chosen.

        'max' is resolved to the platform CPU count (capped at
        CD8SCAPE_MAX_THREADS, default 8), mirroring Julia's logic.
        """
        if self._tmax.isChecked():
            import os
            cap = int(os.environ.get("CD8SCAPE_MAX_THREADS", "8"))
            return min(os.cpu_count() or 1, cap)
        return max(1, self._t.value())


# ---------------------------------------------------------------------------
# CLI argument helpers
# ---------------------------------------------------------------------------

def _inject_suffix(args: List[str], suffix: str) -> List[str]:
    """Return *args* with any existing --suffix value replaced by *suffix*.

    If --suffix is absent, it is appended.  This lets us reuse the user's
    run options for the 'run simulated' step while forcing the suffix to
    match whatever the simulate step produced.
    """
    result: List[str] = []
    replaced, i = False, 0
    while i < len(args):
        if args[i] == "--suffix" and i + 1 < len(args):
            result += ["--suffix", suffix]
            i += 2
            replaced = True
        else:
            result.append(args[i])
            i += 1
    if not replaced:
        result += ["--suffix", suffix]
    return result


# ---------------------------------------------------------------------------
# Helper: scrollable page wrapper
# ---------------------------------------------------------------------------

def _scrolled(inner: QWidget) -> QScrollArea:
    scroll = QScrollArea()
    scroll.setWidgetResizable(True)
    scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
    scroll.setWidget(inner)
    return scroll


# ---------------------------------------------------------------------------
# Page 0 — Setup
# ---------------------------------------------------------------------------

class SetupPage(QWidget):
    def __init__(self, app: "CD8scapeApp"):
        super().__init__()
        self._app = app
        self._install_thread: Optional[_StreamThread] = None

        inner = QWidget()
        inner.setObjectName("page_root")
        vbox = QVBoxLayout(inner)
        vbox.setContentsMargins(40, 32, 40, 32)
        vbox.setSpacing(20)

        # Title
        vbox.addWidget(_label("Setup", "title_page"))
        sub = QLabel("Configure the tools CD8scape needs before running an analysis.")
        sub.setObjectName("lbl_info")
        vbox.addWidget(sub)
        vbox.addWidget(_sep())

        # ── netMHCpan ──────────────────────────────────────────────────────
        nmhc_box = QGroupBox("netMHCpan path")
        nmhc_layout = QVBoxLayout(nmhc_box)
        nmhc_layout.setSpacing(8)

        note = QLabel(
            "netMHCpan must be installed separately (download from "
            "services.healthtech.dtu.dk). Paste the full path to its executable below."
        )
        note.setWordWrap(True)
        note.setObjectName("lbl_info")
        nmhc_layout.addWidget(note)

        path_row = QHBoxLayout()
        self._path_edit = QLineEdit()
        self._path_edit.setPlaceholderText("/path/to/netMHCpan-4.x/netMHCpan")
        saved = read_netmhcpan_path()
        if saved:
            self._path_edit.setText(saved)
        browse_btn = _btn("Browse…", "btn_secondary")
        browse_btn.setFixedWidth(90)
        browse_btn.clicked.connect(self._browse)
        path_row.addWidget(self._path_edit, 1)
        path_row.addWidget(browse_btn)
        nmhc_layout.addLayout(path_row)

        action_row = QHBoxLayout()
        val_btn  = _btn("Validate", "btn_outline")
        save_btn = _btn("Save", "btn_primary")
        val_btn.setFixedWidth(100)
        save_btn.setFixedWidth(100)
        val_btn.clicked.connect(self._validate)
        save_btn.clicked.connect(self._save)
        self._nmhc_status = QLabel("")
        action_row.addWidget(val_btn)
        action_row.addWidget(save_btn)
        action_row.addSpacing(12)
        action_row.addWidget(self._nmhc_status)
        action_row.addStretch()
        nmhc_layout.addLayout(action_row)
        vbox.addWidget(nmhc_box)

        # ── Perl ──────────────────────────────────────────────────────────
        perl_box = QGroupBox("Perl")
        perl_layout = QHBoxLayout(perl_box)
        self._perl_status = QLabel("")
        refresh_btn = _btn("Re-check", "btn_secondary")
        refresh_btn.setFixedWidth(90)
        refresh_btn.clicked.connect(self._check_perl)
        perl_layout.addWidget(self._perl_status, 1)
        perl_layout.addWidget(refresh_btn)
        vbox.addWidget(perl_box)
        self._check_perl()

        # ── Julia dependencies ────────────────────────────────────────────
        julia_box = QGroupBox("Julia dependencies")
        julia_layout = QVBoxLayout(julia_box)
        julia_note = QLabel(
            "Click 'Install / verify' to run CD8scape's setup step. "
            "This downloads Julia packages and verifies the netMHCpan path. "
            "Required once per machine."
        )
        julia_note.setWordWrap(True)
        julia_note.setObjectName("lbl_info")
        julia_layout.addWidget(julia_note)

        install_row = QHBoxLayout()
        self._install_btn = _btn("Install / verify dependencies", "btn_outline")
        self._install_btn.setFixedWidth(260)
        self._install_btn.clicked.connect(self._run_prep)

        self._install_status = QLabel("")
        self._install_status.setObjectName("lbl_info")

        install_row.addWidget(self._install_btn)
        install_row.addSpacing(14)
        install_row.addWidget(self._install_status)
        install_row.addStretch()
        julia_layout.addLayout(install_row)

        # Progress bar (indeterminate while running, hidden otherwise)
        self._install_bar = QProgressBar()
        self._install_bar.setRange(0, 100)
        self._install_bar.setValue(0)
        self._install_bar.setTextVisible(False)
        self._install_bar.setFixedHeight(10)
        self._install_bar.hide()
        julia_layout.addWidget(self._install_bar)

        # Terminal output — hidden by default, revealed by toggle
        self._term_toggle = _btn("▶  Show terminal output", "btn_link")
        self._term_toggle.setCheckable(True)
        self._term_toggle.setChecked(False)
        self._term_toggle.hide()
        julia_layout.addWidget(self._term_toggle)

        self._install_log = QTextEdit()
        self._install_log.setReadOnly(True)
        self._install_log.setFixedHeight(180)
        self._install_log.hide()
        julia_layout.addWidget(self._install_log)

        self._term_toggle.toggled.connect(lambda v: (
            self._install_log.setVisible(v),
            self._term_toggle.setText(
                ("▼" if v else "▶") + "  Show terminal output"
            ),
        ))

        vbox.addWidget(julia_box)

        vbox.addStretch()

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(_scrolled(inner))

    # ── actions ──────────────────────────────────────────────────────────

    def _browse(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Select netMHCpan executable", str(Path.home())
        )
        if path:
            self._path_edit.setText(path)

    def _validate(self) -> None:
        raw = self._path_edit.text().strip()
        result = validate_netmhcpan(raw)
        if result.ok:
            self._nmhc_status.setText(f"✓ {result.message}")
            self._nmhc_status.setObjectName("lbl_ok")
        else:
            self._nmhc_status.setText(f"✗ {result.message}")
            self._nmhc_status.setObjectName("lbl_err")
        self._nmhc_status.style().unpolish(self._nmhc_status)
        self._nmhc_status.style().polish(self._nmhc_status)

    def _save(self) -> None:
        raw = self._path_edit.text().strip()
        if not raw:
            self._nmhc_status.setText("✗ No path entered.")
            self._nmhc_status.setObjectName("lbl_err")
            return
        write_netmhcpan_path(raw)
        result = validate_netmhcpan(raw)
        if result.ok:
            self._nmhc_status.setText("✓ Saved and validated.")
            self._nmhc_status.setObjectName("lbl_ok")
        else:
            self._nmhc_status.setText(f"Saved (but: {result.message})")
            self._nmhc_status.setObjectName("lbl_warn")
        self._nmhc_status.style().unpolish(self._nmhc_status)
        self._nmhc_status.style().polish(self._nmhc_status)

    def _check_perl(self) -> None:
        perl = shutil.which("perl")
        if perl:
            self._perl_status.setText(f"✓ Perl found: {perl}")
            self._perl_status.setObjectName("lbl_ok")
        else:
            self._perl_status.setText(
                "✗ Perl not found on PATH. netMHCpan requires Perl — "
                "install it and reopen this app."
            )
            self._perl_status.setObjectName("lbl_err")
        self._perl_status.style().unpolish(self._perl_status)
        self._perl_status.style().polish(self._perl_status)

    def _run_prep(self) -> None:
        if self._install_thread and self._install_thread.isRunning():
            return
        self._install_log.clear()
        self._install_btn.setEnabled(False)
        self._install_btn.setText("Running…")

        # Show indeterminate progress bar; hide any previous result
        self._install_bar.setRange(0, 0)
        self._install_bar.show()
        self._install_status.setText("Running setup…")
        self._install_status.setObjectName("lbl_info")
        self._install_status.style().unpolish(self._install_status)
        self._install_status.style().polish(self._install_status)

        self._install_thread = _StreamThread([["prep"]], self)
        self._install_thread.line.connect(self._append_log)
        self._install_thread.finished_ok.connect(self._prep_done)
        self._install_thread.start()

    def _append_log(self, text: str) -> None:
        self._install_log.append(text)
        self._install_log.moveCursor(QTextCursor.MoveOperation.End)

    def _prep_done(self, ok: bool) -> None:
        self._install_btn.setEnabled(True)
        self._install_btn.setText("Install / verify dependencies")

        # Switch progress bar to determinate and show result
        self._install_bar.setRange(0, 100)
        if ok:
            self._install_bar.setValue(100)
            self._install_status.setText("✓ Setup complete")
            self._install_status.setObjectName("lbl_ok")
            self._install_log.append("\n✓ Setup complete.")
        else:
            self._install_bar.setValue(0)
            self._install_status.setText("✗ Setup failed — check terminal output")
            self._install_status.setObjectName("lbl_err")
            self._install_log.append("\n✗ Setup failed — see output below.")
            # Auto-expand terminal so the error is visible
            self._term_toggle.setChecked(True)
        self._install_status.style().unpolish(self._install_status)
        self._install_status.style().polish(self._install_status)

        # Reveal the "Show terminal output" toggle
        self._term_toggle.show()


# ---------------------------------------------------------------------------
# Page 1 — Data
# ---------------------------------------------------------------------------

class DataPage(QWidget):
    def __init__(self, app: "CD8scapeApp"):
        super().__init__()
        self._app = app

        inner = QWidget()
        inner.setObjectName("page_root")
        vbox = QVBoxLayout(inner)
        vbox.setContentsMargins(40, 32, 40, 32)
        vbox.setSpacing(20)

        vbox.addWidget(_label("Choose your data folder", "title_page"))
        sub = QLabel(
            "CD8scape reads all input files from a single folder. "
            "Select an example dataset or point to your own data."
        )
        sub.setWordWrap(True)
        sub.setObjectName("lbl_info")
        vbox.addWidget(sub)
        vbox.addWidget(_sep())

        # ── Example datasets ─────────────────────────────────────────────
        ex_box = QGroupBox("Example datasets")
        ex_layout = QVBoxLayout(ex_box)
        ex_layout.setSpacing(8)

        ex_note = QLabel(
            "Use one of the bundled datasets to try out CD8scape without your own data."
        )
        ex_note.setWordWrap(True)
        ex_layout.addWidget(ex_note)

        ex_btn_row = QHBoxLayout()
        for name, folder in EXAMPLE_DATASETS.items():
            btn = _btn(name, "btn_outline")
            btn.setFixedWidth(200)
            btn.clicked.connect(lambda checked, f=folder, n=name: self._use_example(f, n))
            if not folder.is_dir():
                btn.setEnabled(False)
                btn.setToolTip(f"Not found: {folder}")
            ex_btn_row.addWidget(btn)
        ex_btn_row.addStretch()
        ex_layout.addLayout(ex_btn_row)

        self._ex_status = QLabel("")
        self._ex_status.setObjectName("lbl_ok")
        ex_layout.addWidget(self._ex_status)
        vbox.addWidget(ex_box)

        # ── Custom folder ─────────────────────────────────────────────────
        own_box = QGroupBox("Or use your own folder")
        own_layout = QVBoxLayout(own_box)
        own_layout.setSpacing(8)

        folder_row = QHBoxLayout()
        self._folder_edit = QLineEdit()
        self._folder_edit.setPlaceholderText("/path/to/your/data/folder")
        self._folder_edit.textChanged.connect(self._on_path_changed)
        browse_btn = _btn("Browse…", "btn_secondary")
        browse_btn.setFixedWidth(90)
        browse_btn.clicked.connect(self._browse)
        folder_row.addWidget(self._folder_edit, 1)
        folder_row.addWidget(browse_btn)
        own_layout.addLayout(folder_row)

        self._folder_status = QLabel("")
        own_layout.addWidget(self._folder_status)

        # File listing
        self._file_list = QListWidget()
        self._file_list.setFixedHeight(180)
        self._file_list.setSelectionMode(QListWidget.SelectionMode.NoSelection)
        own_layout.addWidget(self._file_list)
        self._file_count_lbl = QLabel("")
        self._file_count_lbl.setObjectName("lbl_info")
        own_layout.addWidget(self._file_count_lbl)

        vbox.addWidget(own_box)
        vbox.addStretch()

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(_scrolled(inner))

    def _use_example(self, folder: Path, name: str) -> None:
        self._app.folder_path = folder
        self._folder_edit.setText(str(folder))
        self._ex_status.setText(f"✓ Using: {name}")
        self._refresh_listing(folder)

    def _browse(self) -> None:
        start = str(self._app.folder_path or Path.home())
        folder = QFileDialog.getExistingDirectory(self, "Select data folder", start)
        if folder:
            self._folder_edit.setText(folder)

    def _on_path_changed(self, text: str) -> None:
        check = resolve_folder(text)
        if check.ok:
            self._app.folder_path = check.path
            self._folder_status.setText(f"✓ {check.path}")
            self._folder_status.setObjectName("lbl_ok")
            self._refresh_listing(check.path)
        else:
            self._app.folder_path = None
            self._folder_status.setText(check.message or "")
            self._folder_status.setObjectName("lbl_err" if text.strip() else "lbl_info")
            self._file_list.clear()
            self._file_count_lbl.clear()
        self._folder_status.style().unpolish(self._folder_status)
        self._folder_status.style().polish(self._folder_status)

    def _refresh_listing(self, folder: Path) -> None:
        self._file_list.clear()
        entries, total = list_folder(folder, limit=50)
        for name in entries:
            self._file_list.addItem(name)
        extra = total - len(entries)
        if extra > 0:
            self._file_count_lbl.setText(f"{total} items  (+{extra} more not shown)")
        else:
            self._file_count_lbl.setText(f"{total} item{'s' if total != 1 else ''}")

    def refresh(self) -> None:
        """Called when the page becomes visible to restore state."""
        if self._app.folder_path:
            self._folder_edit.setText(str(self._app.folder_path))


# ---------------------------------------------------------------------------
# Page 2 — Prepare
# ---------------------------------------------------------------------------

class PreparePage(QWidget):
    def __init__(self, app: "CD8scapeApp"):
        super().__init__()
        self._app = app

        inner = QWidget()
        inner.setObjectName("page_root")
        vbox = QVBoxLayout(inner)
        vbox.setContentsMargins(40, 32, 40, 32)
        vbox.setSpacing(20)

        vbox.addWidget(_label("Prepare data", "title_page"))
        sub = QLabel("What kind of input data do you have?")
        sub.setObjectName("lbl_info")
        vbox.addWidget(sub)
        vbox.addWidget(_sep())

        # ── Choice radio buttons ──────────────────────────────────────────
        choice_box = QGroupBox("Input type")
        choice_layout = QVBoxLayout(choice_box)
        choice_layout.setSpacing(10)

        self._btn_group = QButtonGroup(self)
        self._radios: List[QRadioButton] = []
        for c in PREP_CHOICES:
            rb = QRadioButton(c.label)
            rb.setChecked(c.key == self._app.prep_key)
            self._btn_group.addButton(rb, PREP_CHOICES.index(c))
            self._radios.append(rb)
            choice_layout.addWidget(rb)
        vbox.addWidget(choice_box)

        # Detail / file requirements info box
        self._detail_lbl = QLabel("")
        self._detail_lbl.setObjectName("callout_info")
        self._detail_lbl.setWordWrap(True)
        self._detail_lbl.setTextFormat(Qt.TextFormat.PlainText)
        vbox.addWidget(self._detail_lbl)

        # ── Advanced options (collapsible) ────────────────────────────────
        adv_toggle = _btn("▶  Advanced options", "btn_link")
        adv_toggle.setCheckable(True)
        adv_toggle.setChecked(False)
        vbox.addWidget(adv_toggle)

        self._adv_box = QGroupBox("Advanced options")
        self._adv_box.hide()
        adv_inner = QVBoxLayout(self._adv_box)
        adv_inner.setSpacing(8)

        self._read_opts     = ReadOptionsWidget(self)
        self._simulate_opts = SimulateOptionsWidget(self)
        adv_inner.addWidget(self._read_opts)
        adv_inner.addWidget(self._simulate_opts)
        vbox.addWidget(self._adv_box)

        adv_toggle.toggled.connect(lambda v: (
            self._adv_box.setVisible(v),
            adv_toggle.setText(("▼" if v else "▶") + "  Advanced options"),
        ))

        vbox.addStretch()

        # Connect signals
        self._btn_group.buttonClicked.connect(self._on_choice_changed)
        self._update_detail()

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(_scrolled(inner))

    def _on_choice_changed(self) -> None:
        idx = self._btn_group.checkedId()
        if idx >= 0:
            self._app.prep_key = PREP_CHOICES[idx].key
        self._update_detail()

    def _update_detail(self) -> None:
        key = self._app.prep_key
        # Show file requirements for the selected choice
        req = FILE_REQUIREMENTS.get(key, "")
        self._detail_lbl.setText(req)

        # Toggle which options widget is visible
        is_sim = key == "simulation"
        self._simulate_opts.setVisible(is_sim)
        self._read_opts.setVisible(not is_sim)

    def current_option_args(self) -> List[str]:
        key = self._app.prep_key
        if key == "simulation":
            return self._simulate_opts.args()
        return self._read_opts.args()


# ---------------------------------------------------------------------------
# Page 3 — Run
# ---------------------------------------------------------------------------

class RunPage(QWidget):
    def __init__(self, app: "CD8scapeApp"):
        super().__init__()
        self._app = app

        inner = QWidget()
        inner.setObjectName("page_root")
        vbox = QVBoxLayout(inner)
        vbox.setContentsMargins(40, 32, 40, 32)
        vbox.setSpacing(20)

        vbox.addWidget(_label("Run analysis", "title_page"))
        sub = QLabel("What kind of HLA analysis do you want to run?")
        sub.setObjectName("lbl_info")
        vbox.addWidget(sub)
        vbox.addWidget(_sep())

        # ── Choice radio buttons ──────────────────────────────────────────
        choice_box = QGroupBox("Analysis type")
        choice_layout = QVBoxLayout(choice_box)
        choice_layout.setSpacing(10)

        self._btn_group = QButtonGroup(self)
        self._radios: List[QRadioButton] = []
        for c in RUN_CHOICES:
            rb = QRadioButton(c.label)
            rb.setChecked(c.key == self._app.run_key)
            self._btn_group.addButton(rb, RUN_CHOICES.index(c))
            self._radios.append(rb)
            choice_layout.addWidget(rb)
        vbox.addWidget(choice_box)

        # Detail info box
        self._detail_lbl = QLabel("")
        self._detail_lbl.setObjectName("callout_info")
        self._detail_lbl.setWordWrap(True)
        vbox.addWidget(self._detail_lbl)

        # Alleles file warning — updates based on choice
        self._alleles_warn = QLabel("")
        self._alleles_warn.setObjectName("callout_warn")
        self._alleles_warn.setWordWrap(True)
        vbox.addWidget(self._alleles_warn)

        # ── Advanced options (collapsible) ────────────────────────────────
        adv_toggle = _btn("▶  Advanced options", "btn_link")
        adv_toggle.setCheckable(True)
        adv_toggle.setChecked(False)
        vbox.addWidget(adv_toggle)

        self._adv_box = QGroupBox("Advanced options")
        self._adv_box.hide()
        adv_inner = QVBoxLayout(self._adv_box)

        self._std_opts   = RunOptionsWidget(show_supertype_note=False, parent=self)
        self._super_opts = RunOptionsWidget(show_supertype_note=True,  parent=self)
        adv_inner.addWidget(self._std_opts)
        adv_inner.addWidget(self._super_opts)
        vbox.addWidget(self._adv_box)

        adv_toggle.toggled.connect(lambda v: (
            self._adv_box.setVisible(v),
            adv_toggle.setText(("▼" if v else "▶") + "  Advanced options"),
        ))

        # ── Percentile analysis ───────────────────────────────────────────
        vbox.addWidget(_sep())

        pct_box = QGroupBox("Percentile analysis")
        pct_outer = QVBoxLayout(pct_box)
        pct_outer.setSpacing(10)

        pct_desc = QLabel(
            "Compare observed results against a simulated background to compute "
            "enrichment percentiles.  When enabled, three extra steps run automatically: "
            "simulate → run (simulated) → percentile."
        )
        pct_desc.setWordWrap(True)
        pct_desc.setObjectName("lbl_info")
        pct_outer.addWidget(pct_desc)

        self._pct_enable = QCheckBox("Include percentile analysis in this run")
        pct_outer.addWidget(self._pct_enable)

        # Settings sub-widget (hidden until enabled)
        self._pct_settings = QWidget()
        pct_sv = QVBoxLayout(self._pct_settings)
        pct_sv.setContentsMargins(0, 4, 0, 0)
        pct_sv.setSpacing(12)

        sim_box = QGroupBox("Simulation settings")
        sim_inner = QVBoxLayout(sim_box)
        self._pct_sim_opts = SimulateOptionsWidget(self)
        sim_inner.addWidget(self._pct_sim_opts)
        pct_sv.addWidget(sim_box)

        self._pct_per_allele = QCheckBox(
            "Per-allele percentiles (--per-allele)  —  enable if you also "
            "checked --per-allele in Advanced options above"
        )
        pct_sv.addWidget(self._pct_per_allele)

        pct_outer.addWidget(self._pct_settings)
        self._pct_settings.hide()
        self._pct_enable.toggled.connect(self._pct_settings.setVisible)

        vbox.addWidget(pct_box)
        vbox.addStretch()

        self._btn_group.buttonClicked.connect(self._on_choice_changed)
        self._update_detail()

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(_scrolled(inner))

    def _on_choice_changed(self) -> None:
        idx = self._btn_group.checkedId()
        if idx >= 0:
            self._app.run_key = RUN_CHOICES[idx].key
        self._update_detail()

    def _update_detail(self) -> None:
        key = self._app.run_key
        choice = choice_by_key(RUN_CHOICES, key)
        self._detail_lbl.setText(choice.detail)

        if key == "standard":
            self._alleles_warn.setText(
                "Requires alleles.txt in your data folder.\n"
                "Format: one HLA allele per line, e.g.  HLA-A03:01"
            )
            self._std_opts.setVisible(True)
            self._super_opts.setVisible(False)
        else:
            self._alleles_warn.setText(
                "Requires supertype_panel.csv in your data folder.\n"
                "Format: two columns — Allele, Frequency"
            )
            self._std_opts.setVisible(False)
            self._super_opts.setVisible(True)

    def current_option_args(self) -> List[str]:
        key = self._app.run_key
        if key == "supertype":
            return self._super_opts.args()
        return self._std_opts.args()

    def effective_threads(self) -> int:
        """Thread count currently selected (resolves 'max' to a real number)."""
        key = self._app.run_key
        opts = self._super_opts if key == "supertype" else self._std_opts
        return opts.effective_threads()

    def run_suffix(self) -> str:
        """The --suffix value for the observed run (empty = CD8scape default)."""
        key = self._app.run_key
        opts = self._super_opts if key == "supertype" else self._std_opts
        return opts.suffix()

    def is_per_allele(self) -> bool:
        """True if --per-allele is ticked in the active run options widget."""
        key = self._app.run_key
        opts = self._super_opts if key == "supertype" else self._std_opts
        return opts.per_allele()

    def sim_relative_size(self) -> float:
        """Relative size of the percentile simulation (1.0 = --n 1000 baseline)."""
        return self._pct_sim_opts.sim_relative_size()

    # ── Percentile accessors used by ExecutePage ──────────────────────────

    @property
    def include_percentile(self) -> bool:
        return self._pct_enable.isChecked()

    def sim_option_args(self) -> List[str]:
        return self._pct_sim_opts.args()

    def sim_suffix(self) -> str:
        """The --suffix value chosen for the simulate step (default 'simulated')."""
        args = self._pct_sim_opts.args()
        try:
            return args[args.index("--suffix") + 1]
        except (ValueError, IndexError):
            return "simulated"

    def percentile_option_args(self) -> List[str]:
        return ["--per-allele"] if self._pct_per_allele.isChecked() else []


# ---------------------------------------------------------------------------
# Page 4 — Execute
# ---------------------------------------------------------------------------

class ExecutePage(QWidget):
    def __init__(self, app: "CD8scapeApp"):
        super().__init__()
        self._app = app
        self._thread: Optional[_StreamThread] = None
        self._pre_run_snapshot: Dict[str, float] = {}
        self._n_active_steps: int = 2
        self._steps_done: int = 0
        self._estimator: Optional[ProgressEstimator] = None
        # Snapshotted in prepare() — forwarded to OutputPage after run
        self._snap_run_suffix: str = ""
        self._snap_per_allele: bool = False
        self._snap_include_pct: bool = False
        self._snap_pct_per_allele: bool = False

        vbox = QVBoxLayout(self)
        vbox.setContentsMargins(0, 0, 0, 0)
        vbox.setSpacing(0)

        # ── Summary panel (top) ───────────────────────────────────────────
        summary_widget = QWidget()
        summary_widget.setObjectName("page_root")
        summary_layout = QVBoxLayout(summary_widget)
        summary_layout.setContentsMargins(40, 24, 40, 16)
        summary_layout.setSpacing(8)

        summary_layout.addWidget(_label("Execute", "title_page"))

        self._summary_lbl = QLabel("")
        self._summary_lbl.setObjectName("callout_info")
        self._summary_lbl.setWordWrap(True)
        summary_layout.addWidget(self._summary_lbl)

        run_row = QHBoxLayout()
        self._run_btn = _btn("▶  Run Analysis", "btn_primary")
        self._run_btn.setFixedHeight(40)
        self._run_btn.setFixedWidth(180)
        self._run_btn.clicked.connect(self._start_run)

        self._status_lbl = QLabel("")
        self._status_lbl.setObjectName("lbl_info")

        run_row.addWidget(self._run_btn)
        run_row.addSpacing(16)
        run_row.addWidget(self._status_lbl)
        run_row.addStretch()
        summary_layout.addLayout(run_row)

        vbox.addWidget(summary_widget)
        vbox.addWidget(_sep())

        # ── Tab widget: Progress / Terminal ───────────────────────────────
        self._tabs = QTabWidget()
        self._tabs.setDocumentMode(True)

        # Tab 1 — Progress view
        prog_widget = QWidget()
        prog_widget.setObjectName("page_root")
        prog_layout = QVBoxLayout(prog_widget)
        prog_layout.setContentsMargins(32, 24, 32, 24)
        prog_layout.setSpacing(20)

        # ── Overall progress bar (honest: completed_steps / total_steps) ──
        overall_hdr = QHBoxLayout()
        overall_title = QLabel("Overall progress")
        overall_title.setObjectName("title_section")
        self._progress_lbl = QLabel("0%")
        self._progress_lbl.setObjectName("lbl_info")
        self._progress_lbl.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        overall_hdr.addWidget(overall_title)
        overall_hdr.addWidget(self._progress_lbl)
        prog_layout.addLayout(overall_hdr)

        self._overall_bar = QProgressBar()
        self._overall_bar.setRange(0, 2)      # updated in prepare()
        self._overall_bar.setValue(0)
        self._overall_bar.setTextVisible(False)
        self._overall_bar.setFixedHeight(14)
        prog_layout.addWidget(self._overall_bar)

        self._step_label = QLabel("Ready to run")
        self._step_label.setObjectName("lbl_info")
        prog_layout.addWidget(self._step_label)

        prog_layout.addSpacing(8)
        prog_layout.addWidget(_sep())
        prog_layout.addSpacing(8)

        # ── Per-step status cards (no individual bars — icon + title + status) ──
        # Tuple: (QFrame card, icon_lbl, name_lbl, status_lbl)
        # Cards 2-4 are percentile-only, hidden until percentile is enabled.
        self._step_rows: List[tuple] = []
        step_info = [
            ("Step 1 — Prepare data",      "Reads input files → variants.csv / frames.csv"),
            ("Step 2 — Run (observed)",    "Runs netMHCpan on your observed variants"),
            ("Step 3 — Simulate variants", "Generates synthetic variants for background"),
            ("Step 4 — Run (simulated)",   "Runs netMHCpan on the simulated variants"),
            ("Step 5 — Percentile",        "Computes enrichment percentiles vs. simulated"),
        ]
        for title, subtitle in step_info:
            card = QFrame()
            card.setObjectName("step_card")
            # QFrame clips its styled background correctly; QWidget does not.
            # The #step_card selector avoids cascading into child labels.
            card.setStyleSheet(
                "QFrame#step_card { background: #ffffff; border: 1px solid #e5e5ea;"
                " border-radius: 10px; }"
            )
            card_layout = QVBoxLayout(card)
            card_layout.setContentsMargins(20, 12, 20, 12)
            card_layout.setSpacing(2)

            row = QHBoxLayout()
            icon_lbl = QLabel("○")
            icon_lbl.setFixedWidth(22)
            icon_lbl.setStyleSheet(
                "font-size:18px; color:#a8a8ae; border:none; background:transparent;"
            )
            name_lbl = QLabel(title)
            name_lbl.setObjectName("title_section")
            name_lbl.setStyleSheet("border:none; background:transparent;")
            status_lbl = QLabel("Waiting")
            status_lbl.setObjectName("lbl_info")
            status_lbl.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
            status_lbl.setStyleSheet("border:none; background:transparent;")
            row.addWidget(icon_lbl)
            row.addSpacing(6)
            row.addWidget(name_lbl, 1)
            row.addWidget(status_lbl)
            card_layout.addLayout(row)

            sub_lbl = QLabel(subtitle)
            sub_lbl.setObjectName("lbl_info")
            sub_lbl.setStyleSheet(
                "border:none; background:transparent; padding-left: 28px;"
            )
            card_layout.addWidget(sub_lbl)

            prog_layout.addWidget(card)
            self._step_rows.append((card, icon_lbl, name_lbl, status_lbl))

        # Cards 2-4 are percentile-only; hidden until percentile is enabled
        for card, *_ in self._step_rows[2:]:
            card.hide()

        prog_layout.addStretch()
        self._tabs.addTab(prog_widget, "Progress")

        # Tab 2 — Terminal output (raw log)
        self._log = QTextEdit()
        self._log.setReadOnly(True)
        self._tabs.addTab(self._log, "Terminal output")

        vbox.addWidget(self._tabs, 1)

    # ── Called when step 4 becomes visible ───────────────────────────────

    def prepare(self) -> None:
        """Populate summary + reset progress. Called just before the page is shown.

        If a run is already in progress we skip the reset entirely — navigating
        back to this step mid-run (e.g. via the step bar) must not wipe the
        live progress state.
        """
        if self._thread and self._thread.isRunning():
            return   # don't disturb an active run

        # Collect current option args from wizard pages
        self._app.prep_option_args = self._app.prepare_page.current_option_args()
        self._app.run_option_args  = self._app.run_page.current_option_args()

        prep_choice = choice_by_key(PREP_CHOICES, self._app.prep_key)
        run_choice  = choice_by_key(RUN_CHOICES,  self._app.run_key)
        run_page    = self._app.run_page
        include_pct = run_page.include_percentile

        folder_str = str(self._app.folder_path) if self._app.folder_path else "—"
        summary_lines = [
            f"Folder:   {folder_str}",
            f"Prepare:  {prep_choice.label}",
            f"Run:      {run_choice.label}",
        ]
        if include_pct:
            summary_lines.append(f"Simulate: {run_page.sim_suffix()} suffix")
            summary_lines.append("Percentile: enabled")
        self._summary_lbl.setText("\n".join(summary_lines))

        # Show / hide percentile step cards
        for card, *_ in self._step_rows[2:]:
            card.setVisible(include_pct)

        # Build an estimator that knows the pipeline shape, thread count,
        # and simulation size (so a --n 100 run gets a proportionally smaller
        # slice of the bar than a --p 1.0 run).
        n_threads = self._app.run_page.effective_threads()
        sim_factor = (
            self._app.run_page.sim_relative_size() if include_pct else 1.0
        )
        self._estimator = ProgressEstimator(
            include_percentile=include_pct,
            n_threads=n_threads,
            sim_size_factor=sim_factor,
        )

        # Snapshot run config so output highlighting is correct after the run
        # even if the user changes settings on the Run page mid-run.
        run_page = self._app.run_page
        self._snap_run_suffix    = run_page.run_suffix()
        self._snap_per_allele    = run_page.is_per_allele()
        self._snap_include_pct   = include_pct
        self._snap_pct_per_allele = "--per-allele" in run_page.percentile_option_args()

        # Reset all cards and overall bar
        self._steps_done = 0
        self._n_active_steps = 5 if include_pct else 2
        for i in range(len(self._step_rows)):
            self._set_card_state(i, "○", "#a8a8ae", "Waiting", "lbl_info")
        self._update_overall_bar()

        self._tabs.setCurrentIndex(0)

        self._run_btn.setEnabled(True)
        self._status_lbl.setText("")

    # ── Build argv tails ─────────────────────────────────────────────────

    def _build_steps(self) -> List[List[str]]:
        folder = str(self._app.folder_path)

        # Step 1 — prepare observed (read / read --aa / simulate)
        prep_choice = choice_by_key(PREP_CHOICES, self._app.prep_key)
        prep_args = (
            [prep_choice.cli_command, folder]
            + prep_choice.cli_fixed_args
            + self._app.prep_option_args
        )

        # Step 2 — run observed
        run_choice = choice_by_key(RUN_CHOICES, self._app.run_key)
        run_args = (
            [run_choice.cli_command, folder]
            + run_choice.cli_fixed_args
            + self._app.run_option_args
        )

        steps = [prep_args, run_args]

        if self._app.run_page.include_percentile:
            run_page   = self._app.run_page
            sim_sfx    = run_page.sim_suffix()       # e.g. "simulated"
            sim_opts   = run_page.sim_option_args()  # n/p/seed/suffix/latest
            pct_opts   = run_page.percentile_option_args()

            # Step 3 — simulate
            sim_args = ["simulate", folder] + sim_opts

            # Step 4 — run simulated: same command as step 2, suffix overridden
            run_sim_args = (
                [run_choice.cli_command, folder]
                + run_choice.cli_fixed_args
                + _inject_suffix(self._app.run_option_args, sim_sfx)
            )

            # Step 5 — percentile
            pct_args = ["percentile", folder] + pct_opts

            steps += [sim_args, run_sim_args, pct_args]

        return steps

    # ── Snapshot helpers ─────────────────────────────────────────────────

    def _take_snapshot(self) -> Dict[str, float]:
        folder = self._app.folder_path
        if not folder or not folder.is_dir():
            return {}
        return {
            p.name: p.stat().st_mtime
            for p in folder.iterdir()
            if p.is_file()
        }

    def _new_or_modified(self, before: Dict[str, float]) -> List[Path]:
        folder = self._app.folder_path
        if not folder or not folder.is_dir():
            return []
        after = {
            p.name: p.stat().st_mtime
            for p in folder.iterdir()
            if p.is_file()
        }
        results = []
        for name, mtime in sorted(after.items()):
            if name not in before or before[name] < mtime:
                results.append(folder / name)
        return results

    # ── Run ──────────────────────────────────────────────────────────────

    # ── Step card state helpers ───────────────────────────────────────────

    def _set_card_state(
        self, i: int, icon: str, icon_color: str,
        status_text: str, status_obj: str
    ) -> None:
        """Update a step card's icon and status label."""
        _card, lbl_icon, _, lbl_status = self._step_rows[i]
        lbl_icon.setText(icon)
        lbl_icon.setStyleSheet(
            f"font-size:18px; color:{icon_color}; border:none; background:transparent;"
        )
        lbl_status.setText(status_text)
        lbl_status.setObjectName(status_obj)
        lbl_status.style().unpolish(lbl_status)
        lbl_status.style().polish(lbl_status)

    def _update_overall_bar(self) -> None:
        """Refresh the overall progress bar and its text label.

        The bar uses the ProgressEstimator's weighted overall percentage when
        available (finer granularity than simple step-counting), and falls back
        to the coarse step fraction if the estimator hasn't been created yet.
        """
        done  = self._steps_done
        total = self._n_active_steps

        if self._estimator is not None:
            pct = self._estimator.overall_pct
            # Bar range 0–100 for smooth sub-step updates
            self._overall_bar.setRange(0, 100)
            self._overall_bar.setValue(pct)
        else:
            pct = int(done / total * 100) if total else 0
            self._overall_bar.setRange(0, total)
            self._overall_bar.setValue(done)

        self._progress_lbl.setText(f"{pct}%")
        if done == total and total > 0:
            self._step_label.setText("Complete")
        elif done > 0:
            self._step_label.setText(f"Step {done + 1} of {total}")
        else:
            self._step_label.setText(f"Step 1 of {total}")

    def _update_phase(self, phase: int, rc: int) -> None:
        """Called by _StreamThread.phase_done — update card, estimator, and bar."""
        success = rc == 0
        if self._estimator is not None:
            self._estimator.complete_step(phase, success=success)

        if success:
            self._set_card_state(phase, "✓", "#34c759", "Complete", "lbl_ok")
            self._steps_done += 1
            self._update_overall_bar()
            next_phase = phase + 1
            if next_phase < self._n_active_steps:
                self._set_card_state(next_phase, "◉", "#0071e3", "Running…", "lbl_info")
                if self._estimator is not None:
                    self._estimator.start_step(next_phase)
        else:
            self._set_card_state(phase, "✗", "#ff3b30", "Failed", "lbl_err")
            self._update_overall_bar()

    def _start_run(self) -> None:
        if not self._app.folder_path:
            QMessageBox.warning(self, "No folder", "No folder selected.")
            return
        if self._thread and self._thread.isRunning():
            return

        self._log.clear()
        self._run_btn.setEnabled(False)
        self._status_lbl.setText("Running…")
        self._app.set_nav_enabled(False)

        # Reset all cards, then mark step 0 as running
        self._steps_done = 0
        self._n_active_steps = 5 if self._app.run_page.include_percentile else 2
        for i in range(len(self._step_rows)):
            self._set_card_state(i, "○", "#a8a8ae", "Waiting", "lbl_info")
        self._set_card_state(0, "◉", "#0071e3", "Running…", "lbl_info")
        if self._estimator is not None:
            self._estimator.start_step(0)
        self._update_overall_bar()

        self._pre_run_snapshot = self._take_snapshot()

        steps = self._build_steps()
        self._thread = _StreamThread(steps, self)
        self._thread.line.connect(self._append_log)
        self._thread.phase_done.connect(self._update_phase)
        self._thread.finished_ok.connect(self._run_finished)
        self._thread.start()

    def _append_log(self, text: str) -> None:
        self._log.append(text)
        self._log.moveCursor(QTextCursor.MoveOperation.End)
        # Feed the line to the estimator; update the bar if within-step
        # progress was found (i.e. the overall_pct actually changed).
        if self._estimator is not None:
            prev = self._estimator.overall_pct
            self._estimator.ingest_line(text)
            if self._estimator.overall_pct != prev:
                self._update_overall_bar()

    def _run_finished(self, ok: bool) -> None:
        self._app.set_nav_enabled(True)
        if ok:
            self._status_lbl.setText("✓ Analysis complete")
            self._status_lbl.setObjectName("lbl_ok")
        else:
            self._status_lbl.setText("✗ Failed — see Terminal output tab")
            self._status_lbl.setObjectName("lbl_err")
            # Automatically switch to the terminal tab so the error is visible
            self._tabs.setCurrentIndex(1)
        self._status_lbl.style().unpolish(self._status_lbl)
        self._status_lbl.style().polish(self._status_lbl)

        output_files = self._new_or_modified(self._pre_run_snapshot)

        if ok:
            # Hand off to the Output page and navigate there
            self._app.output_page.populate(
                files=output_files,
                snap_run_suffix=self._snap_run_suffix,
                snap_sim_suffix=(
                    self._app.run_page.sim_suffix()
                    if self._snap_include_pct else ""
                ),
                snap_per_allele=self._snap_per_allele,
                snap_include_pct=self._snap_include_pct,
                snap_pct_per_allele=self._snap_pct_per_allele,
                folder=self._app.folder_path,
            )
            self._app._set_step(self._app._STEP_OUTPUT)

    # ── (output file actions are in OutputPage) ──────────────────────────


# ---------------------------------------------------------------------------
# Page 5 — Output
# ---------------------------------------------------------------------------

class OutputPage(QWidget):
    """Wizard step 5: browse, preview, and manage the files produced by a run."""

    def __init__(self, app: "CD8scapeApp"):
        super().__init__()
        self._app = app
        self._output_files: List[Path] = []
        self._snap_run_suffix: str = ""
        self._snap_per_allele: bool = False
        self._snap_include_pct: bool = False
        self._snap_pct_per_allele: bool = False

        vbox = QVBoxLayout(self)
        vbox.setContentsMargins(0, 0, 0, 0)
        vbox.setSpacing(0)

        # ── Header ────────────────────────────────────────────────────────
        hdr = QWidget()
        hdr.setObjectName("page_root")
        hdr_layout = QVBoxLayout(hdr)
        hdr_layout.setContentsMargins(40, 24, 40, 16)
        hdr_layout.setSpacing(4)
        hdr_layout.addWidget(_label("Output", "title_page"))

        self._subtitle_lbl = QLabel("Run an analysis to see output files here.")
        self._subtitle_lbl.setObjectName("lbl_info")
        hdr_layout.addWidget(self._subtitle_lbl)

        # Variant fates funnel (shown when variant_fates.csv is available)
        self._fates_lbl = QLabel("")
        self._fates_lbl.setObjectName("callout_info")
        self._fates_lbl.setWordWrap(True)
        self._fates_lbl.hide()
        hdr_layout.addWidget(self._fates_lbl)

        vbox.addWidget(hdr)
        vbox.addWidget(_sep())

        # ── Body: file list (left) + preview (right) ──────────────────────
        body = QWidget()
        body.setObjectName("page_root")
        body_h = QHBoxLayout(body)
        body_h.setContentsMargins(40, 20, 40, 20)
        body_h.setSpacing(20)

        # ── Left column: file list + action buttons ───────────────────────
        left = QVBoxLayout()
        left.setSpacing(10)

        files_hdr = _label("Output files", "title_section")
        left.addWidget(files_hdr)

        self._output_list = QListWidget()
        left.addWidget(self._output_list, 1)

        # Action buttons (2×2 grid feels cleaner than a single row)
        btn_grid = QVBoxLayout()
        btn_grid.setSpacing(6)
        row1 = QHBoxLayout()
        self._open_btn = _btn("Open selected", "btn_outline")
        self._open_btn.clicked.connect(self._open_selected)
        self._reveal_btn = _btn("Show in Finder", "btn_secondary")
        self._reveal_btn.clicked.connect(self._reveal_selected)
        row1.addWidget(self._open_btn)
        row1.addWidget(self._reveal_btn)
        row2 = QHBoxLayout()
        self._zip_btn = _btn("Save all as ZIP…", "btn_outline")
        self._zip_btn.clicked.connect(self._save_zip)
        self._delete_btn = _btn("Delete output files", "btn_danger")
        self._delete_btn.clicked.connect(self._delete_outputs)
        row2.addWidget(self._zip_btn)
        row2.addWidget(self._delete_btn)
        btn_grid.addLayout(row1)
        btn_grid.addLayout(row2)
        left.addLayout(btn_grid)

        # ── Right column: preview ─────────────────────────────────────────
        right = QVBoxLayout()
        right.setSpacing(6)
        right.addWidget(_label("File preview", "title_section"))

        # Stacked pane: 0 = table (CSV/TSV), 1 = plain text fallback
        self._preview_stack = QStackedWidget()

        # ── Table view (for CSV / TSV) ────────────────────────────────────
        self._preview_table = QTableWidget()
        self._preview_table.setEditTriggers(
            QTableWidget.EditTrigger.NoEditTriggers
        )
        self._preview_table.setAlternatingRowColors(True)
        self._preview_table.setSelectionBehavior(
            QTableWidget.SelectionBehavior.SelectRows
        )
        self._preview_table.verticalHeader().setVisible(False)
        self._preview_table.horizontalHeader().setStretchLastSection(True)
        self._preview_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.ResizeToContents
        )
        self._preview_table.setShowGrid(True)
        self._preview_table.setStyleSheet(
            "QTableWidget { font-size: 12px; }"
            "QTableWidget::item { padding: 3px 8px; }"
            "QHeaderView::section { font-weight: bold; padding: 4px 8px;"
            " background: #f2f2f7; border: none; border-bottom: 1px solid #d1d1d6; }"
        )
        self._preview_stack.addWidget(self._preview_table)

        # ── Plain-text fallback ───────────────────────────────────────────
        self._preview_text = QTextEdit()
        self._preview_text.setReadOnly(True)
        self._preview_text.setPlaceholderText("Select a file to preview its contents.")
        _mono = QFont(
            "Menlo" if platform.system() == "Darwin"
            else "Consolas" if platform.system() == "Windows"
            else "DejaVu Sans Mono"
        )
        _mono.setPointSize(11)
        self._preview_text.setFont(_mono)
        self._preview_stack.addWidget(self._preview_text)

        self._preview_stack.setCurrentIndex(1)   # start on text pane
        right.addWidget(self._preview_stack, 1)

        body_h.addLayout(left, 2)
        body_h.addLayout(right, 3)
        vbox.addWidget(body, 1)

        # Wire preview on selection change
        self._output_list.currentItemChanged.connect(self._update_preview)

        # Adjust platform label
        if platform.system() == "Windows":
            self._reveal_btn.setText("Show in Explorer")
        elif platform.system() == "Linux":
            self._reveal_btn.setText("Show in Files")

    # ── Public API called by ExecutePage after a run ──────────────────────

    def populate(
        self,
        *,
        files: List[Path],
        snap_run_suffix: str,
        snap_sim_suffix: str,
        snap_per_allele: bool,
        snap_include_pct: bool,
        snap_pct_per_allele: bool,
        folder: Optional[Path],
    ) -> None:
        """Called by ExecutePage._run_finished after a successful pipeline run."""
        self._output_files = files
        self._snap_run_suffix   = snap_run_suffix
        self._snap_sim_suffix   = snap_sim_suffix
        self._snap_per_allele   = snap_per_allele
        self._snap_include_pct  = snap_include_pct
        self._snap_pct_per_allele = snap_pct_per_allele
        self._folder = folder

        self._output_list.clear()
        self._clear_preview()
        self._fates_lbl.hide()
        self._subtitle_lbl.setText("Run complete — results are listed below.")

        self._populate_file_list()
        self._parse_and_show_fates()

    # ── File list population ──────────────────────────────────────────────

    def _primary_output_names(self) -> set:
        sfx = f"_{self._snap_run_suffix}" if self._snap_run_suffix else ""
        names: set = set()
        if self._snap_include_pct:
            if self._snap_pct_per_allele:
                names.add("percentile_per_allele_best_ranks.csv")
            else:
                names.add("percentile_harmonic_mean_best_ranks.csv")
            names.add(f"harmonic_mean_best_ranks{sfx}.csv")
            if self._snap_per_allele:
                names.add(f"per_allele_best_ranks{sfx}.csv")
        else:
            names.add(f"harmonic_mean_best_ranks{sfx}.csv")
            if self._snap_per_allele:
                names.add(f"per_allele_best_ranks{sfx}.csv")
        return names

    def _populate_file_list(self) -> None:
        if not self._output_files:
            return
        primary_names = self._primary_output_names()
        primary   = [p for p in self._output_files if p.name in primary_names]
        secondary = [p for p in self._output_files if p.name not in primary_names]

        self._output_list.clear()
        for p in primary:
            size_kb = p.stat().st_size // 1024
            item = QListWidgetItem(f"★  {p.name}  ({size_kb} KB)")
            item.setData(Qt.ItemDataRole.UserRole, str(p))
            item.setForeground(QColor("#0071e3"))
            font = item.font(); font.setBold(True); item.setFont(font)
            self._output_list.addItem(item)
        if primary and secondary:
            div = QListWidgetItem("── other files ──")
            div.setFlags(Qt.ItemFlag.NoItemFlags)
            div.setForeground(QColor("#a8a8ae"))
            self._output_list.addItem(div)
        for p in secondary:
            size_kb = p.stat().st_size // 1024
            item = QListWidgetItem(f"  {p.name}  ({size_kb} KB)")
            item.setData(Qt.ItemDataRole.UserRole, str(p))
            self._output_list.addItem(item)

        # Auto-select the first primary file so the preview loads immediately
        if self._output_list.count() > 0:
            self._output_list.setCurrentRow(0)

    # ── Variant fates ─────────────────────────────────────────────────────

    @staticmethod
    def _count_fates(path: Path) -> Optional[Tuple[int, int, int, int]]:
        """Parse variant_fates.csv and return (total, in_frame, non_syn, passed).
        Returns None on any error or if the file is empty."""
        import csv as _csv
        try:
            total = in_frame = non_syn = passed = 0
            with open(path, newline="", encoding="utf-8", errors="replace") as fh:
                for row in _csv.DictReader(fh):
                    total += 1
                    ff = (row.get("frame_filter") or "").strip()
                    pf = (row.get("peptide_filter") or "").strip()
                    bf = (row.get("binding_filter") or "").strip()
                    if ff == "In frame":
                        in_frame += 1
                        if pf in ("Non-synonymous", "Stop codon"):
                            non_syn += 1
                            if bf == "Passed":
                                passed += 1
            return (total, in_frame, non_syn, passed) if total > 0 else None
        except Exception:
            return None

    @staticmethod
    def _fates_line(label: str, counts: Tuple[int, int, int, int]) -> str:
        total, in_frame, non_syn, passed = counts
        return (
            f"{label}:  {total} total"
            f"  →  {in_frame} in frame"
            f"  →  {non_syn} non-synonymous"
            f"  →  {passed} passed binding"
        )

    def _resolve_fates_path(self, suffix: str) -> Optional[Path]:
        folder = getattr(self, "_folder", None)
        if not folder or not folder.is_dir():
            return None
        candidates = []
        if suffix:
            candidates.append(folder / f"variant_fates_{suffix}.csv")
        candidates.append(folder / "variant_fates.csv")
        return next((c for c in candidates if c.exists()), None)

    def _parse_and_show_fates(self) -> None:
        lines: List[str] = []

        # Observed variants
        obs_path = self._resolve_fates_path(self._snap_run_suffix)
        if obs_path is not None:
            counts = self._count_fates(obs_path)
            if counts:
                lines.append(self._fates_line("Observed variants", counts))

        # Simulated variants (only when percentile mode was active)
        if self._snap_include_pct:
            sim_sfx = getattr(self, "_snap_sim_suffix", "simulated")
            sim_path = self._resolve_fates_path(sim_sfx)
            if sim_path is not None:
                counts = self._count_fates(sim_path)
                if counts:
                    lines.append(self._fates_line("Simulated variants", counts))

        if lines:
            self._fates_lbl.setText("\n".join(lines))
            self._fates_lbl.show()

    # ── Preview ───────────────────────────────────────────────────────────

    def _clear_preview(self) -> None:
        self._preview_table.clearContents()
        self._preview_table.setRowCount(0)
        self._preview_table.setColumnCount(0)
        self._preview_text.clear()
        self._preview_stack.setCurrentIndex(1)

    def _show_table(self, path: Path) -> None:
        """Render a CSV or TSV file as a QTableWidget."""
        import csv as _csv
        delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
        MAX_ROWS = 500
        try:
            with open(path, newline="", encoding="utf-8", errors="replace") as fh:
                reader = _csv.reader(fh, delimiter=delimiter)
                rows: List[List[str]] = []
                for i, row in enumerate(reader):
                    if i > MAX_ROWS:
                        rows.append([f"… ({i - 1} rows shown, file truncated)"])
                        break
                    rows.append(row)
        except OSError as exc:
            self._preview_text.setPlainText(f"Could not read file:\n{exc}")
            self._preview_stack.setCurrentIndex(1)
            return

        if not rows:
            self._preview_text.setPlainText("(empty file)")
            self._preview_stack.setCurrentIndex(1)
            return

        headers = rows[0]
        data    = rows[1:]

        self._preview_table.setRowCount(len(data))
        self._preview_table.setColumnCount(len(headers))
        self._preview_table.setHorizontalHeaderLabels(headers)

        for r, row in enumerate(data):
            for c, val in enumerate(row):
                item = QTableWidgetItem(val)
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self._preview_table.setItem(r, c, item)

        # Resize columns to content but cap each at 200 px so wide tables don't overflow
        self._preview_table.resizeColumnsToContents()
        for col in range(self._preview_table.columnCount()):
            if self._preview_table.columnWidth(col) > 200:
                self._preview_table.setColumnWidth(col, 200)

        self._preview_stack.setCurrentIndex(0)

    def _update_preview(self, current, _previous) -> None:
        if current is None:
            self._clear_preview()
            return
        path_str = current.data(Qt.ItemDataRole.UserRole)
        if not path_str:
            self._clear_preview()
            return
        p = Path(path_str)
        if not p.exists() or not p.is_file():
            self._preview_text.setPlainText(f"File not found:\n{p}")
            self._preview_stack.setCurrentIndex(1)
            return

        suffix = p.suffix.lower()
        if suffix in {".csv", ".tsv"}:
            self._show_table(p)
        elif suffix in {".pep", ".txt", ".out", ".fasta", ".fa", ".log", ".dat"}:
            try:
                with open(p, encoding="utf-8", errors="replace") as fh:
                    lines: List[str] = []
                    for i, line in enumerate(fh):
                        if i >= 300:
                            lines.append("… (preview truncated at 300 lines)")
                            break
                        lines.append(line.rstrip("\r\n"))
                self._preview_text.setPlainText("\n".join(lines))
                self._preview_stack.setCurrentIndex(1)
            except OSError as exc:
                self._preview_text.setPlainText(f"Could not read file:\n{exc}")
                self._preview_stack.setCurrentIndex(1)
        else:
            self._preview_text.setPlainText(
                f"Binary / unsupported type ({p.suffix}) — open externally to view."
            )
            self._preview_stack.setCurrentIndex(1)

    # ── Actions ───────────────────────────────────────────────────────────

    def _selected_path(self) -> Optional[Path]:
        item = self._output_list.currentItem()
        if item:
            v = item.data(Qt.ItemDataRole.UserRole)
            return Path(v) if v else None
        return None

    def _open_selected(self) -> None:
        p = self._selected_path()
        if not p and self._output_files:
            p = self._output_files[0]
        if p and p.exists():
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(p)))

    def _reveal_selected(self) -> None:
        p = self._selected_path()
        if not p and self._output_files:
            p = self._output_files[0]
        if not p:
            return
        system = platform.system()
        if system == "Darwin":
            subprocess.Popen(["open", "-R", str(p)])
        elif system == "Windows":
            subprocess.Popen(["explorer", "/select,", str(p)])
        else:
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(p.parent)))

    def _save_zip(self) -> None:
        if not self._output_files:
            return
        folder = getattr(self, "_folder", None)
        default_name = (folder.name if folder else "cd8scape") + "_output.zip"
        dest, _ = QFileDialog.getSaveFileName(
            self, "Save output files as ZIP",
            str(Path.home() / default_name), "ZIP archives (*.zip)",
        )
        if not dest:
            return
        try:
            with zipfile.ZipFile(dest, "w", zipfile.ZIP_DEFLATED) as zf:
                for p in self._output_files:
                    if p.exists():
                        zf.write(p, p.name)
            QMessageBox.information(
                self, "Saved", f"Saved {len(self._output_files)} file(s) to:\n{dest}"
            )
        except OSError as exc:
            QMessageBox.critical(self, "Error", f"Could not write ZIP:\n{exc}")

    def _delete_outputs(self) -> None:
        if not self._output_files:
            return
        names = "\n".join(f"  • {p.name}" for p in self._output_files if p.exists())
        reply = QMessageBox.question(
            self, "Delete output files?",
            f"This will permanently delete:\n\n{names}\n\nContinue?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
        )
        if reply != QMessageBox.StandardButton.Yes:
            return
        deleted, failed = 0, 0
        for p in self._output_files:
            if p.exists():
                try:
                    p.unlink(); deleted += 1
                except OSError:
                    failed += 1
        self._output_list.clear()
        self._output_files.clear()
        self._fates_lbl.hide()
        self._clear_preview()
        self._subtitle_lbl.setText("Run an analysis to see output files here.")
        msg = f"Deleted {deleted} file(s)."
        if failed:
            msg += f"\n{failed} file(s) could not be deleted."
        QMessageBox.information(self, "Done", msg)


# ---------------------------------------------------------------------------
# Main window
# ---------------------------------------------------------------------------

class CD8scapeApp(QMainWindow):
    # Step indices
    _STEP_SETUP    = 0
    _STEP_DATA     = 1
    _STEP_PREPARE  = 2
    _STEP_RUN      = 3
    _STEP_EXECUTE  = 4
    _STEP_OUTPUT   = 5

    def __init__(self):
        super().__init__()
        self.setWindowTitle("CD8scape")
        self.setMinimumSize(940, 680)
        self.resize(1040, 760)

        # ── Shared state ──────────────────────────────────────────────────
        self.folder_path: Optional[Path] = None
        self.prep_key: str = PREP_CHOICES[0].key
        self.run_key: str  = RUN_CHOICES[0].key
        self.prep_option_args: List[str] = []
        self.run_option_args:  List[str] = []

        # ── Central widget ────────────────────────────────────────────────
        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        # Step bar
        self.step_bar = StepBar(STEP_NAMES)
        root.addWidget(self.step_bar)
        root.addWidget(_sep())

        # Pages
        self.stack = QStackedWidget()
        self.setup_page   = SetupPage(self)
        self.data_page    = DataPage(self)
        self.prepare_page = PreparePage(self)
        self.run_page     = RunPage(self)
        self.execute_page = ExecutePage(self)
        self.output_page  = OutputPage(self)

        for page in (
            self.setup_page,
            self.data_page,
            self.prepare_page,
            self.run_page,
            self.execute_page,
            self.output_page,
        ):
            self.stack.addWidget(_scrolled(page))
        root.addWidget(self.stack, 1)

        # Nav bar
        root.addWidget(_sep())
        nav = QWidget()
        nav.setObjectName("nav_bar")
        nav_layout = QHBoxLayout(nav)
        nav_layout.setContentsMargins(24, 10, 24, 10)

        self.back_btn = _btn("← Back", "btn_secondary")
        self.back_btn.setFixedWidth(110)
        self.back_btn.clicked.connect(self._go_back)

        self.next_btn = _btn("Next →", "btn_primary")
        self.next_btn.setFixedWidth(160)
        self.next_btn.clicked.connect(self._go_next)

        nav_layout.addWidget(self.back_btn)
        nav_layout.addStretch()
        nav_layout.addWidget(self.next_btn)
        root.addWidget(nav)

        self.step_bar.step_clicked.connect(self._set_step)

        self._current_step = 0
        self._set_step(0)

    # ── Navigation ────────────────────────────────────────────────────────

    def _set_step(self, step: int) -> None:
        self._current_step = step
        self.stack.setCurrentIndex(step)
        self.step_bar.set_step(step)

        self.back_btn.setVisible(step > 0)

        if step == self._STEP_EXECUTE:
            self.next_btn.hide()
            self.execute_page.prepare()
        elif step == self._STEP_OUTPUT:
            # Output is the terminal step; no Next button needed
            self.next_btn.hide()
        else:
            self.next_btn.show()
            if step == self._STEP_RUN:
                self.next_btn.setText("Run Analysis →")
            else:
                self.next_btn.setText("Next →")

        if step == self._STEP_DATA:
            self.data_page.refresh()

    def _go_back(self) -> None:
        if self._current_step > 0:
            self._set_step(self._current_step - 1)

    def _go_next(self) -> None:
        step = self._current_step

        if step == self._STEP_DATA and self.folder_path is None:
            QMessageBox.warning(
                self,
                "No folder selected",
                "Please select a working folder before continuing.\n\n"
                "You can use one of the example datasets or choose your own folder.",
            )
            return

        if step < self._STEP_EXECUTE:
            self._set_step(step + 1)
        # Output step is reached automatically after a run; Next is hidden there

    def set_nav_enabled(self, enabled: bool) -> None:
        """Lock / unlock all navigation (back, next, and step bar) during a run."""
        self.back_btn.setEnabled(enabled)
        self.next_btn.setEnabled(enabled)
        self.step_bar.set_locked(not enabled)

    def closeEvent(self, event) -> None:  # type: ignore[override]
        """Wait for any background threads before letting Qt destroy widgets.

        Without this, Qt emits "QThread: Destroyed while thread is still
        running" when the window is closed mid-run, because the QThread
        objects (owned by child widgets) are torn down before the underlying
        OS thread has exited.
        """
        threads = [
            self.setup_page._install_thread,
            self.execute_page._thread,
        ]
        for t in threads:
            if t is not None and t.isRunning():
                # Give the thread up to 5 s to finish naturally.
                # If it's still going (e.g. a very long netMHCpan run), we
                # accept the close anyway — the process will be reaped by the
                # OS once Python exits.
                t.wait(5000)
        event.accept()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    app = QApplication(sys.argv)
    app.setApplicationName("CD8scape")
    app.setOrganizationName("CD8scape")
    app.setStyleSheet(STYLE)

    # Use Fusion style for consistent cross-platform look
    app.setStyle("Fusion")

    window = CD8scapeApp()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
