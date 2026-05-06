"""
plots_widget.py — Embedded Matplotlib plot viewer for CD8scape outputs.

Provides PlotViewer, a QWidget with sub-tabs that render the same family
of diagnostic plots as the reference notebook (plot_output.ipynb):

  • HMBR        — ancestral vs derived harmonic-mean best rank (Plot A)
  • log₂ FC     — log₂ fold-change in HMBR per mutation (Plot B)
  • Simulated   — per-mutation density of null distribution with observed mark
  • Per-Allele  — per-allele log₂ fold-change scatter (requires --per-allele)
  • Allele Box  — per-allele distribution boxplot (requires --per-allele)

All plots are faceted by protein / reading-frame and sorted by genomic locus.

Requires matplotlib (>=3.7).  If matplotlib is absent the widget degrades
gracefully to a plain "install matplotlib" notice.
"""
from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QFrame,
    QLabel,
    QScrollArea,
    QSizePolicy,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

# ── Optional matplotlib import ────────────────────────────────────────────────
_MPL_OK = False
try:
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.figure import Figure
    from matplotlib.gridspec import GridSpec
    from matplotlib.patches import Patch
    import matplotlib.ticker as mticker
    import numpy as np
    _MPL_OK = True
except ImportError:
    pass

# ── Protein-name → short label (matching plot_a.R / plot_sim.R) ──────────────
FRAME_ABBREVS: Dict[str, str] = {
    "RNA-dependent RNA polymerase": "RdRp",
    "surface glycoprotein":         "S",
    "envelope protein":             "E",
    "membrane glycoprotein":        "M",
    "nucleocapsid phosphoprotein":  "N",
    "ORF3a protein":                "ORF3a",
    "ORF6 protein":                 "ORF6",
    "ORF7a protein":                "ORF7a",
    "ORF7b protein":                "ORF7b",
    "ORF8 protein":                 "ORF8",
    "helicase":                     "helicase",
    "endoRNAse":                    "endoRNAse",
}

# HLA-locus colour palette (approximates R viridis mako begin=0.2, end=0.8)
LOCUS_COLORS: Dict[str, str] = {
    "A": "#2d1160",
    "B": "#2272b5",
    "C": "#6ecdc8",
}
LOCUS_FALLBACK = "#888888"


# ── Helpers ───────────────────────────────────────────────────────────────────

def _safe_float(val: str) -> Optional[float]:
    try:
        f = float(val)
        return None if (math.isnan(f) or math.isinf(f)) else f
    except (ValueError, TypeError):
        return None


def _load_csv(path: Path) -> List[Dict[str, str]]:
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        return [dict(row) for row in csv.DictReader(fh)]


def _find_file(folder: Path, stem: str, suffix: str) -> Optional[Path]:
    """Locate <stem>_<suffix>.csv or <stem>.csv inside folder."""
    if suffix:
        p = folder / f"{stem}_{suffix}.csv"
        if p.exists():
            return p
    p = folder / f"{stem}.csv"
    return p if p.exists() else None


def _viridis(n: int) -> List[tuple]:
    """Return n evenly-spaced colours from the viridis colormap."""
    if not _MPL_OK or n == 0:
        return []
    import matplotlib.cm as cm
    cmap = cm.get_cmap("viridis")
    return [cmap(i / max(n - 1, 1)) for i in range(n)]


def _frame_order_and_colors(
    rows: List[Dict],
) -> Tuple[List[str], Dict[str, tuple]]:
    """
    Return (frame_order, frame_colors) where frame_order is sorted by the
    minimum genomic locus appearing in each frame (matching the notebook).
    """
    min_locus: Dict[str, float] = {}
    for r in rows:
        frame = r.get("Frame", "")
        locus = _safe_float(r.get("Locus", ""))
        if frame and locus is not None:
            if frame not in min_locus or locus < min_locus[frame]:
                min_locus[frame] = locus
    frame_order = sorted(min_locus, key=lambda f: min_locus[f])
    colors = _viridis(len(frame_order))
    frame_colors = dict(zip(frame_order, colors))
    return frame_order, frame_colors


def _mutations_by_frame(
    rows: List[Dict],
    frame_order: List[str],
    value_cols: List[str],
) -> Dict[str, List[Dict]]:
    """Group and sort rows by frame / locus, keeping only rows with ≥1 valid value."""
    grouped: Dict[str, List[Dict]] = {f: [] for f in frame_order}
    for r in rows:
        frame = r.get("Frame", "")
        if frame not in grouped:
            continue
        locus = _safe_float(r.get("Locus", ""))
        if locus is None:
            continue
        if not any(_safe_float(r.get(c, "")) is not None for c in value_cols):
            continue
        grouped[frame].append(r)
    for frame in grouped:
        grouped[frame].sort(key=lambda r: _safe_float(r.get("Locus", "")) or 0)
    return grouped


def _apply_facet_style(ax, col_idx: int, frame: str, y_label: str) -> None:
    """Apply shared facet aesthetics to a subplot axis."""
    abbrev = FRAME_ABBREVS.get(frame, frame)
    ax.set_title(abbrev, fontsize=9, fontweight="bold", pad=3,
                 bbox=dict(facecolor="#f2f2f7", edgecolor="none", pad=2))
    ax.set_facecolor("white")
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
        spine.set_color("#aaaaaa")
    ax.tick_params(axis="both", labelsize=8, length=3)
    ax.grid(axis="y", linewidth=0.3, color="#eeeeee", zorder=0)
    if col_idx == 0:
        ax.set_ylabel(y_label, fontsize=9)
    else:
        ax.set_ylabel("")
        ax.tick_params(axis="y", labelleft=False)


def _make_faceted_figure(
    frame_order: List[str],
    grouped: Dict[str, List[Dict]],
    fig_height: float = 5.0,
    mutation_width: float = 0.55,
    min_fig_width: float = 10.0,
) -> Tuple["Figure", List, List[str]]:
    """
    Build a Figure with one subplot per frame.
    Returns (fig, axes_list, active_frames).
    """
    active = [f for f in frame_order if grouped.get(f)]
    if not active:
        raise ValueError("No data to plot.")
    counts = [len(grouped[f]) for f in active]
    total = sum(counts)
    fig_width = max(min_fig_width, total * mutation_width + len(active) * 0.8 + 1.5)

    fig = Figure(figsize=(fig_width, fig_height), dpi=100)
    fig.patch.set_facecolor("#f5f5f7")
    gs = GridSpec(
        1, len(active),
        figure=fig,
        width_ratios=counts,
        wspace=0.05,
        left=0.06, right=0.99,
        bottom=0.30, top=0.88,
    )
    axes = [fig.add_subplot(gs[0, i]) for i in range(len(active))]
    return fig, axes, active


# ── Individual plot renderers ─────────────────────────────────────────────────

def _render_hmbr(rows: List[Dict]) -> "Figure":
    """
    Plot A: dodged bar chart of HMBR_A (ancestral) vs HMBR_D (derived),
    faceted by Frame, sorted by Locus.  Y-axis is sqrt-scaled.
    """
    frame_order, frame_colors = _frame_order_and_colors(rows)
    grouped = _mutations_by_frame(rows, frame_order, ["HMBR_A", "HMBR_D"])

    fig, axes, active = _make_faceted_figure(frame_order, grouped, fig_height=5.2)

    # Collect all finite HMBR values for global y-range
    all_vals = []
    for r in rows:
        for col in ("HMBR_A", "HMBR_D"):
            v = _safe_float(r.get(col, ""))
            if v is not None:
                all_vals.append(v)
    y_max_lin = max(all_vals) * 1.08 if all_vals else 10.0
    sqrt_max = math.sqrt(y_max_lin)

    bar_w = 0.38
    threshold_sqrt = math.sqrt(2.0)

    for col_idx, (ax, frame) in enumerate(zip(axes, active)):
        mutations = grouped[frame]
        color = frame_colors.get(frame, (0.5, 0.5, 0.5, 1.0))
        n = len(mutations)

        for i, r in enumerate(mutations):
            ha = _safe_float(r.get("HMBR_A", "")) or 0.0
            hd = _safe_float(r.get("HMBR_D", "")) or 0.0
            ax.bar(i - bar_w / 2, math.sqrt(ha), width=bar_w,
                   color=color, alpha=0.9, zorder=2, linewidth=0)
            ax.bar(i + bar_w / 2, math.sqrt(hd), width=bar_w,
                   color=color, alpha=0.45, zorder=2, linewidth=0)

        ax.axhline(threshold_sqrt, linestyle="--", linewidth=0.6,
                   color="#555555", zorder=3)
        ax.set_xlim(-0.65, n - 0.35)
        ax.set_ylim(0, sqrt_max * 1.08)

        ax.set_xticks(range(n))
        ax.set_xticklabels(
            [r.get("Mutation", "") for r in mutations],
            rotation=90, ha="right", va="top", fontsize=7,
        )

        # Custom sqrt-scale y-ticks: show nice linear values
        nice_lin = [v for v in (0, 1, 2, 5, 10, 20, 50, 100, 200)
                    if math.sqrt(v) <= sqrt_max * 1.08]
        ax.set_yticks([math.sqrt(v) for v in nice_lin])
        if col_idx == 0:
            ax.set_yticklabels([str(v) for v in nice_lin], fontsize=8)
        else:
            ax.set_yticklabels([])

        _apply_facet_style(ax, col_idx, frame, "HMBR  (√ scale)")

    # Legend: Ancestral / Derived
    legend_handles = [
        Patch(facecolor="grey", alpha=0.9, label="Ancestral"),
        Patch(facecolor="grey", alpha=0.45, label="Derived"),
    ]
    axes[-1].legend(handles=legend_handles, fontsize=8, loc="upper right",
                    frameon=True, framealpha=0.9, edgecolor="#cccccc")

    fig.suptitle("Harmonic Mean Best Rank — Ancestral vs Derived",
                 fontsize=10, y=0.97)
    return fig


def _render_log2fc(rows: List[Dict]) -> "Figure":
    """
    Plot B: bar chart of log₂ fold-change in HMBR per mutation,
    faceted by Frame.
    """
    frame_order, frame_colors = _frame_order_and_colors(rows)
    grouped = _mutations_by_frame(rows, frame_order, ["log2_foldchange_HMBR"])

    fig, axes, active = _make_faceted_figure(frame_order, grouped, fig_height=5.0)

    all_fc = [_safe_float(r.get("log2_foldchange_HMBR", ""))
              for r in rows if _safe_float(r.get("log2_foldchange_HMBR", "")) is not None]
    y_lim = max(abs(v) for v in all_fc) * 1.1 if all_fc else 3.0

    bar_w = 0.65

    for col_idx, (ax, frame) in enumerate(zip(axes, active)):
        mutations = grouped[frame]
        color = frame_colors.get(frame, (0.5, 0.5, 0.5, 1.0))
        n = len(mutations)

        for i, r in enumerate(mutations):
            fc = _safe_float(r.get("log2_foldchange_HMBR", ""))
            if fc is not None:
                ax.bar(i, fc, width=bar_w, color=color, alpha=0.88,
                       zorder=2, linewidth=0)

        ax.axhline(0, linestyle="--", linewidth=0.6, color="#555555", zorder=3)
        ax.set_xlim(-0.65, n - 0.35)
        ax.set_ylim(-y_lim, y_lim)

        ax.set_xticks(range(n))
        ax.set_xticklabels(
            [r.get("Mutation", "") for r in mutations],
            rotation=90, ha="right", va="top", fontsize=7,
        )
        _apply_facet_style(ax, col_idx, frame,
                           "log₂ Fold Change HMBR")

    fig.suptitle("log₂ Fold Change in Harmonic Mean Best Rank",
                 fontsize=10, y=0.97)
    return fig


def _render_simulated(obs_rows: List[Dict], sim_rows: List[Dict]) -> "Figure":
    """
    Plot sim: for each enriched observed mutation, plot a kernel-density
    estimate of the simulated null distribution and mark the observed value
    with its percentile.  Arranged in a 3-row grid (matching notebook).
    """
    # Filter observed to enriched mutations only
    enriched = []
    for r in obs_rows:
        fc = _safe_float(r.get("log2_foldchange_HMBR", ""))
        ha = _safe_float(r.get("HMBR_A", ""))
        hd = _safe_float(r.get("HMBR_D", ""))
        # Only show enriched (FC > 0) and not trivially high-rank variants
        if fc is not None and fc > 0:
            if not (ha is not None and hd is not None and ha > 2 and hd > 2):
                enriched.append(r)

    if not enriched:
        raise ValueError(
            "No enriched mutations found (log₂FC > 0).\n"
            "Run with percentile analysis enabled to see the simulated background."
        )

    # Simulated log2FC values
    sim_vals = [_safe_float(r.get("log2_foldchange_HMBR", ""))
                for r in sim_rows]
    sim_vals = np.array([v for v in sim_vals if v is not None], dtype=float)
    if len(sim_vals) < 10:
        raise ValueError("Insufficient simulated data for density estimate.")

    # Sort enriched by Locus
    enriched.sort(key=lambda r: _safe_float(r.get("Locus", "")) or 0)

    # Frame → color
    frame_order, frame_colors = _frame_order_and_colors(obs_rows)

    # KDE using Silverman's rule
    bw = 1.06 * float(np.std(sim_vals)) * len(sim_vals) ** (-0.2)
    n_kde = 512
    x_full = np.linspace(sim_vals.min() - 3 * bw, sim_vals.max() + 3 * bw, n_kde)

    def _kde_at(x_pts: np.ndarray) -> np.ndarray:
        diff = x_pts[:, None] - sim_vals[None, :]          # (n_x, n_sim)
        kernel = np.exp(-0.5 * (diff / bw) ** 2)
        return kernel.mean(axis=1) / (bw * math.sqrt(2 * math.pi))

    y_full = _kde_at(x_full)
    global_y_max = float(y_full.max())

    # ECDF for percentile annotation
    sorted_sim = np.sort(sim_vals)

    def _percentile(v: float) -> float:
        return float(np.searchsorted(sorted_sim, v, side="right")) / len(sorted_sim) * 100

    # Layout: up to 3 rows, columns fill as needed
    n_mut = len(enriched)
    n_rows = min(3, n_mut)
    n_cols = math.ceil(n_mut / n_rows)
    fig_width = max(10, n_cols * 3.5)
    fig_height = n_rows * 3.0 + 0.8

    fig = Figure(figsize=(fig_width, fig_height), dpi=100)
    fig.patch.set_facecolor("#f5f5f7")

    for idx, r in enumerate(enriched):
        row_i = idx % n_rows
        col_i = idx // n_rows
        ax = fig.add_subplot(n_rows, n_cols, row_i * n_cols + col_i + 1)

        frame = r.get("Frame", "")
        obs_fc = _safe_float(r.get("log2_foldchange_HMBR", "")) or 0.0
        mutation = r.get("Mutation", "")
        color = frame_colors.get(frame, (0.5, 0.5, 0.5, 1.0))

        # Window the density around the observed value
        left_margin, right_margin = 2.5, 3.5
        x_window = x_full[(x_full >= obs_fc - left_margin) &
                           (x_full <= obs_fc + right_margin)]
        y_window = _kde_at(x_window)

        ax.fill_between(x_window, y_window, alpha=0.75, color=color, linewidth=0)
        ax.axvline(obs_fc, linestyle="--", linewidth=1.0, color="#333333", zorder=3)

        pct = _percentile(obs_fc)
        abbrev = FRAME_ABBREVS.get(frame, frame)
        label = f"{mutation} ({abbrev})\n{obs_fc:.3g}  ({pct:.3g}th pctile)"
        ax.text(obs_fc + 0.12, global_y_max * 0.82, label,
                fontsize=7, ha="left", va="top", color="#222222",
                bbox=dict(facecolor="white", alpha=0.7, edgecolor="none", pad=1))

        ax.set_xlabel("log₂FC HMBR", fontsize=7)
        ax.set_ylabel("Density", fontsize=7) if col_i == 0 else None
        ax.set_facecolor("white")
        ax.tick_params(labelsize=7)
        for spine in ax.spines.values():
            spine.set_linewidth(0.4)
            spine.set_color("#aaaaaa")

    fig.suptitle("Simulated Null Distribution — Observed Values Marked",
                 fontsize=10, y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    return fig


def _render_per_allele_scatter(
    pa_rows: List[Dict],
    frame_order: List[str],
) -> "Figure":
    """
    Per-allele log₂ fold-change scatter, faceted by Frame.
    Each point is one (mutation, allele) pair, coloured by HLA locus.
    """
    # Enrich rows with derived fields
    enriched: List[Dict] = []
    for r in pa_rows:
        fc = _safe_float(r.get("log2_foldchange_BR", ""))
        locus = _safe_float(r.get("Locus", ""))
        mhc = r.get("MHC", "")
        if fc is None or locus is None:
            continue
        allele_short = mhc.replace("HLA-", "").replace("HLA*", "")
        hla_locus = ""
        for ch in mhc:
            if ch in ("A", "B", "C"):
                hla_locus = ch
                break
        enriched.append({**r, "log2_foldchange_BR": fc, "Locus": locus,
                         "allele_short": allele_short, "HLA_locus": hla_locus})

    if not enriched:
        raise ValueError("No valid per-allele data found.")

    grouped = _mutations_by_frame(
        enriched, frame_order, ["log2_foldchange_BR"]
    )

    fig, axes, active = _make_faceted_figure(
        frame_order, grouped, fig_height=5.0, mutation_width=0.60
    )

    all_fc = [r["log2_foldchange_BR"] for r in enriched
              if isinstance(r["log2_foldchange_BR"], float)]
    y_lim = max(abs(v) for v in all_fc) * 1.1 if all_fc else 3.0

    rng = np.random.default_rng(42)

    for col_idx, (ax, frame) in enumerate(zip(axes, active)):
        mutations = grouped[frame]
        # Group points per mutation position
        mut_positions: Dict[str, int] = {}
        mut_list = []
        for r in mutations:
            mut = r.get("Mutation", "")
            if mut not in mut_positions:
                mut_positions[mut] = len(mut_list)
                mut_list.append(mut)

        for r in mutations:
            mut = r.get("Mutation", "")
            fc = r["log2_foldchange_BR"]
            hla_locus = r["HLA_locus"]
            x = mut_positions[mut] + rng.uniform(-0.15, 0.15)
            color = LOCUS_COLORS.get(hla_locus, LOCUS_FALLBACK)
            ax.scatter(x, fc, s=28, color=color, alpha=0.85, zorder=2,
                       linewidths=0)

        ax.axhline(0, linestyle="--", linewidth=0.6, color="#555555", zorder=3)
        n = len(mut_list)
        ax.set_xlim(-0.65, n - 0.35)
        ax.set_ylim(-y_lim, y_lim)
        ax.set_xticks(range(n))
        ax.set_xticklabels(mut_list, rotation=90, ha="right", va="top", fontsize=7)
        _apply_facet_style(ax, col_idx, frame, "log₂FC Best Rank (per allele)")

    # HLA-locus legend
    handles = [Patch(facecolor=LOCUS_COLORS.get(l, LOCUS_FALLBACK), label=l)
               for l in ("A", "B", "C")]
    axes[-1].legend(handles=handles, title="HLA locus", fontsize=8,
                    title_fontsize=8, loc="upper right",
                    frameon=True, framealpha=0.9, edgecolor="#cccccc")

    fig.suptitle("Per-Allele log₂ Fold Change in Best Rank",
                 fontsize=10, y=0.97)
    return fig


def _render_per_allele_box(pa_rows: List[Dict]) -> "Figure":
    """
    Boxplot of per-allele log₂ fold-change grouped by allele,
    ordered by HLA locus then allele name.
    """
    enriched: List[Dict] = []
    for r in pa_rows:
        fc = _safe_float(r.get("log2_foldchange_BR", ""))
        mhc = r.get("MHC", "")
        if fc is None:
            continue
        allele_short = mhc.replace("HLA-", "").replace("HLA*", "")
        hla_locus = ""
        for ch in mhc:
            if ch in ("A", "B", "C"):
                hla_locus = ch
                break
        enriched.append({"fc": fc, "allele_short": allele_short,
                         "HLA_locus": hla_locus})

    if not enriched:
        raise ValueError("No valid per-allele data found.")

    # Allele order: locus then name
    allele_meta = {}
    for r in enriched:
        a = r["allele_short"]
        if a not in allele_meta:
            allele_meta[a] = r["HLA_locus"]
    allele_order = sorted(allele_meta, key=lambda a: (allele_meta[a], a))

    # Collect values per allele
    by_allele: Dict[str, List[float]] = {a: [] for a in allele_order}
    for r in enriched:
        by_allele[r["allele_short"]].append(r["fc"])

    n_alleles = len(allele_order)
    fig_width = max(8, n_alleles * 0.9 + 2)
    fig = Figure(figsize=(fig_width, 5.0), dpi=100)
    fig.patch.set_facecolor("#f5f5f7")
    ax = fig.add_subplot(111)

    rng = np.random.default_rng(42)
    all_fc = [r["fc"] for r in enriched]
    y_lim = max(abs(v) for v in all_fc) * 1.1 if all_fc else 3.0

    for i, allele in enumerate(allele_order):
        vals = np.array(by_allele[allele])
        locus = allele_meta[allele]
        color = LOCUS_COLORS.get(locus, LOCUS_FALLBACK)

        # Boxplot (no fliers — we draw points ourselves)
        bp = ax.boxplot(
            vals, positions=[i], widths=0.5,
            patch_artist=True, showfliers=False,
            boxprops=dict(facecolor="none", color="#444444", linewidth=0.8),
            whiskerprops=dict(color="#444444", linewidth=0.8),
            capprops=dict(color="#444444", linewidth=0.8),
            medianprops=dict(color="#333333", linewidth=1.2),
        )

        # Jittered points
        jitter = rng.uniform(-0.15, 0.15, size=len(vals))
        ax.scatter(i + jitter, vals, s=22, color=color, alpha=0.7,
                   zorder=3, linewidths=0)

    ax.axhline(0, linestyle="--", linewidth=0.6, color="#555555", zorder=2)
    ax.set_xlim(-0.7, n_alleles - 0.3)
    ax.set_ylim(-y_lim, y_lim)
    ax.set_xticks(range(n_alleles))
    ax.set_xticklabels(allele_order, rotation=90, ha="right", va="top", fontsize=8)
    ax.set_ylabel("log₂ Fold Change Best Rank", fontsize=9)
    ax.set_xlabel("Allele", fontsize=9)
    ax.set_facecolor("white")
    ax.tick_params(labelsize=8)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
        spine.set_color("#aaaaaa")
    ax.grid(axis="y", linewidth=0.3, color="#eeeeee", zorder=0)

    # Locus legend
    handles = [Patch(facecolor=LOCUS_COLORS.get(l, LOCUS_FALLBACK), label=l)
               for l in ("A", "B", "C")]
    ax.legend(handles=handles, title="HLA locus", fontsize=8,
              title_fontsize=8, loc="upper right",
              frameon=True, framealpha=0.9, edgecolor="#cccccc")

    fig.suptitle("Per-Allele log₂ Fold Change — Distribution by Allele",
                 fontsize=10, y=0.97)
    fig.subplots_adjust(left=0.10, right=0.97, bottom=0.28, top=0.91)
    return fig


# ── Scrollable canvas widget ──────────────────────────────────────────────────

class _PlotCanvas(QScrollArea):
    """A scrollable area containing one FigureCanvas."""

    def __init__(self, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self.setFrameShape(QFrame.Shape.NoFrame)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self._fig: Optional["Figure"] = None
        self._show_placeholder("No data loaded.")

    def set_figure(self, fig: "Figure") -> None:
        if self._fig is not None:
            try:
                self._fig.clf()
            except Exception:
                pass
        self._fig = fig
        canvas = FigureCanvas(fig)
        canvas.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding
        )
        # Minimum size from figure dimensions so scrollbars appear as needed
        w_px = int(fig.get_figwidth() * fig.dpi)
        h_px = int(fig.get_figheight() * fig.dpi)
        canvas.setMinimumSize(w_px, h_px)
        self.setWidget(canvas)
        self.setWidgetResizable(False)
        canvas.draw()

    def show_error(self, msg: str) -> None:
        lbl = QLabel(msg)
        lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
        lbl.setWordWrap(True)
        lbl.setStyleSheet("color: #ff3b30; font-size: 12px; padding: 24px;")
        self._set_label(lbl)

    def _show_placeholder(self, msg: str) -> None:
        lbl = QLabel(msg)
        lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
        lbl.setStyleSheet("color: #8e8e93; font-size: 12px; padding: 24px;")
        self._set_label(lbl)

    def _set_label(self, lbl: QLabel) -> None:
        self.setWidget(lbl)
        self.setWidgetResizable(True)


# ── Main PlotViewer widget ────────────────────────────────────────────────────

class PlotViewer(QWidget):
    """
    Tabbed plot viewer for CD8scape output files.

    Call load() after each analysis run to populate the plots.
    """

    TAB_HMBR     = 0
    TAB_LOG2FC   = 1
    TAB_SIM      = 2
    TAB_PA_SCAT  = 3
    TAB_PA_BOX   = 4

    def __init__(self, parent: Optional[QWidget] = None):
        super().__init__(parent)
        vbox = QVBoxLayout(self)
        vbox.setContentsMargins(0, 0, 0, 0)
        vbox.setSpacing(0)

        if not _MPL_OK:
            notice = QLabel(
                "matplotlib is not installed.\n\n"
                "Run:  pip install matplotlib\n"
                "then restart CD8scape to enable plot output."
            )
            notice.setAlignment(Qt.AlignmentFlag.AlignCenter)
            notice.setStyleSheet(
                "color: #636366; font-size: 13px; padding: 32px;"
            )
            vbox.addWidget(notice)
            return

        self._tabs = QTabWidget()
        self._tabs.setDocumentMode(True)

        self._canvas_hmbr    = _PlotCanvas()
        self._canvas_log2fc  = _PlotCanvas()
        self._canvas_sim     = _PlotCanvas()
        self._canvas_pa_scat = _PlotCanvas()
        self._canvas_pa_box  = _PlotCanvas()

        self._tabs.addTab(self._canvas_hmbr,    "HMBR")
        self._tabs.addTab(self._canvas_log2fc,  "log₂ FC")
        self._tabs.addTab(self._canvas_sim,     "Simulated")
        self._tabs.addTab(self._canvas_pa_scat, "Per-Allele")
        self._tabs.addTab(self._canvas_pa_box,  "Allele Box")

        vbox.addWidget(self._tabs)
        self._reset_placeholders()

    # ── Public API ────────────────────────────────────────────────────────────

    def load(
        self,
        folder: Optional[Path],
        run_suffix: str = "",
        sim_suffix: str = "simulated",
        has_per_allele: bool = False,
        has_percentile: bool = False,
    ) -> None:
        """
        Locate output CSVs inside *folder* and render all available plots.

        Parameters
        ----------
        folder        : data folder used for this run
        run_suffix    : --suffix value used for the observed run (may be "")
        sim_suffix    : --suffix used for the simulate step (default "simulated")
        has_per_allele: True if --per-allele was requested
        has_percentile: True if percentile analysis was run
        """
        if not _MPL_OK:
            return
        if folder is None or not folder.is_dir():
            self._reset_placeholders("No output folder available.")
            return

        # ── Locate files ──────────────────────────────────────────────────────
        hmbr_path    = _find_file(folder, "harmonic_mean_best_ranks", run_suffix)
        sim_path     = _find_file(folder, "harmonic_mean_best_ranks", sim_suffix) \
                       if has_percentile else None
        pa_path      = _find_file(folder, "per_allele_best_ranks", run_suffix) \
                       if has_per_allele else None

        # ── Load observed HMBR ────────────────────────────────────────────────
        obs_rows: List[Dict] = []
        if hmbr_path and hmbr_path.exists():
            try:
                obs_rows = _load_csv(hmbr_path)
            except Exception as exc:
                self._set_all_error(f"Could not read {hmbr_path.name}:\n{exc}")
                return
        else:
            self._reset_placeholders(
                "harmonic_mean_best_ranks.csv not found.\n"
                "Run an analysis first."
            )
            return

        # Pre-compute frame order from observed data (used by per-allele plots)
        frame_order, _ = _frame_order_and_colors(obs_rows)

        # ── HMBR plot ─────────────────────────────────────────────────────────
        self._render_tab(self._canvas_hmbr, _render_hmbr, obs_rows)

        # ── log₂ FC plot ──────────────────────────────────────────────────────
        self._render_tab(self._canvas_log2fc, _render_log2fc, obs_rows)

        # ── Simulated density plot ────────────────────────────────────────────
        if has_percentile and sim_path and sim_path.exists():
            try:
                sim_rows = _load_csv(sim_path)
            except Exception as exc:
                self._canvas_sim.show_error(
                    f"Could not read {sim_path.name}:\n{exc}"
                )
                sim_rows = []
            if sim_rows:
                self._render_tab(
                    self._canvas_sim, _render_simulated, obs_rows, sim_rows
                )
        else:
            self._canvas_sim._show_placeholder(
                "Simulated background not available.\n"
                "Enable percentile analysis to see this plot."
            )

        # ── Per-allele plots ──────────────────────────────────────────────────
        if has_per_allele and pa_path and pa_path.exists():
            try:
                pa_rows = _load_csv(pa_path)
            except Exception as exc:
                msg = f"Could not read {pa_path.name}:\n{exc}"
                self._canvas_pa_scat.show_error(msg)
                self._canvas_pa_box.show_error(msg)
                pa_rows = []
            if pa_rows:
                self._render_tab(
                    self._canvas_pa_scat,
                    _render_per_allele_scatter,
                    pa_rows, frame_order,
                )
                self._render_tab(self._canvas_pa_box, _render_per_allele_box,
                                 pa_rows)
        else:
            msg = (
                "Per-allele data not available.\n"
                "Enable --per-allele in Run options to see this plot."
            )
            self._canvas_pa_scat._show_placeholder(msg)
            self._canvas_pa_box._show_placeholder(msg)

    # ── Internal helpers ──────────────────────────────────────────────────────

    def _render_tab(
        self,
        canvas: _PlotCanvas,
        renderer,
        *args,
    ) -> None:
        """Call renderer(*args) and push the resulting Figure into canvas."""
        try:
            fig = renderer(*args)
            canvas.set_figure(fig)
        except Exception as exc:
            canvas.show_error(str(exc))

    def _reset_placeholders(self, msg: str = "Run an analysis to see plots.") -> None:
        if not _MPL_OK:
            return
        for canvas in (
            self._canvas_hmbr,
            self._canvas_log2fc,
            self._canvas_sim,
            self._canvas_pa_scat,
            self._canvas_pa_box,
        ):
            canvas._show_placeholder(msg)

    def _set_all_error(self, msg: str) -> None:
        for canvas in (
            self._canvas_hmbr,
            self._canvas_log2fc,
            self._canvas_sim,
            self._canvas_pa_scat,
            self._canvas_pa_box,
        ):
            canvas.show_error(msg)
