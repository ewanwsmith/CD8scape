"""
plots_widget.py — Interactive plot viewer for CD8scape outputs.

Uses pyqtgraph (a PyQt-native library) for all rendering — trackpad scroll,
pinch-to-zoom, and drag-to-pan all work without any custom event handling.

Three conditional plot tabs, mirroring the reference notebook (plot_output.ipynb):

  • Escape Scores         — log₂ fold-change in HMBR per mutation, faceted by
                            reading frame.  Always shown after a run.

  • Percentile Scores     — per-mutation density of the simulated null
                            distribution with the observed value marked.
                            Shown only when percentile analysis was run.

  • Per-Allele Escape     — per-allele log₂ fold-change scatter (Per Mutation
    Scores                  sub-tab) and boxplot by allele (By Allele sub-tab).
                            Shown only when --per-allele was requested.

Saving: each plot panel has a "Save…" button that writes a PNG or SVG with a
standardised filename:
    CD8scape_escape_scores_{folder}_{suffix}.png
    CD8scape_percentile_scores_{folder}_{suffix}.png
    CD8scape_per_allele_escape_scores_{folder}_{suffix}.png

COLUMN_TOOLTIPS is imported by app_qt.py to add hover-over explanations to
the file-preview table in the Files tab.

Requires: pyqtgraph>=0.13, numpy>=1.24
"""
from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QFileDialog,
    QFrame,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

# ── pyqtgraph — lazy import so pip-install mid-session works immediately ───────
_PG_OK = False
pg = None  # type: ignore[assignment]

def _ensure_pg() -> bool:
    """Try to import pyqtgraph (and its exporters).  Safe to call repeatedly."""
    global _PG_OK, pg
    if _PG_OK:
        return True
    try:
        import importlib, importlib.util
        importlib.invalidate_caches()
        import pyqtgraph as _pg
        import pyqtgraph.exporters  # noqa: F401 — registers exporters
        _pg.setConfigOptions(antialias=True, background="#ffffff")
        pg = _pg
        _PG_OK = True
    except ImportError:
        pass
    return _PG_OK

_ensure_pg()  # attempt at startup; succeeds if already installed

# ── Column tooltips exported to app_qt.py ─────────────────────────────────────
COLUMN_TOOLTIPS: Dict[str, str] = {
    "Frame":
        "Reading frame / protein the variant falls in.",
    "Locus":
        "Genomic position of the variant (nucleotide index).",
    "Mutation":
        "Amino-acid change or indel description.",
    "HMBR_A":
        "Harmonic Mean Best Rank — ancestral allele.\n"
        "Lower values mean a peptide is predicted to bind better.",
    "HMBR_D":
        "Harmonic Mean Best Rank — derived (mutant) allele.\n"
        "Lower values mean a peptide is predicted to bind better.",
    "foldchange_HMBR":
        "Derived HMBR ÷ ancestral HMBR.\n"
        ">1 → escape (derived binds worse); <1 → gained binding.",
    "log2_foldchange_HMBR":
        "log₂(derived HMBR / ancestral HMBR).\n"
        "Positive → immune escape (derived binds worse).\n"
        "Negative → increased predicted binding after mutation.",
    "MHC":
        "HLA allele used for this binding prediction.",
    "ELBR_A":
        "Expected Log Best Rank — ancestral allele (per-allele metric).",
    "ELBR_D":
        "Expected Log Best Rank — derived allele (per-allele metric).",
    "foldchange_BR":
        "Derived best rank ÷ ancestral best rank for this allele.\n"
        ">1 → escape; <1 → gained binding.",
    "log2_foldchange_BR":
        "log₂(derived / ancestral) best rank for this allele.\n"
        "Positive → escape; Negative → increased binding.",
    "Percentile":
        "Enrichment percentile relative to the simulated background.\n"
        "Higher → more significant predicted immune escape.",
    "max_escape_allele":
        "The HLA allele showing the largest predicted escape for this mutation.",
    "max_escape_log2_fc":
        "log₂ fold change for the max-escape allele.",
}

# ── Protein-name abbreviations ────────────────────────────────────────────────
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

# Viridis-sampled palette (covers up to 20 reading frames)
_VIRIDIS_HEX = [
    "#440154", "#46085c", "#470d60", "#481769", "#482273",
    "#472d7b", "#453882", "#414287", "#3d4d8a", "#38578c",
    "#33618d", "#2e6b8e", "#29758e", "#257f8e", "#218a8d",
    "#1d948a", "#1a9e86", "#1fa87f", "#2db27d", "#41bc77",
]

# HLA-locus palette (approximates viridis mako)
LOCUS_COLORS: Dict[str, str] = {
    "A": "#2d1160",
    "B": "#2272b5",
    "C": "#6ecdc8",
}
_LOCUS_FALLBACK = "#888888"


# ── Data helpers ──────────────────────────────────────────────────────────────

def _safe_float(val: Any) -> Optional[float]:
    try:
        f = float(val)
        return None if (math.isnan(f) or math.isinf(f)) else f
    except (ValueError, TypeError):
        return None


def _load_csv(path: Path) -> List[Dict[str, str]]:
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        return [dict(row) for row in csv.DictReader(fh)]


def _find_file(folder: Path, stem: str, suffix: str) -> Optional[Path]:
    if suffix:
        p = folder / f"{stem}_{suffix}.csv"
        if p.exists():
            return p
    p = folder / f"{stem}.csv"
    return p if p.exists() else None


def _frame_order_and_colors(
    rows: List[Dict],
) -> Tuple[List[str], Dict[str, str], Dict[str, str]]:
    """
    Return (frame_order, frame_colors, mutation_colors).

    frame_colors   maps frame  → hex color (one viridis colour per frame).
    mutation_colors maps mutation → hex color, populated only when all rows
                   share a single frame — each mutation gets its own colour
                   spread evenly across the full viridis palette.
    When multiple frames are present, mutation_colors is empty and callers
    should colour by frame.
    """
    min_locus: Dict[str, float] = {}
    for r in rows:
        f = r.get("Frame", "")
        locus = _safe_float(r.get("Locus", ""))
        if f and locus is not None:
            if f not in min_locus or locus < min_locus[f]:
                min_locus[f] = locus
    order = sorted(min_locus, key=lambda f: min_locus[f])
    n = len(order)
    palette = _VIRIDIS_HEX[:n] if n <= len(_VIRIDIS_HEX) else (
        _VIRIDIS_HEX * (n // len(_VIRIDIS_HEX) + 1)
    )[:n]
    frame_colors = dict(zip(order, palette))

    # Per-mutation colours when all data is in one frame
    mutation_colors: Dict[str, str] = {}
    if n <= 1:
        muts: List[str] = []
        seen: set = set()
        for r in rows:
            m = r.get("Mutation", "")
            if m and m not in seen:
                seen.add(m)
                muts.append(m)
        nm = len(muts)
        if nm > 0:
            v = _VIRIDIS_HEX
            indices = [round(i * (len(v) - 1) / max(1, nm - 1)) for i in range(nm)]
            mutation_colors = {m: v[idx] for m, idx in zip(muts, indices)}

    return order, frame_colors, mutation_colors


def _mutations_by_frame(
    rows: List[Dict],
    frame_order: List[str],
    value_cols: List[str],
) -> Dict[str, List[Dict]]:
    grouped: Dict[str, List[Dict]] = {f: [] for f in frame_order}
    for r in rows:
        f = r.get("Frame", "")
        if f not in grouped:
            continue
        locus = _safe_float(r.get("Locus", ""))
        if locus is None:
            continue
        if not any(_safe_float(r.get(c, "")) is not None for c in value_cols):
            continue
        grouped[f].append(r)
    for f in grouped:
        grouped[f].sort(key=lambda r: _safe_float(r.get("Locus", "")) or 0)
    return grouped


def _gaussian_kde(data: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Silverman-bandwidth Gaussian KDE."""
    bw = 1.06 * float(np.std(data)) * len(data) ** (-0.2)
    diff = x[:, None] - data[None, :]
    kernel = np.exp(-0.5 * (diff / bw) ** 2)
    return kernel.mean(axis=1) / (bw * math.sqrt(2 * math.pi))


# ── pyqtgraph styling helpers ─────────────────────────────────────────────────

def _hex_to_rgba(h: str, alpha: int = 200) -> Tuple[int, int, int, int]:
    h = h.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return (r, g, b, alpha)


def _style_plot(p: "pg.PlotItem", title: str = "") -> None:
    """Apply shared visual style to a pyqtgraph PlotItem."""
    if title:
        p.setTitle(title, color="#1d1d1f", size="9pt")
    p.getAxis("bottom").setPen(pg.mkPen("#aaaaaa", width=0.5))
    p.getAxis("left").setPen(pg.mkPen("#aaaaaa", width=0.5))
    p.getAxis("bottom").setTextPen(pg.mkPen("#333333"))
    p.getAxis("left").setTextPen(pg.mkPen("#333333"))
    p.showGrid(x=False, y=True, alpha=0.12)
    # Mouse interaction is disabled per-ViewBox; zoom/pan is handled by
    # the outer _ZoomScrollArea handles zoom/pan at the pane level.
    p.getViewBox().setMouseMode(pg.ViewBox.PanMode)


# ── Plot builders ─────────────────────────────────────────────────────────────

def _build_escape_scores(
    rows: List[Dict],
    frame_order: List[str],
    frame_colors: Dict[str, str],
    mutation_colors: Dict[str, str] = {},
) -> "_PassthroughGLW":
    """
    Bar chart: log₂ fold-change in HMBR per mutation, faceted by frame.
    All mutations shown; dashed line at 0.
    """
    grouped = _mutations_by_frame(rows, frame_order, ["log2_foldchange_HMBR"])
    active = [f for f in frame_order if grouped.get(f)]
    if not active:
        raise ValueError("No log₂ fold-change data found.")

    all_fc = [_safe_float(r.get("log2_foldchange_HMBR", ""))
              for r in rows if _safe_float(r.get("log2_foldchange_HMBR", "")) is not None]
    y_lim = max(abs(v) for v in all_fc) * 1.12 if all_fc else 3.0

    glw = _PassthroughGLW()
    glw.setBackground("#f5f5f7")

    first_plot = None
    for col_i, frame in enumerate(active):
        mutations = grouped[frame]
        n = len(mutations)
        abbrev = FRAME_ABBREVS.get(frame, frame)

        p = glw.addPlot(row=0, col=col_i, title=abbrev,
                        axisItems={"bottom": _RotatedAxisItem()})
        _style_plot(p)
        p.getViewBox().setMouseEnabled(x=False, y=False)
        p.setTitle(abbrev, color="#1d1d1f", size="9pt",
                   bold=True, italic=False)

        # Proportional column widths
        glw.ci.layout.setColumnStretchFactor(col_i, max(1, n))

        heights = np.array([
            _safe_float(r.get("log2_foldchange_HMBR", "")) or 0.0
            for r in mutations
        ])
        x = np.arange(n, dtype=float)

        if mutation_colors:
            # Single-frame run: one viridis colour per mutation
            for i, r in enumerate(mutations):
                c_hex = mutation_colors.get(r.get("Mutation", ""), "#440154")
                bar = pg.BarGraphItem(
                    x=np.array([x[i]]), height=np.array([heights[i]]),
                    width=0.65,
                    brush=pg.mkBrush(*_hex_to_rgba(c_hex, 210)),
                    pen=pg.mkPen(None),
                )
                p.addItem(bar)
        else:
            color_hex = frame_colors.get(frame, "#440154")
            p.addItem(pg.BarGraphItem(
                x=x, height=heights, width=0.65,
                brush=pg.mkBrush(*_hex_to_rgba(color_hex, 210)),
                pen=pg.mkPen(None),
            ))

        # y = 0 dashed line
        zero = pg.InfiniteLine(
            pos=0, angle=0,
            pen=pg.mkPen("#555555", width=1,
                         style=Qt.PenStyle.DashLine),
        )
        p.addItem(zero)

        # x-axis: mutation names (rotated via _RotatedAxisItem)
        ticks = [[(i, mutations[i].get("Mutation", "")) for i in range(n)]]
        p.getAxis("bottom").setTicks(ticks)

        p.setXRange(-0.6, n - 0.4, padding=0)
        p.setYRange(-y_lim, y_lim, padding=0)

        if col_i == 0:
            p.setLabel("left", "log₂ FC HMBR", color="#333333", size="8pt")
            first_plot = p
        else:
            p.hideAxis("left")
            if first_plot:
                p.setYLink(first_plot)

    return glw


def _build_percentile_scores(
    obs_rows: List[Dict],
    sim_rows: List[Dict],
    frame_colors: Dict[str, str],
    mutation_colors: Dict[str, str] = {},
) -> "_PassthroughGLW":
    """
    One density panel per mutation (all mutations with valid data),
    showing the simulated null KDE with the observed log₂FC marked.
    """
    sim_vals = np.array(
        [v for r in sim_rows
         if (v := _safe_float(r.get("log2_foldchange_HMBR", ""))) is not None],
        dtype=float,
    )
    if len(sim_vals) < 10:
        raise ValueError("Insufficient simulated data for density estimate.")

    sorted_sim = np.sort(sim_vals)

    def _pctile(v: float) -> float:
        return float(np.searchsorted(sorted_sim, v, side="right")) / len(sorted_sim) * 100

    # All mutations with valid log2FC (no positivity filter)
    candidates = []
    for r in obs_rows:
        fc = _safe_float(r.get("log2_foldchange_HMBR", ""))
        locus = _safe_float(r.get("Locus", ""))
        if fc is not None and locus is not None:
            candidates.append(r)
    candidates.sort(key=lambda r: _safe_float(r.get("Locus", "")) or 0)

    if not candidates:
        raise ValueError("No mutations with valid fold-change data found.")

    n_mut = len(candidates)
    n_cols = max(1, math.ceil(n_mut / 3))
    n_rows = min(3, n_mut)

    # KDE grid on the full sim distribution
    x_full = np.linspace(sim_vals.min() - 1.5, sim_vals.max() + 1.5, 512)
    y_full = _gaussian_kde(sim_vals, x_full)
    y_max = float(y_full.max())

    glw = _PassthroughGLW()
    glw.setBackground("#f5f5f7")

    for idx, r in enumerate(candidates):
        col_i = idx // n_rows
        row_i = idx % n_rows

        frame = r.get("Frame", "")
        obs_fc = _safe_float(r.get("log2_foldchange_HMBR", "")) or 0.0
        mutation = r.get("Mutation", "")
        abbrev = FRAME_ABBREVS.get(frame, frame)
        pct = _pctile(obs_fc)

        title = f"{mutation}  ({abbrev})"
        p = glw.addPlot(row=row_i, col=col_i, title=title)
        _style_plot(p)
        p.getViewBox().setMouseEnabled(x=False, y=False)
        p.setTitle(title, color="#1d1d1f", size="8pt")

        if mutation_colors:
            color_hex = mutation_colors.get(mutation, "#440154")
        else:
            color_hex = frame_colors.get(frame, "#440154")
        rgba = _hex_to_rgba(color_hex, 180)

        # Window the density near the observed value
        left_m, right_m = 2.5, 3.5
        mask = (x_full >= obs_fc - left_m) & (x_full <= obs_fc + right_m)
        x_win = x_full[mask]
        y_win = y_full[mask]

        # Filled density curve
        curve = p.plot(
            x_win, y_win,
            fillLevel=0,
            brush=pg.mkBrush(*rgba),
            pen=pg.mkPen(color_hex, width=1),
        )

        # Observed value vertical line
        obs_line = pg.InfiniteLine(
            pos=obs_fc, angle=90,
            pen=pg.mkPen("#333333", width=1.5,
                         style=Qt.PenStyle.DashLine),
            label=f"{obs_fc:.2g}\n({pct:.1f}th %ile)",
            labelOpts=dict(
                position=0.80,
                color="#1d1d1f",
                movable=False,
                fill=pg.mkBrush(255, 255, 255, 180),
            ),
        )
        p.addItem(obs_line)

        p.setLabel("bottom", "log₂ FC HMBR", color="#555555", size="7pt")
        if col_i == 0:
            p.setLabel("left", "Density", color="#555555", size="7pt")
        else:
            p.hideAxis("left")

    return glw


def _build_per_allele_scatter(
    pa_rows: List[Dict],
    frame_order: List[str],
    frame_colors: Dict[str, str],
    mutation_colors: Dict[str, str] = {},
) -> "_PassthroughGLW":
    """
    Scatter: per-allele log₂ fold-change per mutation, faceted by frame,
    coloured by HLA locus (A / B / C).
    """
    # Enrich rows
    enriched = []
    for r in pa_rows:
        fc = _safe_float(r.get("log2_foldchange_BR", ""))
        locus = _safe_float(r.get("Locus", ""))
        mhc = r.get("MHC", "")
        if fc is None or locus is None:
            continue
        allele_short = mhc.replace("HLA-", "").replace("HLA*", "")
        hla_locus = next((c for c in mhc if c in ("A", "B", "C")), "")
        enriched.append({**r, "log2_foldchange_BR": fc, "Locus": locus,
                         "allele_short": allele_short, "HLA_locus": hla_locus})

    if not enriched:
        raise ValueError("No valid per-allele data found.")

    grouped = _mutations_by_frame(enriched, frame_order, ["log2_foldchange_BR"])
    active = [f for f in frame_order if grouped.get(f)]
    if not active:
        raise ValueError("No per-allele data found for any frame.")

    all_fc = [r["log2_foldchange_BR"] for r in enriched]
    y_lim = max(abs(v) for v in all_fc) * 1.12 if all_fc else 3.0

    rng = np.random.default_rng(42)

    glw = _PassthroughGLW()
    glw.setBackground("#f5f5f7")

    first_plot = None
    for col_i, frame in enumerate(active):
        mutations_data = grouped[frame]
        abbrev = FRAME_ABBREVS.get(frame, frame)

        p = glw.addPlot(row=0, col=col_i,
                        axisItems={"bottom": _RotatedAxisItem()})
        _style_plot(p, title=abbrev)
        p.getViewBox().setMouseEnabled(x=False, y=False)

        n_muts = len({r.get("Mutation", "") for r in mutations_data})
        glw.ci.layout.setColumnStretchFactor(col_i, max(1, n_muts))

        # Group by mutation position
        mut_order: List[str] = []
        mut_positions: Dict[str, int] = {}
        for r in mutations_data:
            m = r.get("Mutation", "")
            if m not in mut_positions:
                mut_positions[m] = len(mut_order)
                mut_order.append(m)

        # Scatter points
        xs, ys, brushes = [], [], []
        for r in mutations_data:
            m = r.get("Mutation", "")
            fc = r["log2_foldchange_BR"]
            hla_locus = r["HLA_locus"]
            x_pos = mut_positions[m] + float(rng.uniform(-0.15, 0.15))
            color = LOCUS_COLORS.get(hla_locus, _LOCUS_FALLBACK)
            xs.append(x_pos)
            ys.append(fc)
            brushes.append(pg.mkBrush(*_hex_to_rgba(color, 200)))

        scatter = pg.ScatterPlotItem(
            x=np.array(xs), y=np.array(ys),
            size=7, pen=pg.mkPen(None), brush=brushes,
        )
        p.addItem(scatter)

        # y = 0 line
        p.addItem(pg.InfiniteLine(
            pos=0, angle=0,
            pen=pg.mkPen("#555555", width=1, style=Qt.PenStyle.DashLine),
        ))

        ticks = [[(i, mut_order[i]) for i in range(len(mut_order))]]
        p.getAxis("bottom").setTicks(ticks)

        p.setXRange(-0.65, len(mut_order) - 0.35, padding=0)
        p.setYRange(-y_lim, y_lim, padding=0)

        if col_i == 0:
            p.setLabel("left", "log₂ FC Best Rank", color="#333333", size="8pt")
            first_plot = p
        else:
            p.hideAxis("left")
            if first_plot:
                p.setYLink(first_plot)

    return glw


def _build_per_allele_box(pa_rows: List[Dict]) -> "pg.GraphicsLayoutWidget":
    """
    Boxplot: distribution of per-allele log₂ FC, grouped by allele,
    ordered by HLA locus then allele name.  Jittered individual points shown.
    """
    enriched = []
    for r in pa_rows:
        fc = _safe_float(r.get("log2_foldchange_BR", ""))
        mhc = r.get("MHC", "")
        if fc is None:
            continue
        allele_short = mhc.replace("HLA-", "").replace("HLA*", "")
        hla_locus = next((c for c in mhc if c in ("A", "B", "C")), "")
        enriched.append({"fc": fc, "allele_short": allele_short,
                         "HLA_locus": hla_locus})

    if not enriched:
        raise ValueError("No valid per-allele data found.")

    allele_meta: Dict[str, str] = {}
    for r in enriched:
        if r["allele_short"] not in allele_meta:
            allele_meta[r["allele_short"]] = r["HLA_locus"]
    allele_order = sorted(allele_meta, key=lambda a: (allele_meta[a], a))

    by_allele: Dict[str, List[float]] = {a: [] for a in allele_order}
    for r in enriched:
        by_allele[r["allele_short"]].append(r["fc"])

    all_fc = [r["fc"] for r in enriched]
    y_lim = max(abs(v) for v in all_fc) * 1.12 if all_fc else 3.0

    glw = _PassthroughGLW()
    glw.setBackground("#f5f5f7")
    p = glw.addPlot(row=0, col=0, axisItems={"bottom": _RotatedAxisItem()})
    _style_plot(p)
    p.getViewBox().setMouseEnabled(x=False, y=False)
    p.setLabel("left", "log₂ FC Best Rank", color="#333333", size="8pt")

    rng = np.random.default_rng(42)
    box_pen = pg.mkPen("#444444", width=0.9)
    median_pen = pg.mkPen("#222222", width=1.8)

    for i, allele in enumerate(allele_order):
        vals = np.array(by_allele[allele])
        if len(vals) == 0:
            continue
        locus = allele_meta[allele]
        color = LOCUS_COLORS.get(locus, _LOCUS_FALLBACK)

        q1, q2, q3 = float(np.percentile(vals, 25)), float(np.percentile(vals, 50)), float(np.percentile(vals, 75))
        iqr = q3 - q1
        w_lo = max(float(vals.min()), q1 - 1.5 * iqr)
        w_hi = min(float(vals.max()), q3 + 1.5 * iqr)

        # Box (Q1 → Q3)
        box = pg.BarGraphItem(
            x=np.array([i]), y0=np.array([q1]), y1=np.array([q3]),
            width=0.45,
            brush=pg.mkBrush(200, 200, 200, 160),
            pen=box_pen,
        )
        p.addItem(box)

        # Median line
        p.plot([i - 0.225, i + 0.225], [q2, q2], pen=median_pen)

        # Whiskers + caps
        for y_end, y_base in [(w_hi, q3), (w_lo, q1)]:
            p.plot([i, i], [y_base, y_end], pen=box_pen)
            p.plot([i - 0.1, i + 0.1], [y_end, y_end], pen=box_pen)

        # Jittered individual points
        jitter = rng.uniform(-0.14, 0.14, size=len(vals))
        scatter = pg.ScatterPlotItem(
            x=np.full(len(vals), i) + jitter,
            y=vals,
            size=5,
            pen=pg.mkPen(None),
            brush=pg.mkBrush(*_hex_to_rgba(color, 190)),
        )
        p.addItem(scatter)

    # y = 0 line
    p.addItem(pg.InfiniteLine(
        pos=0, angle=0,
        pen=pg.mkPen("#555555", width=1, style=Qt.PenStyle.DashLine),
    ))

    ticks = [[(i, allele_order[i]) for i in range(len(allele_order))]]
    p.getAxis("bottom").setTicks(ticks)

    p.setXRange(-0.65, len(allele_order) - 0.35, padding=0)
    p.setYRange(-y_lim, y_lim, padding=0)

    return glw


# ── Font helper ───────────────────────────────────────────────────────────────

def _small_font():
    from PyQt6.QtGui import QFont
    f = QFont()
    f.setPointSize(7)
    return f


# ── Rotated bottom-axis labels ────────────────────────────────────────────────

class _RotatedAxisItem(pg.AxisItem):
    """
    Bottom AxisItem that draws tick labels rotated -90° so they don't overlap.

    pyqtgraph's setStyle() does not accept 'tickTextAngle', so we subclass
    and override drawPicture() to rotate each label ourselves.
    """
    _LABEL_PX = 68  # vertical space reserved for rotated labels

    def __init__(self, **kwargs):
        kwargs["orientation"] = "bottom"
        super().__init__(**kwargs)
        self.setStyle(tickTextHeight=self._LABEL_PX, tickTextOffset=2,
                      tickFont=_small_font())

    def drawPicture(self, p, axisSpec, tickSpecs, textSpecs):
        from PyQt6.QtCore import QRectF
        p.setRenderHint(p.RenderHint.Antialiasing, False)
        p.setRenderHint(p.RenderHint.TextAntialiasing, True)

        # Axis line
        pen, pt1, pt2 = axisSpec
        p.setPen(pen)
        p.drawLine(pt1, pt2)

        # Tick marks
        for pen, pt1, pt2 in tickSpecs:
            p.setPen(pen)
            p.drawLine(pt1, pt2)

        # Rotated labels — pivot around the top-centre of each text rect
        if self.style.get("showValues", True):
            font = self.style.get("tickFont") or self.font()
            p.setFont(font)
            p.setPen(self.textPen())
            for rect, _flags, text in textSpecs:
                p.save()
                p.translate(rect.center().x(), rect.top())
                p.rotate(-90)
                # After -90° rotation the label extends leftward;
                # right-align so the end of the text sits near the tick.
                p.drawText(
                    QRectF(-(self._LABEL_PX - 2), -rect.width() / 2,
                           self._LABEL_PX - 2, rect.width()),
                    Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter,
                    text,
                )
                p.restore()


# ── Scroll + zoom container ───────────────────────────────────────────────────

class _PassthroughGLW(pg.GraphicsLayoutWidget):
    """
    GraphicsLayoutWidget that ignores wheel events so they propagate up
    to the parent _ZoomScrollArea instead of being swallowed by pyqtgraph.
    All ViewBoxes inside should also have setMouseEnabled(False).
    """
    def wheelEvent(self, ev):
        ev.ignore()   # bubble up to the QScrollArea parent


class _ZoomScrollArea(QScrollArea):
    """
    QScrollArea holding a _PassthroughGLW.

    • Two-finger trackpad scroll  → pan  (QScrollArea native behaviour)
    • Ctrl + scroll               → zoom (resize the inner GLW)
    • +/− buttons in _make_panel → zoom via zoom_in() / zoom_out()
    """
    _MIN_ZOOM = 0.25
    _MAX_ZOOM = 5.0

    def __init__(self, glw: "_PassthroughGLW", base_height: int = 420):
        super().__init__()
        self._glw = glw
        self._base_height = base_height
        self._zoom = 1.0
        self._initialised = False
        self.setWidget(glw)
        self.setWidgetResizable(False)
        self.setFrameShape(QFrame.Shape.NoFrame)

    # ── size management ──────────────────────────────────────────
    def _apply_zoom(self, z: float) -> None:
        self._zoom = max(self._MIN_ZOOM, min(self._MAX_ZOOM, z))
        vp_w = max(200, self.viewport().width())
        new_w = max(vp_w, int(vp_w * self._zoom))
        new_h = max(self._base_height, int(self._base_height * self._zoom))
        self._glw.setFixedSize(new_w, new_h)

    def zoom_level(self) -> float:
        return self._zoom

    def zoom_in(self)    -> None: self._apply_zoom(self._zoom * 1.3)
    def zoom_out(self)   -> None: self._apply_zoom(self._zoom / 1.3)
    def zoom_reset(self) -> None: self._apply_zoom(1.0)

    # ── Qt events ───────────────────────────────────────────────
    def showEvent(self, ev):
        super().showEvent(ev)
        if not self._initialised:
            self._initialised = True
            self._apply_zoom(1.0)

    def resizeEvent(self, ev):
        super().resizeEvent(ev)
        if self._initialised:
            self._apply_zoom(self._zoom)   # keep ratio on window resize

    def wheelEvent(self, ev):
        if ev.modifiers() & Qt.KeyboardModifier.ControlModifier:
            factor = 1.13 if ev.angleDelta().y() > 0 else 1.0 / 1.13
            self._apply_zoom(self._zoom * factor)
            ev.accept()
        else:
            super().wheelEvent(ev)         # default: scroll the area


# ── Panel wrapper (plot widget + save button) ─────────────────────────────────

def _make_panel(
    pg_widget: "_PassthroughGLW",
    save_stem: str,
    base_height: int = 420,
) -> QWidget:
    """
    Wrap a _PassthroughGLW in a scroll area with zoom controls and a save button.

    Interaction:
      Scroll (two-finger trackpad)  → pan the plot
      Ctrl + scroll                 → zoom in / out
      +  /  −  buttons              → zoom in / out
    """
    panel = QWidget()
    vbox = QVBoxLayout(panel)
    vbox.setContentsMargins(0, 0, 0, 0)
    vbox.setSpacing(0)

    # ── Controls row ─────────────────────────────────────────────
    ctrl = QHBoxLayout()
    ctrl.setContentsMargins(8, 4, 8, 4)

    hint_lbl = QLabel("Scroll to pan · Ctrl+scroll or ＋/－ to zoom")
    hint_lbl.setObjectName("lbl_info")
    hint_lbl.setStyleSheet("font-size: 11px; color: #8e8e93;")

    zoom_out_btn = QPushButton("－")
    zoom_out_btn.setFixedSize(26, 26)
    zoom_out_btn.setObjectName("btn_outline")

    zoom_lbl = QLabel("100%")
    zoom_lbl.setFixedWidth(42)
    zoom_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
    zoom_lbl.setStyleSheet("font-size: 11px; color: #555;")

    zoom_in_btn = QPushButton("＋")
    zoom_in_btn.setFixedSize(26, 26)
    zoom_in_btn.setObjectName("btn_outline")

    save_btn = QPushButton("Save plot…")
    save_btn.setObjectName("btn_outline")
    save_btn.setFixedHeight(26)

    ctrl.addWidget(hint_lbl)
    ctrl.addStretch()
    ctrl.addWidget(zoom_out_btn)
    ctrl.addWidget(zoom_lbl)
    ctrl.addWidget(zoom_in_btn)
    ctrl.addSpacing(10)
    ctrl.addWidget(save_btn)

    vbox.addLayout(ctrl)

    # ── Scroll area ───────────────────────────────────────────────
    scroll = _ZoomScrollArea(pg_widget, base_height=base_height)
    vbox.addWidget(scroll, 1)

    # Wire zoom buttons; update label on every zoom change
    def _apply(z: float) -> None:
        scroll._apply_zoom(z)
        zoom_lbl.setText(f"{int(round(scroll.zoom_level() * 100))}%")

    zoom_in_btn.clicked.connect(lambda: _apply(scroll.zoom_level() * 1.3))
    zoom_out_btn.clicked.connect(lambda: _apply(scroll.zoom_level() / 1.3))

    # Patch _ZoomScrollArea to also update the label on Ctrl+scroll
    _orig_apply = scroll._apply_zoom

    def _patched_apply(z: float) -> None:
        _orig_apply(z)
        zoom_lbl.setText(f"{int(round(scroll.zoom_level() * 100))}%")

    scroll._apply_zoom = _patched_apply

    # ── Save ──────────────────────────────────────────────────────
    def _save():
        default = Path.home() / f"{save_stem}.png"
        dest, _ = QFileDialog.getSaveFileName(
            panel, "Save plot", str(default),
            "PNG image (*.png);;SVG vector (*.svg)",
        )
        if not dest:
            return
        try:
            if dest.lower().endswith(".svg"):
                exp = pg.exporters.SVGExporter(pg_widget.scene())
            else:
                exp = pg.exporters.ImageExporter(pg_widget.scene())
                exp.parameters()["width"] = 2400
            exp.export(dest)
        except Exception as exc:
            QMessageBox.critical(panel, "Save failed", str(exc))

    save_btn.clicked.connect(_save)
    return panel


def _make_error_panel(msg: str) -> QWidget:
    w = QWidget()
    lbl = QLabel(msg)
    lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
    lbl.setWordWrap(True)
    lbl.setStyleSheet("color: #ff3b30; font-size: 12px; padding: 24px;")
    vbox = QVBoxLayout(w)
    vbox.addWidget(lbl)
    return w


def _make_install_panel() -> QWidget:
    w = QWidget()
    lbl = QLabel(
        "pyqtgraph is not installed.\n\n"
        "Go to the Setup step and click\n"
        '"Install Python plotting library" to enable plots.'
    )
    lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
    lbl.setWordWrap(True)
    lbl.setStyleSheet("color: #636366; font-size: 13px; padding: 32px;")
    vbox = QVBoxLayout(w)
    vbox.addWidget(lbl)
    return w


# ── Main PlotViewer widget ────────────────────────────────────────────────────

class PlotViewer(QWidget):
    """
    Tabbed interactive plot viewer for CD8scape output files.

    Call load() after each run.  Tabs are added/removed dynamically:
      • "Escape Scores"             — always present after a run
      • "Percentile Scores"         — only when has_percentile=True
      • "Per-Allele Escape Scores"  — only when has_per_allele=True
    """

    def __init__(self, parent: Optional[QWidget] = None):
        super().__init__(parent)
        vbox = QVBoxLayout(self)
        vbox.setContentsMargins(0, 0, 0, 0)
        vbox.setSpacing(0)

        if not _PG_OK:
            vbox.addWidget(_make_install_panel())
            return

        self._tabs = QTabWidget()
        self._tabs.setDocumentMode(True)
        vbox.addWidget(self._tabs)

    # ── Public API ────────────────────────────────────────────────────────────

    def load(
        self,
        folder: Optional[Path],
        folder_name: str = "",
        run_suffix: str = "",
        sim_suffix: str = "simulated",
        has_per_allele: bool = False,
        has_percentile: bool = False,
    ) -> None:
        # Re-attempt import in case pyqtgraph was just installed this session
        _ensure_pg()

        if not _PG_OK:
            self._tabs.clear()
            self._tabs.addTab(_make_install_panel(), "Escape Scores")
            return

        self._tabs.clear()

        if folder is None or not folder.is_dir():
            self._tabs.addTab(_make_error_panel("No output folder available."),
                              "Escape Scores")
            return

        fname = folder_name or folder.name

        # Standardised filename stem helpers
        suffix_part = f"_{run_suffix}" if run_suffix else ""
        sim_part    = f"_{sim_suffix}" if sim_suffix else "_simulated"

        def stem(plot_type: str, s: str = suffix_part) -> str:
            return f"CD8scape_{plot_type}_{fname}{s}"

        # ── Load observed HMBR CSV ────────────────────────────────────────
        hmbr_path = _find_file(folder, "harmonic_mean_best_ranks", run_suffix)
        if not (hmbr_path and hmbr_path.exists()):
            self._tabs.addTab(
                _make_error_panel(
                    "harmonic_mean_best_ranks.csv not found.\n"
                    "Run an analysis to generate output files."
                ),
                "Escape Scores",
            )
            return

        try:
            obs_rows = _load_csv(hmbr_path)
        except Exception as exc:
            self._tabs.addTab(
                _make_error_panel(f"Could not read {hmbr_path.name}:\n{exc}"),
                "Escape Scores",
            )
            return

        frame_order, frame_colors, mutation_colors = \
            _frame_order_and_colors(obs_rows)

        # ── Tab 1 — Escape Scores ─────────────────────────────────────────
        try:
            glw = _build_escape_scores(
                obs_rows, frame_order, frame_colors, mutation_colors)
            tab = _make_panel(glw, stem("escape_scores"))
        except Exception as exc:
            tab = _make_error_panel(str(exc))
        self._tabs.addTab(tab, "Escape Scores")

        # ── Tab 2 — Percentile Scores (conditional) ───────────────────────
        if has_percentile:
            sim_path = _find_file(folder, "harmonic_mean_best_ranks", sim_suffix)
            if sim_path and sim_path.exists():
                try:
                    sim_rows = _load_csv(sim_path)
                    glw = _build_percentile_scores(
                        obs_rows, sim_rows, frame_colors, mutation_colors)
                    tab = _make_panel(glw, stem("percentile_scores", sim_part),
                                      base_height=500)
                except Exception as exc:
                    tab = _make_error_panel(str(exc))
            else:
                tab = _make_error_panel(
                    f"Simulated output not found ({sim_suffix}).\n"
                    "Re-run with percentile analysis enabled."
                )
            self._tabs.addTab(tab, "Percentile Scores")

        # ── Tab 3 — Per-Allele Escape Scores (conditional) ────────────────
        if has_per_allele:
            pa_path = _find_file(folder, "per_allele_best_ranks", run_suffix)
            if pa_path and pa_path.exists():
                try:
                    pa_rows = _load_csv(pa_path)
                    # Sub-tab widget
                    pa_tabs = QTabWidget()
                    pa_tabs.setDocumentMode(True)

                    # Per-mutation scatter
                    try:
                        glw_scat = _build_per_allele_scatter(
                            pa_rows, frame_order, frame_colors, mutation_colors
                        )
                        pa_tabs.addTab(
                            _make_panel(glw_scat,
                                        stem("per_allele_escape_scores_per_mutation")),
                            "Per Mutation",
                        )
                    except Exception as exc:
                        pa_tabs.addTab(_make_error_panel(str(exc)), "Per Mutation")

                    # By-allele boxplot
                    try:
                        glw_box = _build_per_allele_box(pa_rows)
                        pa_tabs.addTab(
                            _make_panel(glw_box,
                                        stem("per_allele_escape_scores_by_allele")),
                            "By Allele",
                        )
                    except Exception as exc:
                        pa_tabs.addTab(_make_error_panel(str(exc)), "By Allele")

                    tab = pa_tabs
                except Exception as exc:
                    tab = _make_error_panel(str(exc))
            else:
                tab = _make_error_panel(
                    "per_allele_best_ranks.csv not found.\n"
                    "Re-run with --per-allele enabled."
                )
            self._tabs.addTab(tab, "Per-Allele Escape Scores")

    def reset(self) -> None:
        if _PG_OK:
            self._tabs.clear()
