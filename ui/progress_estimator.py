"""
progress_estimator.py — Live overall-progress estimator for CD8scape pipelines.

Architecture
------------
Each CD8scape run step (run / run_supertype) contains two distinct phases:

  Phase 1 — netMHCpan:
    run_netMHCpan.jl emits within-step progress in two formats:

      Sequential (--t 1):
          "Chunk N/M, batch B/B2 (PPP peptides, AA alleles). PCT% complete."
      Parallel (--t > 1):
          "N/M units complete. PCT%"

    Both formats are emitted with \\r overwrite; Python's universal-newlines
    mode treats \\r as a line terminator so each update arrives as a separate
    yielded string.

  Phase 2 — post-netMHCpan processing:
    Three Julia scripts run sequentially after netMHCpan finishes:

      process_scores.jl     — joins netMHCpan TSV with peptide labels.
                              First output: "Sorting by Locus…"
      process_best_ranks.jl — streams processed_output.csv, computing
                              per-(Locus, MHC, Mutation) best ranks.
                              Progress marker every 10 M rows:
                                "  NM rows scanned (…)"
                              then materialises HMBR and fold changes.
      variant_fates.jl      — traces variants through filter stages.
                              Final line: "variant_fates.csv written …"

    For large datasets the post-processing phase can rival or exceed netMHCpan
    in wall-clock time.  Without accounting for it the progress bar freezes at
    ~99 % for many minutes.

Step-weight model
-----------------
The overall bar is a weighted sum of per-step fractions:

    overall = Σ frac[i] * weight[i]

Weights are calibrated to typical wall-clock proportions (1 thread):

    2-step:  read ≈ 1 %,   run ≈ 99 %
    5-step:  read ≈ 1 %,   run_obs ≈ 48 %,  simulate ≈ 2 %,
             run_sim ≈ 48 %,  percentile ≈ 2 %  (before sim-size adjustment)

Thread-count model
------------------
More threads speed up netMHCpan (the parallel phase) but not post-processing
or I/O-bound steps.  We use sqrt(n_threads) as an effective speedup, capped
at sqrt(32) ≈ 5.7×.  Higher thread counts deflate run-step weights (relative
to fixed-cost steps) and shrink the netMHCpan cap within each run step.

Intra-step split
----------------
Within a run step the fraction is split between the two phases:

    frac[step] = nmhc_progress * NMHC_CAP
                 + post_progress * (POST_CAP - NMHC_CAP)

where:
    NMHC_CAP  = 0.75  — netMHCpan saturates here (not 1.0)
    POST_CAP  = 0.99  — hard ceiling until complete_step() fires

NMHC_CAP is deflated at higher thread counts because fast netMHCpan means
post-processing is a larger fraction of actual wall time:

    effective_NMHC_CAP = NMHC_CAP * (1 - 0.12 * log2(n_threads))
    (minimum 0.55)

Post-processing milestone table
--------------------------------
Ordered milestones map known log-line patterns to positions (0.0–1.0) in the
post-processing window.  Fractions only ever advance (high-watermark).

    process_scores.jl phase:
        "Sorting by Locus"           → 0.08
        "Saving results to"          → 0.18
        "Processing completed"       → 0.28

    process_best_ranks.jl phase:
        "Reading input file:"        → 0.32
        "Streaming CSV to compute"   → 0.36
        "rows scanned"               → nudge (+0.03 each, cap 0.70)
        "Successfully processed"     → 0.72
        "Materialising.*ancestral"   → 0.77
        "Found best ranks for anc"   → 0.81
        "Materialising.*derived"     → 0.84
        "Found best ranks for der"   → 0.87
        "Saved best ranks to"        → 0.91
        "Calculating harmonic"       → 0.93
        "Saved harmonic mean"        → 0.97

    variant_fates.jl phase:
        "variant_fates.*written"     → 1.00
"""
from __future__ import annotations

import math
import os
import re
import sys
from typing import FrozenSet, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Regex patterns
# ---------------------------------------------------------------------------

# netMHCpan: sequential
_RE_SEQ = re.compile(
    r"Chunk\s+(\d+)/(\d+).*?(\d+)\s*%\s*complete",
    re.IGNORECASE,
)

# netMHCpan: parallel
_RE_PAR = re.compile(
    r"(\d+)/(\d+)\s+units\s+complete\.\s*(\d+)\s*%",
    re.IGNORECASE,
)

# post-processing: rows scanned in process_best_ranks
_RE_ROWS = re.compile(r"(\d+)M rows scanned", re.IGNORECASE)

# Post-processing milestone table: (pattern, post_fraction)
# Ordered by expected appearance in the log.  Each entry advances the
# post-processing fraction only if it's larger than the current value.
_POST_MILESTONES: List[Tuple[re.Pattern, float]] = [
    (re.compile(r"Sorting by Locus",                  re.I), 0.08),
    (re.compile(r"Saving results to",                 re.I), 0.18),
    (re.compile(r"Processing completed",              re.I), 0.28),
    (re.compile(r"Reading input file:",               re.I), 0.32),
    (re.compile(r"Streaming CSV to compute",          re.I), 0.36),
    (re.compile(r"Successfully processed",            re.I), 0.72),
    (re.compile(r"Materialising.*ancestral",          re.I), 0.77),
    (re.compile(r"Found best ranks for anc",          re.I), 0.81),
    (re.compile(r"Materialising.*derived",            re.I), 0.84),
    (re.compile(r"Found best ranks for der",          re.I), 0.87),
    (re.compile(r"Saved best ranks to",               re.I), 0.91),
    (re.compile(r"Calculating harmonic",              re.I), 0.93),
    (re.compile(r"Saved harmonic mean",               re.I), 0.97),
    (re.compile(r"variant_fates.*written",            re.I), 1.00),
]

# Fraction of step weight used by netMHCpan at 1 thread.  Post-processing
# gets the remaining (POST_CAP - _NMHC_CAP_BASE) = 0.24 of the step.
_NMHC_CAP_BASE: float = 0.75
_POST_CAP: float = 0.99          # hard ceiling until complete_step() fires
_ROWS_NUDGE: float = 0.03        # per "NM rows scanned" message
_ROWS_STREAMING_CAP: float = 0.70  # streaming can't push past here alone


# ---------------------------------------------------------------------------
# Step-weight tables
# ---------------------------------------------------------------------------

_RUN_STEPS_2: FrozenSet[int] = frozenset({1})
_RUN_STEPS_5: FrozenSet[int] = frozenset({1, 3})

_BASE_2 = [0.01, 0.99]

_BASE_5_RUN_OBS = 0.48
_BASE_5_READ    = 0.01
_BASE_5_PCT     = 0.02


def _make_5step_base(sim_size_factor: float) -> List[float]:
    sf = max(0.01, sim_size_factor)
    w_run_sim  = _BASE_5_RUN_OBS * sf
    w_simulate = min(0.01 + 0.01 * sf, 0.05)
    raw = [_BASE_5_READ, _BASE_5_RUN_OBS, w_simulate, w_run_sim, _BASE_5_PCT]
    total = sum(raw)
    return [w / total for w in raw]


_ENV_WEIGHTS = os.environ.get("CD8SCAPE_PROGRESS_WEIGHTS", "").strip()


def _parse_env_weights(n: int) -> Optional[List[float]]:
    if not _ENV_WEIGHTS:
        return None
    try:
        parts = [float(x.strip()) for x in _ENV_WEIGHTS.split(",")]
        if len(parts) != n:
            return None
        total = sum(parts)
        if total <= 0:
            return None
        return [p / total for p in parts]
    except ValueError:
        return None


def _adjusted_weights(
    base: List[float],
    run_indices: FrozenSet[int],
    n_threads: int,
) -> List[float]:
    speedup = min(n_threads, 32) ** 0.5
    adjusted = [
        w / speedup if i in run_indices else w
        for i, w in enumerate(base)
    ]
    total = sum(adjusted)
    return [w / total for w in adjusted] if total > 0 else list(base)


def _nmhc_cap(n_threads: int) -> float:
    """
    Effective netMHCpan cap within a run step, deflated at higher thread
    counts to reflect that faster netMHCpan means post-processing is a
    larger share of actual wall time.

        cap(1)  = 0.75
        cap(4)  ≈ 0.67
        cap(8)  ≈ 0.60
        cap(16) ≈ 0.53  (minimum 0.55 applied)
    """
    if n_threads <= 1:
        return _NMHC_CAP_BASE
    reduction = 0.12 * math.log2(max(1, n_threads))
    return max(0.55, _NMHC_CAP_BASE - reduction)


# ---------------------------------------------------------------------------
# ProgressEstimator
# ---------------------------------------------------------------------------

class ProgressEstimator:
    """
    Ingests CD8scape log lines and produces a continuously updated overall
    progress estimate (0–100%).

    Parameters
    ----------
    include_percentile : bool
        False → 2 steps [read, run]; True → 5 steps [read, run_obs, sim, run_sim, pct].
    n_threads : int
        The --t value passed to the run command.
    sim_size_factor : float
        Relative simulation size (1.0 = --n 1000 baseline).
        Ignored when include_percentile is False.
    """

    def __init__(
        self,
        include_percentile: bool = False,
        n_threads: int = 1,
        sim_size_factor: float = 1.0,
    ) -> None:
        self._include_pct = include_percentile
        self._n_threads = max(1, int(n_threads))
        self._sim_size_factor = max(0.01, float(sim_size_factor))
        self._n_steps = 5 if include_percentile else 2
        self._run_indices = _RUN_STEPS_5 if include_percentile else _RUN_STEPS_2
        self._nmhc_cap = _nmhc_cap(self._n_threads)

        if include_percentile:
            base = _make_5step_base(self._sim_size_factor)
        else:
            base = list(_BASE_2)

        env_w = _parse_env_weights(self._n_steps)
        if env_w is not None:
            self._weights: List[float] = env_w
        else:
            self._weights = _adjusted_weights(base, self._run_indices, self._n_threads)

        # Per-step overall fraction [0.0, 1.0]
        self._frac: List[float] = [0.0] * self._n_steps

        # Per-step post-processing tracking
        self._in_post: List[bool]  = [False] * self._n_steps  # entered post phase?
        self._post_frac: List[float] = [0.0] * self._n_steps  # post phase fraction [0,1]

        self._current: int = -1

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def start_step(self, step_index: int) -> None:
        if 0 <= step_index < self._n_steps:
            self._current = step_index

    def ingest_line(self, line: str) -> None:
        step = self._current
        if step < 0 or step >= self._n_steps:
            return
        if step not in self._run_indices:
            return

        cap = self._nmhc_cap

        # ── Phase 1: netMHCpan progress ───────────────────────────────────
        if not self._in_post[step]:
            pct = _parse_nmhc_pct(line)
            if pct is not None:
                # Map [0, 100] → [0, cap), never reaching the cap
                self._frac[step] = min(pct / 100.0 * cap, cap * 0.999)
                return

        # ── Phase 2: post-processing milestones ───────────────────────────
        # Check whether this line signals the start of post-processing
        # (any milestone match or a "rows scanned" hit transitions us)
        post_hit = _parse_post_progress(line, self._post_frac[step])
        if post_hit is not None:
            self._in_post[step] = True
            self._post_frac[step] = post_hit
            # Combine: netMHCpan part (at cap) + post-processing part
            post_window = _POST_CAP - cap
            self._frac[step] = min(
                cap + post_hit * post_window,
                _POST_CAP,
            )

    def complete_step(self, step_index: int, success: bool) -> None:
        if 0 <= step_index < self._n_steps:
            if success:
                self._frac[step_index] = 1.0

    # ------------------------------------------------------------------
    # Queries
    # ------------------------------------------------------------------

    @property
    def overall_fraction(self) -> float:
        return min(
            sum(f * w for f, w in zip(self._frac, self._weights)),
            1.0,
        )

    @property
    def overall_pct(self) -> int:
        return int(self.overall_fraction * 100)

    @property
    def n_steps(self) -> int:
        return self._n_steps

    @property
    def n_threads(self) -> int:
        return self._n_threads

    @property
    def include_percentile(self) -> bool:
        return self._include_pct

    def step_weights_pct(self) -> List[int]:
        return [round(w * 100) for w in self._weights]

    def step_fraction(self, step_index: int) -> float:
        if 0 <= step_index < self._n_steps:
            return self._frac[step_index]
        return 0.0

    @property
    def sim_size_factor(self) -> float:
        return self._sim_size_factor

    def summary(self) -> str:
        step_str = ", ".join(
            f"s{i}={'post' if self._in_post[i] else 'nmhc'}@{int(f*100)}%"
            for i, f in enumerate(self._frac)
        )
        return (
            f"ProgressEstimator(overall={self.overall_pct}%"
            f" threads={self._n_threads}"
            f" nmhc_cap={self._nmhc_cap:.0%}"
            f" sim_factor={self._sim_size_factor:.2f}"
            f" | {step_str})"
        )


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _parse_nmhc_pct(line: str) -> Optional[int]:
    """Extract a netMHCpan percentage from a log line; None if not found."""
    m = _RE_SEQ.search(line)
    if m:
        return int(m.group(3))
    m = _RE_PAR.search(line)
    if m:
        return int(m.group(3))
    return None


def _parse_post_progress(line: str, current: float) -> Optional[float]:
    """
    Check *line* against the post-processing milestone table and the
    row-count nudge pattern.

    Returns the new post-processing fraction (0.0–1.0) if this line
    advances progress, else None.  Only advances (high-watermark).
    """
    best = current

    # Named milestones (high-watermark)
    for pattern, frac in _POST_MILESTONES:
        if frac > best and pattern.search(line):
            best = frac
            break  # milestones are ordered; first match is the right one

    # "NM rows scanned" nudge — each occurrence adds _ROWS_NUDGE up to cap
    if _RE_ROWS.search(line):
        nudged = min(current + _ROWS_NUDGE, _ROWS_STREAMING_CAP)
        if nudged > best:
            best = nudged

    return best if best > current else None


def infer_threads_from_args(run_args: List[str], cpu_count: Optional[int] = None) -> int:
    """Extract the effective thread count from a CD8scape argv tail."""
    n_cpu = cpu_count or os.cpu_count() or 1
    cap = max(1, int(os.environ.get("CD8SCAPE_MAX_THREADS", "8")))

    for i, arg in enumerate(run_args):
        if arg in ("--t", "--thread") and i + 1 < len(run_args):
            val = run_args[i + 1].strip().lower()
            if val == "max":
                return min(n_cpu, cap)
            try:
                return max(1, min(int(val), cap))
            except ValueError:
                return 1
    return 1


# ---------------------------------------------------------------------------
# Standalone smoke-test / pipe mode
# ---------------------------------------------------------------------------

def _cli_main() -> None:  # pragma: no cover
    """
    Pipe CD8scape stdout through this script to watch overall progress.

        julia CD8scape.jl run /data --t 4 2>&1 | \\
            python ui/progress_estimator.py --steps 2 --threads 4
    """
    import argparse

    parser = argparse.ArgumentParser(
        description="Pipe CD8scape output to get live overall progress."
    )
    parser.add_argument("--steps",   type=int, choices=[2, 5], default=2)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--quiet",   action="store_true")
    args = parser.parse_args()

    include_pct = args.steps == 5
    est = ProgressEstimator(include_percentile=include_pct, n_threads=args.threads)

    print(
        f"[progress_estimator] {args.steps} steps, {args.threads} thread(s), "
        f"percentile={'yes' if include_pct else 'no'}"
        f" | nmhc_cap={est._nmhc_cap:.0%}"
        f" | weights={est.step_weights_pct()}",
        file=sys.stderr,
    )

    _STEP_MARKERS = [
        (re.compile(r"simulate.*variant",    re.I), 2),
        (re.compile(r"Simulated variant",    re.I), 2),
        (re.compile(r"percentile",           re.I), 4),
        (re.compile(r"Variants written to",  re.I), 0),
        (re.compile(r"Peptides\.pep.*written", re.I), 1),
        (re.compile(r"Using NetMHCpan",      re.I), 1),
    ]

    current_step = 0
    est.start_step(0)
    prev_pct = -1

    try:
        for raw_line in sys.stdin:
            line = raw_line.rstrip("\r\n")
            if not args.quiet:
                print(line)

            for pattern, step_idx in _STEP_MARKERS:
                if step_idx > current_step and pattern.search(line):
                    est.complete_step(current_step, success=True)
                    current_step = step_idx
                    est.start_step(current_step)
                    break

            est.ingest_line(line)

            pct = est.overall_pct
            if pct != prev_pct:
                print(
                    f"\r[{pct:3d}%] {est.summary()}",
                    end="", file=sys.stderr, flush=True,
                )
                prev_pct = pct

    except KeyboardInterrupt:
        pass

    print(file=sys.stderr)
    est.complete_step(current_step, success=True)
    print(f"[progress_estimator] final: {est.overall_pct}%", file=sys.stderr)


if __name__ == "__main__":
    _cli_main()
