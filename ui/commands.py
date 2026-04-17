"""
commands.py — CLI manifest and per-command Streamlit renderers.

This is a *manifest* of CD8scape's CLI surface — nothing more. Each
`render_*_options` function:

    * Displays Streamlit controls that map 1:1 onto CD8scape CLI flags.
    * Returns the list of CLI arguments that should be appended after
      the subcommand name (and folder path, where applicable).

No function in this module knows what CD8scape does with the flags — it
only knows which flags exist, their type, their default, and the
human-readable description taken from CD8scape's own --help. If CD8scape
adds or removes a flag, update this file; nothing else should need to
change.

Help text is paraphrased from CD8scape's README / --help so that users
see the same description Anthropic-agnostically in the UI tooltip as
they would on the command line.
"""

from __future__ import annotations

from typing import List

import streamlit as st


# Commands CD8scape understands, in the same order as in ./CD8scape.jl --help.
COMMANDS: List[str] = [
    "prep",
    "read",
    "simulate",
    "run",
    "run_supertype",
    "percentile",
]


COMMAND_HELP = {
    "prep":          "Set up the environment (install Julia dependencies, "
                     "validate the netMHCpan path and Perl).",
    "read":          "Parse variant + reading-frame input files from the "
                     "folder, producing variants.csv and frames.csv.",
    "simulate":      "Read reading frames, then generate simulated "
                     "single-nucleotide variants per frame.",
    "run":           "Run the peptide-generation + netMHCpan pipeline for an "
                     "individual HLA genotype.",
    "run_supertype": "Run the pipeline for a representative supertype HLA "
                     "panel.",
    "percentile":    "Compute observed log2 fold-change percentiles against "
                     "the simulated distribution.",
}


# -- shared renderers ---------------------------------------------------------


def _suffix_and_latest(key_prefix: str, default_suffix: str = "") -> List[str]:
    """Render the shared --suffix / --latest / --no-latest controls."""
    args: List[str] = []

    suffix = st.text_input(
        "Output suffix (`--suffix`)",
        value=default_suffix,
        key=f"{key_prefix}_suffix",
        placeholder="leave blank for default",
        help=(
            "Append _<name> before the extension of all output files "
            "(e.g. 'run1' produces variants_run1.csv). Leave blank to use "
            "CD8scape's default behaviour."
        ),
    )
    if suffix.strip():
        args += ["--suffix", suffix.strip()]

    latest_choice = st.radio(
        "When multiple input candidates exist",
        options=["latest", "no-latest"],
        index=0,
        format_func=lambda x: {
            "latest":    "Pick the most recent (--latest, CD8scape default)",
            "no-latest": "Error on ambiguity (--no-latest)",
        }[x],
        key=f"{key_prefix}_latest",
        horizontal=True,
        help=(
            "When reading an input file without an explicit suffix and "
            "multiple candidates exist (e.g. frames.csv and "
            "frames_simulated.csv), --latest picks the most recently "
            "modified file; --no-latest raises an error instead."
        ),
    )
    if latest_choice == "no-latest":
        args.append("--no-latest")
    # --latest is CD8scape's default, no need to emit it explicitly.
    return args


# -- per-command renderers ----------------------------------------------------


def render_prep_options() -> List[str]:
    st.caption(
        "`prep` installs Julia dependencies and validates the netMHCpan / "
        "Perl setup. There are no options."
    )
    return []


def render_read_options() -> List[str]:
    args: List[str] = []
    if st.checkbox(
        "Pass `--aa` (variants are amino-acid format, .aa file)",
        value=False,
        key="read_aa",
        help=(
            "Parse amino-acid-level variants from a .aa file instead of "
            "VCF or Samfire trajectories. See ./CD8scape.jl --help for the "
            "file format."
        ),
    ):
        args.append("--aa")
    args += _suffix_and_latest("read")
    return args


def render_simulate_options() -> List[str]:
    args: List[str] = []

    sampling = st.radio(
        "Sampling mode",
        options=["all", "number", "proportion"],
        index=1,  # CD8scape's help text calls out --n as the common path.
        format_func=lambda x: {
            "all":        "Write all variants (no --n / --p)",
            "number":     "Sample a fixed number (--n)",
            "proportion": "Sample a proportion (--p)",
        }[x],
        key="sim_sampling",
        help=(
            "CD8scape: if both --n and --p are provided, --n takes "
            "precedence; if neither is provided, all variants are written."
        ),
    )
    if sampling == "number":
        n = st.number_input(
            "Number of variants (`--n`)",
            min_value=1,
            value=1000,
            step=1,
            key="sim_n",
            help="CD8scape default: 1000 when --n is given without a value.",
        )
        args += ["--n", str(int(n))]
    elif sampling == "proportion":
        p = st.slider(
            "Proportion to sample (`--p`)",
            min_value=0.001,
            max_value=0.999,
            value=0.10,
            step=0.001,
            key="sim_p",
            help="CD8scape default: 0.1 when --p is given without a value.",
        )
        args += ["--p", f"{p:g}"]

    seed = st.number_input(
        "Random seed (`--seed`)",
        min_value=0,
        value=1320,
        step=1,
        key="sim_seed",
        help="RNG seed. CD8scape default: 1320.",
    )
    args += ["--seed", str(int(seed))]

    # CD8scape defaults --suffix to 'simulated' for simulate when omitted,
    # so pre-fill the field to match. Users can override or clear it.
    args += _suffix_and_latest("sim", default_suffix="simulated")
    return args


def _render_run_like_options(key_prefix: str) -> List[str]:
    args: List[str] = []

    use_max = st.checkbox(
        "Use maximum parallel chunks (`--t max`)",
        value=False,
        key=f"{key_prefix}_tmax",
        help=(
            "Run NetMHCpan chunks in parallel up to CD8SCAPE_MAX_THREADS "
            "(default 8). Uncheck to pick a specific number."
        ),
    )
    if use_max:
        args += ["--t", "max"]
    else:
        t = st.number_input(
            "Parallel chunks (`--t`)",
            min_value=1,
            value=1,
            step=1,
            key=f"{key_prefix}_t",
            help=(
                "Max number of chunks to execute in parallel when running "
                "NetMHCpan. CD8scape default: 1."
            ),
        )
        args += ["--t", str(int(t))]

    if st.checkbox(
        "Pass `--per-allele`",
        value=False,
        key=f"{key_prefix}_per_allele",
        help=(
            "Compute log2 fold change per allele in the provided genome; "
            "writes a separate per_allele_best_ranks.csv. Only alleles with "
            "ancestral EL rank ≤ 2% are included."
        ),
    ):
        args.append("--per-allele")

    if st.checkbox(
        "Pass `--verbose`",
        value=False,
        key=f"{key_prefix}_verbose",
        help="Preserve per-allele logs and temp files for debugging.",
    ):
        args.append("--verbose")

    args += _suffix_and_latest(key_prefix)
    return args


def render_run_options() -> List[str]:
    return _render_run_like_options("run")


def render_run_supertype_options() -> List[str]:
    st.caption(
        "Note: supertype-panel alleles are population-frequency surrogates, "
        "not an individual's genotype. CD8scape prints a reminder when "
        "`--per-allele` is used with this command."
    )
    return _render_run_like_options("run_supertype")


def render_percentile_options() -> List[str]:
    args: List[str] = []

    if st.checkbox(
        "Pass `--per-allele`",
        value=False,
        key="pct_per_allele",
        help=(
            "Operate on per-allele fold changes instead of harmonic-mean "
            "best ranks. Percentiles are then computed per allele."
        ),
    ):
        args.append("--per-allele")

    sim_file = st.text_input(
        "Simulated file path (`--s`)",
        value="",
        key="pct_s",
        placeholder="leave blank for CD8scape default",
        help=(
            "Path to the simulated file. CD8scape default: "
            "harmonic_mean_best_ranks_simulated.csv (or "
            "per_allele_best_ranks_simulated.csv when --per-allele is set)."
        ),
    )
    if sim_file.strip():
        args += ["--s", sim_file.strip()]

    obs_file = st.text_input(
        "Observed file path (`--o`)",
        value="",
        key="pct_o",
        placeholder="leave blank for CD8scape default",
        help=(
            "Path to the observed file. CD8scape default: the most recent "
            "matching file in the data folder, excluding _simulated files."
        ),
    )
    if obs_file.strip():
        args += ["--o", obs_file.strip()]

    return args


# -- dispatch -----------------------------------------------------------------


def render_options(command: str) -> List[str]:
    """Render the option controls for `command` and return the CLI argv tail.

    The returned list does NOT include the subcommand name or the folder
    path — only the option arguments.
    """
    if command == "prep":
        return render_prep_options()
    if command == "read":
        return render_read_options()
    if command == "simulate":
        return render_simulate_options()
    if command == "run":
        return render_run_options()
    if command == "run_supertype":
        return render_run_supertype_options()
    if command == "percentile":
        return render_percentile_options()
    raise ValueError(f"Unknown CD8scape command: {command!r}")


def needs_folder(command: str) -> bool:
    """All CD8scape commands except `prep` require a folder path."""
    return command != "prep"
