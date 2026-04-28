"""
workflow.py — Workflow abstraction layer for CD8scape UI.

Provides a higher-level view of CD8scape that maps user-facing concepts
("What kind of data do I have?", "What analysis do I want?") onto the
underlying CLI commands.

Users should never need to know that `read --aa` and `read` are the same
subcommand, or that `run_supertype` exists as a separate entry point. They
choose based on their data and intent; this module resolves the mapping.

CLI mapping
-----------
Prepare step (Step 1):
    "Nucleotide input"     → read
    "Amino acid input"     → read --aa
    "Simulate variants"    → simulate

Run step (Step 2):
    "Standard analysis"    → run
    "Supertype analysis"   → run_supertype

Execution
---------
The two steps always run in sequence: prepare → run. Both are assembled as
standard CD8scape argv tails by the caller (app.py). This module only holds
the declarative mapping — no subprocess logic lives here.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List


@dataclass(frozen=True)
class WorkflowChoice:
    """One selectable option in a workflow step.

    Attributes
    ----------
    key:
        Unique identifier used in st.session_state to remember the selection.
    label:
        Short label shown in the radio button.
    detail:
        One-sentence description shown below the radio (as a caption).
    cli_command:
        CD8scape subcommand to invoke for this choice.
    cli_fixed_args:
        Extra CLI arguments that are always appended for this choice,
        *before* any user-configurable option args.
        Example: ``["--aa"]`` for the amino-acid read path.
    """

    key: str
    label: str
    detail: str
    cli_command: str
    cli_fixed_args: List[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Step 1 — Prepare data
# ---------------------------------------------------------------------------

PREP_CHOICES: List[WorkflowChoice] = [
    WorkflowChoice(
        key="nucleotide",
        label="Nucleotide input",
        detail=(
            "Parse nucleotide-level variants into a variants table. "
            "Needs: a variant file (.vcf / .vcf.gz or single_locus_trajectories…out) "
            "and reading-frame files (sequences.fasta + Consensus.fa for NCBI data, "
            "or Reading_Frames.dat for Samfire data)."
        ),
        cli_command="read",
        cli_fixed_args=[],
    ),
    WorkflowChoice(
        key="amino_acid",
        label="Amino acid input (.aa file)",
        detail=(
            "Read pre-called amino-acid-level variants directly from a .aa file. "
            "Needs: a .aa variant file and reading-frame files "
            "(sequences.fasta + Consensus.fa, or Reading_Frames.dat)."
        ),
        cli_command="read",
        cli_fixed_args=["--aa"],
    ),
    WorkflowChoice(
        key="simulation",
        label="Simulate variants",
        detail=(
            "Generate synthetic single-nucleotide variants for benchmarking. "
            "Needs: reading-frame files only "
            "(sequences.fasta + Consensus.fa, or Reading_Frames.dat). "
            "No real variant file required."
        ),
        cli_command="simulate",
        cli_fixed_args=[],
    ),
]

# ---------------------------------------------------------------------------
# Step 2 — Run analysis
# ---------------------------------------------------------------------------

RUN_CHOICES: List[WorkflowChoice] = [
    WorkflowChoice(
        key="standard",
        label="Standard analysis",
        detail=(
            "Run netMHCpan for an individual HLA genotype. "
            "Reads alleles from alleles.txt (one HLA allele per line, "
            "e.g. HLA-A03:01). Produces harmonic-mean best-rank scores."
        ),
        cli_command="run",
        cli_fixed_args=[],
    ),
    WorkflowChoice(
        key="supertype",
        label="Supertype panel analysis",
        detail=(
            "Run for a representative population supertype HLA panel. "
            "Reads alleles from supertype_panel.csv (columns: Allele, Frequency) "
            "instead of alleles.txt. Use this for population-level analyses."
        ),
        cli_command="run_supertype",
        cli_fixed_args=[],
    ),
]

# ---------------------------------------------------------------------------
# Lookup helpers
# ---------------------------------------------------------------------------


def choice_by_key(choices: List[WorkflowChoice], key: str) -> WorkflowChoice:
    """Return the WorkflowChoice whose key matches *key*.

    Raises ``KeyError`` if no match is found — this is always a programming
    error (the key must come from the same list), so a loud failure is right.
    """
    for c in choices:
        if c.key == key:
            return c
    raise KeyError(
        f"No WorkflowChoice with key {key!r} in the provided list. "
        f"Available keys: {[c.key for c in choices]}"
    )
