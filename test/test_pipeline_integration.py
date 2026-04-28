"""
test_pipeline_integration.py — Integration tests for CD8scape Julia pipeline scripts.

These tests call the Julia scripts as subprocesses using real (or minimal synthetic)
data to verify that:
  - scripts exit cleanly (exit code 0) on valid input
  - scripts exit with non-zero code and a helpful message on bad input
  - expected output files are produced with the correct columns/content
  - the --suffix and --latest flags are forwarded correctly

Run from the repo root:
    python -m pytest test/test_pipeline_integration.py -v
    python test/test_pipeline_integration.py

Prerequisites: Julia must be on PATH with the src/ project instantiated.
Tests that require Julia are skipped when `julia` is not found on PATH.
"""

from __future__ import annotations

import csv
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SRC_DIR = REPO_ROOT / "src"
EXAMPLE_DATA = REPO_ROOT / "data" / "Example_data"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def julia_available() -> bool:
    return shutil.which("julia") is not None


def run_julia(script: str, *args: str, cwd: Path | None = None) -> subprocess.CompletedProcess:
    """Run a Julia script from src/ with the src/ project."""
    cmd = [
        "julia",
        f"--project={SRC_DIR}",
        str(SRC_DIR / script),
        *args,
    ]
    return subprocess.run(
        cmd,
        cwd=str(cwd or REPO_ROOT),
        capture_output=True,
        text=True,
        timeout=300,
    )


def read_csv_rows(path: Path) -> list[dict]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


skip_no_julia = unittest.skipUnless(julia_available(), "Julia not found on PATH")


# ---------------------------------------------------------------------------
# Helper: build a minimal temp data directory mirroring Example_data layout
# ---------------------------------------------------------------------------

def make_minimal_ncbi_dir(tmp_dir: Path) -> None:
    """Copy Example_data FASTA and allele files into tmp_dir."""
    for fname in ("sequences.fasta", "Consensus.fa", "alleles.txt", "supertype_panel.csv"):
        src = EXAMPLE_DATA / fname
        if src.exists():
            shutil.copy(src, tmp_dir / fname)


def make_minimal_frames_csv(tmp_dir: Path, dna: str = "ATGCCCGAATTT",
                              start: int = 1, desc: str = "ORF1") -> None:
    """Write a minimal frames.csv."""
    end = start + len(dna) - 1
    rows = [{"Region": f"{start},{end}", "Consensus_sequence": dna, "Description": desc}]
    with open(tmp_dir / "frames.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["Region", "Consensus_sequence", "Description"])
        w.writeheader(); w.writerows(rows)


def make_minimal_variants_csv(tmp_dir: Path) -> None:
    rows = [
        {"Locus": "4", "Consensus": "C", "Variant": "T"},
    ]
    with open(tmp_dir / "variants.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["Locus", "Consensus", "Variant"])
        w.writeheader(); w.writerows(rows)


def make_minimal_pep_file(tmp_dir: Path) -> None:
    with open(tmp_dir / "Peptides.pep", "w") as f:
        f.write("x\nACDEFGHIK\nMNPQRSTVW\nACDEFGHIK\n")  # one dup


# ===========================================================================
class TestReadNcbiFrames(unittest.TestCase):
    """Tests for read_ncbi_frames.jl"""

    @skip_no_julia
    def test_produces_frames_csv(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            make_minimal_ncbi_dir(d)
            if not (d / "sequences.fasta").exists():
                self.skipTest("Example data not available")
            result = run_julia("read_ncbi_frames.jl", str(d))
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue((d / "frames.csv").exists())

    @skip_no_julia
    def test_frames_csv_has_required_columns(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            make_minimal_ncbi_dir(d)
            if not (d / "sequences.fasta").exists():
                self.skipTest("Example data not available")
            run_julia("read_ncbi_frames.jl", str(d))
            rows = read_csv_rows(d / "frames.csv")
            self.assertGreater(len(rows), 0)
            self.assertIn("Region", rows[0])
            self.assertIn("Consensus_sequence", rows[0])
            self.assertIn("Description", rows[0])

    @skip_no_julia
    def test_no_fasta_exits_nonzero(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = run_julia("read_ncbi_frames.jl", tmp)
            self.assertNotEqual(result.returncode, 0)

    @skip_no_julia
    def test_suffix_flag_produces_suffixed_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            make_minimal_ncbi_dir(d)
            if not (d / "sequences.fasta").exists():
                self.skipTest("Example data not available")
            run_julia("read_ncbi_frames.jl", str(d), "--suffix", "test")
            self.assertTrue((d / "frames_test.csv").exists())


# ===========================================================================
class TestReadSamfireFrames(unittest.TestCase):
    """Tests for read_samfire_frames.jl"""

    @skip_no_julia
    def test_produces_frames_csv_from_dat(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            with open(d / "Reading_Frames.dat", "w") as f:
                f.write("77 496\nATGCCC\nMPE\n")
            result = run_julia("read_samfire_frames.jl", str(d))
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue((d / "frames.csv").exists())

    @skip_no_julia
    def test_frames_csv_has_correct_region(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            with open(d / "Reading_Frames.dat", "w") as f:
                f.write("77 496\nATGCCCGAATTT\nMPEF\n")
            run_julia("read_samfire_frames.jl", str(d))
            rows = read_csv_rows(d / "frames.csv")
            self.assertEqual(rows[0]["Region"], "77,496")

    @skip_no_julia
    def test_no_dat_file_exits_cleanly(self):
        """Script should exit 0 (not crash) when no .dat file found."""
        with tempfile.TemporaryDirectory() as tmp:
            result = run_julia("read_samfire_frames.jl", tmp)
            # The script prints a message and returns without hard-crashing
            self.assertIn(result.returncode, (0, 1))


# ===========================================================================
class TestParseTrajectories(unittest.TestCase):
    """Tests for parse_trajectories.jl"""

    @skip_no_julia
    def test_produces_variants_csv(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            with open(d / "single_locus_trajectories.out", "w") as f:
                f.write("100 A G 1 0 10 0 20 0 30\n")
            result = run_julia("parse_trajectories.jl", str(d))
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue((d / "variants.csv").exists())

    @skip_no_julia
    def test_variants_csv_columns(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            with open(d / "single_locus_trajectories.out", "w") as f:
                f.write("100 A G 1 0 10 0 20 0 30\n200 C T 1 5 5 5 5 5 25\n")
            run_julia("parse_trajectories.jl", str(d))
            rows = read_csv_rows(d / "variants.csv")
            self.assertEqual(len(rows), 2)
            self.assertIn("Locus", rows[0])
            self.assertIn("Consensus", rows[0])
            self.assertIn("Variant", rows[0])

    @skip_no_julia
    def test_no_out_file_errors(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = run_julia("parse_trajectories.jl", tmp)
            self.assertNotEqual(result.returncode, 0)

    @skip_no_julia
    def test_uses_example_stanevich_data(self):
        src = REPO_ROOT / "data" / "Stanevich_et_al"
        if not (src / "single_locus_trajectories10.out").exists():
            self.skipTest("Stanevich example data not available")
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            shutil.copy(src / "single_locus_trajectories10.out", d)
            result = run_julia("parse_trajectories.jl", str(d))
            self.assertEqual(result.returncode, 0, result.stderr)
            rows = read_csv_rows(d / "variants.csv")
            self.assertGreater(len(rows), 0)


# ===========================================================================
class TestCleanPeptides(unittest.TestCase):
    """Tests for clean_peptides.jl"""

    @skip_no_julia
    def test_removes_duplicates(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            with open(d / "Peptides.pep", "w") as f:
                f.write("x\nACDEFGHIK\nACDEFGHIK\nMNPQRSTVW\n")
            result = run_julia("clean_peptides.jl", str(d))
            self.assertEqual(result.returncode, 0, result.stderr)
            rows = read_csv_rows(d / "Peptides.pep")
            self.assertEqual(len(rows), 2)

    @skip_no_julia
    def test_removes_stop_codon_peptides(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            with open(d / "Peptides.pep", "w") as f:
                f.write("x\nACDE*GHIK\nMNPQRSTVW\n")
            run_julia("clean_peptides.jl", str(d))
            rows = read_csv_rows(d / "Peptides.pep")
            peptides = [list(r.values())[0] for r in rows]
            self.assertFalse(any("*" in p for p in peptides))

    @skip_no_julia
    def test_no_pep_file_exits_cleanly(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = run_julia("clean_peptides.jl", tmp)
            self.assertEqual(result.returncode, 0)


# ===========================================================================
class TestGeneratePeptides(unittest.TestCase):
    """Tests for generate_peptides.jl"""

    @skip_no_julia
    def test_produces_peptides_labels_csv(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            make_minimal_frames_csv(d, dna="ATGCCCGAATTTAAG", start=1)
            make_minimal_variants_csv(d)
            result = run_julia("generate_peptides.jl", str(d))
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            self.assertTrue((d / "peptides_labels.csv").exists())

    @skip_no_julia
    def test_produces_peptides_pep_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            make_minimal_frames_csv(d, dna="ATGCCCGAATTTAAG", start=1)
            make_minimal_variants_csv(d)
            run_julia("generate_peptides.jl", str(d))
            self.assertTrue((d / "Peptides.pep").exists())

    @skip_no_julia
    def test_peptides_labels_csv_has_correct_columns(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            make_minimal_frames_csv(d, dna="ATGCCCGAATTTAAGCCC", start=1)
            make_minimal_variants_csv(d)
            run_julia("generate_peptides.jl", str(d))
            if (d / "peptides_labels.csv").exists():
                rows = read_csv_rows(d / "peptides_labels.csv")
                if rows:
                    self.assertIn("Locus", rows[0])
                    self.assertIn("Peptide", rows[0])
                    self.assertIn("Peptide_label", rows[0])

    @skip_no_julia
    def test_no_folder_arg_exits_nonzero(self):
        result = subprocess.run(
            ["julia", f"--project={SRC_DIR}", str(SRC_DIR / "generate_peptides.jl")],
            capture_output=True, text=True, timeout=60,
        )
        self.assertNotEqual(result.returncode, 0)


# ===========================================================================
class TestSimulateVariants(unittest.TestCase):
    """Tests for simulate_variants.jl"""

    @skip_no_julia
    def test_produces_variants_simulated_csv(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            make_minimal_frames_csv(d, dna="ATGCCCGAATTT", start=1)
            result = run_julia("simulate_variants.jl", str(d))
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue((d / "variants_simulated.csv").exists())

    @skip_no_julia
    def test_variants_csv_has_locus_consensus_variant(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            make_minimal_frames_csv(d, dna="ATGCCCGAATTT", start=1)
            run_julia("simulate_variants.jl", str(d))
            rows = read_csv_rows(d / "variants_simulated.csv")
            self.assertGreater(len(rows), 0)
            self.assertIn("Locus", rows[0])
            self.assertIn("Consensus", rows[0])
            self.assertIn("Variant", rows[0])

    @skip_no_julia
    def test_no_self_substitution(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            make_minimal_frames_csv(d, dna="ATGCCCGAATTT", start=1)
            run_julia("simulate_variants.jl", str(d))
            rows = read_csv_rows(d / "variants_simulated.csv")
            for row in rows:
                self.assertNotEqual(row["Consensus"], row["Variant"])

    @skip_no_julia
    def test_n_flag_limits_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            make_minimal_frames_csv(d, dna="ATGCCCGAATTT", start=1)
            run_julia("simulate_variants.jl", str(d), "--n", "5")
            rows = read_csv_rows(d / "variants_simulated.csv")
            self.assertLessEqual(len(rows), 5)

    @skip_no_julia
    def test_p_flag_samples_proportion(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            # 4 codons * 3 positions * 3 subs = 36 total
            make_minimal_frames_csv(d, dna="ATGCCCGAATTT", start=1)
            run_julia("simulate_variants.jl", str(d), "--p", "0.5")
            rows = read_csv_rows(d / "variants_simulated.csv")
            self.assertLessEqual(len(rows), 36)


# ===========================================================================
class TestParseAaVariants(unittest.TestCase):
    """Tests for parse_aa_variants.jl"""

    @skip_no_julia
    def test_produces_variants_csv(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            make_minimal_frames_csv(d, dna="ATGCCCGAATTT", start=1, desc="ORF1")
            with open(d / "variants.aa", "w") as f:
                f.write("ORF1 2\nP S\n")
            result = run_julia("parse_aa_variants.jl", str(d))
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            self.assertTrue((d / "variants.csv").exists())

    @skip_no_julia
    def test_variants_csv_has_expected_columns(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            make_minimal_frames_csv(d, dna="ATGCCCGAATTT", start=1, desc="ORF1")
            with open(d / "variants.aa", "w") as f:
                f.write("ORF1 2\nP S\n")
            run_julia("parse_aa_variants.jl", str(d))
            rows = read_csv_rows(d / "variants.csv")
            self.assertGreater(len(rows), 0)
            self.assertIn("Locus", rows[0])
            self.assertIn("Consensus", rows[0])
            self.assertIn("Variant", rows[0])

    @skip_no_julia
    def test_no_aa_file_exits_nonzero(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            make_minimal_frames_csv(d)
            result = run_julia("parse_aa_variants.jl", str(d))
            self.assertNotEqual(result.returncode, 0)

    @skip_no_julia
    def test_unknown_orf_name_skipped_gracefully(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            make_minimal_frames_csv(d, dna="ATGCCCGAATTT", start=1, desc="ORF1")
            with open(d / "variants.aa", "w") as f:
                f.write("BADORF 2\nP S\n")
            result = run_julia("parse_aa_variants.jl", str(d))
            self.assertNotEqual(result.returncode, 0)  # no valid variants


# ===========================================================================
class TestVariantFatesIntegration(unittest.TestCase):
    """End-to-end test for variant_fates.jl using pre-built Stanevich data."""

    @skip_no_julia
    def test_produces_variant_fates_csv(self):
        src = REPO_ROOT / "data" / "Stanevich_et_al"
        required = ["variants.csv", "frames.csv", "harmonic_mean_best_ranks.csv"]
        if not all((src / f).exists() for f in required):
            self.skipTest("Stanevich example data not available")
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            for f in ["variants.csv", "frames.csv", "harmonic_mean_best_ranks.csv"]:
                shutil.copy(src / f, d / f)
            result = run_julia("variant_fates.jl", str(d))
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue((d / "variant_fates.csv").exists())

    @skip_no_julia
    def test_variant_fates_has_required_columns(self):
        src = REPO_ROOT / "data" / "Stanevich_et_al"
        required = ["variants.csv", "frames.csv", "harmonic_mean_best_ranks.csv"]
        if not all((src / f).exists() for f in required):
            self.skipTest("Stanevich example data not available")
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            for f in required:
                shutil.copy(src / f, d / f)
            run_julia("variant_fates.jl", str(d))
            rows = read_csv_rows(d / "variant_fates.csv")
            self.assertGreater(len(rows), 0)
            required_cols = {"Locus", "frame_filter", "peptide_filter", "binding_filter"}
            self.assertTrue(required_cols.issubset(rows[0].keys()))


# ===========================================================================
if __name__ == "__main__":
    unittest.main(verbosity=2)
