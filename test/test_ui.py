"""
test_ui.py — Unit tests for CD8scape UI modules.

Tests every non-GUI module in ui/:
    paths.py            folder resolution and listing
    setup.py            netMHCpan path read/write/validate, env check
    workflow.py         workflow choice definitions and lookup
    runner.py           Julia detection, command building, repo root
    progress_estimator  thread inference and progress parsing

Run from the repo root:
    python -m pytest test/test_ui.py -v
    python test/test_ui.py
"""

from __future__ import annotations

import os
import sys
import tempfile
import unittest
from pathlib import Path

# ---------------------------------------------------------------------------
# Make ui/ importable regardless of where tests are run from
# ---------------------------------------------------------------------------

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "ui"))

from paths import FolderCheck, list_folder, resolve_folder
from setup import (
    SETTINGS_PATH,
    check_env,
    read_netmhcpan_path,
    validate_netmhcpan,
    write_netmhcpan_path,
)
from workflow import PREP_CHOICES, RUN_CHOICES, WorkflowChoice, choice_by_key
from runner import (
    CD8scapeNotFoundError,
    JuliaNotFoundError,
    build_command,
    find_julia,
    get_repo_root,
)
from progress_estimator import ProgressEstimator, infer_threads_from_args


# ===========================================================================
# paths.py
# ===========================================================================


class TestResolveFolderRejectsInvalid(unittest.TestCase):
    def test_empty_string(self):
        r = resolve_folder("")
        self.assertFalse(r.ok)
        self.assertIn("enter", r.message.lower())

    def test_whitespace_only(self):
        r = resolve_folder("   ")
        self.assertFalse(r.ok)

    def test_nonexistent_path(self):
        r = resolve_folder("/absolutely/does/not/exist/xyz123")
        self.assertFalse(r.ok)
        self.assertIn("exist", r.message.lower())

    def test_file_not_folder(self):
        r = resolve_folder(__file__)
        self.assertFalse(r.ok)
        self.assertIn("file", r.message.lower())


class TestResolveFolderAcceptsValid(unittest.TestCase):
    def test_test_directory_itself(self):
        r = resolve_folder(str(Path(__file__).parent))
        self.assertTrue(r.ok)
        self.assertEqual(r.message, "")

    def test_repo_root(self):
        r = resolve_folder(str(REPO_ROOT))
        self.assertTrue(r.ok)

    def test_tilde_expansion(self):
        r = resolve_folder("~")
        self.assertTrue(r.ok)

    def test_strips_quotes(self):
        # Finder on macOS sometimes wraps drag-and-dropped paths in quotes
        r = resolve_folder(f'"{REPO_ROOT}"')
        self.assertTrue(r.ok)


class TestListFolder(unittest.TestCase):
    def setUp(self):
        self.test_dir = Path(__file__).parent

    def test_returns_nonempty_for_test_dir(self):
        entries, total = list_folder(self.test_dir)
        self.assertGreater(total, 0)
        self.assertGreater(len(entries), 0)

    def test_folders_have_trailing_slash(self):
        # If any subdirectory exists inside, its name ends with /
        entries, _ = list_folder(REPO_ROOT)
        dirs = [e for e in entries if e.endswith("/")]
        self.assertGreater(len(dirs), 0, "Expected at least one subdirectory in repo root")

    def test_limit_caps_returned_entries(self):
        entries, total = list_folder(REPO_ROOT, limit=2)
        self.assertLessEqual(len(entries), 2)
        self.assertGreaterEqual(total, len(entries))

    def test_limit_none_returns_all(self):
        entries, total = list_folder(self.test_dir, limit=None)
        self.assertEqual(len(entries), total)

    def test_known_files_present(self):
        entries, _ = list_folder(self.test_dir, limit=None)
        names = [e.rstrip("/") for e in entries]
        self.assertIn("test_ui.py", names)


# ===========================================================================
# setup.py
# ===========================================================================


class TestValidateNetmhcpan(unittest.TestCase):
    def test_empty_path(self):
        r = validate_netmhcpan("")
        self.assertFalse(r.ok)

    def test_whitespace_path(self):
        r = validate_netmhcpan("   ")
        self.assertFalse(r.ok)

    def test_nonexistent_path(self):
        r = validate_netmhcpan("/nonexistent/netMHCpan")
        self.assertFalse(r.ok)
        self.assertIn("not found", r.message.lower())

    def test_directory_not_file(self):
        r = validate_netmhcpan(str(REPO_ROOT))
        self.assertFalse(r.ok)
        self.assertIn("directory", r.message.lower())

    def test_non_executable_file(self):
        with tempfile.NamedTemporaryFile(delete=False) as f:
            tmp = Path(f.name)
        try:
            tmp.chmod(0o644)  # readable but not executable
            r = validate_netmhcpan(str(tmp))
            self.assertFalse(r.ok)
            self.assertIn("executable", r.message.lower())
        finally:
            tmp.unlink(missing_ok=True)

    def test_executable_file(self):
        with tempfile.NamedTemporaryFile(delete=False) as f:
            tmp = Path(f.name)
        try:
            tmp.chmod(0o755)
            r = validate_netmhcpan(str(tmp))
            self.assertTrue(r.ok)
            self.assertEqual(r.resolved, tmp.resolve())
        finally:
            tmp.unlink(missing_ok=True)


class TestReadWriteNetmhcpanPath(unittest.TestCase):
    """Write/read roundtrip using a temporary settings file."""

    def _with_tmp_settings(self, fn):
        """Run fn(tmp_path) with SETTINGS_PATH temporarily replaced."""
        import setup as _setup
        orig = _setup.SETTINGS_PATH
        with tempfile.NamedTemporaryFile(suffix=".txt", delete=False) as f:
            tmp = Path(f.name)
        try:
            tmp.unlink()  # start with no file
            _setup.SETTINGS_PATH = tmp
            fn(tmp)
        finally:
            _setup.SETTINGS_PATH = orig
            tmp.unlink(missing_ok=True)

    def test_missing_file_returns_none(self):
        def check(_):
            self.assertIsNone(read_netmhcpan_path())
        self._with_tmp_settings(check)

    def test_write_then_read_roundtrip(self):
        def check(_):
            write_netmhcpan_path("/fake/netMHCpan")
            self.assertEqual(read_netmhcpan_path(), "/fake/netMHCpan")
        self._with_tmp_settings(check)

    def test_write_overwrites_existing_value(self):
        def check(_):
            write_netmhcpan_path("/old/path")
            write_netmhcpan_path("/new/path")
            self.assertEqual(read_netmhcpan_path(), "/new/path")
        self._with_tmp_settings(check)

    def test_placeholder_returns_none(self):
        def check(tmp):
            tmp.write_text("NETMHCPAN=full/path/to/netMHCpan\n")
            self.assertIsNone(read_netmhcpan_path())
        self._with_tmp_settings(check)

    def test_comment_lines_ignored(self):
        def check(tmp):
            tmp.write_text("# this is a comment\nNETMHCPAN=/real/path\n")
            self.assertEqual(read_netmhcpan_path(), "/real/path")
        self._with_tmp_settings(check)


class TestCheckEnv(unittest.TestCase):
    def test_returns_env_status(self):
        status = check_env()
        # Perl is almost always available on macOS/Linux
        self.assertIsInstance(status.perl_ok, bool)
        self.assertIsInstance(status.netmhcpan_ok, bool)

    def test_perl_found_on_macos_linux(self):
        import shutil
        if shutil.which("perl"):
            status = check_env()
            self.assertTrue(status.perl_ok)
            self.assertIsNotNone(status.perl_path)

    def test_ready_requires_both(self):
        status = check_env()
        self.assertEqual(status.ready, status.perl_ok and status.netmhcpan_ok)


# ===========================================================================
# workflow.py
# ===========================================================================


class TestWorkflowChoices(unittest.TestCase):
    def test_prep_choices_nonempty(self):
        self.assertGreater(len(PREP_CHOICES), 0)

    def test_run_choices_nonempty(self):
        self.assertGreater(len(RUN_CHOICES), 0)

    def test_all_choices_have_required_fields(self):
        for choice in PREP_CHOICES + RUN_CHOICES:
            self.assertIsInstance(choice.key, str)
            self.assertIsInstance(choice.label, str)
            self.assertIsInstance(choice.detail, str)
            self.assertIsInstance(choice.cli_command, str)
            self.assertIsInstance(choice.cli_fixed_args, list)

    def test_prep_keys_are_unique(self):
        keys = [c.key for c in PREP_CHOICES]
        self.assertEqual(len(keys), len(set(keys)))

    def test_run_keys_are_unique(self):
        keys = [c.key for c in RUN_CHOICES]
        self.assertEqual(len(keys), len(set(keys)))

    def test_expected_prep_commands_present(self):
        commands = {c.cli_command for c in PREP_CHOICES}
        self.assertIn("read", commands)
        self.assertIn("simulate", commands)

    def test_expected_run_commands_present(self):
        commands = {c.cli_command for c in RUN_CHOICES}
        self.assertIn("run", commands)
        self.assertIn("run_supertype", commands)

    def test_amino_acid_choice_has_aa_flag(self):
        aa = next((c for c in PREP_CHOICES if "--aa" in c.cli_fixed_args), None)
        self.assertIsNotNone(aa, "Expected a PREP_CHOICE with --aa fixed arg")


class TestChoiceByKey(unittest.TestCase):
    def test_finds_valid_key(self):
        key = PREP_CHOICES[0].key
        result = choice_by_key(PREP_CHOICES, key)
        self.assertEqual(result.key, key)

    def test_raises_on_invalid_key(self):
        with self.assertRaises(KeyError):
            choice_by_key(PREP_CHOICES, "completely_nonexistent_key")

    def test_finds_each_prep_choice(self):
        for c in PREP_CHOICES:
            found = choice_by_key(PREP_CHOICES, c.key)
            self.assertIs(found, c)

    def test_finds_each_run_choice(self):
        for c in RUN_CHOICES:
            found = choice_by_key(RUN_CHOICES, c.key)
            self.assertIs(found, c)


# ===========================================================================
# runner.py
# ===========================================================================


class TestGetRepoRoot(unittest.TestCase):
    def test_returns_path_containing_cd8scape_jl(self):
        root = get_repo_root()
        self.assertTrue((root / "CD8scape.jl").exists(),
                        f"CD8scape.jl not found at {root}")

    def test_returns_path_object(self):
        self.assertIsInstance(get_repo_root(), Path)


class TestFindJulia(unittest.TestCase):
    def test_returns_string_or_raises(self):
        try:
            julia = find_julia()
            self.assertIsInstance(julia, str)
            self.assertTrue(Path(julia).exists(), f"Julia path does not exist: {julia}")
        except JuliaNotFoundError:
            pass  # acceptable — Julia may not be installed in CI

    def test_error_message_mentions_julia(self):
        import shutil
        if shutil.which("julia") is None:
            with self.assertRaises(JuliaNotFoundError) as ctx:
                find_julia()
            self.assertIn("julia", str(ctx.exception).lower())


class TestBuildCommand(unittest.TestCase):
    @unittest.skipUnless(
        __import__("shutil").which("julia"), "Julia not on PATH"
    )
    def test_returns_list_starting_with_julia(self):
        cmd = build_command(["--help"])
        self.assertIsInstance(cmd, list)
        self.assertTrue(cmd[0].endswith("julia") or "julia" in cmd[0].lower())

    @unittest.skipUnless(
        __import__("shutil").which("julia"), "Julia not on PATH"
    )
    def test_contains_cd8scape_jl(self):
        cmd = build_command(["--help"])
        self.assertTrue(any("CD8scape.jl" in part for part in cmd))

    @unittest.skipUnless(
        __import__("shutil").which("julia"), "Julia not on PATH"
    )
    def test_appends_args(self):
        cmd = build_command(["read", "/some/folder", "--aa"])
        self.assertIn("read", cmd)
        self.assertIn("--aa", cmd)


# ===========================================================================
# progress_estimator.py
# ===========================================================================


class TestInferThreads(unittest.TestCase):
    def test_no_args_returns_one(self):
        self.assertEqual(infer_threads_from_args([]), 1)

    def test_t_flag_with_number(self):
        self.assertEqual(infer_threads_from_args(["--t", "4"]), 4)

    def test_thread_flag_with_number(self):
        self.assertEqual(infer_threads_from_args(["--thread", "8"]), 8)

    def test_t_max_returns_positive(self):
        n = infer_threads_from_args(["--t", "max"])
        self.assertGreater(n, 0)

    def test_unrelated_args_return_one(self):
        self.assertEqual(infer_threads_from_args(["--aa", "--per-allele"]), 1)

    def test_invalid_number_returns_one(self):
        self.assertEqual(infer_threads_from_args(["--t", "notanumber"]), 1)


class TestProgressEstimator(unittest.TestCase):
    def test_initial_progress_is_zero(self):
        est = ProgressEstimator(include_percentile=False, n_threads=1)
        self.assertEqual(est.overall_fraction, 0.0)
        self.assertEqual(est.overall_pct, 0)

    def test_complete_both_steps_gives_full_progress(self):
        est = ProgressEstimator(include_percentile=False, n_threads=1)
        est.start_step(0)
        est.complete_step(0, success=True)
        est.start_step(1)
        est.complete_step(1, success=True)
        self.assertAlmostEqual(est.overall_fraction, 1.0)
        self.assertEqual(est.overall_pct, 100)

    def test_two_step_mode_has_two_steps(self):
        self.assertEqual(ProgressEstimator(include_percentile=False).n_steps, 2)

    def test_five_step_mode_has_five_steps(self):
        self.assertEqual(ProgressEstimator(include_percentile=True).n_steps, 5)

    def test_overall_fraction_bounded(self):
        est = ProgressEstimator(include_percentile=False, n_threads=4)
        self.assertGreaterEqual(est.overall_fraction, 0.0)
        self.assertLessEqual(est.overall_fraction, 1.0)

    def test_weights_sum_to_100_pct(self):
        # round() on floats can produce 99 or 101 for 5-step; allow ±1
        total = sum(ProgressEstimator(include_percentile=False).step_weights_pct())
        self.assertIn(total, (99, 100, 101))

    def test_five_step_weights_sum_to_100_pct(self):
        total = sum(ProgressEstimator(include_percentile=True).step_weights_pct())
        self.assertIn(total, (99, 100, 101))

    def test_ingest_sequential_nmhc_line_advances_progress(self):
        est = ProgressEstimator(include_percentile=False, n_threads=1)
        est.start_step(1)
        before = est.overall_fraction
        est.ingest_line("Chunk 1/2, batch 1/1 (100 peptides, 6 alleles). 50% complete.")
        self.assertGreater(est.overall_fraction, before)

    def test_ingest_parallel_nmhc_line_advances_progress(self):
        est = ProgressEstimator(include_percentile=False, n_threads=4)
        est.start_step(1)
        before = est.overall_fraction
        est.ingest_line("2/4 units complete. 50%")
        self.assertGreater(est.overall_fraction, before)

    def test_ingest_post_processing_milestone_advances_progress(self):
        est = ProgressEstimator(include_percentile=False, n_threads=1)
        est.start_step(1)
        est.ingest_line("Chunk 1/1, batch 1/1 (100 peptides, 6 alleles). 100% complete.")
        before = est.overall_fraction
        est.ingest_line("Sorting by Locus...")
        self.assertGreater(est.overall_fraction, before)

    def test_n_threads_property(self):
        self.assertEqual(ProgressEstimator(n_threads=8).n_threads, 8)

    def test_include_percentile_property(self):
        self.assertTrue(ProgressEstimator(include_percentile=True).include_percentile)
        self.assertFalse(ProgressEstimator(include_percentile=False).include_percentile)


# ===========================================================================
# Runner
# ===========================================================================

if __name__ == "__main__":
    unittest.main(verbosity=2)
