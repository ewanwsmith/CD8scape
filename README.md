# CD8scape

CD8scape runs netMHCpan on genetic variants, either for an individual's HLA genotype or for a representative supertype panel. It predicts whether a mutation weakens CD8+ T-cell recognition — i.e. immune escape.

## Features
- Peptide generation for consensus and variant loci
- MHC binding prediction with netMHCpan
- Output parsing and best-rank calculation
- Harmonic mean best rank (HMBR) and fold change across the allele panel
- Per-allele log2 fold change for every allele in the genome (`--per-allele`)
- Percentile benchmarking against a background variant distribution — real-world data by default, simulated variants as a fallback

## Requirements
- Perl 5
- Julia v1.11+
- [netMHCpan 4.2](https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/) (4.1 also works)
- Python 3.8+ and PyQt6 (desktop app only — PyQt6 installs itself on first launch)

## Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/ewanwsmith/CD8scape.git
   cd CD8scape
   ```
2. Copy `src/settings.txt.example` to `src/settings.txt` and set your netMHCpan path:
   ```bash
   cp src/settings.txt.example src/settings.txt
   # Edit src/settings.txt to set your local NETMHCPAN path
   ```
3. Install Julia dependencies and validate settings:
   ```bash
   ./CD8scape.jl prep
   ```
   This activates the project environment, installs the Julia packages, validates the netMHCpan path, and checks that Perl is available.

## Desktop app

There's a PyQt6 desktop app if you'd rather not use the command line. It's a thin wrapper: every action builds and runs a `./CD8scape.jl <command> [flags…]` call and streams the output live, so file discovery, variant parsing, and error messages all come from CD8scape itself, not from separate GUI code. See `ui/README.md` for the internals.

**Launch:**
- **macOS**: double-click `launch.command` in Finder (or run `python launch.py` from a terminal). A Terminal window stays open while the app runs.
- **Windows**: double-click `launch.bat`.

PyQt6 installs into your Python environment the first time you run the launcher — no manual `pip install`. `src/settings.txt` (the netMHCpan path) is created and validated from inside the app, so if you only use the GUI you can skip steps 2–3 above.

The app is a wizard, one page per pipeline stage:

| Step | Page | What it does |
|---|---|---|
| 1 | **Setup** | Enter your netMHCpan path, check Perl is available, and install Julia dependencies (same as `./CD8scape.jl prep`). Once per machine. |
| 2 | **Data** | Pick an example dataset or point to your own data folder. |
| 3 | **Prepare** | Choose input type (nucleotide VCF/Samfire or amino-acid `.aa`) and parsing options (same as `read`/`simulate`). |
| 4 | **Run** | Pick analysis type — individual genotype, supertype panel, or percentile benchmarking (against a real-world background file, or a simulated one) — plus options like `--per-allele`, `--verbose`, and thread count. |
| 5 | **Execute** | Run the pipeline and watch stdout/stderr live, with an estimated time remaining. |
| 6 | **Output** | Browse and open results, preview a table for the selected file, save everything as a ZIP, or delete the run's output. |

### Troubleshooting

- **PyQt6 won't install**: the launcher installs it into your current Python 3.8+ environment. If that fails, run `pip install PyQt6` yourself and re-launch.
- **netMHCpan path rejected on the Setup page**: the app validates the path exactly as `CD8scape.jl prep` does — point it at the netMHCpan executable, not the folder containing it.
- **macOS "app can't be opened" warning**: right-click `launch.command` → Open the first time to get past Gatekeeper (the script is unsigned).
- **Output page slow to load, or a very long delete confirmation**: this comes from `--verbose` runs, which keep per-allele logs and temp files. The delete dialog summarises large file lists rather than printing every name; use **Show Details…** to see the full list.

## Input Data

CD8scape reads a data folder containing:

- **alleles.txt**: HLA alleles, one per line (e.g. `HLA-A03:01`). Required for `run`.
- **supertype_panel.csv**: columns `Allele` and `Frequency` (and optionally `Locus`). Required for `run_supertype`. Drop it in the data folder to override the project default.
- **Variant file**, one of:
  - `.vcf` or `.vcf.gz` (standard variant call format).
  - `single_locus_trajectories.out` from [Samfire](https://github.com/cjri/samfire) (also matches any `single_locus_trajectories*.out`, or falls back to any `.out` file in the folder).
  - A `.aa` amino-acid variant file, used when `read` is run with `--aa` (see below).
- **Reading frame file**, one of:
  - `sequences.fasta` **and** `consensus.fa` for the NCBI path: `sequences.fasta` gives ORF definitions with coordinate headers (e.g. from [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/)), and `consensus.fa` gives the full reference genome that reading-frame subsequences are cut from.
  - `Reading_Frames.dat` from [Samfire](https://github.com/cjri/samfire).

CD8scape picks the right files automatically. It tries VCF first; if there's no VCF, or parsing fails, it falls back to Samfire trajectories. Pass `--aa` to `read` to read amino-acid variants instead. For reading frames it tries the NCBI path (`sequences.fasta` + `consensus.fa`) first, then falls back to Samfire's `Reading_Frames.dat`.

### Amino-acid variant file format (`.aa`)

The `.aa` format sets variants directly at the amino-acid level — for instance, to reproduce substitutions from a paper when you don't have the underlying sequence data. Pass `--aa` to `read` to use it.

Each variant is two lines:
```
<orf_name> <aa_position>
<ancestral_aa> <derived_aa>
```

- `orf_name` must match a `Description` value in `frames.csv` (set by the reading frame file) exactly.
- `aa_position` is 1-based within the translated protein.
- `ancestral_aa` and `derived_aa` are single-letter amino acid codes.
- Blank lines between records are ignored.

**Example:**
```
Orf3 23
K M

Orf1 45
A T
```

Canonical codons are used for both the ancestral and derived amino acids, regardless of the consensus sequence. So you can force a substitution even where the consensus at that position doesn't encode the stated ancestral amino acid (e.g. reproducing published results without the original data). When a consensus mismatch is overridden, CD8scape prints a warning and updates the frames file in place.

### Example: alleles.txt
```
HLA-A03:01
HLA-A01:01
HLA-B08:01
HLA-B07:02
HLA-C07:02
HLA-C07:01
```

### Folder Structure Examples

Individual genotype run:
```
<your_data_folder>/
    alleles.txt
    variants.vcf                    # or single_locus_trajectories.out
    sequences.fasta                 # NCBI ORF definitions
    consensus.fa                    # full reference genome
```

Supertype panel run:
```
<your_data_folder>/
    supertype_panel.csv
    variants.vcf                    # or single_locus_trajectories.out
    sequences.fasta
    consensus.fa
```

Samfire-only input:
```
<your_data_folder>/
    alleles.txt
    single_locus_trajectories.out
    Reading_Frames.dat
```

Amino-acid variant input (use `read --aa`):
```
<your_data_folder>/
    alleles.txt
    variants.aa                     # amino-acid variant file
    sequences.fasta                 # NCBI ORF definitions (or Reading_Frames.dat)
    consensus.fa                    # full reference genome
```

For Samfire, see the [Samfire GitHub](https://github.com/cjri/samfire) for how to generate the `.out` and `.dat` files.

## Usage
Run every command from the repository root.

### 1. Prepare Environment
```bash
./CD8scape.jl prep
```
`prep` does all the dependency installation and setup. The other commands (`read`, `simulate`, `run`, `run_supertype`, `percentile`) install nothing and assume the environment is already prepared.

### 2. Parse Input Data
```bash
./CD8scape.jl read <folder_path> [--aa] [--suffix <name>] [--latest|--no-latest]
```
Parses variants and reading frames from the data folder into `variants.csv` and `frames.csv`.
- `--aa`: read amino-acid variants from a `.aa` file instead of VCF or Samfire trajectories (see [Amino-acid variant file format](#amino-acid-variant-file-format-aa)).

### 3. Run Pipeline (Individual Genotype)
```bash
./CD8scape.jl run <folder_path> [--t <N|max>|--thread <N|max>] [--per-allele] [--verbose] [--suffix <name>] [--latest|--no-latest]
```
Generates peptides, runs netMHCpan, parses the output, and calculates best ranks and fold changes.
- `--t`/`--thread`: max parallel chunks for netMHCpan (default: 1). `max` uses all available threads up to the safety cap (override with `CD8SCAPE_MAX_THREADS`).
- `--verbose`: keep per-allele logs and temp files for debugging.
- `--per-allele`: compute log2 fold change for every allele in the genome separately. Only alleles with an ancestral EL rank ≤ 2% are included. Writes `per_allele_best_ranks.csv` with columns `Frame`, `Locus`, `Mutation`, `MHC`, `ELBR_A`, `ELBR_D`, `foldchange_BR`, `log2_foldchange_BR`.

### 4. Run Pipeline (Supertype Panel)
```bash
./CD8scape.jl run_supertype <folder_path> [--t <N|max>|--thread <N|max>] [--per-allele] [--verbose] [--suffix <name>] [--latest|--no-latest]
```
Same as `run`, but uses a representative supertype HLA panel.
- `--per-allele` works here too, but panel alleles are population-frequency surrogates, not an individual's genotype — so per-allele results reflect population coverage rather than one person's immunogenicity.

### 5. Simulate Input Data (optional — fallback background)

Percentile benchmarking (step 6) compares your observed variants against a background of comparator variants. Use a **real-world background** where you can: a `harmonic_mean_best_ranks.csv` from running steps 2 and 3/4 on a large panel of naturally occurring variants (e.g. surveillance or consensus data for the pathogen). `simulate` builds a **synthetic background** instead — for when no real-world dataset exists, or when you specifically want an unbiased null where every possible single-nucleotide substitution is equally likely. Skip this step if you already have a real-world background to pass to `percentile --s`.

```bash
./CD8scape.jl simulate <folder_path> [--n <count>] [--p <proportion>] [--seed <int>] [--suffix <name>] [--latest|--no-latest]
```
- Parses reading frames (`frames.csv`) and generates every single-nucleotide variant per reading frame (`variants.csv`).
- Sampling:
   - `--n <count>`: sample an absolute number of variants.
   - `--p <proportion>` (alias `--prop`): sample a proportion in (0,1).
   - If both are given, `--n` wins.
- Defaults: `--n` defaults to `1000` and `--p` to `0.1` when the flag is given without a value; omit both to write all variants. `--seed` sets the RNG seed (default: `1320`).
- After simulating, run step 3 or 4 again on the simulated data (with a `--suffix` like `simulated`) to produce the `harmonic_mean_best_ranks_simulated.csv` background used below.

### 6. Compute Percentiles (Benchmarking)
```bash
./CD8scape.jl percentile <folder_path> [--per-allele] [--s <background_file>] [--o <obs_file>]
```
- Computes observed log2 fold-change percentiles against a background of comparator variants.
- Works on HMBR fold changes (`harmonic_mean_best_ranks.csv`) by default.
- **Background — two options for `--s`:**
  - **Real-world (recommended)**: point `--s` at a `harmonic_mean_best_ranks(_suffix).csv` built by running `read` + `run`/`run_supertype` (steps 2 and 3/4) on naturally occurring variants. This benchmarks against variation actually seen in circulating strains. Skip `simulate` entirely in this mode.
  - **Simulated (fallback)**: run `simulate` + `run`/`run_supertype` (step 5, then 3/4 again) to build `harmonic_mean_best_ranks_simulated.csv` and pass that as `--s`. Use it when there's no real-world dataset, or when you want an exhaustive/sampled set of synthetic substitutions.
- `--per-allele`: work on per-allele fold changes instead (`per_allele_best_ranks.csv`). Percentiles are computed per allele rather than per variant.
- `--s <background_file>`: the background comparator file. Defaults to `harmonic_mean_best_ranks_simulated.csv` (or `per_allele_best_ranks_simulated.csv` with `--per-allele`) — pass a real-world file explicitly to use one.
- `--o <obs_file>`: the observed file (defaults to the most recent match in the data folder, excluding `_simulated`).
- Observed variants are dropped from the background before percentiles are computed.
- Writes `percentile_harmonic_mean_best_ranks.csv` or `percentile_per_allele_best_ranks.csv`. Added columns: `Percentile` (0–100), `Z_i` (the percentile re-expressed as a standard-normal score), and `p_value` (empty for individual variants). Two summary rows are appended:
  - `combined_z`: a parametric combined score (Stouffer's method) testing whether the observed variants sit higher in the background distribution than expected, with a one-tailed p-value.
  - `empirical_p`: a non-parametric permutation p-value that repeatedly draws random sets of variants from the background and compares their mean percentile to the observed mean.

The desktop app exposes the real-world path as a **"Use pre-existing background files"** toggle on the Run page. Turning it on skips the simulate step and the extra run-on-simulated pass, and calls `percentile` straight against the background and observed CSVs you point it at — 3 steps instead of 5.

### Global Options

These are shared across `read`, `simulate`, `run`, and `run_supertype`:

- **`--suffix <name>`**: insert `_<name>` before the extension of every output file. `--suffix foo` gives `variants_foo.csv`, `best_ranks_foo.csv`, `harmonic_mean_best_ranks_foo.csv`, and so on. For `simulate`, the suffix defaults to `simulated`.
- **`--latest`** (default) / **`--no-latest`**: how input files are resolved when there's no `--suffix` and several candidates exist (e.g. `frames.csv` and `frames_simulated.csv`). `--latest` takes the most recently modified file; `--no-latest` errors on ambiguity.

Together these let several analyses (e.g. observed vs. simulated) sit in one data folder without overwriting each other.

## Workflow Summary
1. **prep**: install dependencies.
2. **read**: parse variants and frames from real data.
3. **run/run_supertype**: generate peptides, predict binding, process output, and calculate best ranks and fold changes for your observed data. Add `--per-allele` for per-allele fold changes.
4. **simulate** *(optional, right before percentile)*: generate a synthetic background of single-nucleotide variants, then run step 3 on it (with a `--suffix`). Only needed if you don't already have a real-world background — skip to step 5 if you do.
5. **percentile**: compare observed fold changes to a background — real-world by default (`--s <background_file>`, built by running steps 2–3 on a real variant panel), or the simulated one from step 4. Add `--per-allele` to benchmark per-allele fold changes instead of HMBR.

## How it works

netMHCpan reports an eluted-ligand percentile rank (EL %rank) for each peptide–allele pair — how the peptide's binding score compares against a background of random natural peptides for that allele. Lower ranks mean stronger predicted presentation to CD8+ T cells. CD8scape counts a peptide as a binder at an EL %rank of 2% or better (netMHCpan's weak-binder cut-off).

**Peptide generation.** For each non-synonymous variant, CD8scape enumerates every 8–11-mer peptide whose window covers the mutated residue, in both the ancestral (consensus) and derived (variant) states. Peptides that are synonymous between the two states, or that contain a stop codon, are dropped.

**Best rank.** netMHCpan scores every peptide against each allele. For a given variant, state, and allele, CD8scape keeps the strongest binder — the lowest EL %rank among the peptides covering the mutation. These are the per-allele best ranks `ELBR_A` (ancestral) and `ELBR_D` (derived).

**Harmonic mean best rank (HMBR).** Best ranks are pooled across the allele panel using a harmonic mean, which is dominated by the strongest binders — appropriate because recognition hinges on the single best-presented peptide. An individual genotype weights its alleles equally; the supertype panel weights them by population frequency.

**Fold change.** The escape signal for each variant is the ratio of derived to ancestral HMBR, reported on a log2 scale. A positive log2 fold change means the derived peptide binds more weakly than the ancestral one — predicted escape. Loci where neither state binds (both HMBR above 2) are removed first. With `--per-allele`, the same ratio is computed separately for each allele, restricted to alleles where the ancestral peptide is at least a weak binder.

**Percentile benchmarking.** Observed fold changes are ranked against a background of comparator variants — a real-world panel by default, or a simulated single-nucleotide set — after removing any background entry that matches an observed variant. Each observed variant gets a percentile from 0 to 100: its position in the background distribution. CD8scape also reports two summary measures of whether the observed variants are collectively shifted toward escape: a parametric combined score (`combined_z`, Stouffer's method over the per-variant percentiles) and a non-parametric permutation p-value (`empirical_p`) that repeatedly resamples the background for comparison.

## Output Files
- `variants.csv`, `frames.csv`: parsed input data.
- `Peptides.pep`, `peptides_labels.csv`: generated peptides and labels.
- `netMHCpan_output.tsv`, `processed_output.csv`: raw and processed netMHCpan results.
- `best_ranks.csv`: per-allele best EL ranks for ancestral and derived peptides at each locus.
- `harmonic_mean_best_ranks.csv`: HMBR and log2 fold changes pooled across the panel (`Frame`, `Locus`, `Mutation`, `HMBR_A`, `HMBR_D`, `foldchange_HMBR`, `log2_foldchange_HMBR`).
- `per_allele_best_ranks.csv`: per-allele eluted-ligand best ranks and log2 fold changes for every allele in the genome, filtered to ancestral EL rank ≤ 2% (`Frame`, `Locus`, `Mutation`, `MHC`, `ELBR_A`, `ELBR_D`, `foldchange_BR`, `log2_foldchange_BR`). Written when `--per-allele` is passed to `run` or `run_supertype`.
- `variants_simulated.csv`, `harmonic_mean_best_ranks_simulated.csv`: simulated background variants and HMBR results (from `simulate` + `run`, used as the `--s` fallback background when there's no real-world dataset). A real-world background file has the same shape as `harmonic_mean_best_ranks.csv` but comes from a naturally occurring variant panel — name it whatever you like and pass it via `--s`.
- `per_allele_best_ranks_simulated.csv`: per-allele results for the simulated background (from `simulate` + `run --per-allele`).
- `percentile_harmonic_mean_best_ranks.csv`: observed HMBR with `Percentile`, `Z_i`, and `p_value` columns relative to the background, plus the `combined_z` (parametric combined score) and `empirical_p` (permutation p-value) summary rows.
- `percentile_per_allele_best_ranks.csv`: observed per-allele fold changes with the same columns and summary rows (from `percentile --per-allele`).

## Advanced Configuration

Set these environment variables to override internal defaults for large-panel or resource-constrained runs:

| Variable | Default | Description |
|---|---|---|
| `CD8SCAPE_MAX_THREADS` | 8 | Safety cap on parallel netMHCpan chunks when `--t max` is used. |
| `CD8SCAPE_ALLELE_CHAR_LIMIT` | 1023 | Max character length of the allele string passed to a single netMHCpan call. Alleles are batched past this limit. |
| `CD8SCAPE_ALLELE_COUNT_LIMIT` | 75 | Max number of alleles passed to a single netMHCpan call. Alleles are batched past this limit. |

## Example Data

`data/Example_data/` is a minimal synthetic dataset for testing the pipeline end to end. It has a 53-amino-acid ORF, two amino-acid variants (`.aa` format), a three-allele genotype, and a small supertype panel.

Real-world SARS-CoV-2 data used in development is in `data/Stanevich_et_al/`.

## Citation
If you use CD8scape, please cite the repository and netMHCpan.

The real-world example data (`data/Stanevich_et_al/`) is from:
Stanevich, O.V., Alekseeva, E.I., Sergeeva, M. et al. SARS-CoV-2 escape from cytotoxic T cells during long-term COVID-19. Nat Commun 14, 149 (2023). https://doi.org/10.1038/s41467-022-34033-x


## License
GNU Public License
