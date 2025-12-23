## Quick Start

- Prep environment and external tools:
   - Activates the `src` Julia environment, installs deps, and checks NetMHCpan + Perl.

```bash
./CD8scape.jl prep
```

- Simulate variants from frames.csv (writes `variants.csv`):

```bash
julia --project=src src/simulate_variants.jl /Users/e.smith.5/Desktop/CD8scape_example_copy
# optional flags: --n 1000  |  --p 0.1  |  --seed 42
```

- Generate and clean peptides (writes `peptides_labels.csv` and `Peptides.pep`):

```bash
julia --project=src src/generate_peptides.jl /Users/e.smith.5/Desktop/CD8scape_example_copy
julia --project=src src/clean_peptides.jl   /Users/e.smith.5/Desktop/CD8scape_example_copy
```

- Run NetMHCpan over peptides (writes `netMHCpan_output.tsv`):
   - For specific alleles, ensure `alleles.txt` exists in the folder; NetMHCpan path set in `src/settings.txt`.

```bash
julia --project=src src/run_netMHCpan.jl --folder /Users/e.smith.5/Desktop/CD8scape_example_copy
```

- Run NetMHCpan with global supertype panel:

```bash
julia --project=src src/run_netMHCpan_global.jl --folder /Users/e.smith.5/Desktop/CD8scape_example_copy
```

- Parse a VCF to variants.csv (if using VCF instead of simulation):

```bash
julia --project=src src/parse_vcf.jl /Users/e.smith.5/Desktop/CD8scape_example_copy
```
# CD8scape

CD8scape runs netMHCpan on genetic variants for individual HLA genotypes or representative supertype panels.

## Features
- Automated peptide generation for consensus and variant loci
- MHC binding prediction using netMHCpan
- Robust output parsing and best-rank calculation
- Harmonic mean best rank (HMBR) and fold change analysis
 

## Requirements
- Perl 5
- Julia v1.11+
- [netMHCpan 4.1](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/)

## Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/ewanwsmith/CD8scape.git
   cd CD8scape
   ```
2. Copy `src/settings.txt.example` to `src/settings.txt` and set your netMHCpan path.
   ```bash
   cp src/settings.txt.example src/settings.txt
   # Edit src/settings.txt to set your local NETMHCPAN path
   ```
3. Install Julia and Perl dependencies:
   ```bash
   ./CD8scape.jl prep
   ```

## Input Data

CD8scape expects a folder containing the following files:

- **alleles.txt**: List of HLA alleles, one per line (e.g., `HLA-A01:01`).
- **Variant file**: Either
  - `single_locus_trajectories.out` (from [Samfire](https://github.com/cjri/samfire)), or
  - Any `.out` file with Samfire-compatible format, or
  - `.vcf` or `.vcf.gz` file (standard variant call format).
- **Reading frame file**: Either
  - `Reading_Frames.dat` (from [Samfire](https://github.com/cjri/samfire)), or
  - `sequences.fasta` (reference genome from [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/)).

### Example: alleles.txt
```
HLA-A01:01
HLA-A02:01
HLA-B07:02
HLA-C07:02
```

### Folder Structure Example
```
<your_data_folder>/
    alleles.txt
    single_locus_trajectories.out   # or variants.vcf
    Reading_Frames.dat              # or sequences.fasta
```

CD8scape will automatically detect and use the appropriate files for variant and reading frame parsing. For Samfire, see [Samfire GitHub](https://github.com/cjri/samfire) for details on generating `.out` and `.dat` files.

## Usage
All commands are run from the repository root:

### 1. Prepare Environment
```bash
./CD8scape.jl prep
```
Note: `prep` performs all dependency installation and environment setup. The other commands (`read`, `simulate`, `run`, `run_supertype`) do not install packages and assume the environment is already prepared.

### 2. Parse Input Data
```bash
./CD8scape.jl read <folder_path>
```
- Parses variants and reading frames, outputs `variants.csv` and `frames.csv`.

### 3. Simulate Input Data
```bash
./CD8scape.jl simulate <folder_path> [--n <count>] [--p <proportion>] [--seed <int>]
```
- Parses reading frames (writes `frames.csv`) and generates exhaustive simulated single-nucleotide variants per reading frame (writes `variants.csv`).
- Sampling options:
   - `--n <count>`: sample an absolute number of variants.
   - `--p <proportion>` (alias `--prop`): sample a proportion in (0,1).
   - If both `--n` and `--p` are provided, `--n` takes precedence.
 - Defaults: `--n` defaults to `1000` and `--p` defaults to `0.1` when you choose to sample by count or proportion; no sampling flags writes all variants. `--seed` sets RNG seed (default: `1320`).

### 4. Run Pipeline (Individual Genotype)
```bash
./CD8scape.jl run <folder_path> [--t <N>|--thread <N>] [--verbose]
```
- Generates peptides, runs netMHCpan, parses output, calculates best ranks and fold changes.
- `--t <N|max>`/`--thread <N|max>` runs up to N peptide chunks in parallel (default 1). Use `max` to request the capped maximum. For safety, concurrency is capped by default to half of logical CPUs; if you request above the cap, it will be reduced with a notice. You can adjust the cap with `CD8SCAPE_MAX_THREADS`.
- `--verbose` preserves per-allele logs and temp files for inspection.

### 5. Run Pipeline (Supertype Panel)
```bash
./CD8scape.jl run_supertype <folder_path> [--t <N>|--thread <N>] [--verbose]
```
- As above, but uses a representative supertype HLA panel.
- `--t <N|max>`/`--thread <N|max>` runs up to N peptide chunks in parallel (default 1). Use `max` to request the capped maximum. For safety, concurrency is capped by default to half of logical CPUs; if you request above the cap, it will be reduced with a notice. You can adjust the cap with `CD8SCAPE_MAX_THREADS`.

### 6. Percentile Scoring
Compute the percentile of observed variant HMBR log2 fold-change values against a simulated background distribution.

```bash
./CD8scape.jl percentile <folder_path> [--s <sim_file>] [--o <obs_file>]
```

- **Simulation input**: defaults to `harmonic_mean_best_ranks_simulated.csv` in `<folder_path>`. Override with `--s` (relative to folder or absolute path).
- **Observed input**: defaults to the most recent `harmonic_mean_best_ranks*.csv` in `<folder_path>` whose suffix is not `_simulated`. Override with `--o`.
- **Overlap handling**: simulated entries that exactly match observed variants (by `Frame|Locus|Mutation` → or `Locus|Mutation` → or `Locus`) are removed before scoring.
- **Output**: writes `percentile_harmonic_mean_best_ranks(_<suffix>).csv` in `<folder_path>`, adding a new `Percentile` column (0–100).

### Suffixes and Multi‑Run Folders
- Purpose: Keep outputs from different runs side‑by‑side in the same data folder without clobbering files.
- Flag: `--suffix <name>` inserts `_name` before file extensions for outputs, and is preferred for inputs when present.
- Discovery: If the suffixed input is not present, the tool falls back to the most recent matching file when `--latest` is used (default). Disable with `--no-latest`.

Behavior by stage
- read: Writes `frames_<suffix>.csv` and `variants_<suffix>.csv`.
- simulate: Defaults to `--suffix simulated` when none provided; writes `frames_simulated.csv` and `variants_simulated.csv`.
- run: Prefers `frames_<suffix>.csv` and `variants_<suffix>.csv` if they exist; otherwise falls back to latest `frames*.csv` and `variants*.csv`. Outputs are suffixed: `peptides_labels_<suffix>.csv`, `netMHCpan_output_<suffix>.tsv`, `processed_peptides_<suffix>.csv`, etc.
- run_supertype: Same suffix handling and outputs as `run` but for the supertype panel.

Examples
```bash
# Read with a tag
./CD8scape.jl read /path/to/data --suffix readtag

# Simulate with default suffix 'simulated'
./CD8scape.jl simulate /path/to/data --n 100

# Simulate with a custom suffix
./CD8scape.jl simulate /path/to/data --suffix simtag --n 50

# Run using the simulated tag (produces suffixed outputs)
./CD8scape.jl run /path/to/data --suffix simulated

# Run with a new tag; if frames_<tag>.csv and variants_<tag>.csv don't exist,
# generate_peptides will fall back to the most recent frames*/variants*.
./CD8scape.jl run /path/to/data --suffix fallbacktest

# Run supertype with its own tag (same suffix-first, latest-fallback behavior)
./CD8scape.jl run_supertype /path/to/data --suffix sttag
```

 

## Workflow Summary
1. **prep**: Install dependencies
2. **read**: Parse variants and frames
3. **run/run_supertype**: Generate peptides, predict binding, process output, calculate best ranks and fold changes
 

## Output Files
- `variants.csv`, `frames.csv`: Parsed input data
- `Peptides.pep`, `peptides_labels.csv`: Generated peptides and labels
- `netMHCpan_output.tsv`, `processed_output.csv`: Raw and processed netMHCpan results
- `best_ranks.csv`, `harmonic_mean_best_ranks.csv`: Best ranks and fold change analysis
- `percentile_harmonic_mean_best_ranks.csv`: Observed fold-change annotated with Percentile vs simulated background


## Removed
- The context-based pipeline under `src/context_run` was removed in favor of the simulation workflow.
 

## Citation
If you use CD8scape in your research, please cite the repository and netMHCpan as appropriate.

Example data is from:
Stanevich, O.V., Alekseeva, E.I., Sergeeva, M. et al. SARS-CoV-2 escape from cytotoxic T cells during long-term COVID-19. Nat Commun 14, 149 (2023). https://doi.org/10.1038/s41467-022-34033-x


## License
GNU Public License