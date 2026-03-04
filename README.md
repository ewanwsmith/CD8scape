# CD8scape

CD8scape runs netMHCpan on genetic variants for individual HLA genotypes or representative supertype panels.

## Features
- Automated peptide generation for consensus and variant loci
- MHC binding prediction using netMHCpan
- Robust output parsing and best-rank calculation
- Harmonic mean best rank (HMBR) and fold change analysis
- Max escape allele identification per variant (with biological caveat for supertype panel mode)
- Simulated variant generation for percentile benchmarking

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
2. Copy `src/settings.txt.example` to `src/settings.txt` and set your netMHCpan path:
   ```bash
   cp src/settings.txt.example src/settings.txt
   # Edit src/settings.txt to set your local NETMHCPAN path
   ```
3. Install Julia dependencies and validate settings:
   ```bash
   ./CD8scape.jl prep
   ```
   This activates the project environment, installs required Julia packages, validates the netMHCpan path, and checks that Perl is available.

## Input Data

CD8scape expects a data folder containing the following files:

- **alleles.txt**: List of HLA alleles, one per line (e.g. `HLA-A03:01`). Required for `run`.
- **supertype_panel.csv**: CSV with columns `Allele` and `Frequency` (and optionally `Locus`). Required for `run_supertype`. Can also be placed in the data folder to override the project default.
- **Variant file** (one of the following):
  - `.vcf` or `.vcf.gz` file (standard variant call format).
  - `single_locus_trajectories.out` from [Samfire](https://github.com/cjri/samfire) (also matches any `single_locus_trajectories*.out`, or falls back to any `.out` file in the folder).
  - A `.aa` amino-acid variant file — used when `read` is run with `--aa` (see below).
- **Reading frame file** (one of the following):
  - `sequences.fasta` **and** `consensus.fa` for the NCBI path: `sequences.fasta` provides ORF definitions with coordinate headers (e.g. from [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/)), and `consensus.fa` provides the full reference genome from which reading frame subsequences are extracted.
  - `Reading_Frames.dat` from [Samfire](https://github.com/cjri/samfire).

CD8scape will automatically detect and use the appropriate files for variant and reading frame parsing. VCF files are tried first; if no VCF is found or parsing fails, Samfire trajectory parsing is attempted. When `--aa` is passed to `read`, amino-acid variants are read instead. For reading frames, the NCBI path (`sequences.fasta` + `consensus.fa`) is tried first, falling back to Samfire's `Reading_Frames.dat`.

### Amino-acid variant file format (`.aa`)

The `.aa` format specifies variants directly at the amino-acid level — for example, to replicate substitutions reported in the literature without underlying sequence data. Pass `--aa` to the `read` command to use this format.

Each variant is two lines:
```
<orf_name> <aa_position>
<ancestral_aa> <derived_aa>
```

- `orf_name` must exactly match a `Description` value in `frames.csv` (set by the reading frame file).
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

Canonical codons are always used for both ancestral and derived amino acids, regardless of the actual consensus sequence. This means substitutions can be forced even when the consensus at that position does not encode the specified ancestral amino acid (e.g. when replicating published results without access to the original sequence data). A warning is printed and the frames file is updated in place whenever a consensus mismatch is overridden.

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

For Samfire, see [Samfire GitHub](https://github.com/cjri/samfire) for details on generating `.out` and `.dat` files.

## Usage
All commands are run from the repository root:

### 1. Prepare Environment
```bash
./CD8scape.jl prep
```
Note: `prep` performs all dependency installation and environment setup. The other commands (`read`, `simulate`, `run`, `run_supertype`, `percentile`) do not install packages and assume the environment is already prepared.

### 2. Parse Input Data
```bash
./CD8scape.jl read <folder_path> [--aa] [--suffix <name>] [--latest|--no-latest]
```
Parses variants and reading frames from the data folder, producing `variants.csv` and `frames.csv`.
- `--aa`: read amino-acid variants from a `.aa` file instead of VCF or Samfire trajectories (see [Amino-acid variant file format](#amino-acid-variant-file-format-aa) above).

### 3. Simulate Input Data
```bash
./CD8scape.jl simulate <folder_path> [--n <count>] [--p <proportion>] [--seed <int>] [--suffix <name>] [--latest|--no-latest]
```
- Parses reading frames (writes `frames.csv`) and generates exhaustive simulated single-nucleotide variants per reading frame (writes `variants.csv`).
- Sampling options:
   - `--n <count>`: sample an absolute number of variants.
   - `--p <proportion>` (alias `--prop`): sample a proportion in (0,1).
   - If both `--n` and `--p` are provided, `--n` takes precedence.
- Defaults: `--n` defaults to `1000` and `--p` defaults to `0.1` when the flag is provided without a value; omitting both flags writes all variants. `--seed` sets RNG seed (default: `1320`).

### 4. Run Pipeline (Individual Genotype)
```bash
./CD8scape.jl run <folder_path> [--t <N|max>|--thread <N|max>] [--max-escape] [--verbose] [--suffix <name>] [--latest|--no-latest]
```
- Generates peptides, runs netMHCpan, parses output, calculates best ranks and fold changes.
- `--t`/`--thread`: max parallel chunks for netMHCpan (default: 1). Use `max` to use the safety cap.
- `--verbose`: preserve per-allele logs and temp files for debugging.
- `--max-escape`: compute the single allele from the panel showing the largest predicted escape for each variant. Adds two columns to `harmonic_mean_best_ranks.csv`:
  - `max_escape_allele`: the HLA allele with the highest log2 fold change (ancestral EL rank ≤ 2% required; `missing` if no allele shows genuine escape).
  - `max_escape_log2_fc`: the corresponding log2(EL_Rank_derived / EL_Rank_ancestral) for that allele.

### 5. Run Pipeline (Supertype Panel)
```bash
./CD8scape.jl run_supertype <folder_path> [--t <N|max>|--thread <N|max>] [--max-escape] [--verbose] [--suffix <name>] [--latest|--no-latest]
```
- As above, but uses a representative supertype HLA panel.
- `--max-escape` is available but interpret results with caution: panel alleles are population-frequency surrogates rather than an individual's genotype, so the `max_escape_allele` result may not be biologically meaningful. A warning is printed when this flag is used with `run_supertype`.

### 6. Compute Percentiles (Benchmarking)
```bash
./CD8scape.jl percentile <folder_path> [--s <sim_file>] [--o <obs_file>]
```
- Computes observed HMBR log2 fold-change percentiles relative to a simulated distribution.
- `--s <sim_file>`: path to the simulated HMBR file (defaults to `harmonic_mean_best_ranks_simulated.csv` in the data folder).
- `--o <obs_file>`: path to the observed HMBR file (defaults to the most recent `harmonic_mean_best_ranks*.csv` in the data folder, excluding `_simulated`).
- Observed variants are excluded from the simulated distribution before computing percentiles.
- Writes `percentile_harmonic_mean_best_ranks.csv` (with a `Percentile` column, 0–100) to the data folder.

### Global Options

The following options are shared across `read`, `simulate`, `run`, and `run_supertype`:

- **`--suffix <name>`**: Append `_<name>` before the file extension of all output files. For example, `--suffix foo` produces `variants_foo.csv`, `best_ranks_foo.csv`, `harmonic_mean_best_ranks_foo.csv`, etc. For `simulate`, the suffix defaults to `simulated` when omitted.
- **`--latest`** (default) / **`--no-latest`**: Controls how input files are resolved when no `--suffix` is given and multiple candidates exist (e.g. `frames.csv` and `frames_simulated.csv`). `--latest` picks the most recently modified file; `--no-latest` raises an error on ambiguity.

These options allow multiple independent analyses (e.g. observed vs. simulated) to coexist in the same data folder without overwriting each other.

## Workflow Summary
1. **prep**: Install dependencies
2. **read**: Parse variants and frames from real data
3. **simulate**: Generate simulated single-nucleotide variants from reading frames
4. **run/run_supertype**: Generate peptides, predict binding, process output, calculate best ranks and fold changes
5. **percentile**: Compare observed fold changes to the simulated distribution

## Output Files
- `variants.csv`, `frames.csv`: Parsed input data
- `Peptides.pep`, `peptides_labels.csv`: Generated peptides and labels
- `netMHCpan_output.tsv`, `processed_output.csv`: Raw and processed netMHCpan results
- `best_ranks.csv`: Per-allele best EL ranks for ancestral and derived peptides at each locus.
- `harmonic_mean_best_ranks.csv`: Harmonic mean best ranks (HMBR) and log2 fold changes across the panel. With `--max-escape`, also includes `max_escape_allele` and `max_escape_log2_fc`.
- `variants_simulated.csv`, `harmonic_mean_best_ranks_simulated.csv`: Simulated variant data and HMBR results
- `percentile_harmonic_mean_best_ranks.csv`: Observed HMBR with percentile relative to simulated distribution

## Citation
If you use CD8scape in your research, please cite the repository and netMHCpan as appropriate.

Example data is from:
Stanevich, O.V., Alekseeva, E.I., Sergeeva, M. et al. SARS-CoV-2 escape from cytotoxic T cells during long-term COVID-19. Nat Commun 14, 149 (2023). https://doi.org/10.1038/s41467-022-34033-x


## License
GNU Public License