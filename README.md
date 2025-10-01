# CD8scape

CD8scape runs netMHCpan on genetic variants for individual HLA genotypes or representative supertype panels.

## Features
- Automated peptide generation for consensus and variant loci
- MHC binding prediction using netMHCpan
- Robust output parsing and best-rank calculation
- Harmonic mean best rank (HMBR) and fold change analysis
- Context run for placing output metrics in the context of a simulated random background

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

- **alleles.txt**: List of HLA alleles, one per line (e.g., `HLA-A01:01`). File name is case-insensitive.
- **Variant file**: Either
  - `single_locus_trajectories.out` (from [Samfire](https://github.com/cjri/samfire)), or
  - Any `.out` file with SAMFIRE format, or
  - `.vcf` or `.vcf.gz` file (standard variant call format).
- **Reading frame file**: Either
  - `Reading_Frames.dat` (from [Samfire](https://github.com/cjri/samfire)), or
  - `sequences.fa/fasta` (reference genome from [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/)).
- **consensus.fa/fasta** (required only with NCBI sequences): Reference consensus sequence file containing the full genome sequence. This file is needed when using `sequences.fa/fasta` from NCBI, as it provides the complete reference sequence from which reading frame regions are extracted based on coordinates specified in the sequences file headers.

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
    alleles.txt                     # HLA alleles (case-insensitive filename)
    single_locus_trajectories.out   # or any .out file with SAMFIRE format, or .vcf/.vcf.gz
    Reading_Frames.dat              # or sequences.fa/fasta + consensus.fa/fasta
```

**Note**: When using NCBI format with `sequences.fa/fasta`, you must also provide `consensus.fa/fasta` containing the full reference genome. The sequences file contains coordinate ranges that are extracted from the consensus sequence.

CD8scape uses case-insensitive file discovery, so files can be named with any capitalization (e.g., `alleles.txt`, `Alleles.txt`, `ALLELES.TXT` are all acceptable). Both `.fa` and `.fasta` extensions are supported for sequence files.

CD8scape will automatically detect and use the appropriate files for variant and reading frame parsing. For Samfire, see [Samfire GitHub](https://github.com/cjri/samfire) for details on generating `.out` and `.dat` files.

### Reading Frame Formats

CD8scape supports two reading frame input formats:

1. **Samfire format**: Uses `Reading_Frames.dat` containing pre-extracted reading frame sequences. This is a self-contained format where each reading frame is already provided as a sequence.

2. **NCBI format**: Uses `sequences.fa/fasta` + `consensus.fa/fasta` where:
   - `sequences.fa/fasta` contains headers with coordinate ranges (e.g., `>accession:77..496`)
   - `consensus.fa/fasta` contains the full reference genome sequence
   - CD8scape extracts the specified coordinate ranges from the consensus sequence to generate reading frames

The NCBI format is useful when working with reference genomes from databases like NCBI Virus, where you have coordinate information and need to extract specific regions from a full genome sequence.

## Usage
All commands are run from the repository root:

### 1. Prepare Environment
```bash
./CD8scape.jl prep
```

### 2. Parse Input Data
```bash
./CD8scape.jl read <folder_path>
```
- Parses variants and reading frames, outputs `variants.csv` and `frames.csv`.

### 3. Run Pipeline (Individual Genotype)
```bash
./CD8scape.jl run <folder_path>
```
- Generates peptides, runs netMHCpan, parses output, calculates best ranks and fold changes.

### 4. Run Pipeline (Supertype Panel)
```bash
./CD8scape.jl run_supertype <folder_path>
```
- As above, but uses a representative supertype HLA panel.

### 5. Context Run (Benchmarking)
```bash
./CD8scape.jl context <folder_path> [--supertype] [--n_loci <number_of_loci>] [--seed <random_seed>] [--force]
```
- Runs context-sensitive pipeline to generate simulated loci and peptides, runs netMHCpan, and compares observed fold changes to a distribution of random mutations for the same genome.
- `--n_loci <number_of_loci>` sets the number of simulated loci (default: 1000).
- `--seed <random_seed>` sets the random number seed for reproducible loci selection (default: 1320).
- If `--supertype` is provided, uses the representative supertype HLA panel for predictions.
- If intermediate results exist, resumes from the appropriate step.
- Use `--force` to rerun the full pipeline.

## Workflow Summary
1. **prep**: Install dependencies
2. **read**: Parse variants and frames
3. **run/run_supertype**: Generate peptides, predict binding, process output, calculate best ranks and fold changes
4. **context**: Simulate background, run predictions, compare observed results to context distribution

## Output Files
- `variants.csv`, `frames.csv`: Parsed input data
- `Peptides.pep`, `peptides_labels.csv`: Generated peptides and labels  
- `netMHCpan_output.tsv`, `processed_output.csv`: Raw and processed netMHCpan results
- `processed_peptides.csv`: Processed peptides with scores
- `best_ranks.csv`, `harmonic_mean_best_ranks.csv`: Best ranks and fold change analysis
- Context run outputs:
  - `context_peptides.pep`: Generated context peptides
  - `context_processed_netMHCpan_output.csv`: Processed context results
  - `context_scores.csv`, `context_harmonic_mean_best_ranks.csv`: Context analysis results
  - `harmonic_mean_best_ranks_with_percentile.csv`: Final results with percentile comparisons

## Citation

If you use CD8scape in your research, please cite the repository and netMHCpan as appropriate.

Example data is from:

Stanevich, O.V., Alekseeva, E.I., Sergeeva, M. et al. SARS-CoV-2 escape from cytotoxic T cells during long-term COVID-19. Nat Commun 14, 149 (2023). https://doi.org/10.1038/s41467-022-34033-x

## License
GNU Public License