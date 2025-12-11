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
   - `--p <proportion>` (alias `--prop`): sample a proportion in (0,1].
   - If both `--n` and `--p` are provided, `--n` takes precedence.
 - Defaults: `--n` defaults to `1000` and `--p` defaults to `0.1` when you choose to sample by count or proportion; no sampling flags writes all variants. `--seed` sets RNG seed (default: `1320`).

### 4. Run Pipeline (Individual Genotype)
```bash
./CD8scape.jl run <folder_path>
```
- Generates peptides, runs netMHCpan, parses output, calculates best ranks and fold changes.

### 5. Run Pipeline (Supertype Panel)
```bash
./CD8scape.jl run_supertype <folder_path>
```
- As above, but uses a representative supertype HLA panel.

 

## Workflow Summary
1. **prep**: Install dependencies
2. **read**: Parse variants and frames
3. **run/run_supertype**: Generate peptides, predict binding, process output, calculate best ranks and fold changes
 

## Output Files
- `variants.csv`, `frames.csv`: Parsed input data
- `Peptides.pep`, `peptides_labels.csv`: Generated peptides and labels
- `netMHCpan_output.tsv`, `processed_output.csv`: Raw and processed netMHCpan results
- `best_ranks.csv`, `harmonic_mean_best_ranks.csv`: Best ranks and fold change analysis
 

## Citation

If you use CD8scape in your research, please cite the repository and netMHCpan as appropriate.

Example data is from:

Stanevich, O.V., Alekseeva, E.I., Sergeeva, M. et al. SARS-CoV-2 escape from cytotoxic T cells during long-term COVID-19. Nat Commun 14, 149 (2023). https://doi.org/10.1038/s41467-022-34033-x

## License
GNU Public License