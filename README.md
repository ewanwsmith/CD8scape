# CD8scape

CD8scape is a pipeline for running netMHCpan on genetic variants, supporting both individual HLA genotypes and supertype panels, with context-sensitive analysis for benchmarking.

## Features
- Automated peptide generation for consensus and variant loci
- MHC binding prediction using netMHCpan
- Robust output parsing and best-rank calculation
- Harmonic mean best rank (HMBR) and fold change analysis
- Context run for benchmarking observed results against simulated background

## Requirements
- Perl 5
- Julia v1.11+
- `alleles.txt` file listing HLA alleles (one per line, e.g. `HLA-A01:01`)

## Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/ewanwsmith/CD8scape.git
   cd CD8scape
   ```
2. Download and install netMHCpan, then set its path in `src/settings.txt`.
   - [netMHCpan 4.1 download](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/)
   - Add a line to `src/settings.txt`: `NETMHCPAN=/full/path/to/netMHCpan`
3. Install Julia and Perl dependencies:
   ```bash
   ./CD8scape.jl prep
   ```

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
./CD8scape.jl context <folder_path> [--supertype] [--n_loci <number_of_loci>] [--force]
```
- Runs context-sensitive pipeline to generate simulated loci and peptides, runs netMHCpan, and benchmarks observed fold changes against a background distribution.
- `--n_loci <number_of_loci>` sets the number of simulated loci (default: 1000).
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
- `best_ranks.csv`, `harmonic_mean_best_ranks.csv`: Best ranks and fold change analysis
- `context_scores.csv`, `context_harmonic_mean_best_ranks.csv`, `harmonic_mean_best_ranks_with_percentile.csv`: Context run outputs

## Citation
If you use CD8scape in your research, please cite the repository and netMHCpan as appropriate.

## License
MIT License