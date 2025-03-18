## CD8scape
Runs netMHCpan for variants. 

# Requirements:
- perl 5
- Julia v1.11
- R v4.4.1
- [netMHCpan4.1](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/) installation. A path to this should be added to [settings.txt](src/settings.txt)
- alleles.txt file present in the directory. This should be a list of the desrired HLA alleles for which to run netMHCpan. Each allele should be on a new line, e.g. 
    HLA-A01:01
    HLA-A02:01

# Methods:
- ./CD8scape.jl prep
    1) Installs Julia and R dependencies. 
- ./CD8scape.jl read $folder_path
    1) Looks for a [Samfire](https://github.com/cjri/samfire) single_locus_trajectories.out file in the folder (or, failing this, any .out file) from which to read variants. If this is missing, looks for a .vcf or .vcf.gz file from which to read variants. This is output to a file called variants.csv.
    2) Looks for a [Samfire](https://github.com/cjri/samfire) Reading_Frames.dat file in the folder from which to read open reading frames. If this is missing, looks for a reference reading frames file from [NCBI](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/) (this should be called sequences.fasta). This is output to a file called frames.csv.
- ./CD8scape.jl run $folder_path
    1) Generates consensus and variant peptides of lengths 8, 9, 10, and 11 for each variant locus. This will produce Peptides.pep, which is input to later steps, and peptides_labels.csv, which contains arbitrary labels used to later join consensus-variant peptide pairs. 
    2) Runs netMHCpan's MHC prediction model to generate eluted ligand (EL) scores and associated % ranks for each peptide. This will produce netMHCpan_output.tsv. 
    3) Runs process_output.pl to convert netMHCpan_output.tsv to a more readable processed_output.csv. 
    4) Finds the best rank for each provided HLA type and saves these as best_ranks.csv
    5) Calculates the harmonic mean best rank for each locus across loci for both ancestral and derived states, and then calculates the fold change between these so that a fold change of >1 represents a predicted weakning of presentation. 