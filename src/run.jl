# run_netMHCpan.jl

# Change directory to netMHCpan installation path
cd("/Users/e.smith.5/Documents/PhD/Software/netMHCpan-4.1")

# Path to the file containing alleles, one per line
alleles_file = "/Users/e.smith.5/Documents/PhD/CD8scape/data/RSV_example/alleles.txt"

# Read alleles from file and join them with commas
allele_list = readlines(alleles_file)
alleles = join(allele_list, ",")

# Define other variable options
binding_affinity_option = "-BA"
xlsfile_path = "/Users/e.smith.5/Documents/PhD/CD8scape/data/RSV_example/RSV.tsv"

# Construct the command with `-xls` directly in the command text
cmd = `./netMHCpan -p /Users/e.smith.5/Documents/PhD/CD8scape/data/RSV_example/Peptides.pep \
    $(binding_affinity_option) \
    -xls \
    -a $(alleles) \
    -xlsfile $(xlsfile_path)`

# Run the command
run(cmd)