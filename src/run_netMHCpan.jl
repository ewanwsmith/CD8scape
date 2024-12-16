#!/usr/bin/env julia

# run_netMHCpan.jl

using Base

# Function to parse command-line arguments
function parse_arguments()
    args = Dict()
    for (i, arg) in enumerate(ARGS)
        if arg in ["--folder", "-f"]
            args["folder"] = ARGS[i + 1]
        end
    end
    if !haskey(args, "folder")
        error("Usage: ./run_netMHCpan.jl --folder /path/to/data")
    end
    return args
end

# Main function
function main()
    # Parse the arguments
    args = parse_arguments()
    folder_path = args["folder"]

    # Paths to required files in the specified folder
    alleles_file = joinpath(folder_path, "alleles.txt")
    peptides_file = joinpath(folder_path, "Peptides.pep")
    xlsfile_path = joinpath(folder_path, "netMHCpan_output.tsv")

    # Check that required files exist
    for file in [alleles_file, peptides_file]
        if !isfile(file)
            error("File $file not found. Please ensure it exists in the specified folder.")
        end
    end

    # Change directory to netMHCpan installation path
    cd("/Users/e.smith.5/Documents/PhD/Software/netMHCpan-4.1")

    # Read only the first column from the alleles file
    allele_list = open(alleles_file) do file
        [split(line, r"\s+")[1] for line in readlines(file) if !isempty(line)]
    end

    # Join the alleles with commas
    alleles = join(allele_list, ",")

    # Construct the netMHCpan command
    cmd = `./netMHCpan -p $(peptides_file) \
        -xls \
        -a $(alleles) \
        -xlsfile $(xlsfile_path)`

    # Run the command
    println("Running netMHCpan with command:")
    println(cmd)
    run(cmd)
end

main()