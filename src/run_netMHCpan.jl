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

# Function to get netMHCpan path from settings.txt
function get_netMHCpan_path()
    settings_file = joinpath(@__DIR__, "settings.txt")
    if !isfile(settings_file)
        error("settings.txt not found in the script directory. Please provide the file.")
    end

    # Read the first non-comment line
    open(settings_file) do file
        for line in eachline(file)
            line = strip(line)
            if !isempty(line) && !startswith(line, "#")
                return line
            end
        end
    end
    error("No valid netMHCpan path found in settings.txt.")
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

    # Get netMHCpan installation path
    netMHCpan_path = get_netMHCpan_path()

    # Change directory to netMHCpan installation path
    cd(netMHCpan_path)

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