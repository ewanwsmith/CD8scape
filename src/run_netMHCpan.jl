#!/usr/bin/env julia

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

# Function to get netMHCpan path
function get_netMHCpan_path()
    settings_file = "/Users/e.smith.5/Documents/PhD/CD8scape/src/settings.txt"

    if !isfile(settings_file)
        error("settings.txt not found at $settings_file. Please provide the file.")
    end

    # Read the first non-comment, non-empty line and clean it
    for line in readlines(settings_file)
        line = strip(line)
        if !isempty(line) && !startswith(line, "#")
            path = replace(line, r"['\"]" => "") |> normpath
            println("Failed to parse path to netMHCpan installation : netMHCpan path -> ", path)

            if isdir(path)
                return path  # Return immediately if valid
            else
                error("ERROR: Path '$path' does not exist.")
            end
        end
    end

    # Error if no valid lines are found
    error("ERROR: No valid netMHCpan path found in settings.txt.")
end

# Main function
function main()
    # Parse arguments
    args = parse_arguments()
    folder_path = args["folder"]

    # Get netMHCpan installation path
    netMHCpan_path = get_netMHCpan_path()

    # Change to the netMHCpan directory
    println("Changing directory to -> ", netMHCpan_path)
    cd(netMHCpan_path)
    println("Current working directory -> ", pwd())

    # Paths to required files
    alleles_file = joinpath(folder_path, "alleles.txt")
    peptides_file = joinpath(folder_path, "Peptides.pep")
    xlsfile_path = joinpath(folder_path, "netMHCpan_output.tsv")

    # Check required files
    for file in [alleles_file, peptides_file]
        if !isfile(file)
            error("File $file not found. Please ensure it exists in the specified folder.")
        end
    end

    # Read alleles
    allele_list = open(alleles_file) do file
        [split(line, r"\s+")[1] for line in readlines(file) if !isempty(line)]
    end
    alleles = join(allele_list, ",")

    # Construct and run netMHCpan command
    cmd = `./netMHCpan -p $(peptides_file) -xls -a $(alleles) -xlsfile $(xlsfile_path)`
    println("Running netMHCpan with command:")
    println(cmd)
    run(cmd)
end

main()