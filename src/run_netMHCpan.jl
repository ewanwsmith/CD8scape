#!/usr/bin/env julia
using Serialization
using DataFrames
"""
run_netMHCpan.jl

Runs NetMHCpan for all peptides using provided alleles.

Usage:
    julia run_netMHCpan.jl --folder <folder_path>

Arguments:
    --folder   Path to the folder containing input files.
"""

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

function find_settings_file()
    current_dir = pwd()
    while current_dir != "/"
        candidate = joinpath(current_dir, "src", "settings.txt")
        if isfile(candidate)
            return candidate
        end
        current_dir = dirname(current_dir)
    end
    error("settings.txt not found in any 'src' directory above current working directory.")
end

function get_netMHCpan_executable()
    settings_file = find_settings_file()

    if !isfile(settings_file)
        error("settings.txt not found at $settings_file. Please provide the file.")
    end

    for line in readlines(settings_file)
        s = strip(line)
        if isempty(s) || startswith(s, "#")
            continue
        end
        if occursin('=', s)
            k, v = strip.(split(s, '=', limit=2))
            if uppercase(k) == "NETMHCPAN"
                p = normpath(replace(v, "~" => homedir()))
                println("NetMHCpan path (settings) -> ", p)
                if isfile(p)
                    return p
                else
                    error("ERROR: Path '$p' does not exist.")
                end
            end
        else
            p = normpath(replace(s, "~" => homedir()))
            println("NetMHCpan path (settings) -> ", p)
            if isfile(p)
                return p
            else
                error("ERROR: Path '$p' does not exist.")
            end
        end
    end

    error("ERROR: No valid netMHCpan path found in settings.txt.")
end

function main()
    args = parse_arguments()
    folder_path = args["folder"]

    netMHCpan_exe = get_netMHCpan_executable()
    println("Using NetMHCpan executable: ", netMHCpan_exe)

    # File paths
    alleles_file = joinpath(folder_path, "alleles.txt")
    peptides_file = joinpath(folder_path, "Peptides.pep")
    xlsfile_path = joinpath(folder_path, "netMHCpan_output.tsv")
    cache_file = joinpath(folder_path, "results_cache.jls")

    # Check if required files exist
    for file in [alleles_file, peptides_file]
        if !isfile(file)
            error("File $file not found. Please ensure it exists in the specified folder.")
        end
    end

    # Test write permissions
    test_file = joinpath(folder_path, "test_write.txt")
    try
        open(test_file, "w") do io
            write(io, "test")
        end
        rm(test_file)
    catch e
        error("ERROR: Cannot write to output folder. Check permissions!")
    end

    # Read alleles
    allele_list = open(alleles_file) do file
        [replace(split(line, r"\s+")[1], "*" => "") for line in readlines(file) if !isempty(line)]
    end

    # Read peptides
    peptide_list = open(peptides_file) do file
        [strip(line) for line in readlines(file) if !isempty(strip(line))]
    end

    # Load cache if exists, else initialize empty Dict
    cache = Dict{Tuple{String,String},Any}()
    if isfile(cache_file)
        try
            cache = deserialize(cache_file)
            println("Loaded cache with $(length(cache)) entries from $cache_file")
        catch e
            println("Warning: Could not load cache file, starting with empty cache. Error: $e")
            cache = Dict{Tuple{String,String},Any}()
        end
    else
        println("No cache file found, starting with empty cache.")
    end

    # TODO: Chunking and cache lookup logic will go here
end

main()
