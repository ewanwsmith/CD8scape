#!/usr/bin/env julia
"""
Lightweight status logger. When `overwrite=true`, uses a carriage return
without newline to update the same line; otherwise prints a newline.
"""
function status(msg; overwrite=false)
    if overwrite
        print("\r", msg)
        flush(stdout)
    else
        println(msg)
    end
end
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
                if isfile(p)
                    return p
                else
                    error("ERROR: Path '$p' does not exist.")
                end
            end
        else
            p = normpath(replace(s, "~" => homedir()))
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

    # Chunking logic: run NetMHCpan in batches of 1000 peptides per chunk, per allele
    chunk_size = 1000
    total_peptides = length(peptide_list)
    peptide_chunks = [peptide_list[i:min(i+chunk_size-1, total_peptides)] for i in 1:chunk_size:total_peptides]
    temp_out_files = String[]
    total_chunks = length(peptide_chunks)
    total_alleles = length(allele_list)
    # Precompute chunk sizes and cumulative counts for accurate progress
    chunk_lengths = map(length, peptide_chunks)
    cumulative_lengths = cumsum([0; chunk_lengths[1:end-1]])
    for (chunk_idx, chunk) in enumerate(peptide_chunks)
        chunk_peps = chunk
        if isempty(chunk_peps)
            continue
        end
        temp_pep_file = joinpath(folder_path, "_temp_peptides_$(chunk_idx).pep")
        open(temp_pep_file, "w") do io
            for pep in chunk_peps
                println(io, pep)
            end
        end
        temp_out_files_chunk = String[]
        for (allele_idx, allele) in enumerate(allele_list)
            temp_out_file = joinpath(folder_path, "_temp_netMHCpan_output_$(chunk_idx)_$(allele_idx).tsv")
            push!(temp_out_files_chunk, temp_out_file)
            cmd = Cmd([netMHCpan_exe, "-p", temp_pep_file, "-xls", "-a", allele, "-xlsfile", temp_out_file])
            # Progress based on total peptides across alleles, accounting for variable chunk sizes
            total_work = total_alleles * total_peptides
            processed_before_chunk = total_alleles * cumulative_lengths[chunk_idx]
            processed_in_current = allele_idx * length(chunk_peps)
            percent_done = Int(floor(100 * (processed_before_chunk + processed_in_current) / total_work))
            status("Running chunk $(chunk_idx) / $(total_chunks). Chunk size: $(length(chunk_peps)) peptides. $(percent_done)% complete."; overwrite=true)
            try
                run(pipeline(cmd, stdout=devnull, stderr=devnull))
                if !isfile(temp_out_file)
                    error("ERROR: NetMHCpan did not create the expected output file for chunk/allele.")
                end
            catch e
                println("Error running NetMHCpan on chunk $(chunk_idx), allele $(allele): ", e)
                exit(1)
            end
        end
        append!(temp_out_files, temp_out_files_chunk)
    end
    # Merge all temp output files into final netMHCpan_output.tsv
    # Ensure the progress line ends with a newline before merging message
    print("\n")
    status("Merging chunk outputs into $xlsfile_path ...")
    open(xlsfile_path, "w") do out_io
        for (i, temp_file) in enumerate(temp_out_files)
            open(temp_file, "r") do in_io
                for (j, line) in enumerate(eachline(in_io))
                    # Write header only for first file
                    if i == 1 || (j > 1 && !startswith(line, "#"))
                        println(out_io, line)
                    end
                end
            end
        end
    end
    # Cleanup temp files
    status("Cleaning up temporary files...")
    for f in readdir(folder_path)
        if (startswith(f, "_temp_peptides") && endswith(f, ".pep")) || (startswith(f, "_temp_netMHCpan_output") && endswith(f, ".tsv"))
            temp_path = joinpath(folder_path, f)
            if isfile(temp_path)
                rm(temp_path; force=true)
            end
        end
    end
    println("Temporary files deleted.")
end

main()
