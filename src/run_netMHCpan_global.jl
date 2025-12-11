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
"""
run_netMHCpan_global.jl

Runs NetMHCpan for all peptides using a global supertype panel.

Usage:
    julia run_netMHCpan_global.jl --folder <folder_path>

Arguments:
    --folder   Path to the folder containing input files.
"""

include("env.jl")
using CSV, DataFrames
using Serialization

function parse_arguments()
    args = Dict()
    for (i, arg) in enumerate(ARGS)
        if arg in ["--folder", "-f"]
            args["folder"] = ARGS[i + 1]
        end
    end
    if !haskey(args, "folder")
        error("Usage: ./run_netMHCpan_global.jl --folder /path/to/data")
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

function get_netMHCpan_path()
    # 1) Prefer ENV set by env.jl
    if haskey(ENV, "NETMHCPAN")
        p = ENV["NETMHCPAN"] |> normpath
        # println("NetMHCpan path (ENV) -> ", p)  # Removed as per user request
        if isfile(p) && Base.Filesystem.isexecutable(p)
            return p
        else
            error("ERROR: NETMHCPAN in environment points to a non-existent or non-executable file: '$p'")
        end
    end
    # 2) Otherwise read from settings.txt (supports NETMHCPAN=/full/path or bare path)
    settings_file = find_settings_file()
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
                if isfile(p) && Base.Filesystem.isexecutable(p)
                    return p
                else
                    error("ERROR: Path '$p' does not exist or is not executable.")
                end
            end
        else
            p = normpath(replace(s, "~" => homedir()))
            println("NetMHCpan path (settings) -> ", p)
            if isfile(p) && Base.Filesystem.isexecutable(p)
                return p
            else
                error("ERROR: Path '$p' does not exist or is not executable.")
            end
        end
    end
    error("ERROR: No valid NETMHCPAN file path found in settings.txt.")
end

function main()
    args = parse_arguments()
    folder_path = args["folder"]

    netMHCpan_path = get_netMHCpan_path()
    println("Using NetMHCpan from: ", netMHCpan_path)

    # File paths
    function find_representatives_file(folder_path::AbstractString)
        # 1) Prefer a supertype_panel.csv inside the dataset folder passed via --folder
        local_candidate = abspath(joinpath(folder_path, "supertype_panel.csv"))
        if isfile(local_candidate)
            return local_candidate
        end
        # 2) Fall back to current working directory (if someone runs the script from the dataset folder)
        cwd_candidate = abspath("supertype_panel.csv")
        if isfile(cwd_candidate)
            return cwd_candidate
        end
        # 3) Ascend the directory tree to find src/supertype_panel.csv (project default)
        dir = @__DIR__
        while true
            candidate = joinpath(dir, "src", "supertype_panel.csv")
            if isfile(candidate)
                return candidate
            end
            parent = dirname(dir)
            if parent == dir
                break
            end
            dir = parent
        end
        error("supertype_panel.csv not found. Looked in: $(local_candidate), $(cwd_candidate), and project src/ upwards from $(abspath(@__DIR__)).")
    end
    representatives_file = find_representatives_file(folder_path)
    peptides_file = joinpath(folder_path, "Peptides.pep")
    xlsfile_path = joinpath(folder_path, "netMHCpan_output.tsv")
    cache_file = joinpath(folder_path, "results_cache.jls")

    # Check if required files exist
    if !isfile(peptides_file)
        error("File $peptides_file not found. Please ensure it exists in the specified folder.")
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

    # Read peptides
    peptide_list = open(peptides_file) do file
        [strip(line) for line in readlines(file) if !isempty(strip(line))]
    end


    # TODO: Chunking and cache lookup logic will go here

    # Read alleles from the supertype_panel.csv file (robust to header quirks)
    df = CSV.read(representatives_file, DataFrame; normalizenames=true)
    
    # Normalize column names: strip spaces/NBSP and lowercase
    function _clean_sym(n)
        # Accept Symbol or String column names, normalize to lowercase Symbol
        str = n isa Symbol ? String(n) : String(n)
        # remove regular spaces and non-breaking spaces
        str = replace(str, r"[\s\u00A0]+" => "")
        return Symbol(lowercase(str))
    end
    rename!(df, Dict(n => _clean_sym(n) for n in names(df)))
    
    # Discover the allele and (optional) frequency columns robustly
    function _find_col(df::DataFrame, target::String)
        for n in names(df)
            cleaned = lowercase(replace(String(n), r"[\s\u00A0]+" => ""))
            if cleaned == target
                return n
            end
        end
        return nothing
    end
    allele_col = _find_col(df, "allele")
    if allele_col === nothing
        cols = join(string.(names(df)), ", ")
        error("ERROR: Column 'Allele' not found in $(representatives_file). Columns present: $cols")
    end
    freq_col = _find_col(df, "frequency")
    
    # Optional Frequency filtering if present (accepts string or number)
    if freq_col !== nothing
        df[!, freq_col] = map(x -> try
                x === missing ? missing :
                x isa Number ? float(x) :
                parse(Float64, replace(String(x), r"[\s\u00A0]" => ""))
            catch
                missing
            end, df[!, freq_col])
        df = filter(row -> row[freq_col] !== missing && row[freq_col] != 0.0, df)
    end
    
    # Clean allele strings in the discovered allele column
    clean_allele = function(a)
        s = String(a)
        s = replace(s, r"\u00A0" => "")              # non-breaking spaces
        s = strip(s)
        # Insert colon if four consecutive digits after '*'
        s = replace(s, r"\*(\d{2})(\d{2})" => s"*\1:\2")
        return s
    end
    df[!, allele_col] = map(a -> a === missing ? missing : clean_allele(a), df[!, allele_col])
    
    # Build allele list for netMHCpan (remove '*' as expected by downstream)
    allele_list = [replace(String(a), "*" => "") for a in collect(skipmissing(df[!, allele_col])) if !isempty(strip(String(a)))]
    alleles = join(allele_list, ",")
    netmhcpan_exec = netMHCpan_path

    # Prepare to chunk peptides and check cache

    chunk_size = 1000
    total_peptides = length(peptide_list)
    peptide_chunks = [peptide_list[i:min(i+chunk_size-1, total_peptides)] for i in 1:chunk_size:total_peptides]
    # No cache logic
    # For each chunk, run NetMHCpan only on uncached peptides
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
        # For each allele, run NetMHCpan and update progress
        temp_out_files_chunk = String[]
        for (allele_idx, allele) in enumerate(allele_list)
            temp_out_file = joinpath(folder_path, "_temp_netMHCpan_output_$(chunk_idx)_$(allele_idx).tsv")
            push!(temp_out_files_chunk, temp_out_file)
            cmd = Cmd([netmhcpan_exec, "-p", temp_pep_file, "-xls", "-a", allele, "-xlsfile", temp_out_file])
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
            # Parse NetMHCpan output and update cache (wide format)
            # (Parsing handled below when merging outputs)
        end
        append!(temp_out_files, temp_out_files_chunk)
    end
    # Merge all temp output files into final netmhcpan_output.tsv
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
    status("Merged output written to $xlsfile_path")
    # Cleanup temp files
    status("Cleaning up temporary files...")
    # Remove all temp peptides and NetMHCpan output files matching patterns (including extra underscores/numbers)
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
