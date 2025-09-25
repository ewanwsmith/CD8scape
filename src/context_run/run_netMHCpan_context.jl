#!/usr/bin/env julia
"""
This script runs NetMHCpan for context peptides and manages input/output files.

Purpose:
- Parses arguments and finds NetMHCpan executable.
- Runs NetMHCpan on cleaned context peptides.
- Handles output file management and error checking.

Usage:
    julia run_netMHCpan_context.jl --folder <path_to_folder> --mode panel|supertype

Arguments:
    --folder   Path to the folder containing input files.
    --mode     Run mode: panel or supertype.
"""

using CSV, DataFrames

function parse_arguments()
    args = Dict()
    for (i, arg) in enumerate(ARGS)
        if arg in ["--folder", "-f"]
            args["folder"] = ARGS[i + 1]
        elseif arg == "--mode"
            args["mode"] = ARGS[i + 1]
        end
    end
    if !haskey(args, "folder") || !haskey(args, "mode")
        error("Usage: ./run_netMHCpan_context.jl --folder /path/to/data --mode panel|supertype")
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
    if haskey(ENV, "NETMHCPAN")
        p = ENV["NETMHCPAN"] |> normpath
        println("NetMHCpan path (ENV) -> ", p)
        if isfile(p) && isexecutable(p)
            return p
        else
            error("NETMHCPAN in environment points to a non-existent or non-executable file: '$p'")
        end
    end

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
                if isfile(p) && isexecutable(p)
                    return p
                else
                    error("Path '$p' does not exist or is not executable.")
                end
            end
        else
            p = normpath(replace(s, "~" => homedir()))
            println("NetMHCpan path (settings) -> ", p)
            if isfile(p) && isexecutable(p)
                return p
            else
                error("Path '$p' does not exist or is not executable.")
            end
        end
    end
    error("No valid NETMHCPAN file path found in settings.txt.")
end

function run_netmhcpan(cmd, xlsfile_path)
    println("Running NetMHCpan with command:")
    println(cmd)

    try
        run(cmd)
        if isfile(xlsfile_path)
            println("NetMHCpan output found at $xlsfile_path")
        else
            error("ERROR: NetMHCpan did not create the expected output file.")
        end
        println("NetMHCpan completed successfully.")
    catch e
        println("Error running NetMHCpan: ", e)
        exit(1)
    end
end

function main()
    args = parse_arguments()
    folder = args["folder"]
    mode = args["mode"]
    netmhcpan = get_netMHCpan_path()
    peptides_file = joinpath(folder, "context_peptides.pep")
    if !isfile(peptides_file)
        error("File $peptides_file not found. Please ensure it exists.")
    end
    xlsfile_path = joinpath(folder, "netMHCpan_output.tsv")

    testfile = joinpath(folder, "test_write.txt")
    try
        open(testfile, "w") do io
            write(io, "test")
        end
        rm(testfile)
    catch
        error("ERROR: Cannot write to output folder. Check permissions!")
    end

    # Build allele format map from MHC_pseudo.dat
    settings_file = find_settings_file()
    pseudo_file = joinpath(dirname(settings_file), "MHC_pseudo.dat")
    allele_format_map = Dict{String,String}()
    if isfile(pseudo_file)
        for line in readlines(pseudo_file)
            s = strip(line)
            if isempty(s) || startswith(s, "#")
                continue
            end
            allele = split(s)[1]
            key = replace(replace(allele, "*" => ""), ":" => "")
            allele_format_map[key] = allele
        end
    else
        println("Warning: MHC_pseudo.dat not found at $pseudo_file. Allele format resolution will be skipped.")
    end

    if mode == "panel"
        alleles_file = joinpath(folder, "alleles.txt")
        if !isfile(alleles_file)
            error("File $alleles_file not found.")
        end
        raw_allele_list = open(alleles_file) do file
            [replace(split(line, r"\s+")[1], "*" => "") for line in readlines(file) if !isempty(line)]
        end
        allele_list = String[]
        skipped_alleles = String[]
        for a in raw_allele_list
            key = replace(replace(a, "*" => ""), ":" => "")
            if haskey(allele_format_map, key)
                push!(allele_list, allele_format_map[key])
            else
                push!(skipped_alleles, a)
            end
        end
        if !isempty(skipped_alleles)
            println("Warning: The following alleles were not found in MHC_pseudo.dat and will be skipped:")
            println(join(skipped_alleles, ", "))
        end
        alleles = join(allele_list, ",")


    elseif mode == "supertype"
        function find_representatives_file(folder_path::AbstractString)
            local_candidate = abspath(joinpath(folder_path, "supertype_panel.csv"))
            if isfile(local_candidate)
                return local_candidate
            end
            cwd_candidate = abspath("supertype_panel.csv")
            if isfile(cwd_candidate)
                return cwd_candidate
            end
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
            error("supertype_panel.csv not found.")
        end

        representatives_file = find_representatives_file(folder)
        squished_panel_path = joinpath(folder, "supertype_panel_squished.csv")
        squishing_map_path = joinpath(folder, "squishing_map.csv")
        df = CSV.read(representatives_file, DataFrame; normalizenames=true)

        function _clean_sym(n)
            str = String(n)
            str = replace(str, r"[\s\u00A0]+" => "")
            return Symbol(lowercase(str))
        end
        rename!(df, Dict(n => _clean_sym(n) for n in names(df)))
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
        freq_col = _find_col(df, "frequency")

        # Robust allele normalization
        function normalize_allele(s::AbstractString)
            s2 = String(s)
            s2 = replace(s2, '\u00A0' => ' ')
            s2 = replace(s2, r"\s+" => "")
            s2 = uppercase(s2)
            m = match(r"^(HLA-[A-Z]+)\*(\d{2}):(\d{2})", s2)
            if m !== nothing
                return string(m.captures[1], "*", m.captures[2], ":", m.captures[3])
            end
            m = match(r"^(HLA-[A-Z]+)\*(\d{4})$", s2)
            if m !== nothing
                g = m.captures[2]
                return string(m.captures[1], "*", g[1:2], ":", g[3:4])
            end
            m = match(r"^(HLA-[A-Z]+)(\d{2}):(\d{2})$", s2)
            if m !== nothing
                return string(m.captures[1], "*", m.captures[2], ":", m.captures[3])
            end
            m = match(r"^(HLA-[A-Z]+)(\d{4})$", s2)
            if m !== nothing
                g = m.captures[2]
                return string(m.captures[1], "*", g[1:2], ":", g[3:4])
            end
            return s2
        end
        df[!, allele_col] = map(a -> a === missing ? missing : normalize_allele(a), df[!, allele_col])

        # --- Neighbour squishing logic with normalization ---
        neighbours_path = joinpath(@__DIR__, "..", "allele_neighbours.csv")
        squishing_map = DataFrame(Original=String[], Squished=String[])
        if isfile(neighbours_path)
            neighbours_df = CSV.read(neighbours_path, DataFrame)
            rename!(neighbours_df, Dict(n => _clean_sym(n) for n in names(neighbours_df)))
            norm_allele_map = Dict(normalize_allele(String(row.allele)) => row.neighbour for row in eachrow(neighbours_df))
            squished_alleles = String[]
            for a in df[!, allele_col]
                norm_a = normalize_allele(String(a))
                if haskey(norm_allele_map, norm_a)
                    squished = norm_allele_map[norm_a]
                    push!(squished_alleles, squished)
                    push!(squishing_map, (Original=String(a), Squished=squished))
                else
                    println("Warning: squishing step failed for allele: $a (no neighbour match)")
                    push!(squished_alleles, String(a))
                    push!(squishing_map, (Original=String(a), Squished=String(a)))
                end
            end
            df[!, allele_col] = squished_alleles
            CSV.write(squished_panel_path, df)
            println("Squished supertype panel written to $squished_panel_path")
            CSV.write(squishing_map_path, squishing_map)
            println("Allele squishing map written to $squishing_map_path")
        else
            println("Warning: allele_neighbours.csv not found, skipping neighbour squishing.")
        end

        # Build allele list for netMHCpan (remove '*' as expected by downstream)
        raw_allele_list = [replace(String(a), "*" => "") for a in collect(skipmissing(df[!, allele_col])) if !isempty(strip(String(a)))]
        allele_list = String[]
        skipped_alleles = String[]
        for a in raw_allele_list
            key = replace(replace(a, "*" => ""), ":" => "")
            if haskey(allele_format_map, key)
                push!(allele_list, allele_format_map[key])
            else
                push!(skipped_alleles, a)
            end
        end
        if !isempty(skipped_alleles)
            println("Warning: The following alleles were not found in MHC_pseudo.dat and will be skipped:")
            println(join(skipped_alleles, ", "))
        end
        alleles = join(allele_list, ",")

    else
        error("Unsupported mode: $mode. Use panel or supertype.")
    end

    # Chunk alleles by total character length (≤1024 chars for -a argument)
    max_chars = 1024
    allele_chunks = Vector{Vector{String}}()
    current_chunk = String[]
    current_len = 0
    for allele in allele_list
        add_len = length(allele) + (isempty(current_chunk) ? 0 : 1)
        if current_len + add_len > max_chars
            push!(allele_chunks, current_chunk)
            current_chunk = String[]
            current_len = 0
        end
        push!(current_chunk, allele)
        current_len += add_len
    end
    if !isempty(current_chunk)
        push!(allele_chunks, current_chunk)
    end
    temp_files = String[]

    for (i, chunk) in enumerate(allele_chunks)
        if isempty(chunk)
            println("Skipping empty chunk $(i)")
            continue
        end
        alleles_str = join(chunk, ",")
        chunk_out = joinpath(folder, "netMHCpan_output_chunk$(i).tsv")
        cmd = `$netmhcpan -p $peptides_file -xls -a $alleles_str -xlsfile $chunk_out`
        run_netmhcpan(cmd, chunk_out)
        if isfile(chunk_out)
            push!(temp_files, chunk_out)
        end
    end

    # Join chunk outputs
    println("Joining NetMHCpan outputs...")
    open(xlsfile_path, "w") do out_io
        for (i, temp_file) in enumerate(temp_files)
            open(temp_file, "r") do in_io
                for (j, line) in enumerate(eachline(in_io))
                    # Write header only for first chunk
                    if i == 1 || j > 1
                        write(out_io, line, "\n")
                    end
                end
            end
            rm(temp_file)
        end
    end
    println("NetMHCpan outputs joined to $xlsfile_path")
end

main()