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

    if mode == "panel"
        alleles_file = joinpath(folder, "alleles.txt")
        if !isfile(alleles_file)
            error("File $alleles_file not found.")
        end
        allele_list = open(alleles_file) do file
            [replace(split(line, r"\s+")[1], "*" => "") for line in readlines(file) if !isempty(line)]
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

        clean_allele = function(a)
            s = String(a)
            s = replace(s, r"\u00A0" => "")
            s = strip(s)
            s = replace(s, r"\*(\d{2})(\d{2})" => s"*\1:\2")
            return s
        end
        df[!, allele_col] = map(a -> a === missing ? missing : clean_allele(a), df[!, allele_col])
        allele_list = [replace(String(a), "*" => "") for a in collect(skipmissing(df[!, allele_col])) if !isempty(strip(String(a)))]
        alleles = join(allele_list, ",")

    else
        error("Unsupported mode: $mode. Use panel or supertype.")
    end

    cmd = `$netmhcpan -p $peptides_file -xls -a $alleles -xlsfile $xlsfile_path`
    run_netmhcpan(cmd, xlsfile_path)
end

main()