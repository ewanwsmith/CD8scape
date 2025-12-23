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
include("path_utils.jl")
using CSV, DataFrames
using Serialization

function parse_arguments()
    args = Dict()
    for (i, arg) in enumerate(ARGS)
        if arg in ["--folder", "-f"]
            args["folder"] = ARGS[i + 1]
        elseif arg == "--verbose"
            args["verbose"] = true
        elseif arg == "--t" || arg == "--thread"
            args["threads_raw"] = ARGS[i + 1]
        elseif arg == "--suffix"
            if i + 1 <= length(ARGS) && !startswith(ARGS[i+1], "--")
                args["suffix"] = ARGS[i + 1]
            else
                args["suffix"] = ""
            end
        elseif arg == "--latest"
            args["latest"] = true
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
    verbose = get(args, "verbose", false)
    suffix  = get(args, "suffix", "")
    latest  = get(args, "latest", true)
    requested = get(args, "threads_raw", nothing)
    # Safety cap: default to half of logical CPUs unless overridden by ENV
    cpu_threads = try
        Sys.CPU_THREADS
    catch
        1
    end
    default_cap = max(1, Int(floor(cpu_threads / 2)))
    env_cap = try
        haskey(ENV, "CD8SCAPE_MAX_THREADS") ? max(1, parse(Int, ENV["CD8SCAPE_MAX_THREADS"])) : nothing
    catch
        nothing
    end
    cap = env_cap === nothing ? default_cap : env_cap
    threads = 1
    if requested !== nothing
        v = lowercase(String(requested))
        if v == "max" || v == "cap"
            threads = cap
        else
            try
                threads = max(1, parse(Int, String(requested)))
            catch
                println("Invalid value for --t/--thread: '$(requested)'. Falling back to 1.")
                threads = 1
            end
            if threads > cap
                println("Capping parallel chunks to $(cap) (requested $(threads); logical CPUs=$(cpu_threads)). Set ENV CD8SCAPE_MAX_THREADS to adjust.")
                threads = cap
            end
        end
    end

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
    # Input peptides: ignore suffix; discover with latest
    peptides_file = resolve_read(joinpath(folder_path, "Peptides.pep"); suffix="", latest=latest)
    # Use lowercase filename to match downstream Perl script expectation
    xlsfile_path = resolve_write(joinpath(folder_path, "netmhcpan_output.tsv"); suffix=suffix)
    cache_file   = resolve_write(joinpath(folder_path, "results_cache.jls"); suffix=suffix)

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

    # If all peptides were dropped upstream, note and continue; downstream will handle
    if isempty(peptide_list)
        println("No peptides remaining after filtering (synonymity/stop codons).")
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
    locus_col = _find_col(df, "locus")
    
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

    # Identify alleles with trailing IMGT/HLA suffix letters (N, L, S, C, A, Q, G, P, etc.)
    # Examples: HLA-A*02:01N, HLA-B*07:02G
    function has_suffix(a)
        s = String(a)
        return occursin(r"(:\d{2})(?:[:\d{2}]*)?[A-Z]$", s)
    end
    suffix_mask = map(a -> a === missing ? false : has_suffix(a), df[!, allele_col])
    removed_count = count(suffix_mask)
    if removed_count > 0
        println("Removing $(removed_count) alleles with suffix letters and reweighting panel")
        df = df[.!suffix_mask, :]
        # Reweight frequencies: if locus present, normalize per locus; otherwise global normalize
        if freq_col !== nothing
            if locus_col !== nothing
                # group by locus and renormalize frequency to sum to original per locus sums
                # First compute sums per locus before removal (from original df including removed rows)
                # We reconstruct original sums using current df plus removed rows frequencies if available.
                # Simpler: normalize so that within each locus the remaining frequencies sum to 1.0
                grp = groupby(df, locus_col)
                for g in grp
                    total = sum(skipmissing(g[!, freq_col]))
                    if total > 0
                        g[!, freq_col] .= g[!, freq_col] ./ total
                    end
                end
                df = combine(grp, names(df) .=> identity)
            else
                total = sum(skipmissing(df[!, freq_col]))
                if total > 0
                    df[!, freq_col] .= df[!, freq_col] ./ total
                end
            end
        end
    end
    
    # Filter out only NA/missing/empty entries; do not restrict to NetMHCpan lists
    valid_mask = map(a -> a !== missing && !isempty(strip(String(a))) && lowercase(strip(String(a))) ∉ ("na", "nan", "missing"), df[!, allele_col])
    dropped_invalid = count(.!valid_mask)
    if dropped_invalid > 0
        println("Dropping $(dropped_invalid) invalid/NA/empty allele entries before NetMHCpan run")
        df = df[valid_mask, :]
    end
    # Build allele list for netMHCpan (remove '*' as expected by downstream)
    allele_list = [replace(String(a), "*" => "") for a in collect(df[!, allele_col])]
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
    # Helper to run a single chunk (all alleles) and return its temp outputs
    function run_chunk(chunk_idx::Int, chunk_peps::Vector{<:AbstractString})
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
            cmd = Cmd([netmhcpan_exec, "-p", temp_pep_file, "-xls", "-a", allele, "-xlsfile", temp_out_file])
            # Progress based on total peptides across alleles, accounting for variable chunk sizes
            total_work = total_alleles * total_peptides
            processed_before_chunk = total_alleles * cumulative_lengths[chunk_idx]
            processed_in_current = allele_idx * length(chunk_peps)
            percent_done = Int(floor(100 * (processed_before_chunk + processed_in_current) / total_work))
            status("Running chunk $(chunk_idx) / $(total_chunks). Chunk size: $(length(chunk_peps)) peptides. $(percent_done)% complete."; overwrite=true)
            try
                if verbose
                    allele_log = joinpath(folder_path, "_temp_netMHCpan_log_$(chunk_idx)_$(allele_idx).txt")
                    open(allele_log, "w") do log_io
                        run(pipeline(cmd, stdout=log_io, stderr=log_io))
                    end
                else
                    run(pipeline(cmd, stdout=devnull, stderr=devnull))
                end
                if !isfile(temp_out_file)
                    error("ERROR: NetMHCpan did not create the expected output file for chunk/allele.")
                end
            catch e
                println("\nError running NetMHCpan on chunk $(chunk_idx), allele $(allele): ", e)
                println("Command attempted: ", join(cmd.exec, " "))
                exit(1)
            end
        end
        return temp_out_files_chunk
    end
    # Run chunks with bounded concurrency based on --t/--thread
    sem = Base.Semaphore(threads)
    results = Vector{Vector{String}}(undef, total_chunks)
    tasks = Task[]

    # Shared progress state for parallel reporting
    total_work = total_alleles * total_peptides
    progress_lock = ReentrantLock()
    active_chunks = Set{Int}()
    processed_ref = Ref(0)
    reporter_task = nothing
    if threads > 1 && total_work > 0
        reporter_task = @async begin
            last_len = 0
            while true
                sleep(0.5)
                running = Int[]
                pct = 100
                done = false
                lock(progress_lock) do
                    running = sort(collect(active_chunks))
                    pct = total_work == 0 ? 100 : Int(floor(100 * processed_ref[] / total_work))
                    done = (processed_ref[] >= total_work) || (isempty(running) && processed_ref[] > 0)
                end
                line = "Running chunks: " * (isempty(running) ? "-" : join(running, ", ")) * " of $(total_chunks) total. $(pct)% complete"
                print("\r", line)
                if length(line) < last_len
                    print(" "^(last_len - length(line)))
                end
                flush(stdout)
                last_len = length(line)
                if done
                    break
                end
            end
            println()
        end
    end

    for (chunk_idx, chunk) in enumerate(peptide_chunks)
        chunk_peps = chunk
        if isempty(chunk_peps)
            results[chunk_idx] = String[]
            continue
        end
        t = @async begin
            Base.acquire(sem)
            try
                if threads == 1
                    # Original per-chunk progress output
                    results[chunk_idx] = run_chunk(chunk_idx, chunk_peps)
                else
                    # Parallel mode: suppress per-chunk printing and update global progress
                    lock(progress_lock) do
                        push!(active_chunks, chunk_idx)
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
                        cmd = Cmd([netmhcpan_exec, "-p", temp_pep_file, "-xls", "-a", allele, "-xlsfile", temp_out_file])
                        try
                            if verbose
                                allele_log = joinpath(folder_path, "_temp_netMHCpan_log_$(chunk_idx)_$(allele_idx).txt")
                                open(allele_log, "w") do log_io
                                    run(pipeline(cmd, stdout=log_io, stderr=log_io))
                                end
                            else
                                run(pipeline(cmd, stdout=devnull, stderr=devnull))
                            end
                            if !isfile(temp_out_file)
                                error("ERROR: NetMHCpan did not create the expected output file for chunk/allele.")
                            end
                        catch e
                            println("\nError running NetMHCpan on chunk $(chunk_idx), allele $(allele): ", e)
                            println("Command attempted: ", join(cmd.exec, " "))
                            exit(1)
                        end
                        lock(progress_lock) do
                            processed_ref[] += length(chunk_peps)
                        end
                    end
                    results[chunk_idx] = temp_out_files_chunk
                end
            finally
                if threads > 1 && total_work > 0
                    lock(progress_lock) do
                        delete!(active_chunks, chunk_idx)
                    end
                end
                Base.release(sem)
            end
        end
        push!(tasks, t)
    end

    for t in tasks
        wait(t)
    end
    if reporter_task !== nothing
        wait(reporter_task)
    end
    for r in results
        append!(temp_out_files, r)
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
    # If merged file is empty (no outputs), write a minimal stub matching Perl expectations
    try
        sz = filesize(xlsfile_path)
        if sz == 0
            open(xlsfile_path, "w") do out_io
                # Allele line: tab-separated allele names (without '*')
                println(out_io, join(allele_list, "\t"))
                # Header line: Pos, Peptide, ID, then 4 columns per allele: core, icore, EL-score, EL_Rank
                header = String["Pos", "Peptide", "ID"]
                for _ in allele_list
                    append!(header, ["core", "icore", "EL-score", "EL_Rank"])
                end
                println(out_io, join(header, "\t"))
            end
        end
    catch
        # If filesize check fails, proceed without stub
    end
    status("Merged output written to $xlsfile_path")
    # Cleanup temp files
    status("Cleaning up temporary files...")
    # Remove temp files; when verbose, keep peptides and outputs, only delete logs.
    del_pep = 0; del_out = 0; del_logs = 0
    keep_pep = 0; keep_out = 0; keep_logs = 0
    kept_logs_list = String[]
    for f in readdir(folder_path)
        pep = (startswith(f, "_temp_peptides") && endswith(f, ".pep"))
        out = (startswith(f, "_temp_netMHCpan_output") && endswith(f, ".tsv"))
        log = (startswith(f, "_temp_netMHCpan_log") && endswith(f, ".txt"))
        temp_path = joinpath(folder_path, f)
        if pep
            if verbose
                keep_pep += 1
            elseif isfile(temp_path)
                try; rm(temp_path; force=true); del_pep += 1; catch; end
            end
        elseif out
            if verbose
                keep_out += 1
            elseif isfile(temp_path)
                try; rm(temp_path; force=true); del_out += 1; catch; end
            end
        elseif log
            if verbose
                keep_logs += 1
                push!(kept_logs_list, f)
            elseif isfile(temp_path)
                try; rm(temp_path; force=true); del_logs += 1; catch; end
            end
        end
    end
    if verbose
        shown = min(length(kept_logs_list), 10)
        remaining = length(kept_logs_list) - shown
        shown_list = kept_logs_list[1:shown]
        msg = "Cleanup summary (verbose): kept peptides=$(keep_pep), kept outputs=$(keep_out), kept logs=$(keep_logs). Deleted peptides=$(del_pep), outputs=$(del_out), logs=$(del_logs). Kept logs (showing $(shown)" * (remaining > 0 ? ", +$(remaining) more" : "") * "): "
        println(msg, join(shown_list, ", "))
    end
end
main()
