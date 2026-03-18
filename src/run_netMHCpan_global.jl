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
    cap = get_env_int("CD8SCAPE_MAX_THREADS", default_cap)
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
    catch
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

    # Build allele list for netMHCpan (remove '*') and gather frequencies
    raw_netmhcpan = [replace(String(a), "*" => "") for a in collect(df[!, allele_col])]
    raw_canonical  = [String(a) for a in collect(df[!, allele_col])]
    raw_freqs = freq_col !== nothing ?
        Float64[coalesce(x, 0.0) for x in df[!, freq_col]] :
        fill(1.0 / max(1, nrow(df)), nrow(df))

    # Allele squishing: merge alleles that are identical after * removal (i.e. the same netMHCpan call).
    # Frequencies are summed for merged alleles so downstream weighting remains correct.
    seen_idx  = Dict{String, Int}()
    allele_list      = String[]   # deduplicated netMHCpan strings (no *)
    allele_canonical = String[]   # canonical form with * (for squished panel CSV)
    allele_freq_sum  = Float64[]  # summed frequencies
    squish_originals = Dict{String, Vector{String}}()  # netMHCpan string → all canonical forms

    for (net_a, canon_a, f) in zip(raw_netmhcpan, raw_canonical, raw_freqs)
        if haskey(seen_idx, net_a)
            i = seen_idx[net_a]
            allele_freq_sum[i] += f
            push!(squish_originals[net_a], canon_a)
        else
            push!(allele_list, net_a)
            push!(allele_canonical, canon_a)
            push!(allele_freq_sum, f)
            seen_idx[net_a] = length(allele_list)
            squish_originals[net_a] = [canon_a]
        end
    end

    n_squished_groups = count(v -> length(v) > 1, values(squish_originals))
    if n_squished_groups > 0
        println("Allele squishing: $(n_squished_groups) group(s) of equivalent alleles merged:")
        for net_a in allele_list
            if length(squish_originals[net_a]) > 1
                idx = seen_idx[net_a]
                println("  $(net_a): $(join(squish_originals[net_a], ", ")) → combined frequency $(round(allele_freq_sum[idx]; digits=4))")
            end
        end
        println("Panel reduced from $(length(raw_netmhcpan)) to $(length(allele_list)) unique alleles for netMHCpan.")
    end

    # Save squished panel so process_best_ranks_supertype can use the merged frequencies
    squished_panel_path = joinpath(folder_path, "supertype_panel_squished.csv")
    CSV.write(squished_panel_path, DataFrame(:Allele => allele_canonical, :Frequency => allele_freq_sum))
    println("Squished panel written to $(squished_panel_path)")

    # Save squishing map: one row per original panel allele, showing what it was merged into
    squish_map_df = DataFrame(
        :Original_Allele  => raw_canonical,
        :Squished_Into    => [allele_canonical[seen_idx[net_a]] for net_a in raw_netmhcpan],
        :Was_Squished     => [length(squish_originals[net_a]) > 1 for net_a in raw_netmhcpan],
    )
    squish_map_path = joinpath(folder_path, "squishing_map.csv")
    CSV.write(squish_map_path, squish_map_df)
    println("Squishing map written to $(squish_map_path)")

    # netMHCpan 4.2 enforces two limits on the -a allele argument:
    #   1. A 1024-character PTYPE_LINE limit on the joined allele string.
    #   2. An internal array limit of 78 alleles per call (79+ causes a crash).
    # Build batches greedily so each batch respects both limits.
    allele_char_limit  = get_env_int("CD8SCAPE_ALLELE_CHAR_LIMIT", 1023)
    allele_count_limit = get_env_int("CD8SCAPE_ALLELE_COUNT_LIMIT", 75)
    allele_batches = Vector{Vector{String}}()
    let current = String[], cur_len = 0
        for allele in allele_list
            entry_len = length(allele) + (isempty(current) ? 0 : 1)   # +1 for comma
            if !isempty(current) && (cur_len + entry_len > allele_char_limit || length(current) >= allele_count_limit)
                push!(allele_batches, current)
                current = [allele]
                cur_len = length(allele)
            else
                push!(current, allele)
                cur_len += entry_len
            end
        end
        isempty(current) || push!(allele_batches, current)
    end
    n_batches = length(allele_batches)
    if n_batches > 1
        println("Splitting $(length(allele_list)) alleles into $(n_batches) batches (char limit: $(allele_char_limit), count limit: $(allele_count_limit)).")
    end

    # Chunk peptides
    chunk_size = 1000
    total_peptides = length(peptide_list)
    peptide_chunks = [peptide_list[i:min(i+chunk_size-1, total_peptides)] for i in 1:chunk_size:total_peptides]
    total_chunks = length(peptide_chunks)
    total_alleles = length(allele_list)
    total_units  = total_chunks * n_batches   # each unit = one netMHCpan call

    # Write all peptide chunk files upfront (fast sequential I/O).
    # netMHCpan only reads these files, so multiple concurrent tasks can share a peptide file.
    temp_pep_files = Vector{Union{String,Nothing}}(undef, total_chunks)
    for (chunk_idx, chunk_peps) in enumerate(peptide_chunks)
        if isempty(chunk_peps)
            temp_pep_files[chunk_idx] = nothing
            continue
        end
        f = joinpath(folder_path, "_temp_peptides_$(chunk_idx).pep")
        open(f, "w") do io
            for pep in chunk_peps; println(io, pep); end
        end
        temp_pep_files[chunk_idx] = f
    end

    # results_matrix[chunk_idx, batch_idx] = output file path (nothing if chunk was empty)
    results_matrix = Matrix{Union{String,Nothing}}(nothing, total_chunks, n_batches)

    # Single work unit: one netMHCpan call for (chunk_idx, batch_idx).
    function run_unit(chunk_idx::Int, batch_idx::Int, pep_file::String)
        batch = allele_batches[batch_idx]
        out_file = joinpath(folder_path, "_temp_netMHCpan_output_$(chunk_idx)_$(batch_idx).tsv")
        cmd = Cmd([netMHCpan_path, "-p", pep_file, "-xls", "-a", join(batch, ","), "-xlsfile", out_file])
        try
            if verbose
                log_file = joinpath(folder_path, "_temp_netMHCpan_log_$(chunk_idx)_$(batch_idx).txt")
                open(log_file, "w") do io; run(pipeline(cmd, stdout=io, stderr=io)); end
            else
                run(pipeline(cmd, stdout=devnull, stderr=devnull))
            end
            if !isfile(out_file)
                error("NetMHCpan did not create output for chunk $(chunk_idx), batch $(batch_idx).")
            end
        catch e
            println("\nError on chunk $(chunk_idx), batch $(batch_idx): ", e)
            println("Command: ", join(cmd.exec, " "))
            exit(1)
        end
        return out_file
    end

    if threads == 1
        # Sequential: iterate over chunks then batches, showing progress per unit.
        for chunk_idx in 1:total_chunks
            pep_file = temp_pep_files[chunk_idx]
            pep_file === nothing && continue
            for batch_idx in 1:n_batches
                pct = Int(floor(100 * ((chunk_idx-1)*n_batches + (batch_idx-1)) / max(1, total_units)))
                status("Chunk $(chunk_idx)/$(total_chunks), batch $(batch_idx)/$(n_batches) ($(length(peptide_chunks[chunk_idx])) peptides, $(total_alleles) alleles). $(pct)% complete."; overwrite=true)
                results_matrix[chunk_idx, batch_idx] = run_unit(chunk_idx, batch_idx, pep_file)
            end
        end
    else
        # Parallel: all (chunk, batch) units compete for `threads` concurrent slots.
        # This naturally adapts: many chunks + few batches → threads spread across chunks;
        # few chunks + many batches → threads spread across batches; mixed → both.
        sem = Base.Semaphore(threads)
        progress_lock = ReentrantLock()
        completed_ref = Ref(0)

        reporter_task = @async begin
            last_len = 0
            while true
                sleep(0.5)
                pct, done, count = lock(progress_lock) do
                    c = completed_ref[]
                    p = total_units == 0 ? 100 : Int(floor(100 * c / total_units))
                    p, c >= total_units, c
                end
                line = "$(count)/$(total_units) units complete. $(pct)%"
                print("\r", line)
                if length(line) < last_len; print(" "^(last_len - length(line))); end
                flush(stdout)
                last_len = length(line)
                if done; break; end
            end
            println()
        end

        all_tasks = Task[]
        for chunk_idx in 1:total_chunks
            pep_file = temp_pep_files[chunk_idx]
            pep_file === nothing && continue
            for batch_idx in 1:n_batches
                ci, bi, pf = chunk_idx, batch_idx, pep_file   # capture for closure
                t = @async begin
                    Base.acquire(sem)
                    try
                        out_file = run_unit(ci, bi, pf)
                        lock(progress_lock) do
                            results_matrix[ci, bi] = out_file
                            completed_ref[] += 1
                        end
                    finally
                        Base.release(sem)
                    end
                end
                push!(all_tasks, t)
            end
        end
        for t in all_tasks; wait(t); end
        wait(reporter_task)
    end

    # Merge: for each chunk, horizontally stitch per-batch files (fixed cols from batch 1,
    # allele cols from all batches in order), then vertically concatenate chunks.
    # netMHCpan 4.2 appends "Ave" and "NB" summary columns at the end of each row/header.
    # These are row-level averages (not per-allele); we detect and strip them so the output
    # has exactly 4 cols per allele, which process_output.pl requires.
    print("\n")
    status("Merging chunk outputs into $xlsfile_path ...")

    open(xlsfile_path, "w") do out_io
        # Line 1: allele names (tab-separated, two leading tabs) for process_output.pl
        println(out_io, "\t\t" * join(allele_list, '\t'))

        if isempty(peptide_list)
            header = String["Pos", "Peptide", "ID"]
            for _ in allele_list; append!(header, ["core", "icore", "EL-score", "EL_Rank"]); end
            println(out_io, join(header, '\t'))
        else
            col_header_written = false
            trim_counts = Int[]  # trailing summary cols to drop from each batch (set on first chunk)
            for chunk_idx in 1:total_chunks
                temp_pep_files[chunk_idx] === nothing && continue
                batch_files = [results_matrix[chunk_idx, bi] for bi in 1:n_batches]
                any(f -> f === nothing || !isfile(f), batch_files) && continue
                lines_per_batch = [[split(line, '\t') for line in eachline(f)] for f in batch_files]
                nlines = length(lines_per_batch[1])
                for bl in lines_per_batch
                    length(bl) == nlines || error("Mismatched line counts across batches in chunk $(chunk_idx).")
                end
                if !col_header_written
                    fixed_hdr = lines_per_batch[1][2][1:min(3, length(lines_per_batch[1][2]))]
                    raw_hdr_cols = [bl[2][4:end] for bl in lines_per_batch]
                    # Detect trailing Ave/NB summary cols per batch from the column header
                    trim_counts = [length(c) >= 2 && strip(c[end]) == "NB" && strip(c[end-1]) == "Ave" ? 2 : 0 for c in raw_hdr_cols]
                    allele_cols = [c[1:end-t] for (c, t) in zip(raw_hdr_cols, trim_counts)]
                    println(out_io, join(vcat(fixed_hdr, reduce(vcat, allele_cols)), '\t'))
                    col_header_written = true
                end
                for row_idx in 3:nlines
                    fp = lines_per_batch[1][row_idx]
                    raw_data_cols = [bl[row_idx][4:end] for bl in lines_per_batch]
                    data_allele_cols = [c[1:end-t] for (c, t) in zip(raw_data_cols, trim_counts)]
                    row = vcat(fp[1:min(3,length(fp))], reduce(vcat, data_allele_cols))
                    println(out_io, join(row, '\t'))
                end
            end
            col_header_written || error("Could not construct column header: no valid chunk output files found.")
        end
    end

    status("Merged output written to $xlsfile_path")
    # Cleanup temp files
    status("Cleaning up temporary files...")
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
