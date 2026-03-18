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
include("path_utils.jl")
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
    # threads already determined and capped above (supports --t max)

    netMHCpan_exe = get_netMHCpan_executable()
    println("Using NetMHCpan executable: ", netMHCpan_exe)

    # File paths
    # Inputs: ignore suffix; discover with latest
    alleles_file = resolve_read(joinpath(folder_path, "alleles.txt"); suffix="", latest=latest)
    peptides_file = resolve_read(joinpath(folder_path, "Peptides.pep"); suffix="", latest=latest)
    xlsfile_path = resolve_write(joinpath(folder_path, "netMHCpan_output.tsv"); suffix=suffix)
    cache_file   = resolve_write(joinpath(folder_path, "results_cache.jls"); suffix=suffix)

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

    # Read alleles (remove * for netMHCpan; strip comments/extra fields)
    allele_list = open(alleles_file) do file
        [replace(split(line, r"\s+")[1], "*" => "") for line in readlines(file) if !isempty(strip(line)) && !startswith(strip(line), "#")]
    end

    # Read peptides
    peptide_list = open(peptides_file) do file
        [strip(line) for line in readlines(file) if !isempty(strip(line))]
    end

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

    # Chunk peptides into batches of 1000
    chunk_size = 1000
    total_peptides = length(peptide_list)
    peptide_chunks = [peptide_list[i:min(i+chunk_size-1, total_peptides)] for i in 1:chunk_size:total_peptides]
    total_chunks = length(peptide_chunks)
    total_alleles = length(allele_list)
    total_units  = total_chunks * n_batches   # each unit = one netMHCpan call

    # Write all peptide chunk files upfront (fast sequential I/O).
    # netMHCpan only reads these files, so concurrent batch tasks can safely share them.
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
        cmd = Cmd([netMHCpan_exe, "-p", pep_file, "-xls", "-a", join(batch, ","), "-xlsfile", out_file])
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
            println("Error on chunk $(chunk_idx), batch $(batch_idx): ", e)
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
    # Detect and strip trailing "Ave"/"NB" summary cols that netMHCpan 4.2 appends.
    print("\n")
    status("Merging chunk outputs into $xlsfile_path ...")
    open(xlsfile_path, "w") do out_io
        println(out_io, "\t\t" * join(allele_list, '\t'))

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
        if !col_header_written
            header = String["Pos", "Peptide", "ID"]
            for _ in allele_list; append!(header, ["core", "icore", "EL-score", "EL_Rank"]); end
            println(out_io, join(header, '\t'))
        end
    end
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
