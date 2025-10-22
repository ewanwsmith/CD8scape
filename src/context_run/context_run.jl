#!/usr/bin/env julia
"""
This script orchestrates the full context peptide pipeline.

Purpose:
- Runs all context-specific processing steps in order.
- Handles argument parsing and error checking.
- Passes options to each step as needed.

Usage:
    julia context_run.jl <folder_path> [--n_loci <number_of_loci>] [--mode panel|supertype] [--override] [--seed <random_seed>]

Arguments:
    <folder_path>   Path to the folder containing input files.
    --n_loci        Number of loci to simulate (default: 1000).
    --mode          Run mode: panel or supertype (default: panel).
    --override      Optional flag to override existing outputs.
    --seed          Random number seed (default: 1320).
"""

using CSV, DataFrames, FilePathsBase
import Dates

# Simple directory-based lock helpers (cooperative locking). This is not a full
# cross-platform advisory lock, but it prevents accidental concurrent writers in
# typical usage where multiple orchestrators run on the same filesystem.
const _LOCK_HANDLES = Dict{String, IO}()

# Memoize whether descriptor-based flock is supported in this Julia build.
const _FLOCK_SUPPORTED = Ref{Union{Bool, Nothing}}(nothing)

function flock_supported()
    if _FLOCK_SUPPORTED[] !== nothing
        return _FLOCK_SUPPORTED[]
    end
    # Quick checks: ensure Libc exports flock and some fileno entry exists
    try
        if !isdefined(Libc, :flock)
            _FLOCK_SUPPORTED[] = false
            return false
        end
    catch
        _FLOCK_SUPPORTED[] = false
        return false
    end

    # Try to obtain a file descriptor for a temporary IO and attempt a non-blocking lock.
    tmp = nothing
    try
        tmp = open("/dev/null", "w")
        # attempt fileno using Libc or Base.Libc if available
        fd = try
            Libc.fileno(tmp)
        catch
            try
                Base.Libc.fileno(tmp)
            catch
                nothing
            end
        end
        if fd === nothing
            _FLOCK_SUPPORTED[] = false
            return false
        end
        # Try flock call non-blocking; if not available or fails, consider unsupported
        try
            res = Libc.flock(fd, Libc.LOCK_EX | Libc.LOCK_NB)
            if res == 0
                # release
                try
                    Libc.flock(fd, Libc.LOCK_UN)
                catch
                end
                _FLOCK_SUPPORTED[] = true
            else
                _FLOCK_SUPPORTED[] = false
            end
        catch
            _FLOCK_SUPPORTED[] = false
        end
    catch
        _FLOCK_SUPPORTED[] = false
    finally
        try
            if tmp !== nothing
                close(tmp)
            end
        catch
        end
    end
    return _FLOCK_SUPPORTED[]
end

# Helper to obtain file descriptor in a portable way. Tries Libc then Base.Libc.
function safe_fileno(io::IO)
    # Attempt common locations for fileno
    try
        return Libc.fileno(io)
    catch
        try
            return Base.Libc.fileno(io)
        catch
            throw(ErrorException("fileno-unavailable"))
        end
    end
end

# Acquire an exclusive file lock using flock. Returns true on success, false on timeout.
function lock_acquire(lockfile::AbstractString; timeout_seconds::Int=30)
    # Directory-lock is the default for cross-platform reliability. Only attempt
    # descriptor-based flocking if the user explicitly requests it via
    # CD8S_USE_FLOCK=1 and the runtime indicates flock is supported.
    if !(get(ENV, "CD8S_USE_FLOCK", "0") == "1" && flock_supported())
        lockdir = lockfile * ".dirlock"
        start = time()
        while true
            try
                mkdir(lockdir)
                try
                    open(joinpath(lockdir, "owner.txt"), "w") do out
                        println(out, getpid(), " ", Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS"))
                    end
                catch
                end
                return true
            catch
                if time() - start > timeout_seconds
                    return false
                end
                sleep(0.1)
            end
        end
    end

    # Try to use Libc.flock if present (POSIX systems). If unavailable, fall back to dir lock.
    try
        # open (or create) the lock file for read/write
    io = open(lockfile, "w+")
    fd = safe_fileno(io)
        start = time()
        while true
            # Try non-blocking exclusive lock
            res = try
                Libc.flock(fd, Libc.LOCK_EX | Libc.LOCK_NB)
            catch
                # flock not available or failed; close io and fallback
                close(io)
                throw(ErrorException("flock-unavailable"))
            end
            if res == 0
                # success
                _LOCK_HANDLES[lockfile] = io
                # write diagnostics
                try
                    seek(io, 0)
                    write(io, string(getpid(), " ", Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS"), "\n"))
                    flush(io)
                catch
                end
                return true
            end
            if time() - start > timeout_seconds
                close(io)
                return false
            end
            sleep(0.1)
        end
    catch e
        # Fallback to cooperative directory lock if flock is not available
        if occursin("flock-unavailable", String(e)) || isa(e, ErrorException)
            lockdir = lockfile * ".dirlock"
            start = time()
            while true
                try
                    mkdir(lockdir)
                    try
                        open(joinpath(lockdir, "owner.txt"), "w") do out
                            println(out, getpid(), " ", Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS"))
                        end
                    catch
                    end
                    return true
                catch
                    if time() - start > timeout_seconds
                        return false
                    end
                    sleep(0.1)
                end
            end
        else
            # unknown error; rethrow
            rethrow(e)
        end
    end
end

function lock_release(lockfile::AbstractString)
    # Release flock handle if present, else remove dir lock
    try
        if haskey(_LOCK_HANDLES, lockfile)
            io = _LOCK_HANDLES[lockfile]
            fd = safe_fileno(io)
            try
                Libc.flock(fd, Libc.LOCK_UN)
            catch
            end
            try
                close(io)
            catch
            end
            delete!(_LOCK_HANDLES, lockfile)
            # attempt to remove the lockfile
            try
                rm(lockfile)
            catch
            end
            return true
        else
            # Try directory fallback
            lockdir = lockfile * ".dirlock"
            try
                rm(lockdir; recursive=true)
                return true
            catch e
                println("Warning: failed to remove lockdir $lockdir: $e")
                return false
            end
        end
    catch e
        println("Warning: lock_release encountered error: $e")
        return false
    end
end

# Atomic write CSV: write to temp file in same directory then rename
function atomic_write_csv(path::AbstractString, df::DataFrame)
    dir = dirname(path)
    base = basename(path)
    # Ensure target directory exists
    try
        mkpath(dir)
    catch
        # ignore
    end

    # Create a temporary filename in the same directory to allow an atomic rename.
    # Use process id and timestamp to minimise collisions.
    tmp_path = joinpath(dir, "$(base).tmp.$(getpid()).$(Int(time()*1000))")
    try
        # Write directly to the temporary path, then atomically rename to final path.
        CSV.write(tmp_path, df)
        Base.rename(tmp_path, path)
    catch e
        # try cleanup of tmp_path if it still exists
        try
            if isfile(tmp_path)
                Base.rm(tmp_path)
            end
        catch
        end
        rethrow(e)
    end
end


function main()


    if !isdefined(Main, :ARGS) || isempty(ARGS)
        println("Error: folder_path must be provided as a command-line argument.")
        exit(1)
    end


    folder_path = ARGS[1]
    frames_path = joinpath(folder_path, "frames.csv")

    # Define attempt before any references to it
    attempt = 1

    if !isfile(frames_path)
        println("Error: frames.csv not found in the provided folder path: $folder_path")
        exit(1)
    end
    # Default number of loci
    n_loci = 1000
    mode = "panel"
    override = false
    seed = 1320
    verbose = false

    # Parse additional arguments
    for (i, arg) in enumerate(ARGS)
        if arg == "--n_loci" && i < length(ARGS)
                println("[DIAG] Skipping cache write: context_peptides.pep missing for attempt $attempt.")
            try
                n_loci = parse(Int, ARGS[i + 1])
            catch
                println("Error: Invalid value for --n_loci.")
                exit(1)
            end
        elseif arg == "--mode" && i < length(ARGS)
            mode = ARGS[i + 1]
                println("[DIAG] Skipping cache write: all loci already processed for attempt $attempt.")
            if mode ∉ ["panel", "supertype"]
                println("Error: Invalid value for --mode.")
                exit(1)
            end
                println("[DIAG] Skipping cache write: error filtering duplicates for attempt $attempt.")
        elseif arg == "--verbose" || arg == "-v"
            verbose = true
        elseif arg == "--override"
            override = true
        elseif arg == "--seed" && i < length(ARGS)
            try
                seed = parse(Int, ARGS[i + 1])
            catch
                println("Error: Invalid value for --seed.")
                println("[DIAG] Skipping cache write: error in generation for attempt $attempt.")
                exit(1)
            end
        end
    end

    # ...existing code...
    generate_script = joinpath(@__DIR__, "generate_context_peptides.jl")
                println("[DIAG] Skipping cache write: error reading peptide file for attempt $attempt.")

    # Make max attempts configurable (default increased to 50)
    max_attempts = 200
    # Allow override from CLI --max_attempts
    for (i, arg) in enumerate(ARGS)
        if arg == "--max_attempts" && i < length(ARGS)
            try
                max_attempts = max(1, parse(Int, ARGS[i + 1]))
            catch
                println("Warning: invalid value for --max_attempts; using default $max_attempts")
            end
        end
    end

    attempt = 1
                    println("[DIAG] Skipping cache write: NetMHCpan run failed for attempt $attempt.")
    current_n = n_loci

    # Prepare a persistent cache directory to store per-attempt outputs so we don't re-run
    # previously-generated variants. Each attempt is executed in its own subfolder under
    # `context_cache/attempt_*`. Successful loci are aggregated into master files under
    # the primary folder_path and only loci that pass HMBR filtering are included there.
    cache_dir = joinpath(folder_path, "context_cache")
    mkpath(cache_dir)

                    println("[DIAG] Skipping cache write: NetMHCpan did not produce output for attempt $attempt.")
    master_hm_df = DataFrame()
    master_best_df = DataFrame()

    while attempt <= max_attempts
        println("Context generation attempt $attempt: generating $current_n loci (seed=$seed)")

        attempt_dir = joinpath(cache_dir, "attempt_$(attempt)_seed_$(seed)_n_$(current_n)")

        # If this attempt was already executed earlier, reuse its outputs instead of re-running.
        if isdir(attempt_dir) && isfile(joinpath(attempt_dir, "context_harmonic_mean_best_ranks.csv"))
            println("Found cached results for attempt $attempt at $attempt_dir — reusing outputs")
        else
            mkpath(attempt_dir)
            # Copy alleles.txt only if mode is panel
            if mode == "panel"
                try
                    alleles_candidates = filter(f -> lowercase(basename(f)) == "alleles.txt", readdir(folder_path; join=true))
                    if !isempty(alleles_candidates)
                        alleles_src = first(alleles_candidates)
                        alleles_dest = joinpath(attempt_dir, "alleles.txt")
                        if !isfile(alleles_dest)
                            cp(alleles_src, alleles_dest)
                            println("Copied alleles file to attempt directory: $(alleles_dest)")
                        end
                    else
                        println("Note: no alleles.txt found in main folder ($folder_path); NetMHCpan runner may skip alleles.")
                    end
                catch e
                    println("Warning: could not copy alleles.txt into attempt dir: $e")
                end
            end
            # Always copy supertype_panel.csv if mode is supertype
            if mode == "supertype"
                try
                    supertype_src_candidates = filter(f -> lowercase(basename(f)) == "supertype_panel.csv", readdir(folder_path; join=true))
                    if !isempty(supertype_src_candidates)
                        supertype_src = first(supertype_src_candidates)
                        supertype_dest = joinpath(attempt_dir, "supertype_panel.csv")
                        if !isfile(supertype_dest)
                            cp(supertype_src, supertype_dest)
                            println("Copied supertype_panel.csv to attempt directory: $(supertype_dest)")
                        end
                    else
                        println("Warning: no supertype_panel.csv found in main folder ($folder_path); supertype run may fail.")
                    end
                catch e
                    println("Warning: could not copy supertype_panel.csv into attempt dir: $e")
                end
            end
            # Copy frames.csv into attempt directory so generate script can operate there
            dest_frames = joinpath(attempt_dir, "frames.csv")
            try
                if isfile(dest_frames)
                    println("Note: frames.csv already exists in attempt directory $attempt_dir — reusing existing file")
                else
                    cp(frames_path, dest_frames)
                end
            catch e
                println("Error copying frames.csv to attempt directory: $e")
                exit(1)
            end

            # Run generation in the attempt directory
            local gen_cmd = `julia $generate_script --folder $attempt_dir --n_loci $current_n --seed $seed`
            try
                run(gen_cmd)
            catch e
                println("Error running generation in attempt $attempt: $e")
                # escalate the seed and try next attempt
                attempt += 1
                seed += 1
                current_n = current_n + max(10, Int(ceil(0.5 * current_n)))
                continue
            end

            # After generation we proceed with the rest of the pipeline inside attempt_dir
            # Case-insensitive lookup of context_peptides_labels.csv / context_peptides.pep (fresh from generation)
            labels_file = joinpath(attempt_dir, "context_peptides_labels.csv")
            matches = filter(f -> lowercase(basename(f)) == "context_peptides.pep", readdir(attempt_dir; join=true))
            if isempty(matches)
                println("Warning: context_peptides.pep not found in attempt directory $attempt_dir after generation. Skipping attempt.")
                # don't abort whole pipeline; try again with a new seed/size
                attempt += 1
                seed += 1
                current_n = current_n + max(10, Int(ceil(0.5 * current_n)))
                continue
            end

            # If we have already aggregated passed loci from previous attempts, avoid re-running
            # those loci: load labels CSV and filter out rows where Locus is already present in master_hm_df
            if isfile(labels_file) && nrow(master_hm_df) > 0
                try
                    labdf = CSV.read(labels_file, DataFrame)
                    original_count = nrow(labdf)
                    # keep only loci not already present in master_hm_df
                    labdf = filter(row -> !(row.Locus in master_hm_df.Locus), labdf)
                    filtered_count = nrow(labdf)
                    if filtered_count == 0
                        println("All generated loci in attempt $attempt have already been processed in previous attempts. Skipping attempt.")
                        attempt += 1
                        seed += 1
                        current_n = max(current_n + 10, current_n + Int(ceil(0.1 * current_n)))
                        continue
                    elseif filtered_count < original_count
                        println("Filtered out $(original_count - filtered_count) duplicate loci from attempt $attempt; $filtered_count new loci remain.")
                        # write filtered labels and regenerate the .pep file for downstream steps
                        CSV.write(joinpath(attempt_dir, "context_peptides_labels.filtered.csv"), labdf)
                        open(joinpath(attempt_dir, "context_peptides.pep"), "w") do io
                            for r in eachrow(labdf)
                                println(io, r.Peptide)
                            end
                        end
                        # update matches to point to the (now-overwritten) pep file
                        matches = [joinpath(attempt_dir, "context_peptides.pep")]
                    end
                catch e
                    println("Warning: could not filter generated labels for duplicates in attempt $attempt: $e")
                end
            end

            # Call clean_peptides_context.jl in attempt_dir
            clean_script = joinpath(@__DIR__, "clean_peptides_context.jl")
            try
                run(`julia $clean_script $attempt_dir`)
            catch e
                println("Error running clean_peptides_context.jl in attempt $attempt: $e")
                exit(1)
            end

            netmhc_script = joinpath(@__DIR__, "run_netMHCpan_context.jl")
            local pepfile = first(matches)

            # Peptide-level cache: a master processed CSV keyed by Peptide lives in the main folder_path.
            # If a peptide has already been processed previously, reuse its processed rows rather than
            # re-running NetMHCpan for that sequence. Only missing peptide sequences are run through
            # NetMHCpan for this attempt.
            master_cache_file = joinpath(folder_path, "peptide_netmhc_processed.csv")
            cached_df = DataFrame()
            try
                if isfile(master_cache_file)
                    cached_df = CSV.read(master_cache_file, DataFrame)
                    # Normalize peptide column name so downstream code can rely on :Peptide
                    try
                        pcol = nothing
                        for nm in names(cached_df)
                            try
                                if lowercase(String(nm)) == "peptide"
                                    pcol = nm; break
                                end
                            catch
                            end
                        end
                        if pcol !== nothing && pcol != :Peptide
                            rename!(cached_df, pcol => :Peptide)
                        end
                    catch e
                        println("[DIAG] Warning: could not normalize peptide column name: $e")
                    end
                    println("[DIAG] Loaded peptide NetMHCpan cache with $(nrow(cached_df)) rows from $master_cache_file")
                else
                    println("[DIAG] No peptide NetMHCpan cache found at $master_cache_file; starting with empty cache.")
                    cached_df = DataFrame()
                end
            catch e
                println("Warning: could not read master peptide cache: $e")
                cached_df = DataFrame()
            end

            # Read unique peptides for this attempt
            attempt_peps = String[]
            try
                attempt_peps = readlines(pepfile)
                attempt_peps = unique(filter(x -> !isempty(strip(x)), attempt_peps))
            catch e
                println("Error reading peptide file $pepfile: $e")
                attempt += 1; seed += 1
                current_n = current_n + max(10, Int(ceil(0.5 * current_n)))
                continue
            end

            # Find which peptides are missing from the cache. Be robust to column name types
            cached_peptides = String[]
            peptide_col = nothing
            if nrow(cached_df) > 0
                for nm in names(cached_df)
                    try
                        if lowercase(String(nm)) == "peptide"
                            peptide_col = nm
                            break
                        end
                    catch
                        # ignore coercion errors
                    end
                end
                if peptide_col !== nothing
                    try
                        cached_peptides = unique(String.(cached_df[!, peptide_col]))
                    catch e
                        println("[DIAG] Warning: could not coerce cached peptide column to String: $e")
                        cached_peptides = String[]
                    end
                end
            end
            # Diagnostic: print a sample of attempt_peps and cached_peptides for debugging
            println("[DIAG] Sample attempt_peps (first 5): ", join(first(attempt_peps, min(5, length(attempt_peps))), ", "))
            if !isempty(cached_peptides)
                println("[DIAG] Sample cached_df.Peptide (first 5): ", join(first(cached_peptides, min(5, length(cached_peptides))), ", "))
            else
                println("[DIAG] Sample cached_df.Peptide (first 5): <none found> (available columns: $(names(cached_df)))")
            end
            missing_peptides = setdiff(attempt_peps, cached_peptides)

            # If there are missing peptides, run NetMHCpan only on them and append processed results to master cache
            if !isempty(missing_peptides)
                missing_file = joinpath(attempt_dir, "missing_peptides.pep")
                open(missing_file, "w") do io
                    for p in missing_peptides
                        println(io, p)
                    end
                end

                local netmhc_cmd = `julia $netmhc_script --folder $attempt_dir --peptides $missing_file --mode $mode`
                # Run NetMHCpan for missing peptides
                try
                    run(netmhc_cmd)
                catch e
                    println("Warning: NetMHCpan run failed for missing-peptide set in attempt $attempt: $e")
                    # treat as recoverable: increment attempt parameters and retry
                    attempt += 1
                    seed += 1
                    incr = max(100, Int(ceil(0.2 * current_n)))
                    max_per_attempt = 200_000
                    current_n = min(current_n + incr, max_per_attempt)
                    println("Retrying with seed=$seed and current_n=$current_n (after NetMHCpan failure)")
                    continue
                end

                # Process NetMHCpan output for missing peptides
                netmhcpan_output = joinpath(attempt_dir, "netMHCpan_output.tsv")
                if !isfile(netmhcpan_output)
                    println("Warning: NetMHCpan did not produce output for missing peptides in attempt $attempt. Will retry generation if attempts remain.")
                    attempt += 1; seed += 1
                    incr = max(100, Int(ceil(0.2 * current_n)))
                    max_per_attempt = 200_000
                    current_n = min(current_n + incr, max_per_attempt)
                    continue
                end

                try
                    run(`perl $(joinpath(@__DIR__, "process_output_context.pl")) $netmhcpan_output`)
                catch e
                    println("Error running process_output_context.pl for missing peptides in attempt $attempt: $e")
                    exit(1)
                end

                # Read processed results for missing peptides and append to master cache
                new_processed = joinpath(attempt_dir, "context_processed_netMHCpan_output.csv")
                if isfile(new_processed)
                    try
                        df_new = CSV.read(new_processed, DataFrame)
                        println("[DIAG] Attempt processed output: $(nrow(df_new)) rows, sample Peptide: ", join(first(df_new.Peptide, min(5, nrow(df_new))), ", "))
                        if nrow(cached_df) == 0
                            cached_df = df_new
                        else
                            # Append and deduplicate cached entries conservatively by all columns
                            combined = vcat(cached_df, df_new)
                            cached_df = unique(combined)
                        end
                        println("[DIAG] Master cache after append: $(nrow(cached_df)) rows, sample Peptide: ", join(first(cached_df.Peptide, min(5, nrow(cached_df))), ", "))
                        # Persist updated master cache
                        # Attempt to persist updated master cache with good diagnostics and a safe fallback.
                        lockdir = joinpath(folder_path, ".peptide_cache_lock")
                        got = false
                        try
                            got = lock_acquire(lockdir; timeout_seconds=30)
                        catch e
                            println("Warning: lock_acquire raised an exception: $e")
                            @show catch_backtrace = catch_backtrace()
                            got = false
                        end

                        if !got
                            println("Warning: could not acquire lock for master peptide cache; attempting non-atomic write to $master_cache_file")
                            try
                                CSV.write(master_cache_file, cached_df)
                                println("Appended $(nrow(df_new)) processed rows to master peptide cache (non-atomic) at $master_cache_file")
                            catch e
                                println("ERROR: non-atomic CSV.write failed for $master_cache_file: $e")
                                showerror(stderr, e)
                                println()
                            end
                        else
                            # We have the lock: try atomic write, but on failure fall back to non-atomic write and always release lock.
                            try
                                try
                                    atomic_write_csv(master_cache_file, cached_df)
                                    println("Appended $(nrow(df_new)) processed rows to master peptide cache (atomic) at $master_cache_file")
                                catch e
                                    println("ERROR: atomic_write_csv failed for $master_cache_file: $e")
                                    showerror(stderr, e)
                                    println()
                                    println("Falling back to non-atomic CSV.write for $master_cache_file")
                                    try
                                        CSV.write(master_cache_file, cached_df)
                                        println("Appended $(nrow(df_new)) processed rows to master peptide cache (fallback non-atomic) at $master_cache_file")
                                    catch e2
                                        println("ERROR: fallback non-atomic CSV.write also failed for $master_cache_file: $e2")
                                        showerror(stderr, e2)
                                        println()
                                    end
                                end
                            finally
                                try
                                    lock_release(lockdir)
                                catch e
                                    println("Warning: lock_release raised an exception: $e")
                                end
                            end
                        end

                        # Confirm cache file exists and print row count
                        if isfile(master_cache_file)
                            try
                                tmp_df = CSV.read(master_cache_file, DataFrame)
                                println("[DIAG] Confirmed cache file written: $(nrow(tmp_df)) rows at $master_cache_file")
                            catch e
                                println("[DIAG] WARNING: cache file exists but could not be read: $e")
                            end
                        else
                            println("[DIAG] ERROR: Cache file $master_cache_file was not written!")
                        end
                    catch e
                        println("Warning: could not read processed output for missing peptides: $e")
                    end
                else
                    println("Warning: expected processed CSV for missing peptides not found at $new_processed")
                end
            else
                println("All $(length(attempt_peps)) peptides for attempt $attempt are present in the master cache; skipping NetMHCpan run.")
            end

                # Build attempt-level processed output by selecting cached rows for the peptides in this attempt
                try
                if nrow(cached_df) == 0
                    println("Warning: master peptide cache is empty after NetMHCpan step for attempt $attempt. No processed results available.")
                end
                # Select rows where Peptide is in attempt_peps (robust to column name casing/types)
                attempt_rows = DataFrame()
                if nrow(cached_df) > 0
                    # ensure peptide_col is set (may have been set earlier)
                    if peptide_col === nothing
                        for nm in names(cached_df)
                            try
                                if lowercase(String(nm)) == "peptide"
                                    peptide_col = nm; break
                                end
                            catch
                            end
                        end
                    end
                    if peptide_col !== nothing
                        attempt_rows = filter(row -> row[peptide_col] in attempt_peps, cached_df)
                    else
                        println("[DIAG] WARNING: could not find Peptide column in master cache; producing empty attempt_rows")
                    end
                end
                println("[DIAG] Attempt output for scoring: $(nrow(attempt_rows)) rows, sample Peptide: ", (nrow(attempt_rows)>0 ? join(first(attempt_rows[!, peptide_col], min(5, nrow(attempt_rows))), ", ") : "<none>"))
                CSV.write(joinpath(attempt_dir, "context_processed_netMHCpan_output.csv"), attempt_rows)
                # Diagnostic: report how many rows were written and how many labeled rows will be produced after joining
                written = nrow(attempt_rows)
                labeled_rows = "unknown"
                try
                    # if labels file exists, compute exact joined row count to mirror process_scores_context.jl behavior
                    labf = joinpath(attempt_dir, "context_peptides_labels.csv")
                    if isfile(labf)
                        labdf = CSV.read(labf, DataFrame)
                        # perform a left join like process_scores_context.jl and count resulting rows
                        joined_temp = leftjoin(attempt_rows, labdf, on=["Peptide"])  # may expand rows
                        labeled_rows = nrow(joined_temp)
                    else
                        labeled_rows = 0
                    end
                catch e
                    labeled_rows = "unknown (error: $(e))"
                end
                println("Wrote attempt processed NetMHCpan output with $written rows to attempt directory: $attempt_dir; estimated labeled rows after join: $labeled_rows")
            catch e
                println("Error constructing attempt processed NetMHCpan output for attempt $attempt: $e")
                exit(1)
            end

            # Run context-specific score processing for this attempt
            if !success(`julia --project=. $(joinpath(@__DIR__, "process_scores_context.jl")) $attempt_dir`)
                println("Error running process_scores_context.jl for attempt $attempt")
                exit(1)
            end

            # After scoring, read the produced context_scores.csv to report exact joined row count
            try
                scores_file = joinpath(attempt_dir, "context_scores.csv")
                if isfile(scores_file)
                    scdf = CSV.read(scores_file, DataFrame)
                    println("Scoring produced $(nrow(scdf)) rows in context_scores.csv for attempt $attempt")
                    # If this differs from the number of processed peptide rows written earlier, explain why
                    try
                        if written !== nothing && written != "unknown"
                            if nrow(scdf) != written
                                println("Note: $(nrow(scdf)) rows found in context_scores.csv but $written unique peptide rows were written from NetMHCpan. This can happen because peptide labels may map multiple rows to the same peptide (one peptide -> multiple peptide_labels), causing the join to expand rows.")
                            end
                        end
                    catch
                        # ignore any diagnostic failures
                    end
                else
                    println("Warning: expected scores file $scores_file not found after scoring step for attempt $attempt")
                end
            catch e
                println("Warning: could not read context_scores.csv for attempt $attempt: $e")
            end

            # Run best ranks for this attempt
            best_ranks_context_script = joinpath(@__DIR__, "process_best_ranks_context.jl")
            if mode == "supertype"
                best_ranks_context_cmd = `julia --project=. $best_ranks_context_script $attempt_dir --supertype`
            else
                best_ranks_context_cmd = `julia --project=. $best_ranks_context_script $attempt_dir`
            end
            run(best_ranks_context_cmd)
        end

        # Inspect resulting harmonic mean output for this attempt and merge passed loci into master tables

        attempt_hm_file = joinpath(attempt_dir, "context_harmonic_mean_best_ranks.csv")
        attempt_best_file = joinpath(attempt_dir, "context_best_ranks.csv")
        attempt_passed = 0
        attempt_total = 0
        attempt_filtered = 0
        if isfile(attempt_hm_file)
            try
                df_hm = CSV.read(attempt_hm_file, DataFrame)
                attempt_passed = nrow(df_hm)
                # Try to get the total loci before filtering (if available)
                # This requires reading the unfiltered file if it exists, otherwise estimate from best ranks
                if isfile(attempt_best_file)
                    try
                        df_best = CSV.read(attempt_best_file, DataFrame)
                        attempt_total = length(unique(df_best.Locus))
                    catch
                        attempt_total = attempt_passed
                    end
                else
                    attempt_total = attempt_passed
                end
                attempt_filtered = attempt_total - attempt_passed
                println("Attempt $attempt: $attempt_filtered loci filtered out as non-binding (of $attempt_total total loci in attempt)")

                # If master is empty just take df_hm, otherwise append only new Locus rows
                if nrow(master_hm_df) == 0
                    master_hm_df = df_hm
                else
                    new_rows = filter(row -> !(row.Locus in master_hm_df.Locus), df_hm)
                    if nrow(new_rows) > 0
                        master_hm_df = vcat(master_hm_df, new_rows)
                    end
                end

                # Aggregate corresponding best ranks for passed loci (if available)
                if isfile(attempt_best_file)
                    try
                        df_best = CSV.read(attempt_best_file, DataFrame)
                        # Filter df_best to only those loci that are present in df_hm
                        loci_in_hm = unique(df_hm.Locus)
                        df_best_filtered = filter(row -> row.Locus in loci_in_hm, df_best)
                        if nrow(master_best_df) == 0
                            master_best_df = df_best_filtered
                        else
                            # avoid duplicate rows for the same Peptide_label / Locus
                            append_rows = filter(row -> !(row.Peptide_label in master_best_df.Peptide_label), df_best_filtered)
                            if nrow(append_rows) > 0
                                master_best_df = vcat(master_best_df, append_rows)
                            end
                        end
                    catch e
                        println("Warning: could not read or merge best ranks for attempt $attempt: $e")
                    end
                end

            catch e
                println("Warning: could not read harmonic mean file for attempt $attempt: $e")
            end
        else
            println("No harmonic mean file produced for attempt $attempt at $attempt_dir")
        end

        total_passed = nrow(master_hm_df)
        println("Attempt $attempt result: $attempt_passed loci passed in this attempt; $total_passed loci aggregated total")

        # Cleanup attempt directory to avoid leaving large intermediate files around,
        # unless the user asked for verbose output.
        if !verbose
            try
                # remove the attempt directory recursively
                rm(attempt_dir; recursive=true)
                println("Removed temporary attempt directory: $attempt_dir")
            catch e
                println("Warning: failed to remove attempt directory $attempt_dir: $e")
            end
        else
            println("Verbose mode: preserving attempt directory: $attempt_dir")
        end

        if total_passed >= n_loci
            println("Reached target of $n_loci loci after $attempt attempt(s);")
            println("Note: more loci may have been aggregated than requested; all aggregated loci (>$n_loci) are kept in the output.")
            break
        else
            if attempt == max_attempts
                println("Warning: reached max attempts ($max_attempts) but only $total_passed loci passed. Proceeding with available loci.")
                break
            end
            # Increase the number of loci to generate in the next attempt to try to cover the deficit.
            deficit = n_loci - total_passed
            # Grow by deficit plus 50% overshoot to reduce number of retries
            grow = deficit + Int(ceil(0.5 * deficit))
                # cap per-attempt generation to avoid runaway sizes
                max_per_attempt = 200_000
                candidate = current_n + grow
                current_n = min(max(candidate, current_n + 10), max_per_attempt)
            seed += 1
            println("Not enough loci yet (needed $n_loci). Increasing next generation to $current_n and retrying (seed=$seed)")
            attempt += 1
            continue
        end
    end

    # After attempts complete, write aggregated master results back to primary folder_path so
    # downstream steps / users can find the files. Only loci that passed HMBR filtering will be
    # present in these master files.
    hm_out = joinpath(folder_path, "context_harmonic_mean_best_ranks.csv")
    best_out = joinpath(folder_path, "context_best_ranks.csv")
    if nrow(master_hm_df) > 0
        CSV.write(hm_out, master_hm_df)
        println("Aggregated harmonic mean best ranks written to: $hm_out")
    else
        println("Warning: no loci passed HMBR filtering across attempts; no aggregated harmonic mean file written. Expected path would have been: $hm_out")
    end

    if nrow(master_best_df) > 0
        CSV.write(best_out, master_best_df)
        println("Aggregated best ranks written to: $best_out")
    else
        println("Note: no best-rank rows were aggregated; no aggregated best ranks file written. Expected path would have been: $best_out")
    end

    try
        main()
    catch e
        println("ERROR: ", e)
        Base.show_backtrace(stderr, catch_backtrace())
        rethrow(e)
    end
end

main()
