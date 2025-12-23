#!/usr/bin/env julia
"""
percentile.jl

Compute observed HMBR log2 fold-change percentiles relative to a simulated distribution.

Usage:
    julia src/percentile.jl <folder_path> [--s <sim_file>] [--o <obs_file>]

Defaults:
- Simulation file defaults to harmonic_mean_best_ranks_simulated.csv in <folder_path>.
- Observed file defaults to the most recent harmonic_mean_best_ranks*.csv in <folder_path>
  whose suffix is not exactly _simulated.

Output:
- Writes percentile_harmonic_mean_best_ranks(_suffix).csv in <folder_path>, adding a
  Percentile column (0-100) to the observed rows.
"""

using CSV, DataFrames, StatsBase

# Utility: resolve file path relative to folder if not absolute
function _resolve_path(folder::AbstractString, p::AbstractString)
    if isempty(p)
        return p
    end
    return isabspath(p) ? String(p) : joinpath(folder, p)
end

# Utility: discover latest observed file (excluding suffix == "simulated")
function _discover_observed(folder::AbstractString)::String
    files = readdir(folder; join=true)
    # match harmonic_mean_best_ranks*.csv and exclude *_simulated.csv
    candidates = String[]
    for f in files
        base = basename(f)
        if endswith(base, ".csv") && startswith(base, "harmonic_mean_best_ranks")
            name, ext = splitext(base)
            # Suffix portion after the base name
            suffix_part = replace(name, "harmonic_mean_best_ranks" => "")
            # Exclude exactly _simulated
            if lowercase(suffix_part) == "_simulated"
                continue
            end
            push!(candidates, f)
        end
    end
    if isempty(candidates)
        error("No observed harmonic_mean_best_ranks*.csv found in $(folder) (excluding _simulated).")
    end
    # Choose latest by modification time
    mtimes = map(p -> stat(p).mtime, candidates)
    idx = argmax(mtimes)
    return candidates[idx]
end

# Utility: extract suffix from observed filename for output naming
function _suffix_from_observed(obspath::AbstractString)::String
    base = basename(obspath)
    name, ext = splitext(base)
    suffix_part = replace(name, "harmonic_mean_best_ranks" => "")
    if isempty(suffix_part)
        return ""
    end
    # Drop leading underscore
    if startswith(suffix_part, "_")
        suffix_part = suffix_part[2:end]
    end
    return suffix_part
end

# Utility: build identity key function based on common columns
function _make_key(df::DataFrame)
    cols = names(df)
    if all(in([:Frame, :Locus, :Mutation]), cols)
        return row -> string(row.Frame, "|", row.Locus, "|", row.Mutation)
    elseif all(in([:Locus, :Mutation]), cols)
        return row -> string(row.Locus, "|", row.Mutation)
    elseif :Locus in cols
        return row -> string(row.Locus)
    else
        # Fallback to row number (no overlap filtering possible)
        return row -> string(ROW_NUMBER[])  # will be replaced inline where used
    end
end

# Main program
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 1
        println("Usage: julia src/percentile.jl <folder_path> [--s <sim_file>] [--o <obs_file>]")
        exit(1)
    end

    folder_path = ARGS[1]
    # Wrap flag parsing in a helper to avoid soft-scope warnings
    function _parse_s_o(argv::Vector{String})
        s = ""; o = ""
        i = 1
        while i <= length(argv)
            a = argv[i]
            if a == "--s"
                if i < length(argv)
                    s = argv[i+1]
                    i += 1
                end
            elseif a == "--o"
                if i < length(argv)
                    o = argv[i+1]
                    i += 1
                end
            end
            i += 1
        end
        return s, o
    end
    sim_file_cli, obs_file_cli = _parse_s_o(ARGS[2:end])

    # Resolve simulation path
    sim_path = isempty(sim_file_cli) ? joinpath(folder_path, "harmonic_mean_best_ranks_simulated.csv") : _resolve_path(folder_path, sim_file_cli)
    if !isfile(sim_path)
        error("Simulation file not found: " * sim_path)
    end

    # Resolve observed path (discover latest if not provided)
    obs_path = isempty(obs_file_cli) ? _discover_observed(folder_path) : _resolve_path(folder_path, obs_file_cli)
    if !isfile(obs_path)
        error("Observed file not found: " * obs_path)
    end

    println("Reading simulation: " * sim_path)
    println("Reading observed:   " * obs_path)

    sim_df = CSV.read(sim_path, DataFrame)
    obs_df = CSV.read(obs_path, DataFrame)

    # Ensure log2_foldchange_HMBR exists; compute from HMBR_D/HMBR_A if needed
    hascol(df, sym) = hasproperty(df, sym)
    function ensure_log2!(df::DataFrame)
        if hascol(df, :log2_foldchange_HMBR)
            return
        end
        if hascol(df, :HMBR_A) && hascol(df, :HMBR_D)
            # Try to coerce to Float64
            a = Vector{Union{Missing,Float64}}(undef, nrow(df))
            d = Vector{Union{Missing,Float64}}(undef, nrow(df))
            for (i, row) in enumerate(eachrow(df))
                ai = row[:HMBR_A]
                di = row[:HMBR_D]
                try
                    a[i] = ai === missing ? missing : Float64(ai)
                catch
                    a[i] = missing
                end
                try
                    d[i] = di === missing ? missing : Float64(di)
                catch
                    d[i] = missing
                end
            end
            fold = Vector{Union{Missing,Float64}}(undef, nrow(df))
            for i in 1:nrow(df)
                if a[i] === missing || d[i] === missing || a[i] == 0
                    fold[i] = missing
                else
                    fold[i] = d[i] / a[i]
                end
            end
            log2fc = Vector{Union{Missing,Float64}}(undef, nrow(df))
            for i in 1:nrow(df)
                if fold[i] === missing
                    log2fc[i] = missing
                else
                    log2fc[i] = log2(fold[i])
                end
            end
            df[!, :log2_foldchange_HMBR] = log2fc
        else
            error("File missing 'log2_foldchange_HMBR' and cannot compute from HMBR_A/HMBR_D.")
        end
    end

    ensure_log2!(sim_df)
    ensure_log2!(obs_df)

    # Overlap removal: drop simulated entries that match observed variants
    # Build keys on shared identifier columns
    has_frame = hascol(sim_df, :Frame) && hascol(obs_df, :Frame)
    has_locus = hascol(sim_df, :Locus) && hascol(obs_df, :Locus)
    has_mut   = hascol(sim_df, :Mutation) && hascol(obs_df, :Mutation)
    keyfun_sim = nothing
    keyfun_obs = nothing
    if has_frame && has_locus && has_mut
        keyfun_sim = row -> string(row.Frame, "|", row.Locus, "|", row.Mutation)
        keyfun_obs = row -> string(row.Frame, "|", row.Locus, "|", row.Mutation)
    elseif has_locus && has_mut
        keyfun_sim = row -> string(row.Locus, "|", row.Mutation)
        keyfun_obs = row -> string(row.Locus, "|", row.Mutation)
    elseif has_locus
        keyfun_sim = row -> string(row.Locus)
        keyfun_obs = row -> string(row.Locus)
    else
        # No stable key available; effectively disable overlap filtering
        keyfun_sim = row -> "sim"
        keyfun_obs = row -> "obs"
    end

    obs_keys = Set{String}()
    for r in eachrow(obs_df)
        push!(obs_keys, keyfun_obs(r))
    end

    keep_mask = trues(nrow(sim_df))
    for (i, r) in enumerate(eachrow(sim_df))
        if keyfun_sim(r) in obs_keys
            keep_mask[i] = false
        end
    end
    sim_df_filtered = sim_df[keep_mask, :]

    # Build simulated distribution from available numeric values
    sim_vals = Vector{Float64}()
    for v in sim_df_filtered[!, :log2_foldchange_HMBR]
        if v isa Missing
            continue
        end
        try
            push!(sim_vals, Float64(v))
        catch
            # skip non-parsable
        end
    end

    if isempty(sim_vals)
        println("Warning: Simulated distribution is empty after filtering; Percentile will be missing.")
        obs_df[!, :Percentile] = [missing for _ in 1:nrow(obs_df)]
    else
        F = ecdf(sim_vals)
        # Percentile as proportion <= x times 100
        perc = Float64[]
        for v in obs_df[!, :log2_foldchange_HMBR]
            if v isa Missing
                push!(perc, NaN)
            else
                try
                    push!(perc, 100.0 * F(Float64(v)))
                catch
                    push!(perc, NaN)
                end
            end
        end
        obs_df[!, :Percentile] = perc
    end

    # Write output next to observed, carrying through observed suffix
    out_suffix = _suffix_from_observed(obs_path)
    out_name = isempty(out_suffix) ? "percentile_harmonic_mean_best_ranks.csv" : string("percentile_harmonic_mean_best_ranks_", out_suffix, ".csv")
    out_path = joinpath(folder_path, out_name)
    CSV.write(out_path, obs_df)
    println("Wrote percentiles to " * out_path)
end
