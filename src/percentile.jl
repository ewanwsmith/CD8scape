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


With --per-allele:
- Simulation file defaults to per_allele_best_ranks_simulated.csv in <folder_path>.
- Observed file defaults to the most recent per_allele_best_ranks*.csv in <folder_path>
  whose suffix is not exactly _simulated.
- Computes percentiles of log2_foldchange_BR per (Locus, Mutation, MHC).
- Output: percentile_per_allele_best_ranks(_suffix).csv
"""

using CSV, DataFrames, StatsBase

# Utility: resolve file path relative to folder if not absolute
function _resolve_path(folder::AbstractString, p::AbstractString)
    if isempty(p)
        return p
    end
    return isabspath(p) ? String(p) : joinpath(folder, p)
end

# Utility: discover latest observed file matching base_name (excluding suffix == "simulated")
function _discover_observed(folder::AbstractString, base_name::AbstractString)::String
    files = readdir(folder; join=true)
    candidates = String[]
    for f in files
        base = basename(f)
        if endswith(base, ".csv") && startswith(base, base_name)
            name, ext = splitext(base)
            suffix_part = replace(name, base_name => "")
            if lowercase(suffix_part) == "_simulated"
                continue
            end
            push!(candidates, f)
        end
    end
    if isempty(candidates)
        error("No observed $(base_name)*.csv found in $(folder) (excluding _simulated).")
    end
    mtimes = map(p -> stat(p).mtime, candidates)
    idx = argmax(mtimes)
    return candidates[idx]
end

# Utility: extract suffix from observed filename for output naming
function _suffix_from_observed(obspath::AbstractString, base_name::AbstractString)::String
    base = basename(obspath)
    name, ext = splitext(base)
    suffix_part = replace(name, base_name => "")
    if isempty(suffix_part)
        return ""
    end
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

# Stouffer helpers: probit (Abramowitz & Stegun) and normcdf (rational approx)
function _probit(p::Float64)::Float64
    p <= 0.0 && return -Inf; p >= 1.0 && return Inf
    p < 0.5 && return -_probit(1.0 - p)
    t = sqrt(-2.0 * log(1.0 - p))
    c = (2.515517, 0.802853, 0.010328); d = (1.432788, 0.189269, 0.001308)
    return t - (c[1] + c[2]*t + c[3]*t^2) / (1.0 + d[1]*t + d[2]*t^2 + d[3]*t^3)
end
function _normcdf(z::Float64)::Float64
    z < 0.0 && return 1.0 - _normcdf(-z)
    t = 1.0 / (1.0 + 0.2316419 * z)
    poly = t*(0.319381530 + t*(-0.356563782 + t*(1.781477937 + t*(-1.821255978 + t*1.330274429))))
    return 1.0 - exp(-z^2/2) / sqrt(2π) * poly
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
    per_allele = any(a -> a == "--per-allele", ARGS[2:end])

    # Choose base names and fold-change column based on mode
    if per_allele
        base_name    = "per_allele_best_ranks"
        fc_col       = :log2_foldchange_BR
        out_prefix   = "percentile_per_allele_best_ranks"
    else
        base_name    = "harmonic_mean_best_ranks"
        fc_col       = :log2_foldchange_HMBR
        out_prefix   = "percentile_harmonic_mean_best_ranks"
    end

    # Resolve simulation path
    sim_path = isempty(sim_file_cli) ? joinpath(folder_path, base_name * "_simulated.csv") : _resolve_path(folder_path, sim_file_cli)
    if !isfile(sim_path)
        error("Simulation file not found: " * sim_path)
    end

    # Resolve observed path (discover latest if not provided)
    obs_path = isempty(obs_file_cli) ? _discover_observed(folder_path, base_name) : _resolve_path(folder_path, obs_file_cli)
    if !isfile(obs_path)
        error("Observed file not found: " * obs_path)
    end

    println("Reading simulation: " * sim_path)
    println("Reading observed:   " * obs_path)

    sim_df = CSV.read(sim_path, DataFrame)
    obs_df = CSV.read(obs_path, DataFrame)

    # Ensure the fold-change column exists; for HMBR mode, recompute from A/D if missing
    hascol(df, sym) = hasproperty(df, sym)
    function ensure_log2!(df::DataFrame, col::Symbol)
        if hascol(df, col)
            return
        end
        if col == :log2_foldchange_HMBR && hascol(df, :HMBR_A) && hascol(df, :HMBR_D)
            a = Vector{Union{Missing,Float64}}(undef, nrow(df))
            d = Vector{Union{Missing,Float64}}(undef, nrow(df))
            for (i, row) in enumerate(eachrow(df))
                ai = row[:HMBR_A]; di = row[:HMBR_D]
                try; a[i] = ai === missing ? missing : Float64(ai); catch; a[i] = missing; end
                try; d[i] = di === missing ? missing : Float64(di); catch; d[i] = missing; end
            end
            fold = [a[i] === missing || d[i] === missing || a[i] == 0 ? missing : d[i] / a[i] for i in 1:nrow(df)]
            df[!, :log2_foldchange_HMBR] = [f === missing ? missing : log2(f) for f in fold]
        elseif col == :log2_foldchange_BR && hascol(df, :ELBR_A) && hascol(df, :ELBR_D)
            fold = [a === missing || d === missing || a == 0 ? missing : d / a
                    for (a, d) in zip(df[!, :ELBR_A], df[!, :ELBR_D])]
            df[!, :log2_foldchange_BR] = [f === missing ? missing : log2(f) for f in fold]
        else
            error("File missing '$(col)' and cannot recompute it from available columns.")
        end
    end

    ensure_log2!(sim_df, fc_col)
    ensure_log2!(obs_df, fc_col)

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
    for v in sim_df_filtered[!, fc_col]
        if v isa Missing; continue; end
        try; push!(sim_vals, Float64(v)); catch; end
    end

    if isempty(sim_vals)
        println("Warning: Simulated distribution is empty after filtering; Percentile will be missing.")
        obs_df[!, :Percentile] = [missing for _ in 1:nrow(obs_df)]
    else
        F = ecdf(sim_vals)
        perc = Float64[]
        for v in obs_df[!, fc_col]
            if v isa Missing
                push!(perc, NaN)
            else
                try; push!(perc, 100.0 * F(Float64(v))); catch; push!(perc, NaN); end
            end
        end
        obs_df[!, :Percentile] = perc
    end

    # Stouffer's combined Z: Z_i per variant, COMBINED summary row
    zi_col = Union{Float64, Missing}[]
    valid_z = Float64[]; valid_perc = Float64[]
    for p in obs_df[!, :Percentile]
        if p isa Missing || (p isa Number && (!isfinite(p) || isnan(p)))
            push!(zi_col, missing)
        else
            z = _probit(clamp(Float64(p), 0.01, 99.99) / 100.0)
            push!(zi_col, z)
            push!(valid_z, z)
            push!(valid_perc, Float64(p))
        end
    end
    obs_df[!, :Z_i]     = zi_col
    obs_df[!, :p_value] = Union{Float64, Missing}[missing for _ in 1:nrow(obs_df)]

    if !isempty(valid_z)
        N_z        = length(valid_z)
        Z_combined = sum(valid_z) / sqrt(N_z)
        p_combined = 1.0 - _normcdf(Z_combined)
        combined_vals = Dict{Symbol, Any}(Symbol(col) => missing for col in names(obs_df))
        combined_vals[:Percentile] = sum(valid_perc) / N_z
        combined_vals[:Z_i]        = Z_combined
        combined_vals[:p_value]    = p_combined
        if haskey(combined_vals, :Mutation)
            combined_vals[:Mutation] = "combined_z"
        end
        combined_row = DataFrame(Dict(k => [v] for (k, v) in combined_vals))
        obs_df = vcat(obs_df, combined_row, cols=:union)
        println("Stouffer Z=$(round(Z_combined, digits=4))  p=$(round(p_combined, digits=6))  N=$(N_z)")
    end

    # Empirical permutation test: sample 9999 sets of k from sim_vals, compare mean percentile
    n_perm = 9999
    if !isempty(valid_perc) && !isempty(sim_vals)
        k_obs = length(valid_perc)
        obs_mean_perc = sum(valid_perc) / k_obs
        if length(sim_vals) < k_obs
            println("Warning: sim pool ($(length(sim_vals))) < k ($(k_obs)); skipping empirical_p.")
        else
            local exceed_count = 0
            for _ in 1:n_perm
                samp = sample(sim_vals, k_obs; replace=false)
                samp_mean_perc = sum(100.0 * F(v) for v in samp) / k_obs
                if samp_mean_perc >= obs_mean_perc
                    exceed_count += 1
                end
            end
            p_empirical = exceed_count / n_perm
            emp_vals = Dict{Symbol, Any}(Symbol(col) => missing for col in names(obs_df))
            emp_vals[:Percentile] = obs_mean_perc
            emp_vals[:p_value]    = p_empirical
            if haskey(emp_vals, :Mutation)
                emp_vals[:Mutation] = "empirical_p"
            end
            emp_row = DataFrame(Dict(kk => [vv] for (kk, vv) in emp_vals))
            obs_df = vcat(obs_df, emp_row, cols=:union)
            println("Empirical p=$(round(p_empirical, digits=6))  k=$(k_obs)  N_sim=$(length(sim_vals))  n_perm=$(n_perm)")
        end
    end

    # Write output next to observed, carrying through observed suffix
    out_suffix = _suffix_from_observed(obs_path, base_name)
    out_name = isempty(out_suffix) ? string(out_prefix, ".csv") : string(out_prefix, "_", out_suffix, ".csv")
    out_path = joinpath(folder_path, out_name)
    CSV.write(out_path, obs_df)
    println("Wrote percentiles to " * out_path)
end
