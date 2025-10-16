#!/usr/bin/env julia
"""
process_best_ranks_supertype.jl

Processes best ranks for peptides using a representative supertype HLA panel.

Usage:
    julia process_best_ranks_supertype.jl <folder_path>

Arguments:
    <folder_path>   Path to the folder containing processed peptides.
"""

using DataFrames, CSV, Statistics, StatsBase

# -----------------------------
# Helpers
# -----------------------------

# Weighted harmonic mean for positive values only; returns `missing` if no valid terms
function hmean_weighted(x::AbstractVector, w::AbstractVector)
    @assert length(x) == length(w) "x and w must have the same length"
    keep = .!ismissing.(x) .& .!ismissing.(w) .& (w .> 0) .& (x .> 0)
    if !any(keep)
        return missing
    end
    xk = Float64.(x[keep])
    wk = Float64.(w[keep])
    return sum(wk) / sum(wk ./ xk)
end

# Normalize allele strings to canonical two-field form: HLA-<LOCUS>*XX:YY
# Handles:
# - missing '*' (e.g., HLA-A01:01 -> HLA-A*01:01)
# - missing ':' with 4 digits (e.g., HLA-A*2402 or HLA-A2402 -> HLA-A*24:02)
# - extra fields (e.g., HLA-A*01:01:01:01 -> HLA-A*01:01)
# - NBSPs/whitespace/case noise
function normalize_allele(s::AbstractString)
    s2 = String(s)
    s2 = replace(s2, '\u00A0' => ' ')            # remove NBSP
    s2 = replace(s2, r"\s+" => "")               # remove all spaces
    s2 = uppercase(s2)

    # Case 1: HLA-A*01:01 (possibly with extra fields) -> take first two fields
    m = match(r"^(HLA-[A-Z]+)\*(\d{2}):(\d{2})", s2)
    if m !== nothing
        return string(m.captures[1], "*", m.captures[2], ":", m.captures[3])
    end

    # Case 2: HLA-A*2402 -> insert colon
    m = match(r"^(HLA-[A-Z]+)\*(\d{4})$", s2)
    if m !== nothing
        g = m.captures[2]
        return string(m.captures[1], "*", g[1:2], ":", g[3:4])
    end

    # Case 3: HLA-A01:01 -> add star
    m = match(r"^(HLA-[A-Z]+)(\d{2}):(\d{2})$", s2)
    if m !== nothing
        return string(m.captures[1], "*", m.captures[2], ":", m.captures[3])
    end

    # Case 4: HLA-A2402 -> add star and colon
    m = match(r"^(HLA-[A-Z]+)(\d{4})$", s2)
    if m !== nothing
        g = m.captures[2]
        return string(m.captures[1], "*", g[1:2], ":", g[3:4])
    end

    return s2
end

# Compute weighted harmonic mean for a subgroup using allele names to fetch weights
function hmean_from_group(best_el_ranks::AbstractVector, alleles::AbstractVector, wmap::Dict{String, Float64})
    w = Vector{Union{Missing, Float64}}(undef, length(alleles))
    @inbounds for i in eachindex(alleles)
        a = normalize_allele(String(alleles[i]))
        w[i] = get(wmap, a, missing)
    end
    return hmean_weighted(best_el_ranks, w)
end

# Clean BOM/whitespace/case and normalize to snake_case
clean_colname(c) = begin
    s = String(c)
    s = replace(s, r"^\ufeff" => "")  # strip BOM
    s = strip(s)
    s = lowercase(s)
    s = replace(s, r"[^a-z0-9]+" => "_")  # non-alnum -> _
    s = replace(s, r"_+" => "_")
    s = strip(s, '_')
    Symbol(s)
end

# -----------------------------
# Args & IO
# -----------------------------

if length(ARGS) != 1
    println("Usage: julia process_best_ranks_supertype.jl <folder_path>")
    exit(1)
end

folder_path = ARGS[1]
input_file = joinpath(folder_path, "processed_peptides.csv")
println("Reading input file: $input_file")

# Read the input CSV into a DataFrame with error handling
try
    global df = CSV.read(input_file, DataFrame)
    println("Successfully loaded data with $(nrow(df)) rows")
catch e
    println("Error reading input file: $e")
    exit(1)
end

# -----------------------------
# Compute best ranks per (Locus, MHC) for C and V peptides
# -----------------------------

function find_best_ranks(df, pattern)
    subset = filter(row -> endswith(row.Peptide_label, pattern), df)
    if isempty(subset)
        println("Warning: No peptides found for pattern '$pattern'")
        return DataFrame(Locus = Int[], MHC = String[], Best_EL_Rank = Float64[], Peptide_Type = String[], Description = String[], Sequence = String[])
    end
    grouped = groupby(subset, [:Locus, :MHC])
    best_rows = combine(grouped) do sdf
        idx = argmin(sdf.EL_Rank)
        (; Best_EL_Rank = sdf.EL_Rank[idx],
           Description = sdf.Peptide_label[idx],
           Sequence = sdf.Peptide[idx])
    end
    return best_rows
end

println("Calculating best ranks for consensus peptides (_C)...")
best_C = find_best_ranks(df, "_C")
best_C.Peptide_Type .= "C"
println("Found best ranks for consensus peptides: $(nrow(best_C)) entries")

println("Calculating best ranks for variant peptides (_V)...")
best_V = find_best_ranks(df, "_V")
best_V.Peptide_Type .= "V"
println("Found best ranks for variant peptides: $(nrow(best_V)) entries")

best_ranks = vcat(best_C, best_V)

# Clean descriptions to remove trailing suffix (e.g., _C, _V, _1)
best_ranks.Description = [replace(strip(string(s)), r"_[^_]*$" => "") for s in best_ranks.Description]

# Reorder to put Locus then Description first
best_ranks = select(best_ranks, :Locus, :Description, Not([:Locus, :Description]))

# Compute shared Description root per Locus
function shared_prefix(strings)
    cleaned = [replace(s, r"_([^_]*)$" => "") for s in strings]
    if isempty(cleaned)
        return ""
    end
    prefix = cleaned[1]
    for s in cleaned[2:end]
        minlen = min(length(prefix), length(s))
        i = 1
        while i <= minlen && prefix[i] == s[i]
            i += 1
        end
        prefix = prefix[1:i-1]
        if isempty(prefix)
            break
        end
    end
    prefix = replace(prefix, r"_$" => "")
    return prefix
end

description_roots = combine(groupby(best_ranks, :Locus)) do sdf
    (; Locus = sdf.Locus[1], Description = shared_prefix(sdf.Description))
end

# Save best_ranks.csv
best_ranks_file = joinpath(folder_path, "best_ranks.csv")
CSV.write(best_ranks_file, best_ranks)
println("Saved best ranks to $best_ranks_file")

# -----------------------------
# Weighted HMBR using supertype_panel.csv
# -----------------------------

println("Calculating weighted harmonic mean best ranks (HMBR) for each locus...")
if !isempty(best_ranks)
    freq_path = joinpath(folder_path, "supertype_panel.csv")
    if !isfile(freq_path)
        println("Error: Could not find supertype_panel.csv in $folder_path.")
        exit(1)
    end
    println("Using allele frequency file: $(basename(freq_path))")

    freq_df = CSV.read(freq_path, DataFrame)
    rename!(freq_df, Dict(names(freq_df) .=> clean_colname.(names(freq_df))))

    # Build a lowercase string lookup of current column names
    name_strs = string.(names(freq_df))
    lower_lookup = Dict(lowercase(s) => s for s in name_strs)

    # If both 'allele' and 'frequency' exist, use them directly
    if haskey(lower_lookup, "allele") && haskey(lower_lookup, "frequency")
        allele_sym = Symbol(lower_lookup["allele"])
        freq_sym   = Symbol(lower_lookup["frequency"])
        freq_df = select(freq_df, allele_sym => :allele, freq_sym => :frequency)
        freq_df.frequency = Float64.(freq_df.frequency)
        freq_df.allele = normalize_allele.(string.(freq_df.allele))
        wmap = Dict(String(a) => f for (a, f) in zip(freq_df.allele, freq_df.frequency))

    else
        # Candidates to recognize columns
        allele_candidates = [:allele, :allele_name, :mhc, :mhc_allele]
        freq_candidates   = [:frequency, :freq, :pop_freq, :allele_frequency, :weight, :w]

        allele_col = findfirst(c -> c in names(freq_df), allele_candidates)
        freq_col   = findfirst(c -> c in names(freq_df), freq_candidates)

        # Heuristic fallback if headers are unexpected
        if allele_col === nothing || freq_col === nothing
            cols = names(freq_df)
            # Try to guess allele column: string-like with many values matching allele patterns
            function looks_like_allele(v)
                if !(eltype(v) <: AbstractString)
                    return false
                end
                n = min(50, length(v))
                m = 0
                for i in 1:n
                    s = String(v[i])
                    if occursin(r"[A-Za-z]+[-:]\d", s)
                        m += 1
                    end
                end
                return m >= max(5, ceil(Int, 0.2n))
            end
            allele_guess_idx = findfirst(i -> looks_like_allele(freq_df[!, cols[i]]), eachindex(cols))

            # Try to guess frequency column: numeric and mostly between 0 and 1
            function looks_like_freq(v)
                T = eltype(v)
                if T <: Real
                    n = min(50, length(v))
                    if n == 0
                        return false
                    end
                    vals = collect(skipmissing(v[1:n]))
                    return !isempty(vals) && mean((0 .<= vals) .& (vals .<= 1)) >= 0.7
                else
                    n = min(50, length(v))
                    parsed = Float64[]
                    for i in 1:n
                        x = tryparse(Float64, String(v[i]))
                        if x !== nothing
                            push!(parsed, x)
                        end
                    end
                    return !isempty(parsed) && mean((0 .<= parsed) .& (parsed .<= 1)) >= 0.7
                end
            end
            freq_guess_idx = findfirst(i -> looks_like_freq(freq_df[!, cols[i]]), eachindex(cols))

            if allele_col === nothing && allele_guess_idx !== nothing
                allele_sym = cols[allele_guess_idx]
            end
            if freq_col === nothing && freq_guess_idx !== nothing
                freq_sym = cols[freq_guess_idx]
            end
        end

        if (allele_col === nothing && !(@isdefined allele_sym)) || (freq_col === nothing && !(@isdefined freq_sym))
            println("Error: could not identify allele or frequency columns in $(basename(freq_path)). Available columns: ", join(string.(names(freq_df)), ", "))
            exit(1)
        end

        # Map to :allele and :frequency
        if !(@isdefined allele_sym)
            allele_sym = allele_candidates[allele_col]
        end
        if !(@isdefined freq_sym)
            freq_sym   = freq_candidates[freq_col]
        end

        freq_df = select(freq_df, allele_sym => :allele, freq_sym => :frequency)
        freq_df.frequency = Float64.(freq_df.frequency)
        freq_df.allele = normalize_allele.(string.(freq_df.allele))
        wmap = Dict(String(a) => f for (a, f) in zip(freq_df.allele, freq_df.frequency))
    end

    # === DEBUG WEIGHTS START (remove this whole block when done debugging) ===
    try
        raw_preds = unique(String.(best_ranks.MHC))
        norm_preds = unique(normalize_allele.(raw_preds))
        weight_keys  = collect(keys(wmap))

        # Show a few examples of raw -> normalized
        nshow = min(5, length(raw_preds))
        sample_pairs = [ raw_preds[i] * " -> " * normalize_allele(raw_preds[i]) for i in 1:nshow ]
        println("DEBUG: Sample predicted alleles (raw -> normalized): ", sample_pairs)

        # Show a few sample keys from the weight map
        println("DEBUG: Sample weight keys: ", first(weight_keys, min(5, length(weight_keys))))

        missing_in_panel = sort(setdiff(norm_preds, weight_keys))
        matched = length(norm_preds) - length(missing_in_panel)
        println("DEBUG weights: total alleles used (normalized) = $(length(norm_preds)), matched = $matched, missing = $(length(missing_in_panel))")

        # Write unmatched alleles; keep columns the same length
        debug_file = joinpath(folder_path, "debug_unmatched_alleles.csv")
        CSV.write(debug_file, DataFrame(
            Predicted_Allele = missing_in_panel,
            Closest_Match_Key = fill("", length(missing_in_panel))  # placeholder
        ))
        println("DEBUG weights: wrote unmatched alleles to $debug_file")
    catch e
        println("DEBUG weights: failed to compute/report unmatched alleles: ", e)
    end
    # === DEBUG WEIGHTS END ===

    # Compute weighted harmonic mean per (Locus, Peptide_Type) using per-row allele in best_ranks.MHC
    grouped = groupby(best_ranks, [:Locus, :Peptide_Type])
    agg = combine(grouped) do sdf
        (; HMBR = hmean_from_group(Float64.(sdf.Best_EL_Rank), sdf.MHC, wmap))
    end

    # Pivot to separate C and V (HMBR only)
    pivot_df = unstack(agg, :Peptide_Type, :HMBR)
    rename!(pivot_df, Dict("C" => "HMBR_C", "V" => "HMBR_V"))

    # Identify and report missing values before fold change calculation
    for row in eachrow(pivot_df)
        if ismissing(row.HMBR_C)
            println("Fold change could not be calculated for locus $(row.Locus) due to missing consensus rank or weights.")
        elseif ismissing(row.HMBR_V)
            println("Fold change could not be calculated for locus $(row.Locus) due to missing variant rank or weights.")
        end
    end

    before_filter = nrow(pivot_df)
    pivot_df = filter(row ->
        !ismissing(row.HMBR_C) && !ismissing(row.HMBR_V) &&
        !(row.HMBR_C > 2 && row.HMBR_V > 2),
        pivot_df
    )
    removed_count = before_filter - nrow(pivot_df)
    println("Removed $removed_count loci where both ancestral and derived states were predicted to be non-binding (HMBR > 2)")

    # Fold change V/C
    pivot_df.foldchange_HMBR = pivot_df.HMBR_V ./ pivot_df.HMBR_C

    # Add log2-transformed fold change (safe handling for missing/non-positive values)
    function safe_log2(x)
        if ismissing(x)
            return missing
        end
        xf = tryparse(Float64, string(x))
        if xf === nothing
            return missing
        end
        if xf <= 0
            return missing
        end
        return log2(xf)
    end
    pivot_df.foldchange_HMBR_log2 = [safe_log2(v) for v in pivot_df.foldchange_HMBR]

    # Attach Description and reorder columns (keep original output schema)
    pivot_df = leftjoin(pivot_df, description_roots, on = :Locus)
    pivot_df = select(pivot_df, :Locus, :Description, :HMBR_C, :HMBR_V, :foldchange_HMBR, :foldchange_HMBR_log2)

    # Save results
    harmonic_mean_file = joinpath(folder_path, "harmonic_mean_best_ranks.csv")
    CSV.write(harmonic_mean_file, pivot_df)
    println("Saved weighted harmonic mean best ranks with fold change to $harmonic_mean_file")
else
    println("No valid best rank data available. Skipping harmonic mean calculations.")
end