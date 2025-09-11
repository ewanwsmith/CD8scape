#!/usr/bin/env julia
"""
This script processes best ranks for context peptides, including weighted supertype calculations.

Purpose:
- Reads context scores and identifies best ranks for consensus and variant peptides.
- Calculates harmonic mean best ranks (HMBR) for each locus.
- Supports weighted calculation using allele frequencies for supertype panel runs.
- Outputs best ranks and harmonic mean results to CSV.

Usage:
    julia process_best_ranks_context.jl <folder_path> [--supertype]

Arguments:
    <folder_path>   Path to the folder containing context_scores.csv.
    --supertype     Optional flag to use weighted allele frequencies.
"""

using DataFrames, CSV, Statistics, StatsBase

# Parse command-line arguments for input folder and optional --supertype flag
if length(ARGS) < 1 || length(ARGS) > 2
    println("Usage: julia process_best_ranks_context.jl <folder_path> [--supertype]")
    exit(1)
end

folder_path = ARGS[1]
use_supertype = length(ARGS) == 2 && ARGS[2] == "--supertype"
input_file = joinpath(folder_path, "context_scores.csv")
println("Reading input file: $input_file")

function find_best_ranks(df, suffix)
    # Filter and select best ranks, including HLA
    subdf = filter(row -> !ismissing(row.Peptide_label) && endswith(row.Peptide_label, suffix), df)
    grouped = groupby(subdf, [:Locus, :Peptide_label, :HLA])
    best = combine(grouped, :EL_Rank => minimum => :Best_EL_Rank)
    return best
end

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

function hmean_weighted(x::AbstractVector, w::AbstractVector)
    # Weighted harmonic mean
    return sum(w) / sum(w ./ x)
end

# Robust allele normalization (copied from supertype method)
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

try
    df = CSV.read(input_file, DataFrame; delim=',')
    println("Successfully loaded data with $(nrow(df)) rows")

    # Ensure best_ranks includes Allele column for weighted calculation
    best_C = find_best_ranks(df, "_C")
    best_C.Peptide_Type .= "C"
    println("Found best ranks for consensus peptides: $(nrow(best_C)) entries")

    best_V = find_best_ranks(df, "_V")
    best_V.Peptide_Type .= "V"
    println("Found best ranks for variant peptides: $(nrow(best_V)) entries")

    best_ranks = vcat(best_C, best_V)
    best_ranks.Peptide_label = [replace(strip(string(s)), r"_[^_]*$" => "") for s in best_ranks.Peptide_label]
    # Ensure HLA is included, then rename to Allele for matching
    rename!(best_ranks, Dict(:HLA => :Allele))
    best_ranks = select(best_ranks, :Locus, :Peptide_label, :Allele, :Best_EL_Rank, :Peptide_Type, Not([:Locus, :Peptide_label, :Allele, :Best_EL_Rank, :Peptide_Type]))

    # Ensure description_roots is always defined before harmonic mean calculations
    description_roots = combine(groupby(best_ranks, :Locus)) do sdf
        (; Locus = sdf.Locus[1], Description = shared_prefix(sdf.Peptide_label))
    end

    println("Calculating harmonic mean best ranks (HMBR) for each locus...")
    if !isempty(best_ranks)
        if use_supertype
            freq_path = joinpath(folder_path, "supertype_panel.csv")
            freq_df = CSV.read(freq_path, DataFrame)
            freq_df.allele = normalize_allele.(string.(freq_df.Allele))
            freq_df.frequency = Float64.(freq_df.Frequency)
            wmap = Dict(String(a) => f for (a, f) in zip(freq_df.allele, freq_df.frequency))

            grouped = groupby(best_ranks, [:Locus, :Peptide_Type])
            agg = combine(grouped) do sdf
                alleles = normalize_allele.(string.(sdf.Allele))
                ranks = Float64.(sdf.Best_EL_Rank)
                weights = [get(wmap, String(a), missing) for a in alleles]
                (; HMBR = hmean_weighted(ranks, weights))
            end
            pivot_df = unstack(agg, :Peptide_Type, :HMBR)
            rename!(pivot_df, Dict("C" => "HMBR_C", "V" => "HMBR_V"))
            pivot_df.foldchange_HMBR = pivot_df.HMBR_V ./ pivot_df.HMBR_C
            pivot_df = leftjoin(pivot_df, description_roots, on = :Locus)
            pivot_df = select(pivot_df, :Locus, :Description, Not([:Locus, :Description]))
            harmonic_mean_file = joinpath(folder_path, "context_harmonic_mean_best_ranks.csv")
            CSV.write(harmonic_mean_file, pivot_df)
            println("Saved weighted harmonic mean best ranks with fold change to $harmonic_mean_file")
        else
            # Standard (unweighted) calculation
            pivot_df = unstack(combine(groupby(best_ranks, [:Locus, :Peptide_Type]),
                :Best_EL_Rank => harmmean => :HMBR), :Peptide_Type, :HMBR)
            rename!(pivot_df, Dict("C" => "HMBR_C", "V" => "HMBR_V"))
            pivot_df.foldchange_HMBR = pivot_df.HMBR_V ./ pivot_df.HMBR_C
            pivot_df = leftjoin(pivot_df, description_roots, on = :Locus)
            pivot_df = select(pivot_df, :Locus, :Description, Not([:Locus, :Description]))
            harmonic_mean_file = joinpath(folder_path, "context_harmonic_mean_best_ranks.csv")
            CSV.write(harmonic_mean_file, pivot_df)
            println("Saved harmonic mean best ranks with fold change to $harmonic_mean_file")
        end
    else
        println("No valid best rank data available. Skipping harmonic mean calculations.")
    end

    best_ranks_file = joinpath(folder_path, "context_best_ranks.csv")
    CSV.write(best_ranks_file, best_ranks)
    println("Saved best ranks to $best_ranks_file")
catch e
    println("Error reading input file: $e")
    exit(1)
end

# All code using best_ranks is now inside the try block above. No code below references best_ranks.
