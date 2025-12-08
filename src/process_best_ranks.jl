#!/usr/bin/env julia
"""
process_best_ranks.jl

Processes best ranks for consensus and variant peptides for each locus.

Usage:
    julia process_best_ranks.jl <folder_path>

Arguments:
    <folder_path>   Path to the folder containing processed_peptides.csv.
"""

using DataFrames, CSV, Statistics, StatsBase

# Parse command-line argument for input folder
if length(ARGS) != 1
    println("Usage: julia script.jl <folder_path>")
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

# Function to find minimum EL_Rank by Locus and MHC for given peptide pattern
function find_best_ranks(df, pattern)
    subset = filter(row -> !ismissing(row.Peptide_label) && endswith(row.Peptide_label, pattern), df)
    if isempty(subset)
        println("Warning: No peptides found for pattern '$pattern'")
        return DataFrame(Locus = Int[], MHC = String[], Best_EL_Rank = Float64[], Peptide_Type = String[], Description = String[], Sequence = String[])
    end
    # Coerce EL_Rank column to Float64 if not already
    if !(eltype(subset.EL_Rank) <: AbstractFloat)
        subset.EL_Rank = map(x -> try
            x isa Number ? float(x) : parse(Float64, String(x))
        catch
            NaN
        end, subset.EL_Rank)
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

# Find best ranks separately for consensus (_C) and variant (_V) peptides
println("Calculating best ranks for consensus peptides (_C)...")
best_C = find_best_ranks(df, "_C")
best_C.Peptide_Type .= "C"
println("Found best ranks for consensus peptides: $(nrow(best_C)) entries")

println("Calculating best ranks for variant peptides (_V)...")
best_V = find_best_ranks(df, "_V")
best_V.Peptide_Type .= "V"
println("Found best ranks for variant peptides: $(nrow(best_V)) entries")

# Combine consensus and variant results and save as best_ranks.csv
best_ranks = vcat(best_C, best_V)

# Clean descriptions to remove suffix (e.g., _C, _V, _1, etc.)
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

# Pivot best_ranks to have separate columns for HMBR_C and HMBR_V
println("Calculating harmonic mean best ranks (HMBR) for each locus...")
if !isempty(best_ranks)
    pivot_df = unstack(combine(groupby(best_ranks, [:Locus, :Peptide_Type]),
        :Best_EL_Rank => harmmean => :HMBR), :Peptide_Type, :HMBR)
    # Only rename columns if they exist
    colnames = names(pivot_df)
    rename_pairs = Pair{Symbol,Symbol}[]
    if "C" in colnames
        push!(rename_pairs, Symbol("C") => :HMBR_C)
    end
    if "V" in colnames
        push!(rename_pairs, Symbol("V") => :HMBR_V)
    end
    if !isempty(rename_pairs)
        rename!(pivot_df, rename_pairs...)
    end

    # Warn if no variant data is present
    if !("HMBR_V" in names(pivot_df))
        println("Warning: No variant peptides found. Skipping fold change calculations and HMBR_V output.")
    else
        # Identify and report missing values before fold change calculation
        for row in eachrow(pivot_df)
            if ismissing(row.HMBR_C)
                println("Fold change could not be calculated for locus $(row.Locus) due to missing consensus rank.")
            elseif ismissing(row.HMBR_V)
                println("Fold change could not be calculated for locus $(row.Locus) due to missing variant rank.")
            end
        end

        # Filter out loci where both HMBR_C and HMBR_V are greater than 2 (non-binding)
        before_filter = nrow(pivot_df)
        pivot_df = filter(row -> !ismissing(row.HMBR_C) && !ismissing(row.HMBR_V) && !(row.HMBR_C > 2 && row.HMBR_V > 2), pivot_df)
        removed_count = before_filter - nrow(pivot_df)
        println("Removed $removed_count loci where both ancestral and derived states were predicted to be non-binding (HMBR > 2)")

        # Calculate fold change (Derived (_V) / Ancestral (_C))
        pivot_df.foldchange_HMBR = pivot_df.HMBR_V ./ pivot_df.HMBR_C
    end

    # Compute shared Description root per Locus
    pivot_df = leftjoin(pivot_df, description_roots, on = :Locus)

    # Reorder to put Locus then Description first
    pivot_df = select(pivot_df, :Locus, :Description, Not([:Locus, :Description]))

    # Save harmonic mean results with fold change
    harmonic_mean_file = joinpath(folder_path, "harmonic_mean_best_ranks.csv") 
    CSV.write(harmonic_mean_file, pivot_df)
    println("Saved harmonic mean best ranks to $harmonic_mean_file")
else
    println("No valid best rank data available. Skipping harmonic mean calculations.")
end