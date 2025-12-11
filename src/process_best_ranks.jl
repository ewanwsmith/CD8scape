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
        raw_desc = String(sdf.Peptide_label[idx])
        # Keep mutation info; remove trailing _C/_V, strip numeric suffix (e.g. _17), and replace spaces with underscores
        cleaned_desc = replace(raw_desc, r"_(C|V)$" => "")
        cleaned_desc = replace(cleaned_desc, r"_\d+$" => "")
        cleaned_desc = replace(cleaned_desc, ' ' => '_')
        (; Best_EL_Rank = sdf.EL_Rank[idx],
           Description = cleaned_desc,
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


# Combine consensus and variant results
best_ranks = vcat(best_C, best_V)

# --- Map Locus to protein Description from frames.csv ---
frames_file = joinpath(folder_path, "frames.csv")
if isfile(frames_file)
    frames_df = CSV.read(frames_file, DataFrame)
    # Extract start positions from Region (assume first number is start)
    function region_to_start(region)
        # Handles both "start,end" and "start1,end1;start2,end2" cases
        split(strip(string(region)), ";")[1] |> x -> split(x, ",")[1] |> y -> parse(Int, y)
    end
    # No longer override Description from frames; retain mutation and protein label from Peptide_label
else
    println("Warning: frames.csv not found in $folder_path. Protein labels will not be mapped.")
end

# Reorder to put Locus then Description first
best_ranks = select(best_ranks, :Locus, :Description, Not([:Locus, :Description]))

# Final tidy on Description: remove any trailing underscores and numeric suffixes
best_ranks.Description = replace.(best_ranks.Description, r"_\d+$" => "")
best_ranks.Description = replace.(best_ranks.Description, r"_+$" => "")

# Compute one Description per Locus: prefer current cleaned, otherwise first seen
description_roots = combine(groupby(best_ranks, :Locus)) do sdf
    (; Locus = sdf.Locus[1], Description = sdf.Description[1])
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
    end

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

    # Calculate fold change (Derived (_V) / Ancestral (_C)) for all valid rows
    if ("HMBR_C" in names(pivot_df)) && ("HMBR_V" in names(pivot_df))
        pivot_df.foldchange_HMBR = pivot_df.HMBR_V ./ pivot_df.HMBR_C
        pivot_df.log2_foldchange_HMBR = log2.(pivot_df.foldchange_HMBR)
    else
        println("DEBUG: Cannot calculate foldchange_HMBR or log2_foldchange_HMBR due to missing columns.")
    end

    # Compute shared Description root per Locus
    pivot_df = leftjoin(pivot_df, description_roots, on = :Locus)

    # Keep Description exactly as mapped; ensure no trailing underscores or numeric suffixes
    if :Description in names(pivot_df)
        pivot_df.Description = replace.(pivot_df.Description, r"_\d+$" => "")
        pivot_df.Description = replace.(pivot_df.Description, r"_+$" => "")
    end

    # Calculate fold change and log2 after all joins/cleaning
    if ("HMBR_C" in names(pivot_df)) && ("HMBR_V" in names(pivot_df))
        pivot_df.foldchange_HMBR = pivot_df.HMBR_V ./ pivot_df.HMBR_C
        pivot_df.log2_foldchange_HMBR = log2.(pivot_df.foldchange_HMBR)
    end

    # Prepare final output columns (use already computed values)
    output_cols = [:Locus, :Description, :HMBR_C, :HMBR_V, :foldchange_HMBR, :log2_foldchange_HMBR]
    pivot_df = select(pivot_df, output_cols...)

    # Save harmonic mean results with fold change
    harmonic_mean_file = joinpath(folder_path, "harmonic_mean_best_ranks.csv") 
    CSV.write(harmonic_mean_file, pivot_df)
    println("Saved harmonic mean best ranks to $harmonic_mean_file")

    # Minimal test: write Locus and log2_foldchange_HMBR to a separate CSV for debugging
    if :log2_foldchange_HMBR in names(pivot_df)
        minimal_test_file = joinpath(folder_path, "harmonic_mean_best_ranks_log2_test.csv")
        minimal_df = select(pivot_df, :Locus, :log2_foldchange_HMBR)
        CSV.write(minimal_test_file, minimal_df)
        println("Minimal test CSV written to $minimal_test_file")
    else
    end
else
    println("No valid best rank data available. Skipping harmonic mean calculations.")
end