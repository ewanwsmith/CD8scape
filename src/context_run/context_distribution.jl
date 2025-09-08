#!/usr/bin/env julia
"""
This script compares observed fold changes (HMBR) to a context distribution and calculates percentiles.

Purpose:
- For each observed locus, computes the percentile of its fold change within the context distribution.
- Prints the percentile for each locus.
- Appends the percentile as a new column (foldchange_percentile) in the observed DataFrame and saves the result.

Usage:
    julia context_distribution.jl <folder_path>

Arguments:
    <folder_path>   Path to the folder containing input files.
"""

using CSV, DataFrames, Statistics

# Parse folder path from command line
if length(ARGS) < 1
    println("Usage: julia context_distribution.jl <folder_path>")
    exit(1)
end
folder_path = ARGS[1]

# File paths (relative to folder_path)
context_file = joinpath(folder_path, "context_harmonic_mean_best_ranks.csv")
observed_file = joinpath(folder_path, "harmonic_mean_best_ranks.csv")
output_file = joinpath(folder_path, "harmonic_mean_best_ranks_with_percentile.csv")

# Column name for fold change/HMBR
foldchange_col = "foldchange_HMBR"

# Read context and observed data
context_df = CSV.read(context_file, DataFrame)
observed_df = CSV.read(observed_file, DataFrame)

# Extract context fold changes
context_foldchanges = context_df[:, foldchange_col]

# Function to calculate percentile
function percentile_rank(value, distribution)
    n_below = count(x -> x <= value, distribution)
    return 100 * n_below / length(distribution)
end

# Calculate and append percentiles
percentiles = [percentile_rank(val, context_foldchanges) for val in observed_df[:, foldchange_col]]
observed_df[:, "foldchange_percentile"] = percentiles

# Print results for each locus
for row in eachrow(observed_df)
    locus_desc = haskey(row, "Description") ? row["Description"] : "Unknown"
    println("Locus: $(locus_desc) Fold change = $(row[foldchange_col]), Percentile = $(row["foldchange_percentile"])%")
end

# Save updated DataFrame
CSV.write(output_file, observed_df)
println("Saved updated DataFrame with percentiles to $(output_file)")
