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


# Check if observed file exists
if !isfile(observed_file)
    println("Observed file $(observed_file) missing. Skipping distribution step.")
    exit(0)
end

# Check if context file exists
if !isfile(context_file)
    println("Context distribution file $(context_file) missing. Skipping percentile calculation.")
    observed_df = CSV.read(observed_file, DataFrame)
    CSV.write(output_file, observed_df)
    println("Saved simulated output to $(output_file)")
    exit(0)
end

# Read context and observed data
context_df = CSV.read(context_file, DataFrame)
observed_df = CSV.read(observed_file, DataFrame)

# Extract context fold changes and transform with log2
function safe_log2(x)
    if ismissing(x)
        return missing
    end
    xf = tryparse(Float64, string(x))
    if xf === nothing
        return missing
    end
    if xf <= 0
        println("Warning: non-positive fold change encountered (value=$(x)); treating as missing in log2 transform.")
        return missing
    end
    return log2(xf)
end

context_foldchanges_raw = context_df[:, foldchange_col]
context_foldchanges = [safe_log2(v) for v in context_foldchanges_raw]

# If context is empty, skip percentile step and output observed_df as simulated output
if isempty(context_foldchanges)
    println("Context distribution empty. Skipping percentile calculation.")
    CSV.write(output_file, observed_df)
    println("Saved simulated output to $(output_file)")
    exit(0)
end

# Function to calculate percentile
# Percentiles are computed on the log2-transformed fold changes (inclusive definition)
function percentile_rank_log2(value_raw, distribution_log2)
    # transform observed value
    val = safe_log2(value_raw)
    vals = collect(skipmissing(distribution_log2))
    if isempty(vals) || ismissing(val)
        return NaN
    end
    n_below = count(x -> x <= val, vals)
    return 100 * n_below / length(vals)
end

# Calculate and append percentiles
percentiles = [percentile_rank_log2(val, context_foldchanges) for val in observed_df[:, foldchange_col]]
observed_df[:, "foldchange_percentile"] = percentiles

# Print results for each locus (show raw and log2)
for row in eachrow(observed_df)
    locus_desc = haskey(row, "Description") ? row["Description"] : "Unknown"
    raw = row[foldchange_col]
    lg = safe_log2(raw)
    pct = row["foldchange_percentile"]
    println("Locus: $(locus_desc) Fold change = $(raw), log2 = $(lg), Percentile = $(pct)%")
end

# Save updated DataFrame
CSV.write(output_file, observed_df)
println("Saved updated DataFrame with percentiles to $(output_file)")
