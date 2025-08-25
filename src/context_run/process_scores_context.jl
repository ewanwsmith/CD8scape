#!/usr/bin/env julia
"""
This script joins NetMHCpan output with peptide labels for downstream analysis.

Purpose:
- Reads NetMHCpan prediction results and peptide label files.
- Normalizes column names and checks for required columns.
- Joins data on peptide sequence.
- Outputs annotated scores to CSV.

Usage:
    julia process_scores_context.jl <folder_path>

Arguments:
    <folder_path>   Path to the folder containing NetMHCpan and label files.
"""


using DataFrames, CSV

# Print usage and exit if arguments are incorrect
function print_usage_and_exit()
    println("Usage: julia process_scores_context.jl <folder_path>")
    exit(1)
end

if length(ARGS) != 1
    print_usage_and_exit()
end

# Parse input folder and file paths
folder_path = ARGS[1]
output_file = joinpath(folder_path, "context_scores.csv")
netmhc_file = joinpath(folder_path, "context_processed_netMHCpan_output.csv")
labels_file = joinpath(folder_path, "context_peptides_labels.csv")

println("[process_scores_context] Reading NetMHCpan output: $netmhc_file")
println("[process_scores_context] Reading peptide labels: $labels_file")

try
    # Read NetMHCpan output CSV
    netmhc_df = CSV.read(netmhc_file, DataFrame; delim=',')
    # Read peptide labels CSV
    labels_df = CSV.read(labels_file, DataFrame)

    # Normalize column names to remove leading/trailing spaces and unify case
    rename!(netmhc_df, Dict(c => strip(String(c)) for c in names(netmhc_df)))
    rename!(labels_df, Dict(c => strip(String(c)) for c in names(labels_df)))

    # Check required columns in both dataframes
    required_netmhc = ["Peptide"]
    required_labels = ["Peptide", "Peptide_label"]
    for col in required_netmhc
        if col ∉ names(netmhc_df)
            error("Missing column '$col' in NetMHCpan output.")
        end
    end
    for col in required_labels
        if col ∉ names(labels_df)
            error("Missing column '$col' in peptide labels.")
        end
    end

    # Join NetMHCpan output with peptide labels on Peptide column
    joined_df = leftjoin(netmhc_df, labels_df, on=["Peptide"])
    println("[process_scores_context] Joined columns: ", names(joined_df))

    # Write joined DataFrame to output CSV
    CSV.write(output_file, joined_df)
    println("[process_scores_context] Saved scores to $output_file")
catch e
    # Print error and exit if any step fails
    println("[process_scores_context] Error processing scores: $e")
    exit(1)
end
