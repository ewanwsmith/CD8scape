#!/usr/bin/env julia
"""
This script cleans context peptides for downstream NetMHCpan analysis.

Purpose:
- Reads raw context peptide file.
- Removes empty lines, stop codons, and duplicates.
- Outputs cleaned peptides to CSV.

Usage:
    julia clean_peptides_context.jl <folder_path>

Arguments:
    <folder_path>   Path to the folder containing context_peptides.pep.
"""

using DataFrames, CSV

function print_usage_and_exit()
    println("Usage: julia clean_peptides_context.jl <folder_path>")
    exit(1)
end

if length(ARGS) != 1
    print_usage_and_exit()
end

folder_path = ARGS[1]
input_file = joinpath(folder_path, "context_peptides.pep")
output_file = joinpath(folder_path, "context_cleaned_peptides.csv")

println("[clean_peptides_context] Reading input file: $input_file")

try
    lines = readlines(input_file)
    cleaned = filter(line -> !isempty(strip(line)) && !occursin("*", line), lines) # Example: remove empty and stop codon lines
    cleaned = unique(cleaned) # Remove duplicates
    open(output_file, "w") do f
        for line in cleaned
            write(f, line * "\n")
        end
    end
    println("[clean_peptides_context] Saved cleaned peptides to $output_file")
catch e
    println("[clean_peptides_context] Error cleaning peptides: $e")
    exit(1)
end

# If output DataFrame or CSV uses Description_Root, rename to Description
