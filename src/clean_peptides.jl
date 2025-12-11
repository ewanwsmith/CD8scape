#!/usr/bin/env julia
"""
clean_peptides.jl

Removes empty lines, stop codons, and duplicate peptides from input .pep file.

Usage:
    julia clean_peptides.jl <folder_path>

Arguments:
    <folder_path>   Path to the folder containing Peptides.pep.
"""

using CSV
using DataFrames

# Function to read, filter, and save unique peptides from Peptides.pep
function process_peptides(folder_path::String)
    # Construct the file path
    files = readdir(folder_path)
    peptide_file = findfirst(f -> lowercase(f) == "peptides.pep", lowercase.(files))
    if peptide_file === nothing
        println("Error: No peptide file (peptides.pep) found in folder.")
        return
    end
    file_path = joinpath(folder_path, files[peptide_file])
    
    # Read the Peptides.pep file into a DataFrame
    df = CSV.read(file_path, DataFrame)
    
    # Filter for unique values in the DataFrame (assuming one column is present)
    df_unique = unique(df)
    
    # Remove peptides containing stop codons (denoted by '*')
    df_filtered = filter(row -> !occursin(r"\*", row[1]), df_unique)
    
    # Save the filtered DataFrame back to the same file
    CSV.write(file_path, df_filtered)
    
    println("Unique peptides without stop codons have been filtered and saved to $(file_path)")
end

# Command-line argument for folder path
function main()
    if length(ARGS) != 1
        println("Usage: julia clean_peptides.jl --folder <folder_path>")
        return
    end
    
    folder_path = ARGS[1]
    process_peptides(folder_path)
end

main()