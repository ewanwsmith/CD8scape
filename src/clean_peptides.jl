#!/usr/bin/env julia

using CSV
using DataFrames

# Function to read, filter, and save unique peptides from Peptides.pep
function process_peptides(folder_path::String)
    # Construct the file path
    file_path = joinpath(folder_path, "Peptides.pep")
    
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
        println("Usage: julia process_peptides.jl --folder <folder_path>")
        return
    end
    
    folder_path = ARGS[1]
    process_peptides(folder_path)
end

main()