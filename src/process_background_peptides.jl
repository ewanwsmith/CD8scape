#!/usr/bin/env julia

using DataFrames, CSV
using ArgParse

# Parse command-line arguments
function parse_args()
    parser = ArgParse.ArgumentParser()
    @add_argument parser "--folder" "-f" help="Path to the folder containing the processed_peptides.csv file" required=true
    return parse_args(parser)
end

# Main function to process the file
function main(args)
    folder_path = args["folder"]
    input_file_path = joinpath(folder_path, "processed_peptides.csv")

    # Check if the file exists
    if !isfile(input_file_path)
        println("Error: processed_peptides.csv not found in folder $folder_path")
        return
    end

    # Read the CSV file into a DataFrame
    peptides = CSV.read(input_file_path, DataFrame)

    # Drop the 'Locus' column
    select!(peptides, Not(:Locus))

    # Filter rows where Peptide_label ends with "_A"
    background_peptides = filter(row -> endswith(row[:Peptide_label], "_A"), peptides)

    # Function to extract loci range from Peptide_label
    function extract_loci(peptide_label::String)
        loci_str = split(peptide_label, "_")[1]
        start_end = split(loci_str, "-")
        start = parse(Int, start_end[1])
        stop = parse(Int, start_end[2])
        return collect(start:stop)
    end

    # Apply the extract_loci function to the Peptide_label column and create the Loci column
    background_peptides.Loci = map(extract_loci, background_peptides.Peptide_label)

    # Create a list to store rows
    rows = []

    # Flatten the Loci column and corresponding EL_Rank values into individual rows
    for i in 1:nrow(background_peptides)
        loci_values = background_peptides.Loci[i]
        el_rank_values = repeat([background_peptides.EL_Rank[i]], length(loci_values))  # Repeat the EL_Rank for each locus
        append!(rows, zip(loci_values, el_rank_values, repeat([background_peptides.MHC[i]], length(loci_values))))  # Append locus, EL_Rank, and MHC to the rows list
    end

    # Convert the rows list into a DataFrame
    loci_expanded = DataFrame(rows, [:Loci, :EL_Rank, :MHC])

    # Group by Loci and MHC, then get the minimum EL_Rank for each combination of Loci and MHC
    grouped = groupby(loci_expanded, [:Loci, :MHC])
    min_ranks = combine(grouped, :EL_Rank => minimum)

    # Pivot the DataFrame to a wider format: one column per MHC type for each locus
    wide_df = unstack(min_ranks, :Loci, :MHC, :EL_Rank_minimum)

    # Function to compute the harmonic mean of EL_Rank values across the different MHC types for each locus
    function harmonic_mean(values::Vector{Float64})
        n = length(values)
        return n / sum(1 ./ values)  # Harmonic mean formula
    end

    # Add a new column for the harmonic mean of ranks across all MHC types for each locus
    wide_df.Harmonic_Mean_Best_Rank = map(row -> harmonic_mean(collect(row[2:end])), eachrow(wide_df))

    # Save the result to a CSV file in the same folder
    output_file_path = joinpath(folder_path, "best_ranks_by_locus.csv")
    CSV.write(output_file_path, wide_df)

    println("Output saved to $output_file_path")
end

# Run the script
args = parse_args()
main(args)