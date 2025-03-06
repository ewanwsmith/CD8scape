#!/usr/bin/env julia

using DataFrames, CSV

# Manual argument parsing function
function parse_arguments()
    args = Dict()
    for (i, arg) in enumerate(ARGS)
        if arg in ["--folder", "-f"]
            args["folder"] = ARGS[i + 1]
        end
    end
    if !haskey(args, "folder")
        error("Usage: ./process_scores.jl --folder /path/to/data")
    end
    return args
end

# Process and join files
function process_and_join(folder_path::String)::DataFrame
    mhcpan_path = joinpath(folder_path, "processed_output.csv")
    peptides_path = joinpath(folder_path, "peptides_labels.csv")

    if !(isfile(mhcpan_path) && isfile(peptides_path))
        error("Required files not found. Ensure both processed_output.csv and peptides_labels.csv are present in $folder_path.")
    end

    mhcpan_df = CSV.read(mhcpan_path, DataFrame)
    peptides_df = CSV.read(peptides_path, DataFrame)

    joined_df = leftjoin(peptides_df, mhcpan_df, on="Peptide")
    dropped_rows = joined_df[ismissing.(joined_df.Pos), :]

    if "Locus" in names(dropped_rows)
        println("Dropped rows (Peptide and Locus):")
        println(dropped_rows[:, [:Peptide, :Locus]])
    else
        println("Column 'Locus' not found in dropped rows.")
    end

    return joined_df
end

# Reshape HLA-related data
function reshape_hla_data(df::DataFrame)::DataFrame
    df[:, :MHC] = df[:, :HLA]
    select!(df, Not(:HLA))  # Remove old HLA column
    return df
end

# Main function to run the pipeline
function main()
    args = parse_arguments()
    folder_path = args["folder"]

    joined_df = process_and_join(folder_path)
    reshaped_df = reshape_hla_data(joined_df)

    # Drop rows where MHC is missing
    filtered_df = filter(row -> !ismissing(row.MHC), reshaped_df)

    println("Sorting by Locus...")
    sort!(filtered_df, :Locus)

    output_path = joinpath(folder_path, "processed_peptides.csv")
    println("Saving results to $output_path...")
    CSV.write(output_path, filtered_df)

    println("Processing completed successfully!")
end

main()