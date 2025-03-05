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
    joined_df = dropmissing(joined_df, :Pos)
    select!(joined_df, Not(:Pos))  # Drop the 'Pos' column

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

# Filter out non-binding peptides where EL_rank is >2
function filter_non_binding_peptides(df::DataFrame)::DataFrame
    grouped_counts = Dict()

    for hla in unique(df.MHC)
        hla_df = filter(row -> row.MHC == hla, df)
        total_count = nrow(hla_df)
        dropped_count = sum(hla_df.EL_Rank .> 2)
        percent_dropped = total_count > 0 ? round(dropped_count / total_count * 100, digits=2) : 0.0

        if dropped_count > 0
            grouped_counts[hla] = (dropped_count, percent_dropped)
        end
    end

    # Print count and percentage of dropped peptides for each HLA type
    for (hla, (count, percent)) in grouped_counts
        println("Dropped $count peptides for $hla ($percent% of total) as non-binding.")
    end

    return filter(row -> row.EL_Rank <= 2, df)
end

# Main function to run the pipeline
function main()
    args = parse_arguments()
    folder_path = args["folder"]

    println("Processing and joining files...")
    result_df = process_and_join(folder_path)

    println("Reshaping HLA-related data...")
    reshaped_df = reshape_hla_data(result_df)

    println("Filtering out non-binding peptides...")
    filtered_df = filter_non_binding_peptides(reshaped_df)

    # Sort by Locus before saving
    println("Sorting by Locus...")
    sort!(filtered_df, :Locus)

    # Save the output to filtered_peptides.csv in the input folder
    output_path = joinpath(folder_path, "filtered_peptides.csv")
    println("Saving results to $output_path...")
    CSV.write(output_path, filtered_df)

    println("Processing completed successfully!")
end

main()
