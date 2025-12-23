#!/usr/bin/env julia
"""
process_scores.jl

Joins NetMHCpan output with peptide labels for downstream analysis.

Usage:
    julia process_scores.jl --folder <folder_path> [--suffix <name>] [--latest]

Arguments:
    --folder   Path to the folder containing NetMHCpan and label files.
"""

using DataFrames, CSV
include("path_utils.jl")

# Manual argument parsing function
function parse_arguments()
    args = Dict()
    for (i, arg) in enumerate(ARGS)
        if arg in ["--folder", "-f"]
            args["folder"] = ARGS[i + 1]
        elseif arg == "--suffix"
            if i + 1 <= length(ARGS) && !startswith(ARGS[i+1], "--")
                args["suffix"] = ARGS[i + 1]
            else
                args["suffix"] = ""
            end
        elseif arg == "--latest"
            args["latest"] = true
        end
    end
    if !haskey(args, "folder")
        error("Usage: ./process_scores.jl --folder /path/to/data [--suffix <name>] [--latest]")
    end
    return args
end

# Process and join files
function process_and_join(folder_path::String; suffix::AbstractString="", latest::Bool=false)::DataFrame
    mhcpan_path    = resolve_read(joinpath(folder_path, "processed_output.csv"); suffix=suffix, latest=latest)
    peptides_path  = resolve_read(joinpath(folder_path, "peptides_labels.csv"); suffix=suffix, latest=latest)

    if !(isfile(mhcpan_path) && isfile(peptides_path))
        error("Required files not found. Ensure both processed_output.csv and peptides_labels.csv are present in $folder_path.")
    end

    mhcpan_df = CSV.read(mhcpan_path, DataFrame)
    peptides_df = CSV.read(peptides_path, DataFrame)

    # Perform left join: attach Peptide_label and Locus to mhcpan_df
    joined_df = leftjoin(mhcpan_df, peptides_df, on="Peptide")

    dropped_rows = joined_df[ismissing.(joined_df.Pos), :]

    if "Locus" in names(joined_df)
        if !isempty(dropped_rows)
            grouped_dropped = combine(groupby(dropped_rows, :Locus), nrow => :Dropped_Count)
            total_per_locus = combine(groupby(joined_df, :Locus), nrow => :Total_Count)

            stats = leftjoin(grouped_dropped, total_per_locus, on=:Locus)
            stats[:, :Percentage] = round.((stats.Dropped_Count ./ stats.Total_Count) .* 100, digits=2)

            for row in eachrow(stats)
                println("$(row.Dropped_Count) peptides dropped from locus $(row.Locus) due to stop codon(s) in the peptide ($(row.Percentage)% of total for locus $(row.Locus))")
            end
        end
    else
        println("Column 'Locus' not found after join. This usually means no matching peptides between NetMHCpan output and peptide labels. Skipping locus-based stats.")
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
    suffix  = get(args, "suffix", "")
    latest  = get(args, "latest", true)
    joined_df = process_and_join(folder_path; suffix=suffix, latest=latest)
    reshaped_df = reshape_hla_data(joined_df)

    # Drop rows where MHC is missing
    filtered_df = filter(row -> !ismissing(row.MHC), reshaped_df)

    println("Sorting by Locus...")
    sort!(filtered_df, :Locus)

    output_path = resolve_write(joinpath(folder_path, "processed_peptides.csv"); suffix=suffix)
    println("Saving results to $output_path...")
    CSV.write(output_path, filtered_df)

    println("Processing completed successfully!")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end