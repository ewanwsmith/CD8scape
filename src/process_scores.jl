#!/usr/bin/env julia

# Include environment settings
include("./env.jl")

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
    mhcpan_path = joinpath(folder_path, "netMHCPan_output_processed.csv")
    peptides_path = joinpath(folder_path, "peptides_labels.csv")

    if !(isfile(mhcpan_path) && isfile(peptides_path))
        error("Required files not found. Ensure both netMHCPan_output_processed.csv and peptides_labels.csv are present in $folder_path.")
    end

    mhcpan_df = CSV.read(mhcpan_path, DataFrame)
    peptides_df = CSV.read(peptides_path, DataFrame)

    joined_df = leftjoin(peptides_df, mhcpan_df, on="Peptide")
    dropped_rows = joined_df[ismissing.(joined_df.Pos), :]
    joined_df = dropmissing(joined_df, :Pos)

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
    hla_score_columns = filter(col -> occursin(r"_EL-score$", col), names(df))
    hla_rank_columns = filter(col -> occursin(r"_EL_Rank$", col), names(df))
    hla_types_score = replace.(hla_score_columns, r"_EL-score$" => "")
    hla_types_rank = replace.(hla_rank_columns, r"_EL_Rank$" => "")
    hla_types = intersect(hla_types_score, hla_types_rank)

    long_df = DataFrame()

    for hla in hla_types
        score_col = Symbol("$hla" * "_EL-score")
        rank_col = Symbol("$hla" * "_EL_Rank")

        hla_df = select(df, [:Locus, :Peptide, :Peptide_label, :Pos, :ID], copycols=false)
        hla_df[:, :MHC] = repeat([hla], nrow(hla_df))
        hla_df[:, :EL_score] = df[!, score_col]
        hla_df[:, :EL_rank] = df[!, rank_col]

        append!(long_df, hla_df)
    end

    return long_df
end

# Combine rows based on Peptide_label and MHC
function combine_rows_by_label_and_mhc(df::DataFrame)::DataFrame
    sort!(df, [:MHC, :Peptide_label])
    mhc_levels = unique(df.MHC)
    combined_df = DataFrame()

    for mhc in mhc_levels
        mhc_df = filter(row -> row.MHC == mhc, df)
        c_rows = filter(row -> endswith(row.Peptide_label, "_C"), eachrow(mhc_df))
        v_rows = filter(row -> endswith(row.Peptide_label, "_V"), eachrow(mhc_df))

        c_labels = replace.(getproperty.(c_rows, :Peptide_label), r"_C$" => "")
        v_labels = replace.(getproperty.(v_rows, :Peptide_label), r"_V$" => "")
        common_labels = intersect(c_labels, v_labels)

        for label in common_labels
            c_row = first(filter(row -> startswith(row.Peptide_label, "$label"), c_rows))
            v_row = first(filter(row -> startswith(row.Peptide_label, "$label"), v_rows))

            combined_row = DataFrame(
                Locus = [c_row.Locus],
                Peptide = [c_row.Peptide],
                Peptide_label = [label],
                MHC = [mhc],
                EL_score_C = [c_row.EL_score],
                EL_score_V = [v_row.EL_score],
                EL_rank_C = [c_row.EL_rank],
                EL_rank_V = [v_row.EL_rank]
            )

            append!(combined_df, combined_row)
        end
    end

    return combined_df
end

# Calculate net scores
function calculate_net_scores(df::DataFrame)::DataFrame
    required_columns = ["EL_score_C", "EL_score_V", "EL_rank_C", "EL_rank_V"]
    for col in required_columns
        if !(col in names(df))
            error("Missing required column: $col")
        end
    end

    df[:, :EL_score_net] = df[:, Symbol("EL_score_C")] .- df[:, Symbol("EL_score_V")]
    df[:, :EL_rank_net] = df[:, Symbol("EL_rank_C")] .- df[:, Symbol("EL_rank_V")]

    return df
end

# Main function to run the pipeline
function main()
    args = parse_arguments()
    folder_path = args["folder"]

    println("Processing and joining files...")
    result_df = process_and_join(folder_path)

    println("Reshaping HLA-related data...")
    long_result_df = reshape_hla_data(result_df)

    println("Combining rows by Peptide_label and MHC...")
    combined_result_df = combine_rows_by_label_and_mhc(long_result_df)

    println("Calculating net scores...")
    final_result_df = calculate_net_scores(combined_result_df)

    # Save the output to net_scores.csv in the input folder
    output_path = joinpath(folder_path, "net_scores.csv")
    println("Saving results to $output_path...")
    CSV.write(output_path, final_result_df)

    println("Processing completed successfully!")
end

main()