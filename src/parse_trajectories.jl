#!/usr/bin/env julia

using DataFrames
using CSV
using FilePathsBase

"""
    parse_trajectories(trajectories_file)

Parse the `.out` file, discover all unique time points, build columns for 
each time point, and return a DataFrame. No separate `Times.in` file is needed.
"""
function parse_trajectories(trajectories_file)
    # Storage for raw data (one Dict per "locus" line in .out)
    data_storage = []

    # A set to collect all unique time points across the entire file
    unique_times = Set{Int}()

    open(trajectories_file, "r") do f
        while !eof(f)
            line = readline(f)
            if isempty(strip(line))
                continue  # skip empty lines
            end

            tokens = split(line)
            # tokens structure:
            # 1) Locus
            # 2) Consensus
            # 3) Variant
            # 4) num_timepoints
            # 5..end) (time + A + C + G + T + total_reads) repeated

            locus           = parse(Int, tokens[1])
            consensus       = tokens[2]
            variant         = tokens[3]
            num_timepoints  = parse(Int, tokens[4])
            data_tokens     = tokens[5:end]

            row_dict = Dict{Symbol, Any}(
                :Locus      => locus,
                :Consensus  => consensus,
                :Variant    => variant
            )

            idx = 1
            for _ in 1:num_timepoints
                time         = parse(Int, data_tokens[idx]);   idx += 1
                count_a      = parse(Int, data_tokens[idx]);   idx += 1
                count_c      = parse(Int, data_tokens[idx]);   idx += 1
                count_g      = parse(Int, data_tokens[idx]);   idx += 1
                count_t      = parse(Int, data_tokens[idx]);   idx += 1
                total_reads  = parse(Int, data_tokens[idx]);   idx += 1

                # Record this time in our set of unique time points
                push!(unique_times, time)

                # Store counts in the row_dict for this time
                row_dict[Symbol("Count_A_t$time")]     = count_a
                row_dict[Symbol("Count_C_t$time")]     = count_c
                row_dict[Symbol("Count_G_t$time")]     = count_g
                row_dict[Symbol("Count_T_t$time")]     = count_t
                row_dict[Symbol("Total_Reads_t$time")] = total_reads
            end

            push!(data_storage, row_dict)
        end
    end

    # Sort unique time points so columns appear in ascending time order
    time_points = sort(collect(unique_times))

    # Build the list of columns for the DataFrame
    columns = [:Locus, :Consensus, :Variant]
    for t in time_points
        for base in ["A", "C", "G", "T"]
            push!(columns, Symbol("Count_$(base)_t$(t)"))
        end
        push!(columns, Symbol("Total_Reads_t$(t)"))
    end

    # Initialize the DataFrame with the desired columns
    df = DataFrame()
    for col in columns
        df[!, col] = Any[]
    end

    # Populate the DataFrame row by row
    for row_dict in data_storage
        # Construct a row vector, ensuring every column is filled
        # Default 0 for missing timepoint columns
        row_values = Vector{Any}(undef, length(columns))
        for (i, col) in pairs(columns)
            if haskey(row_dict, col)
                row_values[i] = row_dict[col]
            elseif col in (:Locus, :Consensus, :Variant)
                # Must exist in the row
                row_values[i] = row_dict[col]
            else
                # For time-based columns that don't exist, default to 0
                row_values[i] = 0
            end
        end
        push!(df, row_values)
    end

    return df
end

function main()
    # Check that a folder path was provided
    if length(ARGS) < 1
        println("Usage: julia script.jl <folder_path>")
        return
    end

    folder_path = ARGS[1]

    # Attempt to find single_locus_trajectories.out
    trajectories_file = joinpath(folder_path, "single_locus_trajectories.out")

    if !isfile(trajectories_file)
        # If that file doesn't exist, look for any .out file in the folder
        out_files = filter(f -> endswith(f, ".out"), readdir(folder_path))
        if isempty(out_files)
            error("No 'single_locus_trajectories.out' and no other .out files found in $folder_path")
        else
            trajectories_file = joinpath(folder_path, out_files[1])
            println("No single_locus_trajectories.out found. Using `$(out_files[1])` instead.")
        end
    end

    # Build output file path
    output_file = joinpath(folder_path, "variants.csv")

    # Parse the trajectories (no need for Times.in anymore)
    trajectories_df = parse_trajectories(trajectories_file)

    # Keep only Locus, Consensus, and Variant columns
    final_df = select(trajectories_df, :Locus, :Consensus, :Variant)

    # Save the DataFrame to variants.csv
    CSV.write(output_file, final_df)
    println("DataFrame saved to $output_file")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end