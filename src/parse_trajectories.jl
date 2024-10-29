#!/usr/bin/env julia

include("./env.jl")

using DataFrames, CSV, FilePathsBase

# Function to read Times.in and return a vector of time points
function read_times(times_file)
    times = []
    open(times_file, "r") do f
        for line in eachline(f)
            push!(times, parse(Int, strip(line)))
        end
    end
    return times
end

# Function to parse the .out file and return a DataFrame
function parse_trajectories(trajectories_file, times)
    # Initialize the DataFrame with predefined columns and correct types
    data = DataFrame(Locus = Int[], Consensus = String[], Variant = String[])
    
    for t in times
        for base in ["A", "C", "G", "T"]
            col_name = Symbol("Count_$(base)_t$(t)")
            data[!, col_name] = Int[]
        end
        total_reads_col_name = Symbol("Total_Reads_t$(t)")
        data[!, total_reads_col_name] = Int[]
    end

    open(trajectories_file, "r") do f
        while !eof(f)
            line = readline(f)
            # Skip empty lines
            if isempty(strip(line))
                continue
            end
            tokens = split(line)
            # Header for a locus
            locus = parse(Int, tokens[1])
            consensus = tokens[2]
            variant = tokens[3]
            num_timepoints = parse(Int, tokens[4])
            # Initialize a dictionary to hold counts for each time point
            row = Dict{Symbol, Any}()
            row[:Locus] = locus
            row[:Consensus] = consensus
            row[:Variant] = variant
            # Initialize counts with zeros for all time points
            for t in times
                for base in ["A", "C", "G", "T"]
                    row[Symbol("Count_$(base)_t$(t)")] = 0
                end
                row[Symbol("Total_Reads_t$(t)")] = 0
            end
            # Remaining tokens are the data points
            data_tokens = tokens[5:end]
            idx = 1
            for _ in 1:num_timepoints
                time = parse(Int, data_tokens[idx])
                idx += 1
                count_a = parse(Int, data_tokens[idx])
                idx += 1
                count_c = parse(Int, data_tokens[idx])
                idx += 1
                count_g = parse(Int, data_tokens[idx])
                idx += 1
                count_t = parse(Int, data_tokens[idx])
                idx += 1
                total_reads = parse(Int, data_tokens[idx])
                idx += 1
                # Update counts in row
                row[Symbol("Count_A_t$time")] = count_a
                row[Symbol("Count_C_t$time")] = count_c
                row[Symbol("Count_G_t$time")] = count_g
                row[Symbol("Count_T_t$time")] = count_t
                row[Symbol("Total_Reads_t$time")] = total_reads
            end
            # Now that row has all the columns, we can push it into the DataFrame
            push!(data, row)
        end
    end
    return data
end

# Main function to process command-line arguments and run the script
function main()
    # Check if the required arguments are provided
    if length(ARGS) < 1
        println("Usage: julia script_name.jl path_to_out_file [path_to_Times.in]")
        exit(1)
    end

    # Paths to the input files
    trajectories_file = ARGS[1]

    if length(ARGS) >= 2
        times_file = ARGS[2]
    else
        # Default to Times.in in the same directory as the .out file
        trajectories_path = Path(trajectories_file)
        times_file = string(parent(trajectories_path), "/Times.in")
    end

    # Check if files exist
    if !isfile(trajectories_file)
        println("Error: Trajectories file not found at path: $trajectories_file")
        exit(1)
    end

    if !isfile(times_file)
        println("Error: Times.in file not found at path: $times_file")
        exit(1)
    end

    # Read the times and trajectories
    times = read_times(times_file)
    trajectories_df = parse_trajectories(trajectories_file, times)

    # Save the DataFrame to a CSV file
    output_file = string(parent(Path(trajectories_file)), "/trajectories.csv")
    CSV.write(output_file, trajectories_df)
    println("DataFrame saved to $output_file")
end

# Run the main function
main()