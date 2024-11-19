include("./env.jl")

using DataFrames
using CSV
using FilePathsBase

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
function parse_trajectories(trajectories_file, times_file)
    # Read the times
    times = read_times(times_file)

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

    # Return the parsed DataFrame
    return data
end

# Function to determine the output file path
function get_output_file_path(trajectories_file)
    trajectories_path = Path(trajectories_file)
    output_file = string(parent(trajectories_path), "/trajectories.csv")
    return output_file
end

# Replace with your actual file paths
trajectories_file = "/Users/e.smith.5/Documents/PhD/CD8scape/data/RSV_example/single_locus_trajectories10.out"
times_file = "/Users/e.smith.5/Documents/PhD/CD8scape/data/RSV_example/Times.in"

# Parse the trajectories
trajectories_df = parse_trajectories(trajectories_file, times_file)

# Optionally, save the DataFrame to a CSV file
output_file = get_output_file_path(trajectories_file)
CSV.write(output_file, trajectories_df)
println("DataFrame saved to $output_file")