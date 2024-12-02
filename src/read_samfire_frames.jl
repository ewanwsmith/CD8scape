# Include environment settings
include("./env.jl")

using DataFrames
using CSV
using FilePathsBase

# Removed ArgParse as it's no longer needed
# using ArgParse

# Function to read the .dat file and create the DataFrame with a combined "Region" column
function read_dat_file(filepath::String)::DataFrame
    # Read all lines from the file
    lines = readlines(filepath)
    
    # Initialize vectors to store the data
    regions = String[]
    consensus = String[]
    
    # Total number of lines
    n = length(lines)
    
    # Iterate over the lines in steps of 3 (assuming 3 lines per entry)
    i = 1
    while i <= n
        # Parse the Start and End from the first line
        start_end = split(strip(lines[i]))
        if length(start_end) < 2
            @warn "Line $i does not contain two elements. Skipping."
            i += 1
            continue
        end
        start_val = tryparse(Int, start_end[1])
        end_val = tryparse(Int, start_end[2])
        
        if isnothing(start_val) || isnothing(end_val)
            @warn "Invalid Start or End values at line $i. Skipping."
            i += 1
            continue
        end
        
        # Combine Start and End into "Region" as "Start,End"
        region = "$(start_val),$(end_val)"
        push!(regions, region)
        
        # Read the Consensus_sequence from the next line
        if i + 1 > n
            @warn "Missing Consensus_sequence after line $i. Skipping."
            break
        end
        consensus_seq = strip(lines[i + 1])
        push!(consensus, consensus_seq)
        
        # Skip the third line (protein sequence)
        i += 3
    end
    
    # Create the DataFrame with "Region" and "Consensus_sequence"
    frames_df = DataFrame(
        Region = regions,
        Consensus_sequence = consensus
    )
    
    return frames_df
end

# Function to save DataFrame as CSV in the same directory as the input file
function save_dataframe_as_csv(df::DataFrame, input_filepath::String, output_filename::String="frames.csv")
    # Extract the directory from the input file path
    directory = dirname(input_filepath)
    
    # Construct the output file path
    output_filepath = joinpath(directory, output_filename)
    
    # Write the DataFrame to CSV
    CSV.write(output_filepath, df)
    
    println("DataFrame successfully saved to $output_filepath")
end

# Removed command-line argument parsing function
# function parse_command_line_args()
#     ...
# end

# Main function to orchestrate the processing
function main()
    # Option 1: Use a predefined file path
    # Uncomment the following line and set your file path
    # frames_filepath = "/path/to/your/frames.dat"

    # Option 2: Prompt the user for the file path
    print("Enter the path to the frames .dat file: ")
    frames_filepath = readline()

    # Check if the provided file exists
    if !isfile(frames_filepath)
        @error "The file '$frames_filepath' does not exist."
        return
    end
    
    # Read the .dat file and create the DataFrame
    frames_df = read_dat_file(frames_filepath)
    
    # Display the DataFrame (optional)
    println("Generated DataFrame:")
    println(frames_df)
    
    # Save the DataFrame as frames.csv in the same directory as frames_filepath
    save_dataframe_as_csv(frames_df, frames_filepath)
end

# Execute the main function
main()