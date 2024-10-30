include("./env.jl")

using DataFrames
using CSV

# Define a function to perform the operation with the folder path as an argument
function process_data(folder_path::String)
    # Construct file paths for frames and trajectories CSV files
    frames_path = joinpath(folder_path, "frames.csv")
    trajectories_path = joinpath(folder_path, "trajectories.csv")
    
    # Read the frames and trajectories data
    frames = CSV.read(frames_path, DataFrame)
    trajectories = CSV.read(trajectories_path, DataFrame)
    
    # Create Start and End columns by extracting the first Start and last End values
    frames[!, :Start] = [parse(Int, split(split(region, ";")[1], ",")[1]) for region in frames.Region]
    frames[!, :End] = [parse(Int, split(split(region, ";")[end], ",")[2]) for region in frames.Region]
    
    # Perform a cross join and then filter based on the condition
    result = crossjoin(trajectories, frames)
    result = filter(row -> row.Start <= row.Locus <= row.End, result)

    return result
end

joined = process_data("/Users/e.smith.5/Documents/PhD/CD8scape/data/RSV_example")

# generate a relative locus for use in the sequence
joined[!, :Relative_Locus] = joined.Locus .- joined.Start

