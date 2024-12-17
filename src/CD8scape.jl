#!/usr/bin/env julia

# CD8scape.jl
#
# Usage:
#   ./CD8scape.jl read samfire <folder_path>
#     - runs parse_trajectories.jl and read_samfire_frames.jl
#   ./CD8scape.jl read ncbi <folder_path>
#     - runs parse_vcf.jl and read_ncbi_frames.jl
#
# Make sure you have set executable permissions:
#   chmod +x CD8scape.jl
#
# Also ensure that env.jl, parse_trajectories.jl, read_samfire_frames.jl,
# parse_vcf.jl, and read_ncbi_frames.jl are in the same directory or 
# correctly referenced.

if length(ARGS) < 3
    println("""
Usage:
  ./CD8scape.jl read samfire <folder_path>
  ./CD8scape.jl read ncbi <folder_path>
""")
    exit(1)
end

command = ARGS[1]
subcommand = ARGS[2]
folder_path = ARGS[3]

# Include the environment setup (activates local environment)
include("./env.jl")

if command == "read"
    if subcommand == "samfire"
        # Run SamFire-specific scripts
        run(`julia --project=. ./parse_trajectories.jl $folder_path`)
        run(`julia --project=. ./read_samfire_frames.jl $folder_path`)
        println("SamFire reading complete.")
    elseif subcommand == "ncbi"
        # Run NCBI-specific scripts
        run(`julia --project=. ./parse_vcf.jl $folder_path`)
        run(`julia --project=. ./read_ncbi_frames.jl $folder_path`)
        println("NCBI reading complete.")
    else
        println("Error: Unsupported subcommand '$subcommand'. Choose 'samfire' or 'ncbi'.")
        exit(1)
    end
else
    println("Error: Invalid command '$command'. Use 'read'.")
    exit(1)
end

println("All scripts have run successfully.")