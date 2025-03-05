#!/usr/bin/env julia

###############################################################################
# CD8scape.jl
#
# USAGE:
#   ./CD8scape.jl prep
#   ./CD8scape.jl read <folder_path>
#   ./CD8scape.jl run  <folder_path>
#
# DESCRIPTION:
#   This script provides three main functionalities: "prep", "read", and "run".
#
#   prep:
#     Sets up the environment by running src/env.jl. This command should be run
#     once before using "read" or "run".
#
#   read:
#     Attempts to parse data (trajectories or VCF) and then read frames 
#     (SamFire or NCBI) in order. It tries src/parse_trajectories.jl first, and if that 
#     fails, src/parse_vcf.jl. Then tries src/read_samfire_frames.jl, and if that fails, src/read_ncbi_frames.jl.
#
#   run:
#     Executes the peptide-generation + NetMHCpan pipeline with <folder_path>.
#
# OPTIONS:
#   --help, -h : Show this help message and exit.
#
###############################################################################

function print_help()
    println("""
CD8scape.jl - A tool for running netMHCpan and managing related data.

USAGE:
  ./CD8scape.jl prep
  ./CD8scape.jl read <folder_path>
  ./CD8scape.jl run  <folder_path>

COMMANDS:
  prep    Set up the environment by running src/env.jl.
  read    Attempt to parse trajectories (or VCF) and read frames (SamFire or NCBI).
  run     Run the peptide-generation and NetMHCpan pipeline on parsed data.

OPTIONS:
  --help, -h
      Print this help message and exit.

DESCRIPTION:
  The 'prep' command initializes the environment by executing src/env.jl. This 
  should be done before using the 'read' or 'run' commands to ensure that 
  all necessary dependencies and settings are configured.

  The 'read' command attempts to parse and read data from the specified folder:
    1. src/parse_trajectories.jl is tried first.
       If it fails, src/parse_vcf.jl is tried.
    2. src/read_samfire_frames.jl is then tried.
       If it fails, src/read_ncbi_frames.jl is tried.

  If both attempts fail at either step, the process is aborted with an error message.

  The 'run' command assumes that data parsing/reading is complete and attempts to:
    - Generate peptides with src/generate_peptides.jl
    - Run NetMHCpan predictions with src/run_netMHCpan.jl
    - Process the resulting NetMHCpan output with src/process_output.pl
    - Process the processed scores with src/process_scores.jl
    - Plot those scores with src/plot_scores.jl

EXAMPLES:
  # Setting up the environment
  ./CD8scape.jl prep

  # Parsing data from a folder "my_data"
  ./CD8scape.jl read my_data

  # Running the pipeline on a folder "my_data"
  ./CD8scape.jl run my_data

Ensure all required scripts (src/env.jl, src/parse_trajectories.jl, src/parse_vcf.jl, 
src/read_samfire_frames.jl, src/read_ncbi_frames.jl, src/generate_peptides.jl, src/run_netMHCpan.jl, 
src/process_output.pl, src/process_scores.jl, and src/plot_scores.jl) are in the `src/` directory.
""")
end

# Helper function to run a command safely and silently handle errors.
# We redirect stderr to /dev/null to prevent error logs from child scripts.
function safe_run(cmd::Cmd)
    devnull = open("/dev/null", "w")
    try
        run(pipeline(cmd, stderr=devnull))
        return true
    catch
        return false
    finally
        close(devnull)
    end
end

# If no arguments provided or help requested, print help.
if length(ARGS) == 0 || ARGS[1] == "--help" || ARGS[1] == "-h"
    print_help()
    exit(0)
end

command = ARGS[1]

# Process commands based on user input
if command == "prep"
    if length(ARGS) > 1
        println("Error: 'prep' command does not take any additional arguments.")
        exit(1)
    end

    try
        include("src/env.jl")
        println("Environment setup completed successfully by running src/env.jl.")
    catch e
        println("Error running src/env.jl: $e")
        exit(1)
    end

elseif command == "read"
    if length(ARGS) < 2
        println("Error: Missing folder_path for read command.")
        exit(1)
    end

    folder_path = ARGS[2]

    parse_success = safe_run(`julia --project=. src/parse_trajectories.jl $folder_path`)
    if !parse_success
        println("parse_trajectories.jl failed, attempting parse_vcf.jl ...")
        parse_success = safe_run(`julia --project=. src/parse_vcf.jl $folder_path`)
        if !parse_success
            println("Both src/parse_trajectories.jl and src/parse_vcf.jl failed.")
            exit(1)
        end
    end

    samfire_success = safe_run(`julia --project=. src/read_samfire_frames.jl $folder_path`)
    if !samfire_success
        println("Attempting src/read_ncbi_frames.jl ...")
        ncbi_success = safe_run(`julia --project=. src/read_ncbi_frames.jl $folder_path`)
        if !ncbi_success
            println("Both src/read_samfire_frames.jl and src/read_ncbi_frames.jl failed.")
            exit(1)
        end
    end

    println("Reading stage finished successfully.")

elseif command == "run"
    if length(ARGS) < 2
        println("Error: Missing folder_path for run command.")
        exit(1)
    end

    folder_path = ARGS[2]
    
    netmhcpan_output = joinpath(folder_path, "netmhcpan_output.tsv")
    processed_output = joinpath(folder_path, "processed_output.csv")

    if !isfile(netmhcpan_output)
        println("Error: NetMHCpan output file '$netmhcpan_output' does not exist.")
        exit(1)
    end

    if !safe_run(`julia --project=. src/generate_peptides.jl $folder_path`)
        println("Error running src/generate_peptides.jl")
        exit(1)
    end

    if !safe_run(`julia --project=. src/run_netMHCpan.jl --folder $folder_path`)
        println("Error running src/run_netMHCpan.jl")
        exit(1)
    end

    try
        perl_output = read(`perl src/process_output.pl $netmhcpan_output`, String)
        open(processed_output, "w") do f
            write(f, perl_output)
        end
        println("Successfully ran src/process_output.pl.")
    catch e
        println("Error running src/process_output.pl: $e")
        exit(1)
    end

    if !safe_run(`julia --project=. src/process_scores.jl --folder $folder_path`)
        println("Error running src/process_scores.jl")
        exit(1)
    end

    println("Run stage finished successfully.")

else
    println("Error: Invalid command '$command'.")
    print_help()
    exit(1)
end