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
#     Sets up the environment by running env.jl. This command should be run
#     once before using "read" or "run".
#
#   read:
#     Attempts to parse data (trajectories or VCF) and then read frames 
#     (SamFire or NCBI) in order. It tries parse_trajectories first, and if that 
#     fails, parse_vcf. Then tries read_samfire_frames, and if that fails, read_ncbi_frames.
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
  prep    Set up the environment by running env.jl.
  read    Attempt to parse trajectories (or VCF) and read frames (SamFire or NCBI).
  run     Run the peptide-generation and NetMHCpan pipeline on parsed data.

OPTIONS:
  --help, -h
      Print this help message and exit.

DESCRIPTION:
  The 'prep' command initializes the environment by executing env.jl. This 
  should be done before using the 'read' or 'run' commands to ensure that 
  all necessary dependencies and settings are configured.

  The 'read' command attempts to parse and read data from the specified folder:
    1. parse_trajectories.jl is tried first.
       If it fails, parse_vcf.jl is tried.
    2. read_samfire_frames.jl is then tried.
       If it fails, read_ncbi_frames.jl is tried.

  If both attempts fail at either step, the process is aborted with an error message.

  The 'run' command assumes that data parsing/reading is complete and attempts to:
    - Generate peptides with generate_peptides.jl
    - Run NetMHCpan predictions with run_netMHCpan.jl
    - Process the resulting NetMHCpan output with process_output.pl
    - Process the processed scores with process_scores.jl
    - Plot those scores with plot_scores.jl

EXAMPLES:
  # Setting up the environment
  ./CD8scape.jl prep

  # Parsing data from a folder "my_data"
  ./CD8scape.jl read my_data

  # Running the pipeline on a folder "my_data"
  ./CD8scape.jl run my_data

Ensure all required scripts (env.jl, parse_trajectories.jl, parse_vcf.jl, 
read_samfire_frames.jl, read_ncbi_frames.jl, generate_peptides.jl, run_netMHCpan.jl, 
process_output.pl, process_scores.jl, and plot_scores.jl) are in the same directory or correctly referenced.
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
    # We expect only 1 argument: prep
    if length(ARGS) > 1
        println("Error: 'prep' command does not take any additional arguments.")
        println("Usage: ./CD8scape.jl prep")
        exit(1)
    end

    # Include the environment setup (activates local environment).
    try
        include("./env.jl")
        println("Environment setup completed successfully by running env.jl.")
    catch e
        println("Error running env.jl: $e")
        exit(1)
    end

elseif command == "read"
    # We expect 2 arguments: read, folder_path
    if length(ARGS) < 2
        println("Error: Missing folder_path for read command.")
        println("Usage: ./CD8scape.jl read <folder_path>")
        exit(1)
    end

    folder_path = ARGS[2]

    # 1. Parse trajectories or VCF
    parse_success = safe_run(`julia --project=. ./parse_trajectories.jl $folder_path`)
    if !parse_success
        println("parse_trajectories.jl failed, attempting parse_vcf.jl ...")
        parse_success = safe_run(`julia --project=. ./parse_vcf.jl $folder_path`)
        if !parse_success
            println("Both parse_trajectories.jl and parse_vcf.jl failed.")
            println("Error: Unable to parse trajectory or VCF data.")
            exit(1)
        end
    end

    # 2. Read SamFire frames or NCBI frames
    samfire_success = safe_run(`julia --project=. ./read_samfire_frames.jl $folder_path`)
    if !samfire_success
        println("SamFire frames could not be read (e.g., Reading_Frames.dat not found).")
        println("Attempting to read NCBI frames instead...")
        ncbi_success = safe_run(`julia --project=. ./read_ncbi_frames.jl $folder_path`)
        if !ncbi_success
            println("Both read_samfire_frames.jl and read_ncbi_frames.jl failed.")
            println("Error: Unable to read frames from SamFire or NCBI data.")
            exit(1)
        else
            # If NCBI frames succeeded
            println("Frames DataFrame saved to $(joinpath(folder_path, "frames.csv"))")
        end
    else
        # If SamFire frames succeeded
        println("Frames DataFrame saved to $(joinpath(folder_path, "frames.csv"))")
    end

    println("Reading stage finished successfully.")

elseif command == "run"
    # We expect 2 arguments: run, folder_path
    if length(ARGS) < 2
        println("Error: Missing folder_path for run command.")
        println("Usage: ./CD8scape.jl run <folder_path>")
        exit(1)
    end

    folder_path = ARGS[2]
    
    # Define paths for NetMHCpan output and processed output
    netmhcpan_output = joinpath(folder_path, "netmhcpan_output.tsv")        # Adjust if different
    processed_output = joinpath(folder_path, "processed_output.csv")        # Adjust if different

    # Verify NetMHCpan output exists
    if !isfile(netmhcpan_output)
        println("Error: NetMHCpan output file '$netmhcpan_output' does not exist.")
        exit(1)
    end

    # Run the pipeline steps, checking success at each step

    # 1. Generate Peptides
    if !safe_run(`julia --project=. ./generate_peptides.jl $folder_path`)
        println("Error running generate_peptides.jl")
        exit(1)
    end

    # 2. Run NetMHCpan
    if !safe_run(`julia --project=. ./run_netMHCpan.jl --folder $folder_path`)
        println("Error running run_netMHCpan.jl")
        exit(1)
    end

    # 3. Process NetMHCpan Output with Perl Script by capturing output in Julia
    try
        # Capture the output of the Perl script
        perl_output = read(`perl ./process_output.pl $netmhcpan_output`, String)
        
        # Write the captured output to processed_output.csv
        open(processed_output, "w") do f
            write(f, perl_output)
        end
        
        println("Successfully ran process_output.pl. Output saved to $processed_output")
    catch e
        println("Error running process_output.pl: $e")
        exit(1)
    end

    # 4. Process Scores
    if !safe_run(`julia --project=. ./process_scores.jl --folder $folder_path`)
        println("Error running process_scores.jl")
        exit(1)
    end

    # 5. Plot Scores
    if !safe_run(`julia --project=. ./plot_scores.jl $folder_path`)
        println("Error running plot_scores.jl")
        exit(1)
    end

    println("Pipeline run stage finished successfully.")

else
    # If the command is not recognized, print help and exit.
    println("Error: Invalid command '$command'. Expected 'prep', 'read', or 'run'.")
    print_help()
    exit(1)
end