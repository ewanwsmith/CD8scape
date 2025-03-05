#!/usr/bin/env julia

###############################################################################
# CD8scape.jl
#
# USAGE:
#   ./CD8scape.jl prep
#   ./CD8scape.jl read <folder_path>
#   ./CD8scape.jl run  <folder_path>
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
""")
end

# Helper function to run commands safely and handle errors
function safe_run(cmd::Cmd)
    try
        run(cmd)
        return true
    catch
        return false
    end
end

# Print help if no arguments are provided
if length(ARGS) == 0 || ARGS[1] in ["--help", "-h"]
    print_help()
    exit(0)
end

command = ARGS[1]

# Process "prep" command
if command == "prep"
    if length(ARGS) > 1
        println("Error: 'prep' command does not take any additional arguments.")
        exit(1)
    end

    try
        include("src/env.jl")
        println("Environment setup completed successfully.")
    catch e
        println("Error running src/env.jl: $e")
        exit(1)
    end

# Process "read" command
elseif command == "read"
    if length(ARGS) < 2
        println("Error: Missing folder_path for read command.")
        exit(1)
    end

    folder_path = ARGS[2]

    parse_success = safe_run(`julia --project=. src/parse_trajectories.jl $folder_path`)
    if !parse_success
        parse_success = safe_run(`julia --project=. src/parse_vcf.jl $folder_path`)
        if !parse_success
            println("Error: Both trajectory parsing methods failed.")
            exit(1)
        end
    end

    samfire_success = safe_run(`julia --project=. src/read_samfire_frames.jl $folder_path`)
    if !samfire_success
        ncbi_success = safe_run(`julia --project=. src/read_ncbi_frames.jl $folder_path`)
        if !ncbi_success
            println("Error: Both frame-reading methods failed.")
            exit(1)
        end
    end

    println("Reading stage finished successfully.")

# Process "run" command
elseif command == "run"
    if length(ARGS) < 2
        println("Error: Missing folder_path for run command.")
        exit(1)
    end

    folder_path = ARGS[2]
    netmhcpan_output = joinpath(folder_path, "netmhcpan_output.tsv")
    processed_output = joinpath(folder_path, "processed_output.csv")

    # Step 1: Generate Peptides
    if !safe_run(`julia --project=. src/generate_peptides.jl $folder_path`)
        println("Error running src/generate_peptides.jl")
        exit(1)
    end

    # Step 2: Run NetMHCpan
    if !safe_run(`julia --project=. src/run_netMHCpan.jl --folder $folder_path`)
        println("Error: NetMHCpan did not run successfully.")
        exit(1)
    end

    # Step 3: Verify NetMHCpan Output
    if !isfile(netmhcpan_output)
        println("Error: Expected NetMHCpan output file '$netmhcpan_output' does not exist.")
        exit(1)
    end

    # Step 4: Process Output with Perl Script
    try
        perl_output = read(`perl src/process_output.pl $netmhcpan_output`, String)
        open(processed_output, "w") do f
            write(f, perl_output)
        end
    catch e
        println("Error running src/process_output.pl: $e")
        exit(1)
    end

    # Step 5: Process Scores
    if !safe_run(`julia --project=. src/process_scores.jl --folder $folder_path`)
        println("Error running src/process_scores.jl")
        exit(1)
    end

    # Step 6: Process Best Ranks
    if !safe_run(`julia --project=. src/process_best_ranks.jl $folder_path`)
        println("Error running src/process_best_ranks.jl")
        exit(1)
    end

    println("Run stage finished successfully.")

else
    println("Error: Invalid command '$command'.")
    print_help()
    exit(1)
end