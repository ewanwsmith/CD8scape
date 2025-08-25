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
  ./CD8scape.jl run_supertype  <folder_path>
  ./CD8scape.jl context <folder_path>

COMMANDS:
  prep    Set up the environment by running src/env.jl.
  read    Attempt to parse variants input (samfire trajectories or .vcf) and read frames (SamFire or NCBI).
  run     Run the peptide-generation and NetMHCpan pipeline on parsed data.
  run_supertype Run the peptide-generation and NetMHCpan pipeline on parsed data for a representative supertpe HLA panel.
  context Run the context-sensitive peptide processing pipeline.

OPTIONS:
  --help, -h
      Print this help message and exit.
""")
end

# Helper function to run commands safely
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

    # Try NCBI frame reading first
    ncbi_success = safe_run(`julia --project=. src/read_ncbi_frames.jl $folder_path`)
    frames_csv_path = joinpath(folder_path, "frames.csv")

    if !ncbi_success || !isfile(frames_csv_path)
        println("read_ncbi_frames.jl failed or frames.csv not found. Trying read_samfire_frames.jl instead.")
        samfire_success = safe_run(`julia --project=. src/read_samfire_frames.jl $folder_path`)

        if !samfire_success || !isfile(frames_csv_path)
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

    # Generate Peptides
    if !safe_run(`julia --project=. src/generate_peptides.jl $folder_path`)
        println("Error running src/generate_peptides.jl")
        exit(1)
    end

    # Clean Peptides
    if !safe_run(`julia --project=. src/clean_peptides.jl $folder_path`)
        println("Error running src/clean_peptides.jl")
        exit(1)
    end

    # Run NetMHCpan
    if !safe_run(`julia --project=. src/run_netMHCpan.jl --folder $folder_path`)
        println("Error: NetMHCpan did not run successfully.")
        exit(1)
    end

    # Verify NetMHCpan Output
    if !isfile(netmhcpan_output)
        println("Error: Expected NetMHCpan output file '$netmhcpan_output' does not exist.")
        exit(1)
    end

    # Process Output with Perl Script
    try
        perl_output = read(`perl src/process_output.pl $netmhcpan_output`, String)
        open(processed_output, "w") do f
            write(f, perl_output)
        end
    catch e
        println("Error running src/process_output.pl: $e")
        exit(1)
    end

    # Process Scores
    if !safe_run(`julia --project=. src/process_scores.jl --folder $folder_path`)
        println("Error running src/process_scores.jl")
        exit(1)
    end

    # Process Best Ranks
    if !safe_run(`julia --project=. src/process_best_ranks.jl $folder_path`)
        println("Error running src/process_best_ranks.jl")
        exit(1)
    end

    println("Run stage finished successfully.")

# Process "run_supertype" command
elseif command == "run_supertype"
    if length(ARGS) < 2
        println("Error: Missing folder_path for run_supertype command.")
        exit(1)
    end

    folder_path = ARGS[2]
    netmhcpan_output = joinpath(folder_path, "netmhcpan_output.tsv")
    processed_output = joinpath(folder_path, "processed_output.csv")

    # Generate Peptides
    if !safe_run(`julia --project=. src/generate_peptides.jl $folder_path`)
        println("Error running src/generate_peptides.jl")
        exit(1)
    end

    # Clean Peptides
    if !safe_run(`julia --project=. src/clean_peptides.jl $folder_path`)
        println("Error running src/clean_peptides.jl")
        exit(1)
    end

    # Run NetMHCpan
    if !safe_run(`julia --project=. src/run_netMHCpan_global.jl --folder $folder_path`)
        println("Error: NetMHCpan did not run successfully.")
        exit(1)
    end

    # Verify NetMHCpan Output
    if !isfile(netmhcpan_output)
        println("Error: Expected NetMHCpan output file '$netmhcpan_output' does not exist.")
        exit(1)
    end

    # Process Output with Perl Script
    try
        perl_output = read(`perl src/process_output.pl $netmhcpan_output`, String)
        open(processed_output, "w") do f
            write(f, perl_output)
        end
    catch e
        println("Error running src/process_output.pl: $e")
        exit(1)
    end

    # Process Scores
    if !safe_run(`julia --project=. src/process_scores.jl --folder $folder_path`)
        println("Error running src/process_scores.jl")
        exit(1)
    end

    # Process Best Ranks
    if !safe_run(`julia --project=. src/process_best_ranks.jl $folder_path`)
        println("Error running src/process_best_ranks.jl")
        exit(1)
    end

    println("run_supertype stage finished successfully.")

# Process "context" command
elseif command == "context"

    if length(ARGS) < 2
        println("Error: Missing folder_path for context command.")
        exit(1)
    end

    folder_path = ARGS[2]
    remaining_args = ARGS[3:end]

    # Check for --force
    force = any(x -> x == "--force", remaining_args)

    # Output file to check
    context_output = joinpath(folder_path, "context_scores.csv")

    if force
        println("--force specified: Running full context pipeline (this may take a while!)")
        local cmd = `julia --project=. src/context_run/context_run.jl $folder_path`
        for arg in remaining_args
            cmd = `$cmd $arg`
        end
        if !safe_run(cmd)
            println("Error running src/context_run/context_run.jl")
            exit(1)
        end
        # Always run distribution step after full context run
        println("Running distribution step...")
        if !safe_run(`julia --project=. src/context_run/context_distribution.jl $folder_path`)
            println("Error running context_distribution.jl")
            exit(1)
        end
        println("Context stage finished successfully.")
    elseif isfile(context_output)
        println("Context scores found: Carrying on from context_scores.csv and running remaining steps.")
        # Run process_best_ranks_context.jl (panel or supertype)
        mode = "panel"
        if any(x -> x == "--supertype", remaining_args)
            mode = "supertype"
        end
        best_ranks_cmd = mode == "supertype" ?
            `julia --project=. src/context_run/process_best_ranks_context.jl $folder_path --supertype` :
            `julia --project=. src/context_run/process_best_ranks_context.jl $folder_path`
        if !safe_run(best_ranks_cmd)
            println("Error running process_best_ranks_context.jl")
            exit(1)
        end
        # Run distribution step
        println("Running distribution step...")
        if !safe_run(`julia --project=. src/context_run/context_distribution.jl $folder_path`)
            println("Error running context_distribution.jl")
            exit(1)
        end
        println("Context distribution step finished successfully.")
    else
        println("Context scores not found: Running full context pipeline (this may take a while !)")
        local cmd = `julia --project=. src/context_run/context_run.jl $folder_path`
        for arg in remaining_args
            cmd = `$cmd $arg`
        end
        if !safe_run(cmd)
            println("Error running src/context_run/context_run.jl")
            exit(1)
        end
        # Always run distribution step after full context run
        println("Running distribution step...")
        if !safe_run(`julia --project=. src/context_run/context_distribution.jl $folder_path`)
            println("Error running context_distribution.jl")
            exit(1)
        end
        println("Context stage finished successfully.")
    end

else
    println("Error: Invalid command '$command'.")
    print_help()
    exit(1)
end