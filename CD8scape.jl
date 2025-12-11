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
  

COMMANDS:
  prep    Set up the environment by running src/env.jl.
  read    Attempt to parse variants input (samfire trajectories or .vcf) and read frames (SamFire or NCBI).
  simulate Generate frames and exhaustive amino-acid variant set per frame.
  run     Run the peptide-generation and NetMHCpan pipeline on parsed data.
  run_supertype Run the peptide-generation and NetMHCpan pipeline on parsed data for a representative supertpe HLA panel.
  

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
    # Setup environment
    if !safe_run(`julia --project=. src/env.jl`)
        println("Error running src/env.jl")
        exit(1)
    end
    println("Environment setup finished successfully.")

# Process "read" command
elseif command == "read"
    if length(ARGS) < 2
        println("Error: Missing folder_path for read command.")
        exit(1)
    end
    folder_path = ARGS[2]
    frames_csv_path = joinpath(folder_path, "frames.csv")
    variants_csv_path = joinpath(folder_path, "variants.csv")

    # Read frames (NCBI first, then Samfire fallback)
    ncbi_success = safe_run(`julia --project=. src/read_ncbi_frames.jl $folder_path`)
    if !ncbi_success || !isfile(frames_csv_path)
        println("read_ncbi_frames.jl failed or frames.csv not found. Trying read_samfire_frames.jl instead.")
        samfire_success = safe_run(`julia --project=. src/read_samfire_frames.jl $folder_path`)
        if !samfire_success || !isfile(frames_csv_path)
            println("Error: Both frame-reading methods failed.")
            exit(1)
        end
    end

    # Parse variants: prefer VCF if present, otherwise trajectories
    vcf_files = filter(f -> endswith(f, ".vcf") || endswith(f, ".vcf.gz"), readdir(folder_path; join=true))
    parse_ok = false
    if !isempty(vcf_files)
        parse_ok = safe_run(`julia --project=. src/parse_vcf.jl $folder_path`)
        if !parse_ok
            println("parse_vcf.jl failed; trying parse_trajectories.jl.")
        end
    end
    if !parse_ok
        parse_ok = safe_run(`julia --project=. src/parse_trajectories.jl $folder_path`)
    end
    if !parse_ok || !isfile(variants_csv_path)
        println("Error: Failed to parse variants from VCF or trajectories.")
        exit(1)
    end

    println("Read stage finished successfully.")

# Process "simulate" command
elseif command == "simulate"
    if length(ARGS) < 2
        println("Error: Missing folder_path for simulate command.")
        exit(1)
    end

    folder_path = ARGS[2]
    frames_csv_path = joinpath(folder_path, "frames.csv")
    extra_args = ARGS[3:end]

    # Try NCBI frame reading first
    ncbi_success = safe_run(`julia --project=. src/read_ncbi_frames.jl $folder_path`)
    if !ncbi_success || !isfile(frames_csv_path)
        println("read_ncbi_frames.jl failed or frames.csv not found. Trying read_samfire_frames.jl instead.")
        samfire_success = safe_run(`julia --project=. src/read_samfire_frames.jl $folder_path`)
        if !samfire_success || !isfile(frames_csv_path)
            println("Error: Both frame-reading methods failed.")
            exit(1)
        end
    end

    # Forward extra arguments to simulate_variants.jl
    local sim_cmd = `julia --project=. src/simulate_variants.jl $folder_path`
    for arg in extra_args
        sim_cmd = `$sim_cmd $arg`
    end
    if !safe_run(sim_cmd)
        println("Error running src/simulate_variants.jl")
        exit(1)
    end

    println("Simulate stage finished successfully.")

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

    println("run_supertype method finished successfully.")

  

else
    println("Error: Invalid command '$command'.")
    print_help()
    exit(1)
end