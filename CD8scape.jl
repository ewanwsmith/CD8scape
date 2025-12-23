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
    ./CD8scape.jl simulate <folder_path> [--n <count>] [--p <proportion>] [--seed <int>]
    ./CD8scape.jl run  <folder_path> [--t <N|max>|--thread <N|max>] [--verbose]
    ./CD8scape.jl run_supertype  <folder_path> [--t <N|max>|--thread <N|max>] [--verbose]
  

COMMANDS:
    prep         Set up the environment by running src/env.jl.
    read         Parse variants input (Samfire trajectories or .vcf) and read frames (SamFire or NCBI).
    simulate     Read frames, then generate and optionally sample simulated variants per frame.
    run          Run the peptide-generation and NetMHCpan pipeline on parsed data.
    run_supertype Run the peptide-generation and NetMHCpan pipeline on parsed data for a representative supertype HLA panel.
  

OPTIONS:
  --help, -h
      Print this help message and exit.
  --t <N|max>, --thread <N|max>
      Max number of chunks to execute in parallel when running NetMHCpan (default 1). Use 'max' to use the safety cap.
  --verbose
      Preserve per-allele logs and temp files for debugging.
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

# Include path helpers for suffix/latest handling
include("src/path_utils.jl")

# Extract common flags from extra arguments
function parse_suffix_latest(argv::Vector{String})
    suffix = ""
    latest = true
    i = 1
    while i <= length(argv)
        a = argv[i]
        if a == "--suffix"
            if i + 1 <= length(argv) && !startswith(argv[i+1], "--")
                i += 1
                suffix = argv[i]
            end
        elseif a == "--latest"
            latest = true
        elseif a == "--no-latest"
            latest = false
        end
        i += 1
    end
    return suffix, latest
end

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
    extra_args = ARGS[3:end]
    suffix, latest = parse_suffix_latest(extra_args)
    # Expected output path for frames (do not require existence before running readers)
    frames_csv_path = resolve_write(joinpath(folder_path, "frames.csv"); suffix=suffix)
    variants_csv_path = resolve_write(joinpath(folder_path, "variants.csv"); suffix=suffix)

    # Read frames (NCBI first, then Samfire fallback)
    local read_ncbi_cmd = `julia --project=. src/read_ncbi_frames.jl $folder_path`
    if suffix != ""; read_ncbi_cmd = `$read_ncbi_cmd --suffix $suffix`; end
    if latest; read_ncbi_cmd = `$read_ncbi_cmd --latest`; else read_ncbi_cmd = `$read_ncbi_cmd --no-latest`; end
    ncbi_success = safe_run(read_ncbi_cmd)
    if !ncbi_success || !isfile(frames_csv_path)
        println("read_ncbi_frames.jl failed or frames.csv not found. Trying read_samfire_frames.jl instead.")
        local read_sam_cmd = `julia --project=. src/read_samfire_frames.jl $folder_path`
        if suffix != ""; read_sam_cmd = `$read_sam_cmd --suffix $suffix`; end
        if latest; read_sam_cmd = `$read_sam_cmd --latest`; else read_sam_cmd = `$read_sam_cmd --no-latest`; end
        samfire_success = safe_run(read_sam_cmd)
        if !samfire_success || !isfile(frames_csv_path)
            println("Error: Both frame-reading methods failed.")
            exit(1)
        end
    end

    # Parse variants: prefer VCF if present, otherwise trajectories
    vcf_files = filter(f -> endswith(f, ".vcf") || endswith(f, ".vcf.gz"), readdir(folder_path; join=true))
    parse_ok = false
    if !isempty(vcf_files)
            local parse_vcf_cmd = `julia --project=. src/parse_vcf.jl $folder_path`
            if suffix != ""; parse_vcf_cmd = `$parse_vcf_cmd --suffix $suffix`; end
            parse_ok = safe_run(parse_vcf_cmd)
        if !parse_ok
            println("parse_vcf.jl failed; trying parse_trajectories.jl.")
        end
    end
    if !parse_ok
            local parse_traj_cmd = `julia --project=. src/parse_trajectories.jl $folder_path`
            if suffix != ""; parse_traj_cmd = `$parse_traj_cmd --suffix $suffix`; end
            parse_ok = safe_run(parse_traj_cmd)
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
    extra_args = ARGS[3:end]
    suffix, latest = parse_suffix_latest(extra_args)
    # If no suffix provided, default to 'simulated' so frames/variants share it
    if isempty(suffix)
        suffix = "simulated"
    end
    # Expected output path for frames (created by readers)
    frames_csv_path = resolve_write(joinpath(folder_path, "frames.csv"); suffix=suffix)

    # Try NCBI frame reading first
    local read_ncbi_cmd = `julia --project=. src/read_ncbi_frames.jl $folder_path`
    if suffix != ""; read_ncbi_cmd = `$read_ncbi_cmd --suffix $suffix`; end
    if latest; read_ncbi_cmd = `$read_ncbi_cmd --latest`; else read_ncbi_cmd = `$read_ncbi_cmd --no-latest`; end
    ncbi_success = safe_run(read_ncbi_cmd)
    if !ncbi_success || !isfile(frames_csv_path)
        println("read_ncbi_frames.jl failed or frames.csv not found. Trying read_samfire_frames.jl instead.")
        local read_sam_cmd = `julia --project=. src/read_samfire_frames.jl $folder_path`
        if suffix != ""; read_sam_cmd = `$read_sam_cmd --suffix $suffix`; end
        if latest; read_sam_cmd = `$read_sam_cmd --latest`; else read_sam_cmd = `$read_sam_cmd --no-latest`; end
        samfire_success = safe_run(read_sam_cmd)
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
    extra_args = ARGS[3:end]
    suffix, latest = parse_suffix_latest(extra_args)
    verbose = any(a -> a == "--verbose", extra_args)
    # Extract optional thread count and pass through
    threads_arg = String[]
    for (i, a) in enumerate(extra_args)
        if a == "--t" || a == "--thread"
            if i <= length(extra_args) - 1
                push!(threads_arg, a)
                push!(threads_arg, extra_args[i+1])
            end
        end
    end
    # Expected output for NetMHCpan runner (construct path without requiring existence)
        netmhcpan_output = resolve_write(joinpath(folder_path, "netmhcpan_output.tsv"); suffix=suffix)
        processed_output = resolve_write(joinpath(folder_path, "processed_output.csv"); suffix=suffix)
    skip_marker = joinpath(folder_path, ".cd8scape_skipped")

    # Generate Peptides
        local gen_cmd = `julia --project=. src/generate_peptides.jl $folder_path`
        if suffix != ""; gen_cmd = `$gen_cmd --suffix $suffix`; end
        if latest; gen_cmd = `$gen_cmd --latest`; else gen_cmd = `$gen_cmd --no-latest`; end
        if !safe_run(gen_cmd)
        println("Error running src/generate_peptides.jl")
        exit(1)
    end

    # Clean Peptides
    if !safe_run(`julia --project=. src/clean_peptides.jl $folder_path`)
        println("Error running src/clean_peptides.jl")
        exit(1)
    end

    # If no peptides remain after cleaning, skip gracefully (post-clean stage)
    try
        pep_file = joinpath(folder_path, "Peptides.pep")
        if isfile(pep_file)
            pep_lines = readlines(pep_file)
            if isempty(pep_lines)
                println("Skipping downstream processing: no peptides remain after cleaning.")
                println("Run stage finished successfully (no data).")
                exit(0)
            end
        end
    catch
        # Continue if read fails; subsequent steps will handle
    end

    # Run NetMHCpan
    local run_cmd = `julia --project=. src/run_netMHCpan.jl --folder $folder_path`
    for t in threads_arg
        run_cmd = `$run_cmd $t`
    end
    if verbose
        run_cmd = `$run_cmd --verbose`
    end
    if suffix != ""; run_cmd = `$run_cmd --suffix $suffix`; end
    if latest; run_cmd = `$run_cmd --latest`; else run_cmd = `$run_cmd --no-latest`; end
    if !safe_run(run_cmd)
        println("Error: NetMHCpan did not run successfully.")
        exit(1)
    end
    # If the runner signaled skip, end process gracefully now (post-NetMHCpan)
    if isfile(skip_marker)
        println("Skipping downstream processing due to runner skip marker.")
        # Clean up marker to avoid future confusion
        try; rm(skip_marker; force=true); catch; end
        println("Run stage finished successfully (no data).")
        exit(0)
    end

    # Verify NetMHCpan Output
    if !isfile(netmhcpan_output)
        println("Error: Expected NetMHCpan output file '$netmhcpan_output' does not exist.")
        exit(1)
    end
    # Continue regardless of number of lines; downstream scripts will handle empty outputs

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

    # If processed_output has no data rows, skip remaining stages gracefully
    try
        lines = readlines(processed_output)
        data_rows = length(lines) > 1 ? (length(lines) - 1) : 0
        if data_rows == 0
            println("Skipping downstream processing: processed_output has 0 data rows.")
            println("Run stage finished successfully (no data).")
            exit(0)
        end
    catch
        # If reading fails, continue; subsequent steps may handle errors
    end

    # Process Scores
    local scores_cmd = `julia --project=. src/process_scores.jl --folder $folder_path`
    if suffix != ""; scores_cmd = `$scores_cmd --suffix $suffix`; end
    if latest; scores_cmd = `$scores_cmd --latest`; else scores_cmd = `$scores_cmd --no-latest`; end
    if !safe_run(scores_cmd)
        println("Error running src/process_scores.jl")
        exit(1)
    end

    # Process Best Ranks
    local ranks_cmd = `julia --project=. src/process_best_ranks.jl $folder_path`
    if suffix != ""; ranks_cmd = `$ranks_cmd --suffix $suffix`; end
    if latest; ranks_cmd = `$ranks_cmd --latest`; else ranks_cmd = `$ranks_cmd --no-latest`; end
    if !safe_run(ranks_cmd)
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
    extra_args = ARGS[3:end]
    suffix, latest = parse_suffix_latest(extra_args)
    verbose = any(a -> a == "--verbose", extra_args)
    # Extract optional thread count and pass through
    threads_arg = String[]
    for (i, a) in enumerate(extra_args)
        if a == "--t" || a == "--thread"
            if i <= length(extra_args) - 1
                push!(threads_arg, a)
                push!(threads_arg, extra_args[i+1])
            end
        end
    end
    netmhcpan_output = resolve_write(joinpath(folder_path, "netmhcpan_output.tsv"); suffix=suffix)
    processed_output = resolve_write(joinpath(folder_path, "processed_output.csv"); suffix=suffix)
    # Define skip marker path for graceful early exit
    skip_marker = joinpath(folder_path, ".cd8scape_skipped")

    # Generate Peptides
    local gen_cmd = `julia --project=. src/generate_peptides.jl $folder_path`
    if suffix != ""; gen_cmd = `$gen_cmd --suffix $suffix`; end
    if latest; gen_cmd = `$gen_cmd --latest`; else gen_cmd = `$gen_cmd --no-latest`; end
    if !safe_run(gen_cmd)
        println("Error running src/generate_peptides.jl")
        exit(1)
    end

    # Clean Peptides
    if !safe_run(`julia --project=. src/clean_peptides.jl $folder_path`)
        println("Error running src/clean_peptides.jl")
        exit(1)
    end

    # If no peptides remain after cleaning, skip gracefully (post-clean stage)
    try
        pep_file = joinpath(folder_path, "Peptides.pep")
        if isfile(pep_file)
            pep_lines = readlines(pep_file)
            if isempty(pep_lines)
                println("Skipping downstream processing: no peptides remain after cleaning.")
                println("run_supertype method finished successfully (no data).")
                exit(0)
            end
        end
    catch
        # Continue if read fails; subsequent steps will handle
    end

    # Run NetMHCpan
    local run_cmd = `julia --project=. src/run_netMHCpan_global.jl --folder $folder_path`
    for t in threads_arg
        run_cmd = `$run_cmd $t`
    end
    if verbose
        run_cmd = `$run_cmd --verbose`
    end
    if suffix != ""; run_cmd = `$run_cmd --suffix $suffix`; end
    if latest; run_cmd = `$run_cmd --latest`; else run_cmd = `$run_cmd --no-latest`; end
    if !safe_run(run_cmd)
        println("Error: NetMHCpan did not run successfully.")
        exit(1)
    end
    # If the runner signaled skip, end process gracefully now (post-NetMHCpan)
    if isfile(skip_marker)
        println("Skipping downstream processing due to runner skip marker.")
        # Clean up marker to avoid future confusion
        try; rm(skip_marker; force=true); catch; end
        println("run_supertype method finished successfully (no data).")
        exit(0)
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

    # If processed_output has no data rows, skip remaining stages gracefully
    try
        lines = readlines(processed_output)
        data_rows = length(lines) > 1 ? (length(lines) - 1) : 0
        if data_rows == 0
            println("Skipping downstream processing: processed_output has 0 data rows.")
            println("run_supertype method finished successfully (no data).")
            exit(0)
        end
    catch
        # If reading fails, continue; subsequent steps may handle errors
    end

    # Process Scores
    local scores_cmd = `julia --project=. src/process_scores.jl --folder $folder_path`
    if suffix != ""; scores_cmd = `$scores_cmd --suffix $suffix`; end
    if latest; scores_cmd = `$scores_cmd --latest`; else scores_cmd = `$scores_cmd --no-latest`; end
    if !safe_run(scores_cmd)
        println("Error running src/process_scores.jl")
        exit(1)
    end

    # Process Best Ranks
    local ranks_cmd = `julia --project=. src/process_best_ranks.jl $folder_path`
    if suffix != ""; ranks_cmd = `$ranks_cmd --suffix $suffix`; end
    if latest; ranks_cmd = `$ranks_cmd --latest`; else ranks_cmd = `$ranks_cmd --no-latest`; end
    if !safe_run(ranks_cmd)
        println("Error running src/process_best_ranks.jl")
        exit(1)
    end

    println("run_supertype method finished successfully.")

  

else
    println("Error: Invalid command '$command'.")
    print_help()
    exit(1)
end