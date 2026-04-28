#!/usr/bin/env julia

###############################################################################
# CD8scape.jl
#
# USAGE:
#   ./CD8scape.jl prep
#   ./CD8scape.jl read        <folder_path> [--aa] [--suffix <name>] [--latest|--no-latest]
#   ./CD8scape.jl simulate    <folder_path> [--n <count>] [--p <proportion>] [--seed <int>] [--suffix <name>] [--latest|--no-latest]
#   ./CD8scape.jl run         <folder_path> [--t <N|max>] [--per-allele] [--verbose] [--suffix <name>] [--latest|--no-latest]
#   ./CD8scape.jl run_supertype <folder_path> [--t <N|max>] [--per-allele] [--verbose] [--suffix <name>] [--latest|--no-latest]
#   ./CD8scape.jl percentile  <folder_path> [--s <sim_file>] [--o <obs_file>] [--per-allele]
#
###############################################################################

function print_help()
    println("""
CD8scape.jl - A tool for running netMHCpan and managing related data.

USAGE:
    ./CD8scape.jl prep
    ./CD8scape.jl read <folder_path> [--aa] [--suffix <name>] [--latest|--no-latest]
    ./CD8scape.jl simulate <folder_path> [--n <count>] [--p <proportion>] [--seed <int>] [--suffix <name>] [--latest|--no-latest]
    ./CD8scape.jl run  <folder_path> [--t <N|max>|--thread <N|max>] [--per-allele] [--verbose] [--suffix <name>] [--latest|--no-latest]
    ./CD8scape.jl run_supertype  <folder_path> [--t <N|max>|--thread <N|max>] [--per-allele] [--verbose] [--suffix <name>] [--latest|--no-latest]
    ./CD8scape.jl percentile <folder_path> [--s <sim_file>] [--o <obs_file>] [--per-allele]
  

COMMANDS:
    prep         Set up the environment by running src/env.jl.
    read         Parse variants input (Samfire trajectories, .vcf, or .aa amino-acid variants) and read frames (SamFire or NCBI).
    simulate     Read frames, then generate and optionally sample simulated variants per frame.
    run          Run the peptide-generation and NetMHCpan pipeline on parsed data.
    run_supertype Run the peptide-generation and NetMHCpan pipeline on parsed data for a representative supertype HLA panel.
    percentile  Compute observed HMBR log2 fold-change percentiles relative to simulated distribution.
  

OPTIONS:
  --help, -h
      Print this help message and exit.
  --aa
      Parse amino-acid-level variants from a .aa file instead of VCF or trajectories.
      The .aa file format is two lines per variant:
        <orf_name> <aa_position>
        <ancestral_aa> <derived_aa>
      The orf_name must match a Description in frames.csv; aa_position is 1-based.
  --suffix <name>
      Append _<name> before the file extension of all output files (e.g. --suffix foo
      produces variants_foo.csv, best_ranks_foo.csv, etc.). For 'simulate', defaults to
      'simulated' when omitted.
  --latest (default), --no-latest
      When reading input files without a suffix, if multiple candidates exist (e.g.
      frames.csv and frames_simulated.csv), --latest picks the most recently modified
      file. --no-latest errors on ambiguity instead.
  --t <N|max>, --thread <N|max>
      Max number of chunks to execute in parallel when running NetMHCpan (default 1). Use 'max' to use the safety cap.
  --per-allele
      Perform the log2 fold change calculation per allele in the provided genome,
      writing a separate per_allele_best_ranks.csv file containing one row per
      (Locus, Mutation, allele) with columns:
        Frame              - protein / region label
        MHC                - the allele identifier
        ELBR_A            - best ancestral EL rank for that allele
        ELBR_D            - best derived EL rank for that allele
        foldchange_BR      - ELBR_D / ELBR_A
        log2_foldchange_BR - log2(foldchange_BR)
      Only alleles where ancestral EL rank ≤ 2% are included.
      Available for run, run_supertype, and percentile. For percentile, switches from
      harmonic_mean_best_ranks to per_allele_best_ranks files and uses log2_foldchange_BR
      as the distribution column; output is percentile_per_allele_best_ranks.csv.
      Note: for run_supertype, the panel alleles are population-frequency surrogates
      rather than an individual's genotype.
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
    use_aa = false
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
        elseif a == "--aa"
            use_aa = true
        end
        i += 1
    end
    return suffix, latest, use_aa
end

# Process "prep" command
if command == "prep"
    # Setup environment
    if !safe_run(`julia --startup-file=no --project=. src/env.jl`)
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
    suffix, latest, use_aa = parse_suffix_latest(extra_args)
    # Expected output path for frames (do not require existence before running readers)
    frames_csv_path = resolve_write(joinpath(folder_path, "frames.csv"); suffix=suffix)
    variants_csv_path = resolve_write(joinpath(folder_path, "variants.csv"); suffix=suffix)

    # Read frames (NCBI first, then Samfire fallback)
    local read_ncbi_cmd = `julia --startup-file=no --project=. src/read_ncbi_frames.jl $folder_path`
    if suffix != ""; read_ncbi_cmd = `$read_ncbi_cmd --suffix $suffix`; end
    if latest; read_ncbi_cmd = `$read_ncbi_cmd --latest`; else read_ncbi_cmd = `$read_ncbi_cmd --no-latest`; end
    ncbi_success = safe_run(read_ncbi_cmd)
    if !ncbi_success || !isfile(frames_csv_path)
        println("read_ncbi_frames.jl failed or frames.csv not found. Trying read_samfire_frames.jl instead.")
        local read_sam_cmd = `julia --startup-file=no --project=. src/read_samfire_frames.jl $folder_path`
        if suffix != ""; read_sam_cmd = `$read_sam_cmd --suffix $suffix`; end
        if latest; read_sam_cmd = `$read_sam_cmd --latest`; else read_sam_cmd = `$read_sam_cmd --no-latest`; end
        samfire_success = safe_run(read_sam_cmd)
        if !samfire_success || !isfile(frames_csv_path)
            println("Error: Both frame-reading methods failed.")
            exit(1)
        end
    end

    # Parse variants
    parse_ok = false
    if use_aa
        # --aa flag: parse amino-acid-level variants from a .aa file
        local parse_aa_cmd = `julia --startup-file=no --project=. src/parse_aa_variants.jl $folder_path`
        if suffix != ""; parse_aa_cmd = `$parse_aa_cmd --suffix $suffix`; end
        if latest; parse_aa_cmd = `$parse_aa_cmd --latest`; else parse_aa_cmd = `$parse_aa_cmd --no-latest`; end
        parse_ok = safe_run(parse_aa_cmd)
        if !parse_ok
            println("Error: parse_aa_variants.jl failed.")
        end
    else
        # Prefer VCF if present, otherwise trajectories
        vcf_files = filter(f -> endswith(f, ".vcf") || endswith(f, ".vcf.gz"), readdir(folder_path; join=true))
        if !isempty(vcf_files)
            local parse_vcf_cmd = `julia --startup-file=no --project=. src/parse_vcf.jl $folder_path`
            if suffix != ""; parse_vcf_cmd = `$parse_vcf_cmd --suffix $suffix`; end
            parse_ok = safe_run(parse_vcf_cmd)
            if !parse_ok
                println("parse_vcf.jl failed; trying parse_trajectories.jl.")
            end
        end
        if !parse_ok
            local parse_traj_cmd = `julia --startup-file=no --project=. src/parse_trajectories.jl $folder_path`
            if suffix != ""; parse_traj_cmd = `$parse_traj_cmd --suffix $suffix`; end
            parse_ok = safe_run(parse_traj_cmd)
        end
    end
    if !parse_ok || !isfile(variants_csv_path)
        println("Error: Failed to parse variants.")
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
    suffix, latest, _ = parse_suffix_latest(extra_args)
    # If no suffix provided, default to 'simulated' so frames/variants share it
    if isempty(suffix)
        suffix = "simulated"
    end
    # Expected output path for frames (created by readers)
    frames_csv_path = resolve_write(joinpath(folder_path, "frames.csv"); suffix=suffix)

    # Try NCBI frame reading first
    local read_ncbi_cmd = `julia --startup-file=no --project=. src/read_ncbi_frames.jl $folder_path`
    if suffix != ""; read_ncbi_cmd = `$read_ncbi_cmd --suffix $suffix`; end
    if latest; read_ncbi_cmd = `$read_ncbi_cmd --latest`; else read_ncbi_cmd = `$read_ncbi_cmd --no-latest`; end
    ncbi_success = safe_run(read_ncbi_cmd)
    if !ncbi_success || !isfile(frames_csv_path)
        println("read_ncbi_frames.jl failed or frames.csv not found. Trying read_samfire_frames.jl instead.")
        local read_sam_cmd = `julia --startup-file=no --project=. src/read_samfire_frames.jl $folder_path`
        if suffix != ""; read_sam_cmd = `$read_sam_cmd --suffix $suffix`; end
        if latest; read_sam_cmd = `$read_sam_cmd --latest`; else read_sam_cmd = `$read_sam_cmd --no-latest`; end
        samfire_success = safe_run(read_sam_cmd)
        if !samfire_success || !isfile(frames_csv_path)
            println("Error: Both frame-reading methods failed.")
            exit(1)
        end
    end

    # Forward extra arguments to simulate_variants.jl
    local sim_cmd = `julia --startup-file=no --project=. src/simulate_variants.jl $folder_path`
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
    suffix, latest, _ = parse_suffix_latest(extra_args)
    verbose = any(a -> a == "--verbose", extra_args)
    per_allele = any(a -> a == "--per-allele", extra_args)
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
        local gen_cmd = `julia --startup-file=no --project=. src/generate_peptides.jl $folder_path`
        if suffix != ""; gen_cmd = `$gen_cmd --suffix $suffix`; end
        if latest; gen_cmd = `$gen_cmd --latest`; else gen_cmd = `$gen_cmd --no-latest`; end
        if !safe_run(gen_cmd)
        println("Error running src/generate_peptides.jl")
        exit(1)
    end

    # Clean Peptides
    if !safe_run(`julia --startup-file=no --project=. src/clean_peptides.jl $folder_path`)
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
    local run_cmd = `julia --startup-file=no --project=. src/run_netMHCpan.jl --folder $folder_path`
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
    #
    # NOTE (memory): Stream the Perl script's stdout straight to disk instead of
    # buffering the entire output in RAM via `read(..., String)`. For large
    # NetMHCpan outputs (e.g. ~10 GB) the processed CSV is similarly large, and
    # materialising it as a single Julia string causes an OOM crash.
    try
        run(pipeline(`perl src/process_output.pl $netmhcpan_output`, stdout = processed_output))
    catch e
        println("Error running src/process_output.pl: $e")
        exit(1)
    end

    # If processed_output has no data rows, skip remaining stages gracefully.
    # Avoid readlines() here (also buffers the whole file); instead count newlines
    # while streaming.
    try
        data_rows = 0
        open(processed_output, "r") do io
            # Skip header line if present
            isheader = !eof(io)
            isheader && readline(io)
            while !eof(io)
                readline(io)
                data_rows += 1
                # Early exit: once we know there is data, stop counting.
                data_rows >= 1 && break
            end
        end
        if data_rows == 0
            println("Skipping downstream processing: processed_output has 0 data rows.")
            println("Run stage finished successfully (no data).")
            exit(0)
        end
    catch
        # If reading fails, continue; subsequent steps may handle errors
    end

    # Process Scores
    local scores_cmd = `julia --startup-file=no --project=. src/process_scores.jl --folder $folder_path`
    if suffix != ""; scores_cmd = `$scores_cmd --suffix $suffix`; end
    if latest; scores_cmd = `$scores_cmd --latest`; else scores_cmd = `$scores_cmd --no-latest`; end
    if !safe_run(scores_cmd)
        println("Error running src/process_scores.jl")
        exit(1)
    end

    # Process Best Ranks
    local ranks_cmd = `julia --startup-file=no --project=. src/process_best_ranks.jl $folder_path`
    if suffix != ""; ranks_cmd = `$ranks_cmd --suffix $suffix`; end
    if latest; ranks_cmd = `$ranks_cmd --latest`; else ranks_cmd = `$ranks_cmd --no-latest`; end
    if per_allele; ranks_cmd = `$ranks_cmd --per-allele`; end
    if !safe_run(ranks_cmd)
        println("Error running src/process_best_ranks.jl")
        exit(1)
    end

    # Variant fates — trace each variant through pipeline filter stages
    local fates_cmd = `julia --startup-file=no --project=. src/variant_fates.jl $folder_path`
    if suffix != ""; fates_cmd = `$fates_cmd --suffix $suffix`; end
    if latest; fates_cmd = `$fates_cmd --latest`; else fates_cmd = `$fates_cmd --no-latest`; end
    if !safe_run(fates_cmd)
        println("Warning: src/variant_fates.jl failed — variant fates summary will not be available.")
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
    suffix, latest, _ = parse_suffix_latest(extra_args)
    verbose = any(a -> a == "--verbose", extra_args)
    per_allele = any(a -> a == "--per-allele", extra_args)
    if per_allele
        println("Note: --per-allele is active on a supertype panel run. Supertype panel alleles are population-frequency surrogates, not an individual's genotype.")
    end
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
    local gen_cmd = `julia --startup-file=no --project=. src/generate_peptides.jl $folder_path`
    if suffix != ""; gen_cmd = `$gen_cmd --suffix $suffix`; end
    if latest; gen_cmd = `$gen_cmd --latest`; else gen_cmd = `$gen_cmd --no-latest`; end
    if !safe_run(gen_cmd)
        println("Error running src/generate_peptides.jl")
        exit(1)
    end

    # Clean Peptides
    if !safe_run(`julia --startup-file=no --project=. src/clean_peptides.jl $folder_path`)
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
    local run_cmd = `julia --startup-file=no --project=. src/run_netMHCpan_global.jl --folder $folder_path`
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
    #
    # NOTE (memory): Stream the Perl script's stdout straight to disk instead of
    # buffering the entire output in RAM via `read(..., String)`. For the
    # supertype panel (~2000 HLA alleles) the processed CSV can reach tens of GB,
    # so materialising it as a single Julia string causes an OOM crash.
    try
        run(pipeline(`perl src/process_output.pl $netmhcpan_output`, stdout = processed_output))
    catch e
        println("Error running src/process_output.pl: $e")
        exit(1)
    end

    # If processed_output has no data rows, skip remaining stages gracefully.
    # Avoid readlines() here (also buffers the whole file); stream instead.
    try
        data_rows = 0
        open(processed_output, "r") do io
            isheader = !eof(io)
            isheader && readline(io)
            while !eof(io)
                readline(io)
                data_rows += 1
                # Early exit: once we know there is data, stop counting.
                data_rows >= 1 && break
            end
        end
        if data_rows == 0
            println("Skipping downstream processing: processed_output has 0 data rows.")
            println("run_supertype method finished successfully (no data).")
            exit(0)
        end
    catch
        # If reading fails, continue; subsequent steps may handle errors
    end

    # Process Scores
    local scores_cmd = `julia --startup-file=no --project=. src/process_scores.jl --folder $folder_path`
    if suffix != ""; scores_cmd = `$scores_cmd --suffix $suffix`; end
    if latest; scores_cmd = `$scores_cmd --latest`; else scores_cmd = `$scores_cmd --no-latest`; end
    if !safe_run(scores_cmd)
        println("Error running src/process_scores.jl")
        exit(1)
    end

    # Process Best Ranks
    local ranks_cmd = `julia --startup-file=no --project=. src/process_best_ranks.jl $folder_path`
    if suffix != ""; ranks_cmd = `$ranks_cmd --suffix $suffix`; end
    if latest; ranks_cmd = `$ranks_cmd --latest`; else ranks_cmd = `$ranks_cmd --no-latest`; end
    if per_allele; ranks_cmd = `$ranks_cmd --per-allele`; end
    if !safe_run(ranks_cmd)
        println("Error running src/process_best_ranks.jl")
        exit(1)
    end

    # Variant fates — trace each variant through pipeline filter stages
    local fates_cmd = `julia --startup-file=no --project=. src/variant_fates.jl $folder_path`
    if suffix != ""; fates_cmd = `$fates_cmd --suffix $suffix`; end
    if latest; fates_cmd = `$fates_cmd --latest`; else fates_cmd = `$fates_cmd --no-latest`; end
    if !safe_run(fates_cmd)
        println("Warning: src/variant_fates.jl failed — variant fates summary will not be available.")
    end

    println("run_supertype method finished successfully.")

  

elseif command == "percentile"
    if length(ARGS) < 2
        println("Error: Missing folder_path for percentile command.")
        exit(1)
    end

    folder_path = ARGS[2]
    extra_args = ARGS[3:end]
    # Wrap parsing in a helper to avoid soft-scope warnings
    function _parse_s_o(argv::Vector{String})
        fwd = String[]
        i = 1
        while i <= length(argv)
            a = argv[i]
            if a == "--s" || a == "--o"
                push!(fwd, a)
                if i < length(argv) && !startswith(argv[i+1], "--")
                    push!(fwd, argv[i+1])
                    i += 1
                end
            end
            i += 1
        end
        return fwd
    end
    forward = _parse_s_o(extra_args)
    per_allele = any(a -> a == "--per-allele", extra_args)

    local pct_cmd = `julia --startup-file=no --project=. src/percentile.jl $folder_path`
    for f in forward
        pct_cmd = `$pct_cmd $f`
    end
    if per_allele; pct_cmd = `$pct_cmd --per-allele`; end
    if !safe_run(pct_cmd)
        println("Error running src/percentile.jl")
        exit(1)
    end

    println("Percentile stage finished successfully.")

else
    println("Error: Invalid command '$command'.")
    print_help()
    exit(1)
end