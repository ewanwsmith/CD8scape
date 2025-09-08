#!/usr/bin/env julia
"""
This script orchestrates the full context peptide pipeline.

Purpose:
- Runs all context-specific processing steps in order.
- Handles argument parsing and error checking.
- Passes options to each step as needed.

Usage:
    julia context_run.jl <folder_path> [--n_loci <number_of_loci>] [--mode panel|supertype] [--override]

Arguments:
    <folder_path>   Path to the folder containing input files.
    --n_loci        Number of loci to simulate (default: 1000).
    --mode          Run mode: panel or supertype (default: panel).
    --override      Optional flag to override existing outputs.
"""

using CSV, DataFrames, FilePathsBase


function main()


    if !isdefined(Main, :ARGS) || isempty(ARGS)
        println("Error: folder_path must be provided as a command-line argument.")
        exit(1)
    end

    folder_path = ARGS[1]
    frames_path = joinpath(folder_path, "frames.csv")

    if !isfile(frames_path)
        println("Error: frames.csv not found in the provided folder path: $folder_path")
        exit(1)
    end

    # Default number of loci
    n_loci = 1000
    mode = "panel"
    override = false

    # Parse additional arguments
    for (i, arg) in enumerate(ARGS)
        if arg == "--n_loci" && i < length(ARGS)
            try
                n_loci = parse(Int, ARGS[i + 1])
            catch
                println("Error: Invalid value for --n_loci.")
                exit(1)
            end
        elseif arg == "--mode" && i < length(ARGS)
            mode = ARGS[i + 1]
            if mode ∉ ["panel", "supertype"]
                println("Error: Invalid value for --mode.")
                exit(1)
            end
        elseif arg == "--override"
            override = true
        end
    end

    # Call generate_context_peptides.jl
    generate_script = joinpath(@__DIR__, "generate_context_peptides.jl")
    local cmd = `julia $generate_script --folder $folder_path --n_loci $n_loci`

    run(cmd)

    # Case-insensitive lookup of context_peptides.pep
    matches = filter(f -> lowercase(basename(f)) == "context_peptides.pep", readdir(folder_path; join=true))
    if isempty(matches)
        println("Error: context_peptides.pep not found (case-insensitive) in folder $folder_path")
        exit(1)
    end


    # Call clean_peptides_context.jl
    clean_script = joinpath(@__DIR__, "clean_peptides_context.jl")
    local clean_cmd = `julia $clean_script $folder_path`
    run(clean_cmd)

    # Run NetMHCpan based on mode (no --prefix context_)
    netmhc_script = joinpath(@__DIR__, "run_netMHCpan_context.jl")
    local netmhc_cmd = `julia $netmhc_script --folder $folder_path --peptides $(first(matches)) --mode $mode`
    run(netmhc_cmd)


    # Post-process NetMHCpan output using context-specific Perl script
    netmhcpan_output = joinpath(folder_path, "netMHCpan_output.tsv")
    processed_output = joinpath(folder_path, "context_processed_netMHCpan_output.csv")

    if !isfile(netmhcpan_output)
        println("Error: Expected NetMHCpan output file 'netMHCpan_output.tsv' does not exist.")
        exit(1)
    end

    try
        run(`perl $(joinpath(@__DIR__, "process_output_context.pl")) $netmhcpan_output`)
    catch e
        println("Error running process_output_context.pl: $e")
        exit(1)
    end


    # Run context-specific score processing
    if !success(`julia --project=. $(joinpath(@__DIR__, "process_scores_context.jl")) $folder_path`)
        println("Error running process_scores_context.jl")
        exit(1)
    end

    # Run appropriate best ranks script based on mode
    if mode == "panel"
        best_ranks_context_script = joinpath(@__DIR__, "process_best_ranks_context.jl")
        best_ranks_context_cmd = `julia --project=. $best_ranks_context_script $folder_path`
        run(best_ranks_context_cmd)
    elseif mode == "supertype"
        best_ranks_context_script = joinpath(@__DIR__, "process_best_ranks_context.jl")
        best_ranks_context_cmd = `julia --project=. $best_ranks_context_script $folder_path --supertype`
        run(best_ranks_context_cmd)
    end


end

main()