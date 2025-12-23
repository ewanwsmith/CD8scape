#!/usr/bin/env julia
"""
process_best_ranks.jl

Processes best ranks for ancestral and derived peptides for each locus.

Usage:
    julia process_best_ranks.jl <folder_path> [--suffix <name>] [--latest|--no-latest]

Arguments:
    <folder_path>   Path to the folder containing processed_peptides.csv.
"""

using DataFrames, CSV, Statistics, StatsBase
include("path_utils.jl")

if abspath(PROGRAM_FILE) == @__FILE__
    # Parse positional folder and optional flags
    if length(ARGS) < 1
        println("Usage: julia process_best_ranks.jl <folder_path> [--suffix <name>] [--latest]")
        exit(1)
    end

    folder_path = ARGS[1]
    # Helper to parse flags from the remainder of ARGS to avoid soft-scope warnings
    function _parse_suffix_latest(argv::Vector{String})
        sfx = ""
        lat = true
        i = 1
        while i <= length(argv)
            a = argv[i]
            if a == "--suffix"
                if i + 1 <= length(argv) && !startswith(argv[i+1], "--")
                    i += 1
                    sfx = argv[i]
                end
            elseif a == "--latest"
                lat = true
            elseif a == "--no-latest"
                lat = false
            end
            i += 1
        end
        return sfx, lat
    end
    suffix, latest = _parse_suffix_latest(ARGS[2:end])

    input_file = resolve_read(joinpath(folder_path, "processed_peptides.csv"); suffix=suffix, latest=latest)
    println("Reading input file: $input_file")

    # Read the input CSV into a DataFrame with error handling
    try
        global df = CSV.read(input_file, DataFrame)
        println("Successfully loaded data with $(nrow(df)) rows")
    catch e
        println("Error reading input file: $e")
        exit(1)
    end

# Function to find minimum EL_Rank by Locus, MHC, and AA change for given peptide pattern
function find_best_ranks(df, pattern)
    subset = filter(row -> !ismissing(row.Peptide_label) && endswith(row.Peptide_label, pattern), df)
    if isempty(subset)
        println("Warning: No peptides found for pattern '$pattern'")
        return DataFrame(Locus = Int[], MHC = String[], Change = String[], Best_EL_Rank = Float64[], Peptide_Type = String[], Description = String[], Sequence = String[])
    end
    # Coerce EL_Rank column to Float64 if not already
    if !(eltype(subset.EL_Rank) <: AbstractFloat)
        subset.EL_Rank = map(x -> try
            x isa Number ? float(x) : parse(Float64, String(x))
        catch
            NaN
        end, subset.EL_Rank)
    end
    # Derive AA change key from Peptide_label (e.g., "A17T" from "A17T_ProteinName_3_ D")
    change_keys = replace.(subset.Peptide_label, r"_(C|V|A|D)$" => "")
    change_keys = replace.(change_keys, r"_\d+$" => "")
    change_keys = replace.(change_keys, ' ' => '_')
    subset[!, :Mutation] = String.(first.(split.(change_keys, "_")))

    grouped = groupby(subset, [:Locus, :MHC, :Mutation])
    best_rows = combine(grouped) do sdf
        idx = argmin(sdf.EL_Rank)
          raw_desc = String(sdf.Peptide_label[idx])
          cleaned = replace(raw_desc, r"_(C|V|A|D)$" => "")
          cleaned = replace(cleaned, r"_\d+$" => "")
          cleaned = replace(cleaned, ' ' => '_')
          # Frame is the last part after underscores
          parts = split(cleaned, "_")
          frame = parts[end]
          (; Best_EL_Rank = sdf.EL_Rank[idx],
              Frame = frame,
              Sequence = sdf.Peptide[idx])
    end
    return best_rows
end

# Find best ranks separately for ancestral (_A, legacy _C) and derived (_D, legacy _V) peptides
    println("Calculating best ranks for ancestral peptides (_A)...")
    best_C = find_best_ranks(df, "_A")
    best_C.Peptide_Type .= "A"
    println("Found best ranks for ancestral peptides: $(nrow(best_C)) entries")

    println("Calculating best ranks for derived peptides (_D)...")
    best_V = find_best_ranks(df, "_D")
    best_V.Peptide_Type .= "D"
    println("Found best ranks for derived peptides: $(nrow(best_V)) entries")


# Combine ancestral and derived results
    best_ranks = vcat(best_C, best_V)

# --- Map Locus to protein Description from frames.csv ---
    # Context:
    # Previously, `Frame` was parsed from `Peptide_label` by taking the
    # last underscore-delimited token, which truncated protein names.
    # We now replace `Frame` using frames.csv by matching each row's
    # `Locus` to the Region ranges and taking the corresponding full
    # `Description`. If no match is found, we keep the original value.
    frames_file = resolve_read(joinpath(folder_path, "frames.csv"); suffix="", latest=latest)
    if isfile(frames_file)
        frames_df = CSV.read(frames_file, DataFrame)

        # Parse a Region string like "77,496" or multiple coords "77,496;606,980"
        # into a vector of (start, end) tuples
        function region_to_segments(region_str::AbstractString)::Vector{Tuple{Int,Int}}
            rs = strip(String(region_str))
            isempty(rs) && return Tuple{Int,Int}[]
            segments = Tuple{Int,Int}[]
            for part in split(rs, ';')
                se = split(part, ',')
                if length(se) == 2
                    start_pos = try
                        parse(Int, se[1])
                    catch
                        continue
                    end
                    end_pos = try
                        parse(Int, se[2])
                    catch
                        continue
                    end
                    push!(segments, (start_pos, end_pos))
                end
            end
            return segments
        end

        # Build a lookup table of regions → description
        regions_lookup = [(region_to_segments(fr.Region), String(fr.Description)) for fr in eachrow(frames_df)]

        # Replace Frame (derived from Peptide_label) with full Description matched by Locus.
        # If multiple frames overlap a locus, the first match in frames.csv is used.
        if !isempty(best_ranks)
            best_ranks[!, :Frame] = [
                begin
                    loc = r.Locus
                    # find first frames.csv entry whose any segment contains locus
                    desc = nothing
                    for (segments, d) in regions_lookup
                        for (s,e) in segments
                            if loc >= s && loc <= e
                                desc = d
                                break
                            end
                        end
                        if desc !== nothing
                            break
                        end
                    end
                    desc === nothing ? String(r.Frame) : desc
                end for r in eachrow(best_ranks)
            ]
        end
    else
        println("Warning: frames.csv not found in $folder_path. Protein labels will not be mapped.")
    end

# Reorder: put Frame, Locus, Mutation first
    best_ranks = select(best_ranks, :Frame, :Locus, :Mutation, Not([:Frame, :Locus, :Mutation]))

# Compute one Frame per (Locus, Mutation): prefer current mapped, otherwise first seen
    description_roots = DataFrame(Locus=Int[], Mutation=String[], Frame=String[])
    if !isempty(best_ranks)
        description_roots = combine(groupby(best_ranks, [:Locus, :Mutation])) do sdf
            (; Locus = sdf.Locus[1], Mutation = sdf.Mutation[1], Frame = sdf.Frame[1])
        end
    end

# Save best_ranks.csv
    best_ranks_file = resolve_write(joinpath(folder_path, "best_ranks.csv"); suffix=suffix)
    CSV.write(best_ranks_file, best_ranks)
    println("Saved best ranks to $best_ranks_file")

# Pivot best_ranks to have separate columns for HMBR_A and HMBR_D
    println("Calculating harmonic mean best ranks (HMBR) per Frame/Locus/Mutation...")
    if !isempty(best_ranks)
        pivot_df = unstack(
            combine(groupby(best_ranks, [:Frame, :Locus, :Mutation, :Peptide_Type]),
                    :Best_EL_Rank => harmmean => :HMBR),
            :Peptide_Type, :HMBR)
    # Only rename columns if they exist
    colnames = names(pivot_df)
    rename_pairs = Pair{Symbol,Symbol}[]
    if "A" in colnames
        push!(rename_pairs, Symbol("A") => :HMBR_A)
    end
    if "D" in colnames
        push!(rename_pairs, Symbol("D") => :HMBR_D)
    end
    if !isempty(rename_pairs)
        rename!(pivot_df, rename_pairs...)
    end

    # Warn if no derived data is present
    if !("HMBR_D" in names(pivot_df))
        println("Warning: No derived peptides found. Skipping fold change calculations and HMBR_D output.")
    end

    # Identify and report missing values before fold change calculation (per Locus, Mutation)
    missing_msgs = Set{Tuple{Int,String,String}}()  # (Locus, Mutation, which)
    for row in eachrow(pivot_df)
        locus = row.Locus
        change = get(row, :Mutation, "?")
        if ismissing(get(row, :HMBR_A, missing))
            push!(missing_msgs, (locus, String(change), "ancestral"))
        end
        if ismissing(get(row, :HMBR_D, missing))
            push!(missing_msgs, (locus, String(change), "derived"))
        end
    end
    for (locus, change, which) in missing_msgs
        println("Fold change could not be calculated for locus $(locus) (change $(change)) due to missing $(which) rank.")
    end

    # Filter out loci where both HMBR_A and HMBR_D are greater than 2 (non-binding)
    before_filter = nrow(pivot_df)
    pivot_df = filter(row -> !ismissing(get(row, :HMBR_A, missing)) && !ismissing(get(row, :HMBR_D, missing)) && !(get(row, :HMBR_A, 0.0) > 2 && get(row, :HMBR_D, 0.0) > 2), pivot_df)
    removed_count = before_filter - nrow(pivot_df)
    println("Removed $removed_count loci where both ancestral and derived states were predicted to be non-binding (HMBR > 2)")

    # Calculate fold change (Derived (_D) / Ancestral (_A)) for all valid rows
    if ("HMBR_A" in names(pivot_df)) && ("HMBR_D" in names(pivot_df))
        pivot_df.foldchange_HMBR = pivot_df.HMBR_D ./ pivot_df.HMBR_A
        pivot_df.log2_foldchange_HMBR = log2.(pivot_df.foldchange_HMBR)
    else
        println("DEBUG: Cannot calculate foldchange_HMBR or log2_foldchange_HMBR due to missing columns.")
    end

    # Frame is already present from grouping; no join needed

    # Keep Description exactly as mapped; ensure no trailing underscores or numeric suffixes
    # No further cleanup needed for Frame

    # Calculate fold change and log2 after all joins/cleaning
    if ("HMBR_A" in names(pivot_df)) && ("HMBR_D" in names(pivot_df))
        pivot_df.foldchange_HMBR = pivot_df.HMBR_D ./ pivot_df.HMBR_A
        pivot_df.log2_foldchange_HMBR = log2.(pivot_df.foldchange_HMBR)
    end

    # Prepare final output columns (use already computed values)
    output_cols = [:Frame, :Locus, :Mutation, :HMBR_A, :HMBR_D, :foldchange_HMBR, :log2_foldchange_HMBR]
    pivot_df = select(pivot_df, output_cols...)

    # Save harmonic mean results with fold change
    harmonic_mean_file = resolve_write(joinpath(folder_path, "harmonic_mean_best_ranks.csv"); suffix=suffix)
    CSV.write(harmonic_mean_file, pivot_df)
    println("Saved harmonic mean best ranks to $harmonic_mean_file")

    # Minimal test: write Locus and log2_foldchange_HMBR to a separate CSV for debugging
    if :log2_foldchange_HMBR in names(pivot_df)
        minimal_test_file = resolve_write(joinpath(folder_path, "harmonic_mean_best_ranks_log2_test.csv"); suffix=suffix)
        minimal_df = select(pivot_df, :Locus, :log2_foldchange_HMBR)
        CSV.write(minimal_test_file, minimal_df)
        println("Minimal test CSV written to $minimal_test_file")
    else
    end
    else
        println("No valid best rank data available. Skipping harmonic mean calculations.")
    end
end