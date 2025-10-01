#!/usr/bin/env julia
"""
read_ncbi_frames.jl

Reads NCBI reference frames from sequences.fasta for open reading frame extraction.

Usage:
    julia read_ncbi_frames.jl <folder_path>

Arguments:
    <folder_path>   Path to the folder containing sequences.fasta.
"""

"""
Usage:
    julia script_name.jl <folder_path>

This script:
1) Finds sequences.(fa|fasta) for "reference reading frames".
2) Finds consensus.(fa|fasta) for the "consensus sequence".
3) Parses each header in sequences-file for Region & Description.
4) Reads the entire sequence from the consensus-file.
5) Extracts each Region's subsequence from the consensus, storing it as `Consensus_sequence`.
6) Saves a single `frames.csv` with columns: Region, Consensus_sequence, Description.
"""

using DataFrames
using CSV
using FilePathsBase
import Base.Filesystem: joinpath, abspath, isfile, findfirst

# ──────────────────────────────────────────────────────────────────────────────
# Parse the FASTA header, extracting Region (like "77,496" or "77,496;606,980") 
# and Description (any leftover text).
function parse_header(header::AbstractString)
    # Remove '>'
    header = strip(header, ['>'])

    region_str  = ""
    description = ""

    # If header has "join(" syntax, e.g. >join(X:77..496,Y:606..980) ...
    if occursin("join(", header)
        join_part, rest_of_header = split(header, ')'; limit=2)
        raw_coords = join_part[6:end]  # remove 'join('
        rest_of_header = strip(rest_of_header)

        # raw_coords might look like: "acc1:77..496,acc2:606..980"
        region_parts = split(raw_coords, ',')
        formatted = String[]
        for region in region_parts
            acc_and_coords = split(region, ':')
            coords = length(acc_and_coords) == 2 ? acc_and_coords[2] : acc_and_coords[1]
            if occursin("..", coords)
                start_pos, end_pos = split(coords, "..")
                push!(formatted, "$start_pos,$end_pos")
            else
                push!(formatted, "$coords,$coords")
            end
        end
        region_str = join(formatted, ';')
        description = rest_of_header
    else
        # Typical header: >accession:77..496 Description text...
        parts = split(header, ' '; limit=2)
        accession_coords = split(parts[1], ':')

        if length(accession_coords) == 2
            coords = accession_coords[2]
            if occursin("..", coords)
                start_pos, end_pos = split(coords, "..")
                region_str = "$start_pos,$end_pos"
            else
                region_str = "$coords,$coords"
            end
        else
            region_str = ""
        end

        description = length(parts) > 1 ? parts[2] : ""
    end

    # Clean description
    description = lstrip(description, '|')
    description = replace(description, r"\[.*?\]" => "")
    description = strip(description)

    return region_str, description
end

# ──────────────────────────────────────────────────────────────────────────────
# Reads the FASTA file (headers only). Builds a DataFrame with columns:
#   Region, Description
function read_fasta_metadata(file_path::String)::DataFrame
    df = DataFrame(Region=String[], Description=String[])
    open(file_path, "r") do io
        header = ""
        for line in eachline(io)
            line = strip(line)
            if isempty(line)
                continue
            elseif startswith(line, '>')
                # If there's a previous header, parse it first
                if !isempty(header)
                    reg_str, desc = parse_header(header)
                    push!(df, (reg_str, desc))
                end
                header = line
            else
                # Sequence lines are ignored here
                continue
            end
        end
        # Parse the last header if any
        if !isempty(header)
            reg_str, desc = parse_header(header)
            push!(df, (reg_str, desc))
        end
    end
    return df
end

# ──────────────────────────────────────────────────────────────────────────────
# Reads a single-sequence FASTA file (like `Consensus.fa`), returning the
# entire consensus as a single string.
function read_fasta_sequence(file_path::String)::String
    seq = ""
    open(file_path, "r") do io
        for line in eachline(io)
            line = strip(line)
            if isempty(line) || startswith(line, '>')
                continue
            end
            seq *= line
        end
    end
    return seq
end

# ──────────────────────────────────────────────────────────────────────────────
# Given a region string like "77,496" or multiple coords "77,496;606,980"
# plus the consensus sequence, return the concatenated subsequence(s).
function extract_regions_sequence(region_str::String, consensus::String)::String
    if isempty(region_str)
        return ""
    end

    regions = split(region_str, ';')
    subseqs = String[]
    for region in regions
        start_end = split(region, ',')
        if length(start_end) != 2
            error("Invalid region format: $region")
        end
        start_pos = parse(Int, start_end[1])
        end_pos   = parse(Int, start_end[2])
        if start_pos < 1 || end_pos > length(consensus)
            error("Region out of bounds: $start_pos,$end_pos for consensus length $(length(consensus))")
        end
        push!(subseqs, consensus[start_pos:end_pos])
    end
    return join(subseqs, "")
end

# ──────────────────────────────────────────────────────────────────────────────
function main()
    if length(ARGS) < 1
        println("Usage: julia script_name.jl <folder_path>")
        return
    end

    folder_path = abspath(ARGS[1])

    # Attempt to find sequences.(fa|fasta) - case insensitive
    function find_file_case_insensitive(folder, basename, extensions)
        for ext in extensions
            # Check all possible case combinations
            patterns = [
                lowercase(basename) * "." * lowercase(ext),
                lowercase(basename) * "." * uppercase(ext),
                uppercase(basename) * "." * lowercase(ext),
                uppercase(basename) * "." * uppercase(ext),
                titlecase(basename) * "." * lowercase(ext),
                titlecase(basename) * "." * uppercase(ext)
            ]
            for pattern in patterns
                filepath = joinpath(folder, pattern)
                if isfile(filepath)
                    return filepath
                end
            end
        end
        return nothing
    end
    
    sequences_fa = find_file_case_insensitive(folder_path, "sequences", ["fa", "fasta"])
    if sequences_fa === nothing
        error("No sequences.fa or sequences.fasta file found in $folder_path")
    end

    # Attempt to find consensus.(fa|fasta) - case insensitive
    consensus_fa = find_file_case_insensitive(folder_path, "consensus", ["fa", "fasta"])
    if consensus_fa === nothing
        error("No consensus.fa or consensus.fasta file found in $folder_path")
    end

    println("Using sequences file: $sequences_fa")
    println("Using consensus file: $consensus_fa")

    # 1) Parse "sequences" headers → DataFrame(Region, Description)
    df = read_fasta_metadata(sequences_fa)
    
    # 2) Read entire consensus sequence from "consensus"
    consensus_seq = read_fasta_sequence(consensus_fa)

    # 3) For each row, extract the corresponding subsequence from the consensus
    #    and store it in a new column :Consensus_sequence
    df[!, :Consensus_sequence] = [extract_regions_sequence(df.Region[i], consensus_seq) for i in 1:nrow(df)]

    # We want columns: Region, Consensus_sequence, Description
    # So reorder them as desired
    final_df = select(df, :Region, :Consensus_sequence, :Description)

    # 4) Write final DataFrame to frames.csv
    output_csv = joinpath(folder_path, "frames.csv")
    CSV.write(output_csv, final_df)
    println("DataFrame written to $output_csv with columns: Region, Consensus_sequence, Description")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end