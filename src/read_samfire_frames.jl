#!/usr/bin/env julia
"""
read_samfire_frames.jl

Reads Samfire Reading_Frames.dat for open reading frame extraction.

Usage:
    julia read_samfire_frames.jl <folder_path>

Arguments:
    <folder_path>   Path to the folder containing Reading_Frames.dat.
"""

using DataFrames
using CSV
using FilePathsBase
using Logging
import Base.Filesystem: joinpath, dirname, isfile

"""
Reads `Reading_Frames.dat` where each block is 3 lines:
1) `start end`
2) Nucleotide consensus sequence
3) (ignored)

Generates a DataFrame with columns:
- Region            (e.g. "77,496")
- Consensus_sequence
- Description       (auto-labeled "Frame_1", "Frame_2", etc.)

Returns a DataFrame with rows equal to how many 3-line blocks it successfully parsed.
"""
function read_dat_file(filepath::String)::DataFrame
    lines = readlines(filepath)

    regions         = String[]
    consensus       = String[]
    descriptions    = String[]

    frame_counter   = 1
    i = 1
    n = length(lines)

    while i <= n
        # Ensure we have at least 2 lines for Region + Consensus
        # and 1 extra line to skip. So total 3 lines needed.
        if i + 2 > n
            @warn "Incomplete block starting at line $i. Skipping."
            break
        end

        # ─────────────────────────────────────────────────────
        # (1) Region line: parse two integers, e.g. "77 496"
        start_end = split(strip(lines[i]))
        if length(start_end) < 2
            @warn "Line $i does not contain two integers. Skipping this 3-line block."
            i += 3
            continue
        end

        start_val = tryparse(Int, start_end[1])
        end_val   = tryparse(Int, start_end[2])
        if isnothing(start_val) || isnothing(end_val)
            @warn "Invalid Start/End at line $i. Skipping this 3-line block."
            i += 3
            continue
        end

        region_str = "$(start_val),$(end_val)"
        push!(regions, region_str)

        # ─────────────────────────────────────────────────────
        # (2) Consensus line
        consensus_seq = strip(lines[i + 1])
        push!(consensus, consensus_seq)

        # ─────────────────────────────────────────────────────
        # (3) Skip the third line entirely (the amino-acid line)
        # but create a Description = "Frame_<n>"
        desc_str = "Frame_$(frame_counter)"
        push!(descriptions, desc_str)

        frame_counter += 1
        i += 3
    end

    frames_df = DataFrame(
        Region             = regions,
        Consensus_sequence = consensus,
        Description        = descriptions
    )

    return frames_df
end

"""
Writes the given DataFrame to `frames.csv` in the same directory as the `.dat` file.
"""
function save_dataframe_as_csv(df::DataFrame, input_filepath::String, output_filename::String="frames.csv")
    directory = dirname(input_filepath)
    output_filepath = joinpath(directory, output_filename)
    CSV.write(output_filepath, df)
    println("DataFrame successfully saved to $output_filepath")
end

function main()
    if length(ARGS) < 1
        println("Usage: julia script_name.jl <folder_path>")
        exit(1)
    end

    folder_path = ARGS[1]
    dat_filepath = joinpath(folder_path, "Reading_Frames.dat")

    if !isfile(dat_filepath)
        println("The file 'Reading_Frames.dat' was not found in '$folder_path'. Searching for 'sequences.fa' instead.")
        return  # Exits the function without erroring out
    end

    frames_df = read_dat_file(dat_filepath)
    println("Generated DataFrame:")
    display(frames_df)

    save_dataframe_as_csv(frames_df, dat_filepath, "frames.csv")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end