#!/usr/bin/env julia
"""
parse_vcf.jl

Parses .vcf or .vcf.gz file to extract variant information.

Usage:
    julia parse_vcf.jl <folder_path>

Arguments:
    <folder_path>   Path to the folder containing .vcf file.
"""

using CSV
using DataFrames
using CodecZlib
import Base.Filesystem: splitext, dirname, basename, joinpath, readdir
include("path_utils.jl")

function main()
    # Ensure the user provided a folder path
    if length(ARGS) < 1
        println("Usage: julia parse_vcf.jl <folder_path>")
        return
    end

    folder_path = ARGS[1]
    # Optional suffix for output variants; latest not used here
    suffix = ""
    if length(ARGS) >= 2
        i = 2
        while i <= length(ARGS)
            arg = ARGS[i]
            if arg == "--suffix"
                if i + 1 <= length(ARGS) && !startswith(ARGS[i+1], "--")
                    i += 1
                    suffix = ARGS[i]
                end
            end
            i += 1
        end
    end
    vcf_file_path = find_vcf_file(folder_path)
    println("Found VCF file: $vcf_file_path")

    # Read the VCF file into a DataFrame
    vcf_dataframe = read_vcf_with_csv(vcf_file_path)

    # Keep only POS, REF, ALT columns
    select!(vcf_dataframe, [:POS, :REF, :ALT])

    # Expand multi-allelic ALT (comma-separated) into one row per ALT.
    # Also filter to SNVs (single-base REF and ALT) to match downstream expectations.
    expanded = DataFrame(POS=Int[], REF=String[], ALT=String[])
    for r in eachrow(vcf_dataframe)
        ref = String(r[:REF])
        # ALT may be comma-separated
        alts = split(String(r[:ALT]), ",")
        for alt in alts
            alt_str = String(alt)
            # Keep only SNVs: single-base REF and ALT, A/C/G/T
            if length(ref) == 1 && length(alt_str) == 1 && occursin(r"^[ACGTacgt]$", ref) && occursin(r"^[ACGTacgt]$", alt_str)
                push!(expanded, (Int(r[:POS]), uppercase(ref), uppercase(alt_str)))
            end
        end
    end

    # Rename columns:
    #   POS  -> Locus
    #   REF  -> Consensus
    #   ALT  -> Variant
    rename!(expanded, :POS => :Locus, :REF => :Consensus, :ALT => :Variant)

    # Write out the CSV named 'variants.csv' in the same folder
    output_file_path = resolve_write(joinpath(folder_path, "variants.csv"); suffix=suffix)
    CSV.write(output_file_path, expanded)
    println("Variants written to: $output_file_path")
end

function find_vcf_file(folder_path::String)
    # Searches for a file that ends in .vcf or .vcf.gz in the given folder.
    files = readdir(folder_path)
    for file in files
        lowerfile = lowercase(file)
        if endswith(lowerfile, ".vcf") || endswith(lowerfile, ".vcf.gz")
            return joinpath(folder_path, file)
        end
    end
    error("No .vcf or .vcf.gz file found in $folder_path.")
end

function read_vcf_with_csv(vcf_filename)
    file_io = open_vcf_file(vcf_filename)
    lines = readlines(file_io)
    close(file_io)

    # Find header line index (i.e., line starting with '#CHROM')
    header_index = findfirst(line -> startswith(line, "#CHROM"), lines)
    if header_index === nothing
        error("Header line starting with '#CHROM' not found in the VCF file.")
    end

    # Extract header line and data lines
    header_line = lines[header_index]
    data_lines = lines[(header_index + 1):end]

    data_text = join(data_lines, "\n")
    data_io = IOBuffer(data_text)

    # If the header line looks like "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO..."
    # Then splitting on '\t' after removing '#':
    # ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", ...]
    header = split(header_line[2:end], '\t')
    header = String.(header)  # Convert SubString{String} to String

    df = CSV.read(data_io, DataFrame; delim='\t', header=header)
    return df
end

function open_vcf_file(filename)
    try
        if is_gzipped(filename)
            return GzipDecompressorStream(open(filename, "r"))
        else
            return open(filename, "r")
        end
    catch e
        error("Failed to open file: $filename. Error: $(e)")
    end
end

function is_gzipped(filename)
    ext = get_file_extension(filename)
    if ext == ".gz"
        return true
    end
    open(filename, "r") do io
        magic_bytes = read(io, 2)  # Read first 2 bytes
        return magic_bytes == UInt8[0x1f, 0x8b]
    end
end

function get_file_extension(filename)
    _, ext = splitext(filename)
    return ext
end

# Run main if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end