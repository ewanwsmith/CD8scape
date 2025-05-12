#!/usr/bin/env julia

using CSV
using DataFrames
using CodecZlib
import Base.Filesystem: splitext, dirname, basename, joinpath, readdir

function main()
    # Ensure the user provided a folder path
    if length(ARGS) < 1
        println("Usage: julia script.jl <folder_path>")
        return
    end

    folder_path = ARGS[1]
    vcf_file_path = find_vcf_file(folder_path)
    println("Found VCF file: $vcf_file_path")

    # Read the VCF file into a DataFrame
    vcf_dataframe = read_vcf_with_csv(vcf_file_path)

    # Keep only POS, REF, ALT columns
    select!(vcf_dataframe, [:POS, :REF, :ALT])

    # Rename columns:
    #   POS  -> Locus
    #   REF  -> Consensus
    #   ALT  -> Variant
    rename!(vcf_dataframe, :POS => :Locus, :REF => :Consensus, :ALT => :Variant)

    # Write out the CSV named 'variants.csv' in the same folder
    output_file_path = joinpath(folder_path, "variants.csv")
    CSV.write(output_file_path, vcf_dataframe)
    println("DataFrame saved as $output_file_path")
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