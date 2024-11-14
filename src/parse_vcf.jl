#!/usr/bin/env julia

include("./env.jl")  # Include the env.jl file

using ArgParse
using CSV
using DataFrames
using CodecZlib
import Base.Filesystem: splitext

function parse_command_line()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--in"
        help = "Path to the .vcf or .vcf.gz file"
        arg_type = String
        required = true
    end
    return parse_args(s)
end

function read_vcf_with_csv(vcf_filename)
    # Open the VCF file (compressed or uncompressed)
    file_io = open_vcf_file(vcf_filename)

    # Read lines from the VCF file
    lines = readlines(file_io)
    close(file_io)

    # Find the header line starting with "#CHROM"
    header_index = findfirst(line -> startswith(line, "#CHROM"), lines)
    if header_index === nothing
        error("Header line starting with '#CHROM' not found in the VCF file.")
    end

    # Extract the header line and data lines
    header_line = lines[header_index]
    data_lines = lines[(header_index + 1):end]

    # Prepare the data for CSV parsing
    data_text = join(data_lines, "\n")
    data_io = IOBuffer(data_text)

    # Parse the header to remove the leading '#'
    header = split(header_line[2:end], '\t')
    header = String.(header)  # Convert SubString{String} to String

    # Read the data using CSV.jl
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
    # Check if the file extension is .gz
    ext = get_file_extension(filename)
    if ext == ".gz"
        return true
    end
    # Alternatively, check the magic bytes
    open(filename, "r") do io
        magic_bytes = read(io, 2)  # Read first 2 bytes
        return magic_bytes == UInt8[0x1f, 0x8b]
    end
end

function get_file_extension(filename)
    _, ext = splitext(filename)
    return ext
end

function main()
    args = parse_command_line()
    vcf_file_path = args["in"]

    println("Processing file: $vcf_file_path")

    vcf_dataframe = read_vcf_with_csv(vcf_file_path)
    println("First 5 rows of the DataFrame:")
    display(first(vcf_dataframe, 5))
end

main()