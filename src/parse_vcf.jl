using CSV
using DataFrames
using CodecZlib
import Base.Filesystem: splitext, dirname, basename, joinpath

function read_vcf_with_csv(vcf_filename)
    file_io = open_vcf_file(vcf_filename)
    lines = readlines(file_io)
    close(file_io)

    header_index = findfirst(line -> startswith(line, "#CHROM"), lines)
    if header_index === nothing
        error("Header line starting with '#CHROM' not found in the VCF file.")
    end

    header_line = lines[header_index]
    data_lines = lines[(header_index + 1):end]

    data_text = join(data_lines, "\n")
    data_io = IOBuffer(data_text)

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

function get_base_filename(filepath)
    filename = basename(filepath)
    if endswith(filename, ".vcf.gz")
        base_filename = filename[1:end - length(".vcf.gz")]
    elseif endswith(filename, ".vcf")
        base_filename = filename[1:end - length(".vcf")]
    else
        base_filename, _ = splitext(filename)
    end
    return base_filename
end

function get_output_file_path(vcf_file_path)
    dir = dirname(vcf_file_path)
    base_filename = get_base_filename(vcf_file_path)
    output_filename = base_filename * ".csv"
    return joinpath(dir, output_filename)
end

function process_vcf_to_csv(vcf_file_path)
    println("Processing file: $vcf_file_path")

    vcf_dataframe = read_vcf_with_csv(vcf_file_path)
    println("First 5 rows of the DataFrame:")
    display(first(vcf_dataframe, 5))

    output_file_path = get_output_file_path(vcf_file_path)
    CSV.write(output_file_path, vcf_dataframe)
    println("DataFrame saved to $output_file_path")
end

function main()
    # Set the VCF file path here
    vcf_file_path = "/Users/e.smith.5/Documents/PhD/CD8scape/data/Flu_example/E3i_06R_030_B_18V18332.vcf"  # <-- Replace with your file path

    # Call the processing function
    process_vcf_to_csv(vcf_file_path)
end

main()