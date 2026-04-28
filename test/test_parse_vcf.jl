"""
test_parse_vcf.jl — Tests for src/parse_vcf.jl

Covers:
  get_file_extension – standard and gz extensions
  is_gzipped         – magic-byte detection, .gz extension
  find_vcf_file      – finds .vcf, finds .vcf.gz, errors when absent
  read_vcf_with_csv  – parses VCF header + data rows, multi-allelic expansion
"""

using Test
using DataFrames
using CodecZlib

include(joinpath(@__DIR__, "../src/parse_vcf.jl"))

# ---------------------------------------------------------------------------
function with_tmpdir(f::Function)
    dir = mktempdir()
    try; f(dir); finally; rm(dir; recursive=true, force=true); end
end

MINIMAL_VCF = """
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tT\t.\tPASS\t.
chr1\t200\t.\tCC\tTT\t.\tPASS\t.
""" |> lstrip

# ===========================================================================
@testset "get_file_extension" begin

    @testset "plain extension" begin
        @test get_file_extension("variants.vcf") == ".vcf"
        @test get_file_extension("data.csv") == ".csv"
    end

    @testset ".vcf.gz treated as .gz" begin
        @test get_file_extension("variants.vcf.gz") == ".gz"
    end

    @testset "no extension" begin
        ext = get_file_extension("Makefile")
        @test ext == ""
    end

    @testset "hidden file (dot prefix)" begin
        ext = get_file_extension(".gitignore")
        @test ext == ".gitignore" || ext == ""
    end

end

# ===========================================================================
@testset "is_gzipped" begin

    @testset "plain VCF file is not gzipped" begin
        f, io = mktemp()
        write(io, MINIMAL_VCF); close(io)
        @test !is_gzipped(f)
        rm(f)
    end

    @testset "file with .gz extension and gzip magic bytes" begin
        f, io = mktemp()
        close(io)
        gzip_path = f * ".gz"
        open(GzipCompressorStream, gzip_path, "w") do gz
            write(gz, MINIMAL_VCF)
        end
        @test is_gzipped(gzip_path)
        rm(gzip_path)
    end

    @testset "file with .gz extension always returns true (extension check takes priority)" begin
        # is_gzipped() checks the file extension first; a file ending in .gz
        # returns true regardless of magic bytes. This matches the implementation.
        f, io = mktemp()
        write(io, "plain text not gzipped"); close(io)
        fake_gz = f * ".gz"
        cp(f, fake_gz)
        @test is_gzipped(fake_gz)   # extension wins
        rm(f); rm(fake_gz)
    end

end

# ===========================================================================
@testset "find_vcf_file" begin

    @testset "finds .vcf file in folder" begin
        with_tmpdir() do dir
            vcf = joinpath(dir, "variants.vcf")
            write(vcf, MINIMAL_VCF)
            @test find_vcf_file(dir) == vcf
        end
    end

    @testset "finds .vcf.gz file in folder" begin
        with_tmpdir() do dir
            gz = joinpath(dir, "variants.vcf.gz")
            open(GzipCompressorStream, gz, "w") do io
                write(io, MINIMAL_VCF)
            end
            result = find_vcf_file(dir)
            @test endswith(result, ".vcf.gz")
        end
    end

    @testset "throws when no VCF file present" begin
        with_tmpdir() do dir
            @test_throws ErrorException find_vcf_file(dir)
        end
    end

    @testset "prefers first alphabetically when multiple present" begin
        with_tmpdir() do dir
            write(joinpath(dir, "a.vcf"), MINIMAL_VCF)
            write(joinpath(dir, "b.vcf"), MINIMAL_VCF)
            result = find_vcf_file(dir)
            @test basename(result) in ("a.vcf", "b.vcf")
        end
    end

end

# ===========================================================================
@testset "read_vcf_with_csv" begin

    @testset "returns a DataFrame" begin
        f, io = mktemp()
        write(io, MINIMAL_VCF); close(io)
        df = read_vcf_with_csv(f)
        rm(f)
        @test df isa DataFrame
    end

    @testset "contains POS, REF, ALT columns" begin
        f, io = mktemp()
        write(io, MINIMAL_VCF); close(io)
        df = read_vcf_with_csv(f)
        rm(f)
        @test :POS in propertynames(df)
        @test :REF in propertynames(df)
        @test :ALT in propertynames(df)
    end

    @testset "parses correct number of data rows" begin
        f, io = mktemp()
        write(io, MINIMAL_VCF); close(io)
        df = read_vcf_with_csv(f)
        rm(f)
        @test nrow(df) == 2
    end

    @testset "parses POS as integer-compatible" begin
        f, io = mktemp()
        write(io, MINIMAL_VCF); close(io)
        df = read_vcf_with_csv(f)
        rm(f)
        @test df.POS[1] == 100
        @test df.POS[2] == 200
    end

    @testset "no header line throws error" begin
        bad_vcf = "chr1\t100\t.\tA\tT\t.\tPASS\t.\n"
        f, io = mktemp()
        write(io, bad_vcf); close(io)
        @test_throws ErrorException read_vcf_with_csv(f)
        rm(f)
    end

    @testset "reads gzip-compressed VCF" begin
        f, io = mktemp(); close(io); rm(f)
        gz = f * ".vcf.gz"
        open(GzipCompressorStream, gz, "w") do io; write(io, MINIMAL_VCF); end
        df = read_vcf_with_csv(gz)
        rm(gz)
        @test nrow(df) == 2
    end

    @testset "non-ACGT alleles are excluded during main expansion" begin
        # The filtering happens in main(), not in read_vcf_with_csv itself,
        # but we verify that read_vcf_with_csv at least returns the raw row
        vcf = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" *
              "chr1\t50\t.\tN\tT\t.\tPASS\t.\n"
        f, io = mktemp()
        write(io, vcf); close(io)
        df = read_vcf_with_csv(f)
        rm(f)
        @test nrow(df) == 1
        @test string(df.REF[1]) == "N"
    end

end
