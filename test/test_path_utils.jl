"""
test_path_utils.jl — Comprehensive tests for src/path_utils.jl

Covers:
  with_suffix      – inserting suffix before extension
  _candidates_for  – discovering file candidates in a directory
  discover_path    – single/multi-candidate resolution (latest & strict)
  resolve_read     – read-path helper (suffix vs. discovery)
  resolve_write    – write-path builder (no existence check)
  get_env_int      – ENV integer reader with defaults and clamp
"""

using Test
using Random

include(joinpath(@__DIR__, "../src/path_utils.jl"))

# ---------------------------------------------------------------------------
# Helper: create a named temp dir that is cleaned up after the testset
# ---------------------------------------------------------------------------
function with_tmpdir(f::Function)
    dir = mktempdir()
    try
        f(dir)
    finally
        rm(dir; recursive=true, force=true)
    end
end

# ===========================================================================
@testset "with_suffix" begin

    @testset "empty suffix returns path unchanged" begin
        @test with_suffix("/foo/bar/frames.csv", "") == "/foo/bar/frames.csv"
        @test with_suffix("frames.csv", "") == "frames.csv"
        @test with_suffix("/a/b/c.txt", "") == "/a/b/c.txt"
    end

    @testset "inserts suffix before extension" begin
        @test with_suffix("frames.csv", "simulated") == "frames_simulated.csv"
        @test with_suffix("/data/variants.csv", "v2") == "/data/variants_v2.csv"
        @test with_suffix("results.tsv", "run1") == "results_run1.tsv"
    end

    @testset "works with multi-dot filenames" begin
        @test with_suffix("my.data.csv", "test") == "my.data_test.csv"
    end

    @testset "works with absolute paths" begin
        @test with_suffix("/home/user/data/frames.csv", "abc") == "/home/user/data/frames_abc.csv"
    end

    @testset "returns String type" begin
        @test with_suffix("frames.csv", "x") isa String
        @test with_suffix("frames.csv", "") isa String
    end

    @testset "suffix containing underscores" begin
        @test with_suffix("frames.csv", "run_1") == "frames_run_1.csv"
    end

    @testset "file with no extension" begin
        result = with_suffix("Makefile", "test")
        @test endswith(result, "_test")
    end

end

# ===========================================================================
@testset "discover_path" begin

    @testset "single exact file" begin
        with_tmpdir() do dir
            f = joinpath(dir, "frames.csv")
            write(f, "header\n")
            result = discover_path(f; latest=false)
            @test result == f
        end
    end

    @testset "single suffixed file" begin
        with_tmpdir() do dir
            f = joinpath(dir, "frames_simulated.csv")
            write(f, "header\n")
            base = joinpath(dir, "frames.csv")
            result = discover_path(base; latest=false)
            @test result == f
        end
    end

    @testset "no file throws error" begin
        with_tmpdir() do dir
            base = joinpath(dir, "frames.csv")
            @test_throws ErrorException discover_path(base; latest=false)
        end
    end

    @testset "multiple candidates with latest=true returns most-recent" begin
        with_tmpdir() do dir
            f1 = joinpath(dir, "frames.csv");        write(f1, "old\n")
            sleep(0.05)
            f2 = joinpath(dir, "frames_v2.csv");     write(f2, "newer\n")
            sleep(0.05)
            f3 = joinpath(dir, "frames_v3.csv");     write(f3, "newest\n")
            result = discover_path(joinpath(dir, "frames.csv"); latest=true)
            @test result == f3
        end
    end

    @testset "multiple candidates with latest=false throws error" begin
        with_tmpdir() do dir
            write(joinpath(dir, "frames.csv"),    "a\n")
            write(joinpath(dir, "frames_v2.csv"), "b\n")
            @test_throws ErrorException discover_path(joinpath(dir, "frames.csv"); latest=false)
        end
    end

    @testset "ignores unrelated files in directory" begin
        with_tmpdir() do dir
            write(joinpath(dir, "frames.csv"),    "target\n")
            write(joinpath(dir, "variants.csv"),  "other\n")
            write(joinpath(dir, "README.md"),     "doc\n")
            result = discover_path(joinpath(dir, "frames.csv"); latest=false)
            @test basename(result) == "frames.csv"
        end
    end

end

# ===========================================================================
@testset "resolve_read" begin

    @testset "with suffix — uses suffixed path when it exists" begin
        with_tmpdir() do dir
            f = joinpath(dir, "frames_simulated.csv")
            write(f, "x\n")
            result = resolve_read(joinpath(dir, "frames.csv"); suffix="simulated")
            @test result == f
        end
    end

    @testset "with suffix — throws when suffixed path missing" begin
        with_tmpdir() do dir
            @test_throws ErrorException resolve_read(
                joinpath(dir, "frames.csv"); suffix="ghost")
        end
    end

    @testset "without suffix — discovers existing file" begin
        with_tmpdir() do dir
            f = joinpath(dir, "variants.csv")
            write(f, "x\n")
            result = resolve_read(f; suffix="")
            @test result == f
        end
    end

    @testset "without suffix — auto-discovers latest when latest=true" begin
        with_tmpdir() do dir
            f1 = joinpath(dir, "variants.csv");       write(f1, "old\n")
            sleep(0.05)
            f2 = joinpath(dir, "variants_new.csv");   write(f2, "new\n")
            result = resolve_read(joinpath(dir, "variants.csv"); suffix="", latest=true)
            @test result == f2
        end
    end

    @testset "returns String" begin
        with_tmpdir() do dir
            f = joinpath(dir, "frames.csv")
            write(f, "x\n")
            @test resolve_read(f) isa String
        end
    end

end

# ===========================================================================
@testset "resolve_write" begin

    @testset "empty suffix returns base path unchanged" begin
        @test resolve_write("/data/frames.csv"; suffix="") == "/data/frames.csv"
    end

    @testset "inserts suffix like with_suffix" begin
        r = resolve_write("/data/variants.csv"; suffix="simulated")
        @test r == "/data/variants_simulated.csv"
    end

    @testset "does not require file to exist" begin
        # Should not throw even though file doesn't exist
        result = resolve_write("/nonexistent/dir/file.csv"; suffix="x")
        @test endswith(result, "_x.csv")
    end

    @testset "returns String" begin
        @test resolve_write("/tmp/out.csv") isa String
    end

end

# ===========================================================================
@testset "get_env_int" begin

    @testset "missing key returns default" begin
        key = "CDTEST_MISSING_$(randstring(8))"
        @test get_env_int(key, 4) == 4
        @test get_env_int(key, 1) == 1
        @test get_env_int(key, 100) == 100
    end

    @testset "valid integer from ENV" begin
        key = "CDTEST_INT_$(randstring(8))"
        ENV[key] = "8"
        try
            @test get_env_int(key, 1) == 8
        finally
            delete!(ENV, key)
        end
    end

    @testset "non-numeric ENV value returns default" begin
        key = "CDTEST_BAD_$(randstring(8))"
        ENV[key] = "notanumber"
        try
            @test get_env_int(key, 3) == 3
        finally
            delete!(ENV, key)
        end
    end

    @testset "clamps to min_val" begin
        key = "CDTEST_CLAMP_$(randstring(8))"
        ENV[key] = "0"
        try
            @test get_env_int(key, 1; min_val=1) == 1
        finally
            delete!(ENV, key)
        end
        ENV[key] = "-5"
        try
            @test get_env_int(key, 1; min_val=1) == 1
        finally
            delete!(ENV, key)
        end
    end

    @testset "value above min_val passes through" begin
        key = "CDTEST_ABOVE_$(randstring(8))"
        ENV[key] = "16"
        try
            @test get_env_int(key, 1; min_val=1) == 16
        finally
            delete!(ENV, key)
        end
    end

end
