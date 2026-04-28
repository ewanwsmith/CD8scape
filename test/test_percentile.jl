"""
test_percentile.jl — Tests for helper functions in src/percentile.jl

Covers the module-level helpers:
  _resolve_path       – absolute vs. relative path resolution
  _discover_observed  – finds latest non-simulated CSV
  _suffix_from_observed – extracts suffix from observed filename
  _probit             – inverse normal CDF (Abramowitz & Stegun approx)
  _normcdf            – normal CDF (rational approximation)
"""

using Test

# percentile.jl wraps its main() body in `if abspath(PROGRAM_FILE) == @__FILE__`
# so the helpers at module level are safe to include.
include(joinpath(@__DIR__, "../src/percentile.jl"))

# ---------------------------------------------------------------------------
function with_tmpdir(f::Function)
    dir = mktempdir()
    try; f(dir); finally; rm(dir; recursive=true, force=true); end
end

# ===========================================================================
@testset "_probit" begin

    @testset "p=0.5 gives 0 (median)" begin
        @test abs(_probit(0.5)) < 1e-3
    end

    @testset "p=0.0 gives -Inf" begin
        @test _probit(0.0) == -Inf
    end

    @testset "p=1.0 gives +Inf" begin
        @test _probit(1.0) == Inf
    end

    @testset "symmetric: probit(1-p) == -probit(p)" begin
        for p in [0.1, 0.25, 0.75, 0.9]
            @test abs(_probit(1 - p) + _probit(p)) < 1e-6
        end
    end

    @testset "monotonically increasing" begin
        probs = [0.05, 0.25, 0.5, 0.75, 0.95]
        vals  = [_probit(p) for p in probs]
        @test issorted(vals)
    end

    @testset "p=0.975 ≈ 1.96 (95% CI half)" begin
        @test abs(_probit(0.975) - 1.96) < 0.01
    end

    @testset "p=0.841 ≈ 1.0 (one sigma)" begin
        @test abs(_probit(0.841) - 1.0) < 0.02
    end

    @testset "returns Float64" begin
        @test _probit(0.5) isa Float64
    end

end

# ===========================================================================
@testset "_normcdf" begin

    @testset "z=0 gives 0.5" begin
        @test abs(_normcdf(0.0) - 0.5) < 1e-4
    end

    @testset "large positive z approaches 1" begin
        @test _normcdf(10.0) > 0.999
    end

    @testset "large negative z approaches 0" begin
        @test _normcdf(-10.0) < 0.001
    end

    @testset "symmetric: normcdf(-z) == 1 - normcdf(z)" begin
        for z in [0.5, 1.0, 1.96, 2.576]
            @test abs(_normcdf(-z) - (1.0 - _normcdf(z))) < 1e-6
        end
    end

    @testset "z=1.96 ≈ 0.975" begin
        @test abs(_normcdf(1.96) - 0.975) < 0.001
    end

    @testset "z=1.645 ≈ 0.95" begin
        @test abs(_normcdf(1.645) - 0.95) < 0.002
    end

    @testset "monotonically increasing" begin
        zs = [-3.0, -1.0, 0.0, 1.0, 3.0]
        vals = [_normcdf(z) for z in zs]
        @test issorted(vals)
    end

    @testset "returns Float64" begin
        @test _normcdf(0.0) isa Float64
    end

    @testset "probit and normcdf are inverse functions (roundtrip)" begin
        for p in [0.05, 0.1, 0.5, 0.9, 0.95]
            z = _probit(p)
            p_back = _normcdf(z)
            @test abs(p_back - p) < 1e-4
        end
    end

end

# ===========================================================================
@testset "_resolve_path" begin

    @testset "absolute path returned unchanged" begin
        @test _resolve_path("/some/folder", "/abs/path.csv") == "/abs/path.csv"
    end

    @testset "relative path joined to folder" begin
        result = _resolve_path("/data/folder", "obs.csv")
        @test result == "/data/folder/obs.csv"
    end

    @testset "empty string returns empty string" begin
        @test _resolve_path("/data/folder", "") == ""
    end

end

# ===========================================================================
@testset "_discover_observed" begin

    @testset "finds the only non-simulated CSV" begin
        with_tmpdir() do dir
            write(joinpath(dir, "harmonic_mean_best_ranks.csv"), "a\n")
            write(joinpath(dir, "harmonic_mean_best_ranks_simulated.csv"), "b\n")
            result = _discover_observed(dir, "harmonic_mean_best_ranks")
            @test basename(result) == "harmonic_mean_best_ranks.csv"
        end
    end

    @testset "excludes _simulated suffix" begin
        with_tmpdir() do dir
            write(joinpath(dir, "harmonic_mean_best_ranks_simulated.csv"), "sim\n")
            @test_throws ErrorException _discover_observed(dir, "harmonic_mean_best_ranks")
        end
    end

    @testset "picks latest among multiple non-simulated" begin
        with_tmpdir() do dir
            f1 = joinpath(dir, "harmonic_mean_best_ranks.csv")
            write(f1, "old\n"); sleep(0.05)
            f2 = joinpath(dir, "harmonic_mean_best_ranks_v2.csv")
            write(f2, "new\n")
            result = _discover_observed(dir, "harmonic_mean_best_ranks")
            @test basename(result) == "harmonic_mean_best_ranks_v2.csv"
        end
    end

    @testset "throws when no candidates exist" begin
        with_tmpdir() do dir
            @test_throws ErrorException _discover_observed(dir, "harmonic_mean_best_ranks")
        end
    end

end

# ===========================================================================
@testset "_suffix_from_observed" begin

    @testset "no suffix returns empty string" begin
        @test _suffix_from_observed("/folder/harmonic_mean_best_ranks.csv",
                                    "harmonic_mean_best_ranks") == ""
    end

    @testset "extracts suffix correctly" begin
        @test _suffix_from_observed("/folder/harmonic_mean_best_ranks_myrun.csv",
                                    "harmonic_mean_best_ranks") == "myrun"
    end

    @testset "strips leading underscore from suffix" begin
        result = _suffix_from_observed("/folder/harmonic_mean_best_ranks_v2.csv",
                                       "harmonic_mean_best_ranks")
        @test !startswith(result, "_")
        @test result == "v2"
    end

end
