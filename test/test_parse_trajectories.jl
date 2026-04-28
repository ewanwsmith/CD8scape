"""
test_parse_trajectories.jl — Tests for src/parse_trajectories.jl

Covers:
  parse_trajectories – correct column schema, time-point columns,
                       multi-row files, missing time points filled with 0,
                       whitespace-only / empty lines skipped
"""

using Test
using DataFrames

include(joinpath(@__DIR__, "../src/parse_trajectories.jl"))

# ---------------------------------------------------------------------------
function write_traj(content::String)
    f, io = mktemp()
    write(io, content)
    close(io)
    return f
end

# ===========================================================================
@testset "parse_trajectories" begin

    @testset "single row, one time point" begin
        # Format: locus consensus variant num_timepoints time A C G T total
        f = write_traj("100 A G 1 0 50 0 30 0 80\n")
        df = parse_trajectories(f)
        rm(f)
        @test nrow(df) == 1
        @test df.Locus[1]   == 100
        @test df.Consensus[1] == "A"
        @test df.Variant[1]   == "G"
        @test hasproperty(df, :Count_A_t0)
        @test df.Count_A_t0[1] == 50
        @test df.Total_Reads_t0[1] == 80
    end

    @testset "two rows, same time point" begin
        content = "100 A G 1 10 10 0 20 0 30\n" *
                  "200 C T 1 10  0 5  0 5 10\n"
        f = write_traj(content)
        df = parse_trajectories(f)
        rm(f)
        @test nrow(df) == 2
        @test df.Locus[1] == 100
        @test df.Locus[2] == 200
    end

    @testset "multiple time points produce correct columns" begin
        # Two time points: t0 and t5
        f = write_traj("100 A G 2 0 10 0 20 0 30  5 15 0 25 0 40\n")
        df = parse_trajectories(f)
        rm(f)
        @test hasproperty(df, :Count_A_t0)
        @test hasproperty(df, :Count_A_t5)
        @test hasproperty(df, :Total_Reads_t5)
    end

    @testset "missing time point for a row is filled with 0" begin
        # Row 1 has t0 and t5; Row 2 has only t0
        content = "100 A G 2 0 10 0 20 0 30  5 15 0 25 0 40\n" *
                  "200 C T 1 0  5 0  5 0 10\n"
        f = write_traj(content)
        df = parse_trajectories(f)
        rm(f)
        @test hasproperty(df, :Count_A_t5)
        @test df[df.Locus .== 200, :Count_A_t5][1] == 0
    end

    @testset "empty lines are skipped" begin
        content = "100 A G 1 0 10 0 20 0 30\n\n\n200 C T 1 0 5 0 5 0 10\n"
        f = write_traj(content)
        df = parse_trajectories(f)
        rm(f)
        @test nrow(df) == 2
    end

    @testset "column order: Locus, Consensus, Variant first" begin
        f = write_traj("100 A G 1 0 10 0 20 0 30\n")
        df = parse_trajectories(f)
        rm(f)
        cols = names(df)
        @test cols[1] == "Locus"
        @test cols[2] == "Consensus"
        @test cols[3] == "Variant"
    end

    @testset "time points appear in ascending order" begin
        content = "100 A G 3 10 1 0 2 0 3  5 4 0 5 0 9  0 6 0 7 0 13\n"
        f = write_traj(content)
        df = parse_trajectories(f)
        rm(f)
        time_cols = [c for c in names(df) if startswith(c, "Count_A_t")]
        times = sort([parse(Int, replace(c, r"Count_A_t" => "")) for c in time_cols])
        @test times == sort(times)
    end

    @testset "result is a DataFrame" begin
        f = write_traj("100 A G 1 0 10 0 20 0 30\n")
        df = parse_trajectories(f)
        rm(f)
        @test df isa DataFrame
    end

    @testset "Locus column contains integers" begin
        f = write_traj("999 A T 1 0 1 2 3 4 10\n")
        df = parse_trajectories(f)
        rm(f)
        @test df.Locus[1] isa Integer
    end

end
