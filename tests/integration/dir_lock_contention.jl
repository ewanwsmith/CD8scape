#!/usr/bin/env julia
using Test, CSV, DataFrames
import Random

@testset "dir_lock_contention" begin
    N_WRITERS = 16

    tmpdir = Base.mktempdir()
    master = joinpath(tmpdir, "peptide_netmhc_processed.csv")

    # ensure environment uses directory-lock path
    ENV["CD8S_USE_FLOCK"] = "0"

    writer_script = joinpath(@__DIR__, "dir_lock_writer.jl")

    tasks = Vector{Task}(undef, N_WRITERS)
    failures = String[]
    for i in 1:N_WRITERS
        id = string(i)
        cmd = `julia $writer_script $master $id`
        tasks[i] = @async begin
            try
                run(cmd)
            catch e
                push!(failures, "writer task $id failed: $(e)")
            end
        end
    end

    # wait for tasks to finish (with timeout)
    deadline = time() + 30.0
    for t in tasks
        while !istaskdone(t) && time() < deadline
            sleep(0.05)
        end
        if !istaskdone(t)
            push!(failures, "A writer task did not finish in time")
        end
    end

    @test isfile(master)
    df = CSV.read(master, DataFrame)
    @test nrow(df) == N_WRITERS

    # Fail if any writers recorded an error (show failures for debugging)
    if !isempty(failures)
        @info "Writer failures:" failures
    end
    @test isempty(failures)
end
