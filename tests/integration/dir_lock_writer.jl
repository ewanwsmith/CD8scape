#!/usr/bin/env julia
# Simple writer that appends a CSV row to a shared master file using lock_acquire/lock_release
using CSV, DataFrames
include(joinpath(@__DIR__, "../../src/lock_helpers.jl"))

function append_row(master_file::String, id::String)
    lockpath = joinpath(dirname(master_file), ".peptide_cache_lock")
    got = lock_acquire_test(lockpath; timeout_seconds=30)
    if !got
        println("Writer $id: could not acquire lock")
        return false
    end
    try
        # use a simple atomic append helper to avoid heavy CSV read/concat in parallel
        line = "0,WRITER_$(id),TEST,HLA-A01:01,,,0.0,0.0\n"
        atomic_append_line(master_file, line)
        return true
    finally
        lock_release_test(lockpath)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 2
        println("Usage: dir_lock_writer.jl <master_file> <id>")
        exit(1)
    end
    master = ARGS[1]
    id = ARGS[2]
    ok = append_row(master, id)
    println("Writer $id done -> $ok")
end
