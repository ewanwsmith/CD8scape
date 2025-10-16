#!/usr/bin/env julia
using Dates

# Lightweight locking helpers used by tests. Keep dependencies minimal.
const _LOCK_HANDLES_TEST = Dict{String, IO}()

# Try to detect flock support (simple, non-invasive)
function flock_supported_test()
    try
        if !isdefined(Libc, :flock)
            return false
        end
    catch
        return false
    end
    tmp = nothing
    try
        tmp = open("/dev/null", "w")
        fd = nothing
        try
            fd = Libc.fileno(tmp)
        catch
            try
                fd = Base.Libc.fileno(tmp)
            catch
                fd = nothing
            end
        end
        if fd === nothing
            return false
        end
        try
            res = Libc.flock(fd, Libc.LOCK_EX | Libc.LOCK_NB)
            if res == 0
                try Libc.flock(fd, Libc.LOCK_UN) catch end
                return true
            else
                return false
            end
        catch
            return false
        end
    catch
        return false
    finally
        try if tmp !== nothing close(tmp) end catch end
    end
end

function lock_acquire_test(lockfile::AbstractString; timeout_seconds::Int=30)
    # default to directory lock for portability
    lockdir = lockfile * ".dirlock"
    start = time()
    while true
        try
            mkdir(lockdir)
            try
                open(joinpath(lockdir, "owner.txt"), "w") do out
                    println(out, getpid(), " ", Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS"))
                end
            catch
            end
            return true
        catch
            if time() - start > timeout_seconds
                return false
            end
            sleep(0.05)
        end
    end
end

function lock_release_test(lockfile::AbstractString)
    lockdir = lockfile * ".dirlock"
    try
        Base.rm(lockdir; recursive=true)
        return true
    catch
        return false
    end
end

# Atomic append a single CSV line using mktemp+rename. Line should include newline.
function atomic_append_line(master::AbstractString, line::AbstractString)
    dir = dirname(master)
    mkpath(dir)
    # Under directory-lock, it's safe to open the file in append mode and write
    need_header = !isfile(master)
    open(master, "a") do out
        if need_header
            println(out, "Pos,Peptide,ID,HLA,core,icore,EL-score,EL_Rank")
        end
        print(out, line)
    end
    return true
end

export lock_acquire_test, lock_release_test, atomic_append_line, flock_supported_test
