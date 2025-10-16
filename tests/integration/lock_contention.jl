using Test, CSV, DataFrames, Random

function make_writer_script(path::AbstractString)
        open(path, "w") do io
                script = raw"""#!/usr/bin/env bash
set -euo pipefail
master="$1"
pep="$2"
id="$3"
done="$4"
lockdir="${master}.dirlock"
# Acquire directory lock (mkdir is atomic on POSIX)
start=$(date +%s.%N)
while ! mkdir "$lockdir" 2>/dev/null; do
    now=$(date +%s.%N)
    elapsed=$(echo "$now - $start" | bc -l)
    if (( $(echo "$elapsed > 10" | bc -l) )); then echo "lock timeout" >&2; exit 2; fi
    sleep 0.01
done
echo "$id $(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$lockdir/owner.txt"
dir=$(dirname "$master")
base=$(basename "$master")
tmp=$(mktemp -p "$dir" .${base}.tmp.XXXXXX)
# write header and append existing rows (excluding an existing header if present)
printf 'Peptide\n' > "$tmp"
if [ -f "$master" ]; then
    tail -n +2 "$master" >> "$tmp" || true
fi
echo "$pep" >> "$tmp"
# dedupe while preserving header
(head -n1 "$tmp"; tail -n +2 "$tmp" | sort -u) > "${tmp}.dedup"
mv -f "${tmp}.dedup" "$tmp"
mv -f "$tmp" "$master"
rm -rf "$lockdir"
touch "$done"
"""
                println(io, script)
        end
        chmod(path, 0o755)
end

@testset "master cache lock crash-simulation (Julia writers)" begin
    tmp = Base.mktempdir()
    master = joinpath(tmp, "peptide_netmhc_processed.csv")

    jwriter = joinpath(tmp, "jwriter.jl")
    open(jwriter, "w") do io
        script = raw"""#!/usr/bin/env julia
    using CSV, DataFrames
    import Dates

    # portable fileno helper
    function safe_fileno(io::IO)
        try
            return Base.Libc.fileno(io)
        catch
            try
                return Libc.fileno(io)
            catch
                throw(ErrorException("fileno-unavailable"))
            end
        end
    end

    function lock_acquire(lockfile::AbstractString; timeout_seconds::Int=30)
    # Try to use libc flock via ccall (so child process doesn't need Libc import).
    LOCK_SH = 1
    LOCK_EX = 2
    LOCK_NB = 4
    LOCK_UN = 8
    try
    io = open(lockfile, "w+")
    fd = safe_fileno(io)
        start = time()
        while true
            res = try
                # try common libc names (Linux/macOS)
                ccall((:flock, "libc"), Cint, (Cint, Cint), fd, LOCK_EX | LOCK_NB)
            catch
                try
                    ccall((:flock, "libSystem.B.dylib"), Cint, (Cint, Cint), fd, LOCK_EX | LOCK_NB)
                catch
                    close(io)
                    throw(ErrorException("flock-unavailable"))
                end
            end
            if res == 0
                return (true, io)
            end
            if time() - start > timeout_seconds
                close(io)
                return (false, nothing)
            end
            sleep(0.01)
        end
    catch
        # fallback: directory-based lock with stale-lock detection
        lockdir = lockfile * ".dirlock"
        start = time()
        while true
            try
                mkdir(lockdir)
                open(joinpath(lockdir, "owner.txt"), "w") do out
                    println(out, getpid(), " ", Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS"))
                end
                return (true, nothing)
            catch
                # If the lockdir exists, try to detect whether owner is still running.
                try
                    ownerfile = joinpath(lockdir, "owner.txt")
                    if isfile(ownerfile)
                        txt = read(ownerfile, String)
                        parts = split(strip(txt))
                        if !isempty(parts)
                            ownerpid = try parse(Int, parts[1]) catch; nothing end
                            if ownerpid !== nothing
                                # check if process exists (signal 0)
                                exists = try
                                    kill(ownerpid, 0)
                                    true
                                catch
                                    false
                                end
                                if !exists
                                    # stale lock: remove and retry immediately
                                    try Base.rm(lockdir; force=true, recursive=true) catch end
                                    continue
                                end
                            end
                        end
                    end
                catch
                end
                if time() - start > timeout_seconds
                    return (false, nothing)
                end
                sleep(0.01)
            end
        end
    end
end

function lock_release(lockfile::AbstractString, handle)
    try
        if handle !== nothing
            # handle is IO returned by open(lockfile)
            fd = safe_fileno(handle)
            try
                Base.Libc.flock(fd, Base.Libc.LOCK_UN)
            catch
            end
            try close(handle) catch end
            try Base.rm(lockfile) catch end
            return true
        else
            lockdir = lockfile * ".dirlock"
            try Base.rm(lockdir; recursive=true) catch end
            return true
        end
    catch
        return false
    end
end

function atomic_write_csv(path::AbstractString, lines::Vector{String})
    dir = dirname(path)
    base = basename(path)
    tmp = joinpath(dir, "." * base * ".tmp" * string(rand(UInt32)))
    open(tmp, "w") do o
        println(o, "Peptide")
        for l in lines
            println(o, l)
        end
    end
    # debug: report tmp used
    try
        println("atomic_write_csv: wrote tmp=", tmp)
        println("atomic_write_csv: lines=", join(lines, ","))
    catch
    end
    mv(tmp, path; force=true)
end

if length(ARGS) < 4
    println("Usage: jwriter.jl <master_csv> <peptide> <id> <mode>")
    exit(1)
end
master = ARGS[1]
pep = ARGS[2]
id = ARGS[3]
mode = ARGS[4]

lockfile = joinpath(dirname(master), ".peptide_cache_lock")
got, handle = lock_acquire(lockfile; timeout_seconds=10)
if !got
    println("JWriter $id: could not acquire lock")
    exit(2)
end
if handle !== nothing
    println("JWriter $id: acquired flock")
else
    println("JWriter $id: acquired dirlock")
end

lines = String[]
    # Append-only write: create header if missing and append peptide as a new row.
    try
        if !isfile(master)
            open(master, "w") do o
                println(o, "Peptide")
            end
        end
        open(master, "a") do o
            println(o, pep)
        end
        try println("JWriter $id: appended pep=", pep) catch end
    catch e
        try println(stderr, "append-write error: ", e) catch end
        exit(3)
    end

if mode == "crash"
    # simulate abrupt exit without calling lock_release; process exit should release resources
    exit(0)
end

lock_release(lockfile, handle)
println("JWriter $id: done")
"""
        println(io, script)
    end
    chmod(jwriter, 0o755)

    # Spawn several writers, some of which will crash
    peptides = ["J_A","J_B","J_C","J_D","J_E"]
    modes = ["ok","crash","ok","crash","ok"]
    done_markers = String[]
    for i in 1:length(peptides)
        done = joinpath(tmp, "jwriter_done_$(i).txt")
        push!(done_markers, done)
        out = joinpath(tmp, "jwriter_out_$(i).txt")
        err = joinpath(tmp, "jwriter_err_$(i).txt")
    cmd = `$(Sys.BINDIR)/julia $jwriter $master $(peptides[i]) $(i) $(modes[i])`
        @async begin
            open(out, "w") do o
                open(err, "w") do e
                    try
                        run(pipeline(cmd; stdout=o, stderr=e))
                    catch e
                        write(e, "ERROR: $e\n")
                    end
                end
            end
            touch(done)
        end
    end

    # wait for writers to finish
    deadline = time() + 15.0
    while time() < deadline
        if all(isfile, done_markers)
            break
        end
        sleep(0.05)
    end
    if !all(isfile, done_markers)
        println("Timeout waiting for jwriters. tmp: $tmp")
        try
            for f in readdir(tmp; join=true)
                println("  -> ", f)
            end
        catch
        end
        for i in 1:length(peptides)
            out = joinpath(tmp, "jwriter_out_$(i).txt")
            err = joinpath(tmp, "jwriter_err_$(i).txt")
            println("---- jwriter $(i) stderr ----")
            if isfile(err)
                try println(read(err, String)) catch end
            else
                println("(no stderr)")
            end
            println("---- jwriter $(i) stdout ----")
            if isfile(out)
                try println(read(out, String)) catch end
            else
                println("(no stdout)")
            end
        end
    end

        # On failure, print helpful debug logs. Otherwise keep output minimal.
        if !(all(isfile, done_markers) && isfile(master))
            println("--- debug: tmp dir: $tmp ---")
            try
                for f in readdir(tmp; join=true)
                    println("  -> ", f)
                end
            catch
            end
            println("--- master file content ---")
            if isfile(master)
                try
                    println(read(master, String))
                catch
                    println("(could not read master)")
                end
            else
                println("(master missing)")
            end
            for i in 1:length(peptides)
                out = joinpath(tmp, "jwriter_out_$(i).txt")
                err = joinpath(tmp, "jwriter_err_$(i).txt")
                println("---- jwriter $(i) stderr ----")
                if isfile(err)
                    try println(read(err, String)) catch end
                else
                    println("(no stderr)")
                end
                println("---- jwriter $(i) stdout ----")
                if isfile(out)
                    try println(read(out, String)) catch end
                else
                    println("(no stdout)")
                end
            end
        end

        @test all(isfile, done_markers)
        @test isfile(master)
        df = CSV.read(master, DataFrame)
        found = sort(unique(String.(df.Peptide)))
        expected = sort(unique(peptides))
        @test found == expected

    Base.rm(tmp; recursive=true)
end

@testset "master cache lock contention" begin
    tmp = mktempdir()
    master = joinpath(tmp, "peptide_netmhc_processed.csv")

    writer = joinpath(tmp, "writer.jl")
    make_writer_script(writer)

    # define peptides (include a duplicate to test deduplication)
    peptides = ["PEP_A", "PEP_B", "PEP_C", "PEP_D", "PEP_A", "PEP_E", "PEP_F", "PEP_G"]

    # spawn background writers that contend on the same master cache
    done_markers = String[]
    tasks = Task[]
    for (i, pep) in enumerate(peptides)
        done = joinpath(tmp, "writer_done_$(i).txt")
        push!(done_markers, done)
        writer_stdout = joinpath(tmp, "writer_stdout_$(i).txt")
        writer_stderr = joinpath(tmp, "writer_stderr_$(i).txt")
        # open files for redirection
        outio = open(writer_stdout, "w")
        errio = open(writer_stderr, "w")
    cmd = `bash $writer $master $pep $i $done`
        t = @async begin
            try
                run(pipeline(cmd; stdout=outio, stderr=errio))
            catch e
                # record failure to stderr file
                write(errio, "ERROR: $e\n")
            finally
                try close(outio) catch end
                try close(errio) catch end
            end
        end
        push!(tasks, t)
    end

    # give processes a short moment to start
    sleep(0.05)

    # wait for all done markers with a timeout
    deadline = time() + 30.0
    while time() < deadline
        if all(isfile, done_markers)
            break
        end
        sleep(0.05)
    end

    if !all(isfile, done_markers)
        println("Timeout waiting for writers. tmp dir: $tmp")
        try
            for f in readdir(tmp; join=true)
                println("  -> ", f)
            end
        catch
        end
        # print per-writer stderr/stdout for debugging
        for (i, _) in enumerate(peptides)
            sfile = joinpath(tmp, "writer_stderr_$(i).txt")
            ofile = joinpath(tmp, "writer_stdout_$(i).txt")
            println("---- writer $(i) stderr ----")
            if isfile(sfile)
                try println(read(sfile, String)) catch end
            else
                println("(no stderr file yet)")
            end
            println("---- writer $(i) stdout ----")
            if isfile(ofile)
                try println(read(ofile, String)) catch end
            else
                println("(no stdout file yet)")
            end
        end
        @test false
    end

    # read master cache and verify deduplicated peptides present
    @test isfile(master)
    df = CSV.read(master, DataFrame)
    found = sort(unique(String.(df.Peptide)))
    expected = sort(unique(peptides))
    @test found == expected

    # cleanup
    Base.rm(tmp; recursive=true)
end
