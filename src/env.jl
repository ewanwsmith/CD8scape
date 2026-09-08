using Pkg
using TOML

# ---------------------------------------------------------------------------
# Project environment
# ---------------------------------------------------------------------------
# Activate the *repository root* explicitly, so prep works from any cwd.
#
# This must be the root, not src/: CD8scape.jl launches every worker script
# with `--project=.` from the repo root, so the root environment is the one the
# pipeline actually runs in. Activating src/ here meant `prep` installed and
# instantiated a second, parallel environment that nothing else ever used.
const PROJECT_DIR = normpath(joinpath(@__DIR__, ".."))

"""
    stale_manifest(dir) -> Union{String,Nothing}

Path of `dir`'s Manifest.toml if it was resolved by a different Julia minor
series, otherwise `nothing`.

A manifest written by e.g. Julia 1.13 records stdlibs that do not exist in
1.11 (JuliaSyntaxHighlighting, ...). Pkg then aborts with "Could not locate
the source code for the X package" and cannot recover on its own.
"""
function stale_manifest(dir::AbstractString)
    path = joinpath(dir, "Manifest.toml")
    isfile(path) || return nothing
    recorded = try
        get(TOML.parsefile(path), "julia_version", nothing)
    catch
        return nothing   # unparseable — let Pkg report it in its own words
    end
    recorded === nothing && return nothing
    v = try
        VersionNumber(recorded)
    catch
        return nothing
    end
    return (v.major, v.minor) == (VERSION.major, VERSION.minor) ? nothing : path
end

# A manifest from another Julia series cannot be instantiated here. Move it
# aside (rather than deleting it) and let Pkg re-resolve from Project.toml.
let stale = stale_manifest(PROJECT_DIR)
    if stale !== nothing
        backup = stale * ".incompatible"
        mv(stale, backup; force = true)
        @warn "Manifest.toml was resolved by a different Julia version; re-resolving from Project.toml. Previous file kept at $backup."
    end
end

Pkg.activate(PROJECT_DIR, io = devnull)

# Ensure the environment has a valid Project/Manifest
try
    Pkg.resolve(; io = devnull)
catch
    # Ignore resolve issues; we'll instantiate/add below
end

import Base.Filesystem: isexecutable

"""
Collect dependency names from top-level `Project.toml` and `src/Project.toml`.
Falls back gracefully if a file is missing.
"""
function collect_declared_deps()
    dep_names = Set{String}()
    path = joinpath(PROJECT_DIR, "Project.toml")
    if isfile(path)
        try
            tbl = TOML.parsefile(path)
            if haskey(tbl, "deps")
                for (name, _) in tbl["deps"]
                    push!(dep_names, String(name))
                end
            end
        catch
            # ignore parse errors and fall back to the minimal set below
        end
    end
    # If no files found or empty, use minimal required set for this pipeline
    if isempty(dep_names)
        dep_names = Set([
            "DataFrames",
            "CSV",
            "FilePathsBase",
            "ArgParse",
            "CodecZlib",
            "Statistics",
            "StatsBase",
        ])
    end
    return collect(dep_names)
end

installed = try
    keys(Pkg.dependencies())
catch
    # If the manifest is missing or unusable, instantiate and retry. Both steps
    # must be guarded: a broken manifest makes instantiate itself throw, and an
    # uncaught error here aborts prep before any package can be added.
    try
        Pkg.instantiate(; io = devnull)
    catch err
        @warn "Could not instantiate the existing environment; adding packages from scratch." exception = err
    end
    try
        keys(Pkg.dependencies())
    catch
        # Fallback: treat as empty set
        Base.KeySet{String, Dict{String, Pkg.Types.PackageEntry}}(Dict{String, Pkg.Types.PackageEntry}())
    end
end

required_packages = collect_declared_deps()
for pkg in required_packages
    if !(pkg in installed)
        Pkg.add(pkg; io = devnull)
    end
end

# Instantiate the environment (fetches deps listed in Project/Manifest.toml)
Pkg.instantiate(io = devnull)

println("Julia dependencies installed and environment instantiated.")

# Check for netMHCpan executable
settings_path = joinpath(@__DIR__, "settings.txt")

function get_netmhcpan_path(settings_file::String)
    # Read NETMHCPAN from settings.txt
    # Accept either "NETMHCPAN=/full/path/to/netMHCpan" or a single bare line with the path.
    netmhcpan_path = nothing
    for line in eachline(settings_file)
        s = strip(line)
        if isempty(s) || startswith(s, "#")
            continue
        end
        if occursin('=', s)
            k, v = strip.(split(s, '=', limit=2))
            if uppercase(k) == "NETMHCPAN"
                netmhcpan_path = v
                break
            end
        else
            # Bare path fallback
            netmhcpan_path = s
            break
        end
    end
    if netmhcpan_path === nothing
        error("NETMHCPAN not set in settings.txt. Add a line like: NETMHCPAN=/full/path/to/netMHCpan")
    end
    # Expand ~ and make absolute relative to settings file dir if needed
    netmhcpan_path = replace(netmhcpan_path, "~" => homedir())
    if !isabspath(netmhcpan_path)
        netmhcpan_path = abspath(joinpath(dirname(settings_file), netmhcpan_path))
    end
    return netmhcpan_path
end

function check_netmhcpan(exec_path::AbstractString)
    if !isfile(exec_path)
        error("netMHCpan executable not found at $exec_path")
    elseif !Base.Filesystem.isexecutable(exec_path)
        error("netMHCpan found at $exec_path, but it is not executable")
    end
    return exec_path
end

# ---------------------------------------------------------------------------
# Python UI dependencies
# ---------------------------------------------------------------------------

const PYTHON_UI_PACKAGES = [
    ("PyQt6",        "PyQt6>=6.4"),
    ("pyqtgraph",    "pyqtgraph>=0.13"),
    ("numpy",        "numpy>=1.24"),
    ("setproctitle", "setproctitle"),
]

# Test/dev dependencies — installed by prep but not required at runtime
const PYTHON_DEV_PACKAGES = [
    ("pytest", "pytest"),
]

"""
Find a Python 3.8+ executable, mirroring the priority order in launch.command.
Returns the path as a String, or nothing if not found.
"""
function find_python()
    candidates = String[]

    # 1. Active conda environment
    conda_prefix = get(ENV, "CONDA_PREFIX", "")
    if !isempty(conda_prefix)
        push!(candidates, joinpath(conda_prefix, "bin", "python3"))
    end

    # 2. Common conda / Homebrew / system locations
    append!(candidates, [
        "/opt/anaconda3/bin/python3",
        "/opt/miniconda3/bin/python3",
        joinpath(homedir(), "anaconda3", "bin", "python3"),
        joinpath(homedir(), "miniconda3", "bin", "python3"),
        "/opt/homebrew/bin/python3",
        "/usr/local/bin/python3",
        "python3",
    ])

    for c in candidates
        resolved = Sys.which(c)
        resolved === nothing && !isabspath(c) && continue
        path = resolved !== nothing ? resolved : c
        isfile(path) || continue
        try
            ok = readchomp(`$path -c "import sys; print(sys.version_info >= (3,8))"`)
            ok == "True" && return path
        catch
            continue
        end
    end
    return nothing
end

"""
Ensure each required Python UI package is importable, installing via pip if not.
Skips silently if no Python 3.8+ is found (UI is optional for the pipeline).
"""
function ensure_python_ui_deps()
    python = find_python()
    if python === nothing
        @warn "Python 3.8+ not found — skipping UI dependency check. The Julia pipeline will still work."
        return
    end

    for (mod, spec) in vcat(PYTHON_UI_PACKAGES, PYTHON_DEV_PACKAGES)
        importable = try
            ret = run(ignorestatus(`$python -c "import $mod"`))
            success(ret)
        catch
            false
        end

        if !importable
            println("Python: $mod not found — installing $spec ...")
            try
                run(`$python -m pip install $spec --break-system-packages`)
                println("Python: $mod installed.")
            catch e
                @warn "Could not install $spec: $e"
            end
        else
            println("Python: $mod OK")
        end
    end
end

ensure_python_ui_deps()

# Check for Perl installation
function check_perl_installed()
    perl_path = Sys.which("perl")
    if perl_path === nothing
        error("Perl is not installed or not found in PATH.")
    end
end

try
    netmhcpan_path = get_netmhcpan_path(settings_path)
    netmhcpan_exec = check_netmhcpan(netmhcpan_path)
    ENV["NETMHCPAN"] = netmhcpan_exec
    println("netMHCpan executable found and is executable: ", netmhcpan_exec)

    check_perl_installed()
    println("Perl installation found.")

    println("Prep stage finished successfully")
catch e
    println("Error during preparation stage: ", e)
    rethrow(e)
end