using Pkg
using TOML

# Activate the src environment explicitly (so prep works from any cwd)
Pkg.activate(@__DIR__, io = devnull)

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
    # Candidates: project root and src env
    root_proj = abspath(joinpath(@__DIR__, "..", "Project.toml"))
    src_proj  = abspath(joinpath(@__DIR__, "Project.toml"))
    for path in (root_proj, src_proj)
        if isfile(path)
            try
                tbl = TOML.parsefile(path)
                if haskey(tbl, "deps")
                    for (name, _) in tbl["deps"]
                        push!(dep_names, String(name))
                    end
                end
            catch
                # ignore parse errors and continue
            end
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
    # If manifest is missing, instantiate then retry
    Pkg.instantiate(; io = devnull)
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