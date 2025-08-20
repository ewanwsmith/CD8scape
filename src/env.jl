using Pkg

# Activate the environment in the current directory silently
Pkg.activate(".", io = devnull)

import Base.Filesystem: isexecutable

# Ensure that Project.toml exists and includes needed packages
# Only needed the first time you run the script or when updating dependencies
required_packages = [
    "DataFrames",
    "CSV",
    "FilePathsBase",
    "ArgParse",
    "CodecZlib",
    "Statistics",
    "StatsBase",
]

installed = keys(Pkg.dependencies())

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