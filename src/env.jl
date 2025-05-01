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
    for line in eachline(settings_file)
        line = strip(line)
        if isempty(line) || startswith(line, "#")
            continue
        end
        return line  # Assume the first valid line is the path
    end
    error("netMHCpan path not found in settings.txt")
end

function check_netmhcpan(path::AbstractString)
    netmhcpan_exec = joinpath(path, "netMHCpan")
    if !isfile(netmhcpan_exec)
        error("netMHCpan executable not found at $netmhcpan_exec")
    elseif !Base.Filesystem.isexecutable(netmhcpan_exec)
        error("netMHCpan found at $netmhcpan_exec, but it is not executable")
    end
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
    check_netmhcpan(netmhcpan_path)
    println("netMHCpan executable found and is executable.")

    check_perl_installed()
    println("Perl installation found.")

    println("Prep stage finished successfully")
catch e
    println("Error during preparation stage: ", e)
    rethrow(e)
end