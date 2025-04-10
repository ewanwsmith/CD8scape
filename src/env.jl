using Pkg

# Activate the environment in the current directory silently
Pkg.activate(".", io = devnull)

# Install RCall first, silently
Pkg.add("RCall"; io = devnull)

# Load RCall after installation
using RCall

# List of other Julia packages to install silently
packages = [
    "DataFrames",
    "CSV",
    "FilePathsBase",
    "ArgParse",
    "CodecZlib",
    "Statistics",
    "StatsBase",
]

# Add remaining Julia packages without output
for pkg in packages
    Pkg.add(pkg; io = devnull)
end

println("Julia dependencies installed.")