using Pkg

# Activate the environment in the current directory silently
Pkg.activate(".", io = devnull)

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